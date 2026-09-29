'''
create_redmapper_dp2_specz.py

Spec-z training sample for redMaPPer on DP2, from DESI DR1. Script version of
redmapper_dp2_specz.ipynb and redmapper_dp2_specz_red.ipynb.

    download  RSP only: pull "DESI DR1 lite within Rubin DP2 footprint" from the Photo-z Server
    select    clean the spec-z, cut to the geometry mask, pick the training pixels, write the spec file
    red       match to the catalog, keep the red-sequence ones, write the red spec file + training tracts

At NERSC there is no download, select reads DESI's own zall-pix-iron.fits (what the
Photo-z Server product was cut from) straight off CFS.
'''

import sys; import os
import glob
import zipfile
import argparse

import numpy as np
import pandas as pd
import healpy as hp
import healsparse as hsp
import fitsio

from astropy.coordinates import SkyCoord
import astropy.units as u

import redmapper_dp2_common as common
import create_redmapper_dp2_catalogs as rm_cats

'''
Initial definitions
'''
## DESI DR1 redshift catalog at NERSC
DESI_NERSC = '/global/cfs/cdirs/desi/public/dr1/spectro/redux/iron/zcatalog/v1/zall-pix-iron.fits'
DESI_PRODUCT = 'DESI DR1 lite within Rubin DP2 footprint'
DESI_COLUMNS = ['TARGETID', 'SURVEY', 'PROGRAM', 'TARGET_RA', 'TARGET_DEC',
                'Z', 'ZERR', 'ZWARN', 'DELTACHI2', 'ZCAT_PRIMARY', 'SPECTYPE']

## rows per read of the 28M row NERSC file
CHUNK = 2_000_000

## number of DP2 tracts with an object table, for labeling partial red runs
N_TRACTS_TOTAL = 2174


'''
file names
'''
def clean_file(specz_dir, bands, vis_min):
    '''
    Full cleaned in-footprint table, for plots and checks.
    '''
    return os.path.join(specz_dir, f'desi_dr1_clean_{bands}_vis{vis_min}.parquet')


def spec_file(specz_dir, bands, vis_min, red=False, partial=False):
    '''
    The spec file redmapper reads (ra, dec, z, z_err).
    '''
    kind = 'red_' if red else ''
    tag = '_partial' if partial else ''
    return os.path.join(specz_dir, f'dp2_specz_desi_dr1_{kind}{bands}_vis{vis_min}{tag}.fit')


def hpix_file(specz_dir, bands, vis_min, nside_train):
    '''
    RING pixels of the training region.
    '''
    return os.path.join(specz_dir, f'training_hpix_{bands}_vis{vis_min}_nside{nside_train}.txt')


'''
download (RSP)
'''
def run_download(args):
    '''
    Download step: the Photo-z Server's DESI DR1 cut to DP2, into ``--specz-dir``.

    Needs a token in ~/.pzserver_token (make one at pzserver.linea.org.br).

    Parameters
    ----------
    args: `argparse.Namespace`
       Parsed arguments.
    '''
    from pzserver import PzServer

    with open(os.path.expanduser('~/.pzserver_token')) as f:
        pz = PzServer(token=f.read().strip())

    ## finding the product by its title
    products = pz.get_products_list()
    if isinstance(products, pd.DataFrame):
        products = products.to_dict('records')
    matches = [p for p in products if any(DESI_PRODUCT in str(v) for v in p.values())]
    if len(matches) != 1:
        sys.exit(f'{len(matches)} products match "{DESI_PRODUCT}", expected 1')

    ## only downloading once
    if len(glob.glob(f'{args.specz_dir}/**/*.parquet', recursive=True)) == 0:
        pz.download_product(matches[0]['internal_name'], save_in=args.specz_dir)

    ## the download can come as a zip
    for zip_path in glob.glob(f'{args.specz_dir}/*.zip'):
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(args.specz_dir)
        print(f'unzipped {zip_path}')
    print('Done!')


'''
select
'''
def find_desi_file(args):
    '''
    DESI file to read: ``--desi-file``, else a downloaded parquet in ``--specz-dir``, else the NERSC one.

    Parameters
    ----------
    args: `argparse.Namespace`
       Parsed arguments.

    Returns
    -------
    path: `str`
       File path.
    '''
    if args.desi_file is not None:
        return args.desi_file
    downloaded = glob.glob(f'{args.specz_dir}/**/desi_dr1_in_rubin_dp2*.parquet', recursive=True)
    if len(downloaded) > 0:
        return downloaded[0]
    if os.path.exists(DESI_NERSC):
        return DESI_NERSC
    sys.exit('no DESI file: run the download step (RSP) or pass --desi-file')


def clean_chunk(desi, geometry):
    '''
    Quality + footprint cuts on one chunk of DESI.

    Parameters
    ----------
    desi: `pd.DataFrame`
       DESI rows with ``DESI_COLUMNS``.
    geometry: `healsparse.HealSparseMap`
       Fracgood geometry mask.

    Returns
    -------
    specz: `pd.DataFrame`
       Rows passing every cut.
    counts: `list` [`int`]
       Rows left after each cut, in the order the log prints them.
    '''
    ## string columns can come back as bytes or padded
    for col in ['SURVEY', 'PROGRAM', 'SPECTYPE']:
        if desi[col].map(type).eq(bytes).any():
            desi[col] = desi[col].str.decode('utf-8')
        desi[col] = desi[col].astype(str).str.strip()

    ## one redshift per target, no fit warnings, galaxies only
    primary = desi['ZCAT_PRIMARY'].astype(bool).values
    no_warn = (desi['ZWARN'] == 0).values
    is_galaxy = (desi['SPECTYPE'] == 'GALAXY').values
    good = primary & no_warn & is_galaxy

    in_footprint = np.zeros(len(desi), dtype=bool)
    in_footprint[good] = geometry.get_values_pos(desi['TARGET_RA'].values[good],
                                                 desi['TARGET_DEC'].values[good], lonlat=True) > 0

    counts = [len(desi), primary.sum(), (primary & no_warn).sum(), good.sum(), in_footprint.sum()]
    return desi[in_footprint].reset_index(drop=True), counts


def read_clean_desi(path, geometry):
    '''
    Read DESI and apply the cuts, chunked for the full NERSC file.

    Parameters
    ----------
    path: `str`
       Parquet (Photo-z Server download) or FITS (zall-pix-iron.fits).
    geometry: `healsparse.HealSparseMap`
       Fracgood geometry mask.

    Returns
    -------
    specz: `pd.DataFrame`
       Clean in-footprint galaxy redshifts.
    '''
    names = ['all rows', '+ ZCAT_PRIMARY', '+ ZWARN == 0', '+ GALAXY', '+ inside footprint']
    total = np.zeros(len(names), dtype=int)
    kept = []

    if path.endswith('.parquet'):
        specz, counts = clean_chunk(pd.read_parquet(path, columns=DESI_COLUMNS), geometry)
        kept.append(specz)
        total += counts
    else:
        with fitsio.FITS(path) as f:
            n_rows = f[1].get_nrows()
            for start in range(0, n_rows, CHUNK):
                rows = f[1][DESI_COLUMNS][start:start + CHUNK]
                desi = pd.DataFrame({c: rows[c].astype(rows[c].dtype.newbyteorder('='))
                                     for c in DESI_COLUMNS})
                specz, counts = clean_chunk(desi, geometry)
                kept.append(specz)
                total += counts
                print(f'  [{min(start + CHUNK, n_rows)}/{n_rows}] {total[-1]} in footprint so far', flush=True)

    for name, count in zip(names, total):
        print(f'{name:<20} {count}')
    return pd.concat(kept, ignore_index=True)


def z_bins(args):
    '''
    Redshift bin edges for the per-bin counts.
    '''
    return np.arange(args.zmin, args.zmax + 0.001, args.dz)


def training_cells(specz, geometry, depth, args):
    '''
    Area, depth, and spec-z counts for every nside_train RING pixel on the footprint.

    Parameters
    ----------
    specz: `pd.DataFrame`
       Clean in-footprint spec-z.
    geometry: `healsparse.HealSparseMap`
       Fracgood geometry mask.
    depth: `healsparse.HealSparseMap`
       Depth map.
    args: `argparse.Namespace`
       Parsed arguments.

    Returns
    -------
    cells: `pd.DataFrame`
       One row per candidate pixel.
    bin_counts: `pd.DataFrame`
       Usable spec-z per pixel (rows) and redshift bin (columns).
    survey_median: `float`
       Median limmag over the footprint.
    '''
    nside = geometry.nside_sparse
    geom_pix = geometry.valid_pixels
    geom_ra, geom_dec = common.pixel_centers(nside, geom_pix)
    geom_hpix = hp.ang2pix(args.nside_train, geom_ra, geom_dec, lonlat=True, nest=False)   ## RING, for the config
    limmag = depth.get_values_pix(geom_pix)['limmag'].astype(np.float64)   ## FITS is big-endian, pandas won't index it

    ## footprint area and depth per pixel
    cells = pd.DataFrame({'hpix': geom_hpix, 'limmag': limmag})
    cells = cells.groupby('hpix').agg(n_pix=('limmag', 'size'), median_limmag=('limmag', 'median'))
    cells['area'] = cells['n_pix'] * common.pix_area(nside)
    cells['filled'] = cells['area'] / hp.nside2pixarea(args.nside_train, degrees=True)

    ## spec-z per pixel, overall and per program
    specz_hpix = hp.ang2pix(args.nside_train, specz['TARGET_RA'].values, specz['TARGET_DEC'].values,
                            lonlat=True, nest=False)
    cells['n_specz'] = pd.Series(specz_hpix).value_counts().reindex(cells.index, fill_value=0)
    for program in sorted(specz['PROGRAM'].unique()):
        this_program = (specz['PROGRAM'] == program).values
        cells[f'n_{program}'] = pd.Series(specz_hpix[this_program]).value_counts().reindex(cells.index, fill_value=0)
    cells['center_ra'], cells['center_dec'] = hp.pix2ang(args.nside_train, cells.index.values, lonlat=True, nest=False)

    ## spec-z per redshift bin, counting only what redmapper will use
    z_edges = z_bins(args)
    n_bins = len(z_edges) - 1
    usable = (specz['ZERR'] < args.zerr_max).values
    usable_bin = np.digitize(specz['Z'].values[usable], z_edges) - 1
    in_range = (usable_bin >= 0) & (usable_bin < n_bins)
    bin_counts = pd.crosstab(specz_hpix[usable][in_range], usable_bin[in_range])
    bin_counts = bin_counts.reindex(index=cells.index, columns=range(n_bins), fill_value=0)

    return cells, bin_counts, float(np.median(limmag))


def choose_region(cells, bin_counts, survey_median, args):
    '''
    Grow a connected ~``--target-area`` region of typical-depth pixels.

    Starts at the pixel whose weakest redshift bin is best, then keeps adding the
    touching pixel that raises the weakest bin the most (ties go to more spec-z).
    Deep fields (like COSMOS) are skipped by only taking pixels within
    ``--depth-tol`` mag of the survey median depth.

    Parameters
    ----------
    cells: `pd.DataFrame`
       From `training_cells`.
    bin_counts: `pd.DataFrame`
       From `training_cells`.
    survey_median: `float`
       Median limmag over the footprint.
    args: `argparse.Namespace`
       Parsed arguments.

    Returns
    -------
    region: `list` [`int`]
       RING pixels.
    '''
    typical = (cells['median_limmag'] - survey_median).abs() < args.depth_tol
    candidates = cells[typical & (cells['n_specz'] > 0)]
    print(f'{len(candidates)} candidate pixels')
    if len(candidates) == 0:
        sys.exit('no typical-depth pixels with spec-z, pass --training-hpix by hand')

    def weakest_bin(hpix_list):
        return bin_counts.loc[hpix_list].sum().min()

    region = [int(bin_counts.loc[candidates.index].min(axis=1).idxmax())]
    while cells.loc[region, 'area'].sum() < args.target_area:
        touching = set()
        for h in region:
            touching.update(hp.get_all_neighbours(args.nside_train, h, nest=False))
        touching = [int(h) for h in touching if h in candidates.index and h not in region]
        if len(touching) == 0:
            print('ran out of touching candidates')
            break
        best = max(touching, key=lambda h: (weakest_bin(region + [h]), cells.loc[h, 'n_specz']))
        region.append(best)
    return region


def run_select(args):
    '''
    Select step: clean DESI, cut to the mask, pick the training pixels, write the spec file.

    The spec file doesn't depend on which training pixels were picked, only on the
    geometry mask. redmapper only uses spec-z inside the training pixels (plus
    border), matches them to our galaxies itself, and picks out the red ones
    itself, so this hands it every good galaxy redshift in the footprint.

    Parameters
    ----------
    args: `argparse.Namespace`
       Parsed arguments.
    '''
    geometry = hsp.HealSparseMap.read(common.mask_file(args.map_dir, args.bands, args.vis_min, args.nside))
    depth = hsp.HealSparseMap.read(common.depth_file(args.map_dir, args.bands, args.vis_min, args.nside))
    print(f'geometry mask: {geometry.get_valid_area(degrees=True):.1f} deg^2')

    ## reading and cleaning DESI
    desi_path = find_desi_file(args)
    print(f'reading {desi_path}')
    specz = read_clean_desi(desi_path, geometry)
    if len(specz) == 0:
        sys.exit('no spec-z inside the footprint')
    print(specz['PROGRAM'].value_counts().to_string(), '\n')

    ## how much the redmapper z_err cut removes
    print(f'{np.mean(specz["ZERR"] < args.zerr_max) * 100:.1f}% of in-footprint spec-z pass z_err < {args.zerr_max}')

    ## candidate pixels and the training region
    cells, bin_counts, survey_median = training_cells(specz, geometry, depth, args)
    print(f'survey median limmag: {survey_median:.2f}')
    print(cells.sort_values('n_specz', ascending=False).head(20).round(2).to_string(), '\n')

    if args.training_hpix is not None:
        region = [int(h) for h in args.training_hpix.split(',')]
        print('training pixels from --training-hpix')
    else:
        region = choose_region(cells, bin_counts, survey_median, args)
    region = [h for h in region if h in cells.index]
    print(cells.loc[region, ['center_ra', 'center_dec', 'area', 'filled', 'median_limmag', 'n_specz']].round(2).to_string())
    print(f'total: {cells.loc[region, "area"].sum():.1f} deg^2, {cells.loc[region, "n_specz"].sum()} spec-z\n')

    ## per-bin check for the region, split by program
    z_edges = z_bins(args)
    n_bins = len(z_edges) - 1
    usable = specz[specz['ZERR'] < args.zerr_max]
    usable_hpix = hp.ang2pix(args.nside_train, usable['TARGET_RA'].values, usable['TARGET_DEC'].values,
                             lonlat=True, nest=False)
    usable_bin = np.digitize(usable['Z'].values, z_edges) - 1
    region_rows = np.isin(usable_hpix, region) & (usable_bin >= 0) & (usable_bin < n_bins)
    per_bin = pd.crosstab(usable_bin[region_rows], usable['PROGRAM'].values[region_rows])
    per_bin = per_bin.reindex(range(n_bins), fill_value=0)
    per_bin.insert(0, 'z_lo', z_edges[:-1].round(2))
    per_bin['total'] = bin_counts.loc[region].sum().values
    per_bin['enough'] = per_bin['total'] >= args.per_bin_min
    per_bin.index.name = 'bin'
    per_bin.columns.name = None
    print(per_bin.to_string())
    print(f'weakest bin: {per_bin["total"].min()} (z_lo = {per_bin.loc[per_bin["total"].idxmin(), "z_lo"]})\n')

    ## writing the region, the cleaned table, and the redmapper spec file
    os.makedirs(args.specz_dir, exist_ok=True)
    path = hpix_file(args.specz_dir, args.bands, args.vis_min, args.nside_train)
    np.savetxt(path, region, fmt='%d')
    print(f'wrote training pixels (RING, nside {args.nside_train}) to {path}')

    specz.to_parquet(clean_file(args.specz_dir, args.bands, args.vis_min), index=False)

    spec_out = np.zeros(len(specz), dtype=[('ra', 'f8'), ('dec', 'f8'), ('z', 'f4'), ('z_err', 'f4')])
    spec_out['ra'] = specz['TARGET_RA'].values
    spec_out['dec'] = specz['TARGET_DEC'].values
    spec_out['z'] = specz['Z'].values
    spec_out['z_err'] = specz['ZERR'].values
    path = spec_file(args.specz_dir, args.bands, args.vis_min)
    fitsio.write(path, spec_out, clobber=True)
    print(f'wrote {len(spec_out)} spec-z to {path}')

    ## reading it back the way redmapper does
    check = fitsio.read(path, ext=1, lower=True, trim_strings=True)
    print(f'{len(check)} rows, z {check["z"].min():.3f} to {check["z"].max():.3f}, '
          f'all finite: {all(np.isfinite(check[c]).all() for c in check.dtype.names)}')
    print('Done!')


'''
red
'''
def load_template():
    '''
    redMaPPer's own starting red-sequence colors, bc03_colors_lsst.fit.

    Returns
    -------
    template_color: `function`
       template_color(blue, red, z) -> blue - red color of the red sequence at z.
    '''
    import redmapper
    path = os.path.join(os.path.dirname(redmapper.__file__), 'data', 'initcolors', 'bc03_colors_lsst.fit')
    template = fitsio.read(path, ext=1, lower=True)
    template_bands = fitsio.read_header(path, ext=1)['BANDS'].split(',')
    print(f'{path}: bands {template_bands}, z {template["z"].min():.2f} to {template["z"].max():.2f}')

    def template_color(blue, red, z):
        ## COLOR column k is band k minus band k+1
        k = template_bands.index(blue)
        assert template_bands[k + 1] == red
        return np.interp(z, template['z'], template['color'][:, k])
    return template_color


def match_tracts(specz, flux_files, args, bands, ref_i):
    '''
    Match spec-z to catalog galaxies, one tract file at a time.

    Each tract goes through the same cuts as ``build``, then each spec-z gets its
    nearest galaxy within ``--match-radius``. Also finds the tracts with objects in
    the training pixels.

    Parameters
    ----------
    specz: `np.ndarray`
       Spec file rows (ra, dec, z, z_err).
    flux_files: `list` [`str`]
       Tract flux files from extract.
    args: `argparse.Namespace`
       Parsed arguments, carries the build's model choices and mag cut.
    bands: `list` [`str`]
       Bands in the catalog.
    ref_i: `int`
       Index of the reference band in ``bands``.

    Returns
    -------
    matched: `pd.DataFrame`
       One row per matched spec-z, with the galaxy's refmag and colors.
    training_tracts: `list` [`int`]
       Tracts with objects in the training pixels.
    '''
    lo_b, lo_r = [bands.index(b) for b in args.low_color.split('-')]
    hi_b, hi_r = [bands.index(b) for b in args.high_color.split('-')]
    radius = args.match_radius / 3600.0

    spec_coords = SkyCoord(specz['ra'] * u.deg, specz['dec'] * u.deg)
    matched = []
    training_tracts = []

    for n, path in enumerate(flux_files):
        flux = fitsio.read(path)
        tract = int(os.path.basename(path).split('_')[1].replace('tract', ''))

        ## noting tracts that reach into the training pixels
        flux_hpix = hp.ang2pix(args.nside_train, flux['ra'], flux['dec'], lonlat=True, nest=False)
        if np.isin(flux_hpix, args.region).any():
            training_tracts.append(tract)

        gals, _ = rm_cats.fluxes_to_galaxies(flux, args, ref_i, None, bands)
        if len(gals) == 0:
            continue

        ## spec-z near this tract (handles tracts that cross RA = 0)
        ra_g = gals['ra'].copy()
        ra_s = specz['ra'].copy()
        if ra_g.max() - ra_g.min() > 180:
            ra_g = np.where(ra_g > 180, ra_g - 360, ra_g)
            ra_s = np.where(ra_s > 180, ra_s - 360, ra_s)
        near = ((ra_s > ra_g.min() - 0.01) & (ra_s < ra_g.max() + 0.01) &
                (specz['dec'] > gals['dec'].min() - 0.01) & (specz['dec'] < gals['dec'].max() + 0.01))
        if near.sum() == 0:
            continue

        gal_coords = SkyCoord(gals['ra'] * u.deg, gals['dec'] * u.deg)
        idx, sep, _ = spec_coords[near].match_to_catalog_sky(gal_coords)
        good = sep.deg < radius
        if good.sum() == 0:
            continue

        s = specz[near][good]
        g = gals[idx[good]]
        matched.append(pd.DataFrame({
            'ra': s['ra'], 'dec': s['dec'], 'z': s['z'], 'z_err': s['z_err'],
            'refmag': g['refmag'],
            'low_color': g['mag'][:, lo_b] - g['mag'][:, lo_r],
            'high_color': g['mag'][:, hi_b] - g['mag'][:, hi_r],
            'low_color_err': np.hypot(g['mag_err'][:, lo_b], g['mag_err'][:, lo_r]),
            'high_color_err': np.hypot(g['mag_err'][:, hi_b], g['mag_err'][:, hi_r]),
            'sep_arcsec': sep.arcsec[good],
        }))

        if (n + 1) % 200 == 0:
            print(f'  [{n+1}/{len(flux_files)}] {sum(len(m) for m in matched)} matched so far', flush=True)

    if len(matched) == 0:
        return pd.DataFrame(), training_tracts
    matched = pd.concat(matched, ignore_index=True)

    ## a spec-z can land in two tract files near tract edges; keep one
    before = len(matched)
    matched = matched.sort_values('sep_arcsec').drop_duplicates(subset=['ra', 'dec', 'z']).reset_index(drop=True)
    print(f'{len(matched)} spec-z matched to a catalog galaxy ({before - len(matched)} duplicates dropped), '
          f'median separation {matched["sep_arcsec"].median():.2f}"')
    return matched, training_tracts


def gauss(x, mu, sigma):
    return np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (np.sqrt(2 * np.pi) * sigma)


def fit_red(matched, args):
    '''
    Red vs. blue in each redshift bin, filling ``matched['p_red']``.

    Same fit as redmapper's calibration/selectspecred.py: trim to the 5th-99th
    percentile of the offset from the template, then
    ``EcgmmFitter.fit([0.2], [-0.5, 0.0], [0.2, 0.05], offset=2.0)``. The second
    Gaussian (starting on the template) is the red one. Each galaxy then gets a red
    probability from the two Gaussians, each widened by its color error.

    Parameters
    ----------
    matched: `pd.DataFrame`
       From `match_tracts`, with ``delta`` and ``color_err`` filled in.
    args: `argparse.Namespace`
       Parsed arguments.

    Returns
    -------
    fits: `pd.DataFrame`
       Fit results per bin.
    '''
    try:
        from redmapper.fitters import EcgmmFitter
    except ImportError:
        from redmapper.calibration.selectspecred import EcgmmFitter

    matched['p_red'] = np.nan
    fits = []
    z_edges = z_bins(args)
    for lo, hi in zip(z_edges[:-1], z_edges[1:]):
        in_bin = (matched['z'] >= lo) & (matched['z'] < hi)
        delta = matched.loc[in_bin, 'delta'].values
        err = matched.loc[in_bin, 'color_err'].values
        if len(delta) < 20:
            print(f'z {lo:.2f}-{hi:.2f}: only {len(delta)} galaxies, skipping')
            continue

        p5, p99 = np.percentile(delta, [5, 99])
        trim = (delta > p5) & (delta < p99)
        try:
            wt, mu, sigma = EcgmmFitter(delta[trim], err[trim]).fit([0.2], [-0.5, 0.0], [0.2, 0.05], offset=2.0)
        except Exception as e:
            print(f'z {lo:.2f}-{hi:.2f}: fit failed ({e})')
            continue

        wt = np.atleast_1d(wt)
        w_blue, w_red = (wt[0], wt[1]) if len(wt) > 1 else (wt[0], 1 - wt[0])

        like_blue = w_blue * gauss(delta, mu[0], np.hypot(sigma[0], err))
        like_red = w_red * gauss(delta, mu[1], np.hypot(sigma[1], err))
        matched.loc[in_bin, 'p_red'] = like_red / (like_red + like_blue)

        fits.append({'z_lo': round(lo, 2), 'color': matched.loc[in_bin, 'color_name'].iloc[0],
                     'n': len(delta), 'w_red': w_red, 'mu_blue': mu[0], 'mu_red': mu[1],
                     'sigma_blue': sigma[0], 'sigma_red': sigma[1]})
    return pd.DataFrame(fits)


def run_red(args, bands, ref_i):
    '''
    Red step: keep the red-sequence spec-z, write the red spec file and the training tracts.

    The red/blue split is on the color offset from redmapper's bc03 template:
    ``--low-color`` below ``--switch-z`` and ``--high-color`` above (redmapper
    how-to; with no g band the low-z bins lean on r-i). Runs on whatever tracts
    extract has finished and labels the output ``_partial`` until all are done.

    Parameters
    ----------
    args: `argparse.Namespace`
       Parsed arguments.
    bands: `list` [`str`]
       Bands in the catalog.
    ref_i: `int`
       Index of the reference band in ``bands``.
    '''
    ## spec-z, training pixels, and template
    specz = fitsio.read(spec_file(args.specz_dir, args.bands, args.vis_min), lower=True)
    specz = specz.astype(specz.dtype.newbyteorder('='))   ## FITS is big-endian, pandas won't sort it
    specz = specz[specz['z_err'] < args.zerr_max]   ## redmapper drops the rest anyway
    print(f'{len(specz)} spec-z')

    if args.training_hpix is not None:
        args.region = [int(h) for h in args.training_hpix.split(',')]
    else:
        args.region = np.loadtxt(hpix_file(args.specz_dir, args.bands, args.vis_min, args.nside_train),
                                 dtype=int, ndmin=1).tolist()
    print(f'training pixels (RING, nside {args.nside_train}): {args.region}')
    template_color = load_template()

    ## which tracts are done
    flux_files = sorted(glob.glob(common.flux_file(args.cat_dir, '*', args.bands, args.vis_min)))
    n_empty = len(glob.glob(common.flux_file(args.cat_dir, '*', args.bands, args.vis_min) + '.empty'))
    n_done = len(flux_files) + n_empty
    partial = n_done < args.n_tracts
    print(f'{n_done} of {args.n_tracts} tracts extracted' + (' (partial run)' if partial else ''))
    print(f'catalog cuts: refmag {args.refmag_model} {args.ref_band} < {args.mag_cut}, colors {args.color_model}')

    ## matching
    matched, training_tracts = match_tracts(specz, flux_files, args, bands, ref_i)

    tag = '_partial' if partial else ''
    path = os.path.join(args.specz_dir, f'training_tracts_{args.bands}_vis{args.vis_min}{tag}.txt')
    np.savetxt(path, training_tracts, fmt='%d')
    print(f'wrote {len(training_tracts)} training tracts to {path}')
    if len(matched) == 0:
        sys.exit('no spec-z matched a catalog galaxy, nothing to fit')

    ## color offset from the template
    low_z = (matched['z'] < args.switch_z).values
    lo_b, lo_r = args.low_color.split('-')
    hi_b, hi_r = args.high_color.split('-')
    matched['color_name'] = np.where(low_z, args.low_color, args.high_color)
    matched['color'] = np.where(low_z, matched['low_color'], matched['high_color'])
    matched['color_err'] = np.where(low_z, matched['low_color_err'], matched['high_color_err'])
    matched['template'] = np.where(low_z, template_color(lo_b, lo_r, matched['z']),
                                   template_color(hi_b, hi_r, matched['z']))
    matched['delta'] = matched['color'] - matched['template']

    ## red vs. blue per bin. the red gaussian should be the narrow one near offset 0,
    ## if a bin's red and blue come out swapped or on top of each other don't trust that bin
    fits = fit_red(matched, args)
    print(fits.round(3).to_string(index=False), '\n')

    ## red spec-z per bin
    matched['red'] = matched['p_red'] > args.p_red   ## more likely red than blue
    in_training = np.isin(hp.ang2pix(args.nside_train, matched['ra'].values, matched['dec'].values,
                                     lonlat=True, nest=False), args.region)
    z_edges = z_bins(args)
    per_bin = pd.DataFrame({
        'z_lo': z_edges[:-1].round(2),
        'color': np.where(z_edges[:-1] < args.switch_z, args.low_color, args.high_color),
        'matched': np.histogram(matched['z'], bins=z_edges)[0],
        'red': np.histogram(matched['z'][matched['red']], bins=z_edges)[0],
        'red, training': np.histogram(matched['z'][matched['red'] & in_training], bins=z_edges)[0],
    })
    print(per_bin.to_string(index=False))
    print(f'\n{matched["red"].sum()} red spec-z total, {(matched["red"] & in_training).sum()} in the training region'
          + ('   (partial: extract not finished)' if partial else ''))

    ## writing the red spec file, same four columns redmapper reads
    red = matched[matched['red']]
    spec_out = np.zeros(len(red), dtype=[('ra', 'f8'), ('dec', 'f8'), ('z', 'f4'), ('z_err', 'f4')])
    spec_out['ra'] = red['ra'].values
    spec_out['dec'] = red['dec'].values
    spec_out['z'] = red['z'].values
    spec_out['z_err'] = red['z_err'].values
    path = spec_file(args.specz_dir, args.bands, args.vis_min, red=True, partial=partial)
    fitsio.write(path, spec_out, clobber=True)

    ## keeping the whole matched table too, for plots and for checking the selection later
    matched.to_parquet(os.path.join(args.specz_dir, f'dp2_specz_matched_{args.bands}_vis{args.vis_min}{tag}.parquet'),
                       index=False)
    print(f'wrote {len(spec_out)} red spec-z to {path}')
    print('Done!')


def parse_args():
    '''
    Command line arguments.

    Returns
    -------
    args: `argparse.Namespace`
       Parsed arguments.
    '''
    p = argparse.ArgumentParser(description='DESI DR1 spec-z training sample for redmapper on DP2.')
    p.add_argument('step', choices=['download', 'select', 'red'])

    p.add_argument('--bands', default='riz', help='must match the maps and catalog (default riz)')
    p.add_argument('--vis-min', type=int, default=3, help='must match the maps (default 3)')
    p.add_argument('--nside', type=int, default=4096, help='must match the maps (default 4096)')
    p.add_argument('--map-dir', default=os.path.expanduser('~/dp2_redmapper/maps'))
    p.add_argument('--cat-dir', default=os.path.expanduser('~/dp2_redmapper/catalogs'))
    p.add_argument('--specz-dir', default=os.path.expanduser('~/dp2_redmapper/specz'))

    ## select
    p.add_argument('--desi-file', default=None,
                   help=f'DESI parquet or FITS (default: download in --specz-dir, else {DESI_NERSC})')
    p.add_argument('--nside-train', type=int, default=8,
                   help='RING nside of the training pixels (default 8, ~54 deg^2 each)')
    p.add_argument('--target-area', type=float, default=200.0, help='training area in deg^2 (default 200)')
    p.add_argument('--depth-tol', type=float, default=0.15,
                   help='candidate pixels within this many mag of the survey median depth (default 0.15)')
    p.add_argument('--training-hpix', default=None,
                   help='comma separated RING pixels, to set the training region by hand')
    p.add_argument('--zerr-max', type=float, default=0.001, help='calib_spec_max_zerr in redmapper (default 0.001)')
    p.add_argument('--zmin', type=float, default=0.10, help='per-bin counts from here (default 0.10)')
    p.add_argument('--zmax', type=float, default=0.90, help='per-bin counts up to here (default 0.90)')
    p.add_argument('--dz', type=float, default=0.05, help='per-bin width (default 0.05, cosmoDC2 setup)')
    p.add_argument('--per-bin-min', type=int, default=40, help='spec-z wanted per bin (default 40, cosmoDC2)')

    ## red
    p.add_argument('--ref-band', default='z', help='reference band (default z)')
    p.add_argument('--refmag-model', choices=rm_cats.MODELS, default='cmodel', help='same as build (default cmodel)')
    p.add_argument('--color-model', choices=rm_cats.MODELS, default='sersic', help='same as build (default sersic)')
    p.add_argument('--mag-cut', type=float, default=24.0, help='same as build (default 24)')
    p.add_argument('--match-radius', type=float, default=3.0,
                   help='arcsec, from calibration/selectspecred.py (default 3)')
    p.add_argument('--switch-z', type=float, default=0.7,
                   help='use --low-color below this and --high-color above (default 0.7)')
    p.add_argument('--low-color', default='r-i', help='color below --switch-z (default r-i)')
    p.add_argument('--high-color', default='i-z', help='color above --switch-z (default i-z)')
    p.add_argument('--p-red', type=float, default=0.5, help='keep p_red above this (default 0.5)')
    p.add_argument('--n-tracts', type=int, default=N_TRACTS_TOTAL,
                   help=f'tracts in a full extract, below this the output is _partial (default {N_TRACTS_TOTAL})')

    return p.parse_args()


def main():
    args = parse_args()
    bands = list(args.bands)
    if args.ref_band not in bands:
        sys.exit(f'reference band {args.ref_band} is not in {args.bands}')
    ref_i = bands.index(args.ref_band)

    if args.step == 'download':
        run_download(args)
    elif args.step == 'select':
        run_select(args)
    else:
        run_red(args, bands, ref_i)


if __name__ == '__main__':
    main()
