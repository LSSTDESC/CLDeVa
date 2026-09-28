'''
create_redmapper_dp2_catalogs.py
'''

import sys; import os
import glob
import argparse
import multiprocessing as mp

import numpy as np
import healsparse as hsp
import fitsio

from astropy.coordinates import SkyCoord
import astropy.units as u

import redmapper_dp2_common as common
from dustmaps.sfd import SFDQuery

import redmapper

'''
Initial definitions
'''
zeropoint = 22.5  # used this in cosmoDC2

## no GAaP rn
MODELS = ['cmodel', 'sersic']

## flux columns per model
FLUX_COLS = {'cmodel': ('{b}_cModelFlux', '{b}_cModelFluxErr'),
             'sersic': ('{b}_sersicFlux', '{b}_sersicFluxErr')}

EXT_COL = 'griz_model_extendedness'
EXT_CUT = 0.5 #

## GalaxyCatalogMaker._check_galaxies refuses mag >= 90, and mag_err >= 90
MAG_LIMIT = 90.0

START_METHOD = 'fork'


'''
columns
'''
def flux_columns(model, bands):
    '''
    Flux and flux error column names for one model.

    Parameters
    ----------
    model: `str`
       Photometry model, 'cmodel' or 'sersic'.
    bands: `list` [`str`]
       Bands to build names for.

    Returns
    -------
    fluxes: `list` [`str`]
       Flux column per band.
    errs: `list` [`str`]
       Flux error column per band.
    '''
    fpat, epat = FLUX_COLS[model]
    return [fpat.format(b=b) for b in bands], [epat.format(b=b) for b in bands]


def get_columns(bands):
    '''
    Everything extract needs, for all models at once.

    Parameters
    ----------
    bands: `list` [`str`]
       Bands to read.

    Returns
    -------
    cols: `list` [`str`]
       Object table columns to ask the butler for.
    '''
    cols = ['objectId', 'coord_ra', 'coord_dec', EXT_COL,
            'sersic_no_data_flag', 'sersic_unknown_flag']
    for model in MODELS:
        fluxes, errs = flux_columns(model, bands)
        cols += fluxes + errs
    for b in bands:
        cols += [f'{b}_pixelFlags_saturatedCenter', f'{b}_pixelFlags_nodata']
    return cols


def flux_dtype(n_bands):
    '''
    What each tract file holds: fluxes in nJy, one (n_bands) array per model.

    Parameters
    ----------
    n_bands: `int`
       Number of bands.

    Returns
    -------
    dtype: `list` [`tuple`]
       Numpy dtype spec for the per-tract flux array.
    '''
    dtype = [('id', 'i8'), ('ra', 'f8'), ('dec', 'f8'), ('ebv', 'f4'),
             ('sersic_ok', 'i2')]
    for model in MODELS:
        dtype += [(f'{model}_flux', 'f4', n_bands), (f'{model}_flux_err', 'f4', n_bands)]
    return dtype


'''
extract
'''
def select_objects(d, skymap, tract, mask, bands):
    '''
    Cuts that don't depend on the photometry model.

    Footprint, tract ownership, extendedness and pixel flags.

    Parameters
    ----------
    d: `astropy.table.Table`
       One tract of the object table.
    skymap: `lsst.skymap.BaseSkyMap`
       Skymap, for tract ownership.
    tract: `int`
       Tract id.
    mask: `healsparse.HealSparseMap`
       Cell-resolution boolean geometry.
    bands: `list` [`str`]
       Bands whose pixel flags must be clean.

    Returns
    -------
    selected: `np.ndarray`
       Boolean, True for objects that survive every cut.
    funnel: `dict` [`str`, `int`]
       How many objects survived each stage, for the log.
    '''
    ra = common.get_col(d, 'coord_ra', np.nan).astype(float)
    dec = common.get_col(d, 'coord_dec', np.nan).astype(float)
    ok = np.isfinite(ra) & np.isfinite(dec)

    in_footprint = np.zeros(len(ra), dtype=bool)
    in_footprint[ok] = mask.get_values_pos(ra[ok], dec[ok], lonlat=True) > 0

    in_tract = np.zeros(len(ra), dtype=bool)
    in_tract[in_footprint] = common.in_tract_inner(skymap, tract, ra[in_footprint], dec[in_footprint])

    model_extension = common.get_col(d, EXT_COL, np.nan).astype(float) > EXT_CUT

    clean_pixels = np.ones(len(ra), dtype=bool)
    for b in bands:
        clean_pixels &= ~common.get_col(d, f'{b}_pixelFlags_saturatedCenter', True).astype(bool)
        clean_pixels &= ~common.get_col(d, f'{b}_pixelFlags_nodata', True).astype(bool)

    selected = in_footprint & in_tract & model_extension & clean_pixels
    return selected, {'objects': len(ra), 'footprint': in_footprint.sum(), 'tract': in_tract.sum(),
                      'extended': (in_tract & model_extension).sum(), 'clean': selected.sum()}


def extract_tract(butler, skymap, tract, mask, bands, ref_i, extract_mag_cut, sfd):
    '''
    One tract of the object table -> one flux array with every model.

    Parameters
    ----------
    butler: `lsst.daf.butler.Butler`
       Butler to read from.
    skymap: `lsst.skymap.BaseSkyMap`
       Skymap.
    tract: `int`
       Tract id.
    mask: `healsparse.HealSparseMap`
       Cell-resolution boolean geometry.
    bands: `list` [`str`]
       Bands to extract.
    ref_i: `int`
       Index of the reference band in ``bands``.
    extract_mag_cut: `float`
       Keep objects brighter than this in the reference band in any model.
    sfd: `dustmaps.sfd.SFDQuery`
       Dust map, queried for E(B-V) per object.

    Returns
    -------
    out: `np.ndarray`
       Structured flux array for this tract.
    funnel: `dict` [`str`, `int`]
       How many objects survived each stage, for the log.
    '''
    d = common.read_tract(butler, tract, get_columns(bands))
    selected, funnel = select_objects(d, skymap, tract, mask, bands)

    ## keeping anything bright enough in any model
    bright_any = np.zeros(len(d), dtype=bool)
    for model in MODELS:
        fluxes, _ = flux_columns(model, bands)
        ref_mag = common.flux_to_mag(common.get_col(d, fluxes[ref_i], np.nan))
        bright_any |= np.isfinite(ref_mag) & (ref_mag < extract_mag_cut)
    keep = selected & bright_any
    funnel['bright'] = keep.sum()

    out = np.zeros(keep.sum(), dtype=flux_dtype(len(bands)))
    out['id'] = common.get_col(d, 'objectId', -1)[keep]
    out['ra'] = common.get_col(d, 'coord_ra', np.nan)[keep]
    out['dec'] = common.get_col(d, 'coord_dec', np.nan)[keep]

    sersic_ok = (~common.get_col(d, 'sersic_no_data_flag', True).astype(bool)
                 & ~common.get_col(d, 'sersic_unknown_flag', True).astype(bool))
    out['sersic_ok'] = sersic_ok[keep]

    for model in MODELS:
        fluxes, errs = flux_columns(model, bands)
        out[f'{model}_flux'] = np.column_stack([common.get_col(d, f, np.nan)[keep] for f in fluxes])
        out[f'{model}_flux_err'] = np.column_stack([common.get_col(d, e, np.nan)[keep] for e in errs])

    if len(out) > 0:
        out['ebv'] = sfd(SkyCoord(out['ra'] * u.deg, out['dec'] * u.deg))

    return out, funnel


## setting up one butler, mask and dust map per worker
_worker = {}


def init_worker(args, bands, ref_i):
    '''
    Set up one extract worker.

    Parameters
    ----------
    args: `argparse.Namespace`
       Parsed arguments.
    bands: `list` [`str`]
       Bands to extract.
    ref_i: `int`
       Index of the reference band in ``bands``.
    '''
    _worker['butler'], _worker['skymap'] = common.open_butler()
    _worker['mask'] = hsp.HealSparseMap.read(mask_path(args, fine=True))
    _worker['sfd'] = SFDQuery()
    _worker['args'] = args
    _worker['bands'] = bands
    _worker['ref_i'] = ref_i


def extract_worker(tract):
    '''
    Pool entry point: extract one tract with this worker's butler.

    Parameters
    ----------
    tract: `int`
       Tract id.

    Returns
    -------
    message: `str`
       Progress line for the log.
    '''
    return extract_and_save(_worker['butler'], _worker['skymap'], _worker['mask'],
                            _worker['sfd'], tract, _worker['args'], _worker['bands'],
                            _worker['ref_i'])


def extract_and_save(butler, skymap, mask, sfd, tract, args, bands, ref_i):
    '''
    Extract one tract and write its flux file.

    Writes to a temp

    Parameters
    ----------
    butler: `lsst.daf.butler.Butler`
       Butler to read from.
    skymap: `lsst.skymap.BaseSkyMap`
       Skymap.
    mask: `healsparse.HealSparseMap`
       Cell-resolution boolean geometry.
    sfd: `dustmaps.sfd.SFDQuery`
       Dust map.
    tract: `int`
       Tract id.
    args: `argparse.Namespace`
       Parsed arguments.
    bands: `list` [`str`]
       Bands to extract.
    ref_i: `int`
       Index of the reference band in ``bands``.

    Returns
    -------
    message: `str`
       Progress line for the log.
    '''
    path = common.flux_file(args.out_dir, tract, args.bands, args.vis_min)

    try:
        out, funnel = extract_tract(butler, skymap, tract, mask, bands, ref_i,
                                    args.extract_mag_cut, sfd)
    except Exception as err:
        ## gcs throws the odd error
        return f'tract {tract} FAILED: {err}'

    ## clearing out old files
    for old_file in [path, path + '.empty']:
        if os.path.exists(old_file):
            os.remove(old_file)

    ## marking the empty tracts so they get skipped next time
    if len(out) == 0:
        open(path + '.empty', 'w').close()
    else:
        tmp = path + '.tmp'
        fitsio.write(tmp, out)
        os.replace(tmp, path)

    return f'tract {tract}: ' + ' -> '.join(f'{k} {v}' for k, v in funnel.items())


def tract_done(out_dir, tract, bands, vis_min):
    '''
    Has this tract already been extracted?

    Parameters
    ----------
    out_dir: `str`
       Catalog directory.
    tract: `int`
       Tract id.
    bands: `str`
       Band string, e.g. 'riz'.
    vis_min: `int`
       Visits required per band in the geometry.

    Returns
    -------
    done: `bool`
       True if there's a flux file or an empty marker.
    '''
    path = common.flux_file(out_dir, tract, bands, vis_min)
    return os.path.exists(path) or os.path.exists(path + '.empty')


def run_extract(args, bands, ref_i):
    '''
    Extract step: one flux file per tract.

    Tracts that already have a file are skipped unless ``--clobber``. Like the
    maps, this is nearly all butler latency, so workers barely touch the CPU.

    Parameters
    ----------
    args: `argparse.Namespace`
       Parsed arguments.
    bands: `list` [`str`]
       Bands to extract.
    ref_i: `int`
       Index of the reference band in ``bands``.
    '''
    ## opening the butler and getting the tract list
    butler, skymap = common.open_butler()
    tracts = common.get_tract_list(args.tracts, args.tract_list, butler)
    os.makedirs(os.path.join(args.out_dir, 'tract_fluxes'), exist_ok=True)

    ## skipping tracts that are already done
    todo = [int(t) for t in tracts
            if args.clobber or not tract_done(args.out_dir, t, args.bands, args.vis_min)]
    print(f'extract: ref {args.ref_band} < {args.extract_mag_cut} in any model')
    print(f'{len(todo)} tracts to do, {len(tracts) - len(todo)} already done, '
          f'{args.n_procs} process(es)')

    ## extracting one tract at a time, or with a pool of workers
    if args.n_procs == 1:
        mask = hsp.HealSparseMap.read(mask_path(args, fine=True))
        sfd = SFDQuery()
        for n, tract in enumerate(todo):
            print(f'[{n+1}/{len(todo)}] '
                  + extract_and_save(butler, skymap, mask, sfd, tract, args, bands, ref_i),
                  flush=True)
    else:
        with mp.get_context(START_METHOD).Pool(args.n_procs, initializer=init_worker,
                                               initargs=(args, bands, ref_i)) as pool:
            for n, message in enumerate(pool.imap_unordered(extract_worker, todo)):
                print(f'[{n+1}/{len(todo)}] {message}', flush=True)

    print('Done!')


'''
build
'''
def redmapper_valid(ra, dec, mags, mag_errs, refmag, refmag_err):
    '''
    Drop galaxies GalaxyCatalogMaker._check_galaxies would reject.

    It refuses the whole catalog if any object has a non-finite mag or error, a
    value >= 90, a value <= 0, or out-of-range coordinates.

    Parameters
    ----------
    ra: `np.ndarray`
       Right ascension, degrees.
    dec: `np.ndarray`
       Declination, degrees.
    mags: `np.ndarray`
       Color magnitudes, shape (n_gal, n_band).
    mag_errs: `np.ndarray`
       Color magnitude errors, shape (n_gal, n_band).
    refmag: `np.ndarray`
       Total magnitude in the reference band.
    refmag_err: `np.ndarray`
       Error on refmag.

    Returns
    -------
    keep: `np.ndarray`
       Boolean, True for galaxies redmapper will accept.
    dropped: `dict` [`str`, `int`]
       How many failed each rule, for the log.
    '''
    ## finding non-finite mags and errors
    finite = (np.all(np.isfinite(mags), axis=1) & np.all(np.isfinite(mag_errs), axis=1)
              & np.isfinite(refmag) & np.isfinite(refmag_err))

    ## checking the mag, err and coordinate ranges
    with np.errstate(invalid='ignore'):
        below_90 = np.all(mags < MAG_LIMIT, axis=1) & np.all(mag_errs < MAG_LIMIT, axis=1)
        positive = np.all(mags > 0, axis=1) & np.all(mag_errs > 0, axis=1)
        on_sky = (ra >= 0) & (ra <= 360) & (dec >= -90) & (dec <= 90)

    keep = finite & below_90 & positive & on_sky
    dropped = {'non-finite': int((~finite).sum()),
               'mag/err >= 90': int((finite & ~below_90).sum()),
               'mag/err <= 0': int((finite & below_90 & ~positive).sum()),
               'off-sky ra/dec': int((finite & below_90 & positive & ~on_sky).sum())}
    return keep, dropped


def fluxes_to_galaxies(flux, args, ref_i, dust_coeffs, bands):
    '''
    One tract's flux array -> redmapper galaxy array, plus how many each cut dropped.

    Parameters
    ----------
    flux: `np.ndarray`
       Structured flux array from the extract step.
    args: `argparse.Namespace`
       Parsed arguments, for the model choices and the magnitude cut.
    ref_i: `int`
       Index of the reference band in ``bands``.
    dust_coeffs: `dict` [`str`, `float`] or `None`
       A_band/E(B-V) per band, or None to leave magnitudes reddened.
    bands: `list` [`str`]
       Bands in the catalog.

    Returns
    -------
    galaxies: `np.ndarray`
       Structured array in redmapper's galaxy format.
    dropped: `dict` [`str`, `int`]
       How many failed each cut, for the log.
    '''
    ## converting fluxes to mags
    with np.errstate(invalid='ignore', divide='ignore'):
        mags = common.flux_to_mag(flux[f'{args.color_model}_flux'])
        mag_errs = common.flux_to_magerr(flux[f'{args.color_model}_flux'], flux[f'{args.color_model}_flux_err'])
        refmag = common.flux_to_mag(flux[f'{args.refmag_model}_flux'][:, ref_i])
        refmag_err = common.flux_to_magerr(flux[f'{args.refmag_model}_flux'][:, ref_i],
                                           flux[f'{args.refmag_model}_flux_err'][:, ref_i])

    if dust_coeffs is not None:
        extinction = np.column_stack([dust_coeffs[b] * flux['ebv'] for b in bands])
        mags = mags - extinction
        refmag = refmag - extinction[:, ref_i]

    ## cutting on refmag and sersic convergence
    dropped = {}
    with np.errstate(invalid='ignore'):
        bright = np.isfinite(refmag) & (refmag < args.mag_cut)
    dropped[f'refmag >= {args.mag_cut}'] = int((~bright).sum())

    if 'sersic' in (args.refmag_model, args.color_model):
        converged = flux['sersic_ok'].astype(bool)
        dropped['sersic not converged'] = int((bright & ~converged).sum())
        bright &= converged
        
    ## dropping anything redmapper would reject
    valid, rule_drops = redmapper_valid(flux['ra'][bright], flux['dec'][bright], mags[bright],
                                        mag_errs[bright], refmag[bright], refmag_err[bright])
    dropped.update(rule_drops)
    keep = bright.copy()
    keep[bright] = valid

    ## filling in the galaxy array
    gal_dtype = [('id', 'i8'), ('ra', 'f8'), ('dec', 'f8'),
                 ('refmag', 'f4'), ('refmag_err', 'f4'),
                 ('mag', 'f4', len(bands)), ('mag_err', 'f4', len(bands)),
                 ('ebv', 'f4')]
    galaxies = np.zeros(keep.sum(), dtype=gal_dtype)
    galaxies['id'] = flux['id'][keep]
    galaxies['ra'] = flux['ra'][keep]
    galaxies['dec'] = flux['dec'][keep]
    galaxies['refmag'] = refmag[keep]
    galaxies['refmag_err'] = refmag_err[keep]
    galaxies['mag'] = mags[keep]
    galaxies['mag_err'] = mag_errs[keep]
    galaxies['ebv'] = flux['ebv'][keep]

    return galaxies, dropped


def make_info_dict(bands, ref_i, mag_cut, area):
    '''
    Header info for GalaxyCatalogMaker.

    Parameters
    ----------
    bands: `list` [`str`]
       Bands in the catalog.
    ref_i: `int`
       Index of the reference band in ``bands``.
    mag_cut: `float`
       Faint limit in the reference band.
    area: `float`
       Footprint area in square degrees.

    Returns
    -------
    info_dict: `dict`
       Keywords for the catalog master table.
    '''
    info_dict = {'LIM_REF': mag_cut,
                 'REF_IND': ref_i,
                 'AREA': area,
                 'NMAG': len(bands),
                 'MODE': 'LSST',
                 'ZP': zeropoint}
    for i, b in enumerate(bands):
        info_dict[f'{b.upper()}_IND'] = i
    return info_dict


def mask_area(mask):
    '''
    Footprint area from the geometry mask, fracgood weighted.

    Parameters
    ----------
    mask: `healsparse.HealSparseMap`
       Float fracgood mask.

    Returns
    -------
    area: `float`
       Area in square degrees.
    '''
    pixels = mask.valid_pixels
    return float(mask.get_values_pix(pixels).sum()) * common.pix_area(mask.nside_sparse)


def run_build(args, bands, ref_i):

    '''
    Build step: every tract flux file -> one redmapper galaxy catalog.

    Parameters
    ----------
    args: `argparse.Namespace`
       Parsed arguments.
    bands: `list` [`str`]
       Bands in the catalog.
    ref_i: `int`
       Index of the reference band in ``bands``.
    '''
    if args.mag_cut > args.extract_mag_cut:
        sys.exit(f'--mag-cut {args.mag_cut} is fainter than --extract-mag-cut {args.extract_mag_cut}, '
                 f'the tract files don\'t have those galaxies')

    dust_coeffs = common.parse_dust_coeffs(args.dust_coeffs, bands)
    base = common.catalog_base(args.out_dir, args.refmag_model, args.color_model,
                               args.bands, args.vis_min, dered=dust_coeffs is not None)

    ## checking for old catalog files first
    old_files = glob.glob(base + '_*.fit')
    if len(old_files) > 0:
        sys.exit(f'{len(old_files)} old catalog files at {base}_*.fit, move or delete them first')

    files = sorted(glob.glob(common.flux_file(args.out_dir, '*', args.bands, args.vis_min)))
    if len(files) == 0:
        sys.exit('no tract flux files, run extract first')

    ## getting the area and setting up the catalog maker
    area = mask_area(hsp.HealSparseMap.read(mask_path(args)))   ## fracgood weighted
    print(f'build: refmag {args.refmag_model} {args.ref_band} < {args.mag_cut}, colors {args.color_model}, '
          f'dust {"corrected with " + str(dust_coeffs) if dust_coeffs else "NOT corrected"}')
    print(f'{len(files)} tracts, area {area:.1f} deg^2 -> {base}')

    maker = redmapper.GalaxyCatalogMaker(base, make_info_dict(bands, ref_i, args.mag_cut, area))

    ## looping through the tract files and appending galaxies
    n_total = 0
    n_input = 0
    dropped_total = {}
    all_ids = []
    for n, path in enumerate(files):
        flux = fitsio.read(path)
        galaxies, dropped = fluxes_to_galaxies(flux, args, ref_i, dust_coeffs, bands)
        n_input += len(flux)
        for rule, count in dropped.items():
            dropped_total[rule] = dropped_total.get(rule, 0) + count

        all_ids.append(galaxies['id'])
        if len(galaxies) > 0:
            maker.append_galaxies(galaxies)
            n_total += len(galaxies)

        if (n + 1) % 100 == 0:
            print(f'  [{n+1}/{len(files)}] {n_total} galaxies so far')

    ## checking the ids are unique across tracts too
    all_ids = np.concatenate(all_ids)
    n_dup = len(all_ids) - len(np.unique(all_ids))
    if n_dup > 0:
        sys.exit(f'{n_dup} duplicate ids across tracts, not finalizing')

    maker.finalize_catalog()

    print(f'\n{n_input} objects in the tract files')
    for rule, count in dropped_total.items():
        print(f'  dropped {count} for {rule}')
    print(f'{n_total} galaxies -> {base}')
    print('Done!')


def mask_path(args, fine=False):
    '''
    Path to the geometry the catalog cuts on.

    Parameters
    ----------
    args: `argparse.Namespace`
       Parsed arguments.
    fine: `bool`, optional
       True for the cell-resolution boolean map, False for the fracgood map at
       ``args.nside``.

    Returns
    -------
    path: `str`
       File path.
    '''
    if fine:
        path = common.bitpacked_file(args.map_dir, args.bands, args.vis_min, args.nside_sparse)
    else:
        path = common.mask_file(args.map_dir, args.bands, args.vis_min, args.nside)
    if not os.path.exists(path):
        sys.exit(f'no geometry mask at {path}, run create_redmapper_dp2_maps.py first')
    return path


def parse_args():
    '''
    Command line arguments.

    Returns
    -------
    args: `argparse.Namespace`
       Parsed arguments.
    '''
    p = argparse.ArgumentParser(description='Photometric galaxy catalog for redmapper from DP2 (RSP).')
    p.add_argument('step', choices=['extract', 'build'])

    p.add_argument('--bands', default='riz', help='must match the maps (default riz)')
    p.add_argument('--ref-band', default='z', help='reference band (default z)')
    p.add_argument('--vis-min', type=int, default=3, help='must match the maps (default 3)')
    p.add_argument('--nside', type=int, default=4096, help='must match the maps (default 4096)')
    p.add_argument('--nside-sparse', type=int, default=131072, help='must match the maps (default 131072)')
    p.add_argument('--map-dir', default=os.path.expanduser('~/dp2_redmapper/maps'))
    p.add_argument('--out-dir', default=os.path.expanduser('~/dp2_redmapper/catalogs'))

    ## extract
    p.add_argument('--tracts', default=None, help='comma separated, e.g. 10188,8228')
    p.add_argument('--tract-list', default=None, help='text file, one tract per line (default: all tracts)')
    p.add_argument('--extract-mag-cut', type=float, default=24.0,
                   help='keep objects brighter than this in the reference band in any model (default 24)')
    p.add_argument('--clobber', action='store_true', help='redo tracts that already have files')
    p.add_argument('--n-procs', type=int, default=1,
                   help='worker processes for extract, each with its own butler (default 1)')

    ## build
    p.add_argument('--refmag-model', choices=MODELS, default='cmodel',
                   help='total magnitude for refmag (default cmodel)')
    p.add_argument('--color-model', choices=MODELS, default='sersic',
                   help='color-optimized magnitude for mag (default sersic, since GAaP is out)')
    p.add_argument('--mag-cut', type=float, default=24.0,
                   help='refmag cut (default 24, the bright_check from the intro notebook)')
    p.add_argument('--dust-coeffs', default=None,
                   help="A_band/E(B-V) per band, e.g. 'r=..,i=..,z=..' (no default). same values as the maps")

    return p.parse_args()


def main():
    ## reading args and finding the reference band
    args = parse_args()
    bands = list(args.bands)
    if args.ref_band not in bands:
        sys.exit(f'reference band {args.ref_band} is not in {args.bands}')
    ref_i = bands.index(args.ref_band)

    ## running the step asked for
    if args.step == 'extract':
        run_extract(args, bands, ref_i)
    else:
        run_build(args, bands, ref_i)


if __name__ == '__main__':
    main()
