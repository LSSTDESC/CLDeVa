'''
create_redmapper_dp2_specz.py

spec-z training sample for redmapper on DP2 from DESI DR1 (desc_rm, nersc/env_redmapper.sh)

select: every clean DESI galaxy spec-z in our footprint
red: match those to our galaxies, fit red + blue gaussians to the color offset from the redmapper template in z bins, keep the red ones (p_red > 0.5)

python create_redmapper_dp2_specz.py select --bands griz --vis-min 2 --star 16 --map-dir .../maps --flux-dir .../fluxes --specz-dir .../specz

python create_redmapper_dp2_specz.py red --bands griz --vis-min 2 --star 16 --map-dir .../maps --flux-dir .../fluxes --specz-dir .../specz

the red fit runs over z 0 - 1.1 for every setup (above 1.1 the red gaussian doesn't do well)
'''

import sys; import os
import glob
import argparse

import numpy as np
import pandas as pd
import healsparse as hsp
import fitsio

from astropy.coordinates import SkyCoord
import astropy.units as u

import redmapper
from redmapper.fitters import EcgmmFitter

import redmapper_dp2_common as common
from create_redmapper_dp2_catalogs import fluxes_to_galaxies

'''
Initial definitions
'''
DESI_FILE = '/global/cfs/cdirs/desi/public/dr1/spectro/redux/iron/zcatalog/v1/zall-pix-iron.fits'
DESI_COLUMNS = ['TARGET_RA', 'TARGET_DEC', 'Z', 'ZERR', 'ZWARN', 'ZCAT_PRIMARY', 'SPECTYPE']
CHUNK = 2_000_000   ## the file is ~28M rows

ZERR_MAX = 0.001
MATCH_RADIUS = 3.0  ## arcsec
DZ = 0.05 ## z bins for the red/blue fit
ZMIN = 0.0
ZMAX = 1.1
P_RED_CUT = 0.5

## which color to fit at each z (depending on if GRIZ or RIZ)
COLORS_GRIZ = [(0.0, 'g', 'r'), (0.35, 'r', 'i'), (0.7, 'i', 'z')]
COLORS_RIZ = [(0.0, 'r', 'i'), (0.7, 'i', 'z')]

## red/blue fit
DELTA_RANGE = (-1.0, 0.4)

## blue broad + below the red sequence, red narrow + near the template
FIT_OFFSET = 2.0
FIT_BOUNDS = [(0.01, 0.99),                            # blue weight
              (-1.0 + FIT_OFFSET, -0.1 + FIT_OFFSET),  # mu blue
              (-0.2 + FIT_OFFSET, 0.15 + FIT_OFFSET),  # mu red (-0.2 since the dust correction makes colors ~0.04 bluer)
              (0.05, 0.5),                             # sigma blue
              (0.01, 0.08)]                            # sigma red

spec_dtype = [('ra', 'f8'), ('dec', 'f8'), ('z', 'f4'), ('z_err', 'f4')]


'''
select
'''
def select(name, map_dir, specz_dir):
    ## DESI galaxies inside our footprint
    geometry = hsp.HealSparseMap.read(common.geometry_file(map_dir, name))

    kept = []
    with fitsio.FITS(DESI_FILE) as f:
        n_rows = f[1].get_nrows()
        for start in range(0, n_rows, CHUNK):
            d = f[1][DESI_COLUMNS][start:start + CHUNK]
            spectype = np.char.strip(d['SPECTYPE'].astype(str))
            d = d[d['ZCAT_PRIMARY'].astype(bool) & (d['ZWARN'] == 0) & (spectype == 'GALAXY')]
            d = d[geometry.get_values_pos(d['TARGET_RA'], d['TARGET_DEC'], lonlat=True) > 0]

            spec = np.zeros(len(d), dtype=spec_dtype)
            spec['ra'] = d['TARGET_RA']
            spec['dec'] = d['TARGET_DEC']
            spec['z'] = d['Z']
            spec['z_err'] = d['ZERR']
            kept.append(spec)
            print(f'  [{min(start + CHUNK, n_rows)}/{n_rows}] {sum(len(k) for k in kept)} in footprint', flush=True)

    spec = np.concatenate(kept)
    os.makedirs(specz_dir, exist_ok=True)
    fitsio.write(common.specz_file(specz_dir, name), spec, clobber=True)
    print(f'wrote {len(spec)} spec-z, z {spec["z"].min():.2f} to {spec["z"].max():.2f}')
    print('Done!')


'''
red
'''
def match(spec, name, bands, map_dir, flux_dir):
    ## nearest of our galaxies within MATCH_RADIUS one tract file at a time
    cells = hsp.HealSparseMap.read(common.cells_file(map_dir, name))
    spec_coords = SkyCoord(spec['ra'] * u.deg, spec['dec'] * u.deg)
    files = sorted(glob.glob(common.flux_file(flux_dir, '*')))

    matched = []
    for n, path in enumerate(files):
        gals = fluxes_to_galaxies(fitsio.read(path), bands) ## same cuts as the catalog
        gals = gals[cells.get_values_pos(gals['ra'], gals['dec'], lonlat=True)]
        if len(gals) == 0:
            continue

        ## only the spec-z near this tract; tracts across ra = 0 get wrapped
        ra_g = gals['ra'].copy()
        ra_s = spec['ra'].copy()
        if ra_g.max() - ra_g.min() > 180:
            ra_g[ra_g > 180] -= 360
            ra_s[ra_s > 180] -= 360
        near = ((ra_s > ra_g.min() - 0.01) & (ra_s < ra_g.max() + 0.01) &
                (spec['dec'] > gals['dec'].min() - 0.01) & (spec['dec'] < gals['dec'].max() + 0.01))
        if near.sum() == 0:
            continue

        idx, sep, _ = spec_coords[near].match_to_catalog_sky(SkyCoord(gals['ra'] * u.deg, gals['dec'] * u.deg))
        good = sep.arcsec < MATCH_RADIUS
        if good.sum() == 0:
            continue

        s = spec[near][good]
        g = gals[idx[good]]
        df = pd.DataFrame({'ra': s['ra'], 'dec': s['dec'], 'z': s['z'], 'z_err': s['z_err'],
                           'refmag': g['refmag'], 'sep': sep.arcsec[good]})
        for k, b in enumerate(bands):
            df[f'mag_{b}'] = g['mag'][:, k]
            df[f'mag_err_{b}'] = g['mag_err'][:, k]
        matched.append(df)

        if (n + 1) % 200 == 0:
            print(f'  [{n+1}/{len(files)}] {sum(len(m) for m in matched)} matched', flush=True)

    matched = pd.concat(matched, ignore_index=True)
    ## the same spec-z can match in two tract files at the edges, keep the closer one
    matched = matched.sort_values('sep').drop_duplicates(subset=['ra', 'dec', 'z']).reset_index(drop=True)
    print(f'{len(matched)} spec-z matched, median sep {matched["sep"].median():.2f}"')
    return matched


def gauss(x, mu, sigma):
    return np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (np.sqrt(2 * np.pi) * sigma)


def red(name, bands, map_dir, flux_dir, specz_dir):
    spec = fitsio.read(common.specz_file(specz_dir, name))
    spec = spec.astype(spec.dtype.newbyteorder('='))
    spec = spec[spec['z_err'] < ZERR_MAX]
    print(f'{len(spec)} spec-z with z_err < {ZERR_MAX}')

    matched = match(spec, name, bands, map_dir, flux_dir)

    ## redmapper's starting red sequence
    path = os.path.join(os.path.dirname(redmapper.__file__), 'data', 'initcolors', 'bc03_colors_lsst.fit')
    template = fitsio.read(path, ext=1, lower=True)
    template_bands = fitsio.read_header(path, ext=1)['BANDS'].split(',')

    ## color + its error and template color for each galaxy, picking the color by z
    colors = COLORS_GRIZ if 'g' in bands else COLORS_RIZ
    z = matched['z'].values
    matched['color'] = np.full(len(z), np.nan, dtype=np.float32) ## float32 like the mags
    matched['color_err'] = np.full(len(z), np.nan, dtype=np.float32)
    matched['template'] = np.nan
    for z_start, blue, red_band in colors:
        sel = z >= z_start
        matched.loc[sel, 'color'] = (matched[f'mag_{blue}'] - matched[f'mag_{red_band}'])[sel]
        matched.loc[sel, 'color_err'] = np.hypot(matched[f'mag_err_{blue}'], matched[f'mag_err_{red_band}'])[sel]
        matched.loc[sel, 'template'] = np.interp(z[sel], template['z'], template['color'][:, template_bands.index(blue)])
    matched['delta'] = matched['color'] - matched['template']

    ## red vs blue in each z bin
    matched['p_red'] = np.nan
    z_edges = np.arange(ZMIN, ZMAX + 0.001, DZ)
    for lo, hi in zip(z_edges[:-1], z_edges[1:]):
        in_bin = (z >= lo) & (z < hi)
        delta = matched['delta'].values[in_bin]
        err = matched['color_err'].values[in_bin]
        if len(delta) < 20:
            print(f'z {lo:.2f}: only {len(delta)}, skipping')
            continue

        ok = (delta > DELTA_RANGE[0]) & (delta < DELTA_RANGE[1])
        wt, mu, sigma = EcgmmFitter(delta[ok], err[ok]).fit([0.4], [-0.3, 0.0], [0.15, 0.03],
                                                             bounds=FIT_BOUNDS, offset=FIT_OFFSET)

        like_blue = wt[0] * gauss(delta, mu[0], np.hypot(sigma[0], err))
        like_red = wt[1] * gauss(delta, mu[1], np.hypot(sigma[1], err))
        p_red = like_red / (like_red + like_blue)
        p_red[~ok] = 0.0
        matched.loc[in_bin, 'p_red'] = p_red
        print(f'z {lo:.2f}: n {len(delta)}, red fraction {wt[1]:.2f}, mu blue/red {mu[0]:.3f}/{mu[1]:.3f}, '
              f'sigma blue/red {sigma[0]:.3f}/{sigma[1]:.3f}')

    reds = matched[matched['p_red'] > P_RED_CUT]
    out = np.zeros(len(reds), dtype=spec_dtype)
    for c in ['ra', 'dec', 'z', 'z_err']:
        out[c] = reds[c].values
    fitsio.write(common.specz_file(specz_dir, name, red=True), out, clobber=True)
    print(f'wrote {len(out)} red spec-z')

    ## the whole matched table, for plots
    matched.to_parquet(os.path.join(specz_dir, f'dp2_specz_matched_{name}.parquet'), index=False)
    print('Done!')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('step', choices=['select', 'red'])
    p.add_argument('--bands', default='griz')
    p.add_argument('--vis-min', type=int, default=2)
    p.add_argument('--star', type=int, default=None, help='star mask G cut (12, 15, 16), leave out for none')
    p.add_argument('--map-dir', required=True)
    p.add_argument('--flux-dir', required=True)
    p.add_argument('--specz-dir', required=True)
    args = p.parse_args()

    name = common.run_name(args.bands, args.vis_min, args.star)
    if args.step == 'select':
        select(name, args.map_dir, args.specz_dir)
    else:
        red(name, list(args.bands), args.map_dir, args.flux_dir, args.specz_dir)
