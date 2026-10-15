'''
create_redmapper_dp2_catalogs.py

galaxy catalog for redmapper on DP2, two steps:

extract (LSST stack, nersc/env_butler.sh): object table -> one griz flux file per tract, done once for every setup.
    cuts on the biggest footprint (riz, 2 visits, no star mask) so every setup can be built from the same files
    python create_redmapper_dp2_catalogs.py extract --cells .../riz_vis2/maps/dp2_geometry_riz_vis2_bsouth20_nostar_nside131072_bitpacked.hs --flux-dir .../fluxes

build (desc_rm, nersc/env_redmapper.sh): flux files -> dust corrected mags -> redmapper's pixel files + master table
    python create_redmapper_dp2_catalogs.py build --bands griz --vis-min 2 --star 16 --map-dir .../maps --cat-dir .../catalogs --flux-dir .../fluxes
'''

import sys; import os
import glob
import argparse
import multiprocessing as mp

import numpy as np
import healpy as hp
import healsparse as hsp
import fitsio

from astropy.coordinates import SkyCoord
import astropy.units as u
from dustmaps.csfd import CSFDQuery

import redmapper_dp2_common as common

## only build needs redmapper, and the LSST stack (extract) doesn't have it
try:
    import redmapper
except ImportError:
    redmapper = None

'''
Initial definitions
'''
FLUX_BANDS = 'griz' ## what the flux files have. riz setups just use the last three
REF_BAND = 'z' ## shallowest band
MAG_CUT = 24.0 ## refmag cut, after dust correction
EXTRACT_MAG_CUT = 24.5 
ZEROPOINT = 22.5 ## same as cosmoDC2
NSIDE = 4096

MODELS = ['cmodel', 'sersic']
FLUX_COLS = {'cmodel': '{b}_cModelFlux', 'sersic': '{b}_sersicFlux'}

EXT_COL = 'griz_model_extendedness'


'''
extract
'''
def start_worker(cells_path):
    ## setting up workers for parallelization
    ## need to make butler, skymap, cells, and the dust map global for this
    global butler, skymap, cells, csfd
    butler, skymap = common.open_butler()
    cells = hsp.HealSparseMap.read(cells_path)
    csfd = CSFDQuery(map_fname=common.CSFD_DIR + 'csfd_ebv.fits', mask_fname=common.CSFD_DIR + 'mask.fits')


def extract_tract(tract):
    ## one tract of the object table -> flux file with both models and all of griz
    path = common.flux_file(flux_dir, tract)

    ## skips but doesn't quit tracts if there's an error
    try:
        cols = ['objectId', 'coord_ra', 'coord_dec', EXT_COL, 'sersic_no_data_flag', 'sersic_unknown_flag']
        for model in MODELS:
            cols += [FLUX_COLS[model].format(b=b) for b in FLUX_BANDS] + [FLUX_COLS[model].format(b=b) + 'Err' for b in FLUX_BANDS]
        for b in FLUX_BANDS:
            cols += [f'{b}_pixelFlags_saturatedCenter', f'{b}_pixelFlags_nodata']
        d = butler.get('object', skymap=common.SKYMAP, tract=tract, parameters={'columns': cols})
    except Exception as err:
        #rerunning picks these tracts up again
        return f'tract {tract} FAILED: {err}'

    ra = common.get_col(d, 'coord_ra', np.nan).astype(float)
    dec = common.get_col(d, 'coord_dec', np.nan).astype(float)

    ## in cells and owned by this tract
    ok = np.isfinite(ra) & np.isfinite(dec)
    keep = np.zeros(len(ra), dtype=bool)
    keep[ok] = cells.get_values_pos(ra[ok], dec[ok], lonlat=True) > 0
    keep[keep] = common.owned_by_tract(skymap, tract, ra[keep], dec[keep])

    ## extended
    keep &= common.get_col(d, EXT_COL, np.nan).astype(float) > 0.5

    ## brighter than the (loose) cut in the reference band in either model
    bright = np.zeros(len(ra), dtype=bool)
    for model in MODELS:
        ref_mag = common.flux_to_mag(common.get_col(d, FLUX_COLS[model].format(b=REF_BAND), np.nan))
        bright |= np.isfinite(ref_mag) & (ref_mag < EXTRACT_MAG_CUT)
    keep &= bright

    ## setting up rm formatting. bad_pix is saturated or no data center pixel, per band (build only checks the bands it uses)
    nb = len(FLUX_BANDS)
    dtype = [('id', 'i8'), ('ra', 'f8'), ('dec', 'f8'), ('ebv', 'f4'), ('sersic_ok', 'i2'), ('bad_pix', 'i2', nb)]
    for model in MODELS:
        dtype += [(f'{model}_flux', 'f4', nb), (f'{model}_flux_err', 'f4', nb)]
    out = np.zeros(keep.sum(), dtype=dtype)
    out['id'] = common.get_col(d, 'objectId', -1)[keep]
    out['ra'] = ra[keep]
    out['dec'] = dec[keep]
    sersic_ok = ~common.get_col(d, 'sersic_no_data_flag', True).astype(bool) & ~common.get_col(d, 'sersic_unknown_flag', True).astype(bool)
    out['sersic_ok'] = sersic_ok[keep]
    out['bad_pix'] = np.column_stack([(common.get_col(d, f'{b}_pixelFlags_saturatedCenter', True).astype(bool)
                                       | common.get_col(d, f'{b}_pixelFlags_nodata', True).astype(bool))[keep] for b in FLUX_BANDS])
    for model in MODELS:
        out[f'{model}_flux'] = np.column_stack([common.get_col(d, FLUX_COLS[model].format(b=b), np.nan)[keep] for b in FLUX_BANDS])
        out[f'{model}_flux_err'] = np.column_stack([common.get_col(d, FLUX_COLS[model].format(b=b) + 'Err', np.nan)[keep] for b in FLUX_BANDS])

    ## an empty marker file for tracts with nothing in them
    if len(out) == 0:
        open(path + '.empty', 'w').close()
    else:
        out['ebv'] = 0.86 * csfd(SkyCoord(out['ra'] * u.deg, out['dec'] * u.deg))   ## E(B-V) the dust correction uses
        fitsio.write(path, out, clobber=True)
    return f'tract {tract}: {len(ra)} objects -> {len(out)} kept'


'''
build the files
'''
def fluxes_to_galaxies(flux, bands):
    ## one tract's fluxes -> dust corrected redmapper galaxies in these bands. also used by the specz script
    idx = [FLUX_BANDS.index(b) for b in bands]
    ref = FLUX_BANDS.index(REF_BAND)
    with np.errstate(invalid='ignore', divide='ignore'):
        mag = common.flux_to_mag(flux['sersic_flux'][:, idx])
        mag_err = common.flux_to_magerr(flux['sersic_flux'][:, idx], flux['sersic_flux_err'][:, idx])
        refmag = common.flux_to_mag(flux['cmodel_flux'][:, ref])
        refmag_err = common.flux_to_magerr(flux['cmodel_flux'][:, ref], flux['cmodel_flux_err'][:, ref])

    ## dust correction
    for k, b in enumerate(bands):
        mag[:, k] -= common.DUST_COEFFS[b] * flux['ebv']
    refmag -= common.DUST_COEFFS[REF_BAND] * flux['ebv']

    ## checking Sersic, refmag (after the dust correction) and the pixel flags in our bands
    with np.errstate(invalid='ignore'):
        keep = np.isfinite(refmag) & (refmag < MAG_CUT) & flux['sersic_ok'].astype(bool)
        keep &= np.all(flux['bad_pix'][:, idx] == 0, axis=1)

        ## redmapper refuses the whole catalog if any galaxy has a non-finite, >= 90 or <= 0 mag or error
        keep &= np.all(np.isfinite(mag), axis=1) & np.all(np.isfinite(mag_err), axis=1) & np.isfinite(refmag_err)
        keep &= np.all(mag < 90, axis=1) & np.all(mag_err < 90, axis=1)
        keep &= np.all(mag > 0, axis=1) & np.all(mag_err > 0, axis=1)
        keep &= (flux['ra'] >= 0) & (flux['ra'] <= 360) & (flux['dec'] >= -90) & (flux['dec'] <= 90)

    galaxies = np.zeros(keep.sum(), dtype=[('id', 'i8'), ('ra', 'f8'), ('dec', 'f8'), ('refmag', 'f4'), ('refmag_err', 'f4'),
                                           ('mag', 'f4', len(bands)), ('mag_err', 'f4', len(bands)), ('ebv', 'f4')])
    galaxies['id'] = flux['id'][keep]
    galaxies['ra'] = flux['ra'][keep]
    galaxies['dec'] = flux['dec'][keep]
    galaxies['refmag'] = refmag[keep]
    galaxies['refmag_err'] = refmag_err[keep]
    galaxies['mag'] = mag[keep]
    galaxies['mag_err'] = mag_err[keep]
    galaxies['ebv'] = flux['ebv'][keep]
    return galaxies


def build(name, bands, map_dir, cat_dir, flux_dir):
    os.makedirs(cat_dir, exist_ok=True)
    base = common.catalog_base(cat_dir, name)
    if len(glob.glob(base + '_*.fit')) > 0:
        sys.exit(f'there is already a catalog at {base}_*.fit, delete it first (redmapper would add to it)')

    ## the cells for this setup + star mask (the flux files were cut on the biggest footprint)
    cells = hsp.HealSparseMap.read(common.cells_file(map_dir, name))

    ## area from the fracgood mask
    mask = hsp.HealSparseMap.read(common.geometry_file(map_dir, name))
    area = float(mask.get_values_pix(mask.valid_pixels).sum()) * hp.nside2pixarea(NSIDE, degrees=True)

    info = {'LIM_REF': MAG_CUT, 'REF_IND': bands.index(REF_BAND), 'AREA': area, 'NMAG': len(bands), 'MODE': 'LSST', 'ZP': ZEROPOINT}
    for i, b in enumerate(bands):
        info[f'{b.upper()}_IND'] = i
    maker = redmapper.GalaxyCatalogMaker(base, info)

    files = sorted(glob.glob(common.flux_file(flux_dir, '*')))
    print(f'{len(files)} tract files, {area:.1f} deg^2 -> {base}')
    n_gal = 0
    for n, path in enumerate(files):
        galaxies = fluxes_to_galaxies(fitsio.read(path), bands)
        galaxies = galaxies[cells.get_values_pos(galaxies['ra'], galaxies['dec'], lonlat=True)]
        if len(galaxies) > 0:
            maker.append_galaxies(galaxies)
            n_gal += len(galaxies)
        if (n + 1) % 200 == 0:
            print(f'  [{n+1}/{len(files)}] {n_gal} galaxies', flush=True)

    maker.finalize_catalog()
    print(f'{n_gal} galaxies')
    print('Done!')


if __name__ == '__main__':
    ## for running/setting args for different runs
    p = argparse.ArgumentParser()
    p.add_argument('step', choices=['extract', 'build'])
    p.add_argument('--flux-dir', required=True)
    p.add_argument('--cells', help='extract: the cell map to cut on')
    p.add_argument('--tracts', default=None, help='extract: e.g. 10188,8228 for testing (default all)')
    p.add_argument('--n-procs', type=int, default=4)
    p.add_argument('--bands', default='griz', help='build')
    p.add_argument('--vis-min', type=int, default=2, help='build')
    p.add_argument('--star', type=int, default=None, help='build: star mask G cut (12, 15, 16), leave out for none')
    p.add_argument('--map-dir', help='build')
    p.add_argument('--cat-dir', help='build')
    args = p.parse_args()

    flux_dir = args.flux_dir

    if args.step == 'extract':
        os.makedirs(flux_dir, exist_ok=True)

        ## skipping tracts that are already done (delete the files to redo them)
        butler, skymap = common.open_butler()
        tracts = common.get_tracts(butler, args.tracts)
        todo = [t for t in tracts if not os.path.exists(common.flux_file(flux_dir, t))
                and not os.path.exists(common.flux_file(flux_dir, t) + '.empty')]
        print(f'extracting {len(todo)} tracts ({len(tracts) - len(todo)} already done), cut on {args.cells}')
        with mp.get_context('fork').Pool(args.n_procs, initializer=start_worker, initargs=(args.cells,)) as pool:
            for n, message in enumerate(pool.imap_unordered(extract_tract, todo)):
                print(f'[{n+1}/{len(todo)}] {message}', flush=True)
        print('Done!')
    else:
        name = common.run_name(args.bands, args.vis_min, args.star)
        build(name, list(args.bands), args.map_dir, args.cat_dir, flux_dir)
