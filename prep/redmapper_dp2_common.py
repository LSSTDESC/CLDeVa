'''
redmapper_dp2_common.py

stuff all the DP2 redmapper scripts use: butler, flux -> mag, file names
'''

import os

import numpy as np

## the butler only exists in the LSST stack (maps + extract)
## build and specz run in desc_rm without it
try:
    from lsst.daf.butler import Butler
except ImportError:
    Butler = None

SKYMAP = 'lsst_cells_v2'
NJY_ZP = 2.5 * np.log10(3631e9)   ## nJy -> AB


'''
butler
'''
def open_butler():
    butler = Butler('dp2', collections=['dp2'])
    skymap = butler.get('skyMap', skymap=SKYMAP)
    return butler, skymap


def get_tracts(butler, tracts=None):
    ## the tracts given like '10188,8228', otherwise every tract with an object table
    if tracts is not None:
        return [int(t) for t in tracts.split(',')]
    refs = butler.registry.queryDatasets('object', skymap=SKYMAP)
    return sorted({int(ref.dataId['tract']) for ref in refs})


def owned_by_tract(skymap, tract, ra, dec):
    ## tracts overlap, so only keep what the skymap says belongs to this tract
    owner = skymap.findTractIdArray(np.asarray(ra, dtype=float), np.asarray(dec, dtype=float), degrees=True)
    return owner == tract


def get_col(table, name, fill):
    ## object table columns can be masked, fill those in
    col = table[name]
    if hasattr(col, 'filled'):
        col = col.filled(fill)
    return np.asarray(col)


'''
magnitudes
'''
def flux_to_mag(flux):
    ## NaN where the flux isn't positive
    flux = np.asarray(flux, dtype=float)
    mag = np.full(flux.shape, np.nan)
    good = flux > 0
    mag[good] = -2.5 * np.log10(flux[good]) + NJY_ZP
    return mag


def flux_to_magerr(flux, flux_err):
    ## 2.5/ln(10) / SNR, NaN where the flux isn't positive
    flux = np.asarray(flux, dtype=float)
    flux_err = np.asarray(flux_err, dtype=float)
    err = np.full(flux.shape, np.nan)
    good = flux > 0
    err[good] = (2.5 / np.log(10)) * flux_err[good] / flux[good]
    return err


'''
file names. name is the run, e.g. 'griz_vis2_bsouth20_star16'
'''
## dust correction: A_band = coefficient * E(B-V), with E(B-V) = 0.86 x CSFD (dp2utils DustCorrector)
DUST_COEFFS = {'g': 3.6234, 'r': 2.5488, 'i': 1.8842, 'z': 1.4921}
CSFD_DIR = '/global/cfs/cdirs/lsst/utilities/data-registry/lsst_desc_working/user/s-sayan/dust_maps/csfd/'   ## DataRegistry id 145


def run_name(bands, vis_min, star=None):
    ## always b < -20 now. star is the Gaia G cut of the star mask (12, 15, 16) or None for no star mask
    star = 'nostar' if star is None else f'star{star}'
    return f'{bands}_vis{vis_min}_bsouth20_{star}'


def geometry_file(map_dir, name, nside=4096):
    return os.path.join(map_dir, f'dp2_geometry_{name}_nside{nside}.hs')


def cells_file(map_dir, name, nside_sparse=131072):
    return os.path.join(map_dir, f'dp2_geometry_{name}_nside{nside_sparse}_bitpacked.hs')


def depth_file(map_dir, name, band='min', nside=4096):
    ## 'min' is the shallowest band, otherwise that band. dust corrected
    band = '' if band == 'min' else f'{band}band_'
    return os.path.join(map_dir, f'dp2_depth_{band}{name}_nside{nside}.hs')


def cell_counts_file(map_dir, tract, bands):
    return os.path.join(map_dir, 'cell_counts', f'cellcounts_tract{tract}_{bands}.npz')


def star_mask_file(star_dir, star):
    ## True = masked, Gaia stars brighter than G = star
    return os.path.join(star_dir, f'dp2_gaia_star_mask_g{star}_nside131072.hs')


def flux_file(flux_dir, tract):
    ## one set of flux files (griz) for every setup
    return os.path.join(flux_dir, f'flux_tract{tract}.fits')


def catalog_base(cat_dir, name):
    ## redmapper adds _<hpix>.fit and _master_table.fit to this
    return os.path.join(cat_dir, f'dp2_photo_refcmodel_colsersic_dered_{name}')


def specz_file(specz_dir, name, red=False):
    red = 'red_' if red else ''
    return os.path.join(specz_dir, f'dp2_specz_desi_dr1_{red}{name}.fit')
