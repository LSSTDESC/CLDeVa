'''
redmapper_dp2_common.py

Shared pieces for the DP2 redMaPPer inputs: butler access, magnitudes, healpix helpers, and output filenames
'''

import os

import numpy as np
import healpy as hp

from lsst.daf.butler import Butler
from lsst.geom import SpherePoint, degrees

'''
Initial definitions
'''
REPO = 'dp2'
COLLECTION = 'dp2'
SKYMAP = 'lsst_cells_v2'

## nJy -> AB
njy_ab_zeropoint = 2.5 * np.log10(3631e9)

## the maps give a 5 sigma point source depth, redmapper wants a 10 sigma galaxy
## depth. going 5 -> 10 sigma is 2.5*log10(2) = 0.753 mag by definition, and we
## measured 0.71-0.76 in griz on tract 10188, so no extra extended-source penalty
SIG_OFFSET = 2.5 * np.log10(10.0 / 5.0)


'''
butler
'''
def open_butler():
    '''
    Open the DP2 butler and its skymap, same setup as CL_DP2_Intro_IO.ipynb.

    Returns
    -------
    butler: `lsst.daf.butler.Butler`
       Butler on the dp2 repo, dp2 collection.
    skymap: `lsst.skymap.BaseSkyMap`
       The lsst_cells_v2 skymap.
    '''
    butler = Butler(REPO, collections=[COLLECTION])
    skymap = butler.get('skyMap', skymap=SKYMAP)
    return butler, skymap


def list_tracts(butler):
    '''
    Every tract that has an object table.

    Parameters
    ----------
    butler: `lsst.daf.butler.Butler`
       Butler to query.

    Returns
    -------
    tracts: `np.ndarray`
       Sorted tract ids.
    '''
    refs = butler.registry.queryDatasets('object', skymap=SKYMAP)
    return np.array(sorted({int(ref.dataId['tract']) for ref in refs}))


def read_tract(butler, tract, columns):
    '''
    Read one tract of the object table, only the columns asked for.

    The table is ~1250 columns wide, so this matters a lot.

    Parameters
    ----------
    butler: `lsst.daf.butler.Butler`
       Butler to read from.
    tract: `int`
       Tract id.
    columns: `list` [`str`]
       Column names to read.

    Returns
    -------
    table: `astropy.table.Table`
       The requested columns for that tract.
    '''
    return butler.get('object', skymap=SKYMAP, tract=int(tract),
                      parameters={'columns': list(columns)})


def get_col(table, name, fill):
    '''
    Plain numpy array out of a (possibly masked) table column.

    Parameters
    ----------
    table: `astropy.table.Table`
       Table to read from.
    name: `str`
       Column name.
    fill: `object`
       Value to substitute where the column is masked.

    Returns
    -------
    values: `np.ndarray`
       Unmasked column values.
    '''
    col = table[name]
    if hasattr(col, 'filled'):
        col = col.filled(fill)
    return np.asarray(col)


def get_tract_list(args_tracts, args_tract_list, butler):
    '''
    Tracts to run on: --tracts 10188,8228, --tract-list file.txt, or everything

    Parameters
    ----------
    args_tracts: `str` or `None`
       Comma separated tract ids.
    args_tract_list: `str` or `None`
       Path to a text file with one tract per line.
    butler: `lsst.daf.butler.Butler`
       Butler, used only when neither argument is given.

    Returns
    -------
    tracts: `np.ndarray`
       Tract ids to process.
    '''
    if args_tracts is not None:
        return np.array([int(t) for t in args_tracts.split(',')])
    if args_tract_list is not None:
        return np.loadtxt(args_tract_list, dtype=int, ndmin=1)
    return list_tracts(butler)


'''
tract overlaps
'''
def in_tract_inner(skymap, tract, ra, dec):
    '''
    True where a position belongs to this tract according to the skymap

    Parameters
    ----------
    skymap: `lsst.skymap.BaseSkyMap`
       Skymap that defines tract ownership.
    tract: `int`
       Tract id being processed.
    ra: `np.ndarray`
       Right ascension, degrees.
    dec: `np.ndarray`
       Declination, degrees.

    Returns
    -------
    owned: `np.ndarray`
       Boolean, True where this tract owns the position.
    '''
    ra = np.asarray(ra, dtype=float)
    dec = np.asarray(dec, dtype=float)

    if hasattr(skymap, 'findTractIdArray'):
        owner = skymap.findTractIdArray(ra, dec, degrees=True)
    else:
        owner = np.array([skymap.findTract(SpherePoint(r, d, degrees)).getId()
                          for r, d in zip(ra, dec)])

    return np.asarray(owner) == int(tract)


'''
magnitudes
'''
def flux_to_mag(flux):
    '''
    Convert nanojansky flux to AB magnitude.

    Parameters
    ----------
    flux: `np.ndarray`
       Flux in nJy.

    Returns
    -------
    mag: `np.ndarray`
       AB magnitude, NaN where the flux is not positive.
    '''
    flux = np.asarray(flux, dtype=float)
    mag = np.full(flux.shape, np.nan)
    positive = flux > 0
    mag[positive] = -2.5 * np.log10(flux[positive]) + njy_ab_zeropoint
    return mag


def flux_to_magerr(flux, flux_err):
    '''
    Magnitude error from flux and its error 2.5/ln(10) / SNR

    Parameters
    ----------
    flux: `np.ndarray`
       Flux in nJy.
    flux_err: `np.ndarray`
       Flux error in nJy.

    Returns
    -------
    mag_err: `np.ndarray`
       Magnitude error, NaN where the flux is not positive.
    '''
    flux = np.asarray(flux, dtype=float)
    flux_err = np.asarray(flux_err, dtype=float)
    err = np.full(flux.shape, np.nan)
    positive = flux > 0
    err[positive] = (2.5 / np.log(10)) * flux_err[positive] / flux[positive]
    return err


'''
healpix
'''
def pix_area(nside):
    '''
    Area of one healpix pixel.

    Parameters
    ----------
    nside: `int`
       Healpix nside.

    Returns
    -------
    area: `float`
       Pixel area in square degrees.
    '''
    return hp.nside2pixarea(nside, degrees=True)


def pixel_centers(nside, pixels):
    '''
    Sky positions of NEST pixel centers (healsparse is NEST).

    Parameters
    ----------
    nside: `int`
       Healpix nside.
    pixels: `np.ndarray`
       NEST pixel indices.

    Returns
    -------
    ra: `np.ndarray`
       Right ascension, degrees.
    dec: `np.ndarray`
       Declination, degrees.
    '''
    return hp.pix2ang(nside, pixels, lonlat=True, nest=True)


def nest_children(parents, nside_parent, nside_child):
    '''
    Every NEST child pixel of a set of coarse NEST pixels

    NEST indexing

    Parameters
    ----------
    parents: `np.ndarray`
       Coarse NEST pixel indices.
    nside_parent: `int`
       Nside of the coarse pixels.
    nside_child: `int`
       Nside of the fine pixels, a power of two times ``nside_parent``.

    Returns
    -------
    children: `np.ndarray`
       Fine NEST pixel indices covering the parents.
    '''
    ratio = (nside_child // nside_parent) ** 2
    parents = np.asarray(parents, dtype=np.int64)
    return (parents[:, None] * ratio + np.arange(ratio, dtype=np.int64)).ravel()


def ring_hpix(pixels, nside, nside_config):
    '''
    RING pixels at ``nside_config`` that touch the footprint

    redmapper wants RING here, healsparse is NEST

    Parameters
    ----------
    pixels: `np.ndarray`
       Footprint pixels, NEST at ``nside``.
    nside: `int`
       Nside of ``pixels``.
    nside_config: `int`
       Nside of the RING list redmapper wants.

    Returns
    -------
    hpix: `np.ndarray`
       Sorted unique RING pixel indices.
    '''
    ra, dec = pixel_centers(nside, pixels)
    return np.unique(hp.ang2pix(nside_config, ra, dec, lonlat=True, nest=False))


'''
filenames
'''
def counts_file(out_dir, tract, bands, nside):
    '''
    Object-based per-tract counts

    Parameters
    ----------
    out_dir: `str`
       Map directory.
    tract: `int` or `str`
       Tract id, or a glob pattern.
    bands: `str`
       Band string, e.g. 'riz'.
    nside: `int`
       Map nside.

    Returns
    -------
    path: `str`
       File path.
    '''
    return os.path.join(out_dir, 'counts', f'counts_tract{tract}_{bands}_nside{nside}.npz')


def cell_counts_file(out_dir, tract, bands):
    '''
    Per-tract visit counts per coadd cell

    Parameters
    ----------
    out_dir: `str`
       Map directory.
    tract: `int` or `str`
       Tract id, or a glob pattern.
    bands: `str`
       Band string, e.g. 'riz'.

    Returns
    -------
    path: `str`
       File path.
    '''
    return os.path.join(out_dir, 'cell_counts', f'cellcounts_tract{tract}_{bands}.npz')


def mask_file(out_dir, bands, vis_min, nside):
    '''
    The fracgood geometry mask redmapper reads

    Parameters
    ----------
    out_dir: `str`
       Map directory.
    bands: `str`
       Band string, e.g. 'riz'.
    vis_min: `int`
       Visits required per band.
    nside: `int`
       Map nside.

    Returns
    -------
    path: `str`
       File path.
    '''
    return os.path.join(out_dir, f'dp2_geometry_{bands}_vis{vis_min}_nside{nside}.hs')


def bitpacked_file(out_dir, bands, vis_min, nside_sparse):
    '''
    The cell-resolution geometry as bitpacked boolean

    Parameters
    ----------
    out_dir: `str`
       Map directory.
    bands: `str`
       Band string, e.g. 'riz'.
    vis_min: `int`
       Visits required per band.
    nside_sparse: `int`
       Cell-resolution nside.

    Returns
    -------
    path: `str`
       File path.
    '''
    return os.path.join(out_dir, f'dp2_geometry_{bands}_vis{vis_min}_nside{nside_sparse}_bitpacked.hs')


def depth_file(out_dir, bands, vis_min, nside, depth_band='min', dered=False):
    '''
    Depth map path.

    Parameters
    ----------
    out_dir: `str`
       Map directory.
    bands: `str`
       Band string, e.g. 'riz'.
    vis_min: `int`
       Visits required per band.
    nside: `int`
       Map nside.
    depth_band: `str`, optional
       'min' for the shallowest of the bands (the original file name), otherwise
       that one band.
    dered: `bool`, optional
       Whether the depth was dust corrected.

    Returns
    -------
    path: `str`
       File path.
    '''
    band = '' if depth_band == 'min' else f'{depth_band}band_'
    dust = '_dered' if dered else ''
    return os.path.join(out_dir, f'dp2_depth_{band}{bands}_vis{vis_min}{dust}_nside{nside}.hs')


def flux_file(out_dir, tract, bands, vis_min):
    '''
    Per-tract extracted fluxes.

    Parameters
    ----------
    out_dir: `str`
       Catalog directory.
    tract: `int` or `str`
       Tract id, or a glob pattern.
    bands: `str`
       Band string, e.g. 'riz'.
    vis_min: `int`
       Visits required per band in the geometry the fluxes were cut on.

    Returns
    -------
    path: `str`
       File path.
    '''
    return os.path.join(out_dir, 'tract_fluxes', f'flux_tract{tract}_{bands}_vis{vis_min}.fits')


def catalog_base(out_dir, ref_model, color_model, bands, vis_min, dered=False):
    '''
    Stem of the built galaxy catalog, e.g. dp2_photo_refcmodel_colsersic_riz_vis3.

    Parameters
    ----------
    out_dir: `str`
       Catalog directory.
    ref_model: `str`
       Photometry model used for refmag.
    color_model: `str`
       Photometry model used for the colors.
    bands: `str`
       Band string, e.g. 'riz'.
    vis_min: `int`
       Visits required per band in the geometry.
    dered: `bool`, optional
       Whether the magnitudes were dust corrected.

    Returns
    -------
    base: `str`
       Path stem, without the _NNNNNNN.fit suffix GalaxyCatalogMaker adds.
    '''
    dust = '_dered' if dered else ''
    return os.path.join(out_dir, f'dp2_photo_ref{ref_model}_col{color_model}{dust}_{bands}_vis{vis_min}')


'''
dust
'''
def parse_dust_coeffs(text, bands):
    '''
    Parse ``"r=...,i=...,z=..."`` into ``{band: A_band / E(B-V)}``.

    There's deliberately no default: the coefficients have to come from a real
    source for the LSST filters, not from a guess.

    Parameters
    ----------
    text: `str` or `None`
       Comma separated band=coefficient pairs, or None for no dust correction.
    bands: `list` [`str`]
       Bands that need a coefficient.

    Returns
    -------
    coeffs: `dict` [`str`, `float`] or `None`
       Coefficient per band, or None if ``text`` was None.

    Raises
    ------
    ValueError
       If any band is missing a coefficient.
    '''
    if text is None:
        return None
    coeffs = {}
    for item in text.split(','):
        band, value = item.split('=')
        coeffs[band.strip()] = float(value)
    missing = [b for b in bands if b not in coeffs]
    if missing:
        raise ValueError(f'--dust-coeffs is missing bands {missing}')
    return coeffs
