'''
create_redmapper_dp2_maps.py

Geometry mask + depth map for redMaPPer on DP2 on the RSP.
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
from dustmaps.sfd import SFDQuery

import lsst.geom as geom
from lsst.skymap import Index2D

import redmapper_dp2_common as common

'''
Initial definitions
'''
MAGLIM_MAP = 'deepCoadd_psf_maglim_consolidated_map_weighted_mean'
EXPTIME_MAP = 'deepCoadd_exposure_time_consolidated_map_sum'

## cells: patches are a 22x22 grid, the inner (non-overlapping) region is 1-20
N_CELL = 22
CELL_INNER = (1, 21)

## a visit counts for a cell if its unmasked_fraction is above this
WEIGHT_CUT = 0.0 #from shear team's cut

## shear team drops this detector
BAD_DETECTOR = 122

## start method for the stage 1 workers
START_METHOD = 'fork'

depth_dtype = np.dtype([('exptime', 'f4'), ('limmag', 'f4'),
                        ('m50', 'f4'), ('fracgood', 'f4')])


'''
stage 1: visit counts per cell
'''
def coadd_refs(butler, tract, bands):
    '''
    The deep_coadd refs of one tract, keyed by (patch, band).

    One registry query per tract, so the loop below never has to ask whether a
    coadd exists -- if it's not in here, that patch/band wasn't made.

    Parameters
    ----------
    butler: `lsst.daf.butler.Butler`
       Butler to query.
    tract: `int`
       Tract id.
    bands: `list` [`str`]
       Bands to keep.

    Returns
    -------
    refs: `dict` [(`int`, `str`), `lsst.daf.butler.DatasetRef`]
       Ref per (patch, band).
    '''
    refs = butler.query_datasets('deep_coadd', skymap=common.SKYMAP, tract=int(tract))
    return {(int(ref.dataId['patch']), str(ref.dataId['band'])): ref
            for ref in refs if str(ref.dataId['band']) in bands}


def cell_visit_counts(prov, weight_cut):
    '''
    Visit counts per cell for one patch in one band.

    A contribution counts if it covered more than ``weight_cut`` of the cell with
    good pixels and didn't come from the bad detector. One bincount over the
    flattened cell index instead of 484 mask-and-sums.

    Parameters
    ----------
    prov: `lsst.cell_coadds.CoaddProvenance`
       Provenance component of one deep_coadd.
    weight_cut: `float`
       Minimum unmasked_fraction, exclusive, matching the shear group's ``>``.

    Returns
    -------
    counts: `np.ndarray`
       (22, 22) visit count per cell.
    '''
    contrib = prov.contributions
    good = np.asarray(contrib['unmasked_fraction']) > weight_cut
    good &= np.asarray(contrib['detector']) != BAD_DETECTOR

    flat = np.asarray(contrib['cell_i'])[good] * N_CELL + np.asarray(contrib['cell_j'])[good]
    counts = np.bincount(flat, minlength=N_CELL * N_CELL)
    return counts.reshape(N_CELL, N_CELL).astype(np.int16)


def n_patches(tract_info):
    '''
    How many patches a tract has.

    Parameters
    ----------
    tract_info: `lsst.skymap.TractInfo`
       Tract to measure.

    Returns
    -------
    n: `int`
       Patch count, 100 for lsst_cells_v2.
    '''
    n = tract_info.num_patches
    return int(n[0]) * int(n[1])


def tract_cell_counts(butler, skymap, tract, bands, weight_cut):
    '''
    Visit counts for every patch and band of one tract.

    Parameters
    ----------
    butler: `lsst.daf.butler.Butler`
       Butler to read provenance from.
    skymap: `lsst.skymap.BaseSkyMap`
       Skymap.
    tract: `int`
       Tract id.
    bands: `list` [`str`]
       Bands required in the footprint.
    weight_cut: `float`
       Minimum unmasked_fraction, exclusive.

    Returns
    -------
    patches: `np.ndarray`
       Patch ids that had at least one coadd, shape (n_patch,).
    counts: `np.ndarray`
       Visit counts, shape (n_patch, n_band, 22, 22).
    '''
    tract_info = skymap[int(tract)]
    refs = coadd_refs(butler, tract, bands)

    ## looping through patches and bands to get visit counts
    patches, counts = [], []
    for patch in range(n_patches(tract_info)):
        per_band = np.zeros((len(bands), N_CELL, N_CELL), dtype=np.int16)
        found = False
        for b_i, band in enumerate(bands):
            ref = refs.get((patch, band))
            if ref is None:
                continue ## no coadd in this band -> 0 visits -> cell fails
            prov = butler.get(ref.makeComponentRef('provenance'))
            per_band[b_i] = cell_visit_counts(prov, weight_cut)
            found = True
        if found:
            patches.append(patch)
            counts.append(per_band)
    ## returning empty arrays if no patch had a coadd
    if len(patches) == 0:
        return (np.zeros(0, dtype=np.int32),
                np.zeros((0, len(bands), N_CELL, N_CELL), dtype=np.int16))
    
    return np.array(patches, dtype=np.int32), np.array(counts, dtype=np.int16)


## setting up one butler per worker
_worker = {}


def init_worker(bands, weight_cut, out_dir):
    '''
    Open this worker's own butler and remember the run settings.

    Parameters
    ----------
    bands: `list` [`str`]
       Bands required in the footprint.
    weight_cut: `float`
       Minimum unmasked_fraction, exclusive.
    out_dir: `str`
       Map directory.
    '''
    _worker['butler'], _worker['skymap'] = common.open_butler()
    _worker['bands'] = bands
    _worker['weight_cut'] = weight_cut
    _worker['out_dir'] = out_dir


def count_worker(tract):
    '''
    Pool entry point: count one tract with this worker's butler.

    Parameters
    ----------
    tract: `int`
       Tract id.

    Returns
    -------
    message: `str`
       Progress line for the log.
    '''
    return count_and_save(_worker['butler'], _worker['skymap'], tract, _worker['bands'],
                          _worker['weight_cut'], _worker['out_dir'])


def count_and_save(butler, skymap, tract, bands, weight_cut, out_dir):
    '''
    Count one tract and write its file.

    Writes to a temporary name and renames, so a killed run never leaves a
    half-written .npz behind for the next one to trip over.

    Parameters
    ----------
    butler: `lsst.daf.butler.Butler`
       Butler to read provenance from.
    skymap: `lsst.skymap.BaseSkyMap`
       Skymap.
    tract: `int`
       Tract id.
    bands: `list` [`str`]
       Bands required in the footprint.
    weight_cut: `float`
       Minimum unmasked_fraction, exclusive.
    out_dir: `str`
       Map directory.

    Returns
    -------
    message: `str`
       Progress line for the log.
    '''
    path = common.cell_counts_file(out_dir, tract, ''.join(bands))

    try:
        patches, counts = tract_cell_counts(butler, skymap, tract, bands, weight_cut)
    except Exception as err:
        ## skipping a bad tract
        return f'tract {tract} FAILED: {err}'

    tmp = path + '.tmp.npz'
    np.savez(tmp, patches=patches, counts=counts, bands=np.array(bands),
             weight_cut=weight_cut, bad_detector=BAD_DETECTOR)
    os.replace(tmp, path)

    ## counting up, no cut applied
    inner = slice(*CELL_INNER)
    worst = counts.min(axis=1)[:, inner, inner] if len(patches) else np.zeros(1)
    return (f'tract {tract}: {len(patches)} patches, '
            f'{int((worst > 0).sum())} inner cells with visits in every band')


def run_cell_counts(butler, skymap, tracts, bands, weight_cut, out_dir, clobber, n_procs=1):
    '''
    Stage 1: save visit counts per cell for every tract.

    Reading provenance is nearly all butler latency, so workers barely contend
    for CPU and the speedup is close to linear in ``n_procs``.

    Parameters
    ----------
    butler: `lsst.daf.butler.Butler`
       Butler to read provenance from, used only when ``n_procs`` is 1.
    skymap: `lsst.skymap.BaseSkyMap`
       Skymap, used only when ``n_procs`` is 1.
    tracts: `np.ndarray`
       Tract ids to count.
    bands: `list` [`str`]
       Bands required in the footprint.
    weight_cut: `float`
       Minimum unmasked_fraction, exclusive.
    out_dir: `str`
       Map directory; files land in ``out_dir/cell_counts``.
    clobber: `bool`
       Recount tracts that already have a file.
    n_procs: `int`, optional
       Worker processes, each with its own butler.
    '''
    os.makedirs(os.path.join(out_dir, 'cell_counts'), exist_ok=True)

    ## skipping tracts that already have a file
    todo = [int(t) for t in tracts
            if clobber or not os.path.exists(common.cell_counts_file(out_dir, t, ''.join(bands)))]
    print(f'{len(todo)} tracts to count, {len(tracts) - len(todo)} already done, {n_procs} process(es)')

    ## counting one tract at a time
    if n_procs == 1:
        for n, tract in enumerate(todo):
            print(f'[{n+1}/{len(todo)}] '
                  + count_and_save(butler, skymap, tract, bands, weight_cut, out_dir), flush=True)
        return

    ## counting with a pool of workers
    with mp.get_context(START_METHOD).Pool(n_procs, initializer=init_worker,
                                           initargs=(bands, weight_cut, out_dir)) as pool:
        for n, message in enumerate(pool.imap_unordered(count_worker, todo)):
            print(f'[{n+1}/{len(todo)}] {message}', flush=True)


def load_cell_counts(out_dir, bands, weight_cut):
    '''
    Find the stage 1 files and check they were counted at the cut we want.

    Parameters
    ----------
    out_dir: `str`
       Map directory.
    bands: `list` [`str`]
       Bands required in the footprint.
    weight_cut: `float`
       Cut stage 2 expects the counts to have been made with.

    Returns
    -------
    files: `list` [`str`]
       Sorted cell count file paths.
    '''
    files = sorted(glob.glob(common.cell_counts_file(out_dir, '*', ''.join(bands))))
    if len(files) == 0:
        sys.exit(f'no cell count files in {out_dir}/cell_counts for {"".join(bands)}')

    ## checking the counts were made at the cut asked for
    with np.load(files[0]) as data:
        stored = float(data['weight_cut'])
    if stored != weight_cut:
        sys.exit(f'{files[0]} was counted at weight_cut {stored}, not {weight_cut}. '
                 f'rerun stage 1 with --clobber, or pass --weight-cut {stored}')

    print(f'loaded cell counts for {len(files)} tracts, weight_cut {stored}')
    return files


'''
good cells -> sky rectangles -> healsparse
'''
def cell_inner_bbox(patch_info, i, j):
    '''
    Pixel bounding box of one cell's inner region in tract coordinates

    Parameters
    ----------
    patch_info: `lsst.skymap.PatchInfo`
       Patch the cell belongs to.
    i: `int`
       Cell index along x, matching provenance ``cell_i``.
    j: `int`
       Cell index along y, matching provenance ``cell_j``.

    Returns
    -------
    bbox: `lsst.geom.Box2I`
       Inner bounding box in tract pixel coordinates.
    '''
    return patch_info.getCellInfo(Index2D(x=i, y=j)).inner_bbox


def owned_cells(skymap, tract_info, tract, patches):
    '''
    Which inner cells in this tract to prevent overlap

    Cell centers come out of the patch inner bbox arithmetically -- cell (i, j)
    
    starts at patch_min + ((i-1), (j-1))*cell_size

    Parameters
    ----------
    skymap: `lsst.skymap.BaseSkyMap`
       Skymap that defines tract ownership.
    tract_info: `lsst.skymap.TractInfo`
       Tract being rendered.
    tract: `int`
       Tract id.
    patches: `np.ndarray`
       Patch ids in the same order as the counts array.

    Returns
    -------
    owned: `np.ndarray`
       Boolean, shape (n_patch, 20, 20), indexed by (patch, i-1, j-1).
    '''
    ## setting up the cell center offsets
    n_inner = CELL_INNER[1] - CELL_INNER[0]
    offsets = np.arange(n_inner) + 0.5

    ## getting cell centers in tract pixels for every patch
    x = np.zeros((len(patches), n_inner, n_inner))
    y = np.zeros((len(patches), n_inner, n_inner))
    for p, patch in enumerate(patches):
        patch_info = tract_info[int(patch)]
        bbox = cell_inner_bbox(patch_info, CELL_INNER[0], CELL_INNER[0])
        origin = patch_info.getInnerBBox().getMin()
        x[p] = origin.getX() + offsets[:, None] * bbox.getWidth()
        y[p] = origin.getY() + offsets[None, :] * bbox.getHeight()

    ra, dec = tract_info.wcs.pixelToSkyArray(x.ravel(), y.ravel(), degrees=True)
    return common.in_tract_inner(skymap, tract, ra, dec).reshape(x.shape)


def block_min(worst):
    '''
    Worst visit count over each cell's 2x2 block

    A cell passes only if it and its i-1, j-1 neighbours all have enough visits

    Parameters
    ----------
    worst: `np.ndarray`
       Worst-band visit count per cell, shape (n_patch, 22, 22).

    Returns
    -------
    block: `np.ndarray`
       Minimum over the block, shape (n_patch, 20, 20), indexed by (patch, i-1, j-1).
    '''
    lo, hi = CELL_INNER[0] - 1, CELL_INNER[1] - 1
    return np.minimum.reduce([worst[:, lo:hi, lo:hi], worst[:, lo + 1:hi + 1, lo:hi],
                              worst[:, lo:hi, lo + 1:hi + 1], worst[:, lo + 1:hi + 1, lo + 1:hi + 1]])


def cell_runs(good):
    '''
    Merge neighbouring good cells in each row into runs for efficiency

    Parameters
    ----------
    good: `np.ndarray`
       Boolean per cell, shape (22, 22).

    Returns
    -------
    runs: `list` [(`int`, `int`, `int`)]
       (i, j_start, j_end) per run, inclusive.
    '''
    ## looping through rows to find runs of good cells
    runs = []
    for i in range(*CELL_INNER):
        j = CELL_INNER[0]
        while j < CELL_INNER[1]:
            ## skipping bad cells
            if not good[i, j]:
                j += 1
                continue
            j_start = j
            while j < CELL_INNER[1] and good[i, j]:
                j += 1
            runs.append((i, j_start, j - 1))
    return runs


def set_polygons(geometry, polygons, nside_sparse, chunk=200):
    '''
    Turn rectangles into pixels and set them True

    Parameters
    ----------
    geometry: `healsparse.HealSparseMap`
       Boolean map to fill.
    polygons: `list` [`healsparse.geom.Polygon`]
       Rectangles to render.
    nside_sparse: `int`
       Cell-resolution nside.
    chunk: `int`, optional
       Polygons expanded per pass.
    '''
    ## looping through the polygons in chunks
    for start in range(0, len(polygons), chunk):
        ## getting the pixels under each polygon
        pixels = [poly.get_pixels(nside=nside_sparse) for poly in polygons[start:start + chunk]]

        ## setting them all True at once
        geometry.update_values_pix(np.unique(np.concatenate(pixels)), True)


def clear_pixels(geometry, parents, nside, nside_sparse, chunk=5000):
    '''
    Set every fine pixel under these coarse pixels to False

    Parameters
    ----------
    geometry: `healsparse.HealSparseMap`
       Boolean map to punch holes in.
    parents: `np.ndarray`
       Coarse NEST pixels to clear.
    nside: `int`
       Nside of ``parents``.
    nside_sparse: `int`
       Cell-resolution nside.
    chunk: `int`, optional
       Parents expanded per pass.
    '''
    for start in range(0, len(parents), chunk):
        ## getting the fine pixels under this chunk
        children = common.nest_children(parents[start:start + chunk], nside, nside_sparse)
        geometry.update_values_pix(children, False)


def run_polygon(tract_wcs, patch_info, i, j_start, j_end):
    '''
    One run of cells -> a healsparse polygon on the sky

    Parameters
    ----------
    tract_wcs: `lsst.afw.geom.SkyWcs`
       Tract WCS.
    patch_info: `lsst.skymap.PatchInfo`
       Patch the run belongs to.
    i: `int`
       Cell row.
    j_start: `int`
       First cell in the run.
    j_end: `int`
       Last cell in the run, inclusive.

    Returns
    -------
    polygon: `healsparse.geom.Polygon`
       The run as a sky polygon.
    '''
    ## getting the pixel box from the first to the last cell
    box = geom.Box2D(cell_inner_bbox(patch_info, i, j_start))
    box.include(geom.Box2D(cell_inner_bbox(patch_info, i, j_end)))

    ## converting the corners to ra/dec
    corners = [tract_wcs.pixelToSky(corner) for corner in box.getCorners()]
    ra = np.array([c.getRa().asDegrees() for c in corners])
    dec = np.array([c.getDec().asDegrees() for c in corners])

    return hsp.geom.Polygon(ra=ra, dec=dec, value=True)


def build_geometry(skymap, files, bands, vis_min, nside_sparse, nside_cover,
                   use_block=True, clip_tract=True):
    '''
    Render every good cell of every tract into one bitpacked boolean map.

    Parameters
    ----------
    skymap: `lsst.skymap.BaseSkyMap`
       Skymap.
    files: `list` [`str`]
       Stage 1 cell count files.
    bands: `list` [`str`]
       Bands required in the footprint.
    vis_min: `int`
       Visits required per band, over the 2x2 block.
    nside_sparse: `int`
       Cell-resolution nside.
    nside_cover: `int`
       Coverage nside.
    use_block: `bool`, optional
       Apply the shear group's 2x2 block test rather than testing the cell alone.
    clip_tract: `bool`, optional
       Keep only cells this tract owns.

    Returns
    -------
    geometry: `healsparse.HealSparseMap`
       Boolean cell-resolution footprint.
    '''
    ## setting up the empty map
    geometry = hsp.HealSparseMap.make_empty(nside_cover, nside_sparse, bool, bit_packed=True)
    inner = slice(*CELL_INNER)

    ## looping through tracts to render the good cells
    n_runs = 0
    n_block = 0
    n_owned = 0
    for n, path in enumerate(files):
        ## getting the tract id and its counts
        tract = int(os.path.basename(path).split('_')[1].replace('tract', ''))
        with np.load(path) as data:
            patches, counts = data['patches'], data['counts']
        if len(patches) == 0:
            continue

        tract_info = skymap[tract]
        tract_wcs = tract_info.wcs

        ## worst band per cell, then the block test, then which tract
        worst = counts.min(axis=1)
        # testing the 2x2 block or the cell alone
        if use_block:
            passes = block_min(worst) >= vis_min
        else:
            passes = worst[:, slice(*CELL_INNER), slice(*CELL_INNER)] >= vis_min
        keep = passes
        if clip_tract:
            keep = passes & owned_cells(skymap, tract_info, tract, patches)
        n_block += int(passes.sum())
        n_owned += int(keep.sum())
        
        # turning the kept cells into polygons
        polygons = []
        for p, patch in enumerate(patches):
            if not keep[p].any():
                continue
            good = np.zeros((N_CELL, N_CELL), dtype=bool)
            good[inner, inner] = keep[p]
            patch_info = tract_info[int(patch)]
            for i, j_start, j_end in cell_runs(good):
                polygons.append(run_polygon(tract_wcs, patch_info, i, j_start, j_end))

        if len(polygons) > 0:
            set_polygons(geometry, polygons, nside_sparse)
            n_runs += len(polygons)

        if (n + 1) % 100 == 0:
            print(f'  [{n+1}/{len(files)}] {n_runs} cell runs rendered')

    print(f'cells passing the {vis_min} visit test        {n_block}')
    print(f'  of those, kept after the tract clip       {n_owned}')
    print(f'rendered {n_runs} cell runs from {len(files)} tracts')
    return geometry


def fracgood_map(geometry, nside):
    '''
    Degrade the cell-resolution geometry to a float fracgood per pixel.

    The fraction of each nside pixel inside the footprint, which is exactly what
    redMaPPer wants the mask to mean.

    Parameters
    ----------
    geometry: `healsparse.HealSparseMap`
       Boolean cell-resolution footprint.
    nside: `int`
       Output nside.

    Returns
    -------
    frac: `healsparse.HealSparseMap`
       Float fracgood map.
    '''
    return geometry.fracdet_map(nside)


'''
foregrounds and depth
'''
def foreground_cut(pixels, nside, ebv_max, b_min):
    '''
    Dust and galactic latitude cut at pixel centers.

    Parameters
    ----------
    pixels: `np.ndarray`
       NEST pixels at ``nside``.
    nside: `int`
       Map nside.
    ebv_max: `float`
       Maximum SFD E(B-V).
    b_min: `float`
       Minimum |galactic b| in degrees.

    Returns
    -------
    keep: `np.ndarray`
       Boolean, True for pixels that survive.
    ebv: `np.ndarray`
       SFD E(B-V) at every input pixel.
    '''
    ## getting pixel centers
    ra, dec = common.pixel_centers(nside, pixels)
    coords = SkyCoord(ra * u.deg, dec * u.deg)

    ## getting E(B-V) and galactic latitude
    ebv = SFDQuery()(coords)
    gal_b = coords.galactic.b.deg

    return (ebv < ebv_max) & (np.abs(gal_b) > b_min), ebv


def load_depth_maps(butler, bands, nside):
    '''
    Maglim and exposure time per band, degraded on read (never at native nside).

    Parameters
    ----------
    butler: `lsst.daf.butler.Butler`
       Butler to read from.
    bands: `list` [`str`]
       Bands to load.
    nside: `int`
       Nside to degrade to.

    Returns
    -------
    maglim: `dict` [`str`, `healsparse.HealSparseMap`]
       5 sigma point source depth per band.
    exptime: `dict` [`str`, `healsparse.HealSparseMap`]
       Summed exposure time per band.
    '''
    maglim, exptime = {}, {}
    for b in bands:
        maglim[b] = butler.get(MAGLIM_MAP, band=b, skymap=common.SKYMAP,
                               parameters={'degrade_nside': nside})
        exptime[b] = butler.get(EXPTIME_MAP, band=b, skymap=common.SKYMAP,
                                parameters={'degrade_nside': nside})
        print(f'  loaded {b} maps')
    return maglim, exptime


def depth_values(pixels, maglim, exptime, bands):
    '''
    Per-band 10 sigma galaxy depth and exposure time at each pixel.

    Parameters
    ----------
    pixels: `np.ndarray`
       NEST pixels to sample.
    maglim: `dict` [`str`, `healsparse.HealSparseMap`]
       5 sigma point source depth per band.
    exptime: `dict` [`str`, `healsparse.HealSparseMap`]
       Summed exposure time per band.
    bands: `list` [`str`]
       Bands to sample.

    Returns
    -------
    limmag: `dict` [`str`, `np.ndarray`]
       10 sigma galaxy depth per band.
    exp: `dict` [`str`, `np.ndarray`]
       Exposure time per band.
    has_depth: `np.ndarray`
       Boolean, True where every band has both.
    '''
    ## getting depth and exptime per band
    limmag, exp = {}, {}
    for b in bands:
        m = maglim[b].get_values_pix(pixels)
        t = exptime[b].get_values_pix(pixels)
        limmag[b] = np.where(m > 0, m, np.nan) - common.SIG_OFFSET
        exp[b] = np.where(t > 0, t, np.nan)

    ## checking every band has both
    has_depth = np.ones(len(pixels), dtype=bool)
    for b in bands:
        has_depth &= np.isfinite(limmag[b]) & np.isfinite(exp[b])

    return limmag, exp, has_depth


def depth_product(limmag, exp, bands, depth_band):
    '''
    Pick which depth goes in the file.

    'min' is the shallowest band (a color is only as good as its worse half),
    anything else is that one band.

    Parameters
    ----------
    limmag: `dict` [`str`, `np.ndarray`]
       10 sigma galaxy depth per band.
    exp: `dict` [`str`, `np.ndarray`]
       Exposure time per band.
    bands: `list` [`str`]
       Bands in the footprint.
    depth_band: `str`
       'min' or a band name.

    Returns
    -------
    limmag: `np.ndarray`
       Chosen depth per pixel.
    exp: `np.ndarray`
       Chosen exposure time per pixel.
    '''
    ## picking the shallowest band or the one asked for
    if depth_band == 'min':
        return (np.min([limmag[b] for b in bands], axis=0),
                np.min([exp[b] for b in bands], axis=0))
    return limmag[depth_band], exp[depth_band]


'''
writing
'''
def check_outputs(args, bands, overwrite):
    '''
    Refuse to overwrite existing maps before any work is done.

    Same guard the catalog build uses -- the old object-based maps have exactly
    these names, so a rerun would quietly replace them.

    Parameters
    ----------
    args: `argparse.Namespace`
       Parsed arguments.
    bands: `list` [`str`]
       Bands in the footprint.
    overwrite: `bool`
       Allow replacing the files.
    '''
    ## getting the mask paths
    paths = [common.bitpacked_file(args.out_dir, args.bands, args.vis_min, args.nside_sparse),
             common.mask_file(args.out_dir, args.bands, args.vis_min, args.nside)]
    
    ## adding a depth path per depth band
    for depth_band in args.depth_bands.split(','):
        paths.append(common.depth_file(args.out_dir, args.bands, args.vis_min, args.nside,
                                       depth_band, dered=args.dust_coeffs is not None))

    ## stopping if any already exist
    existing = [p for p in paths if os.path.exists(p)]
    if len(existing) > 0 and not overwrite:
        sys.exit(f'{len(existing)} map files already exist, e.g. {existing[0]}. '
                 f'move them or pass --overwrite')


def write_mask(pixels, fracgood, nside, nside_cover, out_path):
    '''
    Write the mask redMaPPer reads: fracgood per pixel, hp.UNSEEN outside.

    Parameters
    ----------
    pixels: `np.ndarray`
       NEST pixels in the footprint.
    fracgood: `np.ndarray`
       Covered fraction per pixel.
    nside: `int`
       Map nside.
    nside_cover: `int`
       Coverage nside.
    out_path: `str`
       Where to write.
    '''
    ## making and writing the mask
    mask = hsp.HealSparseMap.make_empty(nside_cover, nside, np.float32, sentinel=hp.UNSEEN)
    mask[pixels] = fracgood.astype(np.float32)
    mask.write(out_path, clobber=True)
    print(f'wrote geometry mask to {out_path}')
    print(f'  {np.sum(fracgood) * common.pix_area(nside):.1f} deg^2 (fracgood weighted), '
          f'median fracgood {np.median(fracgood):.2f}')


def write_depth(pixels, limmag, exp, fracgood, nside, nside_cover, bands, nsig, zp, out_path):
    '''
    Write the depth map in redMaPPer's format, same header keywords as the how-to.

    Parameters
    ----------
    pixels: `np.ndarray`
       NEST pixels in the footprint.
    limmag: `np.ndarray`
       10 sigma galaxy depth per pixel.
    exp: `np.ndarray`
       Exposure time per pixel.
    fracgood: `np.ndarray`
       Covered fraction per pixel.
    nside: `int`
       Map nside.
    nside_cover: `int`
       Coverage nside.
    bands: `str`
       Band string for the comment card.
    nsig: `float`
       S/N that limmag corresponds to.
    zp: `float`
       redMaPPer reference zeropoint.
    out_path: `str`
       Where to write.
    '''
    ## filling in the depth values
    values = np.zeros(len(pixels), dtype=depth_dtype)
    values['exptime'] = exp
    values['limmag'] = limmag
    values['m50'] = limmag          ## placeholder, same as cosmodc2
    values['fracgood'] = fracgood

    ## making and writing the map
    depth = hsp.HealSparseMap.make_empty(nside_coverage=nside_cover, nside_sparse=nside,
                                         dtype=depth_dtype, primary='limmag',
                                         sentinel=hp.UNSEEN)
    depth[pixels] = values
    depth.write(out_path, clobber=True)

    ## adding the header keywords
    with fitsio.FITS(out_path, 'rw') as ff:
        hdr = fitsio.FITSHDR()
        hdr['ZP'] = zp
        hdr['NSIG'] = nsig
        hdr['NBAND'] = 1  ## not used, but required in formatting
        hdr['W'] = 0.0
        hdr['EFF'] = 1.0
        hdr['COMMENT'] = f'dp2 depth, {bands}, {nsig} sigma galaxy'
        ff[0].write_keys(hdr)

    print(f'wrote depth map to {out_path}')
    print(f'  median limmag {np.median(limmag):.2f}, '
          f'10:90 percentile {np.percentile(limmag, 10):.2f}:{np.percentile(limmag, 90):.2f}')


def parse_args():
    '''
    Command line arguments.

    Returns
    -------
    args: `argparse.Namespace`
       Parsed arguments.
    '''
    p = argparse.ArgumentParser(description='Visit-based geometry mask + depth map for redmapper on DP2 (RSP).')

    p.add_argument('--tracts', default=None, help='comma separated, e.g. 10188,8228')
    p.add_argument('--tract-list', default=None, help='text file, one tract per line (default: all tracts)')
    p.add_argument('--bands', default='riz', help='bands required in the footprint (default riz)')
    p.add_argument('--vis-min', type=int, default=3,
                   help='visits required in every band over the 2x2 block (default 3, shear group)')
    p.add_argument('--weight-cut', type=float, default=WEIGHT_CUT,
                   help='unmasked_fraction a visit must exceed to count (default 0.0, shear group)')
    p.add_argument('--no-block', action='store_true',
                   help='test the cell alone instead of its 2x2 block (NOT the shear group rule)')
    p.add_argument('--no-tract-clip', action='store_true',
                   help='keep cells in the tract overlap borders instead of only owned cells')

    p.add_argument('--nside', type=int, default=4096, help='mask + depth NSIDE (default 4096, what DES used)')
    p.add_argument('--nside-cover', type=int, default=32, help='coverage NSIDE (default 32, shear group)')
    p.add_argument('--nside-sparse', type=int, default=131072,
                   help='cell-resolution NSIDE (default 131072, shear group)')

    p.add_argument('--ebv-max', type=float, default=0.2, help='max SFD E(B-V) (default 0.2, DES)')
    p.add_argument('--b-min', type=float, default=20.0, help='min |galactic b| in degrees (default 20)')
    p.add_argument('--fracgood-min', type=float, default=0.0,
                   help='drop pixels below this fracgood (default 0, keep every pixel the footprint touches)')

    p.add_argument('--depth-bands', default='min',
                   help="comma separated depth maps to write: 'min' (shallowest) and/or band names (default min)")
    p.add_argument('--dust-coeffs', default=None,
                   help="A_band/E(B-V) per band, e.g. 'r=..,i=..,z=..'. same values as the catalog build (no default)")
    p.add_argument('--zp', type=float, default=22.5, help='redmapper reference zeropoint')
    p.add_argument('--nsig', type=float, default=10.0, help='S/N at limmag')
    p.add_argument('--nside-config', type=int, default=2, help='NSIDE for the RING hpix list printed at the end')

    p.add_argument('--n-procs', type=int, default=1,
                   help='worker processes for stage 1, each with its own butler (default 1)')
    p.add_argument('--skip-cells', action='store_true', help='just rebuild the maps from saved cell counts')
    p.add_argument('--clobber', action='store_true', help='recount tracts that already have a cell count file')
    p.add_argument('--overwrite', action='store_true', help='replace existing mask and depth files')
    p.add_argument('--out-dir', default=os.path.expanduser('~/dp2_redmapper/maps'))

    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    bands = list(args.bands)
    pa = common.pix_area(args.nside)

    check_outputs(args, bands, args.overwrite)

    butler, skymap = common.open_butler()

    ## stage 1: counting visits per cell
    if not args.skip_cells:
        tracts = common.get_tract_list(args.tracts, args.tract_list, butler)
        print(f'reading provenance for {len(tracts)} tracts, bands {args.bands}, '
              f'unmasked_fraction > {args.weight_cut}, detector {BAD_DETECTOR} dropped')
        run_cell_counts(butler, skymap, tracts, bands, args.weight_cut, args.out_dir, args.clobber,
                        n_procs=args.n_procs)

    files = load_cell_counts(args.out_dir, bands, args.weight_cut)

    ## stage 2: rendering the good cells
    block = 'the cell alone' if args.no_block else 'the 2x2 block'
    print(f'\nrendering cells with >= {args.vis_min} visits in every band over {block}')
    geometry = build_geometry(skymap, files, bands, args.vis_min, args.nside_sparse, args.nside_cover,
                              use_block=not args.no_block, clip_tract=not args.no_tract_clip)
    print(f'\n{"cells with enough visits":<30} {geometry.get_valid_area(degrees=True):8.1f} deg^2')

    ## degrading to fracgood
    frac = fracgood_map(geometry, args.nside)
    pixels = frac.valid_pixels
    fracgood = frac.get_values_pix(pixels)
    keep = fracgood > args.fracgood_min
    pixels, fracgood = pixels[keep], fracgood[keep]

    ## cutting on dust and galactic latitude
    fg_keep, ebv = foreground_cut(pixels, args.nside, args.ebv_max, args.b_min)
    pixels, fracgood, ebv = pixels[fg_keep], fracgood[fg_keep], ebv[fg_keep]
    print(f'{f"+ E(B-V) < {args.ebv_max}, |b| > {args.b_min}":<30} '
          f'{np.sum(fracgood) * pa:8.1f} deg^2')

    ## loading depth and dropping pixels without it
    print('\nloading survey property maps')
    maglim, exptime = load_depth_maps(butler, bands, args.nside)
    limmag, exp, has_depth = depth_values(pixels, maglim, exptime, bands)
    pixels, fracgood, ebv = pixels[has_depth], fracgood[has_depth], ebv[has_depth]
    limmag = {b: v[has_depth] for b, v in limmag.items()}
    exp = {b: v[has_depth] for b, v in exp.items()}
    final_area = np.sum(fracgood) * pa
    print(f'{"+ depth in every band":<30} {final_area:8.1f} deg^2\n')

    ## dust correcting the depth if asked
    dust_coeffs = common.parse_dust_coeffs(args.dust_coeffs, bands)
    if dust_coeffs is not None:
        limmag = {b: limmag[b] - dust_coeffs[b] * ebv for b in bands}
        print(f'dust corrected depth with {dust_coeffs}')
    else:
        print('depth is NOT dust corrected')

    ## dropping the cut pixels from the fine map
    dropped = np.setdiff1d(frac.valid_pixels, pixels)
    if len(dropped) > 0:
        clear_pixels(geometry, dropped, args.nside, args.nside_sparse)
        print(f'removed {len(dropped)} pixels from the cell-resolution map '
              f'({geometry.get_valid_area(degrees=True):.1f} deg^2 left)')

    ## writing the cell-resolution map
    fine_path = common.bitpacked_file(args.out_dir, args.bands, args.vis_min, args.nside_sparse)
    geometry.write(fine_path, clobber=True)
    print(f'wrote cell-resolution geometry to {fine_path}')

    ## writing the mask
    write_mask(pixels, fracgood, args.nside, args.nside_cover,
               common.mask_file(args.out_dir, args.bands, args.vis_min, args.nside))

    ## writing each depth map
    for depth_band in args.depth_bands.split(','):
        if depth_band != 'min' and depth_band not in bands:
            sys.exit(f'depth band {depth_band} is not in {args.bands}')
        product_limmag, product_exp = depth_product(limmag, exp, bands, depth_band)
        print(f'\ndepth: {depth_band}')
        write_depth(pixels, product_limmag, product_exp, fracgood, args.nside, args.nside_cover,
                    args.bands if depth_band == 'min' else depth_band, args.nsig, args.zp,
                    common.depth_file(args.out_dir, args.bands, args.vis_min, args.nside, depth_band,
                                      dered=dust_coeffs is not None))

    ## getting the hpix list for the config
    hpix = common.ring_hpix(pixels, args.nside, args.nside_config)
    print(f'\nFinal area: {final_area:.1f} deg^2 (fracgood weighted)')
    print(f'config hpix (RING, nside {args.nside_config}): {hpix.tolist()}')
    print('Done!')


if __name__ == '__main__':
    main()
