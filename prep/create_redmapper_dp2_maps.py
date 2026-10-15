'''
create_redmapper_dp2_maps.py

geometry mask + depth maps for redmapper on DP2
needs the LSST stack (butler); run with nersc/env_butler.sh

step 1: count visits per coadd cell for every tract (slow, saved per tract so it can be restarted)
step 2: cells with enough visits -> healsparse map -> fracgood at nside 4096 -> dust/b/depth cuts (b < -20 only)

depth is dust corrected (0.86 x CSFD). writes one version without a star mask, then one for each
Gaia star mask (G < 12, 15, 16, made with make_star_masks.py)

python create_redmapper_dp2_maps.py --bands griz --vis-min 2 --n-procs 8 --out-dir $SCRATCH/dp2_redmapper/v2/griz_vis2/maps --star-dir $SCRATCH/dp2_redmapper/v2/star_masks
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
from dustmaps.csfd import CSFDQuery

import lsst.geom as geom
from lsst.skymap import Index2D

import redmapper_dp2_common as common

'''
Initial definitions
'''
NSIDE = 4096          ## mask + depth, same as DES
NSIDE_COVER = 32      ## shear group's number
NSIDE_SPARSE = 131072 ## cell resolution, shear group

EBV_MAX = 0.2   ## SFD, same as DES
B_MIN = 20.0    ## galactic latitude cut, degrees. keeping only b < -20

STARS = [12, 15, 16]   ## Gaia G cuts of the star masks. each one masks more than the one before

## shear group's rule: a visit counts if it covers any of the cell (unmasked_fraction > 0), and detector 122 is dropped
WEIGHT_CUT = 0.0
BAD_DETECTOR = 122

## cells: a patch is 22x22 cells, the not overlapping ones are 1-20
N_CELL = 22

MAGLIM_MAP = 'deepCoadd_psf_maglim_consolidated_map_weighted_mean'
EXPTIME_MAP = 'deepCoadd_exposure_time_consolidated_map_sum'
SIG_OFFSET = 2.5 * np.log10(10.0 / 5.0) ## 5 sigma point source -> 10 sigma galaxy

ZP = 22.5
NSIG = 10.0

depth_dtype = np.dtype([('exptime', 'f4'), ('limmag', 'f4'), ('m50', 'f4'), ('fracgood', 'f4')])


'''
step 1: visit counts per cell
'''
def start_worker():
    ## setting up workers for parallelization
    ## need to make butler, skymap for this
    global butler, skymap
    butler, skymap = common.open_butler()


def count_tract(tract):
    ## visit counts per cell for every patch and band of one tract, saved to a npz
    try:
        tract_info = skymap[tract]
        n = tract_info.num_patches

        ## gets all coadds in tract
        refs = butler.query_datasets('deep_coadd', skymap=common.SKYMAP, tract=tract)
        refs = {(int(r.dataId['patch']), str(r.dataId['band'])): r for r in refs}

        patches = []
        counts = []
        for patch in range(int(n[0]) * int(n[1])):
            per_band = np.zeros((len(bands), N_CELL, N_CELL), dtype=np.int16)
            found = False
            for k, band in enumerate(bands):
                if (patch, band) not in refs:
                    continue   ## no coadd in this band -> 0 visits
                prov = butler.get(refs[(patch, band)].makeComponentRef('provenance')).contributions

                ## count the good visits in each cell
                good = (np.asarray(prov['unmasked_fraction']) > WEIGHT_CUT) & (np.asarray(prov['detector']) != BAD_DETECTOR)
                cell = np.asarray(prov['cell_i'])[good] * N_CELL + np.asarray(prov['cell_j'])[good]
                per_band[k] = np.bincount(cell, minlength=N_CELL * N_CELL).reshape(N_CELL, N_CELL)
                found = True
            if found:
                patches.append(patch)
                counts.append(per_band)
    except Exception as err:
        ## in case there's an error
        return f'tract {tract} FAILED: {err}'

    counts = np.array(counts, dtype=np.int16).reshape(-1, len(bands), N_CELL, N_CELL)
    np.savez(common.cell_counts_file(out_dir, tract, ''.join(bands)), patches=np.array(patches, dtype=np.int32), counts=counts)
    return f'tract {tract}: {len(patches)} patches'


'''
step 2: good cells -> healsparse map
'''
def cell_box(patch_info, i, j):
    return patch_info.getCellInfo(Index2D(x=i, y=j)).inner_bbox


def owned_cells(tract_info, tract, patches):
    ## getting which of the 20x20 inner cells of each 
    #patch belong to this tract (cell centers -> skymap)
    centers = np.arange(20) + 0.5
    x = np.zeros((len(patches), 20, 20))
    y = np.zeros((len(patches), 20, 20))
    for p, patch in enumerate(patches):
        patch_info = tract_info[int(patch)]
        box = cell_box(patch_info, 1, 1)
        origin = patch_info.getInnerBBox().getMin()
        x[p] = origin.getX() + centers[:, None] * box.getWidth()
        y[p] = origin.getY() + centers[None, :] * box.getHeight()

    ra, dec = tract_info.wcs.pixelToSkyArray(x.ravel(), y.ravel(), degrees=True)
    return common.owned_by_tract(skymap, tract, ra, dec).reshape(x.shape)


def tract_polygons(path, vis_min):
    ## good cells of one tract -> sky rectangles
    tract = int(os.path.basename(path).split('_')[1].replace('tract', ''))
    data = np.load(path)
    patches, counts = data['patches'], data['counts']
    if len(patches) == 0:
        return []
    tract_info = skymap[tract]

    ## shear group's 2x2 block: a cell passes if it and its i-1, j-1 neighbours all have enough visits
    ## in their worst band. block[:, i-1, j-1] is for cell (i, j)
    worst = counts.min(axis=1)
    block = np.minimum.reduce([worst[:, 0:20, 0:20], worst[:, 1:21, 0:20],
                               worst[:, 0:20, 1:21], worst[:, 1:21, 1:21]])
    keep = (block >= vis_min) & owned_cells(tract_info, tract, patches)

    polygons = []
    for p, patch in enumerate(patches):
        patch_info = tract_info[int(patch)]
        good = np.zeros((N_CELL, N_CELL), dtype=bool)
        good[1:21, 1:21] = keep[p]

        ## good cells next to each other in a row -> one rectangle, way fewer polygons
        for i in range(1, 21):
            j = 1
            while j < 21:
                if not good[i, j]:
                    j += 1
                    continue
                j_start = j
                while j < 21 and good[i, j]:
                    j += 1
                box = geom.Box2D(cell_box(patch_info, i, j_start))
                box.include(geom.Box2D(cell_box(patch_info, i, j - 1)))
                corners = [tract_info.wcs.pixelToSky(c) for c in box.getCorners()]
                ra = [c.getRa().asDegrees() for c in corners]
                dec = [c.getDec().asDegrees() for c in corners]
                polygons.append(hsp.geom.Polygon(ra=ra, dec=dec, value=True))
    return polygons


def clear_stars(cells, star_path):
    ## set every cell inside the star mask to False, one coverage pixel at a time
    cov = np.where(cells.coverage_mask)[0]
    star = hsp.HealSparseMap.read(star_path, pixels=cov)
    for cp in cov:
        pix = cells.valid_pixels_single_covpix(cp)
        cells.update_values_pix(pix[star.get_values_pix(pix)], False)


def clear_pixels(cells, pixels):
    ## set every cell-resolution pixel under these nside 4096 pixels to False, in chunks (1024 each)
    ratio = (NSIDE_SPARSE // NSIDE) ** 2
    for start in range(0, len(pixels), 5000):
        parents = pixels[start:start + 5000].astype(np.int64)
        children = (parents[:, None] * ratio + np.arange(ratio, dtype=np.int64)).ravel()
        cells.update_values_pix(children, False)


'''
writing
'''
def write_maps(name, pixels, fracgood, limmag, exptime, cells):
    ## cell-resolution map
    cells.write(common.cells_file(out_dir, name), clobber=True)

    ## fracgood mask, UNSEEN outside
    mask = hsp.HealSparseMap.make_empty(NSIDE_COVER, NSIDE, np.float32, sentinel=hp.UNSEEN)
    mask[pixels] = fracgood.astype(np.float32)
    mask.write(common.geometry_file(out_dir, name), clobber=True)

    ## depth maps, shallowest band ('min') and each band
    for band in ['min'] + bands:
        values = np.zeros(len(pixels), dtype=depth_dtype)
        if band == 'min':
            values['limmag'] = np.min([limmag[b] for b in bands], axis=0)
            values['exptime'] = np.min([exptime[b] for b in bands], axis=0)
        else:
            values['limmag'] = limmag[band]
            values['exptime'] = exptime[band]
        values['m50'] = values['limmag']   ## placeholder, same as cosmoDC2
        values['fracgood'] = fracgood

        path = common.depth_file(out_dir, name, band)
        depth = hsp.HealSparseMap.make_empty(nside_coverage=NSIDE_COVER, nside_sparse=NSIDE, dtype=depth_dtype,
                                             primary='limmag', sentinel=hp.UNSEEN)
        depth[pixels] = values
        depth.write(path, clobber=True)

        ## header keywords redmapper wants
        hdr = fitsio.FITSHDR()
        hdr['ZP'] = ZP
        hdr['NSIG'] = NSIG
        hdr['NBAND'] = 1 ## not used but has to be there
        hdr['W'] = 0.0
        hdr['EFF'] = 1.0
        hdr['COMMENT'] = f'dp2 depth, {"".join(bands) if band == "min" else band}, {NSIG} sigma galaxy, dust corrected'
        with fitsio.FITS(path, 'rw') as f:
            f[0].write_keys(hdr)

    ## the hpix list for the redmapper config (RING, nside 2)
    ra, dec = hp.pix2ang(NSIDE, pixels, lonlat=True, nest=True)
    hpix = np.unique(hp.ang2pix(2, ra, dec, lonlat=True))
    print(f'{name}: {np.sum(fracgood) * hp.nside2pixarea(NSIDE, degrees=True):.1f} deg^2, '
          f'median depth {np.median(np.min([limmag[b] for b in bands], axis=0)):.2f}, config hpix {hpix.tolist()}')


'''
main
'''
p = argparse.ArgumentParser()
p.add_argument('--bands', default='griz')
p.add_argument('--vis-min', type=int, default=2, help='visits needed in every band over the 2x2 block')
p.add_argument('--tracts', default=None, help='e.g. 10188,8228 for testing (default all)')
p.add_argument('--n-procs', type=int, default=4)
p.add_argument('--out-dir', required=True)
p.add_argument('--star-dir', required=True, help='where make_star_masks.py put the star masks')
args = p.parse_args()

bands = list(args.bands)
out_dir = args.out_dir

if __name__ == '__main__':
    os.makedirs(os.path.join(out_dir, 'cell_counts'), exist_ok=True)
    butler, skymap = common.open_butler()

    ## step 1, skipping tracts that already have a file (delete them to recount)
    tracts = common.get_tracts(butler, args.tracts)
    todo = [t for t in tracts if not os.path.exists(common.cell_counts_file(out_dir, t, args.bands))]
    print(f'counting visits for {len(todo)} tracts ({len(tracts) - len(todo)} already done), bands {args.bands}')
    with mp.get_context('fork').Pool(args.n_procs, initializer=start_worker) as pool:
        for n, message in enumerate(pool.imap_unordered(count_tract, todo)):
            print(f'[{n+1}/{len(todo)}] {message}', flush=True)

    ## step 2, every good cell into one map
    files = sorted(glob.glob(common.cell_counts_file(out_dir, '*', args.bands)))
    print(f'\nmaking the cell map from {len(files)} tracts, >= {args.vis_min} visits in {args.bands}')
    cells = hsp.HealSparseMap.make_empty(NSIDE_COVER, NSIDE_SPARSE, bool, bit_packed=True)
    for n, path in enumerate(files):
        polygons = tract_polygons(path, args.vis_min)
        if len(polygons) > 0:
            pix = np.unique(np.concatenate([poly.get_pixels(nside=NSIDE_SPARSE) for poly in polygons]))
            cells.update_values_pix(pix, True)
        if (n + 1) % 200 == 0:
            print(f'  [{n+1}/{len(files)}]', flush=True)
    print(f'cells with enough visits: {cells.get_valid_area(degrees=True):.1f} deg^2')

    ## fraction of each nside 4096 pixel that's covered
    version = [int(v) for v in hsp.__version__.split('.')[:2]]
    assert version >= [1, 15], 'need healsparse >= 1.15, see nersc/env_butler.sh'
    frac = cells.fracdet_map(NSIDE)
    all_pixels = frac.valid_pixels
    fracgood = frac.get_values_pix(all_pixels)
    pixels = all_pixels[fracgood > 0]
    fracgood = fracgood[fracgood > 0]

    ## dust and galactic latitude at the pixel centers
    ra, dec = hp.pix2ang(NSIDE, pixels, lonlat=True, nest=True)
    coords = SkyCoord(ra * u.deg, dec * u.deg)
    ebv = SFDQuery()(coords)
    gal_b = coords.galactic.b.deg

    ## depth in every band, then degraded to nside 4096 
    limmag = {}
    exptime = {}
    has_depth = np.ones(len(pixels), dtype=bool)
    for b in bands:
        m = butler.get(MAGLIM_MAP, band=b, skymap=common.SKYMAP, parameters={'degrade_nside': NSIDE}).get_values_pix(pixels)
        t = butler.get(EXPTIME_MAP, band=b, skymap=common.SKYMAP, parameters={'degrade_nside': NSIDE}).get_values_pix(pixels)
        limmag[b] = np.where(m > 0, m, np.nan) - SIG_OFFSET
        exptime[b] = np.where(t > 0, t, np.nan)
        has_depth &= np.isfinite(limmag[b]) & np.isfinite(exptime[b])

    ## dust correcting the depth, E(B-V) = 0.86 x CSFD
    csfd = CSFDQuery(map_fname=common.CSFD_DIR + 'csfd_ebv.fits', mask_fname=common.CSFD_DIR + 'mask.fits')
    ebv_csfd = 0.86 * csfd(coords)
    for b in bands:
        limmag[b] = limmag[b] - common.DUST_COEFFS[b] * ebv_csfd

    ## footprint cuts: SFD E(B-V) < 0.2 (like DES), b < -20, depth in every band
    keep = (ebv < EBV_MAX) & (gal_b < -B_MIN) & has_depth
    pixels = pixels[keep]
    fracgood = fracgood[keep]
    limmag = {b: limmag[b][keep] for b in bands}
    exptime = {b: exptime[b][keep] for b in bands}
    clear_pixels(cells, np.setdiff1d(all_pixels, pixels))
    write_maps(common.run_name(args.bands, args.vis_min), pixels, fracgood, limmag, exptime, cells)

    ## then each star mask, taking more cells out of the same map every time
    for star in STARS:
        clear_stars(cells, common.star_mask_file(args.star_dir, star))
        frac = cells.fracdet_map(NSIDE).get_values_pix(pixels)
        sel = frac > 0   ## pixels that are completely masked drop out
        write_maps(common.run_name(args.bands, args.vis_min, star), pixels[sel], frac[sel],
                   {b: limmag[b][sel] for b in bands}, {b: exptime[b][sel] for b in bands}, cells)
    print('Done!')
