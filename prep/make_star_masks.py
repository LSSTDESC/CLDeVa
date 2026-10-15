'''
make_star_masks.py

bright star masks with only the brighter Gaia stars (G < 12, 15, 16)
needs the LSST stack (the Gaia stars come from the_monster refcat in the butler)

python make_star_masks.py --out-dir $SCRATCH/dp2_redmapper/v2/star_masks
'''

import sys; import os
import argparse
import multiprocessing as mp

import numpy as np
import healpy as hp
import healsparse as hsp
import lsst.sphgeom as sg
from lsst.daf.butler import Butler
from astropy.coordinates import SkyCoord

import redmapper_dp2_common as common

'''
Initial definitions
'''
## the registered star mask, used for which part of the sky (coverage) and to check against
SHEAR_MASK = '/global/cfs/cdirs/lsst/utilities/data-registry/lsst_desc_working/user/pladuca/.gen_paths/bright_star_mask_1.2.0/DP2_gaia_mask_131k.hs'
REFCAT = 'the_monster_20250219'
REFCAT_COLLECTION = 'refcats/DM-49042/the_monster_20250219'

STARS = [12, 15, 16] ## G cuts

## circle radius vs G, median radius of the circles in the shear team's mask around stars in each 1 magnitude bin
G_TABLE = [5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5]
R_TABLE = [348.9, 252.0, 163.3, 124.7, 95.0, 61.6, 44.5, 32.1, 23.2, 16.8, 13.5, 10.9]   ## arcsec


def start_worker():
    global butler
    butler = Butler('dp2', collections=[REFCAT_COLLECTION])


def read_shard(shard):
    ## Gaia stars brighter than the faintest cut in one refcat shard
    cat = butler.get(REFCAT, htm7=int(shard)) 
    flux = cat['phot_g_mean_flux']
    ok = flux > 0
    g = -2.5 * np.log10(flux[ok]) + 31.4 ## nJy -> AB
    bright = g < max(STARS)
    return np.degrees(cat['coord_ra'][ok][bright]), np.degrees(cat['coord_dec'][ok][bright]), g[bright]


def get_stars(cov, n_procs):
    ## every Gaia star brighter than the faintest cut in these pixels
    htm = sg.HtmPixelization(7) #required for butler
    sub = (cov[:, None] * (1024 // 32) ** 2 + np.arange((1024 // 32) ** 2)).ravel()
    ra, dec = hp.pix2ang(1024, sub, nest=True, lonlat=True)
    shards = np.unique([htm.index(sg.UnitVector3d(sg.LonLat.fromDegrees(r, d))) for r, d in zip(ra, dec)])
    print(f'{len(shards)} refcat shards', flush=True)

    ## reading the shards is all butler latency, so a pool of workers
    ra_all, dec_all, g_all = [], [], []
    with mp.get_context('fork').Pool(n_procs, initializer=start_worker) as pool:
        for n, (r, d, g) in enumerate(pool.imap_unordered(read_shard, shards, chunksize=20)):
            ra_all.append(r)
            dec_all.append(d)
            g_all.append(g)
            if (n + 1) % 2000 == 0:
                print(f'  [{n+1}/{len(shards)}]', flush=True)
    ra, dec, g = np.concatenate(ra_all), np.concatenate(dec_all), np.concatenate(g_all)

    ## shards overlap the edges, only keep stars in our coverage pixels
    inside = np.isin(hp.ang2pix(32, ra, dec, lonlat=True, nest=True), cov)
    return ra[inside], dec[inside], g[inside]


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--out-dir', required=True)
    p.add_argument('--n-procs', type=int, default=8)
    args = p.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    ## coverage pixels of the shear team mask for comp
    cov = np.where(hsp.HealSparseCoverage.read(SHEAR_MASK).coverage_mask)[0]
    ra, dec = hp.pix2ang(32, cov, nest=True, lonlat=True)
    cov = cov[SkyCoord(ra, dec, unit='deg').galactic.b.deg < -18]
    ra, dec, g = get_stars(cov, args.n_procs)
    radius = np.interp(g, G_TABLE, R_TABLE) / 3600.0   ## degrees
    print(f'{len(g)} stars with G < {max(STARS)} in {len(cov)} coverage pixels')

    ## building our mask
    mask = hsp.HealSparseMap.make_empty(32, 131072, bool, bit_packed=True)
    g_lo = -99
    for star in STARS:
        these = np.where((g >= g_lo) & (g < star))[0]
        for start in range(0, len(these), 2000):
            pix = [hsp.geom.Circle(ra=ra[k], dec=dec[k], radius=radius[k], value=True).get_pixels(nside=131072)
                   for k in these[start:start + 2000]]
            mask.update_values_pix(np.unique(np.concatenate(pix)), True)
        path = common.star_mask_file(args.out_dir, star)
        mask.write(path, clobber=True)
        print(f'G < {star}: {len(these)} more stars ({np.sum(g < star)} total), {mask.get_valid_area(degrees=True):.1f} deg^2 masked -> {path}', flush=True)
        g_lo = star

    ## check against shear mask 
    shear = hsp.HealSparseMap.read(SHEAR_MASK, pixels=cov[::20])
    ours = hsp.HealSparseMap.read(common.star_mask_file(args.out_dir, max(STARS)), pixels=cov[::20])
    pix = np.concatenate([ours.valid_pixels_single_covpix(cp) for cp in cov[::20]])
    print(f'check on {len(cov[::20])} coverage pixels: {shear.get_values_pix(pix).mean() * 100:.1f}% of our G < {max(STARS)} '
          f'masked pixels are masked in the shear team mask too')
    print('Done!')
