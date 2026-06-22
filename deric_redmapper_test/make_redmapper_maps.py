'''
make_redmapper_maps.py 

Builds the two healsparse map for redmapper from cosmoDC2

Survey geometry mask as cosmodc2_geometry_mask.hs
depth map as cosmodc2_depth_map.hs
'''

import sys; import os

os.environ["GCR_CONFIG_SOURCE"] = "files"
import GCRCatalogs

import numpy as np
import healpy as hp
import healsparse as hsp
import fitsio
import argparse

depth_dtype = np.dtype([("exptime", "f4"), ("limmag", "f4"),
    ("m50", "f4"), ("fracgood", "f4")])

def load_galaxy_positions(cat, zmin, zmax, mag_cut, ref_band="z"):
    '''
    return ra,dec arrays for all galaxies passing the mag cut
    '''
    quantities = ["ra", "dec"]
    filters = [f"mag_{ref_band}_lsst <= {mag_cut}", f"redshift_true >= {zmin}", f"redshift_true <= {zmax}"]

    data = cat.get_quantities(quantities, filters=filters)
    ra = np.array(data["ra"], dtype=np.float64)
    dec = np.array(data["dec"], dtype=np.float64)
    print(f"loaded {len(ra):,} galaxies") 

    return ra, dec

def radec_to_healpix(ra, dec, nside, nest=True):
    '''
    convert ra and dec arrays to healpix pixel indeces
    nside tells resolution of pixels. Larger nside means finer/smaller pixels. 
    - sparse gives actual resolution of the map, so more for position guiding
    - cover divides the sky into coarse "coverage" pixels and only stores the fine pixels that actually have data

    nest is nesting pixeling order (ring or nested)
    - Ring ordering numbers pixels in horizontal rings around the sphere
    - Nested ordering groups pixels into a hierarchical tree structure where nearby pixels on the sky have nearby indices
    '''
    theta = np.radians(90.0-dec) 
    phi = np.radians(ra) 
    return hp.ang2pix(nside, theta, phi, nest=nest)

def build_geometry_mask(ra, dec, nside_sparse, nside_cover, out_path): 
    '''
    mark all healpix pixels populate by galaxies as "valid"
    '''
    pixels = radec_to_healpix(ra, dec, nside=nside_sparse, nest=True)
    unique_pixels = np.unique(pixels) 

    sparse_map = hsp.HealSparseMap.make_empty(nside_coverage=nside_cover, nside_sparse=nside_sparse, dtype=np.uint8, sentinel=0)
    sparse_map[unique_pixels] = np.uint8(1)

    sparse_map.write(out_path, clobber=True)
    print(f"write geometry mask to {out_path}")

    #area estimate
    pix_area_deg2 = (180.0 / np.pi)**2 * 4*np.pi/(12*nside_sparse**2)
    area_est = len(unique_pixels)*pix_area_deg2
    print(f"est footprint area: {area_est:.1f} deg^2")
    print("cosmodc2 is 440 deg^2") 

    return sparse_map, unique_pixels

def build_depth_mask(unique_pixels, nside_sparse, nside_cover, limmag, exptime, nsig, zp, out_path): 
    '''
    build a uniform depth map over the footprint for rm as healsparse

    cosmodc2: all pixels get same limmag, exptime, m50
    dp2: replace with per-pixel values from Butler deepCoadd metadata
    '''
    n = len(unique_pixels) 
    depth_vals = np.zeros(n, dtype=depth_dtype) 
    depth_vals["exptime"] = np.float32(exptime) 
    depth_vals["limmag"] = np.float32(limmag)
    depth_vals["m50"] = np.float32(limmag) 
    depth_vals["fracgood"] = np.float32(1.0) 

    sparse_map = hsp.HealSparseMap.make_empty(nside_coverage=nside_cover, nside_sparse=nside_sparse, 
                                              dtype=depth_dtype, primary="limmag", sentinel=hp.UNSEEN)

    sparse_map[unique_pixels] = depth_vals

    sparse_map.write(out_path, clobber=True)

    # fits file to add required header keywords
    with fitsio.FITS(out_path, "rw") as ff: #"fits file"
        hdr = fitsio.FITSHDR()
        hdr["ZP"] = zp
        hdr["NSIG"] = nsig
        hdr["NBAND"] = 1 #not used, but required in formatting
        hdr["W"] = 0.0
        hdr["EFF"] = 1.0
        hdr["LIMMAG"] = limmag
        hdr["EXPTIME"] = exptime
        hdr["COMMENT"] = "cosmodc2 depth mask uniform depth" 

        ff[0].write_keys(hdr)

    print(f"limmag={limmag}, exptime={exptime}s, nsig={nsig}")
    print(f"write depth mask to {out_path}")

    return sparse_map

def parse_args():
    p = argparse.ArgumentParser(description="Build redMaPPer geometry + depth healsparse masks from cosmoDC2.")
    
    p.add_argument("--catalog", default="cosmoDC2_v1.1.4_image")
    p.add_argument("--nside-sparse", type=int,   default=32768, help="HealSparse fine NSIDE (default 32768 ~ 0.1 arcmin pixels)")
    p.add_argument("--nside-cover",  type=int,   default=32, help="HealSparse coverage NSIDE (default 32)")
    p.add_argument("--zmin", type=float, default=0.10)
    p.add_argument("--zmax", type=float, default=0.90)
    p.add_argument("--mag-cut", type=float, default=26.0)
    p.add_argument("--limmag", type=float, default=26.0, help="Uniform 5σ limiting magnitude for depth map")
    p.add_argument("--zp", type=float, default=22.5, help="Reference zeropoint (default 22.5 AB)")
    p.add_argument("--nsig", type=float, default=10.0, help="S/N at limmag (default 10)")
    p.add_argument("--exptime", type=float, default=3000.0, help="Nominal effective exposure time in seconds")
    p.add_argument("--out-dir", default="./")
    
    return p.parse_args()

def main():
    args=parse_args()
    os.makedirs(args.out_dir, exist_ok=True) 

    print(f"loading catalog: {args.catalog}") 
    cat = GCRCatalogs.load_catalog(args.catalog) 

    ra, dec = load_galaxy_positions(cat, args.zmin, args.zmax, args.mag_cut) 

    geom_path = os.path.join(args.out_dir, "cosmodc2_geometry_mask.hs") 
    depth_path = os.path.join(args.out_dir, "cosmodc2_depth_mask.hs") 

    _, unique_pixels = build_geometry_mask(ra, dec, args.nside_sparse, args.nside_cover, geom_path) 

    build_depth_mask(unique_pixels, args.nside_sparse, args.nside_cover, limmag=args.limmag, exptime=args.exptime, nsig=args.nsig, zp=args.zp, out_path=depth_path) 

    print(f"Geometry mask: {geom_path}")
    print(f"Depth mask: {depth_path}")
    print("Done!")

if __name__ == "__main__":
    main()