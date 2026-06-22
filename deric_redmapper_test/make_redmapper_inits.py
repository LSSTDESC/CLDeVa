'''
make_redmapper_inits.py

This will build the specz training catalog, the photometric galaxy catalog, and the initial red-sequence template:
    cosmodc2_specz.fits
    cosmodc2_photo.fits
    cosmodc2_redsequence_init.fit

Runs on the "desc_rm" conda environment I set up
'''

import sys; import os

os.environ["GCR_CONFIG_SOURCE"] = "files"
import GCRCatalogs

import numpy as np
import argparse
from astropy.io import fits
from astropy.table import Table
import redmapper

'''
Initial definitions
'''
# BANDS = ['u', 'g', 'r', 'i', 'z', 'y']
# n_bands = len(BANDS)
# reference_band = 'z'; ref_i = 4
BANDS = ['g', 'r', 'i', 'z', 'y']
n_bands = len(BANDS)
reference_band = 'z'; ref_i = 3

zeropoint = 22.5 #LSST standard
area = 440.0 #cosmoDC2 footprint

z_err = 0.0002
mag_err = 0.01 #placeholder
ebv_val = 0.0 #no dust?

dz_specz = 0.05 #specz redshift bin width
dz_rs = 0.02 #red sequence bin width

'''
soecz subsample for training
'''
def build_specz_catalog(cat, zmin, zmax, n_per_bin, seed, out_path): 
    '''
    Uses central galaxies from GCR catalog
    '''

    quants = ["ra", "dec", "redshift_true", "is_central", "halo_mass"]
    filters = ["is_central == True", "halo_mass > 5e13", f"redshift_true >= {zmin}", f"redshift_true <= {zmax}"]

    data = cat.get_quantities(quants, filters=filters)
    ra = np.array(data["ra"]); dec = np.array(data["dec"])
    z = np.array(data["redshift_true"])
    print(f"{len(z)} cent gals in z range")

    rng = np.random.default_rng(seed)
    bins = np.arange(zmin, zmax + dz_specz, dz_specz)
    keep = []

    # loop through
    for z_low, z_high in zip(bins[:-1], bins[1:]):
        i = np.where((z >= z_low) & (z < z_high))[0]
        if len(i) == 0:
            print("empty bin")
            continue

        sel = rng.choice(i, size=min(len(i), n_per_bin), replace=False)
        keep.append(sel)
        print(f"{len(sel)} in [{z_low}, {z_high}]")

    keep = np.concatenate(keep)

    ra_s, dec_s, z_s = ra[keep], dec[keep], z[keep]
    zerr = np.full(len(keep), z_err, dtype=np.float32)

    tab = Table()
    tab["ra"] = ra_s.astype(np.float64)
    tab["dec"] = dec_s.astype(np.float64)
    tab["z"] = z_s.astype(np.float32)
    tab["z_err"] = zerr

    tab.write(out_path, format='fits', overwrite=True)
    print(f"wrote file to {out_path}")
    
    return tab

'''
photometric catalog
'''

def build_photometric_catalog(cat, zmin, zmax, mag_cut, out_path):
    mag_cols = [f"mag_{b}_lsst" for b in BANDS]
    quants = ["ra", "dec", "redshift_true"] + mag_cols #want ztrue since its sim
    filters = [f"mag_{reference_band}_lsst <= {mag_cut}", f"redshift_true >= {zmin}", f"redshift_true <= {zmax}",]

    data = cat.get_quantities(quants, filters=filters) 
    n = len(data["ra"])

    print(f"Loaded {n:,} galaxies with {reference_band}-band leq {mag_cut}")

    ra = np.array(data["ra"], dtype=np.float64)
    dec = np.array(data["dec"], dtype=np.float64) 

    mags = np.zeros((n, n_bands), dtype=np.float32)
    mag_errs = np.zeros((n, n_bands), dtype=np.float32)

    for i, b in enumerate(BANDS): 
        m = np.array(data[f"mag_{b}_lsst"], dtype=np.float32)
        mags[:, i] = m

        # placeholder SNR error estimate
        mag_errs[:, i] = np.clip(0.01* 10**(0.4*(m-mag_cut+1.0)), 0.005, 0.5).astype(np.float32)

    refmag = mags[:, ref_i]
    refmag_err = mag_errs[:, ref_i]
    ebv = np.full(n, ebv_val, dtype=np.float32)

    gal_id = np.arange(n, dtype=np.int64) # unique integer ID per galaxy

    # match redMaPPer expectations for fits file
    '''
    cols = [fits.Column(name="RA", format="D",  array=ra),
        fits.Column(name="DEC", format="D",  array=dec),
        fits.Column(name="MAG", format=f"{n_bands}E", array=mags),
        fits.Column(name="MAG_ERR", format=f"{n_bands}E", array=mag_errs),
        fits.Column(name="REFMAG", format="E", array=refmag),
        fits.Column(name="REFMAG_ERR", format="E", array=refmag_err),
        fits.Column(name="EBV", format="E",  array=ebv)]

    hdr = fits.Header()
    hdr["ZP"] = zeropoint
    hdr["NMAG"] = n_bands
    hdr["MODE"] = "LSST"
    hdr["LIM_REF"] = mag_cut
    hdr["REF_IND"] = ref_i
    hdr["AREA"] = area
    hdr["G_IND"] = 1
    hdr["R_IND"] = 2
    hdr["I_IND"] = 3
    hdr["Z_IND"] = 4
    hdr["Y_IND"] = 5
    hdr["EBV_REF"] = 0.0
    hdr["COMMENT"] = "cosmoDC2 photometric catalog for redMaPPer (EBV=0, no MW dust)"

    hdu_primary = fits.PrimaryHDU(header=hdr)
    hdu_table = fits.BinTableHDU.from_columns(cols, header=hdr)
    hdul = fits.HDUList([hdu_primary, hdu_table])
    hdul.writeto(out_path, overwrite=True)
    print(f"Wrote {n:,} rows to {out_path}")
    '''

    gal_dtype = [("id", "i8"), ("ra", "f8"), ("dec", "f8"), 
                 ("refmag", "f4"), ("refmag_err", "f4"), ("mag", "f4", n_bands), ("mag_err", "f4", n_bands), 
                 ("ebv", "f4"), ("ztrue", "f4")] # sim-only, for comparison
    
    galaxies = np.zeros(n, dtype=gal_dtype)
    galaxies["id"] = gal_id
    galaxies["ra"] = ra
    galaxies["dec"] = dec
    galaxies["refmag"] = refmag
    galaxies["refmag_err"] = refmag_err
    galaxies["mag"] = mags
    galaxies["mag_err"] = mag_errs
    galaxies["ebv"] = ebv
    galaxies["ztrue"] = np.array(data["redshift_true"], dtype=np.float32)

    info_dict = {"LIM_REF": mag_cut,
        "REF_IND": ref_i,
        "AREA": area,
        "NMAG": n_bands,
        "MODE": "LSST",
        "ZP": zeropoint,
        "G_IND": 0, "R_IND": 1, "I_IND": 2, "Z_IND": 3, "Y_IND": 4}

    maker = redmapper.GalaxyCatalogMaker(out_path, info_dict)
    maker.append_galaxies(galaxies)
    maker.finalize_catalog()
    print(f"Wrote healpix-split catalog with base '{out_path}' ({n:,} galaxies)") #had to change this from fits to healpix file
    return ra, dec, mags, refmag

## red-sequence template for initial guesses
def build_redsequence_template(cat, zmin, zmax, out_path):
    '''
    Build redsequence as function of redshift
    '''

    mag_cols = [f"mag_true_{b}_lsst" for b in BANDS]
    quants = ["redshift_true", "is_central"] + mag_cols

    filters = ["is_central == True", f"redshift_true >= {zmin}", f"redshift_true <= {zmax}"]

    data = cat.get_quantities(quants, filters=filters) 
    z = np.array(data["redshift_true"]) 
    mags = np.column_stack([np.array(data[c]) for c in mag_cols])

    # select on g-r
    # gr_obs = mags[:,1] - mags[:, 2]
    gr_obs = mags[:, 0] - mags[:, 1] # g-r selection
    red = gr_obs > 0.5
    bins = np.arange(zmin, zmax+dz_rs, dz_rs) 
    z_out, gr_out, ri_out, iz_out, zy_out = [], [], [], [], []

    # getting median red-seq colors per redshift band dz_rs
    for zlo, zhi in zip(bins[:-1], bins[1:]):
        mask = (z >= zlo) & (z < zhi) & red
        if mask.sum() < 5:
            continue
        
        zmid = 0.5 * (zlo + zhi)
        m = mags[mask]
        
        # gr = np.median(m[:, 1] - m[:, 2])
        # ri = np.median(m[:, 2] - m[:, 3])
        # iz = np.median(m[:, 3] - m[:, 4])
        # zy = np.median(m[:, 4] - m[:, 5])
        gr = np.median(m[:, 0] - m[:, 1]) # g-r
        ri = np.median(m[:, 1] - m[:, 2]) # r-i
        iz = np.median(m[:, 2] - m[:, 3]) # i-z
        zy = np.median(m[:, 3] - m[:, 4]) # z-y
        
        z_out.append(zmid)
        gr_out.append(gr)
        ri_out.append(ri)
        iz_out.append(iz)
        zy_out.append(zy)

    #write table
    header = ("# cosmoDC2 initial color template for redMaPPer\n")
    with open(out_path, "w") as f: 
        f.write(header) 
        for vals in zip(z_out, gr_out, ri_out, iz_out, zy_out):
            f.write(f"  {vals[0]:.4f}  {vals[1]:.4f}  {vals[2]:.4f}  {vals[3]:.4f}  {vals[4]:.4f}\n")

    return np.array(z_out), np.array(gr_out), np.array(ri_out), np.array(iz_out), np.array(zy_out)

def parse_args():
    p = argparse.ArgumentParser(description="Build redMaPPer inputs (spec-z, photometric, red-sequence) from cosmoDC2")
    p.add_argument("--catalog", default="cosmoDC2_v1.1.4_image")
    p.add_argument("--n-per-bin", type=int,   default=25,
                   help="Spec-z centrals per dz=0.05 bin (default 25)")
    p.add_argument("--zmin", type=float, default=0.10)
    p.add_argument("--zmax", type=float, default=0.90)
    p.add_argument("--mag-cut", type=float, default=26.0, help="z-band limiting magnitude for photometric catalog")
    p.add_argument("--out-dir", default="./",help="Output directory for all files")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--specz-only", action="store_true", help="Only rebuild the spec-z catalog; skip photo + template")
    return p.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print(f"Loading catalog: {args.catalog}")
    cat = GCRCatalogs.load_catalog(args.catalog)

    specz_path = os.path.join(args.out_dir, "cosmodc2_specz.fits")
    photo_path = os.path.join(args.out_dir, "cosmodc2_photo") #used to have .fits
    rs_path = os.path.join(args.out_dir, "cosmodc2_redsequence_init.fit")

    build_specz_catalog(cat, args.zmin, args.zmax, args.n_per_bin, args.seed, specz_path)
    if not args.specz_only:
        build_photometric_catalog(cat, args.zmin, args.zmax, args.mag_cut, photo_path)
        build_redsequence_template(cat, args.zmin, args.zmax, rs_path)
    else:
        print("--specz-only: skipped photometric catalog and red-sequence template")

    print("\n All inputs written")
    print(f" Spec-z training: {specz_path}")
    print(f" Photometric catalog: {photo_path}")
    print(f" Red-sequence init: {rs_path}")
    print()


if __name__ == "__main__":
    main()