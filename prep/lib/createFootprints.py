
from astropy.table import Table
from astropy.io import fits
import numpy as np
import healpy as hp
import healsparse as hs
import tables_io
import os, sys, yaml, glob

src = __file__.replace(os.path.basename(__file__), '')

with open(sys.argv[1]) as fstream :
    cfg = yaml.safe_load(fstream)

footDir = os.path.join(cfg['outpath']['base'], cfg['name'], cfg['outpath']['footprint'])
Nside_cov = cfg['outpath']['hp_nside']

## get mask info
if 'masks' in cfg.keys() :
    mask_nsides = [cfg['masks'][mask]['hp_nside'] for mask in cfg['masks'].keys()]
    mask_files  = [cfg['masks'][mask]['file'] for mask in cfg['masks'].keys()]

    mask_map = hs.HealSparseMap.make_empty(Nside_cov, max(mask_nsides), bool)
    for i in range(len(cfg['masks'].keys())) :
        tmp_mask_map = hs.HealSparseMap.make_empty(Nside_cov, mask_nsides[i], bool)
        tmp_mask_pixels = np.load(mask_files[i])
        tmp_mask_map.update_values_pix(tmp_mask_pixels, np.ones(len(tmp_mask_pixels)).astype(bool), nest=True)
        if mask_nsides[i] != max(mask_nsides) :
            tmp_mask_map.upgrade(max(mask_nsides))
        mask_map |= tmp_mask_map
else :
    mask_map = hs.HealSparseMap.make_empty(Nside_cov, Nside_cov, bool)


## get footprint info (mask is True to keep -- ringed)
ftp_map = hs.HealSparseMap.make_empty(Nside_cov, cfg['footprint']['hp_nside'], bool)
ftp_map_pixels = np.load(cfg['footprint']['file'])
ftp_map.update_values_pix(ftp_map_pixels, np.ones(len(ftp_map_pixels)).astype(bool), nest=True)


## go with higher resolution of the two
Nside = max(ftp_map.nside_sparse, mask_map.nside_sparse)

## bring both maps to same resolution
if Nside > ftp_map.nside_sparse :
    ftp_map = ftp_map.upgrade(Nside)
elif Nside > mask_map.nside_sparse :
    mask_map  = mask_map.upgrade(Nside)

## create masked footprint
maskedFootprint_map = ftp_map
maskedFootprint_map &= mask_map

## export footprint in chunks of coverage pixels
covPixels = np.where(maskedFootprint_map.coverage_map.astype(bool))[0]
for covPix in covPixels :
    print(covPix)

    covPix_map = maskedFootprint_map.get_single_covpix_map(covPix)

    pixTable = Table()
    pixTable["pixel"] = covPix_map.valid_pixels
    pixTable["ra"], pixTable["dec"] = covPix_map.valid_pixels_pos(lonlat=True)
    pixTable["signal"] = np.ones(pixTable["ra"].size)

    footFile = os.path.join(footDir, f"{covPix}_footprint.fits")
    pixTable.write(footFile, overwrite=True)
    del pixTable

