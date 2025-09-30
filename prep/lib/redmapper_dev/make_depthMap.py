import healsparse as hs
import os, sys, yaml
import healpy as hp
import numpy as np
import tables_io

cfgFile = sys.argv[1]
with open(cfgFile) as fstream :
    cfg = yaml.safe_load(fstream)

## read in specz data just to know what HEALPix pixels are in footprint
specz = tables_io.read(os.path.join(cfg['oFiles']['base'], cfg['oFiles']['specz']))
pixels = np.unique(
        hp.ang2pix(
            cfg['depthMap']['nside_sparse'],
            specz['ra'],
            specz['dec'],
            nest=cfg['depthMap']['nest'],
            lonlat=True))
del specz

## make the data type for the map
depth_dtype = list(cfg['depthMap']['dtype'].items())

## initialize the map
depth_map = hs.HealSparseMap.make_empty(
        cfg['depthMap']['nside_cov'],
        cfg['depthMap']['nside_sparse'],
        depth_dtype,
        primary='fracgood')
depth_map[pixels] = np.zeros(len(pixels), dtype=depth_dtype)


## fill each map array (using arbitrary values for now)
for key in cfg['depthMap']['initVals'].keys() :
    depth_map[key][depth_map.valid_pixels] += cfg['depthMap']['initVals'][key]

## add metadata to 
depth_map.metadata = cfg['depthMap']['metaData']

depth_map.write(os.path.join(cfg['oFiles']['base'], cfg['oFiles']['depth']), clobber=False)
