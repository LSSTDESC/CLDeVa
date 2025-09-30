import healsparse as hs
import os, sys, yaml
import healpy as hp
import numpy as np
import tables_io

cfgFile = sys.argv[1]
with open(cfgFile) as fstream :
    cfg = yaml.safe_load(fstream)

specz = tables_io.read(os.path.join(cfg['oFiles']['base'], cfg['oFiles']['specz']))
pixels = np.unique(
        hp.ang2pix(
            cfg['geoMask']['nside_sparse'],
            specz['ra'],
            specz['dec'],
            nest=cfg['geoMask']['nest'],
            lonlat=True))
del specz

geoMap = hs.HealSparseMap.make_empty(
        cfg['geoMask']['nside_cov'],
        cfg['geoMask']['nside_sparse'],
        bool,)

geoMap.update_values_pix(pixels, np.ones(len(pixels)).astype(bool), nest=cfg['geoMask']['nest'])

geoMap.write(os.path.join(cfg['oFiles']['base'], cfg['oFiles']['geo']), clobber=False)
