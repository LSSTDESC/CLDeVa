import healsparse as hs
import healpy as hp
import numpy as np
import tables_io


specz = tables_io.read('/pbs/throng/lsst/users/rsolomon/redmapper/data_prep/input_data/specz/roman_rubin_specz.fits')

pixels = np.unique(hp.ang2pix(4096, specz['ra'], specz['dec'], nest=False, lonlat=True))

hsmap = hs.HealSparseMap.make_empty(64, 4096, bool,)

hsmap.update_values_pix(pixels, np.ones(len(pixels)).astype(bool), nest=False)

hsmap.write('/pbs/throng/lsst/users/rsolomon/redmapper/data_prep/masks/roman_rubin_geoMask.hs', clobber=False)
