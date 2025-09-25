import healsparse as hs
import healpy as hp
import numpy as np
import tables_io


## read in specz data just to know what HEALPix pixels are in footprint
specz = tables_io.read('/pbs/throng/lsst/users/rsolomon/redmapper/data_prep/input_data/specz/roman_rubin_specz.fits')
pixels = np.unique(hp.ang2pix(4096, specz['ra'], specz['dec'], nest=False, lonlat=True))


## make the data type for the map
depth_dtype = [('exptime', 'f4'),  # effective exposure time
               ('limmag', 'f4'),   # limited magnitude at nsig (in header)
               ('m50', 'f4'),      # Should be same as LIMMAG (for now)
               ('fracgood', 'f4')] # fraction of good coverage (see mask above)


## initialize the map
depth_map = hs.HealSparseMap.make_empty(64, 4096, depth_dtype, primary='fracgood')
depth_map[pixels] = np.zeros(len(pixels), dtype=depth_dtype)


## fill each map array (using arbitrary values for now)
depth_map['exptime'][depth_map.valid_pixels] += 30
depth_map['limmag'][depth_map.valid_pixels] += 30
depth_map['m50'][depth_map.valid_pixels] += 30
depth_map['fracgood'][depth_map.valid_pixels] += 1


## add metadata to 
metadata = {'ZP': 31.4,
            'NSIG': 10.0,
            'NBAND': 1,
            'W': 0.0,
            'EFF':1.0,
           }
depth_map.metadata = metadata

depth_map.write('/pbs/throng/lsst/users/rsolomon/redmapper/data_prep/masks/roman_rubin_depthMask.hs', clobber=False)
