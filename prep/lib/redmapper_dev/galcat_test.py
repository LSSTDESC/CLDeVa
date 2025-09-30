import redmapper
import healpy as hp
import tables_io, os, sys, yaml
import numpy as np
from astropy.table import Table, vstack
from tqdm import tqdm

config_file = sys.argv[1]
with open(config_file) as fstream :
    cfg = yaml.safe_load(fstream)


iFiles = [os.path.join(cfg['iFiles']['base'], fileName)
        for fileName in os.listdir(cfg['iFiles']['base'])]
oFiles = os.path.join(cfg['oFiles']['base'], cfg['oFiles']['galcat'])
# files will show up as '/path/to/files/foo_XXXXX.fits' (where XXXXX is a
# healpix number) and the mail galaxy file will be
# '/path/to/files/foo_master_table.fits'

info_dict = cfg['photCat']['info_dict']
photDType = list(cfg['photCat']['dtype'].items())
for i in np.where([p[0] in ['mag','mag_err'] for p in photDType])[0] :
    photDType[i] = photDType[i] + (info_dict['NMAG'],)

maker = redmapper.GalaxyCatalogMaker(oFiles, info_dict, nside=cfg['photCat']['nside'])

catKeys = cfg['photCat']['catKeys']

magNames = [cfg['magFmt'].replace('*', band) for band in cfg['bands']]
for iFile in tqdm(iFiles) :
    og = tables_io.read(iFile)
    og = vstack([Table(og[fkey]) for fkey in og.keys() if (fkey != 'metaData')])
    galaxies = Table()
    for key in catKeys.keys() :
        galaxies[key] = og[catKeys[key]]
    
    galaxies['refmag_err'] = np.ones_like(galaxies['refmag']) * info_dict['ZP'] / galaxies['refmag']

    galaxies['mag'] = np.array([og[magName] for magName in magNames]).T
    galaxies['mag_err'] = np.ones_like(galaxies['mag']) * info_dict['ZP'] / galaxies['mag']

    galaxies['ebv'] = np.zeros_like(galaxies['refmag'])
    
    maker.append_galaxies(galaxies.as_array())

maker.finalize_catalog()
