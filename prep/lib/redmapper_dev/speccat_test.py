import redmapper
import healpy as hp
import tables_io, os
import numpy as np
from astropy.table import Table, vstack
from tqdm import tqdm

input_base = '/sps/lsst/groups/desc/shared/xgal/roman-rubin/roman_rubin_2023_v1.1.3/'
input_files = [os.path.join(input_base, fileName) for fileName in os.listdir(input_base)]
# files will show up as '/path/to/files/foo_XXXXX.fits' (where XXXXX is a
# healpix number) and the mail galaxy file will be
# '/path/to/files/foo_master_table.fits'


spec_dtype = [('id', 'i8'),
              ('ra', 'f8'),        # right ascension (degrees)
              ('dec', 'f8'),       # declination (degrees)
              ('z', 'f4'),         # spectrosopic redshift
              ('z_err', 'f4')]     # error on spec-z


info_dict = {}
info_dict['LIM_REF'] = 24.6
info_dict['REF_IND'] = 4
info_dict['AREA'] = hp.nside2pixarea(64, degrees=True)
info_dict['NMAG'] = 6
info_dict['MODE'] = 'LSST'
info_dict['ZP'] = 31.4
info_dict['U_IND'] = 0 # u-band index
info_dict['G_IND'] = 1 # g-band index
info_dict['R_IND'] = 2 # r-band index
info_dict['I_IND'] = 3 # i-band index
info_dict['Z_IND'] = 4 # z-band index
info_dict['Y_IND'] = 5 # y-band index

pkeys = {
        'id': 'galaxy_id',
        'ra': 'ra',
        'dec': 'dec',
        'z': 'redshift',
        #'z_err': 
        'refmag': 'LSST_obs_z',
        }

bands = ['u','g','r','i','z','y']
galaxies = []
for input_file in tqdm(input_files[:2]) :
    # insert code to translate to file format
    og = tables_io.read(input_file)
    og = vstack([Table(og[fkey]) for fkey in og.keys() if (fkey != 'metaData')])
    galaxies.append(Table())

    is_center = (og['um_source_galaxy_upid'] == -1)
    for pkey in pkeys.keys() :
        galaxies[-1][pkey] = og[pkeys[pkey]][is_center]
    
    galaxies[-1]['mag'] = np.array([og[f'LSST_obs_{band}'][is_center] for band in bands]).T
    galaxies[-1]['mag_err'] = np.ones_like(galaxies[-1]['mag']) * info_dict['ZP'] / galaxies[-1]['mag']
    galaxies[-1]['refmag_err'] = np.ones_like(galaxies[-1]['refmag']) * info_dict['ZP'] / galaxies[-1]['refmag']

    galaxies[-1]['z_err'] = 0.0001 * (1 + galaxies[-1]['z'])
    galaxies[-1]['ebv'] = np.zeros_like(galaxies[-1]['refmag'])
    del og

galaxies_vstacked = vstack(galaxies)

## NOW GRAB ~40 GALAXIES PER REDSHIFT BIN
zbins = np.linspace(0,1.5,301)
indices = np.digitize(galaxies_vstacked['z'], bins=zbins, right=True)
zselection = np.zeros_like(galaxies_vstacked['z']).astype(bool)
for i in range(len(zbins)-1) :
    in_bin = (indices == i)
    if sum(in_bin) > (2 * 40) :
        zselection[in_bin] |= np.random.choice([True, False], sum(in_bin), [40/sum(in_bin), (1-40/sum(in_bin))])

galaxies_vstacked = galaxies_vstacked[zselection]

tables_io.write(galaxies_vstacked, './input_data/specz/roman_rubin_specz.fits', fmt='fits')
