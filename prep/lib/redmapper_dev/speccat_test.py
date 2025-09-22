import redmapper
import healpy as hp
import tables_io, os
import numpy as np
from astropy.table import Table, vstack
from tqdm import tqdm

filename_base = 'input_data/specz/'
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

#redmapper_dtype = [('id', 'i8'),             # galaxy id number (unique)
#                   ('ra', 'f8'),             # right ascension (degrees)
#                   ('dec', 'f8'),            # declination (degrees)
#                   ('refmag', 'f4'),         # total magnitude in reference band
#                   ('refmag_err', 'f4'),     # error in total reference mag
#                   ('mag', 'f4', 6),      # mag array
#                   ('mag_err', 'f4', 6),  # magnitude error array
#                   ('ebv', 'f4'),            # E(B-V) (systematics checking)
#                   ('ztrue', 'f4'),]          # ztrue if from a simulated catalog
#                   #('m200', 'f4'),           # m200 of halo if from a simulated catalog
#                   #('central', 'i2'),        # central? 1 if yes (if from sims)
#                   #('halo_id', 'i8')]        # halo_id if from a simulated catalog

info_dict = {}
info_dict['LIM_REF'] = 24.6
info_dict['REF_IND'] = 4
info_dict['AREA'] = hp.nside2pixarea(64, degrees=True)
info_dict['NMAG'] = 6
info_dict['MODE'] = 'LSST'
info_dict['ZP'] = 31.4
#info_dict['B'] = b_array # if magnitudes are actually luptitudes
info_dict['U_IND'] = 0 # u-band index
info_dict['G_IND'] = 1 # g-band index
info_dict['R_IND'] = 2 # r-band index
info_dict['I_IND'] = 3 # i-band index
info_dict['Z_IND'] = 4 # z-band index
info_dict['Y_IND'] = 5 # y-band index
# (etc, for the rest of the bands)

maker = redmapper.GalaxyCatalogMaker(filename_base, info_dict, nside=int(64*2))

pkeys = {
        'id': 'galaxy_id',
        'ra': 'ra',
        'dec': 'dec',
        'z': 'redshift',
        #'z_err': 
        'refmag': 'LSST_obs_z',
        #'refmag_err': ,
        #'mag': ,
        #'mag_err': ,
        #'ebv': ,
        }

bands = ['u','g','r','i','z','y']
for input_file in tqdm(input_files) :
    # insert code to translate to file format
    og = tables_io.read(input_file)
    og = vstack([Table(og[fkey]) for fkey in og.keys() if (fkey != 'metaData')])
    galaxies = Table()

    is_center = (og['um_source_galaxy_upid'] == -1)
    ## NOW GRAB ~40 GALAXIES PER REDSHIFT BIN
    zbins = np.linspace(0,1.5,151)
    indices = np.digitize(og['z'][is_center], bins=zbins, right=True)
    zselection = np.zeros_like(og['z'][is_center]).astype(bool)
    for i in range(len(zbins)-1) :
        in_bin = (indices == i)
        zselection[in_bin] |= np.random.choice([False,True], sum(in_bin))
    for pkey in pkeys.keys() :
        galaxies[pkey] = og[pkeys[pkey]][is_center][zselection]
    
    galaxies['mag'] = np.array([og[f'LSST_obs_{band}'][is_center][zselection] for band in bands]).T
    galaxies['mag_err'] = np.ones_like(galaxies['mag']) * info_dict['ZP'] / galaxies['mag']
    galaxies['refmag_err'] = np.ones_like(galaxies['refmag']) * info_dict['ZP'] / galaxies['refmag']

    galaxies['z_err'] = 0.001 * (1 + galaxies['z'])
    galaxies['ebv'] = np.zeros_like(galaxies['refmag'])


    maker.append_galaxies(galaxies.as_array())

maker.finalize_catalog()
