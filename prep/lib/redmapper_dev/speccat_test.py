import redmapper
import healpy as hp
import tables_io, os, sys, yaml
import numpy as np
from astropy.table import Table, vstack
from tqdm import tqdm


config_file = sys.argv[1]
with open(config_file) as fstream :
    cfg = yaml.safe_load(fstream)


input_base = cfg['iFiles']['base']
input_files = [os.path.join(input_base, fileName) for fileName in os.listdir(input_base)]

spec_dtype = list(cfg['specCat']['dtype'].items())

catKeys = cfg['specCat']['catKeys']

bands = cfg['bands']
magNames = [cfg['magFmt'].replace('*', band) for band in bands]

galaxies = []
for input_file in tqdm(input_files) :
    # insert code to translate to file format
    og = tables_io.read(input_file)
    og = vstack([Table(og[fkey]) for fkey in og.keys() if (fkey != 'metaData')])
    galaxies.append(Table())

    is_center = (og[catKeys['is_center']] == cfg['specCat']['isCenterFlag'])
    is_massive = (og[catKeys['halo_mass']] > float(cfg['specCat']['massCut']))

    for key in catKeys.keys() :
        galaxies[-1][key] = og[catKeys[key]][is_center & is_massive]
    
    galaxies[-1]['mag'] = np.array([og[mag][is_center & is_massive] for mag in magNames]).T
    galaxies[-1]['mag_err'] = np.ones_like(galaxies[-1]['mag']) * 27. / galaxies[-1]['mag']
    galaxies[-1]['refmag_err'] = np.ones_like(galaxies[-1]['refmag']) * 27. / galaxies[-1]['refmag']

    galaxies[-1]['z_err'] = 0.0001 * (1 + galaxies[-1]['z'])
    galaxies[-1]['ebv'] = np.zeros_like(galaxies[-1]['refmag'])

    del og

galaxies_vstacked = vstack(galaxies)

## NOW GRAB ~40 GALAXIES PER REDSHIFT BIN
zbins = np.linspace(
        cfg['specCat']['zbins']['min'],
        cfg['specCat']['zbins']['max'],
        cfg['specCat']['zbins']['nbins'])

ngals_bin = cfg['specCat']['zDensity']['ngals'] / cfg['specCat']['zDensity']['width']
ngals_bin = int(ngals_bin * np.diff(zbins)[0])

indices = np.digitize(galaxies_vstacked['z'], bins=zbins, right=True)
zselection = np.array(np.zeros_like(galaxies_vstacked['z']).astype(bool))
for i in range(len(zbins)-1) :
    in_bin = (indices == i)
    if sum(in_bin) > ngals_bin :
        random_selection = np.zeros_like(zselection[in_bin]).astype(bool)
        random_selection[:ngals_bin] = True
        np.random.shuffle(random_selection)
        zselection[in_bin] |= random_selection

galaxies_vstacked = galaxies_vstacked[zselection]

tables_io.write(galaxies_vstacked, os.path.join(cfg['oFiles']['base'], cfg['oFiles']['specz']))
