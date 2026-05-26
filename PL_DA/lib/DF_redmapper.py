import DF_utils as utils
import redmapper
from astropy.table import Table, vstack
import tables_io, os, sys, yaml
import numpy as np
from tqdm import tqdm



class DF_redmapper :

    def __init__(self, config) :
        if isinstance(config, dict) :
            self.config = config
        elif isinstance(config, str) :
            try :
                self.config = utils._yaml2dict(config)
            except ValueError :
                print(f"\tFILE DOES NOT EXIST: {config}")
                sys.exit()
        self.spec_dtype = list(self.config['specCat']['dtype'].items())
        self.spec_cols = self.config['specCat']['catKeys']

        self.phot_dict = self.config['photCat']['info_dict']
        self.phot_dtype = list(self.config['photCat']['dtype'].items())
        self.phot_cols = self.config['photCat']['catKeys']
        
        self.bands = self.config['bands']
        self.mag_names = [self.config['magFmt'].replace('*', band) for band in self.bands]


    def get_spec_catalog(self) :
        inPath = self.config['iFiles']['base']
        inFiles = [os.path.join(inPath, fileName) for fileName in os.listdir(inPath)]
        
        galaxies = []
        for inFile in tqdm(inFiles) :
            galaxies.append(Table())
            inData = tables_io.read(inFile)
            inData = vstack([Table(inData[fkey]) for fkey in inData.keys() if (fkey != 'metaData')]) 

            is_center = (inData[self.spec_cols['is_center']] == self.config['specCat']['isCenterFlag'])
            is_massive = (inData[self.spec_cols['halo_mass']] > float(self.config['specCat']['massCut']))

            for col in self.spec_cols.keys() :
                galaxies[-1][col] = inData[self.spec_cols[col]][is_center & is_massive]

            galaxies[-1]['mag'] = np.array([inData[mag][is_center & is_massive] for mag in self.mag_names]).T
            galaxies[-1]['mag_err'] = np.ones_like(galaxies[-1]['mag']) * 27. / galaxies[-1]['mag']
            galaxies[-1]['refmag_err'] = np.ones_like(galaxies[-1]['refmag']) * 27. / galaxies[-1]['refmag']
            galaxies[-1]['z_err'] = 0.0001 * (1 + galaxies[-1]['z'])
            galaxies[-1]['ebv'] = np.zeros_like(galaxies[-1]['refmag'])

            del inData

        galaxies = vstack(galaxies)

        zbins = np.linspace(
                self.config['specCat']['zbins']['min'],
                self.config['specCat']['zbins']['max'],
                self.config['specCat']['zbins']['nbins'])

        ngals_bin = self.config['specCat']['zDensity']['ngals'] / self.config['specCat']['zDensity']['width']
        ngals_bin = int(ngals_bin * np.diff(zbins)[0])

        indices = np.digitize(galaxies['z'], bins=zbins, right=True)
        zselection = np.array(np.zeros_like(galaxies['z']).astype(bool))

        for i in range(len(zbins) - 1) :
            in_bin = (indices == i)
            if sum(in_bin) > ngals_bin :
                random_selection = np.zeros_like(zselection[in_bin]).astype(bool)
                random_selection[:ngals_bin] = True
                np.random.shuffle(random_selection)
                zselection[in_bin] |= random_selection

        self.spec_sample = galaxies[zselection]


    def write_spec_catalog(self) :
        if hasattr(self, 'spec_sample'):
            tables_io.write(self.spec_sample, os.path.join(self.config['oFiles']['base'], self.config['oFiles']['specz']))
        else :
            self.get_spec_catalog()
            tables_io.write(self.spec_sample, os.path.join(self.config['oFiles']['base'], self.config['oFiles']['specz']))


    def write_phot_catalog(self) :
        inPath = self.config['iFiles']['base']
        inFiles = [os.path.join(inPath, fileName) for fileName in os.listdir(inPath)]
        outFiles = os.path.join(self.config['oFiles']['base'], self.config['oFiles']['galcat'])
        # files will show up as '/path/to/files/foo_XXXXX.fits' (where XXXXX is a
        # healpix number) and the main galaxy file will be
        # '/path/to/files/foo_master_table.fits'
        
        for i in np.where([p[0] in ['mag', 'mag_err'] for p in self.phot_dtype])[0] :
            self.phot_dtype[i] = self.phot_dtype[i] + (info_dict['NMAG'],)

        maker = redmapper.GalaxyCatalogMaker(outFiles, info_dict, nside=self.config['photCat']['nside'])

        for inFile in tqdm(inFiles) :
            inData = tables_io.read(inFile)
            inData = vstack([Table(inData[fkey]) for fkey in inData.keys() if (fkey != 'metaData')])

            galaxies = Table()

            for col in self.phot_cols.keys() :
                galaxies[col] = inData[self.phot_cols[col]]

            galaxies['refmag_err'] = np.ones_like(galaxies['refmag']) * self.phot_dict['ZP'] / galaxies['refmag']
            galaxies['mag'] = np.array([inData[magName] for magName in self.mag_names]).T
            galaxies['mag_err'] = np.ones_like(galaxies['mag']) * self.phot_dict['ZP'] / galaxies['mag']
            galaxies['ebv'] = np.zeros_like(galaxies['refmag'])

            maker.append_galaxies(galaxies.as_array())

        maker.finalize_catalog()
