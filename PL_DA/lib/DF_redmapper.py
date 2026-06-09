import DF_utils as utils
import redmapper
from astropy.table import Table, vstack
import tables_io, os, sys, yaml, warnings
import numpy as np
from tqdm import tqdm
import GCRCatalogs as gcr



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

    def _load_from_gcr(self) :
        if not hasattr(self, 'gcr_DF') :
            warnings.filterwarnings('ignore', category=RuntimeWarning)
            os.environ['GCR_CONFIG_SOURCE'] = 'files'
            self.gcr_DF = gcr.load_catalog(self.config['gcr_name'])

        filters = ()
        native_filters = ()
        if 'gcr_filters' in self.config.keys() :
            for fltr in self.config['gcr_filters'].keys() :
                for constraint in self.config['gcr_filters'][fltr] :
                    filters = (*filters, f"{fltr}{constraint}")

        gcr_cols = np.unique(
                [self.config['specCat']['catKeys'][key] for key in self.config['specCat']['catKeys'].keys()]
                + [self.config['photCat']['catKeys'][key] for key in self.config['photCat']['catKeys'].keys()]
                + [self.config['magFmt'].replace('*', band) for band in self.config['bands']]).tolist()

        self.gcr_cat = self.gcr_DF.get_quantities(
                gcr_cols,
                filters=filters,
                native_filters=native_filters)

    
    def get_spec_catalog(self) :
        if self.config['DA_type'].upper() == 'FILE' :
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
        elif self.config['DA_type'].upper() == 'GCR' :
            if not hasattr(self, 'gcr_cat') :
                self._load_from_gcr()
            galaxies = Table()
            for col in self.config['specCat']['catKeys'].keys() :
                galaxies[col] = self.gcr_cat[self.config['specCat']['catKeys'][col]]
            

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
        outFiles = os.path.join(self.config['oFiles']['base'], self.config['oFiles']['galcat'])
        # files will show up as '/path/to/files/foo_XXXXX.fits' (where XXXXX is a
        # healpix number) and the main galaxy file will be
        # '/path/to/files/foo_master_table.fits'
        
        for i in np.where([p[0] in ['mag', 'mag_err'] for p in self.phot_dtype])[0] :
            self.phot_dtype[i] = self.phot_dtype[i] + (self.phot_dict['NMAG'],)

        maker = redmapper.GalaxyCatalogMaker(outFiles, self.phot_dict, nside=self.config['photCat']['nside'])

        if self.config['DA_type'].upper() == 'FILE' :
            inPath = self.config['iFiles']['base']
            inFiles = [os.path.join(inPath, fileName) for fileName in os.listdir(inPath)]
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
        elif self.config['DA_type'].upper() == 'GCR' :
            if not hasattr(self, 'gcr_cat') :
                self._load_from_gcr()
            galaxies = Table()
            for col in self.config['photCat']['catKeys'].keys() :
                galaxies[col] = self.gcr_cat[self.config['photCat']['catKeys'][col]]
            galaxies['refmag_err'] = np.ones_like(galaxies['refmag']) * 0.01
            galaxies['mag'] = np.array([self.gcr_cat[magName] for magName in self.mag_names]).T
            galaxies['mag_err'] = np.ones_like(galaxies['mag']) * 0.01
            galaxies['ebv'] = np.zeros_like(galaxies['refmag'])
            
            maker.append_galaxies(galaxies.as_array())

        maker.finalize_catalog()
