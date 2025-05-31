
from astropy.table import Table, join
import numpy as np
import os, sys, yaml, glob


with open(sys.argv[1]) as fstream :
	cfg = yaml.safe_load(fstream)

path = os.path.join(cfg['outpath']['base'], cfg['name'], cfg['outpath']['galaxies'])

pixtract_files = glob.glob(path+'*')
pixels = np.unique([os.path.basename(pt_f).split('_')[0] for pt_f in pixtract_files])

for pix in pixels :
    print(f'{pix}', end=', ')
    
    tracts = glob.glob(f'{path}{pix}_*')
    t = Table.read(tracts[0], format='fits')
    for tract in tracts[1:] :
    	tmp = Table.read(tract, format='fits')
    	t = join(t, tmp, join_type='outer')
    	
    pixFile = f'{path}{pix}.fits'
    t.write(pixFile, overwrite=True)
    [os.remove(f'{tract}') for tract in tracts]
    
    print('DONE')
    del t
