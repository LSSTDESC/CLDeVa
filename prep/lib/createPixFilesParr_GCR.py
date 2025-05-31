
from astropy.table import Table
from scipy.stats import norm
import numpy as np
import healpy as hp
import healsparse as hs
import os, sys, yaml, warnings
import tables_io
import GCRCatalogs as gcr
from utils import get_pixels


with open(sys.argv[1]) as fstream :
	cfg = yaml.safe_load(fstream)

pixel = int(sys.argv[2])     ## nested pixel ID

data = Table()

gal_keys = cfg['gal_cat']['keys']
bands = cfg['gal_cat']['filter_bands']
Nside_cov = cfg['outpath']['hp_nside']
if 'masks' in cfg.keys() :
    mask_nsides = max([cfg['masks'][mask]['hp_nside'] for mask in cfg['masks'].keys()])
else :
    mask_nsides = 32
Nside_sparse = max(cfg['footprint']['hp_nside'], mask_nsides)

## get galaxy data from GCR
warnings.filterwarnings("ignore", category=RuntimeWarning)
cat_name = cfg['gal_cat']['gcr_name']
print(f'Getting GCR catalog: {cat_name}', end='...')
gcr_cat = gcr.load_catalog(cat_name)


filters = ()
native_filters = ()
if 'gcr_filters' in cfg['gal_cat'] :
    filters += tuple(cfg['gal_cat']['gcr_filters'])

## GCR catalogs, if with healpix information, is saved with nside=32 and ringed indexing.
## Convert pixel to nside=32 and ringed to retrieve GCR data.
if Nside_cov > 32 :
    gcr_hsMap = hs.HealSparseMap.make_empty(32, Nside_cov, bool)
    gcr_hsMap.update_values_pix(np.array([pixel,]), np.ones(1).astype(bool), nest=True)
    gcr_pixel = np.where(gcr_hsMap.coverage_map.astype(bool))[0][0]
    gcr_pixel = hp.nest2ring(32, gcr_pixel)
elif Nside_cov == 32 :
    gcr_pixel = hp.nest2ring(32, pixel)

if gcr_cat.has_quantity('cosmodc2_hp_truth') :
    filters += (f'cosmodc2_hp_truth == {gcr_pixel}',)
if 'healpix_pixel' in gcr_cat.native_filter_quantities :
    native_filters += (f'healpix_pixel == {gcr_pixel}',)

gcr = gcr_cat.get_quantities(
        [gal_keys[key] for key in gal_keys.keys()],
        filters=filters,
        native_filters=native_filters)

for key in gal_keys.keys() :
    data[key] = gcr[gal_keys[key]]
del gcr_cat, gcr
print('DONE')


## apply footprint
ftp_map = hs.HealSparseMap.make_empty(Nside_cov, Nside_sparse, bool)
ftp_file = os.path.join(cfg['outpath']['base'], cfg['name'], cfg['outpath']['footprint'], f"{pixel}_footprint.fits")
ftp_pixels = np.array(tables_io.read(ftp_file)['pixel'].tolist())
ftp_map.update_values_pix(ftp_pixels, np.ones(len(ftp_pixels)).astype(bool), nest=True)

galPix = hp.ang2pix(Nside_sparse, data['ra'], data['dec'], nest=True, lonlat=True)
data_in_ftp = np.isin(galPix, ftp_map.valid_pixels)

if sum(data_in_ftp) == 0 :
    sys.exit("EXIT CODE: healpix completely masked out")
data = data[data_in_ftp]


## get redshifts if not in GCR
if 'redshift' not in gal_keys.keys() :
    print('Getting redshifts', end='...')
    ## PZ files also stored with same healpix designation
    zFile = os.path.join(cfg['redshift_cat']['inpath'], cfg['redshift_cat']['fname'].replace('TMP', str(gcr_pixel)))
    if (suffix:=zFile.split('.')[-1]) == 'zip' :
        fmt = zFile.split('.')[-2]
    else :
        fmt = suffix
    zcols = [cfg['redshift_cat']['keys'][k] for k in cfg['redshift_cat']['keys']]
    zData = tables_io.read(zFile, columns=zcols, fmt=fmt)
    
    ## ensure identical ordering between redshift and GCR catalogs
    _, order_gals, order_pzs = np.intersect1d(data['id'], zData['id'], assume_unique=True, return_indices=True)
    data = data[order_gals]
    zData = zData.loc[order_pzs]
    if len(order_pzs) == 0 :
        sys.exit("EXIT CODE: no photo-z's were found in this healpix")
    
    data['redshift'] = np.hstack(zData[cfg['redshift_cat']['keys']['redshift']])
    if 'redshift_err' in cfg['redshift_cat']['keys'] :
        data['redshift_err'] = np.hstack(zData[cfg['redshift_cat']['keys']['redshift_err']])
    del zData
    print('DONE')

## if redshift smearing is required
if cfg['redshift_cat']['ztype'] == 'smearz' :
    print('Smearing redshifts', end='...')
    sigz_params = cfg['redshift_cat']['sigz']
    sigz_poly = lambda x : sigz_params[0] + x * sigz_params[1] + x**2 * sigz_params[2] + x**3 * sigz_params[3]
    noise = norm.rvs(0, sigz_poly(data['redshift']), len(data['redshift']))
    data['redshift'] += noise * ( 1 + data['redshift'])
    print('DONE')


## set non-detects to mag = 99.0
print('Replacing Non-Detections', end='...')
non_detects = {}
for band in bands :
    non_detects[band] = np.isnan(data[f"mag_{band}"]) | np.isinf(data[f"mag_{band}"])

    data[f"mag_{band}"][non_detects[band]] = 99.0
print('DONE')


## add colors
print('Adding colors', end='...')
for (band_,_band) in zip(bands[:-1],bands[1:]) :
    data[f"{band_}{_band}"] = np.array(data[f"mag_{band_}"] - data[f"mag_{_band}"])

    data[f"{band_}{_band}"][non_detects[band_] | non_detects[_band]] = 99.0
print('DONE')


## apply selection mask
mask = np.ones(len(data['id'])).astype(bool)
if True :
    print('Applying detection mask', end='...')
    ## mask if galaxy is fainter than 10-year depth in more than N bands
    mask_bands = []
    for band in bands :
        mask_bands.append(data[f"mag_{band}"] <= cfg['mag_limits']['survey_depth']['depths'][f"mag_{band}"])
    try :
        Npass = cfg['mag_limits']['survey_depth']['Npass']
    except :
        Npass = 3
    mask &= (sum([mask_band.astype(int) for mask_band in mask_bands]) >= Npass)

    if 'sample_cut' in cfg['mag_limits'].keys() :
        for key in cfg['mag_limits']['sample_cut'].keys() :
            mask &= (data[key] <= cfg['mag_limits']['sample_cut'][key])

    print('DONE')
if True :
    print('Applying fractured halo mask', end='...')
    ## mask if a cosmoDC2 fractured halo member
    badMembers = np.load('./lib/fractured_halos/badMembers.npy')
    mask &= np.isin(data['id'], badMembers, invert=True)
    print('DONE')
if 'err_threshold' in cfg['redshift_cat'] :
    print('Applying redshift error threshold mask', end='...')
    ## mask galaxies if their redshift_err > threshold * (1+z)
    mask &= (data['redshift_err'] < cfg['redshift_cat']['err_threshold'] * (1 + data['redshift']))
data = data[mask]


print(f'Getting relevant pixels in pixel.{pixel}', end='...')
pixels = hp.ang2pix(Nside_cov, data['ra'], data['dec'], nest=True, lonlat=True)
print('DONE')


galDir = os.path.join(cfg['outpath']['base'], cfg['name'], cfg['outpath']['galaxies'])
print(f'pixel.{pixel}:', end='\t')
for pix in np.unique(pixels) :
    print(f'{pix}', end=', ')
    
    inpix = (pixels == pix)
    data_inpix = data[inpix]
    
    if not data_inpix :
    	print(f'Error in {pix}')
    	continue
   
    if (Nside_cov < cfg['outpath']['hp_nside']) :
        pixFile = f'{galDir}{pix}_{pixel}.fits'
    else :
        pixFile = f'{galDir}{pix}.fits'
    data_inpix.write(pixFile, overwrite=True)
    
    print('DONE')
    del data_inpix
