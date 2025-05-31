
from astropy.table import Table
from astropy.io import fits
from scipy.stats import norm
import numpy as np
import healpy as hp
import os, sys, yaml
import tables_io


with open(sys.argv[1]) as fstream :
	cfg = yaml.safe_load(fstream)

pixel = sys.argv[2]

galPath = cfg['gal_cat']['inpath']
zPath = cfg['redshift_cat']['inpath']

galFile = os.path.join(cfg['gal_cat']['inpath'], cfg['gal_cat']['fname'].replace('TMP', str(pixel)))
zFile = os.path.join(cfg['redshift_cat']['inpath'], cfg['redshift_cat']['fname'].replace('TMP', str(pixel)))

outDir = os.path.join(cfg['outpath']['base'], cfg['name'], cfg['outpath']['galaxies'])


names = ['id', 'ra', 'dec', 'mag_u', 'mag_g', 'mag_r', 'mag_i', 'mag_z', 'mag_y', 'redshift']

badMembers = np.load('./lib/badMembers.npy')

def getData(pixel, gdata, zdata) :
    pixcut = (pixels == pixel)
    
    ## match redshift file IDs assuming sorted via galaxy_id
    ## g_in_z are the galaxies in the pixel+cuts with redshift counterparts
    ## z_in_g are the redshifts with counterpart galaxies in the pixel+cuts
    g_in_z = np.isin(gdata[cfg['gal_cat']['keys']['id']][pixcut], zdata[cfg['redshift_cat']['keys']['id']],
            assume_unique=True)
    z_in_g = np.isin(zdata[cfg['redshift_cat']['keys']['id']], gdata[cfg['gal_cat']['keys']['id']][pixcut],
            assume_unique=True)
    
    ## remove BO masked regions
    if 'mask' in cfg.keys() :
        hp_mask = tables_io.read(cfg['mask']['path'])['PIXEL']
        nside_mask = cfg['mask']['hp_nside']
        g_pixels4mask = hp.ang2pix(nside_mask,
            gdata[cfg['gal_cat']['keys']['ra']][pixcut],
            gdata[cfg['gal_cat']['keys']['dec']][pixcut],
            nest=True, lonlat=True)
        g_in_z &= np.isin(g_pixels4mask, hp_mask, invert=True)
    z_in_g &= np.isin(zdata[cfg['redshift_cat']['keys']['id']], gdata[cfg['gal_cat']['keys']['id']][pixcut][g_in_z])
    
    
    ## remove members of fractured halos
    g_in_z &= np.isin(gdata[cfg['gal_cat']['keys']['id']][pixcut], badMembers, invert=True)
    z_in_g &= np.isin(zdata[cfg['redshift_cat']['keys']['id']], badMembers, invert=True)
    
    zcut = g_in_z
    if cfg['redshift_cat']['ztype'] == 'sigz' :
        redshift_npix = zdata[cfg['redshift_cat']['keys']['redshift']].values[z_in_g]
    elif cfg['redshift_cat']['ztype'] == 'pz' :
        #redshift_npix = zdata[cfg['redshift_cat']['keys']['redshift']].flatten()[z_in_g]
        redshift_npix = np.hstack(zdata[cfg['redshift_cat']['keys']['redshift']].values)[z_in_g]
    
    
    gal_id_npix = gdata[cfg['gal_cat']['keys']['id']][pixcut][zcut]
    ra_npix     = gdata[cfg['gal_cat']['keys']['ra']][pixcut][zcut]
    dec_npix    = gdata[cfg['gal_cat']['keys']['dec']][pixcut][zcut]
    mag_u_npix  = gdata[cfg['gal_cat']['keys']['mag_u']][pixcut][zcut]
    mag_g_npix  = gdata[cfg['gal_cat']['keys']['mag_g']][pixcut][zcut]
    mag_r_npix  = gdata[cfg['gal_cat']['keys']['mag_r']][pixcut][zcut]
    mag_i_npix  = gdata[cfg['gal_cat']['keys']['mag_i']][pixcut][zcut]
    mag_z_npix  = gdata[cfg['gal_cat']['keys']['mag_z']][pixcut][zcut]
    mag_y_npix  = gdata[cfg['gal_cat']['keys']['mag_y']][pixcut][zcut]
    
    
    gal_npix = Table([gal_id_npix, ra_npix, dec_npix, mag_u_npix, mag_g_npix, mag_r_npix, mag_i_npix, mag_z_npix, mag_y_npix, redshift_npix],
            names=names)
    
    return gal_npix



galcols = [cfg['gal_cat']['keys'][k] for k in cfg['gal_cat']['keys']]

if cfg['pixels'] == 'dc2' :
    galData = tables_io.read(galFile, columns=galcols, fmt=galFile.split('.')[-1], tType=tables_io.types.PD_DATAFRAME)['photometry']
else :
    galData = tables_io.read(galFile, columns=galcols, fmt=galFile.split('.')[-1])


if cfg['redshift_cat']['ztype'] == 'sigz' :
	zcols = [cfg['redshift_cat']['keys'][k] for k in cfg['redshift_cat']['keys']]
	zData = tables_io.read(zFile, columns=zcols, fmt=zFile.split('.')[-1])
	sigz = norm.rvs(0, cfg['redshift_cat']['sigz'], len(zData[cfg['redshift_cat']['keys']['redshift']]))
	zData[cfg['redshift_cat']['keys']['redshift']] += sigz * ( 1 + zData[cfg['redshift_cat']['keys']['redshift']])
elif cfg['redshift_cat']['ztype'] == 'pz' :
	zData = tables_io.read(zFile)


## apply cuts common among all runs and those particular to the current one
cut = np.ones(galData[cfg['gal_cat']['keys']['id']].shape, dtype=bool)

for pc in ['particular', 'common'] :
	for mag_col in cfg['cuts'][pc].keys() :
		if mag_col == 'color' :
			for key in cfg['cuts'][pc][mag_col].keys() :
				cpair = key.split('_')
				tmp = ( galData[cfg['gal_cat']['keys'][f"mag_{cpair[0]}"]] - galData[cfg['gal_cat']['keys'][f"mag_{cpair[1]}"]] )
				for arg in cfg['cuts'][pc][mag_col][key] :
					cut &= eval(arg)
		elif mag_col == 'mag' :
			mask_bands = []
			for band in cfg['cuts'][pc][mag_col].keys() :
				tmp = galData[cfg['gal_cat']['keys'][band]]
				tmpMask = np.ones_like(tmp).astype(bool)
				for arg in cfg['cuts'][pc][mag_col][band] :
					tmpMask &= eval(arg)
				mask_bands.append(tmpMask)
			if pc == 'particular' :
				cut &= np.array([np.all(g) for g in np.vstack(mask_bands).T])
			elif pc == 'common' :
				cut &= (sum([mb.astype(int) for mb in mask_bands]) >= 3)



galData = galData[cut]


## FIND WHICH PIXELS ARE IN THE SKYAREA
print(f'Getting relevant pixels in pixel.{pixel} ...', end='\t')
Nside = 64
pixels = hp.ang2pix(Nside, galData['ra'], galData['dec'], nest=True, lonlat=True)
print('DONE')


print(f'pixel.{pixel}:', end='\t')
for pix in np.unique(pixels) :
	print(f'{pix}', end=', ')

	gal = getData(pix, galData, zData)

	if not gal :
		print(f'#!{pix}!#')
		continue
	
	pixFile = f'{outDir}{pix}_{pixel}.fits'
	gal.write(pixFile, overwrite=True)

	print('DONE')
	del gal
