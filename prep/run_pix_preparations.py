import os, sys, yaml
import shutil
import healsparse as hs
import healpy as hp
import numpy as np
from lib.utils import slurm_submit, run_argparse, get_pixels

args = run_argparse()

with open(args.cfgFile) as fstream :
    cfg = yaml.safe_load(fstream)

createFootprint = ('createFootprint' in args.processes)
createPix = ('createPix' in args.processes)
combinePix = False

galpath = os.path.join(cfg['outpath']['base'], cfg['name'], cfg['outpath']['galaxies'])
footpath = os.path.join(cfg['outpath']['base'], cfg['name'], cfg['outpath']['footprint'])


if createFootprint :
    if os.path.exists(footpath) :
        shutil.rmtree(footpath)
    os.makedirs(footpath)

    job_id = slurm_submit('createFootprints', configFile=args.cfgFile)
else :
    job_id = None


if createPix :
    if os.path.exists(galpath) :
        shutil.rmtree(galpath)
    os.makedirs(galpath)

    pixels = get_pixels(cfg['pixels'], cfg['outpath']['hp_nside'])

    job_ids = []
    for pixel in pixels :
        tmp = slurm_submit('createPixFilesParr_GCR', configFile=args.cfgFile, pixel=pixel, gap=False, dep=job_id)
        job_ids.append(tmp)
    job_id = ','.join(job_ids)
else :
    job_id = None


if combinePix :
    job_id = slurm_submit('combinePixFiles', configFile=args.cfgFile, dep=job_id)
else :
    job_id = None

