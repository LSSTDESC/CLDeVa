import time
import subprocess
import os, sys, yaml
import shutil
import glob
import argparse as argparse
import numpy as np


def create_slurm_script(task, configFile, pixel):
    with open(configFile) as fstream:
        cfg = yaml.safe_load(fstream)

    print ('...task = ', task)

    scr = __file__.replace('utils.py', '')

    logPath = os.path.join(cfg['admin']['slurm']['logPath'].replace('TMP', scr+'../'), cfg['name'])
    scriptPath = os.path.join(cfg['admin']['slurm']['scriptPath'].replace('TMP', scr+'../'), cfg['name'])

    if pixel is None :
        logFile = os.path.join(logPath, cfg['admin']['slurm']['logFile'][task])
        script = os.path.join(scriptPath, cfg['admin']['slurm']['scriptFile'][task])
    else :
        logFile = os.path.join(logPath, f"{pixel}.out")
        script = os.path.join(scriptPath, f"{pixel}.out")

    if not os.path.exists(logPath) :
        os.makedirs(logPath)

    if not os.path.exists(scriptPath) :
        os.makedirs(scriptPath)

    f = open(f"{script}", "w")
    f.write("#!/bin/sh\n")
    f.write(f"#SBATCH --nodes={cfg['admin']['slurm']['Nnodes']}\n")
    f.write(f"#SBATCH --job-name={task}\n")
    f.write(f"#SBATCH --time={cfg['admin']['slurm']['time'][task]}\n")
    f.write(f"#SBATCH --partition=lsst,htc\n")
    f.write(f"#SBATCH --ntasks=2\n")
    f.write(f"#SBATCH --output={logFile}\n")
    f.write(f"#SBATCH --mem={cfg['admin']['slurm']['memory'][task]}GB\n")
    if pixel is None :
        f.write(f"python -u {scr}{task}.py {configFile}\n")
    else :
        f.write(f"python -u {scr}{task}.py {configFile} {pixel}\n")
    f.close()
    return script


def slurm_submit(task, configFile, pixel=None, dep=None, gap=True):

    if (dep is not None) and (gap == True) :
        time.sleep(3)

    script = create_slurm_script(task, configFile, pixel)
    if dep is not None:
        cmd = f"sbatch --depend=afterok:{dep} {script}"
    else:
        cmd = f"sbatch {script}"

    res = subprocess.run(cmd, shell=True, capture_output=True)
    job_id = str(res.stdout).split("Submitted batch job ")[1].split("\\")[0]
    return job_id


def run_argparse() :
    CLI = argparse.ArgumentParser()
    CLI.add_argument('cfgFile', type=str)
    CLI.add_argument('-p', '--processes', nargs='*', type=str, default=['createFootprint', 'createPix', 'combinePix'])

    return CLI.parse_args()


def get_pixels(pixels=None, nside=None) :
    if type(pixels) is str :
        if nside is None :
            raise ValueError('nside is required')
        pixFile = f"./lib/hpix_footprints/nside_{nside}/{pixels}.npy"
        pixels = np.load(pixFile)
    return pixels

