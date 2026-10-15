#!/bin/bash
#SBATCH --account=m1727
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --job-name=dp2_stars
#SBATCH --output=%x_%j.log

## Gaia star masks with only the brighter stars (G < 12, 15, 16)
OUT=${OUT:-$SCRATCH/dp2_redmapper/v2}
source ~/CLDeVa/prep/nersc/env_butler.sh
cd ~/CLDeVa/prep
python -u make_star_masks.py --out-dir $OUT/star_masks
