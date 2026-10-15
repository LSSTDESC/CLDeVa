#!/bin/bash
#SBATCH --account=m1727
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --job-name=dp2_extract
#SBATCH --output=%x_%j.log

## griz flux files for every setup, cut on the biggest footprint (riz, 2 visits, no star mask)
OUT=${OUT:-$SCRATCH/dp2_redmapper/v2}
source ~/CLDeVa/prep/nersc/env_butler.sh
cd ~/CLDeVa/prep
python -u create_redmapper_dp2_catalogs.py extract --n-procs 8 --flux-dir $OUT/fluxes \
    --cells $OUT/riz_vis2/maps/dp2_geometry_riz_vis2_bsouth20_nostar_nside131072_bitpacked.hs $TRACTS
