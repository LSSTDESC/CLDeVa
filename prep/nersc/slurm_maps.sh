#!/bin/bash
#SBATCH --account=m1727
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --job-name=dp2_maps
#SBATCH --output=%x_%j.log

## geometry + depth maps for one setup, no star mask and each star mask
## sbatch slurm_maps.sh griz 2
## the visit counts are already in $OUT/<bands>_vis<n>/maps/cell_counts (linked), so this only makes the maps
BANDS=$1
VIS=$2
OUT=${OUT:-$SCRATCH/dp2_redmapper/v2}
source ~/CLDeVa/prep/nersc/env_butler.sh
cd ~/CLDeVa/prep
python -u create_redmapper_dp2_maps.py --bands $BANDS --vis-min $VIS --n-procs 8 \
    --out-dir $OUT/${BANDS}_vis${VIS}/maps --star-dir $OUT/star_masks $TRACTS
