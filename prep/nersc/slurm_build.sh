#!/bin/bash
#SBATCH --account=m1727
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --job-name=dp2_build
#SBATCH --output=%x_%j.log

## galaxy catalog + spec-z for one setup, no star mask and each star mask
## sbatch slurm_build.sh griz 2
BANDS=$1
VIS=$2
OUT=${OUT:-$SCRATCH/dp2_redmapper/v2}
source ~/CLDeVa/prep/nersc/env_redmapper.sh
cd ~/CLDeVa/prep
DIR=$OUT/${BANDS}_vis${VIS}
for STAR in none 12 15 16; do
    if [ $STAR == none ]; then S=""; else S="--star $STAR"; fi
    echo "star mask: $STAR at $(date)"
    python -u create_redmapper_dp2_catalogs.py build --bands $BANDS --vis-min $VIS $S --map-dir $DIR/maps --cat-dir $DIR/catalogs --flux-dir $OUT/fluxes
    python -u create_redmapper_dp2_specz.py select --bands $BANDS --vis-min $VIS $S --map-dir $DIR/maps --flux-dir $OUT/fluxes --specz-dir $DIR/specz
    python -u create_redmapper_dp2_specz.py red --bands $BANDS --vis-min $VIS $S --map-dir $DIR/maps --flux-dir $OUT/fluxes --specz-dir $DIR/specz
done
