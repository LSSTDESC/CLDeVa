#!/bin/bash
#SBATCH --job-name=dp2_rm_specz
#SBATCH --account=m1727
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --output=dp2_rm_specz_%j.log

## spec-z training sample (desc_rm conda env): DESI DR1 straight off CFS, cut to the geometry
## mask, then matched to the extracted tracts for the red-sequence spec file
## needs the maps and the tract flux files from slurm_dp2_butler.sh, no butler here
##   SELECT_ARGS / RED_ARGS pass extra flags, e.g. SELECT_ARGS="--training-hpix 398,430,367,461,462"
PREP_DIR=${PREP_DIR:?set PREP_DIR to CLDeVa/prep}
OUT_DIR=${OUT_DIR:-$SCRATCH/dp2_redmapper}

source ${PREP_DIR}/nersc/env_redmapper.sh
set -e
cd ${PREP_DIR}

echo "specz step on $(hostname) at $(date), reading ${OUT_DIR}"
DIRS="--map-dir ${OUT_DIR}/maps --cat-dir ${OUT_DIR}/catalogs --specz-dir ${OUT_DIR}/specz"
python -u create_redmapper_dp2_specz.py select ${DIRS} ${SELECT_ARGS}
python -u create_redmapper_dp2_specz.py red ${DIRS} ${RED_ARGS}
echo "finished at $(date)"
