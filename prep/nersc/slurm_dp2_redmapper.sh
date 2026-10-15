#!/bin/bash
#SBATCH --job-name=dp2_rm_build
#SBATCH --account=m1727
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --output=dp2_rm_build_%j.log

## step 2 of 2, redmapper side (desc_rm conda env): tract flux files -> redmapper galaxy catalog
## no butler here, it only reads what slurm_dp2_butler.sh wrote
##   BUILD_ARGS passes extra flags, e.g. BUILD_ARGS="--mag-cut 23.5"
PREP_DIR=${PREP_DIR:?set PREP_DIR to CLDeVa/prep}
OUT_DIR=${OUT_DIR:-$SCRATCH/dp2_redmapper}

source ${PREP_DIR}/nersc/env_redmapper.sh
set -e
cd ${PREP_DIR}

echo "build step on $(hostname) at $(date), reading ${OUT_DIR}"
python -u create_redmapper_dp2_catalogs.py build \
    --map-dir ${OUT_DIR}/maps --out-dir ${OUT_DIR}/catalogs ${BUILD_ARGS}
echo "finished at $(date)"
