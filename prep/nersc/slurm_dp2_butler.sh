#!/bin/bash
#SBATCH --job-name=dp2_rm_butler
#SBATCH --account=m1727
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --output=dp2_rm_butler_%j.log

## step 1 of 2, butler side (LSST stack): geometry mask + depth maps, then the per-tract flux files
## both steps skip tracts already on disk, so resubmitting picks up where it stopped
## submit with submit_dp2.sh, or set PREP_DIR yourself (sbatch copies this script, so it can't find itself)
##   MAPS_ARGS / EXTRACT_ARGS pass extra flags, e.g. EXTRACT_ARGS="--tracts 10188"
PREP_DIR=${PREP_DIR:?set PREP_DIR to CLDeVa/prep}
OUT_DIR=${OUT_DIR:-$SCRATCH/dp2_redmapper}
N_PROCS=${N_PROCS:-4}   ## butler bound, more workers mostly just hits the registry harder

source ${PREP_DIR}/nersc/env_butler.sh
set -e
cd ${PREP_DIR}

echo "butler step on $(hostname) at $(date), writing to ${OUT_DIR}"
python -u create_redmapper_dp2_maps.py --depth-bands min,r,i,z --n-procs ${N_PROCS} \
    --out-dir ${OUT_DIR}/maps ${MAPS_ARGS}
python -u create_redmapper_dp2_catalogs.py extract --n-procs ${N_PROCS} \
    --map-dir ${OUT_DIR}/maps --out-dir ${OUT_DIR}/catalogs ${EXTRACT_ARGS}
echo "finished at $(date)"
