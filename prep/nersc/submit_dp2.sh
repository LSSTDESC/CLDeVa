#!/bin/bash
## submits the butler job, then the redmapper build job to start only if it succeeds
##   ./submit_dp2.sh            full DP2
##   ./submit_dp2.sh --test     one tract (10188) on the debug queue, into $SCRATCH/dp2_redmapper_test
## OUT_DIR, N_PROCS, MAPS_ARGS, EXTRACT_ARGS, BUILD_ARGS are passed through to the jobs
set -e
NERSC_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export PREP_DIR=$(dirname ${NERSC_DIR})

BUTLER_OPTS=""
RM_OPTS=""
if [ "$1" == "--test" ]; then
    export OUT_DIR=${OUT_DIR:-$SCRATCH/dp2_redmapper_test}
    export MAPS_ARGS="--tracts 10188 ${MAPS_ARGS}"
    export EXTRACT_ARGS="--tracts 10188 ${EXTRACT_ARGS}"
    BUTLER_OPTS="--qos=debug --time=00:30:00"
    RM_OPTS="--qos=debug --time=00:15:00"
fi

cd ${NERSC_DIR}
butler_id=$(sbatch --parsable --export=ALL ${BUTLER_OPTS} slurm_dp2_butler.sh)
rm_id=$(sbatch --parsable --export=ALL --dependency=afterok:${butler_id} ${RM_OPTS} slurm_dp2_redmapper.sh)
echo "butler job ${butler_id} -> redmapper build job ${rm_id} (afterok), logs in ${NERSC_DIR}"
