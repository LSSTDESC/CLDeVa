#!/bin/bash
## conda env with redmapper (no LSST stack in here)
RM_ENV=${RM_ENV:-desc_rm}
module load conda
conda activate ${RM_ENV}
export OMP_NUM_THREADS=1
