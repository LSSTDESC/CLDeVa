#!/bin/bash
## LSST stack from cvmfs + the NERSC butler repo index, so Butler('dp2') resolves
## v30.0.8 is the release DP2 was processed with. needs ~/.lsst/db-auth.yaml for the registry
LSST_STACK=${LSST_STACK:-/cvmfs/sw.lsst.eu/almalinux-x86_64/lsst_distrib/v30.0.8}
source ${LSST_STACK}/loadLSST.bash
setup lsst_distrib
export DAF_BUTLER_REPOSITORY_INDEX=/global/cfs/cdirs/lsst/production/gen3/shared/data-repos.yaml
export OMP_NUM_THREADS=1
