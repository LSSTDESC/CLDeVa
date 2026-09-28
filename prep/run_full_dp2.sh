#!/bin/bash
## full DP2 inputs: geometry + depth, then the photometric catalog
## shear-group visit-based geometry (weight_cut 0.0, 2x2 block, vis_min 3)
## both slow steps are butler bound, so 4 workers each
set -e
cd ~/redmapper

python -u create_redmapper_dp2_maps.py --depth-bands min,r,i,z --n-procs 4
python -u create_redmapper_dp2_catalogs.py extract --n-procs 4
python -u create_redmapper_dp2_catalogs.py build
