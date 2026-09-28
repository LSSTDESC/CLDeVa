#!/bin/bash
## finish the catalog after the maps: extract the tracts that are still missing,
## then build. extract now survives a bad tract, so rerun until nothing fails
set -e
cd ~/redmapper

python -u create_redmapper_dp2_catalogs.py extract --n-procs 4
python -u create_redmapper_dp2_catalogs.py build
