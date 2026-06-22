#!/bin/bash
#SBATCH --job-name=rm_cal
#SBATCH --account=m1727
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --output=rm_cal_%j.log

# Environment
module load conda
conda activate desc_rm

export GCR_CONFIG_SOURCE="files"

cd /global/u1/d/dericj/cdc2_rmpr

# Run the calibration
echo "Starting redMaPPer calibration at $(date)"
echo "Host: $(hostname)"

redmapper_calibrate.py -c cal.yml

echo "Finished at $(date)"
