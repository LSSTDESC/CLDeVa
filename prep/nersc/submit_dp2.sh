#!/bin/bash
## everything for riz 2 visits, riz 3 visits and griz 2 visits, all b < -20, each with no star mask
## and with the G < 12, 15, 16 star masks. writes to $SCRATCH/dp2_redmapper/v2
## the visit counts are done already, link them in first:  <setup>/maps/cell_counts -> the old riz / griz counts
cd ~/CLDeVa/prep/nersc

stars=$(sbatch --parsable slurm_star_masks.sh)
riz2=$(sbatch --parsable -d afterok:$stars slurm_maps.sh riz 2)
riz3=$(sbatch --parsable -d afterok:$stars slurm_maps.sh riz 3)
griz2=$(sbatch --parsable -d afterok:$stars slurm_maps.sh griz 2)
fluxes=$(sbatch --parsable -d afterok:$riz2 slurm_extract.sh)
sbatch -d afterok:$fluxes:$riz2 slurm_build.sh riz 2
sbatch -d afterok:$fluxes:$riz3 slurm_build.sh riz 3
sbatch -d afterok:$fluxes:$griz2 slurm_build.sh griz 2
