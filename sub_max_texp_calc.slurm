#!/bin/bash -l
#SBATCH --partition=cpu
#SBATCH --account=bcrv-delta-cpu
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=leggi014@umn.edu
#SBATCH -e error.uvex
#SBATCH -o out.uvex

#conda activate uvex-followup-env

## $1 is params file
## $2 is the downselected allsky CSV. At minimum needs to have 'simulation_id', 'area_(90)', and 'distmean' (or 'dist(90)') columns.
## $3 is the uvex-followup directory

## slurm doesn't actually know where the script lives
FOLLOWUP_DIR=$3

python3 -u ${FOLLOWUP_DIR}/max-texp-by-sky-loc.py $1 $2

#conda deactivate

