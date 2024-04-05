#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=crisw015@umn.edu
#SBATCH -p amdsmall
#SBATCH -e error.uvex
#SBATCH -o out.uvex
cd /home/vuk/crisw015
source .bashrc
cd /home/vuk/crisw015/UVEX/uvex-followup/
conda activate uvenv-etc-update

## $1 is params file
## $2 is the downselected allsky CSV. At minimum needs to have 'simulation_id', 'area_(90)', and 'distmean' (or 'dist(90)') columns.

python3 -u /home/vuk/crisw015/UVEX/uvex-followup/max-texp-by-sky-loc.py $1 $2

conda deactivate

