#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=30gb
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=crisw015@umn.edu
#SBATCH -p amdsmall
#SBATCH -e error.uvex
#SBATCH -o out.uvex
cd /home/vuk/crisw015
source .bashrc
cd /home/vuk/crisw015/UVEX/uvex-followup-etc-update/
conda activate uvenv-etc-update


## $1 is skymap directory
## $2 is the schedule output directory
## $3 is the downselected allsky CSV, with exposure times. Needs to have 'event_id' and 't_exp (ks)' columns.

SKYMAP_DIR=$1
SCHED_DIR=$2
#TEXP_FILE=$3

## field1 is event number, field2 is exposure time in ks
tail -n +2 $3 | while IFS=, read -r field1 field2; do
	num=$field1
	texp=$field2
	echo "Processing event ${num}..."
	echo "Exposure time is ${texp} ks."
	dorado-scheduling "${SKYMAP_DIR}/${num}.fits" -o "${SCHED_DIR}/${num}.ecsv" -j --mission=uvex --exptime="${texp} ks" --duration="3 hr" --roll-step='360 deg' --skygrid-method=sinusoidal --skygrid-step='10 deg2' --nside=128 --delay='10 yr'
done

conda deactivate

