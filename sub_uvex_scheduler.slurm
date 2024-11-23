#!/bin/bash -l
#SBATCH --partition=cpu
#SBATCH --account=bcrv-delta-cpu
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --mem=30gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=leggi014@umn.edu
#SBATCH -e error.uvex
#SBATCH -o out.uvex

#conda activate uvex-followup-env

## $1 is skymap directory
## $2 is the schedule output directory
## $3 is the downselected allsky CSV, with exposure times. Needs to have 'event_id' and 't_exp (ks)' columns.

PARAMS_FILE=$1
#TEXP_FILE=$2

## extract relevant lines
OUTDIR=`grep "^save_directory="  ${PARAMS_FILE} | python3 -c "print(input().split('=')[1])"`
SKYMAP_DIR=`grep "^obs_scenario="  ${PARAMS_FILE} | python3 -c "print(input().split('=')[1]+'/allsky/')"`
DURATION=`grep "^tiling_time="  ${PARAMS_FILE} | python3 -c "print(input().split('=')[1])"`

SCHED_DIR="${OUTDIR}/schedules/"


## field1 is event number, field2 is exposure time in ks
tail -n +2 $2 | while IFS=, read -r field1 field2; do
    echo "duration is ${DURATION}"
	num=$field1
	texp=$field2
	echo "Processing event ${num}..."
	echo "Exposure time is ${texp} ks."
	m4opt schedule "${SKYMAP_DIR}/${num}.fits"  "${SCHED_DIR}/${num}.ecsv" --mission=uvex --exptime="${texp} ks"  --deadline="24hour" --nside=128 --timelimit="20min" --deadline="24hour"
done

#conda deactivate

