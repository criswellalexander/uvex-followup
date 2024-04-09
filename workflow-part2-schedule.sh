#!/bin/bash

## runs the second part of the UVEX observing scenario simulation workflow
## this script will:
## - load the outputs of followup-workflow-part1.sh (batched allsky files with t_exp calculated)
## - run texp_cut_and_batch.py to ensure minimum 500s exposures and format/batch for the scheduler
## - run sub_max_texp_calc.slurm on each batch to determine UVEX ToO schedules

# SKYMAPS_DIR=$1
# OUTDIR=$2
# BAND=$3

PARAMS_FILE="$(dirname $(readlink -e $1))/$(basename $1)"

## this grabs the save_directory line from the params file
OUTDIR=`grep "^save_directory="  ${PARAMS_FILE} | python3 -c "print(input().split('=')[1])"`

## make directories
mkdir "${OUTDIR}/schedules/"

## get absolute path to the uvex-followup directory
FOLLOWUP_DIR=`dirname -- "$( readlink -f -- "$0"; )";`

echo "Preprocessing..."
python3 ${FOLLOWUP_DIR}/texp_cut_and_batch.py $PARAMS_FILE

echo "t_exp processing, formatting, and rebatching complete. Beginning scheduler submission..."

for FILE in ${OUTDIR}/texp_sched/*
do
    echo "Submitting scheduler for batch file ${FILE}..."
    sbatch ${FOLLOWUP_DIR}/sub_uvex_scheduler.slurm $PARAMS_FILE $FILE
done

echo "Done! All jobs submitted."
