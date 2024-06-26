#!/bin/bash

## runs the first part of the UVEX observing scenario simulation workflow
## this script will:
## - create a new head directory [outdir] [3]
## - run localization_cut_and_batch.py to downselect and batch the skymaps in the desired allsky file [1] with associated skymaps in [2]
## - run sub_max_texp_calc.slurm on each batch to get exposure times for each localization and additionally downselect to <10ks exposures

## make sure path to params file is absolute
PARAMS_FILE="$(dirname $(readlink -e $1))/$(basename $1)"

## this grabs the save_directory line from the params file
OUTDIR=`grep "^save_directory="  ${PARAMS_FILE} | python3 -c "print(input().split('=')[1])"`
## make directories
mkdir $OUTDIR
mkdir "${OUTDIR}/texp_out"

## copy the params file to the output directory for future reference
cp $PARAMS_FILE "${OUTDIR}/params.ini"

## get absolute path to the uvex-followup directory
FOLLOWUP_DIR="$(dirname -- "$( readlink -f -- "$0"; )";)"

python3 ${FOLLOWUP_DIR}/localization_cut_and_batch.py $PARAMS_FILE

for BATCH_FILE in ${OUTDIR}/batches/*
do
    echo "Submitting exposure time calculation for batch file ${BATCH_FILE}..."
    sbatch ${FOLLOWUP_DIR}/sub_max_texp_calc.slurm $PARAMS_FILE $BATCH_FILE $FOLLOWUP_DIR
done

echo "Done! All jobs submitted."
