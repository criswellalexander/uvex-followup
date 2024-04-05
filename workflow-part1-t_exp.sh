#!/bin/bash

## runs the first part of the UVEX observing scenario simulation workflow
## this script will:
## - create a new head directory [outdir] [3]
## - run localization_cut_and_batch.py to downselect and batch the skymaps in the desired allsky file [1] with associated skymaps in [2]
## - run sub_max_texp_calc.slurm on each batch to get exposure times for each localization and additionally downselect to <10ks exposures

PARAMS_FILE=$1

## this grabs the save_diretory line from the params file
OUTDIR=`grep "^save_directory="  ${PARAMS_FILE} | python3 -c "print(input().split('=')[1])"`
## make directories
mkdir $OUTDIR
mkdir "${OUTDIR}/texp_out"

python3 /home/vuk/crisw015/UVEX/uvex-followup-etc-update/localization_cut_and_batch.py $PARAMS_FILE

for FILE in ${OUTDIR}/batches/*
do
    echo "Submitting exposure time calculation for batch file ${BATCH_FILE}..."
    sbatch /home/vuk/crisw015/UVEX/uvex-followup-etc-update/sub_max_texp_calc.slurm $PARAMS_FILE $BATCH_FILE
done

echo "Done! All jobs submitted."
