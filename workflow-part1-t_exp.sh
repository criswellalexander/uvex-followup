#!/bin/bash

## runs the first part of the UVEX observing scenario simulation workflow
## this script will:
## - create a new head directory [outdir] [3]
## - run localization_cut_and_batch.py to downselect and batch the skymaps in the desired allsky file [1] with associated skymaps in [2]
## - run sub_max_texp_calc.slurm on each batch to get exposure times for each localization and additionally downselect to <10ks exposures

ALLSKY_FILE=$1
SKYMAPS_DIR=$2
OUTDIR=$3
BAND=$4
MAG_AB=$5
AREA_CUT=$6

## make directories
mkdir $OUTDIR
mkdir "${OUTDIR}/texp_out"

python3 /home/vuk/crisw015/UVEX/uvex-followup-etc-update/localization_cut_and_batch.py "${ALLSKY_FILE}" "${OUTDIR}" --N_batch 40 --max_area $AREA_CUT

for FILE in ${OUTDIR}/batches/*
do
    echo "Submitting exposure time calculation for batch file ${FILE}..."
    sbatch /home/vuk/crisw015/UVEX/uvex-followup-etc-update/sub_max_texp_calc.slurm $SKYMAPS_DIR "${OUTDIR}/texp_out/" $FILE $BAND $MAG_AB
done

echo "Done! All jobs submitted."
