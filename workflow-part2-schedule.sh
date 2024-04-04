#!/bin/bash

## runs the second part of the UVEX observing scenario simulation workflow
## this script will:
## - load the outputs of followup-workflow-part1.sh (batched allsky files with t_exp calculated)
## - run texp_cut_and_batch.py to ensure minimum 500s exposures and format/batch for the scheduler
## - run sub_max_texp_calc.slurm on each batch to determine UVEX ToO schedules

SKYMAPS_DIR=$1
OUTDIR=$2
BAND=$3

## make directories
mkdir "${OUTDIR}/schedules/"

echo "Preprocessing..."
python3 /home/vuk/crisw015/UVEX/uvex-followup-etc-update/texp_cut_and_batch.py "${OUTDIR}/texp_out/" "${OUTDIR}" 40 40 $BAND

echo "t_exp processing, formatting, and rebatching complete. Beginning scheduler submission..."

for FILE in ${OUTDIR}/texp_sched/*
do
    echo "Submitting scheduler for batch file ${FILE}..."
    sbatch /home/vuk/crisw015/UVEX/uvex-followup-etc-update/sub_uvex_scheduler.slurm $SKYMAPS_DIR "${OUTDIR}/schedules/" $FILE
done

echo "Done! All jobs submitted."
