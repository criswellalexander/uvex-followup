#!/bin/bash

## runs the third part of the UVEX observing scenario simulation workflow
## this script will:
## - use all previously computed elements to compute the UVEX follow-up coverage for each event

# ALLSKY_FILE=$1
# SKYMAPS_DIR=$2
# OUTDIR=$3
# BAND=$4
# MAG=$5
# Rate adjustments
# Base is 0.7955 for O6 sim rate; 1 for O5 sim rate, assuming GWTC-2 measurement of 320/Mpc3/yr
# Switching to GWTC-3 median BNS rate of 210/Mpc3/yr takes an additional adjustment of 0.65625 (O6 w/ GWTC-3 rate overall adjustment is 0.52205)
# Switching to GWTC-3 PDB model median rate of 170/Mpc3/yr takes an adjustment of 0.53125 for O5, 0.42261 for O6
# RATE_ADJ=0.53125
# ALLSKY_SCHED="${OUTDIR}/allsky_sched_full.txt"
# SCHED_DIR="${OUTDIR}/schedules/"
# ALLSKY_COV="${OUTDIR}/allsky_coverage_full.txt"
echo "Computing UVEX coverage for all events..."

python3 /home/vuk/crisw015/UVEX/uvex-followup-etc-update//compute_tiling.py $PARAMS_FILE

echo "Getting statistics and making plots..."
python3 /home/vuk/crisw015/UVEX/uvex-followup-etc-update//make-coverage-plots.py $PARAMS_FILE

## package for easy download
# zip "${OUTDIR}/results.zip" "${OUTDIR}/uvex_*"

echo "Done!"
