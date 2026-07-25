#!/bin/bash

set -euo pipefail

cd /global/home/hpc6245/experiment2_large_scale_gamma_B_recovery_HPC

module --force purge
module load StdEnv/2020
module load r/4.1.2

export R_LIBS_USER="$HOME/packages-R4.1"

mkdir -p data results/array_output figures logs notes

echo "Generating the fixed piecewise-B dataset..."
Rscript code/01_Generate_Fixed_Piecewise_Data.R

echo "Submitting the nine MIF2 array tasks..."
ARRAY_OUTPUT=$(sbatch hpc/01_submit_mif2_array.sh)
ARRAY_JOB_ID=$(echo "$ARRAY_OUTPUT" | awk '{print $4}')

if [ -z "$ARRAY_JOB_ID" ]; then
  echo "Could not parse the array job ID."
  exit 1
fi

echo "Array job ID: $ARRAY_JOB_ID"
echo "Submitting combine-and-plot job after successful completion of all array tasks..."

FINISH_OUTPUT=$(sbatch --dependency=afterok:${ARRAY_JOB_ID} hpc/02_combine_and_plot.sh)
FINISH_JOB_ID=$(echo "$FINISH_OUTPUT" | awk '{print $4}')

echo "Finish job ID: $FINISH_JOB_ID"
echo
echo "Experiment submitted."
echo "Check status with:"
echo "  squeue -u $USER"
echo
echo "Results will be written under:"
echo "  results/"
echo "Figures will be written under:"
echo "  figures/"
echo "Logs will be written under:"
echo "  logs/"
