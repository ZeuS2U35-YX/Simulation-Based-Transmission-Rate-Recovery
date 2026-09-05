#!/bin/bash
#SBATCH --job-name=exp5_bspline600
# Task 1 is the separately validated production pilot.
# The submission wrapper supplies the array range explicitly.
#SBATCH --time=08:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/bspline/slurm-%A_%a.out
#SBATCH --error=logs/bspline/slurm-%A_%a.err

set -euo pipefail
: "${SLURM_ARRAY_TASK_ID:?Submit with an explicit --array range.}"
module --force purge
module load StdEnv/2020
module load r/4.1.2
export R_LIBS_USER="$HOME/packages-R4.1"
cd "$SLURM_SUBMIT_DIR"

Rscript --vanilla code/02_fit_bspline_B.R \
  "$SLURM_ARRAY_TASK_ID" \
  ../experiment_4_nmif600_model_comparison/shared_data \
  results/bspline \
  results/paired_input_manifest.csv
