#!/bin/bash
#SBATCH --job-name=exp4_const600
# The submission wrapper supplies the array range explicitly.
#SBATCH --time=24:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/constant/slurm-%A_%a.out
#SBATCH --error=logs/constant/slurm-%A_%a.err

set -euo pipefail
: "${SLURM_ARRAY_TASK_ID:?Submit with an explicit --array range.}"
module --force purge
module load StdEnv/2020
module load r/4.1.2
export R_LIBS_USER="$HOME/packages-R4.1"
cd "$SLURM_SUBMIT_DIR"

Rscript code/03_run_constant_task.R \
  "$SLURM_ARRAY_TASK_ID" \
  shared_data \
  results_raw/constant \
  results/paramlist.csv
