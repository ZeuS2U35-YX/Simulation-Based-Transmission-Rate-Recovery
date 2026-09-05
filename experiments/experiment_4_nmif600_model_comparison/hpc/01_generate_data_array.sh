#!/bin/bash
#SBATCH --job-name=exp4_data
# The submission wrapper supplies the array range explicitly.
#SBATCH --time=01:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/data/slurm-%A_%a.out
#SBATCH --error=logs/data/slurm-%A_%a.err

set -euo pipefail
: "${SLURM_ARRAY_TASK_ID:?Submit with an explicit --array range.}"
module --force purge
module load StdEnv/2020
module load r/4.1.2
export R_LIBS_USER="$HOME/packages-R4.1"
cd "$SLURM_SUBMIT_DIR"

Rscript code/01_generate_shared_data_task.R \
  "$SLURM_ARRAY_TASK_ID" \
  shared_data \
  results/paramlist.csv
