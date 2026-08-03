#!/bin/bash
#SBATCH --job-name=simstudy200
#SBATCH --array=1-200%20
#SBATCH --time=06:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/slurm-simstudy200-%A_%a.out
#SBATCH --error=logs/slurm-simstudy200-%A_%a.err

set -euo pipefail

# Slurm normally starts the job in the directory from which sbatch was run.
# Submit this script from the Experiment 3 directory.
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

module --force purge
module load StdEnv/2020
module load r/4.1.2

export R_LIBS_USER="$HOME/packages-R4.1"

mkdir -p logs Results

Rscript code/02_run_hpc_task.R \
  "${SLURM_ARRAY_TASK_ID}" \
  "Results" \
  "results/paramlist.csv"
