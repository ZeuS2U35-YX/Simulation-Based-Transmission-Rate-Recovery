#!/bin/bash
#SBATCH --job-name=exp2_mif2
#SBATCH --array=1-9
#SBATCH --time=12:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/slurm-mif2-%A_%a.out
#SBATCH --error=logs/slurm-mif2-%A_%a.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:?Submit this script from the Experiment 2 directory}"

module --force purge
module load StdEnv/2020
module load r/4.1.2

export R_LIBS_USER="$HOME/packages-R4.1"

mkdir -p logs
mkdir -p results/array_output

Rscript code/02_Run_MIF2_Array.R "${SLURM_ARRAY_TASK_ID}"
