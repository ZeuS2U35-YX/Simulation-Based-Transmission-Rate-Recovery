#!/bin/bash
#SBATCH --job-name=exp4_post
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/postprocess/slurm-%j.out
#SBATCH --error=logs/postprocess/slurm-%j.err

set -euo pipefail
module --force purge
module load StdEnv/2020
module load r/4.1.2
export R_LIBS_USER="$HOME/packages-R4.1"
cd "$SLURM_SUBMIT_DIR"

Rscript code/04_combine_results.R \
  gamma results_raw/gamma results/combined/gamma 1:200

Rscript code/04_combine_results.R \
  constant results_raw/constant results/combined/constant 1:200

Rscript code/10_regenerate_filtering_mean_B_paths.R \
  results/combined/gamma/combined_B_filtering_means.csv \
  results/combined/gamma/filtering_mean_provenance.csv \
  4 1:200 shared_data

Rscript code/07_generate_task1_comparison_figures.R
Rscript code/08_generate_task117_comparison_figures.R

Rscript code/05_compare_models.R \
  results/combined/gamma \
  results/combined/constant \
  results/comparison \
  figures/comparison

Rscript code/06_make_convergence_diagnostics.R \
  results/combined/gamma \
  results/combined/constant \
  figures/convergence

bash hpc/package_outputs.sh
