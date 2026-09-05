#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

if [[ "${EXP5_CONFIRM_200_TASKS:-}" != "YES" ]]; then
  printf '%s\n' \
    'This runs the single 200-task Experiment 5 batch locally.' \
    'Set EXP5_CONFIRM_200_TASKS=YES to confirm, or use hpc/submit_all.sh.'
  exit 2
fi

Rscript --vanilla code/01_validate_paired_inputs.R

for task_id in $(seq 1 200); do
  Rscript --vanilla code/02_fit_bspline_B.R "$task_id"
done

Rscript --vanilla code/03_combine_results.R
Rscript --vanilla code/04_compare_with_gamma.R
Rscript --vanilla code/05_compare_three_models.R
