#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

while IFS= read -r -d '' script; do
  bash -n "$script"
done < <(find . -type f -name '*.sh' -print0)

check_submission_guard() {
  local script="$1"
  local confirmation_variable="$2"
  local status

  if env "${confirmation_variable}=" bash "$script" >/dev/null 2>&1; then
    printf 'Expected guarded submission wrapper to refuse: %s\n' "$script" >&2
    return 1
  else
    status=$?
  fi
  if [[ "$status" -ne 2 ]]; then
    printf 'Unexpected refusal status %s from %s\n' "$status" "$script" >&2
    return 1
  fi
  bash "$script" --dry-run >/dev/null
}

check_submission_guard \
  experiments/experiment_4_nmif600_model_comparison/hpc/submit_pilot.sh \
  EXP4_CONFIRM_PILOT_SUBMIT
check_submission_guard \
  experiments/experiment_4_nmif600_model_comparison/hpc/submit_all.sh \
  EXP4_CONFIRM_FULL_SUBMIT
check_submission_guard \
  experiments/experiment_5_bspline_B_recovery/hpc/submit_task1_pilot.sh \
  EXP5_CONFIRM_TASK1_PILOT
check_submission_guard \
  experiments/experiment_5_bspline_B_recovery/hpc/submit_all.sh \
  EXP5_CONFIRM_SUBMIT

LC_ALL=C Rscript --vanilla scripts/validate_release.R
