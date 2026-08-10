#!/bin/bash
# Summarize Slurm resource use for the five-task pilot.

set -euo pipefail
cd "$(dirname "$0")/.."

if [[ ! -f results/pilot_job_ids.env ]]; then
  echo "Missing results/pilot_job_ids.env. Recreate it from the pilot submission output first." >&2
  exit 1
fi

# shellcheck disable=SC1091
source results/pilot_job_ids.env

required=(DATA_JOB_ID GAMMA_JOB_ID CONSTANT_JOB_ID POST_JOB_ID)
for name in "${required[@]}"; do
  if [[ -z "${!name:-}" ]]; then
    echo "Missing $name in results/pilot_job_ids.env" >&2
    exit 1
  fi
done

usage_file="results/pilot_resource_usage.psv"
summary_file="results/pilot_resource_summary.txt"

sacct -n -P \
  -j "${DATA_JOB_ID},${GAMMA_JOB_ID},${CONSTANT_JOB_ID},${POST_JOB_ID}" \
  --format=JobIDRaw,JobName,State,Elapsed,ElapsedRaw,MaxRSS,AllocCPUS,ExitCode \
  > "$usage_file"

max_elapsed_for_job() {
  local job_id="$1"
  awk -F'|' -v id="$job_id" '
    $1 ~ ("^" id "_[0-9]+$") && $5 ~ /^[0-9]+$/ {
      if ($5 > max) { max=$5; elapsed=$4; state=$3 }
    }
    END {
      if (max > 0) print elapsed; else print "NA"
    }
  ' "$usage_file"
}

count_state_for_job() {
  local job_id="$1"
  local state="$2"
  awk -F'|' -v id="$job_id" -v st="$state" '
    $1 ~ ("^" id "_[0-9]+$") && $3 == st { n++ }
    END { print n+0 }
  ' "$usage_file"
}

{
  echo "Experiment 4 pilot resource summary"
  echo "Generated: $(date)"
  echo
  echo "Gamma array job: $GAMMA_JOB_ID"
  echo "  completed tasks: $(count_state_for_job "$GAMMA_JOB_ID" COMPLETED)"
  echo "  longest elapsed: $(max_elapsed_for_job "$GAMMA_JOB_ID")"
  echo
  echo "Constant-B array job: $CONSTANT_JOB_ID"
  echo "  completed tasks: $(count_state_for_job "$CONSTANT_JOB_ID" COMPLETED)"
  echo "  longest elapsed: $(max_elapsed_for_job "$CONSTANT_JOB_ID")"
  echo
  echo "Raw Slurm accounting table: $usage_file"
  echo "Inspect MaxRSS, TIMEOUT, FAILED, OUT_OF_MEMORY, CANCELLED, or NODE_FAIL there before full submission."
} > "$summary_file"

cat "$summary_file"
