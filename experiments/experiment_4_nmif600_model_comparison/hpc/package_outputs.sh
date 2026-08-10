#!/bin/bash
# Create one archive that can be downloaded from the HPC after post-processing.

set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p downloads

archive="downloads/experiment_4_nmif600_model_comparison_outputs.tar.gz"
manifest="downloads/experiment_4_nmif600_model_comparison_manifest.txt"

find README.md .gitignore config code hpc \
  shared_data results_raw results figures logs \
  -type f \
  ! -path '*/.gitkeep' \
  ! -name '*.tmp*' \
  | sort > "$manifest"

tar -czf "$archive" \
  --exclude='*/.gitkeep' \
  --exclude='*.tmp*' \
  README.md .gitignore config code hpc \
  shared_data results_raw results figures logs

printf 'Created %s\nCreated %s\n' "$archive" "$manifest"
