#!/usr/bin/env bash

set -euo pipefail
export LC_ALL=C

usage() {
  cat <<'USAGE'
Build the full replication ZIP while preserving the repository directory layout.

Usage:
  scripts/build_full_replication_bundle.sh \
    --source-checkout PATH \
    --data-source PATH \
    [--output-dir PATH] \
    [--version VERSION] \
    [--source-ref REF] \
    [--source-commit SHA] \
    [--skip-production-cardinality]

VERSION must have the form vX.Y.Z or vX.Y.Z-suffix. The source checkout may
be a Git working tree or a GitHub source archive. For a Git working tree, only
files tracked at HEAD are exported. The data source must contain the ignored
full-data directories at their original repository-relative paths.

--skip-production-cardinality is only for small test fixtures. Never use it to
build a release artifact. Existing ZIPs are never overwritten.
USAGE
}

die() {
  echo "$1" >&2
  exit "${2:-1}"
}

require_option_value() {
  if [[ $# -lt 2 || "$2" == --* ]]; then
    die "Missing value for option: $1" 2
  fi
}

contains_line_break() {
  [[ "$1" == *$'\n'* || "$1" == *$'\r'* ]]
}

audit_package_tree() {
  local root="$1"
  local special_path=""
  local candidate=""
  local relative_path=""
  local base_name=""

  special_path="$(find "$root" ! -type f ! -type d -print -quit)"
  if [[ -n "$special_path" ]]; then
    die "Refusing to package a symbolic link or special file: $special_path"
  fi

  while IFS= read -r -d '' candidate; do
    relative_path="${candidate#"$root/"}"
    if contains_line_break "$relative_path"; then
      die "Refusing to package a path containing CR or LF characters."
    fi
    if [[ "$relative_path" == *\\* ]]; then
      die "Refusing to package a path containing a backslash: $relative_path"
    fi

    base_name="${relative_path##*/}"
    case "$base_name" in
      .git|.env|.env.*|.Renviron|.Renviron.*|.netrc|_netrc|.git-credentials|\
      .npmrc|.pypirc|id_rsa|id_rsa.*|id_ed25519|id_ed25519.*|\
      *.pem|*.key|*.p12|*.pfx|*.kdbx|credentials.json|credentials.*.json|\
      service-account*.json|service_account*.json)
        die "Refusing to package a credential-like file: $relative_path"
        ;;
    esac

    case "/$relative_path/" in
      */.git/*|*/.ssh/*|*/.aws/*|*/.gnupg/*)
        die "Refusing to package a credential-like directory: $relative_path"
        ;;
    esac
  done < <(find "$root" -mindepth 1 -print0)
}

check_production_cardinality() {
  local root="$1"
  local group_paths=(
    "experiments/experiment_4_nmif600_model_comparison/shared_data"
    "experiments/experiment_4_nmif600_model_comparison/results_raw/gamma"
    "experiments/experiment_4_nmif600_model_comparison/results_raw/constant"
    "experiments/experiment_5_bspline_B_recovery/results_raw/bspline"
  )
  local group_labels=(
    "Experiment 4 shared inputs"
    "Experiment 4 Gamma-noise results"
    "Experiment 4 constant-B results"
    "Experiment 5 B-spline results"
  )
  local group_index=0
  local task_number=0
  local task_name=""
  local task_path=""
  local summary_path=""

  for group_index in 0 1 2 3; do
    for ((task_number=1; task_number<=200; task_number++)); do
      printf -v task_name 'task_%03d' "$task_number"
      task_path="$root/${group_paths[$group_index]}/$task_name"
      if [[ ! -d "$task_path" || ! -f "$task_path/COMPLETE" ]]; then
        die "Production cardinality check failed for ${group_labels[$group_index]}: missing $task_name/COMPLETE"
      fi
    done
  done

  summary_path="$root/experiments/experiment_5_bspline_B_recovery/results/comparison_three_models/three_model_overall_summary.csv"
  if [[ ! -f "$summary_path" ]]; then
    die "Production cardinality check failed: missing three-model summary."
  fi

  if ! awk -F, '
    function unquote(value) {
      sub(/\r$/, "", value)
      if (value ~ /^".*"$/) {
        value = substr(value, 2, length(value) - 2)
      }
      return value
    }
    NR == 1 { next }
    NF > 1 {
      model = unquote($1)
      n_tasks = unquote($3)
      rows++
      if (n_tasks == "200") {
        seen[model]++
      }
    }
    END {
      if (rows != 3 ||
          seen["gamma_noise"] != 1 ||
          seen["bspline_B"] != 1 ||
          seen["constant_B"] != 1) {
        exit 1
      }
    }
  ' "$summary_path"; then
    die "Production cardinality check failed: the summary must contain exactly gamma_noise, bspline_B, and constant_B with n_tasks=200."
  fi
}

source_checkout=""
data_source=""
output_dir="dist"
version="v1.1.0"
source_ref=""
source_commit=""
skip_production_cardinality=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --source-checkout)
      require_option_value "$@"
      source_checkout="$2"
      shift 2
      ;;
    --data-source)
      require_option_value "$@"
      data_source="$2"
      shift 2
      ;;
    --output-dir)
      require_option_value "$@"
      output_dir="$2"
      shift 2
      ;;
    --version)
      require_option_value "$@"
      version="$2"
      shift 2
      ;;
    --source-ref)
      require_option_value "$@"
      source_ref="$2"
      shift 2
      ;;
    --source-commit)
      require_option_value "$@"
      source_commit="$2"
      shift 2
      ;;
    --skip-production-cardinality)
      skip_production_cardinality=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -z "$source_checkout" || -z "$data_source" ]]; then
  usage >&2
  exit 2
fi

if [[ ! "$version" =~ ^v[0-9]+\.[0-9]+\.[0-9]+(-[A-Za-z0-9][A-Za-z0-9._-]*)?$ ]]; then
  die "Invalid version string: $version (expected vX.Y.Z or vX.Y.Z-suffix)." 2
fi

if contains_line_break "$source_ref"; then
  die "--source-ref must not contain CR or LF characters." 2
fi
if contains_line_break "$source_commit"; then
  die "--source-commit must not contain CR or LF characters." 2
fi

for required_command in awk cat date find git mkdir mktemp mv rm rsync shasum sort tar tr unzip wc zip; do
  if ! command -v "$required_command" >/dev/null 2>&1; then
    die "Required command not found: $required_command"
  fi
done

if [[ ! -f "$source_checkout/README.md" ]]; then
  die "Not a repository checkout or source archive: $source_checkout"
fi

required_full_data=(
  "experiments/experiment_4_nmif600_model_comparison/shared_data"
  "experiments/experiment_4_nmif600_model_comparison/results_raw"
  "experiments/experiment_5_bspline_B_recovery/results_raw"
)

for relative_path in "${required_full_data[@]}"; do
  if [[ ! -d "$data_source/$relative_path" ]]; then
    die "Missing required full-data directory: $data_source/$relative_path"
  fi
  unexpected_data_path="$(find "$data_source/$relative_path" ! -type f ! -type d -print -quit)"
  if [[ -n "$unexpected_data_path" ]]; then
    die "Refusing to copy a symbolic link or special data file: $unexpected_data_path"
  fi
done

source_checkout="$(cd "$source_checkout" && pwd -P)"
data_source="$(cd "$data_source" && pwd -P)"
mkdir -p "$output_dir"
output_dir="$(cd "$output_dir" && pwd -P)"

is_git_checkout=0
if [[ -e "$source_checkout/.git" || -L "$source_checkout/.git" ]]; then
  if [[ -L "$source_checkout/.git" ]]; then
    die "Refusing a source checkout whose .git entry is a symbolic link."
  fi
  git_top_level="$(git -C "$source_checkout" rev-parse --show-toplevel 2>/dev/null || true)"
  if [[ -z "$git_top_level" ]]; then
    die "The source has a .git entry but is not a valid Git working tree."
  fi
  git_top_level="$(cd "$git_top_level" && pwd -P)"
  if [[ "$git_top_level" != "$source_checkout" ]]; then
    die "--source-checkout must point to the Git working-tree root."
  fi
  is_git_checkout=1
  if [[ -z "$source_commit" ]]; then
    source_commit="$(git -C "$source_checkout" rev-parse HEAD)"
  fi
  if [[ -z "$source_ref" ]]; then
    source_ref="$(git -C "$source_checkout" symbolic-ref --short -q HEAD 2>/dev/null || true)"
  fi
fi

if contains_line_break "$source_ref" || contains_line_break "$source_commit"; then
  die "Resolved source metadata contains CR or LF characters."
fi

package_root_name="Simulation-Based-Transmission-Rate-Recovery-${version}-full-replication"
zip_name="${package_root_name}.zip"
final_zip="$output_dir/$zip_name"
final_checksum="$output_dir/$zip_name.sha256"

if [[ -e "$final_zip" || -L "$final_zip" || -e "$final_checksum" || -L "$final_checksum" ]]; then
  die "Refusing to overwrite an existing release artifact: $final_zip or $final_checksum"
fi

staging_parent=""
output_temp_dir=""
cleanup() {
  if [[ -n "$staging_parent" && -d "$staging_parent" ]]; then
    rm -rf -- "$staging_parent"
  fi
  if [[ -n "$output_temp_dir" && -d "$output_temp_dir" ]]; then
    rm -rf -- "$output_temp_dir"
  fi
}
trap cleanup EXIT

staging_parent="$(mktemp -d "${TMPDIR:-/tmp}/transmission-rate-recovery.XXXXXX")"
package_root="$staging_parent/$package_root_name"
tracked_source="$staging_parent/tracked-source"
mkdir -p "$package_root"

if [[ "$is_git_checkout" -eq 1 ]]; then
  mkdir -p "$tracked_source"
  git -C "$source_checkout" archive --format=tar HEAD \
    | tar -xf - -C "$tracked_source"
else
  tracked_source="$source_checkout"
fi

rsync -a \
  --exclude '/.git' \
  --exclude '/experiments/experiment_4_nmif600_model_comparison/shared_data/' \
  --exclude '/experiments/experiment_4_nmif600_model_comparison/results_raw/' \
  --exclude '/experiments/experiment_5_bspline_B_recovery/results_raw/' \
  --exclude '.DS_Store' \
  --exclude '__MACOSX/' \
  --exclude '.Rhistory' \
  --exclude '.RData' \
  --exclude '.Rproj.user/' \
  --exclude '.env' \
  --exclude '.env.*' \
  --exclude '.Renviron' \
  --exclude '.Renviron.*' \
  --exclude '*.pem' \
  --exclude '*.key' \
  --exclude '__pycache__/' \
  --exclude '*.pyc' \
  --exclude '*.out' \
  --exclude '*.err' \
  --exclude '*~' \
  --exclude '.#*' \
  --exclude 'renv/library/' \
  --exclude 'renv/staging/' \
  --exclude 'logs/' \
  --exclude 'downloads/' \
  --exclude 'tmp/' \
  --exclude 'experiments/experiment_4_nmif600_model_comparison/results/pilot_job_ids.env' \
  --exclude 'experiments/experiment_4_nmif600_model_comparison/results/pilot_submission.txt' \
  --exclude 'experiments/experiment_4_nmif600_model_comparison/results/full_submission.txt' \
  --exclude '*.zip' \
  --exclude '*.tar.gz' \
  "$tracked_source/" "$package_root/"

for relative_path in "${required_full_data[@]}"; do
  mkdir -p "$package_root/$relative_path"
  rsync -a \
    --exclude '.DS_Store' \
    --exclude '*.tmp' \
    --exclude '*.log' \
    "$data_source/$relative_path/" "$package_root/$relative_path/"
done

audit_package_tree "$package_root"

if [[ "$skip_production_cardinality" -eq 1 ]]; then
  cardinality_status="skipped-for-test-fixture"
  echo "WARNING: production cardinality checks were skipped for a test fixture." >&2
else
  check_production_cardinality "$package_root"
  cardinality_status="passed"
fi

cat >"$package_root/PACKAGE_METADATA.txt" <<METADATA
package_name=$package_root_name
package_version=$version
source_repository=https://github.com/ZeuS2U35-YX/Simulation-Based-Transmission-Rate-Recovery
source_ref=${source_ref:-not-recorded}
source_commit=${source_commit:-not-recorded}
created_utc=$(date -u '+%Y-%m-%dT%H:%M:%SZ')
layout=single-top-level-directory; repository-relative paths preserved
production_cardinality_check=$cardinality_status
METADATA

manifest_path="$package_root/MANIFEST.csv"
printf 'path,size_bytes,sha256,role,provenance\n' >"$manifest_path"

classify_role() {
  case "$1" in
    experiments/*/shared_data/*) echo "accepted simulated input" ;;
    experiments/*/results_raw/*) echo "raw task-level result" ;;
    experiments/*/results/*) echo "validated or post-processed result" ;;
    experiments/*/figures/*) echo "figure or figure source data" ;;
    experiments/*/code/*|shared_code/*|scripts/*) echo "analysis or workflow code" ;;
    experiments/*/hpc/*) echo "HPC submission workflow" ;;
    */renv.lock|*/renv/activate.R|*/.Rprofile) echo "software environment" ;;
    *.md|*.txt|CITATION.cff|LICENSE) echo "documentation or metadata" ;;
    *) echo "repository file" ;;
  esac
}

classify_provenance() {
  case "$1" in
    PACKAGE_METADATA.txt)
      echo "package build"
      ;;
    experiments/experiment_4_nmif600_model_comparison/shared_data/*|\
    experiments/experiment_4_nmif600_model_comparison/results_raw/*|\
    experiments/experiment_5_bspline_B_recovery/results_raw/*)
      echo "author archival workspace"
      ;;
    *)
      echo "Git source checkout"
      ;;
  esac
}

while IFS= read -r file_path; do
  relative_path="${file_path#"$package_root/"}"
  [[ "$relative_path" == "MANIFEST.csv" || "$relative_path" == "SHA256SUMS" ]] && continue
  file_size="$(wc -c <"$file_path" | tr -d ' ')"
  file_hash="$(shasum -a 256 "$file_path" | awk '{print $1}')"
  escaped_path="${relative_path//\"/\"\"}"
  file_role="$(classify_role "$relative_path")"
  file_provenance="$(classify_provenance "$relative_path")"
  printf '"%s",%s,%s,"%s","%s"\n' \
    "$escaped_path" "$file_size" "$file_hash" "$file_role" "$file_provenance" \
    >>"$manifest_path"
done < <(find "$package_root" -type f | sort)

(
  cd "$package_root"
  : >SHA256SUMS
  while IFS= read -r relative_path; do
    file_hash="$(shasum -a 256 "$relative_path" | awk '{print $1}')"
    printf '%s  %s\n' "$file_hash" "$relative_path" >>SHA256SUMS
  done < <(find . -type f ! -path './SHA256SUMS' | sort)
  shasum -a 256 -c SHA256SUMS >/dev/null
)

output_temp_dir="$(mktemp -d "$output_dir/.full-replication-output.XXXXXX")"
temporary_zip="$output_temp_dir/$zip_name"
temporary_checksum="$output_temp_dir/$zip_name.sha256"

(
  cd "$staging_parent"
  zip -q -r "$temporary_zip" "$package_root_name"
)
zip -T "$temporary_zip" >/dev/null

archive_hash="$(shasum -a 256 "$temporary_zip" | awk '{print $1}')"
archive_size="$(wc -c <"$temporary_zip" | tr -d ' ')"
printf '%s  %s\n' "$archive_hash" "$zip_name" >"$temporary_checksum"

mv "$temporary_zip" "$final_zip"
mv "$temporary_checksum" "$final_checksum"

echo "Created: $final_zip"
echo "Checksum file: $final_checksum"
echo "Size (bytes): $archive_size"
echo "SHA-256: $archive_hash"
