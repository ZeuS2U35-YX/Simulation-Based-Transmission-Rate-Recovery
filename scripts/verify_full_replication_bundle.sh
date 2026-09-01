#!/usr/bin/env bash

set -euo pipefail
export LC_ALL=C

usage() {
  cat <<'USAGE'
Usage:
  scripts/verify_full_replication_bundle.sh \
    [--skip-production-cardinality] PATH_TO_UNPACKED_PACKAGE

--skip-production-cardinality is only for deliberately incomplete test fixtures.
USAGE
}

die() {
  echo "$1" >&2
  exit "${2:-1}"
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
    die "Package contains a symbolic link or special file: $special_path"
  fi

  while IFS= read -r -d '' candidate; do
    relative_path="${candidate#"$root/"}"
    if contains_line_break "$relative_path"; then
      die "Package contains a path with CR or LF characters."
    fi
    if [[ "$relative_path" == *\\* ]]; then
      die "Package contains a path with a backslash: $relative_path"
    fi

    base_name="${relative_path##*/}"
    case "$base_name" in
      .git|.env|.env.*|.Renviron|.Renviron.*|.netrc|_netrc|.git-credentials|\
      .npmrc|.pypirc|id_rsa|id_rsa.*|id_ed25519|id_ed25519.*|\
      *.pem|*.key|*.p12|*.pfx|*.kdbx|credentials.json|credentials.*.json|\
      service-account*.json|service_account*.json)
        die "Package contains a credential-like file: $relative_path"
        ;;
    esac

    case "/$relative_path/" in
      */.git/*|*/.ssh/*|*/.aws/*|*/.gnupg/*)
        die "Package contains a credential-like directory: $relative_path"
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

skip_production_cardinality=0
package_root=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --skip-production-cardinality)
      skip_production_cardinality=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --*)
      die "Unknown option: $1" 2
      ;;
    *)
      if [[ -n "$package_root" ]]; then
        die "Only one unpacked package path may be supplied." 2
      fi
      package_root="$1"
      shift
      ;;
  esac
done

if [[ -z "$package_root" ]]; then
  usage >&2
  exit 2
fi

for required_command in awk cmp find mktemp rm shasum sort; do
  if ! command -v "$required_command" >/dev/null 2>&1; then
    die "Required command not found: $required_command"
  fi
done

if [[ -L "$package_root" || ! -d "$package_root" ]]; then
  die "Not an unpacked package directory: $package_root"
fi
package_root="$(cd "$package_root" && pwd -P)"

audit_package_tree "$package_root"

required_files=(
  "README_FIRST.md"
  "MANIFEST.csv"
  "PACKAGE_METADATA.txt"
  "SHA256SUMS"
  "experiments/experiment_5_bspline_B_recovery/results/comparison_three_models/three_model_overall_summary.csv"
)
required_directories=(
  "experiments/experiment_4_nmif600_model_comparison/shared_data"
  "experiments/experiment_4_nmif600_model_comparison/results_raw"
  "experiments/experiment_5_bspline_B_recovery/results_raw"
)

for relative_path in "${required_files[@]}"; do
  if [[ ! -f "$package_root/$relative_path" ]]; then
    die "Required package file is missing: $relative_path"
  fi
done
for relative_path in "${required_directories[@]}"; do
  if [[ ! -d "$package_root/$relative_path" ]]; then
    die "Required package directory is missing: $relative_path"
  fi
done

verification_tmp="$(mktemp -d "${TMPDIR:-/tmp}/verify-transmission-rate-recovery.XXXXXX")"
cleanup() {
  if [[ -n "${verification_tmp:-}" && -d "$verification_tmp" ]]; then
    rm -rf -- "$verification_tmp"
  fi
}
trap cleanup EXIT

actual_files="$verification_tmp/actual-files.txt"
listed_files="$verification_tmp/listed-files.txt"
listed_files_sorted="$verification_tmp/listed-files-sorted.txt"

(
  cd "$package_root"
  find . -type f ! -path './SHA256SUMS' | sort >"$actual_files"
)
: >"$listed_files"

while IFS= read -r checksum_line || [[ -n "$checksum_line" ]]; do
  checksum_hash="${checksum_line%% *}"
  if [[ ${#checksum_hash} -ne 64 || "$checksum_hash" == *[!0-9A-Fa-f]* ]]; then
    die "Malformed SHA256SUMS hash entry."
  fi

  checksum_prefix="${checksum_hash}  "
  case "$checksum_line" in
    "$checksum_prefix"*) checksum_path="${checksum_line#"$checksum_prefix"}" ;;
    *) die "Malformed SHA256SUMS line; expected two spaces before the path." ;;
  esac

  if contains_line_break "$checksum_path" || [[ "$checksum_path" == *\\* ]]; then
    die "Unsafe path encoding in SHA256SUMS."
  fi
  case "$checksum_path" in
    ./*) ;;
    *) die "SHA256SUMS path must be relative and begin with ./: $checksum_path" ;;
  esac

  checksum_relative="${checksum_path#./}"
  if [[ -z "$checksum_relative" ]]; then
    die "SHA256SUMS contains an empty path."
  fi
  case "/$checksum_relative/" in
    */../*|*/./*) die "SHA256SUMS contains a path traversal component: $checksum_path" ;;
  esac

  printf '%s\n' "$checksum_path" >>"$listed_files"
done <"$package_root/SHA256SUMS"

sort "$listed_files" >"$listed_files_sorted"
if ! cmp -s "$actual_files" "$listed_files_sorted"; then
  die "SHA256SUMS file list does not exactly match the package's regular files."
fi

(
  cd "$package_root"
  shasum -a 256 -c SHA256SUMS >/dev/null
)

if [[ "$skip_production_cardinality" -eq 1 ]]; then
  echo "WARNING: production cardinality checks were skipped for a test fixture." >&2
else
  check_production_cardinality "$package_root"
fi

echo "Replication bundle verification passed: $package_root"
