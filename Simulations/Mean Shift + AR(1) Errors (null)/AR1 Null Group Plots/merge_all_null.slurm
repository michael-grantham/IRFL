#!/bin/bash
#SBATCH --job-name=run_merge_all_null
#SBATCH --output=run_merge_all_null.%j.out
#SBATCH --error=run_merge_all_null.%j.err
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

set -euo pipefail

# --- Submit from the "AR1 Group Plots Null" directory ---
ROOT_DIR="${SLURM_SUBMIT_DIR}"
PARENT_DIR="$(dirname "$ROOT_DIR")"

echo "Submit dir           : $ROOT_DIR"
echo "Parent (one level up): $PARENT_DIR"
echo

# --- Null-case folders to process ---
FOLDERS=(
  "mean_shift_ar1_null_0.0"
  "mean_shift_ar1_null_0.3"
  "mean_shift_ar1_null_0.7"
  "mean_shift_ar1_null_neg0.3"
  "mean_shift_ar1_null_neg0.7"
)

declare -a SUBMITTED_JOBS=()
declare -a SKIPPED_FOLDERS=()

submit_merge_job () {
  local target_dir="$1"
  local script_path="$target_dir/merge.slurm"

  if [[ ! -d "$target_dir" ]]; then
    echo "WARN: Missing directory: $target_dir"
    SKIPPED_FOLDERS+=("$(basename "$target_dir") (no dir)")
    return 0
  fi
  if [[ ! -f "$script_path" ]]; then
    echo "WARN: No merge.slurm in: $target_dir"
    SKIPPED_FOLDERS+=("$(basename "$target_dir") (no merge.slurm)")
    return 0
  fi

  echo "Submitting from: $target_dir"
  {  # run in current shell (preserves arrays)
    cd "$target_dir"
    out="$(sbatch ./merge.slurm 2>&1 || true)"
    echo "$out"
    if grep -q "Submitted batch job" <<<"$out"; then
      job_id="$(awk '/Submitted batch job/{print $4}' <<<"$out")"
      SUBMITTED_JOBS+=("$(basename "$target_dir") -> $job_id")
    else
      SKIPPED_FOLDERS+=("$(basename "$target_dir") (submit failed)")
    fi
  }
}

# --- Main loop ---
for f in "${FOLDERS[@]}"; do
  submit_merge_job "$PARENT_DIR/$f"
  echo
done

# --- Summary report ---
echo "==== Summary ===="
if ((${#SUBMITTED_JOBS[@]})); then
  echo "Submitted jobs:"
  printf '  %s\n' "${SUBMITTED_JOBS[@]}"
else
  echo "No jobs submitted."
fi

if ((${#SKIPPED_FOLDERS[@]})); then
  echo
  echo "Skipped:"
  printf '  %s\n' "${SKIPPED_FOLDERS[@]}"
fi
echo "================="
