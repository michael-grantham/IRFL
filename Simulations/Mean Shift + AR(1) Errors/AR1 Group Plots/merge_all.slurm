#!/bin/bash
#SBATCH --job-name=run_merge_all
#SBATCH --output=run_merge_all.%j.out
#SBATCH --error=run_merge_all.%j.err
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

set -euo pipefail

# Submit from the "AR1 Group Plots" directory (or any child of the parent).
ROOT_DIR="${SLURM_SUBMIT_DIR}"
PARENT_DIR="$(dirname "$ROOT_DIR")"

echo "Submit dir           : $ROOT_DIR"
echo "Parent (one level up): $PARENT_DIR"
echo

# Sibling folders to process
FOLDERS=(
  "mean_shift_ar1_0.0"
  "mean_shift_ar1_0.3"
  "mean_shift_ar1_0.7"
  "mean_shift_ar1_neg0.3"
  "mean_shift_ar1_neg0.7"
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

  # Submit merge.slurm with the working directory set to the folder
  # (avoids cd and ensures relative paths inside merge.slurm work).
  echo "Submitting: $script_path"
  sbatch_out="$(sbatch --chdir="$target_dir" "$script_path" 2>&1 || true)"
  echo "$sbatch_out"

  # Extract the job ID if submission succeeded
  if grep -q "Submitted batch job" <<<"$sbatch_out"; then
    job_id="$(awk '/Submitted batch job/{print $4}' <<<"$sbatch_out")"
    SUBMITTED_JOBS+=("$(basename "$target_dir") -> $job_id")
  else
    echo "WARN: Failed to submit for $target_dir"
    SKIPPED_FOLDERS+=("$(basename "$target_dir") (submit failed)")
  fi
}

# Iterate and submit merge jobs
for f in "${FOLDERS[@]}"; do
  submit_merge_job "$PARENT_DIR/$f"
  echo
done

# Summary
echo "==== Summary ===="
if ((${#SUBMITTED_JOBS[@]})); then
  echo "Submitted jobs:"
  for line in "${SUBMITTED_JOBS[@]}"; do
    echo "  $line"
  done
else
  echo "No jobs submitted."
fi

if ((${#SKIPPED_FOLDERS[@]})); then
  echo
  echo "Skipped:"
  printf '  %s\n' "${SKIPPED_FOLDERS[@]}"
fi
echo "================="
