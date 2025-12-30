#!/bin/bash
#SBATCH --job-name=group_plots_null
#SBATCH --output=group_plots_null.%j.out
#SBATCH --error=group_plots_null.%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -euo pipefail

# --- Submit this from the "AR1 Group Plots (null)" folder ---
AR1_DIR="${SLURM_SUBMIT_DIR}"
PARENT_DIR="$(dirname "$AR1_DIR")"

echo "AR1 Group Plots (null) dir : $AR1_DIR"
echo "Parent (one level up)      : $PARENT_DIR"
echo

# Folders to search (siblings of AR1 Group Plots (null))
FOLDERS=(
  "mean_shift_ar1_null_0.0"
  "mean_shift_ar1_null_0.3"
  "mean_shift_ar1_null_0.7"
  "mean_shift_ar1_null_neg0.3"
  "mean_shift_ar1_null_neg0.7"
)

# Function: copy newest final_out*.csv from a folder into AR1_DIR
copy_newest_final_out() {
  local src_dir="$1"
  if [[ ! -d "$src_dir" ]]; then
    echo "WARN: Missing folder: $src_dir"
    return 0
  fi

  # Find newest matching file (by modification time)
  local newest
  newest="$(find "$src_dir" -maxdepth 1 -type f -name 'final_out*.csv' -printf '%T@ %p\n' \
            | sort -nr | head -n1 | cut -d' ' -f2- || true)"

  if [[ -z "${newest:-}" ]]; then
    echo "WARN: No final_out*.csv in $src_dir"
    return 0
  fi

  echo "Copying: $newest -> $AR1_DIR/"
  cp -f "$newest" "$AR1_DIR"/
}

# --- Copy phase ---
for f in "${FOLDERS[@]}"; do
  src_path="$PARENT_DIR/$f"
  copy_newest_final_out "$src_path"
done

echo
echo "Copied files now present in AR1 Group Plots (null):"
ls -1 "$AR1_DIR"/final_out*.csv || echo "No final_out*.csv files found after copy."
echo

# --- Run the plotting script ---
module purge >/dev/null 2>&1 || true
module load R/4.2 >/dev/null 2>&1 || { echo "ERROR: Could not load R/4.2 module"; exit 1; }

cd "$AR1_DIR"
echo "Running AR1_Null_Group_Plots.R in: $(pwd)"
Rscript --vanilla AR1_Null_Group_Plots.R
echo "Done."
