#!/usr/bin/env bash
# jobstat.sh — summarize job completion by numeric folders
# A folder counts as COMPLETE if it contains at least one .csv file (default extension)

EXT="${1:-csv}"
completed=0
remaining=0
incomplete_list=()

# Find all numeric directories in the current path, sorted numerically
folders=$(find . -maxdepth 1 -type d -printf "%f\n" | grep -E '^[0-9]+$' | sort -n)

if [[ -z "$folders" ]]; then
  echo "No numeric folders found."
  exit 0
fi

for folder in $folders; do
  count=$(find "$folder" -maxdepth 1 -type f -name "*.${EXT}" | wc -l)
  if [[ $count -gt 0 ]]; then
    ((completed++))
  else
    ((remaining++))
    incomplete_list+=("$folder")
  fi
done

echo "Completed: $completed"
echo "Remaining: $remaining"

if (( remaining > 0 )); then
  echo -n "Incomplete folders: "
  (IFS=", "; echo "${incomplete_list[*]}")
fi
