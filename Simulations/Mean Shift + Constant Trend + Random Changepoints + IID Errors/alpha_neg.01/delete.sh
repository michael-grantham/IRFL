#!/bin/bash

###########################################
# delete.sh — Clean up simulation outputs
# Deletes:
#   1) Directories named only with digits (e.g., 01, 12, 305)
#   2) Files ending in .csv, .out, .Rout, or .pdf, or .Rout
###########################################

# --- Delete numbered directories ---
for folder_name in [0-9]*; do
  if [ -d "$folder_name" ] && [[ "$folder_name" =~ ^[0-9]+$ ]]; then
    rm -rf "$folder_name"
    echo "Deleted directory: '$folder_name'"
  fi
done

# --- Delete matching files ---
for ext in csv out pdf Rout; do
  for file in *.$ext; do
    # If no files match, the glob stays literal ('*.csv'), so skip it
    if [ -e "$file" ]; then
      rm -f "$file"
      echo "Deleted file: '$file'"
    fi
  done
done
