#!/bin/bash

# Loop from 1 to 60
for i in {1..60}; do
  folder_name=$(printf "%02d" $i)
  if [ -d "$folder_name" ]; then
    rm -rf $folder_name
  fi
 mkdir $folder_name

#Copy Rcode and Slurm files to folders
cp trend_shift.R $folder_name;
cp run.slurm $folder_name;
cd $folder_name;

#Run simulation under folder
sbatch run.slurm;
cd ..
done;
