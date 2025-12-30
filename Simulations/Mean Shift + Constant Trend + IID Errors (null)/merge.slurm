#!/bin/bash
#SBATCH --time=01:15:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=2048       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=merge-out
 
module load R/4.3
R CMD BATCH --vanilla merge.R
