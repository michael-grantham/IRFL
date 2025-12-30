#!/bin/bash
#SBATCH --time=13:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=2048       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=meanshift-problem
 
module load R/4.3
R CMD BATCH --vanilla trend_shift.R
