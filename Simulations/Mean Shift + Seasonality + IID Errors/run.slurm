#!/bin/bash
#SBATCH --time=11:15:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=1024       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=meanshift-problem

module load R/4.2
R CMD BATCH --vanilla mean_shift_seasonality.R