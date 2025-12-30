# Mean Shift AR(1) Simulation Workflow

This project runs a set of AR(1) mean-shift simulations across multiple parameter folders, merges the results, and then produces group-level graphics.

---

## 1. Prerequisite: Set Up the Private R Library (r_lib)

Before running any simulations or plotting scripts, you must configure your private R library. Without this step, the scripts will fail to load required packages on the cluster.

### Step 1 — Create the Private Library Directory

Choose a persistent location in your workspace or project directory, for example:

mkdir -p /work/<username>/r_lib

Replace <username> with your actual cluster username or project folder.

### Step 2 — Install Required Packages

Install all required R packages into this directory. The key custom packages must be installed from local tarballs:

install.packages("/work/<username>/Segmentor3IsBack_1.0.tar.gz",
                 repos = NULL, type = "source",
                 lib = "/work/<username>/r_lib")

install.packages("/work/<username>/AR1seg_1.0.tar.gz",
                 repos = NULL, type = "source",
                 lib = "/work/<username>/r_lib")

Then install the supporting CRAN packages:

install.packages(c("genlasso", "Matrix", "dplyr", "changepoint", "wbs"),
                 lib = "/work/<username>/r_lib",
                 repos = "https://cloud.r-project.org")

### Step 3 — Add the Library Path in All Scripts

At the top of every R script (simulation, merging, or plotting), ensure the following line appears:

.libPaths(c("/work/<username>/r_lib", .libPaths()))

This ensures R will load all packages from your private library before checking system locations.

Once that line is present, packages can be loaded normally:

library(Segmentor3IsBack)
library(AR1seg)
library(genlasso)

✅ Important: Do not skip this setup. The cluster’s default R environment may not have these packages available, and jobs will fail to start without a properly configured .libPaths().

---

## 2. Run the Simulations

Each parameter setting has its own folder, for example:

mean_shift_ar1_0.0/
mean_shift_ar1_0.3/
mean_shift_ar1_0.7/
mean_shift_ar1_neg0.3/
mean_shift_ar1_neg0.7/

Inside each folder there is a script named run.sh that performs the simulation.

To run them, navigate into each folder and execute:

sh run.sh

This will generate simulation output files (e.g., final_out*.csv) inside that same folder.

Repeat this step for all parameter folders.

---

## 3. Merge All Simulation Results

After every simulation has finished, navigate to the group folder:

cd "AR1 Group Plots/"

Inside this directory, submit the merge job:

sbatch merge_all.slurm

This script visits all sibling simulation folders and runs their respective merge.slurm jobs to combine results.

Once complete, you will see merged output files appear in each simulation directory and/or in the group folder, depending on configuration.

---

## 4. Generate Group Plots

After the merges finish, stay in the AR1 Group Plots/ directory and run:

sbatch group_plot.slurm

This step calls your R plotting pipeline (e.g., AR1_Group_Plots.R) to create the aggregated graphics summarizing all simulation results.

The resulting plots will be saved in the AR1 Group Plots/ directory (typically as .png or .pdf files).

---

## 5. Summary of Workflow

Stage | Location | Command | Description
------|-----------|----------|-------------
Library setup | /work/<username>/ | — | Create and populate r_lib before running anything
Simulation | Each mean_shift_ar1_* folder | sh run.sh | Runs a single simulation scenario
Merge | AR1 Group Plots/ | sbatch merge_all.slurm | Merges all simulation outputs
Plot | AR1 Group Plots/ | sbatch group_plot.slurm | Produces combined visualizations

---

## 6. Notes

- All jobs are managed through Slurm.
- Run all scripts within their respective directories.
- Check job progress with squeue -u $USER and monitor logs with less or tail -f.
- Ensure .slurm scripts have execute permissions (chmod +x filename.slurm).
- To debug interactively, run .slurm scripts directly with bash instead of sbatch.
- Always verify that .libPaths() is set correctly in your R environment before submission.

---

### Example Complete Run

# Run each simulation
cd mean_shift_ar1_0.0 && sh run.sh && cd ..
cd mean_shift_ar1_0.3 && sh run.sh && cd ..
cd mean_shift_ar1_0.7 && sh run.sh && cd ..
cd mean_shift_ar1_neg0.3 && sh run.sh && cd ..
cd mean_shift_ar1_neg0.7 && sh run.sh && cd ..

# Merge and plot
cd "AR1 Group Plots"
sbatch merge_all.slurm
sbatch group_plot.slurm

---

End of README
