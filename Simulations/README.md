# Simulation - Chapter 4

This README describes the top-level workflow for running, monitoring, cleaning, merging, and plotting results found in Chapter 4 of the dissertation across all simulation projects in this repository.

The simulation infrastructure in this repository is designed primarily for execution on a high-performance computing (HPC) cluster managed by Slurm, where large parameter grids are evaluated in parallel using distributed batch jobs.

## 1. Execution Environment Assumptions

### Primary (Intended) Environment

- Slurm workload manager (e.g., sbatch, job arrays, batch scripts)
- Access to multiple compute nodes and sufficient aggregate memory
- Non-interactive batch execution of large simulation grids

The provided shell scripts and `.slurm` files assume:
- A shared filesystem visible to all compute nodes
- Independent simulation jobs writing intermediate results to disk
- Post-hoc aggregation and plotting performed via separate batch jobs

### Manual / Single-Node Execution (Optional)

All simulations can, in principle, be reproduced without a distributed computing environment, for example on a single workstation or laptop. However:

- Wall-clock runtime grows rapidly without parallelization, often super-linearly in grid size
- Peak memory usage may exceed typical single-node limits
- Long-running jobs are more susceptible to interruption or resource exhaustion
- Intermediate result files may accumulate faster than they can be merged or post-processed

In practice, reproducing full simulation grids outside an HPC environment may require:
- Reducing or subsampling parameter grids
- Converting Slurm batch scripts to sequential execution
- Careful monitoring of memory usage and disk I/O

## 2. Running Simulations (Slurm / HPC Workflow)

Most simulation directories (e.g., Mean Shift + IID Errors, Mean Shift + Trend, Seasonal Mean Shift, etc.) follow a standard structure.

Inside each directory, you will typically find:

- start.sh — submits all simulation jobs
- jobstat.sh — reports progress (completed vs. remaining jobs)
- delete.sh — removes all simulation outputs, logs, and intermediate files
- merge.slurm — aggregates results once simulations finish
- group_plot.slurm (optional) — generates summary plots

### Typical Commands

Submit all simulation jobs:

    sh start.sh

Monitor progress:

    sh jobstat.sh

Remove all previous outputs and logs:

    sh delete.sh

## 3. Merging and Plotting Results

After all simulation jobs have completed, results can be merged and visualized.

Merge per-job outputs into a single aggregated dataset:

    sbatch merge.slurm

Depending on the project, this step may also generate plots automatically.

If plotting is handled separately, run:

    sbatch group_plot.slurm

This produces summary figures (e.g., PNG or PDF files) based on the merged results.

## 4. Projects with Custom or Nested Structures

The following simulation groups use multi-level directory layouts, nested parameter grids, or additional aggregation logic:

- Mean Shift + AR(1) Errors
- Mean Shift + AR(1) Errors (null)
- Mean Shift + Constant Trend + Random Changepoints

For these projects, consult the local README.md files inside each directory. Those documents describe hierarchical job submission strategies, multi-stage merging across parameter subfolders, and project-specific plotting workflows.

## 5. Operational Notes

- All simulations are designed to be run under Slurm.
- Execute scripts from within their respective project directories.
- Progress monitoring is handled exclusively via jobstat.sh.
- Always run sh delete.sh before re-running a simulation batch to avoid mixing results.
- For debugging or exploratory testing, .slurm files may be executed directly with bash, but this bypasses scheduler-based resource allocation and isolation.

## 6. Example High-Level Workflow

    cd "Mean Shift + IID Errors"
    sh start.sh
    sh jobstat.sh
    sh delete.sh
    sbatch merge.slurm
    sbatch group_plot.slurm

For Mean Shift + AR(1) Errors, Mean Shift + AR(1) Errors (null), and Mean Shift + Constant Trend + Random Changepoints, refer to the project-specific READMEs for precise instructions.

End of Top-Level README
