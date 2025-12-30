# Installing Segmentor3IsBack and AR1seg on HCC

This README explains how to set up and install the Segmentor3IsBack and AR1seg R packages on the Holland Computing Center (HCC). Both packages are archived on CRAN and must be installed manually from their tarballs.

---

## Step 1 — Download the Tarballs

Download the archived `.tar.gz` files from CRAN:

Segmentor3IsBack: https://cran.r-project.org/src/contrib/Archive/Segmentor3IsBack/Segmentor3IsBack_2.0.tar.gz  
AR1seg: https://cran.r-project.org/src/contrib/Archive/AR1seg/AR1seg_1.0.tar.gz

Save both files locally.

---

## Step 2 — Upload to HCC

Upload the tarballs to your personal R library directory. For me, it is 

/work/xueheng/mgrantham/r_lib

---

## Step 3 — Load R on HCC

SSH into HCC and load the R module:

module load R

Confirm that R loaded correctly:

R --version

---

## Step 4 — Install the Packages

Start R from the terminal:

R

Then inside the R console, run:

lib_path <- "/work/xueheng/mgrantham/r_lib"  
dir.create(lib_path, showWarnings = FALSE, recursive = TRUE)  
.libPaths(lib_path)

install.packages(file.path(lib_path, "Segmentor3IsBack_2.0.tar.gz"), repos = NULL, type = "source", lib = lib_path)  
install.packages(file.path(lib_path, "AR1seg_1.0.tar.gz"), repos = NULL, type = "source", lib = lib_path)

---

## Step 5 — Verify Installation

library(Segmentor3IsBack, lib.loc = "/work/xueheng/mgrantham/r_lib")  
library(AR1seg, lib.loc = "/work/xueheng/mgrantham/r_lib")  

packageVersion("Segmentor3IsBack")  
packageVersion("AR1seg")

If version numbers (e.g. 2.0 and 1.0) appear, installation succeeded.

---

## Step 6 — Make Library Path Persistent

To avoid resetting paths each time, edit (or create) your ~/.Rprofile file and add:

.libPaths(c("/work/xueheng/mgrantham/r_lib", .libPaths()))

This ensures R automatically checks your personal library whenever it starts.

---

## Step 7 — Add to the Top of Every Script

At the very top of all your R scripts, include this line:

.libPaths(c("/work/xueheng/mgrantham/r_lib", .libPaths()))

This ensures scripts find your installed libraries regardless of which compute node is used.

---

## Directory Layout Example

/work/xueheng/mgrantham/  
├── r_lib/  
│   ├── Segmentor3IsBack_2.0.tar.gz  
│   ├── AR1seg_1.0.tar.gz  
│   ├── Segmentor3IsBack/  
│   ├── AR1seg/  
│   └── README.md  
└── AR(1) Errors/  
    ├── mean_shift_ar1_0.0/  
    ├── AR1 Group Plots
    └── ...

---

## Troubleshooting

• Permission denied: ensure installation targets /work/xueheng/mgrantham/r_lib  
• Missing dependency: install it using install.packages("pkgname", lib="/work/xueheng/mgrantham/r_lib")  
• Modules not found: make sure `module load R` has been executed before running R

---

Maintained by: Michael Grantham  
HCC Path: /work/xueheng/mgrantham/r_lib  
Last updated: October 2025