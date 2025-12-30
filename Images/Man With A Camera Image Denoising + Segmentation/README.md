# Man with a Camera — Image Denoising + Segmentation (IRFL)

This repository contains code for **image denoising and segmentation** using IRFL, demonstrated on the classic *Man with a Camera* test image.

The purpose of this repository is to support **reproducibility of figures and algorithmic behavior** for denoising, iterative refinement, and edge/segmentation visualization.

---

## Data description

This repository includes the following image files:

- `camera_original_128x128.png` — clean reference image  
- `camera_heavily_noised_128x128.png` — synthetically corrupted version used as input  

These images are standard benchmark-style test images commonly used in image processing and signal denoising demonstrations.

No external datasets are required for this script.

---

## Reproducibility

All denoised outputs, iteration grids, and edge overlay figures can be reproduced directly using the code in this repository with the provided PNG files.

The scripts specify all parameters, iteration counts, and plotting logic needed to regenerate the results exactly.

---


