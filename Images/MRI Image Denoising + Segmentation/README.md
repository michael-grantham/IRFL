# MRI Brain Scan Image Denoising + Segmentation (IRFL)

This repository contains code supporting the dissertation and associated manuscripts on the **Iteratively Reweighted Fused Lasso (IRFL)**.

The purpose of this repository is to enable reproducibility of the figures and experimental results presented in the dissertation.

---

## Dataset attribution

Some figures use the **Brain MRI Images for Brain Tumor Detection** dataset originally released on Kaggle by Navoneel Chakrabarty.

The dataset consists of labeled 2D brain MRI slices indicating the presence or absence of a tumor.  
Class labels follow the original Kaggle naming convention:

- Tumor present: files named `Y1.jpg`, `Y2.jpg`, …
- No tumor: files named `1 no.jpeg`, `2 no.jpeg`, …

The dataset is publicly available at:

https://www.kaggle.com/datasets/navoneel/brain-mri-images-for-brain-tumor-detection

This dataset is used solely for methodological illustration and benchmarking.  
It is not intended for clinical or diagnostic validation.

---

## Reproducibility

All figures involving brain MRI images that appear in the dissertation can be reproduced using the code in this repository **once the dataset is obtained from Kaggle** and placed in the expected directory structure.

Exact preprocessing steps, random seeds, and plotting routines are provided to ensure reproducibility of the reported images and results.

---

## Notes

- The dataset is **not redistributed** in this repository.
- Users must obtain the dataset directly from Kaggle and comply with its terms of use.
- Results should be interpreted as algorithmic demonstrations rather than medical conclusions.
