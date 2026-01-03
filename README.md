# 🧠 Single-Cell RNA-seq Analysis of Human Glioblastoma 

This repository provides a **complete end-to-end single-cell RNA sequencing (scRNA-seq) analysis pipeline** using **Seurat (R)** for **human glioblastoma brain tissue** generated with **10X Genomics CytAssist FFPE technology**.

The pipeline processes raw 10X `.h5` count matrices to:
- Perform quality control and filtering
- Identify transcriptionally distinct cell populations
- Detect cluster-specific marker genes
- Annotate biological cell types
- Characterize tumor microenvironment heterogeneity

---

##  Dataset Overview

| Attribute | Description |
|---------|------------|
| Platform | 10X Genomics CytAssist (FFPE) |
| Tissue | Human glioblastoma brain tumor |
| Data type | Single-cell RNA-seq |
| Input format | `.h5` raw feature-barcode matrix |
| Species | *Homo sapiens* |

---

##  Objectives

- Identify **cellular heterogeneity** in glioblastoma
- Distinguish **tumor, immune, neuronal, and stromal populations**
- Extract **cluster-specific marker genes**
- Generate an **annotated single-cell atlas** of glioblastoma tissue

---

##  Software Requirements

### R Environment
- **R ≥ 4.2**
- **Seurat (v4 / v5)**

### Required R Packages
```r
install.packages(c(
  "Seurat",
  "tidyverse",
  "dplyr",
  "ggplot2",
  "patchwork",
  "hdf5r"
))
