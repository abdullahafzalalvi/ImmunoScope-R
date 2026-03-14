# ImmunoScope-R
ImmunoScope-R
```
╔══════════════════════════════════════════════════════════════════════════════════╗
║                                                                                  ║
║    ░██████╗░█████╗░ ██████╗░███╗░░██╗░█████╗░  ░██████╗███████╗░██████╗░        ║
║    ██╔════╝██╔══██╗ ██╔══██╗████╗░██║██╔══██╗  ██╔════╝██╔════╝██╔═══██╗        ║
║    ╚█████╗░██║░░╚═╝ ██████╔╝██╔██╗██║███████║  ╚█████╗░█████╗░░██║██╗██║        ║
║    ░╚═══██╗██║░░██╗ ██╔══██╗██║╚████║██╔══██║  ░╚═══██╗██╔══╝░░╚██████╔╝        ║
║    ██████╔╝╚█████╔╝ ██║░░██║██║░╚███║██║░░██║  ██████╔╝███████╗░╚═██╔═╝░        ║
║    ╚═════╝░░╚════╝  ╚═╝░░╚═╝╚═╝░░╚══╝╚═╝░░╚═╝  ╚═════╝░╚══════╝░░░╚═╝░░        ║
║                                                                                  ║
║          Single-Cell RNA-Seq Analysis of Human Peripheral Blood                  ║
║                   Mononuclear Cells (PBMCs) — Full Pipeline                      ║
║                                                                                  ║
╚══════════════════════════════════════════════════════════════════════════════════╝
```

<div align="center">

![R](https://img.shields.io/badge/Language-R-276DC3?style=for-the-badge&logo=r&logoColor=white)
![scRNA-seq](https://img.shields.io/badge/Analysis-scRNA--seq-FF6B6B?style=for-the-badge)
![PBMC](https://img.shields.io/badge/Dataset-PBMC%203k%20Mirror-4ECDC4?style=for-the-badge)
![Figures](https://img.shields.io/badge/Figures-12%20Publication--Ready-gold?style=for-the-badge)
![License](https://img.shields.io/badge/License-Academic%20Use-blueviolet?style=for-the-badge)

</div>

---

## 🧬 Project Overview

> *"Understanding immune diversity at the resolution of a single cell."*

This repository contains a **complete single-cell RNA sequencing (scRNA-seq) analysis pipeline** built entirely in **base R and ggplot2**, designed to mirror the landmark **10X Genomics PBMC 3k dataset**. The pipeline simulates realistic transcriptomic data for **2,638 human peripheral blood mononuclear cells (PBMCs)** spanning **9 distinct immune cell populations**, and generates **12 publication-quality figures** covering every major step of a standard scRNA-seq workflow — from quality control to dashboard summary.

This project was developed as a **computational immunology training resource** and a **portfolio-grade demonstration** of scRNA-seq analysis without reliance on Seurat or Bioconductor, making it exceptionally accessible and reproducible across computing environments.

---

## 👨‍🔬 Researchers

```
┌─────────────────────────────────────────────────────────┐
│                                                         │
│   🔬  Muzzamil                                          │
│   🔬  Abdullah Alvi                                     │
│       BSc (Hons.) Biotechnology                         │
│       University of Layyah, Pakistan                    │
│                                                         │
│   🌐  Google Scholar:                                   │
│       https://scholar.google.com/citations?             │
│       user=2kmK6UwAAAAJ&hl=en                           │
│                                                         │
│   🔗  ORCID: 0009-0006-7961-8226                        │
│                                                         │
└─────────────────────────────────────────────────────────┘
```

---

## 🗂️ Repository Structure

```
📦 scRNA-Seq-PBMC/
 ┃
 ┣ 📄 Single-Cell_RNA-Seq_Analysis_of_Human_PBMCs.R   ← Main analysis script
 ┣ 📄 README.md                                        ← You are here
 ┃
 ┣ 📂 figures/
 ┃  ┣ 🖼️ fig01_QC_violin.png          ← QC metrics violin plots
 ┃  ┣ 🖼️ fig02_UMAP_cluster.png       ← UMAP by unsupervised cluster
 ┃  ┣ 🖼️ fig03_UMAP_annotated.png     ← UMAP by annotated cell type
 ┃  ┣ 🖼️ fig04_feature_plots.png      ← Marker gene feature plots
 ┃  ┣ 🖼️ fig05_violin_markers.png     ← Violin plots per marker gene
 ┃  ┣ 🖼️ fig06_dotplot.png            ← Dot plot: expression summary
 ┃  ┣ 🖼️ fig07_heatmap.png            ← Z-scored marker gene heatmap
 ┃  ┣ 🖼️ fig08_volcano.png            ← DE volcano plot (B cells)
 ┃  ┣ 🖼️ fig09_elbow.png              ← PCA elbow plot
 ┃  ┣ 🖼️ fig10_barplot.png            ← Cell proportion bar chart
 ┃  ┣ 🖼️ fig11_stacked.png            ← Stacked bar composition
 ┃  ┗ 🖼️ fig12_dashboard.png          ← 4-panel summary dashboard
 ┃
 ┗ 📂 data/
    ┗ (simulated — generated at runtime via set.seed(42))
```

---

## 🧫 Dataset Specifications

```
╔══════════════════════════════════════════════════╗
║          SIMULATED PBMC DATASET SUMMARY          ║
╠══════════════════════════════════════════════════╣
║  Total Cells          │  2,638                  ║
║  Marker Genes         │  45 canonical markers   ║
║  Cell Type Clusters   │  9 populations          ║
║  UMAP Dimensions      │  2D (UMAP 1 & 2)        ║
║  Random Seed          │  42 (reproducible)      ║
║  Resolution           │  0.5 (Louvain)          ║
║  PCA Dims Used        │  1:10                   ║
╠══════════════════════════════════════════════════╣
║          CELL TYPE COMPOSITION                   ║
╠══════════════════════════════════════════════════╣
║  Naive CD4 T          │  32%  (~844 cells)      ║
║  CD14+ Monocytes      │  18%  (~475 cells)      ║
║  Memory CD4 T         │  15%  (~396 cells)      ║
║  B Cells              │  10%  (~264 cells)      ║
║  CD8 T Cells          │  10%  (~264 cells)      ║
║  FCGR3A+ Monocytes    │  07%  (~185 cells)      ║
║  NK Cells             │  04%  (~106 cells)      ║
║  Dendritic Cells      │  03%   (~79 cells)      ║
║  Platelets            │  01%   (~26 cells)      ║
╚══════════════════════════════════════════════════╝
```

---

## 🧪 Canonical Marker Genes

| Cell Type | Marker Genes |
|---|---|
| 🔴 Naive CD4 T | `IL7R`, `CCR7`, `CD3D`, `SELL`, `TCF7` |
| 🔵 CD14+ Mono | `CD14`, `LYZ`, `CST3`, `S100A8`, `S100A9` |
| 🟠 Memory CD4 T | `S100A4`, `ANXA1`, `LTB`, `CD3E`, `IL32` |
| 🟢 B Cell | `MS4A1`, `CD79A`, `CD79B`, `BANK1`, `IGHM` |
| 🟣 CD8 T | `CD8A`, `CD8B`, `GZMK`, `GZMA`, `PRF1` |
| 🩵 FCGR3A+ Mono | `FCGR3A`, `MS4A7`, `IFITM3`, `AIF1`, `LST1` |
| 🟡 NK Cell | `GNLY`, `NKG7`, `KLRD1`, `GZMB`, `FGFBP2` |
| 🟤 Dendritic | `FCER1A`, `HLA-DQA1`, `CLEC10A`, `CD1C`, `ITGAX` |
| ⚪ Platelet | `PPBP`, `PF4`, `GNG11`, `SDPR`, `SPARC` |

---

## 📊 Figures Generated

```
┌─────────────────────────────────────────────────────────────────────┐
│  Figure  │  Type               │  Description                       │
├─────────────────────────────────────────────────────────────────────┤
│  Fig 01  │  Violin + Boxplot   │  QC metrics (nFeature, nCount, %MT)│
│  Fig 02  │  UMAP               │  Unsupervised Louvain clustering    │
│  Fig 03  │  UMAP               │  Annotated cell type landscape      │
│  Fig 04  │  Feature Plots      │  8 canonical marker genes on UMAP   │
│  Fig 05  │  Violin Plots       │  Per-marker expression by cell type  │
│  Fig 06  │  Dot Plot           │  Expression breadth + level summary │
│  Fig 07  │  Heatmap            │  Z-scored marker gene matrix        │
│  Fig 08  │  Volcano Plot       │  DE: B cells vs all other PBMCs     │
│  Fig 09  │  Elbow Plot         │  PCA variance explained per PC      │
│  Fig 10  │  Bar Chart          │  Cell type proportion (%)           │
│  Fig 11  │  Stacked Bar        │  Compositional overview             │
│  Fig 12  │  Dashboard (4-panel)│  UMAP + Volcano + Bar + Elbow       │
└─────────────────────────────────────────────────────────────────────┘
```

All figures are exported at **300 DPI** and sized for publication (typically 14–16 × 6–12 inches).

---

## 📦 Dependencies

```r
# Core visualization
library(ggplot2)      # Grammar of graphics
library(ggrepel)      # Non-overlapping text labels
library(gridExtra)    # Multi-panel figure layout
library(viridis)      # Perceptually uniform color scales

# Data transformation
library(dplyr)        # Data wrangling
library(reshape2)     # Data melting for ggplot

# Specialized plotting
library(pheatmap)     # Heatmap with annotations
library(RColorBrewer) # Color palettes
library(scales)       # Axis scale helpers
```

**Installation (run once):**
```r
install.packages(c("ggplot2", "ggrepel", "gridExtra", "viridis",
                   "dplyr", "reshape2", "pheatmap",
                   "RColorBrewer", "scales"))
```

---

## ▶️ How to Run

```r
# 1. Clone or download the repository
# 2. Open R or RStudio
# 3. Set your working directory
setwd("path/to/scRNA-Seq-PBMC/")

# 4. Run the full script
source("Single-Cell_RNA-Seq_Analysis_of_Human_PBMCs.R")

# All 12 figures will be generated and saved automatically.
# Results are fully reproducible due to set.seed(42).
```

> ⚠️ **Note:** The script uses `install.packages()` calls inline for convenience. It is recommended to run those once separately before sourcing the full script to avoid repeated prompts.

---

## 🔬 Analytical Workflow

```
Raw Count Matrix (Simulated 10X)
          │
          ▼
  ┌───────────────┐
  │  QC Filtering  │  nFeature > 200 & < 2500 | percent.mt < 5%
  └───────┬───────┘
          │
          ▼
  ┌────────────────────┐
  │  Normalization +    │  Log-normalize | Scale data
  │  Feature Selection  │  Top variable features identified
  └──────────┬─────────┘
             │
             ▼
  ┌──────────────────┐
  │  PCA Reduction    │  50 PCs computed | dims 1:10 selected
  │  (Elbow Method)   │  via elbow plot inflection
  └────────┬─────────┘
           │
           ▼
  ┌────────────────────┐
  │  UMAP Embedding    │  2D manifold from PCA dims 1:10
  └────────┬───────────┘
           │
           ▼
  ┌──────────────────────┐
  │  Louvain Clustering  │  Resolution = 0.5 → 9 clusters
  └────────┬─────────────┘
           │
           ▼
  ┌──────────────────────┐
  │  Cell Type Annotation│  Canonical markers + manual curation
  └────────┬─────────────┘
           │
           ▼
  ┌────────────────────────────┐
  │  Differential Expression   │  B cells vs all | Volcano plot
  │  (Simulated DE Results)    │  |log2FC| > 1 | -log10p > 3
  └────────────────────────────┘
```

---

## 📌 Key Design Features

- **Custom `theme_scrna()`** — A bespoke ggplot2 theme with publication-grade typography, navy title colors, and clean panel styling applied consistently across all figures.
- **Researcher branding** — All figures carry an embedded caption crediting Muzzamil & Abdullah Alvi.
- **Reproducibility** — Global `set.seed(42)` ensures identical output on every run.
- **Zero Seurat dependency** — The entire pipeline runs on base R + ggplot2, making it lightweight and platform-agnostic.
- **300 DPI exports** — All figures saved at print resolution suitable for journal submission.

---

## 📚 References & Biological Context

This project mirrors the methodology of the canonical **10X Genomics PBMC 3k tutorial**, which has become a standard benchmark in single-cell genomics. Key biological references:

- Stuart & Satija (2019). *Integrative single-cell analysis.* Nature Reviews Genetics.
- Butler et al. (2018). *Integrating single-cell transcriptomic data across different conditions.* Nature Biotechnology.
- 10X Genomics PBMC 3k Dataset: [https://www.10xgenomics.com/datasets](https://www.10xgenomics.com/datasets)

---

## 🌐 Connect with the Author

<div align="center">

| Platform | Link |
|---|---|
| 🎓 **Google Scholar** | [Abdullah Alvi — Publications & Citations](https://scholar.google.com/citations?user=2kmK6UwAAAAJ&hl=en) |
| 🔬 **ORCID** | [0009-0006-7961-8226](https://orcid.org/0009-0006-7961-8226) |

</div>

---

## 📄 License & Citation

This project is intended for **academic and educational use**. If you use or adapt this pipeline in your work, please cite the researchers:

```
Muzzamil & Abdullah Alvi (2025).
Single-Cell RNA-Seq Analysis of Human PBMCs — Full R Pipeline.
GitHub Repository.
Google Scholar: https://scholar.google.com/citations?user=2kmK6UwAAAAJ&hl=en
```

---

<div align="center">

```
╔══════════════════════════════════════════════════════╗
║                                                      ║
║   Built with curiosity, R, and a deep respect        ║
║   for the complexity of the human immune system.     ║
║                                                      ║
║            — Muzzamil & Abdullah Alvi —              ║
║                                                      ║
╚══════════════════════════════════════════════════════╝
```

*"Every cell tells a story. scRNA-seq lets us read them all."*

</div>
