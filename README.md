# TCGA-BRCA RNA-seq Differential Expression Pipeline (Tumor vs Normal)
A complete, beginner-friendly R pipeline for identifying differentially expressed genes (DEGs) between tumor and normal breast tissue using **TCGA-BRCA** RNA-seq data, accessed via `recount3` and analyzed with `DESeq2` — extended with a **survival analysis module** linking the top DEGs to patient clinical outcomes.

## Project Overview
- **Data type:** Bulk RNA-seq gene expression counts (Illumina platform)
- **Source:** TCGA-BRCA (The Cancer Genome Atlas – Breast Cancer), accessed via the `recount3` R package
- **Method:** DESeq2 differential expression analysis, Tumor vs Normal, extended with Kaplan-Meier and Cox proportional hazards survival analysis on top DEGs
- **DEG criteria:** `padj < 0.05` **and** `|log2FoldChange| > 1`
- **Outputs:** DEG results table (CSV) + 6 QC/analysis plots + survival curves and a Cox regression summary table for candidate prognostic genes

## Repository Structure
```
├── scripts/
│   └── TCGA_BRCA_DESeq2_pipeline.R              # Full annotated DEG pipeline script
├── survival_analysis/
│   ├── survival_analysis.R                      # Kaplan-Meier + Cox regression on top DEGs
│   └── plots/                                   # KM survival curves (generated on run)
├── docs/
│   ├── BRCA_RNAseq_Overview_Presentation.pptx   # Project overview slides
│   ├── DESeq2_Line_by_Line_Teaching_Notes.pptx  # Step-by-step teaching notes (slides)
│   └── DESeq2_Line_by_Line_Teaching_Notes.md    # Same notes, as Markdown (renders on GitHub)
├── results/                                     # Pipeline outputs (CSVs + plots)
└── README.md
```

## Results
Run on the full TCGA-BRCA dataset: **1,249 samples** (1,135 Tumor + 114 Normal) × **63,856 genes** (~39k after low-expression filtering) → **10,683 DEGs** (7,190 UP, 3,493 DOWN in Tumor vs Normal).

**Quality control — before vs after VST normalization:**
| Before | After |
|---|---|
| ![Raw counts boxplot](results/boxplot_raw_TCGA.png) | ![VST-normalized boxplot](results/boxplot_vst_TCGA.png) |

**Volcano plot** — all filtered genes, UP (red) / DOWN (blue) / not significant (grey):
![Volcano plot](results/volcano_TCGA.png)

**PCA** — Tumor vs Normal separation:
![PCA plot](results/PCA_TCGA.png)

**MA plot** — fold change vs mean expression:
![MA plot](results/MA_TCGA.png)

**Heatmaps** — top DEGs, expression centered per gene:
| Top 25 DEGs (samples ordered by condition) | Top 30 DEGs (hierarchically clustered) |
|---|---|
| ![Heatmap top 25](results/heatmap_top25_TCGA.png) | ![Heatmap top 30 clustered](results/heatmap_TCGA.png) |

Full result tables: [`DEG_results_TCGA_SIGNIFICANT_ONLY.csv`](results/DEG_results_TCGA_SIGNIFICANT_ONLY.csv) (10,683 DEGs only) and [`DEG_results_TCGA.csv`](results/DEG_results_TCGA.csv) (all ~39k filtered genes).

## Survival Analysis Extension
Differential expression tells you *which* genes differ between tumor and normal tissue — it doesn't tell you whether that difference actually matters for patients. This extension takes the top DEGs from the analysis above and asks: **is this gene's expression level associated with how long patients survive?**

**What it does:**
1. Pulls TCGA-BRCA clinical/survival metadata via `recount3` (same source, same samples as the DEG pipeline — no separate download needed).
2. Takes the top 10 most significant DEGs (by adjusted p-value) from the results table above.
3. For each gene, splits patients into **High vs. Low expression** groups (median split) and runs:
   - **Kaplan-Meier estimation** — visualizes survival curves for the two groups.
   - **Log-rank test** — tests whether the two curves differ significantly.
   - **Cox proportional hazards regression** — uses the *continuous* expression value (not just High/Low) to estimate a hazard ratio per unit increase in expression, which is the more statistically rigorous result.

**Why this matters:** a gene can be strongly differentially expressed between tumor and normal tissue and still have no relationship to patient outcome — and conversely, a modestly differentially expressed gene can be a strong prognostic signal. Running this step is what starts to turn a differential expression analysis into a **biomarker discovery** analysis.

**Outputs:**
- `survival_analysis/plots/KM_<gene>.png` — Kaplan-Meier survival curves, generated only for genes with log-rank p < 0.05
- `survival_analysis/survival_analysis_results.csv` — combined results table: log-rank p-value, Cox hazard ratio, and Cox p-value for every gene tested

**Important caveat:** this is a discovery-level analysis on a single cohort. Any candidate prognostic gene identified here would need to be validated on an independent cohort before being considered clinically meaningful — this pipeline is a starting point for hypothesis generation, not a clinical claim.

See [`survival_analysis/README.md`](survival_analysis/README.md) for full methodology, requirements, and run instructions specific to this module.

## How to Run

**1. DEG pipeline:**
1. Open `scripts/TCGA_BRCA_DESeq2_pipeline.R` in RStudio.
2. Set your working directory to wherever you want results saved (`Session > Set Working Directory > Choose Directory`).
3. Run the whole script (`Source` button, or `Ctrl+A` then `Ctrl+Enter`).
4. The script takes ~15–30 minutes (data download + DESeq2 run). Do not close RStudio while it runs.
5. A `results/` folder will be created containing the DEG table (CSV) and plots (PNG).

**2. Survival analysis (run after step 1, since it reads the DEG results):**
1. Open `survival_analysis/survival_analysis.R` in RStudio, with your working directory set to the `survival_analysis/` folder.
2. Run the script (`Source`, or `Ctrl+A` then `Ctrl+Enter`).
3. Outputs are saved to `survival_analysis/plots/` and `survival_analysis/survival_analysis_results.csv`.

### Requirements
- R (≥ 4.2 recommended)
- **DEG pipeline:** `recount3`, `DESeq2`, `ggplot2`, `pheatmap`, `ggrepel` (install instructions are included at the top of the script)
- **Survival analysis:** `recount3`, `survival`, `survminer`, `dplyr`, `readr`, `readxl`, `janitor`

## Documentation
- **Overview presentation** ([`docs/BRCA_RNAseq_Overview_Presentation.pptx`](docs/BRCA_RNAseq_Overview_Presentation.pptx)) — background on NGS, TCGA-BRCA, and the analysis workflow, aimed at someone new to RNA-seq.
- **Line-by-line teaching notes** ([`docs/DESeq2_Line_by_Line_Teaching_Notes.md`](docs/DESeq2_Line_by_Line_Teaching_Notes.md)) — walks through all 17 steps of the pipeline script in detail, explaining the biological and statistical reasoning behind each block of code. Also available as slides.
- **Survival analysis methodology** ([`survival_analysis/README.md`](survival_analysis/README.md)) — detailed write-up of the Kaplan-Meier/Cox regression methodology and how it connects back to the DEG pipeline.

