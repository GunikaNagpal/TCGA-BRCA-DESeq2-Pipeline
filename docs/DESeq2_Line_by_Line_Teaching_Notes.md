# TCGA-BRCA DESeq2 Pipeline — Complete Line-by-Line Code Notes

Tumor vs Normal | `recount3` + `DESeq2`

**Results at a glance:** 1,249 samples · 63,856 genes → ~39k after filtering · 10,683 DEGs · 6 output plots

---

## Pipeline Overview — 17 Steps at a Glance

| Step | What happens |
|---|---|
| 0–1 | Install & load packages (`BiocManager`, `recount3`, `DESeq2`, `ggplot2`, `pheatmap`, `ggrepel`) |
| 2–3 | Find TCGA-BRCA in recount3 → download RSE object (1,256 samples × 63,856 genes) |
| 4 | `transform_counts()` → convert coverage to TRUE raw integer counts |
| 4b | Boxplot BEFORE normalization (QC baseline) |
| 5–6 | Decode TCGA barcodes (chars 14–15) → separate Tumor (01) vs Normal (11) → build `coldata` |
| 7–8 | Create `DESeqDataSet` object → filter low-expression genes (~63k → ~39k genes) |
| 9 | `DESeq()` — normalization + dispersion estimation + Negative Binomial GLM + Wald test |
| 9b | Boxplot AFTER normalization (VST) — confirm boxes now align |
| 10–11 | Extract results → apply dual filter: `padj < 0.05` AND `\|log2FC\| > 1` |
| 12 | Save DEG CSV (gene, log2FC, padj, UP/DOWN) |
| 13–16 | Generate 4 plots: Volcano, PCA, Heatmap (top 25 DEGs), MA plot |

---

## Section 0 — Key Concepts You Must Know

- **TCGA-BRCA** — The Cancer Genome Atlas – Breast Cancer dataset. 1,249 real patient samples with RNA-seq data from both tumor and normal breast tissue.
- **RNA-seq** — measures which genes are active and how much. Output is a matrix: rows = ~20,000 genes, columns = samples, values = read counts.
- **Raw count** — a whole integer (0, 47, 832...) showing how many RNA reads mapped to a gene. DESeq2 requires integers, not FPKM/TPM.
- **log2FC** — `log2(Tumor / Normal)`. `log2FC = 1` means 2× higher in Tumor; `log2FC = -1` means 2× lower in Tumor. Threshold used here: `|log2FC| > 1`.
- **padj** — Benjamini-Hochberg adjusted p-value. Controls the false discovery rate when testing ~20,000 genes at once. `padj < 0.05` = significant.
- **DEG** — Differentially Expressed Gene. Requires **both** `padj < 0.05` **and** `|log2FC| > 1`. UP = higher in Tumor, DOWN = lower in Tumor.

---

## Step 0 & 1 — Install & Load Packages

**Step 0: Install (run ONCE ever)**
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("recount3", update = FALSE, ask = FALSE)
BiocManager::install("DESeq2",   update = FALSE, ask = FALSE)
install.packages(c("ggplot2", "pheatmap", "ggrepel"))
```
- `if (!requireNamespace(...))` checks IF BiocManager is not already installed.
- `BiocManager::install()` — the `::` means "use `install()` from the `BiocManager` package"; `update=F, ask=F` keeps it silent.
- `ggplot2`/`pheatmap`/`ggrepel` are plotting packages installed from CRAN.

**Step 1: Load (run every session)**
```r
library(recount3)   # enables available_projects(), create_rse()
library(DESeq2)     # enables DESeq(), results(), vst(), plotMA()
library(ggplot2)    # enables ggplot(), geom_point(), ggsave()
library(pheatmap)   # enables pheatmap()
library(ggrepel)    # enables geom_text_repel() for non-overlapping labels
```

**Key difference:** `install.packages()` downloads the tool (once). `library()` picks up the tool to use it today (every session).

**Why BiocManager?** DESeq2 and recount3 live on Bioconductor, a separate repository from CRAN. BiocManager is the bridge that installs Bioconductor packages with correct version matching.

---

## Step 2 & 3 — Find & Download TCGA-BRCA Data

**Step 2: Find the project**
```r
human_projects <- available_projects()
brca_project <- subset(human_projects,
                        project == "BRCA" &
                        project_home == "data_sources/tcga")
print(brca_project)
```
- `available_projects()` downloads the catalogue of all recount3 datasets.
- `subset()` filters rows: keep rows where `project == "BRCA"` **and** `project_home == "data_sources/tcga"` (so it's TCGA's BRCA data, not GTEx or another source).
- `print(brca_project)` shows the one matching row, confirming we found the right dataset.

**Step 3: Download data (takes 5–15 minutes)**
```r
rse <- create_rse(brca_project)
nrow(rse)   # 63,856 genes
ncol(rse)   # 1,256 samples
```
`create_rse()` downloads all expression data into an RSE (RangedSummarizedExperiment) object — a structured container holding (1) the count matrix, (2) gene annotations (chromosome, name), and (3) sample metadata (patient ID, tumor stage), all automatically aligned.

---

## Step 4 — Convert to True Raw Integer Counts

**Why this step exists:** recount3 internally stores data as *coverage* (base-level read depth sums), not standard gene-level counts. We must convert to integers before DESeq2 can use them.

```r
assay(rse, "counts") <- transform_counts(rse)
counts <- assay(rse, "counts")
print(counts[1:3, 1:3])   # sanity check: values must be whole integers
```

| Raw counts (integers) | FPKM / TPM (decimals) |
|---|---|
| 0, 47, 832, 10304... | 0.93, 14.72, 847.3... |
| DESeq2's Negative Binomial model works correctly | Violates DESeq2 assumptions — wrong results |

**Note:** GSE183947 (GEO) gives FPKM; recount3 gives true raw counts — use recount3.

---

## Step 4b — QC Boxplot (Before Normalization)

```r
log_counts_raw <- log2(counts + 1)
png("results/boxplot_raw_TCGA.png", width = 1200, height = 600)
boxplot(log_counts_raw[, 1:50], las = 2, cex.axis = 0.5, col = "lightblue")
dev.off()
```
- `log2(counts + 1)` — the `+1` avoids `log2(0) = -Infinity`.
- `[, 1:50]` plots the first 50 samples only; `las=2` rotates axis labels; `cex.axis=0.5` halves the font size.
- `png()` + `dev.off()` always go together — `png()` opens the file, `dev.off()` closes and saves it. Skipping `dev.off()` leaves an empty file.

**What to look for BEFORE normalization:** boxes at different heights are normal (each sample has different sequencing depth) — look for similar box *shapes*, not heights. A dramatically different box is a suspect sample.

**What to look for AFTER normalization (Step 9b):** boxes at roughly the same height confirms normalization worked. If medians still vary hugely, suspect a batch effect.

---

## Step 5 — Separate Tumor & Normal via TCGA Barcodes

Example barcode: `TCGA-A1-A0SK-01A-11R-A084-07`

| Segment | Meaning |
|---|---|
| `TCGA` | Always TCGA |
| `A1` | Study/institution |
| `A0SK` | Patient ID |
| `01` (chars 14–15) | **Sample type — this is what we extract** |
| `A / 11R / ...` | Other info (ignored) |

Sample type codes: `01` = Primary Solid Tumor (keep), `11` = Solid Tissue Normal (keep), `06` = Metastatic (exclude).

```r
meta <- colData(rse)
sample_type <- substr(meta$tcga.tcga_barcode, 14, 15)
keep_samples <- sample_type %in% c("01", "11")
counts <- counts[, keep_samples]
```
- `colData(rse)` extracts the sample metadata table.
- `substr(string, 14, 15)` extracts characters 14–15 from every barcode.
- `%in%` is a membership check — TRUE if `sample_type` is `"01"` or `"11"`.
- `[, keep_samples]` filters columns, keeping Tumor + Normal and dropping Metastatic.

**Result:** 1,135 Tumor (01) + 114 Normal (11) = 1,249 samples. 7 Metastatic (06) samples removed.

---

## Step 6 & 7 — Build colData & Create DESeq2 Object

**Step 6: Build the sample group table**
```r
condition <- ifelse(sample_type == "11", "Normal", "Tumor")
coldata <- data.frame(row.names = colnames(counts),
                       Condition = factor(condition, levels = c("Normal", "Tumor")))
```
- `ifelse()` is a vectorized if-else: `"11"` maps to `"Normal"`, everything else maps to `"Tumor"`.
- `factor(..., levels = c("Normal","Tumor"))` — level *order* matters: Normal first = reference group.

**Why Normal must come first:** DESeq2 computes `log2FC = second level / first level`. With Normal first and Tumor second, a positive log2FC means higher expression in Tumor. Reversing the order would flip every sign.

**Step 7: Create the DESeqDataSet object**
```r
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData   = coldata,
                               design    = ~ Condition)
```
- Bundles the count matrix and sample metadata into DESeq2's native object.
- `design = ~ Condition` — the formula operator `~` means "model expression as a function of Condition." A more complex model (e.g., with batch correction) would be `~ Batch + Condition`.

---

## Step 8 — Filter Genes | Step 9 — Run DESeq2

**Step 8: Remove low-expression genes**
```r
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
```
- `counts(dds) >= 10` gives TRUE/FALSE per cell; `rowSums()` counts TRUEs per gene; `>= 3` keeps genes with 10+ reads in at least 3 samples.
- Result: ~63k genes → ~39k genes. Genes with near-zero counts almost everywhere can't give reliable statistics and just slow things down.

**Step 9: Run DESeq2 (~10–20 minutes)**
```r
dds <- DESeq(dds)
```
One function performs, internally:
1. **Size factors** — estimate sequencing depth per sample to normalize for library-size differences.
2. **Dispersion** — estimate gene-wise variability; high-variance genes need stronger evidence to be called significant.
3. **NB GLM** — fit a Negative Binomial Generalized Linear Model to each gene's counts.
4. **Wald test** — test H0: log2FC = 0 for each gene, producing a p-value.
5. **BH correction** — Benjamini-Hochberg adjusts all p-values for multiple testing → `padj`.

With 1,249 samples, DESeq2 has enormous statistical power — even tiny 5% changes can reach `padj < 0.0001`. This is exactly why the `|log2FC| > 1` biological filter is added in Step 11.

---

## Step 10 & 11 — Extract Results & Apply the DEG Filter

**Step 10: Extract results table**
```r
res <- results(dds, contrast = c("Condition", "Tumor", "Normal"), alpha = 0.05)
```
`contrast = c(column, numerator, denominator)` → `log2FC = log2(Tumor/Normal)`, so positive = higher in Tumor. Output columns: `baseMean` (average expression), `log2FoldChange`, `lfcSE` (standard error), `stat`, `pvalue`, `padj` — one row per gene.

**Step 11: Apply the dual filter — the key step**
```r
res_df <- as.data.frame(res)
res_df$status <- "NS"
res_df$status[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange > 1]  <- "UP"
res_df$status[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "DOWN"
print(table(res_df$status))
```
- Every gene starts labeled `"NS"` (Not Significant).
- `!is.na(padj) & padj<0.05 & log2FC>1` → both conditions true → `"UP"`.
- Same logic with `log2FC < -1` → `"DOWN"`.
- `table(res_df$status)` gives the final UP/DOWN/NS counts — this is the DEG count you report.

---

## Step 12 — Save CSV | Step 13 — Volcano Plot

**Step 12: Save DEGs to CSV**
```r
res_ordered <- res_df[order(res_df$padj), ]
deg_only <- subset(res_ordered, status %in% c("UP", "DOWN"))
write.csv(deg_only, "results/DEG_results_TCGA_SIGNIFICANT_ONLY.csv", row.names = FALSE)
```
- `order()` sorts by `padj` ascending (most significant first).
- `subset(..., status %in% c("UP","DOWN"))` keeps only DEGs, dropping all `"NS"` genes.
- Output: 10,683 rows, columns `gene | log2FoldChange | padj | status`, sorted by significance.

**Step 13: Volcano plot**
```r
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("UP" = "red", "DOWN" = "blue", "NS" = "grey70")) +
  geom_text_repel(data = top10, aes(label = gene)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")
```
- x = fold change; y = `-log10(padj)`, which flips the scale so a *small* padj shows *high* on the plot; color by UP/DOWN/NS.
- `alpha = 0.5` makes points 50% transparent since tens of thousands of dots would otherwise overlap.
- `geom_text_repel()` labels only the top 10 most significant genes without overlapping text.
- Dashed lines mark the `|log2FC| > 1` and `padj = 0.05` boundaries.

Named "volcano" because the shape resembles an eruption: upper-right = UP DEGs (red), upper-left = DOWN DEGs (blue), bottom-center = NS (grey). Since `-log10(0.001) = 3` and `-log10(0.05) = 1.3`, a smaller padj sits higher on the plot — i.e., more significant.

---

## Step 14, 15, 16 — PCA, Heatmap & MA Plot

**Step 14 — PCA plot**
```r
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "Condition") + ggtitle("PCA: TCGA-BRCA Tumor vs Normal") + theme_bw()
ggsave("results/PCA_TCGA.png")
```
PCA collapses ~20,000 gene values per sample down to 2 numbers (PC1, PC2). A good result shows the Tumor cluster separating from Normal, with some overlap expected since breast tumors span multiple subtypes (Luminal A/B, HER2+, Basal).

**Step 15 — Heatmap**
```r
top_n <- 25
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col = anno, cluster_rows = TRUE, cluster_cols = FALSE)
```
Rows = genes, columns = samples. `mat - rowMeans(mat)` centers each gene around zero so the heatmap shows relative up/down patterns. `cluster_rows=TRUE` groups genes with similar expression patterns; `cluster_cols=FALSE` keeps samples grouped by Tumor/Normal rather than re-ordering them.

**Step 16 — MA plot**
```r
png("results/MA_TCGA.png")
plotMA(res, ylim = c(-5, 5), main = "MA Plot: TCGA-BRCA")
dev.off()
```
x = average expression across all samples (log scale); y = log2 fold change. Blue dots = significant genes. A healthy MA plot shows blue dots scattered above and below y = 0 in a funnel shape (more scatter at low expression, tighter at high expression).

**Note:** all plots use VST-transformed values for visualization, but DESeq2 always uses the raw integer counts internally for the actual statistics — VST is for human eyes only.

---

## Appendix — R Syntax Quick Reference

| Syntax | Meaning |
|---|---|
| `<-` | Assignment. `x <- 5` stores 5 in variable `x` |
| `#` / `##` | Comment — R ignores everything after `#` on a line |
| `$` | Column access, e.g. `meta$barcode` |
| `[ ]` | Indexing, `matrix[rows, cols]`. Blank = all. `[, 1:50]` = all rows, columns 1–50 |
| `::` | Package prefix, e.g. `BiocManager::install()` |
| `c()` | Combine into a vector, e.g. `c(1,2,3)` |
| `%in%` | Membership check, e.g. `x %in% c('01','11')` |
| `!is.na()` | Not missing — TRUE for non-NA values |
| `==` | Equality check (single `=` assigns instead) |
| `&` / `\|` | AND / OR (vectorized) |
| `paste0()` | Concatenate strings with no separator |
| `nrow()` / `ncol()` | Row count / column count |
| `head(x, n)` | First `n` elements/rows of `x` |
| `order(x)` | Indices that would sort `x` ascending |
| `rowMeans()` / `rowSums()` | Mean / sum of each row in a matrix |
| `substr(str, 14, 15)` | Extract characters 14–15 from a string |

---

## Pipeline Complete — Results Folder Contents

| File | Description |
|---|---|
| `DEG_results_TCGA_SIGNIFICANT_ONLY.csv` | 10,683 DEGs — `gene`, `log2FC`, `padj`, `UP/DOWN` status, sorted by padj |
| `DEG_results_TCGA.csv` | Full results table for all filtered genes (~39k rows), including non-significant |
| `boxplot_raw_TCGA.png` | Raw counts before normalization — QC baseline (lightblue boxes) |
| `boxplot_vst_TCGA.png` | VST counts after normalization — boxes should align (lightcoral) |
| `volcano_TCGA.png` | All filtered genes; UP = red, DOWN = blue, NS = grey; top genes labeled |
| `PCA_TCGA.png` | Tumor vs Normal separation check — quality control plot |
| `heatmap_top25_TCGA.png` | Top 25 DEGs, centered expression, samples ordered by condition |
| `heatmap_TCGA.png` | Top 30 DEGs with hierarchical clustering of both genes and samples |
| `MA_TCGA.png` | Average expression vs fold change — checks for normalization bias |
