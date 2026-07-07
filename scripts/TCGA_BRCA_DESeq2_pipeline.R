## ============================================================================
##  TCGA-BRCA RNA-seq Differential Expression Analysis — COMPLETE PIPELINE
##  Tumor vs Normal | recount3 + DESeq2
##
##  DATA TYPE: Transcriptomic data (bulk RNA-seq gene expression counts)
##             Source: TCGA-BRCA (The Cancer Genome Atlas - Breast Cancer)
##             Platform: Illumina RNA-seq, quantified as raw read counts
##             Access: via recount3 R package
##
##  KEY DEFINITIONS:
##    log2 Fold Change (log2FC): Measures how much a gene's expression
##      differs between Tumor and Normal on a log2 scale.
##      e.g. log2FC = 1 means 2x higher in Tumor; log2FC = -1 means 2x lower.
##    Adjusted p-value (padj): Statistical significance after Benjamini-
##      Hochberg correction for multiple testing across all genes.
##      padj < 0.05 means there is less than 5% chance the result is by chance.
##
##  DEG CRITERIA (both must be met):
##    (a) padj < 0.05       — statistically significant
##    (b) |log2FC| > 1      — at least 2-fold change (biologically meaningful)
##    Only genes meeting BOTH criteria are reported as DEGs.
##
##  HOW TO USE THIS SCRIPT (for someone new to R/programming):
##   1. Open RStudio.
##   2. Open this file (File > Open File > select this .R file).
##   3. Set your working directory to a folder where you want results saved.
##        You can do this via: Session > Set Working Directory > Choose Directory
##   4. Run the WHOLE script: click "Source" button (top-right of editor),
##        OR select all (Ctrl+A) and press Ctrl+Enter.
##   5. The script will take 15-30 minutes total (download + DESeq2 run).
##        Do not close RStudio while it is running.
##   6. When finished, a folder called "results" will appear in your
##        working directory, containing the CSV file and 4 plots (PNG images).
##
##  WHAT THIS SCRIPT DOES (big picture):
##   - Downloads real breast cancer gene expression data (TCGA-BRCA) using recount3
##   - Separates Tumor samples from Normal (healthy) samples
##   - Runs DESeq2 to find which genes are significantly different
##     between Tumor and Normal
##   - Applies BOTH standard filters used in published papers:
##        (a) padj < 0.05  (statistically significant)
##        (b) |log2FoldChange| > 1  (at least 2-fold change — biologically meaningful)
##   - Saves a results table (CSV) and 6 plots (boxplot before/after + volcano + PCA + heatmap + MA)
## ============================================================================


## ----------------------------------------------------------------------
## STEP 0: INSTALL REQUIRED PACKAGES (run this section ONLY ONCE EVER)
## ----------------------------------------------------------------------
## If you have already installed these before, you can skip this block.
## It is safe to run again — it will just skip anything already installed.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("recount3", update = FALSE, ask = FALSE)   # gives us TCGA data
BiocManager::install("DESeq2",   update = FALSE, ask = FALSE)   # for differential expression
install.packages("ggplot2",  repos = "https://cran.rstudio.com")   # for volcano + PCA plots
install.packages("pheatmap", repos = "https://cran.rstudio.com")   # for heatmap
install.packages("ggrepel",  repos = "https://cran.rstudio.com")   # for gene labels on volcano


## ----------------------------------------------------------------------
## STEP 1: LOAD LIBRARIES (run this every time you start a new R session)
## ----------------------------------------------------------------------

library(recount3)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)

cat("All packages loaded successfully!\n")


## ----------------------------------------------------------------------
## STEP 2: FIND THE TCGA BREAST CANCER (BRCA) PROJECT IN recount3
## ----------------------------------------------------------------------
## available_projects() lists every dataset recount3 has.
## We filter ("subset") it down to just the BRCA (breast cancer) project
## from the TCGA data source.

cat("Looking up available projects... (this may take 1-2 minutes)\n")

human_projects <- available_projects()

# If RStudio asks: "...does not exist, create directory? (yes/no):"
# type:  yes
# and press Enter. This is just recount3 creating a folder to cache
# downloaded files so future runs are faster.

brca_project <- subset(human_projects,
                        project      == "BRCA" &
                        project_home == "data_sources/tcga")

print(brca_project)
# Expected: project = BRCA, n_samples = 1256


## ----------------------------------------------------------------------
## STEP 3: DOWNLOAD THE RAW GENE EXPRESSION DATA
## ----------------------------------------------------------------------
## This downloads ~1,256 samples x 63,856 genes worth of count data.
## Takes 5-15 minutes depending on your internet connection.

cat("Downloading TCGA-BRCA data... please wait (5-15 mins)\n")

rse <- create_rse(brca_project)

cat("Download complete!\n")
cat("Dimensions: ", nrow(rse), "genes x", ncol(rse), "samples\n")


## ----------------------------------------------------------------------
## STEP 4: CONVERT TO TRUE RAW INTEGER COUNTS
## ----------------------------------------------------------------------
## recount3 stores data in a special "coverage" format internally.
## transform_counts() converts it into normal gene-level read counts
## (whole numbers / integers), which is exactly what DESeq2 requires.

assay(rse, "counts") <- transform_counts(rse)
counts <- assay(rse, "counts")

cat("Genes:", nrow(counts), "| Samples:", ncol(counts), "\n")

# Sanity check: values must be whole numbers (0, 1, 2, 847...)
# NOT decimals (0.93, 1.47...). If you see decimals, something is wrong.
cat("First few values (should be whole numbers):\n")
print(counts[1:3, 1:3])


## ----------------------------------------------------------------------
## STEP 4b: BOXPLOT OF RAW COUNTS — BEFORE NORMALIZATION (QC CHECK)
## ----------------------------------------------------------------------
## We plot log2(raw counts + 1) across a sample of 50 columns BEFORE
## running DESeq2. Each box = one sample. We expect boxes at different
## heights (normal — different sequencing depths) but with roughly
## similar shapes. A dramatically different box may signal a QC issue.
##
## The +1 avoids log2(0) = -Infinity errors for genes with zero counts.
## We plot only 50 columns to keep the figure readable — with 1,249
## samples, showing all at once would be unreadable.

dir.create("results", showWarnings = FALSE)

log_counts_raw <- log2(counts + 1)

png("results/boxplot_raw_TCGA.png", width = 1200, height = 600)
boxplot(log_counts_raw[, 1:50],
        las     = 2,
        cex.axis = 0.5,
        main    = "Raw Counts (log2) — Before Normalization",
        ylab    = "log2(counts + 1)",
        col     = "lightblue")
dev.off()

cat("Raw counts boxplot saved: results/boxplot_raw_TCGA.png\n")


## ----------------------------------------------------------------------
## STEP 5: SEPARATE TUMOR AND NORMAL SAMPLES
## ----------------------------------------------------------------------
## Every TCGA sample barcode encodes a "sample type" in characters 14-15:
##    01 = Primary Solid Tumor   (cancer tissue)
##    06 = Metastatic Tumor      (cancer that has spread -- we EXCLUDE this)
##    11 = Solid Tissue Normal   (healthy tissue)
##
## We keep ONLY 01 (Tumor) and 11 (Normal) for a clean two-group comparison.

meta <- colData(rse)
sample_type <- substr(meta$tcga.tcga_barcode, 14, 15)

cat("All sample types found in the data:\n")
print(table(sample_type))

# Keep only Tumor (01) and Normal (11); drop everything else (e.g. 06)
keep_samples <- sample_type %in% c("01", "11")
counts      <- counts[, keep_samples]
meta        <- meta[keep_samples, ]
sample_type <- sample_type[keep_samples]

cat("\nAfter filtering:\n")
cat("  Tumor samples :", sum(sample_type == "01"), "\n")
cat("  Normal samples:", sum(sample_type == "11"), "\n")


## ----------------------------------------------------------------------
## STEP 6: BUILD THE SAMPLE INFORMATION TABLE (for DESeq2)
## ----------------------------------------------------------------------
## DESeq2 needs an explicit table telling it which sample is Tumor and
## which is Normal. "Normal" is set as the reference level, so positive
## fold changes mean "higher in Tumor".

condition <- ifelse(sample_type == "11", "Normal", "Tumor")

coldata <- data.frame(
  row.names = colnames(counts),
  Condition = factor(condition, levels = c("Normal", "Tumor"))
)

cat("Sample group sizes:\n")
print(table(coldata$Condition))


## ----------------------------------------------------------------------
## STEP 7: CREATE THE DESeq2 OBJECT
## ----------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ Condition
)


## ----------------------------------------------------------------------
## STEP 8: FILTER OUT VERY LOW-EXPRESSION GENES (removes noise)
## ----------------------------------------------------------------------
## Keep a gene only if at least 3 samples have 10+ reads for it.
## This removes genes that are essentially "silent" everywhere and
## cannot give us reliable statistics.

keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]

cat("Genes remaining after filtering:", nrow(dds), "\n")


## ----------------------------------------------------------------------
## STEP 9: RUN DESeq2 — THE CORE STATISTICAL ANALYSIS
## ----------------------------------------------------------------------
## This single function does normalisation, models each gene's
## variability, and statistically tests Tumor vs Normal for every gene.
## With ~1,250 samples this takes roughly 10-20 minutes. Be patient!

cat("Running DESeq2... this will take 10-20 minutes. Please wait...\n")

dds <- DESeq(dds)

cat("DESeq2 finished!\n")


## ----------------------------------------------------------------------
## STEP 9b: BOXPLOT OF VST-NORMALIZED COUNTS — AFTER NORMALIZATION (QC)
## ----------------------------------------------------------------------
## After DESeq2, we apply VST (Variance Stabilizing Transform) to get
## normalized, log-scale values for visualization. blind = FALSE tells
## VST to use the known group structure (Tumor / Normal) when estimating
## dispersion — this is recommended after DESeq2 has already been run.
##
## What we now expect compared to the BEFORE plot:
##   - Boxes at roughly the same height: normalization corrected for
##     sequencing depth differences.
##   - Similar medians across all samples: if some samples still look
##     very different, there may be a batch effect or QC issue.
##
## Comparing BEFORE (Step 4b) and AFTER (this step) side by side is the
## standard way to confirm that normalization worked correctly.

vsd_qc <- vst(dds, blind = FALSE)

png("results/boxplot_vst_TCGA.png", width = 1200, height = 600)
boxplot(assay(vsd_qc)[, 1:50],
        las     = 2,
        cex.axis = 0.5,
        main    = "VST-Normalized Counts — After Normalization",
        ylab    = "VST expression",
        col     = "lightcoral")
dev.off()

cat("VST normalized boxplot saved: results/boxplot_vst_TCGA.png\n")


## ----------------------------------------------------------------------
## STEP 10: EXTRACT RESULTS
## ----------------------------------------------------------------------
## contrast = c("Condition", "Tumor", "Normal") means:
##     log2FoldChange = expression in Tumor  vs  expression in Normal
##     Positive value = higher in Tumor.   Negative value = lower in Tumor.

res <- results(dds,
               contrast = c("Condition", "Tumor", "Normal"),
               alpha    = 0.05)

cat("\n--- Raw summary (padj < 0.05 only, NO fold-change filter yet) ---\n")

summary(res)

## NOTE: This summary alone is NOT the final DEG count!
## With ~1,250 samples, DESeq2 has huge statistical power, so even
## TINY (biologically meaningless) fold changes become "significant".
## We MUST also apply a fold-change filter -- done in Step 11 below.


## ----------------------------------------------------------------------
## STEP 11: APPLY THE STANDARD TWO-PART DEG FILTER  *** (THE KEY FIX) ***
## ----------------------------------------------------------------------
## A gene is only called a real "DEG" (Differentially Expressed Gene) if
## BOTH of these are true (this is the standard used in published papers):
##
##    (a) padj < 0.05            --> statistically significant
##    (b) |log2FoldChange| > 1   --> at least a 2-fold change (biologically meaningful)
##
## Genes meeting both criteria are labelled UP or DOWN.
## Everything else is labelled NS ("Not Significant" / not a real DEG).

res_df <- as.data.frame(res)
res_df$gene   <- rownames(res_df)

# Start everyone as "NS" (not significant)
res_df$status <- "NS"

# Mark UP only if BOTH conditions are met
res_df$status[ !is.na(res_df$padj) &
               res_df$padj < 0.05 &
               res_df$log2FoldChange >  1 ] <- "UP"

# Mark DOWN only if BOTH conditions are met
res_df$status[ !is.na(res_df$padj) &
               res_df$padj < 0.05 &
               res_df$log2FoldChange < -1 ] <- "DOWN"

cat("\n--- FINAL DEG counts (padj < 0.05 AND |log2FC| > 1) ---\n")
print(table(res_df$status))
## This is the number you should actually report in your report/poster,
## NOT the raw padj-only summary from Step 10.


## ----------------------------------------------------------------------
## STEP 12: SAVE RESULTS TO CSV
## ----------------------------------------------------------------------
## Per reporting standards, we report ONLY significant DEGs (UP + DOWN).
## Columns reported: gene, log2FoldChange, padj, status
## NS (non-significant) genes are excluded from the output.

# Order by significance (smallest padj = most significant, first)
res_ordered <- res_df[order(res_df$padj), ]

dir.create("results", showWarnings = FALSE)

# Keep ONLY significant DEGs, and ONLY the key columns
deg_only <- subset(res_ordered, status %in% c("UP", "DOWN"))
deg_report <- deg_only[, c("gene", "log2FoldChange", "padj", "status")]

write.csv(deg_report, "results/DEG_results_TCGA_SIGNIFICANT_ONLY.csv", row.names = FALSE)

cat("\nSignificant DEGs only saved to:\n")
cat("  results/DEG_results_TCGA_SIGNIFICANT_ONLY.csv\n")
cat("  Columns: gene | log2FoldChange | padj | status\n")
cat("  (", nrow(deg_report), "genes total: ",
    sum(deg_report$status == "UP"),   "UP, ",
    sum(deg_report$status == "DOWN"), "DOWN )\n")


## ----------------------------------------------------------------------
## STEP 13: VOLCANO PLOT
## ----------------------------------------------------------------------
## Every gene = one dot.
##   x-axis = log2 fold change (left = down in Tumor, right = up in Tumor)
##   y-axis = -log10(padj)     (higher = more statistically significant)
##   Red    = significantly UP in Tumor
##   Blue   = significantly DOWN in Tumor
##   Grey   = not significant / not a real DEG

# Pick top 10 most significant genes (by padj) to label on the plot
top10 <- head(res_ordered[order(res_ordered$padj), ], 10)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("UP" = "red", "DOWN" = "blue", "NS" = "grey70")) +
  geom_text_repel(data = top10, aes(label = gene), size = 3, max.overlaps = 15) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volcano Plot: TCGA-BRCA Tumor vs Normal",
       subtitle = paste0("UP = ", sum(res_df$status == "UP"),
                          " | DOWN = ", sum(res_df$status == "DOWN"),
                          " | NS = ", sum(res_df$status == "NS")),
       x = "log2 Fold Change", y = "-log10(adjusted p-value)") +
  theme_bw()

ggsave("results/volcano_TCGA.png", width = 8, height = 6, dpi = 150)
cat("Volcano plot saved!\n")


## ----------------------------------------------------------------------
## STEP 14: PCA PLOT (Quality Control)
## ----------------------------------------------------------------------
## Reduces all gene expression data to 2 axes (PC1, PC2) to visually
## check whether Tumor and Normal samples separate as expected.
##
## NOTE: TCGA-BRCA tumors include several molecular subtypes (Luminal A/B,
## HER2+, Basal/TNBC), so the Tumor group is biologically heterogeneous
## and may spread widely / partially overlap with Normal in PCA. This is
## EXPECTED and consistent with published TCGA-BRCA analyses -- it does
## NOT mean the analysis is wrong.

vsd <- vst(dds, blind = FALSE)
# VST = Variance Stabilising Transform, blind = FALSE uses known group
# structure (Tumor / Normal) for dispersion estimation — recommended
# when groups are known. Used ONLY for PCA/heatmap visualisation;
# DESeq2 itself always uses the original raw counts, never VST data.

plotPCA(vsd, intgroup = "Condition") +
  ggtitle("PCA: TCGA-BRCA Tumor vs Normal") +
  theme_bw()

ggsave("results/PCA_TCGA.png", width = 7, height = 5, dpi = 150)
cat("PCA plot saved!\n")


## ----------------------------------------------------------------------
## STEP 15: HEATMAP OF TOP DEGs
## ----------------------------------------------------------------------
## Shows expression patterns of the most significant genes across all samples.
## Each row = one gene (y-axis only, gene names shown).
## Each column = one sample (x-axis labels REMOVED to avoid clutter).
## Red = higher than average expression, Blue = lower than average.
##
## AVAILABLE OPTIONS — change the number below to switch:
##   top_n <- 10    (Top 10 most significant genes)
##   top_n <- 20    (Top 20)
##   top_n <- 25    (Top 25)
##   top_n <- 50    (Top 50)
##
## Currently set to 25 as default. Change top_n to any of the above.

top_n <- 25   # <--- CHANGE THIS to 10, 20, 25, or 50 as needed

top_genes <- head(rownames(res_ordered[!is.na(res_ordered$padj), ]), top_n)

mat  <- assay(vsd)[top_genes, ]
mat  <- mat - rowMeans(mat)   # centre each gene around its own average

anno <- as.data.frame(colData(vsd)[, "Condition", drop = FALSE])

pheatmap(mat,
         annotation_col  = anno,
         show_rownames   = TRUE,    # show gene names on y-axis
         show_colnames   = FALSE,   # NO x-axis sample labels (too cluttered)
         cluster_cols    = FALSE,   # NO column dendrogram (removed as requested)
         cluster_rows    = TRUE,    # keep row dendrogram (gene clustering on y-axis)
         fontsize_row    = 7,
         main            = paste0("Top ", top_n, " DEGs - TCGA-BRCA"),
         filename        = paste0("results/heatmap_top", top_n, "_TCGA.png"),
         width = 8, height = 9)

cat("Heatmap saved! (Top", top_n, "DEGs)\n")
cat("  File: results/heatmap_top", top_n, "_TCGA.png\n")
cat("  To change: set top_n <- 10, 20, 25, or 50 at the start of Step 15\n")


## ----------------------------------------------------------------------
## STEP 16: MA PLOT (Technical Quality Check)
## ----------------------------------------------------------------------
## x-axis = average expression level of each gene
## y-axis = log fold change
## Blue dots = statistically significant genes (padj < 0.05)
## A healthy MA plot has blue dots spread both above and below y = 0,
## with no obvious one-sided bias.

png("results/MA_TCGA.png", width = 800, height = 600)
plotMA(res, ylim = c(-5, 5), main = "MA Plot: TCGA-BRCA Tumor vs Normal")
dev.off()

cat("MA plot saved!\n")


## ----------------------------------------------------------------------
## DONE!
## ----------------------------------------------------------------------

cat("\n=====================================================\n")
cat(" ALL DONE! Check the 'results' folder for:\n")
cat("   - DEG_results_TCGA_SIGNIFICANT_ONLY.csv\n")
cat("       (significant DEGs only: gene | log2FoldChange | padj | status)\n")
cat("   - boxplot_raw_TCGA.png    (raw counts BEFORE normalization — QC check)\n")
cat("   - boxplot_vst_TCGA.png    (VST counts AFTER normalization — QC confirmation)\n")
cat("   - volcano_TCGA.png\n")
cat("   - PCA_TCGA.png\n")
cat("   - heatmap_top", top_n, "_TCGA.png  (change top_n in Step 15 for other sizes)\n")
cat("   - MA_TCGA.png\n")
cat("=====================================================\n")

cat("\nFINAL DEG COUNT (padj < 0.05 AND |log2FC| > 1):\n")
print(table(res_df$status))
