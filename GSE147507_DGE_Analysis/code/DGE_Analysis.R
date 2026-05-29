
#Differential Gene Expression Analysis using DESeq2

#===============================================================================
# 1. PACKAGE INSTALLATION & LOADING
#===============================================================================

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("Biobase", quietly = TRUE)) BiocManager::install("Biobase")
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("ashr", quietly = TRUE)) BiocManager::install("ashr")
if (!requireNamespace("pheatmap", quietly = TRUE)) BiocManager::install("pheatmap")
if (!requireNamespace("ggrepel", quietly = TRUE)) BiocManager::install("ggrepel")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")

library("BiocManager")
library("GEOquery")
library("dplyr")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("ggrepel")
library("DESeq2")
library("readr")
library("stringr")
library("here")

#===============================================================================
# 2. DATA LOADING & CLEANING
#===============================================================================

# Expression Data
COUNTS_raw <- read_tsv(here("datasets", "GSE147507_RawReadCounts_Human.tsv"), show_col_types = FALSE)
COUNTS <- (as.data.frame(COUNTS_raw))
rownames(COUNTS) <- COUNTS$...1
COUNTS <- COUNTS[, -1]
COUNTS <- round(COUNTS)
COUNTS <- as.matrix(COUNTS)
summary(COUNTS[, 1:3])
cat("Missing values in expression data:", sum(is.na(COUNTS)), "\n")

# Metadata
gse_parsed <- getGEO(filename = here("datasets", "GSE147507-GPL18573_series_matrix.txt"))
META <- pData(gse_parsed)
unique(META$characteristics_ch1)
META <- META %>%
  mutate(`cell line:ch1` = case_when(
    str_detect(title, "NHBE") ~ "NHBE",
    str_detect(title, "A549-ACE2") ~ "A549_ACE2",
    str_detect(title, "A549") ~ "A549",
    str_detect(title, "Calu3") ~ "Calu3",
    TRUE ~ "Other"
  ),
  `treatment:ch1`= case_when(
    str_detect(title, "SARS-CoV-2") ~ "SARS_CoV_2",
    str_detect(title, "Mock") ~ "Mock",
    str_detect(title, "IAV") ~ "IAV",
    str_detect(title, "RSV") ~ "RSV",
    TRUE ~ "Other"
  )
  ) %>%
  mutate(group = as.factor(paste(`cell line:ch1`, `treatment:ch1`, sep = "_")))
unique(META$group)
META <- as.data.frame(META)
rownames(META) <- META$title
META <- META[colnames(COUNTS), ]
if (all(colnames(COUNTS) == rownames(META))) {cat("Success! Names and order match perfectly.\n")} else {stop("Mismatch still exists. Double-check your column names.")}

# ==============================================================================
# 3. CELL LINE CONFIGURATION
# ==============================================================================

cell_line_config <- list(
  list(
    id            = "NHBE",
    display_name  = "NHBE",
    treated_group = "NHBE_SARS_CoV_2",
    control_group = "NHBE_Mock",
    dist_palette  = "Blues",
    heat_palette  = "RdYlBu"
  ),
  list(
    id            = "A549",
    display_name  = "A549",
    treated_group = "A549_SARS_CoV_2",
    control_group = "A549_Mock",
    dist_palette  = "YlOrRd",
    heat_palette  = "PiYG"
  ),
  list(
    id            = "A549_ACE2",
    display_name  = "A549-ACE2",
    treated_group = "A549_ACE2_SARS_CoV_2",
    control_group = "A549_ACE2_Mock",
    dist_palette  = "Greens",
    heat_palette  = "BrBG"
  ),
  list(
    id            = "Calu3",
    display_name  = "Calu-3",
    treated_group = "Calu3_SARS_CoV_2",
    control_group = "Calu3_Mock",
    dist_palette  = "Purples",
    heat_palette  = "PRGn"
  )
)

#===============================================================================
# 4. DESeq2 IMPLEMENTATION
#===============================================================================

# Creating DESeq2 Dataset
dds <- DESeqDataSetFromMatrix(countData = COUNTS, colData = META, design = ~ group)
cat("Dimensions of before Filtering:", dim(dds), "\n")

# Filtering low count Genes
threshold <- 10
dds <- dds[ rowMeans(counts(dds)) >= threshold,]
cat("Dimensions after Filtering:", dim(dds), "\n")

# DESeq2 Analysis
set.seed(42)
prdds <- DESeq(dds)

# Normalization
norm_counts <- counts(prdds, normalized = TRUE)
norm_counts <- as.data.frame(norm_counts)

# Transformation
vsd <- varianceStabilizingTransformation(prdds, blind = FALSE)

# ==============================================================================
# 6. QC VISUALIZATIONS
# ==============================================================================

# Scatter Plots: log2 vs VST
par(mfrow=c(1, 2))
lims <- c(-2, 20)
plot(log2(counts(prdds, normalized = TRUE)[,1:2] + 1), pch=16, cex=0.3, main="log2(x + 1)", xlim=lims, ylim=lims)
plot(assay(vsd)[,1:2], pch=16, cex=0.3, main="VST", xlim=lims, ylim=lims)

# Histograms: raw vs VST
par(mfrow=c(1, 2))
hist(counts(prdds))
hist(assay(vsd))

# Sample-to-sample distance heatmaps (one per cell line)
plot_sample_dist <- function(vsd, META, cell_line_id, display_name, dist_palette) {
  samples <- rownames(META[META$`cell line:ch1` == cell_line_id, ])
  vsd_sub <- vsd[, samples]
  sample_dist <- dist(t(assay(vsd_sub)))
  dist_matrix <- as.matrix(sample_dist)
  colors <- colorRampPalette(rev(brewer.pal(9, dist_palette)))(255)
  pheatmap(dist_matrix,
         clustering_distance_rows = sample_dist,
         clustering_distance_cols = sample_dist,
         col = colors,
         fontsize_row = 8,
         fontsize_col = 8,
         main = paste0("Sample Distance Matrix (", display_name," Cells Only)")
  )
}
for (cfg in cell_line_config) {
  plot_sample_dist(vsd, META, cfg$id, cfg$display_name, cfg$dist_palette)
}

#PCA Plot
pca_data <- plotPCA(vsd, intgroup = c("group"), returnData = TRUE)
ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, aes(color = group)) +
  xlab(paste0("PC1: ", round(attr(pca_data, "percentVar")[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(attr(pca_data, "percentVar")[2] * 100), "% variance")) +
  theme_minimal() +
  ggtitle("PCA by Cell Line and Treatment")

#Dispersion Plot
par(mfrow=c(1, 1))
plotDispEsts(prdds, main = "Dispersion Plot",
             genecol="gray20", fitcol="red", finalcol="dodgerblue3")

#===============================================================================
# 7. DIFFERENTIAL GENE EXPRESSION (DGE) ANALYSIS & VISUALIZATIONS PER-CELL-LINE
#===============================================================================

# Generating DESeq2 Results, Applying Shrinkage & Ordering Data by adjusted p-value
get_dge_results <- function(prdds, treated, control) {
  res_bs <- results(prdds, contrast = c("group", treated, control), alpha = 0.05, lfcThreshold = 1)
  res    <- lfcShrink(prdds, contrast = c("group", treated, control), res = res_bs, type = "ashr")
  res    <- na.omit(res)
  res    <- res[order(res$padj), ]
  return(res)
}

res <- lapply(cell_line_config, function(cfg) {
  get_dge_results(prdds, cfg$treated_group, cfg$control_group)
})

names(res) <- sapply(cell_line_config, function(cfg) cfg$id)

# MA Plots Comparison
par(mfrow=c(2, 2))
for (cfg in cell_line_config) {
  DESeq2::plotMA(
    res[[cfg$id]],
              main = paste0("Differential Expression: ", cfg$display_name," - Sars-CoV-2 vs Mock"),
              ylim=c(-10,10), cex=0.5,
              colNonSig=adjustcolor("gray20", alpha.f = 0.5),
              colSig=adjustcolor("dodgerblue3", alpha.f = 0.5))
  abline(h = c(-1, 1), col = "#ff0000", lwd = 1, lty = 2)
  }

# Assigning Gene Status
for (cfg in cell_line_config) {
  tmp <- as.data.frame(res[[cfg$id]])
  tmp$status <- ifelse(tmp$padj < 0.05 & abs(tmp$log2FoldChange) > 1,
                       ifelse(tmp$log2FoldChange > 0, "Upregulated", "Downregulated"),
                       "Non-Significant")
  res[[cfg$id]] <- tmp 
}

# Generating Volcano Plots
volcano_colors <- c("Upregulated" = "#4DAF4A", "Downregulated" = "#E41A1C", "Non-Significant" = "#377EB8")
for (cfg in cell_line_config) {
  tmp <- as.data.frame(res[[cfg$id]])
  top_labels <- head(tmp, 15)
  print(
    ggplot(as.data.frame(tmp),
      aes(x = log2FoldChange, y = -log10(padj), color = status)) +
      geom_point(size = 1, alpha = 0.7) +
      scale_color_manual(values = volcano_colors) +
      theme_minimal() +
      ggtitle(paste0("Volcano Plot of Differentially Expressed Genes: ", cfg$display_name," - SARS-CoV-2 vs Mock")) +
      xlab("log2 Fold Change") +
      ylab("-log10 (Adjusted p-value)") +
      theme(legend.title = element_blank()) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
      geom_text_repel(data = as.data.frame(top_labels),
                      aes(label = rownames(top_labels)),
                      size = 2.5, max.overlaps = 20)
  )
  res[[cfg$id]] <- tmp 
}

# Extracting Significant Genes
vsd_subsets <- list()
sig_genes   <- list()
up_genes    <- list()
down_genes  <- list()

for (cfg in cell_line_config) {
    samples_to_keep <- (rownames(META[META$group %in% c(cfg$treated_group, cfg$control_group), ]))
    vsd_subsets[[cfg$id]] <- (vsd[, samples_to_keep])
    
    sig <- as.data.frame(res[[cfg$id]][res[[cfg$id]]$padj < 0.05 & abs(res[[cfg$id]]$log2FoldChange) > 1, ])
    sig_genes[[cfg$id]] <- sig
    up_genes[[cfg$id]] <- subset(sig, log2FoldChange > 0)
    down_genes[[cfg$id]] <- subset(sig, log2FoldChange < 0)
    }

# Creating Heatmaps of the Top 10 Up-regulated and Down-Regulated Genes
for (cfg in cell_line_config) {
  heatmap_colors <- colorRampPalette(rev(brewer.pal(name = cfg$heat_palette, n = 11)))(255)
  top_up <- head(up_genes[[cfg$id]], 10)
  if (nrow(top_up) > 1) {
    top_up_exp <- assay(vsd_subsets[[cfg$id]])[rownames(top_up), ]
    pheatmap(top_up_exp,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           scale = "row",
           show_colnames = TRUE,
           col=heatmap_colors,
           main = paste0("Top 10 Up-Regulated ", cfg$display_name," Genes: SARS-CoV-2 vs Mock"))
  }
}

for (cfg in cell_line_config) {
  heatmap_colors <- colorRampPalette(rev(brewer.pal(name = cfg$heat_palette, n = 11)))(255)
  top_down <- head(down_genes[[cfg$id]], 10)
  if (nrow(top_down) > 1) {
    top_down_exp <- assay(vsd_subsets[[cfg$id]])[rownames(top_down), ]
    pheatmap(top_down_exp,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             scale = "row",
             show_colnames = TRUE,
             col=heatmap_colors,
             main = paste0("Top 10 Down-Regulated ", cfg$display_name," Genes: SARS-CoV-2 vs Mock"))
  }
}

# Creating DGE Analysis Files
if (!file.exists("results")) {dir.create("results")}

for (cfg in cell_line_config) {
  sig_genes[[cfg$id]]$status <- ifelse(sig_genes[[cfg$id]]$log2FoldChange > 0, "Upregulated", "Downregulated")
  sig_genes[[cfg$id]]$gene_id <- rownames(sig_genes[[cfg$id]])
  sig_genes[[cfg$id]] <- sig_genes[[cfg$id]][, c("gene_id","baseMean","log2FoldChange","lfcSE","pvalue","padj","status")]
  write.table(sig_genes[[cfg$id]], file = paste0("results/significant_DGE_across_", cfg$display_name,".csv"), sep = ",", row.names = FALSE, quote = FALSE)
}

#===============================================================================
# 8. SESSION INFO
#===============================================================================

sessionInfo()
# R version 4.6.0 (2026-04-24 ucrt)
# Bioconductor version 3.23 (BiocManager 1.30.27)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 26200)
# Matrix products: default
# LAPACK version 3.12.1

#===============================================================================
