
#Differential Gene Expression Analysis using DESeq2 

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(c("GEOquery","Biobase"))

if (!requireNamespace("DESeq2", quietly = TRUE))
BiocManager::install("DESeq2")

if (!requireNamespace("ashr", quietly = TRUE))
BiocManager::install("ashr")

if (!requireNamespace("pheatmap", quietly = TRUE))
BiocManager::install("pheatmap")

if (!requireNamespace("ggrepel", quietly = TRUE))
BiocManager::install("ggrepel")

if (!requireNamespace("readr", quietly = TRUE))
install.packages("readr")

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

setwd("path/to/GSE147507_DGE_Analysis")

#Extracting & Cleaning Expression Data
COUNTS_raw <- read_tsv("datasets/GSE147507_RawReadCounts_Human.tsv", show_col_types = FALSE)
COUNTS <- (as.data.frame(COUNTS_raw))
rownames(COUNTS) <- COUNTS$...1
COUNTS <- COUNTS[, -1]

#Rounding off & Converting to Matrix
COUNTS <- round(COUNTS)
COUNTS <- as.matrix(COUNTS)
summary(COUNTS[, 1:3])
cat("Missing values in expression data:", sum(is.na(COUNTS)), "\n")

#Extracting Metadata
gse_parsed <- getGEO(filename = "datasets/GSE147507-GPL18573_series_matrix.txt")
META <- pData(gse_parsed)
unique(META$characteristics_ch1)

#Cleaning Metadata
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

#Listing Unique Sample Conditions
unique(META$group)

#Wrangling Metadata
META <- as.data.frame(META)
rownames(META) <- META$title
META <- META[colnames(COUNTS), ]
if (all(colnames(COUNTS) == rownames(META))) {cat("Success! Names and order match perfectly.\n")} else {stop("Mismatch still exists. Double-check your column names.")}

#Creating DESeq2 Dataset
dds <- DESeqDataSetFromMatrix(countData = COUNTS, colData = META, design = ~ group)
cat("Dimensions of before Filtering:", dim(dds), "\n")

#Filtering low count Genes
threshold <- 10
dds <- dds[ rowMeans(counts(dds)) >= threshold,]
cat("Dimensions after Filtering:", dim(dds), "\n")

#DESeq2 Analysis
set.seed(42)
prdds <- DESeq(dds)

#Normalization
norm_counts <- counts(prdds, normalized = TRUE)
norm_counts <- as.data.frame(norm_counts)

#Transformation
vsd <- varianceStabilizingTransformation(prdds, blind = FALSE)

#Scatter Plots Comparison
par(mfrow=c(1, 2))
lims <- c(-2, 20)
plot(log2(counts(prdds, normalized = TRUE)[,1:2] + 1), pch=16, cex=0.3, main="log2(x + 1)", xlim=lims, ylim=lims)
plot(assay(vsd)[,1:2], pch=16, cex=0.3, main="VST", xlim=lims, ylim=lims)

#Histograms Comparison
par(mfrow=c(1, 2))
hist(counts(prdds))
hist(assay(vsd))

#Sample-to-Sample Distances (Heatmaps)
nhbe_samples <- rownames(META[META$`cell line:ch1` == "NHBE", ])
vsd_nhbe <- vsd[, nhbe_samples]
sample_dist_nhbe <- dist(t(assay(vsd_nhbe)))
sample_dist_matrix_nhbe <- as.matrix(sample_dist_nhbe)
colors_nhbe <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sample_dist_matrix_nhbe,
         clustering_distance_rows = sample_dist_nhbe,
         clustering_distance_cols = sample_dist_nhbe,
         col = colors_nhbe,
         fontsize_row = 8,
         fontsize_col = 8,
         main = "Sample Distance Matrix (NHBE Cells Only)")

a549_samples <- rownames(META[META$`cell line:ch1` == "A549", ])
vsd_a549 <- vsd[, a549_samples]
sample_dist_a549 <- dist(t(assay(vsd_a549)))
sample_dist_matrix_a549 <- as.matrix(sample_dist_a549)
colors_a549 <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")))(255)
pheatmap(sample_dist_matrix_a549,
         clustering_distance_rows = sample_dist_a549,
         clustering_distance_cols = sample_dist_a549,
         col = colors_a549,
         fontsize_row = 8,
         fontsize_col = 8,
         main = "Sample Distance Matrix (A549 Cells Only)")

a549_ace2_samples <- rownames(META[META$`cell line:ch1` == "A549_ACE2", ])
vsd_a549_ace2 <- vsd[, a549_ace2_samples]
sample_dist_a549_ace2 <- dist(t(assay(vsd_a549_ace2)))
sample_dist_matrix_a549_ace2 <- as.matrix(sample_dist_a549_ace2)
colors_a549_ace2 <- colorRampPalette(rev(brewer.pal(9, "Greens")))(255)
pheatmap(sample_dist_matrix_a549_ace2,
         clustering_distance_rows = sample_dist_a549_ace2,
         clustering_distance_cols = sample_dist_a549_ace2,
         col = colors_a549_ace2,
         fontsize_row = 8,
         fontsize_col = 8,
         main = "Sample Distance Matrix (A549-ACE2 Cells Only)")

calu3_samples <- rownames(META[META$`cell line:ch1` == "Calu3", ])
vsd_calu3 <- vsd[, calu3_samples]
sample_dist_calu3 <- dist(t(assay(vsd_calu3)))
sample_dist_matrix_calu3 <- as.matrix(sample_dist_calu3)
colors_calu3 <- colorRampPalette(rev(brewer.pal(9, "Purples")))(255)
pheatmap(sample_dist_matrix_calu3,
         clustering_distance_rows = sample_dist_calu3,
         clustering_distance_cols = sample_dist_calu3,
         col = colors_calu3,
         fontsize_row = 8,
         fontsize_col = 8,
         main = "Sample Distance Matrix (Calu-3 Cells Only)")

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

#Generating DESeq2 Results & Ordering Data by adjusted p-value
res_nhbe_bs <- results(prdds, contrast = c("group", "NHBE_SARS_CoV_2", "NHBE_Mock"), alpha = 0.05, lfcThreshold = 1)
res_nhbe <- lfcShrink(prdds, contrast = c("group", "NHBE_SARS_CoV_2", "NHBE_Mock"), res = res_nhbe_bs, type = "ashr")
res_nhbe <- na.omit(res_nhbe)
res_nhbe_ordered <- res_nhbe[order(res_nhbe$padj),]

res_a549_bs <- results(prdds, contrast = c("group", "A549_SARS_CoV_2", "A549_Mock"), alpha = 0.05, lfcThreshold = 1)
res_a549 <- lfcShrink(prdds, contrast = c("group", "A549_SARS_CoV_2", "A549_Mock"), res = res_a549_bs, type = "ashr")
res_a549 <- na.omit(res_a549)
res_a549_ordered <- res_a549[order(res_a549$padj),]

res_a549_ace2_bs <- results(prdds, contrast = c("group", "A549_ACE2_SARS_CoV_2", "A549_ACE2_Mock"), alpha = 0.05, lfcThreshold = 1)
res_a549_ace2 <- lfcShrink(prdds, contrast = c("group", "A549_ACE2_SARS_CoV_2", "A549_ACE2_Mock"), res = res_a549_ace2_bs, type = "ashr")
res_a549_ace2 <- na.omit(res_a549_ace2)
res_a549_ace2_ordered <- res_a549_ace2[order(res_a549_ace2$padj),]

res_calu3_bs <- results(prdds, contrast = c("group", "Calu3_SARS_CoV_2", "Calu3_Mock"), alpha = 0.05, lfcThreshold = 1)
res_calu3 <- lfcShrink(prdds, contrast = c("group", "Calu3_SARS_CoV_2", "Calu3_Mock"), res = res_calu3_bs, type = "ashr")
res_calu3 <- na.omit(res_calu3)
res_calu3_ordered <- res_calu3[order(res_calu3$padj),]

#MA Plots Comparison
par(mfrow=c(2, 2))
DESeq2::plotMA(res_nhbe_ordered,
               main="Differential Expression: NHBE - Sars-CoV-2 vs Mock",
               ylim=c(-10,10), cex=0.5,
               colNonSig=adjustcolor("gray20", alpha.f = 0.5),
               colSig=adjustcolor("dodgerblue3", alpha.f = 0.5))
abline(h = c(-1, 1), col = "#ff0000", lwd = 1, lty = 2)

DESeq2::plotMA(res_a549_ordered,
               main="Differential Expression: A549 - Sars-CoV-2 vs Mock",
               ylim=c(-10,10), cex=0.5,
               colNonSig=adjustcolor("gray20", alpha.f = 0.5),
               colSig=adjustcolor("dodgerblue3", alpha.f = 0.5))
abline(h = c(-1, 1), col = "#ff0000", lwd = 1, lty = 2)

DESeq2::plotMA(res_a549_ace2_ordered,
               main="Differential Expression: A549_ACE2 - Sars-CoV-2 vs Mock",
               ylim=c(-10,10), cex=0.5,
               colNonSig=adjustcolor("gray20", alpha.f = 0.5),
               colSig=adjustcolor("dodgerblue3", alpha.f = 0.5))
abline(h = c(-1, 1), col = "#ff0000", lwd = 1, lty = 2)

DESeq2::plotMA(res_calu3_ordered,
               main="Differential Expression: Calu3 - Sars-CoV-2 vs Mock",
               ylim=c(-10,10), cex=0.5,
               colNonSig=adjustcolor("gray20", alpha.f = 0.5),
               colSig=adjustcolor("dodgerblue3", alpha.f = 0.5))
abline(h = c(-1, 1), col = "#ff0000", lwd = 1, lty = 2)

#Assigning Gene Status for Volcano Plots
res_nhbe_ordered$status <- ifelse(res_nhbe_ordered$padj < 0.05 & abs(res_nhbe_ordered$log2FoldChange) > 1,
                                  ifelse(res_nhbe_ordered$log2FoldChange > 0, "Upregulated","Downregulated"),
                                  "Non-Significant")

res_a549_ordered$status <- ifelse(res_a549_ordered$padj < 0.05 & abs(res_a549_ordered$log2FoldChange) > 1,
                                  ifelse(res_a549_ordered$log2FoldChange > 0, "Upregulated", "Downregulated"),
                                  "Non-Significant")

res_a549_ace2_ordered$status <- ifelse(res_a549_ace2_ordered$padj < 0.05 & abs(res_a549_ace2_ordered$log2FoldChange) > 1,
                                       ifelse(res_a549_ace2_ordered$log2FoldChange > 0, "Upregulated", "Downregulated"),
                                       "Non-Significant")

res_calu3_ordered$status <- ifelse(res_calu3_ordered$padj < 0.05 & abs(res_calu3_ordered$log2FoldChange) > 1,
                                   ifelse(res_calu3_ordered$log2FoldChange > 0, "Upregulated", "Downregulated"),
                                   "Non-Significant")

#Volcano Plots
volcano_colors <- c("Upregulated" = "#4DAF4A", "Downregulated" = "#E41A1C", "Non-Significant" = "#377EB8")

top_labels_nhbe <- head(res_nhbe_ordered, 15)
top_labels_a549 <- head(res_a549_ordered, 15)
top_labels_a549_ace2 <- head(res_a549_ace2_ordered, 15)
top_labels_calu3 <- head(res_calu3_ordered, 15)

ggplot(as.data.frame(res_nhbe_ordered),
  aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_manual(values = volcano_colors) +
  theme_minimal() +
  ggtitle("Volcano Plot of Differentially Expressed Genes: NHBE - SARS-CoV-2 vs Mock") +
  xlab("log2 Fold Change") +
  ylab("-log10 (Adjusted p-value)") +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_text_repel(data = as.data.frame(top_labels_nhbe),
                  aes(label = rownames(top_labels_nhbe)),
                  size = 2.5, max.overlaps = 20)

ggplot(as.data.frame(res_a549_ordered),
       aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_manual(values = volcano_colors) +
  theme_minimal() +
  ggtitle("Volcano Plot of Differentially Expressed Genes: A549 - SARS-CoV-2 vs Mock") +
  xlab("log2 Fold Change") +
  ylab("-log10 (Adjusted p-value)") +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_text_repel(data = as.data.frame(top_labels_a549),
                  aes(label = rownames(top_labels_a549)),
                  size = 2.5, max.overlaps = 20)

ggplot(as.data.frame(res_a549_ace2_ordered),
       aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_manual(values = volcano_colors) +
  theme_minimal() +
  ggtitle("Volcano Plot of Differentially Expressed Genes: A549-ACE2 - SARS-CoV-2 vs Mock") +
  xlab("log2 Fold Change") +
  ylab("-log10 (Adjusted p-value)") +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_text_repel(data = as.data.frame(top_labels_a549_ace2),
                  aes(label = rownames(top_labels_a549_ace2)),
                  size = 2.5, max.overlaps = 20)

ggplot(as.data.frame(res_calu3_ordered),
       aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_manual(values = volcano_colors) +
  theme_minimal() +
  ggtitle("Volcano Plot of Differentially Expressed Genes: Calu-3 - SARS-CoV-2 vs Mock") +
  xlab("log2 Fold Change") +
  ylab("-log10 (Adjusted p-value)") +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_text_repel(data = as.data.frame(top_labels_calu3),
                  aes(label = rownames(top_labels_calu3)),
                  size = 2.5, max.overlaps = 20)

#Extracting Significant Genes
samples_to_keep_nhbe <- rownames(META[META$group %in% c("NHBE_Mock", "NHBE_SARS_CoV_2"), ])
vsd_nhbe_subset <- vsd[, samples_to_keep_nhbe]

samples_to_keep_a549 <- rownames(META[META$group %in% c("A549_Mock", "A549_SARS_CoV_2"), ])
vsd_a549_subset <- vsd[, samples_to_keep_a549]

samples_to_keep_a549_ace2 <- rownames(META[META$group %in% c("A549_ACE2_Mock", "A549_ACE2_SARS_CoV_2"), ])
vsd_a549_ace2_subset <- vsd[, samples_to_keep_a549_ace2]

samples_to_keep_calu3 <- rownames(META[META$group %in% c("Calu3_Mock", "Calu3_SARS_CoV_2"), ])
vsd_calu3_subset <- vsd[, samples_to_keep_calu3]

sig_genes_nhbe <- as.data.frame(res_nhbe_ordered[res_nhbe_ordered$padj < 0.05 & abs(res_nhbe_ordered$log2FoldChange) > 1, ])
up_genes_nhbe <- subset(sig_genes_nhbe, log2FoldChange > 0)
down_genes_nhbe <- subset(sig_genes_nhbe, log2FoldChange < 0)

sig_genes_a549 <- as.data.frame(res_a549_ordered[res_a549_ordered$padj < 0.05 & abs(res_a549_ordered$log2FoldChange) > 1, ])
up_genes_a549 <- subset(sig_genes_a549, log2FoldChange > 0)
down_genes_a549 <- subset(sig_genes_a549, log2FoldChange < 0)

sig_genes_a549_ace2 <- as.data.frame(res_a549_ace2_ordered[res_a549_ace2_ordered$padj < 0.05 & abs(res_a549_ace2_ordered$log2FoldChange) > 1, ])
up_genes_a549_ace2 <- subset(sig_genes_a549_ace2, log2FoldChange > 0)
down_genes_a549_ace2 <- subset(sig_genes_a549_ace2, log2FoldChange < 0)

sig_genes_calu3 <- as.data.frame(res_calu3_ordered[res_calu3_ordered$padj < 0.05 & abs(res_calu3_ordered$log2FoldChange) > 1, ])
up_genes_calu3 <- subset(sig_genes_calu3, log2FoldChange > 0)
down_genes_calu3 <- subset(sig_genes_calu3, log2FoldChange < 0)

#Creating Heatmaps of the Top 10 Up-regulated and Down-Regulated Genes
heatmap_colors_nhbe <- colorRampPalette(rev(brewer.pal(name = "RdYlBu", n = 11)))(255)

top_up_nhbe <- head(up_genes_nhbe, 10)
if (nrow(top_up_nhbe) > 1) {
top_up_nhbe_exp <- assay(vsd_nhbe_subset)[rownames(top_up_nhbe), ]
pheatmap(top_up_nhbe_exp,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         show_colnames = TRUE,
         col=heatmap_colors_nhbe,
         main = "Top 10 Up-Regulated NHBE Genes: SARS-CoV-2 vs Mock")}

top_down_nhbe <- head(down_genes_nhbe, 10)
if (nrow(top_down_nhbe) > 1) {
top_down_nhbe_exp <- assay(vsd_nhbe_subset)[rownames(top_down_nhbe), ]
pheatmap(top_down_nhbe_exp,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         show_colnames = TRUE,
         col=heatmap_colors_nhbe,
         main = "Top 10 Down-Regulated NHBE Genes: SARS-CoV-2 vs Mock")}

heatmap_colors_a549 <- colorRampPalette(rev(brewer.pal(name = "PiYG", n = 11)))(255)

top_up_a549 <- head(up_genes_a549, 10)
if (nrow(top_up_a549) > 1) {
  top_up_a549_exp <- assay(vsd_a549_subset)[rownames(top_up_a549), ]
  pheatmap(top_up_a549_exp,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           scale = "row",
           show_colnames = TRUE,
           col=heatmap_colors_a549,
           main = "Top 10 Up-Regulated A549 Genes: SARS-CoV-2 vs Mock")}

top_down_a549 <- head(down_genes_a549, 10)
if (nrow(top_down_a549) > 1) {
  top_down_a549_exp <- assay(vsd_a549_subset)[rownames(top_down_a549), ]
  pheatmap(top_down_a549_exp,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           scale = "row",
           show_colnames = TRUE,
           col=heatmap_colors_a549,
           main = "Top 10 Down-Regulated A549 Genes: SARS-CoV-2 vs Mock")}

heatmap_colors_a549_ace2 <- colorRampPalette(rev(brewer.pal(name = "BrBG", n = 11)))(255)

top_up_a549_ace2 <- head(up_genes_a549_ace2, 10)
if (nrow(top_up_a549_ace2) > 1) {
  top_up_a549_ace2_exp <- assay(vsd_a549_ace2_subset)[rownames(top_up_a549_ace2), ]
  pheatmap(top_up_a549_ace2_exp,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           scale = "row",
           show_colnames = TRUE,
           col=heatmap_colors_a549_ace2,
           main = "Top 10 Up-Regulated A549-ACE2 Genes: SARS-CoV-2 vs Mock")}

top_down_a549_ace2 <- head(down_genes_a549_ace2, 10)
if (nrow(top_down_a549_ace2) > 1) {
  top_down_a549_ace2_exp <- assay(vsd_a549_ace2_subset)[rownames(top_down_a549_ace2), ]
  pheatmap(top_down_a549_ace2_exp,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           scale = "row",
           show_colnames = TRUE,
           col=heatmap_colors_a549_ace2,
           main = "Top 10 Down-Regulated A549-ACE2 Genes: SARS-CoV-2 vs Mock")}

heatmap_colors_calu3 <- colorRampPalette(rev(brewer.pal(name = "PRGn", n = 11)))(255)

top_up_calu3 <- head(up_genes_calu3, 10)
if (nrow(top_up_calu3) > 1) {
  top_up_calu3_exp <- assay(vsd_calu3_subset)[rownames(top_up_calu3), ]
  pheatmap(top_up_calu3_exp,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           scale = "row",
           show_colnames = TRUE,
           col=heatmap_colors_calu3,
           main = "Top 10 Up-Regulated Calu-3 Genes: SARS-CoV-2 vs Mock")}

top_down_calu3 <- head(down_genes_calu3, 10)
if (nrow(top_down_calu3) > 1) {
  top_down_calu3_exp <- assay(vsd_calu3_subset)[rownames(top_down_calu3), ]
  pheatmap(top_down_calu3_exp,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           scale = "row",
           show_colnames = TRUE,
           col=heatmap_colors_calu3,
           main = "Top 10 Down-Regulated Calu-3 Genes: SARS-CoV-2 vs Mock")}

#Creating DEG Analysis Files
if (!file.exists("result")) {dir.create("result")}

sig_genes_nhbe$status <- ifelse(sig_genes_nhbe$log2FoldChange > 0, "Upregulated", "Downregulated")
sig_genes_nhbe$gene_id <- rownames(sig_genes_nhbe)
sig_genes_nhbe <- sig_genes_nhbe[, c("gene_id","baseMean","log2FoldChange","lfcSE","pvalue","padj","status")]
write.table(sig_genes_nhbe, file = "result/significant_DE_NHBE_genes.csv", sep = ",", row.names = FALSE, quote = FALSE)

sig_genes_a549$status <- ifelse(sig_genes_a549$log2FoldChange > 0, "Upregulated", "Downregulated")
sig_genes_a549$gene_id <- rownames(sig_genes_a549)
sig_genes_a549 <- sig_genes_a549[, c("gene_id","baseMean","log2FoldChange","lfcSE","pvalue","padj","status")]
write.table(sig_genes_a549, file = "result/significant_DE_A549_genes.csv", sep = ",", row.names = FALSE, quote = FALSE)

sig_genes_a549_ace2$status <- ifelse(sig_genes_a549_ace2$log2FoldChange > 0, "Upregulated", "Downregulated")
sig_genes_a549_ace2$gene_id <- rownames(sig_genes_a549_ace2)
sig_genes_a549_ace2 <- sig_genes_a549_ace2[, c("gene_id","baseMean","log2FoldChange","lfcSE","pvalue","padj","status")]
write.table(sig_genes_a549_ace2, file = "result/significant_DE_A549_ACE2_genes.csv", sep = ",", row.names = FALSE, quote = FALSE)

sig_genes_calu3$status <- ifelse(sig_genes_calu3$log2FoldChange > 0, "Upregulated", "Downregulated")
sig_genes_calu3$gene_id <- rownames(sig_genes_calu3)
sig_genes_calu3 <- sig_genes_calu3[, c("gene_id","baseMean","log2FoldChange","lfcSE","pvalue","padj","status")]
write.table(sig_genes_calu3, file = "result/significant_DE_Calu3_genes.csv", sep = ",", row.names = FALSE, quote = FALSE)































