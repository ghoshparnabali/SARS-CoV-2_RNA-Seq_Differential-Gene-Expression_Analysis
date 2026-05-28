# SARS-CoV-2_RNA-Seq_Differential-Gene-Expression_Analysis_Across_Multiple_Cell-lines
This repository contains an end-to-end RNA-Seq DGE Analysis pipeline built in RStudio to analyze host transcriptional responses to SARS-CoV-2 across multiple cell lines. Using raw read counts from the GSE147507 benchmark study, this pipeline handles complex multi-condition experimental matrices and implements rigorous quality-control protocols.
## Key Technical Implementations & Results:
- **Differential Expression:** Executed via DESeq2 with variance-stabilized transformations applied for accurate PCA clustering and sample-to-sample distance mapping.
- **Statistical Rigor:** Implemented Adaptive Shrinkage (ashr) estimators to minimize and appropriately manage statistical noise inherent in low-count genes.
- **Data Visualization:** Engineered dynamically scaled Volcano plots and hierarchical clustering heatmaps (ggplot2, pheatmap) to accurately visualize biological outliers.
## Summary of Results:
- **Clear Infection Patterns:**  "The PCA plot reveals distinct macro-clustering based on cell line, indicating that the virus does not override the fundamental biological signature of the host cell. But, while the infected cells cluster closely with their respective healthy (mock) controls, SARS-CoV-2 infection induces clear behavioral shifts within individual cell lines, forming micro-clusters separate from the healthy cells."
- **Strong Immune Reaction:** The heatmaps and significant gene datasets reveal a massive spike in genes that trigger inflammation and immune defenses (such as ISG15 and CXCL1).
- **Identifying Key Markers & Cell-Specific Nuance:** The volcano plots successfully capture extreme genetic reactions to the virus. For example, the CACNA1C gene showed massive activation in NHBE and A549-ACE2 cell lines, but was sharply downregulated in standard A549 cells, underscoring the complex, targeted nature of the host response.
