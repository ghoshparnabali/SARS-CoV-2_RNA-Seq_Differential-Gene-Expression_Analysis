# SARS-CoV-2-RNA-Seq-Differential-Gene-Expression-DEG-Analysis
This repository contains an end-to-end DEG Analysis pipeline built in RStudio to analyze host transcriptional responses to SARS-CoV-2 across multiple cell lines. Utilizing raw read counts from the GSE147507 benchmark study, this pipeline handles complex multi-condition experimental matrices and executes rigorous quality control protocols.
## Key Technical Implementations & Results:
- **Differential Expression:** Executed via DESeq2 with variance-stabilized transformations applied for accurate PCA clustering and sample-to-sample distance mapping.
- **Statistical Rigor:** Implemented Adaptive Shrinkage (ashr) estimators to minimize and appropriately manage statistical noise inherent in low-count genes.
- **Data Visualization:** Engineered dynamically scaled Volcano plots and hierarchical clustering heatmaps (ggplot2, pheatmap) to accurately visualize biological outliers.
## Summary of Results:
- **Clear Infection Patterns:** The PCA models demonstrated that SARS-CoV-2 drastically alters cell behavior, cleanly grouping and separating infected cells from healthy (mock) cells.
- **Strong Immune Reaction:** The heatmaps and significant gene datasets revealed a massive spike in genes that trigger inflammation and immune defenses (such as ISG15 and CXCL1).
- **Identifying Key Markers & Cell-Specific Nuance:** The volcano plots successfully captured extreme genetic reactions to the virus. For example, the CACNA1C gene showed massive activation in NHBE and A549-ACE2 cell lines, but was sharply downregulated in standard A549 cells, underscoring the complex, targeted nature of the host response.
