install.packages(c(
  "Seurat", "SeuratDisk", "SoupX", "DoubletFinder",
  "ggplot2", "patchwork", "cowplot", "ComplexHeatmap",
  "ggrepel", "circlize", "janitor", "here"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2", "edgeR", "limma", "nebula", "MAST",
  "fgsea", "msigdbr", "slingshot", "monocle3", "CellChat"
))

install.packages("glmGamPoi")


