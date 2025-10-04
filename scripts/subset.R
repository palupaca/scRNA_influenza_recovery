library(hdf5r)
library(SeuratDisk)
library(Seurat)
library(devtools)
library(data.table)
#library(DoubletFinder)
library(stringr)
library(dplyr)
library(tidyverse)


#check immune population

head(meta_fill
     )

#check group name for subsetting
unique(meta_fill$lineage)
unique(meta_fill$celltype)
unique(meta_fill$subtype)

# pick immune quickly via lineage or celltype
immune_keep <- c(
  "T_lymphocyte","B_lymphocyte","NK_cell","Plasma_cell",
  "Neutrophil","iMON","cMON","pMON",
  "aMAC","iMAC",
  "cDC1","cDC2","maDC"
)

immune <- subset(merged_filtered, subset = celltype %in% immune_keep) #7239

#quick sanity check

set.seed(46)
DefaultAssay(immune) <- "RNA"

# normalize → HVGs → scale → PCA
immune <- NormalizeData(immune, verbose = FALSE)
immune <- FindVariableFeatures(immune, nfeatures = 3000, verbose = FALSE)
immune <- ScaleData(immune, verbose = FALSE)
immune <- RunPCA(immune, npcs = 50, verbose = FALSE)

ElbowPlot(immune, ndims = 50) #40 looks good

dims_use <- 1:40
immune <- FindNeighbors(immune, dims = dims_use, verbose = FALSE)
immune <- FindClusters(immune, resolution = 0.4, verbose = FALSE)
immune <- RunUMAP(immune, dims = dims_use, verbose = FALSE)

# by author celltype (primary sanity check)
DimPlot(immune, group.by = "celltype", label = TRUE, repel = TRUE) +
  ggtitle("Immune UMAP by author celltype")

# by phase to make sure time isn’t borking structure
DimPlot(immune, group.by = "phase") + ggtitle("Immune UMAP by phase")

# split by phase (check for batch drifts across time)
DimPlot(immune, group.by = "celltype", split.by = "phase", ncol = 3) +
  ggtitle("Immune UMAP split by phase")

# by time point
DimPlot(immune, group.by = "Sacrifice_Timepoint", shuffle = TRUE) +
  ggtitle("Immune UMAP by recovery timepoint")

# by sample
DimPlot(immune, group.by = "orig.ident", shuffle = TRUE) +
  ggtitle("Immune UMAP by sample")

# table to see power for pseudobulk
xtabs(~ celltype + Sacrifice_Timepoint, data = immune@meta.data)

# per-sample cell counts per type 
per_sample_counts <- immune@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  group_by(orig.ident, Sacrifice_Timepoint, celltype) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  arrange(celltype, Sacrifice_Timepoint, desc(n_cells))

per_sample_counts %>% head()

write.csv(per_sample_counts, file = "./data_processed/immunecells_per_sample_counts.csv", row.names = F)

#save immune cell seurat object

save(immune, file="/scratch/qlp9135/scRNA_influenza_recovery/data_processed/immune.RData")

#subset T cell for pseudobulk

load("./data_processed/immune.RData")

#add recovery phase labels
phase_map <- c(`6_dpi`="Acute", `11_dpi`="Acute",
               `19_dpi`="EarlyRepair", `25_dpi`="EarlyRepair",
               `42_dpi`="Remodeling", `90_dpi`="Convalescent",
               `3_days`="Control")
immune$phase2 <- dplyr::recode(immune$Sacrifice_Timepoint, !!!phase_map, .default = immune$phase)

meta = immune@meta.data

#filter for T cells
#check T cell metadata
tmeta <- immune@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  filter(celltype == "T_lymphocyte") %>%
  count(orig.ident, phase2, name = "n_cells")

min_cells_T <- 50
keep_samples_T <- tmeta %>%
  filter(n_cells >= min_cells_T) %>%
  pull(orig.ident) %>%
  unique()

T_only <- subset(immune, subset = celltype == "T_lymphocyte" & orig.ident %in% keep_samples_T) #kept all samples
#1541 cells in total

#helper function
make_pseudobulk <- function(seu, celltype_name, min_cells = 50) {
  counts <- GetAssayData(seu, slot="counts")
  md <- seu@meta.data %>%
    mutate(cell = rownames(.)) %>%
    filter(celltype == celltype_name) %>%
    mutate(group_id = paste(orig.ident, celltype, sep="||"))
  tab <- md %>% count(group_id, name="n_cells") %>% filter(n_cells >= min_cells)
  md <- md %>% filter(group_id %in% tab$group_id)
  if (nrow(md) == 0) stop("No groups left after filtering.")
  groups <- unique(md$group_id)
  G <- Matrix(0, nrow = ncol(counts), ncol = length(groups), sparse = TRUE,
              dimnames = list(colnames(counts), groups))
  G[cbind(match(md$cell, rownames(G)), match(md$group_id, colnames(G)))] <- 1L
  pb <- counts %*% G
  meta <- md %>% distinct(group_id, celltype, phase2, Sacrifice_Timepoint, Tamoxifen_Timepoint) 
  rownames(meta) = NULL 
  meta = meta %>%
    tibble::column_to_rownames("group_id") %>%
    as.data.frame()
  list(counts = pb, coldata = meta)
}

pb <- make_pseudobulk(T_only, "T_lymphocyte", min_cells = min_cells_T)  
cts <- round(pb$counts) #Make the count matrix integer
depth <- Matrix::colSums(cts)
cts <- cts[, depth >= 5e4, drop = FALSE]
coldata <- pb$coldata[colnames(cts), , drop = FALSE]
coldata$phase2 <- factor(coldata$phase2, levels = c("Control","Acute","EarlyRepair"))


#removing the control samples -> more straight forward DEseq2 design

sel <- coldata$phase2 %in% c("Acute","EarlyRepair")
cts2 <- cts[, sel, drop=FALSE]
col2 <- droplevels(coldata[sel, ])

# Two-level factor with explicit reference
col2$phase2 <- relevel(factor(col2$phase2), ref = "Acute")


table(col2$phase2) #4 samples per group


library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = cts2,
                              colData   = col2,
                              design    = ~ phase2)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds, fitType = "parametric")

#quick PCA

vsd <- vst(dds, blind=FALSE)


PCA_1 = plotPCA(vsd, intgroup=c("phase2")) #PRETTY!!
PCA_2= plotPCA(vsd, intgroup=c("Sacrifice_Timepoint")) 

#more distinct differences in acute phase (between 6_dpi and 11_dpi)


library(patchwork)
res_ER_vs_A <- results(dds,
                         contrast = c("phase2","EarlyRepair","Acute"),
                       alpha = 0.05, tidy = T)
res_ER_vs_A <- subset(res_ER_vs_A, padj<0.05) #234



summary(res_ER_vs_A)




