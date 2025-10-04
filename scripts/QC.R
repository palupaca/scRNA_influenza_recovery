
renv::activate()
renv::deativate()

setwd("/scratch/qlp9135/scRNA_influenza_recovery")

library(hdf5r)
library(SeuratDisk)
library(Seurat)
library(devtools)
library(data.table)
#install_github("lhe17/nebula")
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)
library(DoubletFinder)
library(stringr)
library(dplyr)
library(tidyverse)

#check meta data
meta = as.data.frame(fread("./data_meta/GSE262927_CellMetaData.csv"))
head(meta)

sample.info = as.data.frame(fread("./data_meta/SampleInfo.csv"))
sample.info = na.omit(sample.info)

sample.fil = sample.info %>%
  filter(Genotype == "Ki67Cre") %>%
  mutate(Time_Value = as.numeric(str_extract(Sacrifice_Timepoint, "\\d+"))) %>%
  mutate(filename = paste0(ID, ".h5"))

#filter ssample for first round of analysis 
sample.firstround = sample.fil %>%
  filter(Time_Value <42)



#sample check
file <- "data_raw/GSE262927_RAW/GSM8181583_EEM-scRNA-113.h5"
counts <- Read10X_h5(file)

obj <- CreateSeuratObject(counts, project = "InfluenzaRecovery", min.cells = 3, min.features = 200) #initial filtering
obj



#load first round samples
objs <- lapply(1:nrow(sample.firstround), function(i) {
  m <- sample.firstround[i,]
  counts <- Read10X_h5(file.path("./data_raw/GSE262927_RAW", m$filename))
  seurat_obj <- CreateSeuratObject(counts, project = "InfluenzaRecovery",
                                   min.cells = 3, min.features = 200)
  seurat_obj$Genotype  <- m$Genotype 
  seurat_obj$Condition <- m$Condition
  seurat_obj$Tamoxifen_Timepoint <- m$Tamoxifen_Timepoint
  seurat_obj$Sacrifice_Timepoint <- m$Sacrifice_Timepoint
  seurat_obj
})

# Merge into one
merged <- merge(objs[[1]], y = objs[-1], add.cell.ids = sample.firstround$ID)

#quick QC
# % mitochondrial genes (mouse → "^mt-")

merged <- NormalizeData(merged)
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^mt-")
merged$percent.mt <- merged@meta.data$percent.mt / 100

# Add number of genes per UMI for each cell to metadata
merged$log10GenesPerUMI <- log10(merged$nFeature_RNA) / log10(merged$nCount_RNA)

colnames(merged@meta.data)
summary(merged$percent.mt)

# Create metadata dataframe
metadata <- merged@meta.data
metadata$cells <- rownames(metadata)


# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Visualize the number of cell counts per sample
metadata %>% 
  ggplot(aes(x=Sacrifice_Timepoint, fill=Sacrifice_Timepoint)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=Sacrifice_Timepoint, x=nUMI, fill= Sacrifice_Timepoint)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

#nUMI looks good (mostly >1000)
ggplot(metadata, aes(x = Sacrifice_Timepoint, y = nUMI, fill = Sacrifice_Timepoint)) +
  geom_violin(trim = FALSE) +
  scale_y_log10() +
  theme_classic()

# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=Sacrifice_Timepoint, x=nGene, fill= Sacrifice_Timepoint)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)


library(scales)
options(scipen = 999)
# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=Sacrifice_Timepoint, x=percent.mt, fill=Sacrifice_Timepoint)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

#need to match meta data orig.ident to pull annotation

# Filter out low quality cells
filtered_seurat <- subset(x = merged, 
                          subset= (nCount_RNA >= 500) & 
                            (nFeature_RNA >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (percent.mt < 0.20))

#check cell number differences
merged
filtered_seurat #no cell filter (count table may be pre-filtered)

# save merged data 
save(merged, file="data_processed/merged.RData")


filtered_seurat <- NormalizeData(filtered_seurat)
filtered_seurat <- FindVariableFeatures(filtered_seurat, selection.method = "vst", nfeatures = 2000)
filtered_seurat <- ScaleData(filtered_seurat)
filtered_seurat <- RunPCA(filtered_seurat, npcs = 30)
filtered_seurat <- RunUMAP(filtered_seurat, dims = 1:30)
filtered_seurat <- FindNeighbors(filtered_seurat, dims = 1:30)
filtered_seurat <- FindClusters(filtered_seurat, resolution = 0.5)


library(DoubletFinder)
#filter out doublet
nExp <- round(0.05 * ncol(merged))   # try 0.05 first

sweep.res <- paramSweep(merged, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Pick the pK with best BCmetric (look at bcmvn output)
pK_optimal <- bcmvn$pK[which.max(bcmvn$BCmetric)]

merged <- doubletFinder_v3(merged,
                              PCs = 1:30,
                              pN = 0.25,
                              pK = as.numeric(as.character(pK_optimal)),
                              nExp = nExp,
                              reuse.pANN = FALSE,
                              sct = FALSE)




