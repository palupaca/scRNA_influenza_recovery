Preprocessing mouse influenza recovery scRNAseq samples
================

``` r
knitr::opts_chunk$set(echo = TRUE)
library(hdf5r)
library(SeuratDisk)
```

    ## The legacy packages maptools, rgdal, and rgeos, underpinning the sp package,
    ## which was just loaded, will retire in October 2023.
    ## Please refer to R-spatial evolution reports for details, especially
    ## https://r-spatial.org/r/2023/05/15/evolution4.html.
    ## It may be desirable to make the sf package available;
    ## package maintainers should consider adding sf to Suggests:.
    ## The sp package is now running under evolution status 2
    ##      (status 2 uses the sf package in place of rgdal)

    ## Registered S3 method overwritten by 'SeuratDisk':
    ##   method            from  
    ##   as.sparse.H5Group Seurat

``` r
library(Seurat)
```

    ## Loading required package: SeuratObject

    ## Loading required package: sp

    ## 
    ## Attaching package: 'SeuratObject'

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, saveRDS

    ## Loading Seurat v5 beta version 
    ## To maintain compatibility with previous workflows, new Seurat objects will use the previous object structure by default
    ## To use new Seurat v5 assays: Please run: options(Seurat.object.assay.version = 'v5')

``` r
library(devtools)
```

    ## Loading required package: usethis

``` r
library(data.table)
#library(DoubletFinder)
library(stringr)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:data.table':
    ## 
    ##     between, first, last

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:dbplyr':
    ## 
    ##     ident, sql

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyverse)
```

    ## Warning in system("timedatectl", intern = TRUE): running command 'timedatectl'
    ## had status 1

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ forcats   1.0.0     ✔ readr     2.1.4
    ## ✔ ggplot2   4.0.0     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::between()     masks data.table::between()
    ## ✖ dplyr::filter()      masks stats::filter()
    ## ✖ dplyr::first()       masks data.table::first()
    ## ✖ purrr::flatten_df()  masks hdf5r::flatten_df()
    ## ✖ lubridate::hour()    masks data.table::hour()
    ## ✖ dplyr::ident()       masks dbplyr::ident()
    ## ✖ lubridate::isoweek() masks data.table::isoweek()
    ## ✖ dplyr::lag()         masks stats::lag()
    ## ✖ dplyr::last()        masks data.table::last()
    ## ✖ lubridate::mday()    masks data.table::mday()
    ## ✖ lubridate::minute()  masks data.table::minute()
    ## ✖ lubridate::month()   masks data.table::month()
    ## ✖ lubridate::quarter() masks data.table::quarter()
    ## ✖ lubridate::second()  masks data.table::second()
    ## ✖ dplyr::sql()         masks dbplyr::sql()
    ## ✖ purrr::transpose()   masks data.table::transpose()
    ## ✖ lubridate::wday()    masks data.table::wday()
    ## ✖ lubridate::week()    masks data.table::week()
    ## ✖ lubridate::yday()    masks data.table::yday()
    ## ✖ lubridate::year()    masks data.table::year()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
#check meta data
meta = as.data.frame(fread("/scratch/qlp9135/scRNA_influenza_recovery/data_meta/GSE262927_CellMetaData.csv"))
head(meta)
```

    ##   cellbarcode               orig.ident nCount_RNA nFeature_RNA percent.mito
    ## 1          NA scRNA-113_H1N1_02to06dpi       4677         2014    2.6085097
    ## 2          NA scRNA-113_H1N1_02to06dpi       4231         1764    1.7237308
    ## 3          NA scRNA-113_H1N1_02to06dpi       8980         2992    4.4362909
    ## 4          NA scRNA-113_H1N1_02to06dpi       4058         1300    2.9078364
    ## 5          NA scRNA-113_H1N1_02to06dpi      15304         4390    1.7294264
    ## 6          NA scRNA-113_H1N1_02to06dpi        809          464    0.3708282
    ##   phase tamoxifen_start_day sacrifice_day sex experimental_group   trace_call
    ## 1   G2M                   2             6   f   H1N1_tam02_sac06 Not_detected
    ## 2   G2M                   2             6   f   H1N1_tam02_sac06     Untraced
    ## 3    G1                   2             6   f   H1N1_tam02_sac06     Untraced
    ## 4    G1                   2             6   f   H1N1_tam02_sac06 Not_detected
    ## 5    G1                   2             6   f   H1N1_tam02_sac06     Untraced
    ## 6   G2M                   2             6   f   H1N1_tam02_sac06 Not_detected
    ##       lineage     celltype      subtype               cb
    ## 1 Endothelium         CAP1         CAP1 AAACCCAAGTCAACAA
    ## 2    Lymphoid T_lymphocyte T_lymphocyte AAACCCAGTTAATCGC
    ## 3 Endothelium         CAP1         CAP1 AAACCCATCACGGACC
    ## 4  Mesenchyme          AF1          AF1 AAACCCATCTTCCCGA
    ## 5 Endothelium         CAP1         CAP1 AAACGAAAGTAGCATA
    ## 6  Mesenchyme          AF1          AF1 AAACGAATCGAGTTGT

``` r
#add unique barcode column for merging
meta$cell = str_c(meta$orig.ident, meta$cb, sep = "_")
head
```

    ## function (x, ...) 
    ## UseMethod("head")
    ## <bytecode: 0x3252140>
    ## <environment: namespace:utils>

``` r
#check GEO sample info
sample.info = as.data.frame(fread("/scratch/qlp9135/scRNA_influenza_recovery/data_meta/SampleInfo.csv"))
sample.info = na.omit(sample.info)

sample.fil = sample.info %>%
  filter(Genotype == "Ki67Cre") %>%
  mutate(Time_Value = as.numeric(str_extract(Sacrifice_Timepoint, "\\d+"))) %>%
  mutate(filename = paste0(ID, ".h5"))

#filter ssample for first round of analysis 
sample.firstround = sample.fil %>%
  filter(Time_Value <42)

library(data.table)

sample.firstround = as.data.frame(fread("/scratch/qlp9135/scRNA_influenza_recovery/data_meta/sample_firstround.csv"))


head(sample.firstround)
```

    ##   V1         ID                       Sample_Name Genotype      Condition
    ## 1  1 GSM8181583     Ki67Cre_H1N1_tam02_sac06_rep1  Ki67Cre H1N1_Infection
    ## 2  2 GSM8181584     Ki67Cre_H1N1_tam07_sac11_rep1  Ki67Cre H1N1_Infection
    ## 3  3 GSM8181585     Ki67Cre_H1N1_tam14_sac19_rep1  Ki67Cre H1N1_Infection
    ## 4  4 GSM8181586     Ki67Cre_H1N1_tam21_sac25_rep1  Ki67Cre H1N1_Infection
    ## 5  5 GSM8181592     Ki67Cre_H1N1_tam02_sac06_rep2  Ki67Cre H1N1_Infection
    ## 6  6 GSM8181593 Ki67Cre_H1N1_homeostasis_d03_rep1  Ki67Cre    Homeostasis
    ##   Tamoxifen_Timepoint Sacrifice_Timepoint Replicate Time_Value      filename
    ## 1               2_dpi               6_dpi         1          6 GSM8181583.h5
    ## 2               7_dpi              11_dpi         1         11 GSM8181584.h5
    ## 3              14_dpi              19_dpi         1         19 GSM8181585.h5
    ## 4              21_dpi              25_dpi         1         25 GSM8181586.h5
    ## 5               2_dpi               6_dpi         2          6 GSM8181592.h5
    ## 6                  no              3_days         1          3 GSM8181593.h5
    ##                      batch
    ## 1 scRNA-113_H1N1_02to06dpi
    ## 2 scRNA-116_H1N1_07to11dpi
    ## 3 scRNA-117_H1N1_14to19dpi
    ## 4 scRNA-118_H1N1_21to25dpi
    ## 5  scRNA-162_H1N1_02to6dpi
    ## 6 scRNA-166_homeostasis_D3

``` r
#load first round samples
objs <- lapply(1:nrow(sample.firstround), function(i) {
  m <- sample.firstround[i,]
  counts <- Read10X_h5(file.path("/scratch/qlp9135/scRNA_influenza_recovery/data_raw/GSE262927_RAW", m$filename))
  seurat_obj <- CreateSeuratObject(counts, project = "InfluenzaRecovery",
                                   min.cells = 3, min.features = 200)
  seurat_obj$Genotype  <- m$Genotype 
  seurat_obj$Condition <- m$Condition
  seurat_obj$Tamoxifen_Timepoint <- m$Tamoxifen_Timepoint
  seurat_obj$Sacrifice_Timepoint <- m$Sacrifice_Timepoint
  seurat_obj
})
```

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')
    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')
    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')
    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')
    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')
    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')
    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')
    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes
    ## ('-')

``` r
# Merge into one
merged <- merge(objs[[1]], y = objs[-1], add.cell.ids = sample.firstround$batch
                )
```

``` r
#merge meta data (clsuter assignment) provided by the authors
seurat_cells <- tibble(cell = colnames(merged)) 

meta_join <- seurat_cells %>%
  left_join(meta, by = c("cell"))


keep <- meta_join %>% 
  filter(!is.na(phase)) %>%
  select(-2)  %>%
 pull(cell) %>%
  unique()

length(keep)
```

    ## [1] 43773

``` r
length(colnames(merged)) 
```

    ## [1] 55062

``` r
#Subset the Seurat object to matched cells only
merged_filtered <- subset(merged, cells = keep)

merged_filtered #43773 total cells
```

    ## An object of class Seurat 
    ## 26400 features across 43773 samples within 1 assay 
    ## Active assay: RNA (26400 features, 0 variable features)
    ##  2 layers present: counts, data

``` r
#add meta data back to seurat object
meta_fill <- meta_join %>% 
  filter(!is.na(phase)) %>%
  select(-2) %>%
  tibble::column_to_rownames("cell")

merged_filtered <- AddMetaData(merged_filtered, meta_fill)

head(merged_filtered@meta.data)
```

    ##                                                         orig.ident nCount_RNA
    ## scRNA-113_H1N1_02to06dpi_AAACCCAAGTCAACAA scRNA-113_H1N1_02to06dpi       4677
    ## scRNA-113_H1N1_02to06dpi_AAACCCAGTTAATCGC scRNA-113_H1N1_02to06dpi       4231
    ## scRNA-113_H1N1_02to06dpi_AAACCCATCACGGACC scRNA-113_H1N1_02to06dpi       8980
    ## scRNA-113_H1N1_02to06dpi_AAACCCATCTTCCCGA scRNA-113_H1N1_02to06dpi       4058
    ## scRNA-113_H1N1_02to06dpi_AAACGAAAGTAGCATA scRNA-113_H1N1_02to06dpi      15304
    ## scRNA-113_H1N1_02to06dpi_AAACGAATCGAGTTGT scRNA-113_H1N1_02to06dpi        809
    ##                                           nFeature_RNA Genotype      Condition
    ## scRNA-113_H1N1_02to06dpi_AAACCCAAGTCAACAA         2014  Ki67Cre H1N1_Infection
    ## scRNA-113_H1N1_02to06dpi_AAACCCAGTTAATCGC         1764  Ki67Cre H1N1_Infection
    ## scRNA-113_H1N1_02to06dpi_AAACCCATCACGGACC         2992  Ki67Cre H1N1_Infection
    ## scRNA-113_H1N1_02to06dpi_AAACCCATCTTCCCGA         1300  Ki67Cre H1N1_Infection
    ## scRNA-113_H1N1_02to06dpi_AAACGAAAGTAGCATA         4390  Ki67Cre H1N1_Infection
    ## scRNA-113_H1N1_02to06dpi_AAACGAATCGAGTTGT          464  Ki67Cre H1N1_Infection
    ##                                           Tamoxifen_Timepoint
    ## scRNA-113_H1N1_02to06dpi_AAACCCAAGTCAACAA               2_dpi
    ## scRNA-113_H1N1_02to06dpi_AAACCCAGTTAATCGC               2_dpi
    ## scRNA-113_H1N1_02to06dpi_AAACCCATCACGGACC               2_dpi
    ## scRNA-113_H1N1_02to06dpi_AAACCCATCTTCCCGA               2_dpi
    ## scRNA-113_H1N1_02to06dpi_AAACGAAAGTAGCATA               2_dpi
    ## scRNA-113_H1N1_02to06dpi_AAACGAATCGAGTTGT               2_dpi
    ##                                           Sacrifice_Timepoint percent.mito
    ## scRNA-113_H1N1_02to06dpi_AAACCCAAGTCAACAA               6_dpi    2.6085097
    ## scRNA-113_H1N1_02to06dpi_AAACCCAGTTAATCGC               6_dpi    1.7237308
    ## scRNA-113_H1N1_02to06dpi_AAACCCATCACGGACC               6_dpi    4.4362909
    ## scRNA-113_H1N1_02to06dpi_AAACCCATCTTCCCGA               6_dpi    2.9078364
    ## scRNA-113_H1N1_02to06dpi_AAACGAAAGTAGCATA               6_dpi    1.7294264
    ## scRNA-113_H1N1_02to06dpi_AAACGAATCGAGTTGT               6_dpi    0.3708282
    ##                                           phase tamoxifen_start_day
    ## scRNA-113_H1N1_02to06dpi_AAACCCAAGTCAACAA   G2M                   2
    ## scRNA-113_H1N1_02to06dpi_AAACCCAGTTAATCGC   G2M                   2
    ## scRNA-113_H1N1_02to06dpi_AAACCCATCACGGACC    G1                   2
    ## scRNA-113_H1N1_02to06dpi_AAACCCATCTTCCCGA    G1                   2
    ## scRNA-113_H1N1_02to06dpi_AAACGAAAGTAGCATA    G1                   2
    ## scRNA-113_H1N1_02to06dpi_AAACGAATCGAGTTGT   G2M                   2
    ##                                           sacrifice_day sex experimental_group
    ## scRNA-113_H1N1_02to06dpi_AAACCCAAGTCAACAA             6   f   H1N1_tam02_sac06
    ## scRNA-113_H1N1_02to06dpi_AAACCCAGTTAATCGC             6   f   H1N1_tam02_sac06
    ## scRNA-113_H1N1_02to06dpi_AAACCCATCACGGACC             6   f   H1N1_tam02_sac06
    ## scRNA-113_H1N1_02to06dpi_AAACCCATCTTCCCGA             6   f   H1N1_tam02_sac06
    ## scRNA-113_H1N1_02to06dpi_AAACGAAAGTAGCATA             6   f   H1N1_tam02_sac06
    ## scRNA-113_H1N1_02to06dpi_AAACGAATCGAGTTGT             6   f   H1N1_tam02_sac06
    ##                                             trace_call     lineage     celltype
    ## scRNA-113_H1N1_02to06dpi_AAACCCAAGTCAACAA Not_detected Endothelium         CAP1
    ## scRNA-113_H1N1_02to06dpi_AAACCCAGTTAATCGC     Untraced    Lymphoid T_lymphocyte
    ## scRNA-113_H1N1_02to06dpi_AAACCCATCACGGACC     Untraced Endothelium         CAP1
    ## scRNA-113_H1N1_02to06dpi_AAACCCATCTTCCCGA Not_detected  Mesenchyme          AF1
    ## scRNA-113_H1N1_02to06dpi_AAACGAAAGTAGCATA     Untraced Endothelium         CAP1
    ## scRNA-113_H1N1_02to06dpi_AAACGAATCGAGTTGT Not_detected  Mesenchyme          AF1
    ##                                                subtype               cb
    ## scRNA-113_H1N1_02to06dpi_AAACCCAAGTCAACAA         CAP1 AAACCCAAGTCAACAA
    ## scRNA-113_H1N1_02to06dpi_AAACCCAGTTAATCGC T_lymphocyte AAACCCAGTTAATCGC
    ## scRNA-113_H1N1_02to06dpi_AAACCCATCACGGACC         CAP1 AAACCCATCACGGACC
    ## scRNA-113_H1N1_02to06dpi_AAACCCATCTTCCCGA          AF1 AAACCCATCTTCCCGA
    ## scRNA-113_H1N1_02to06dpi_AAACGAAAGTAGCATA         CAP1 AAACGAAAGTAGCATA
    ## scRNA-113_H1N1_02to06dpi_AAACGAATCGAGTTGT          AF1 AAACGAATCGAGTTGT

``` r
# save filtered merged data 
save(merged_filtered, file="/scratch/qlp9135/scRNA_influenza_recovery/data_processed/merged_filtered.RData")
```
