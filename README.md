# scRNA_influenza_recovery

Analysis of lung immune cell dynamics during influenza infection and recovery.  
Data: GEO accession [GSE262927].  
Pipeline: R + Seurat, pseudobulk DESeq2, fgsea, trajectories (slingshot), cell-cell communication (CellChat).

## Project structure
- `data_raw/`: original GEO .h5 files
- `data_meta/`: sample metadata (parsed GSM → SRR mapping)
- `data_processed/`: saved Seurat objects
- `scripts/`: modular analysis scripts
- `notebooks/`: Quarto notebooks for exploratory and final figures
- `reports/`: polished analysis report

## Reproducibility
- R version: 4.2.3
- Dependency management: renv
- Install packages: `Rscript install_packages.R`
