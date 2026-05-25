# IEO RNA-seq CLL Multi-Pipeline Analysis

## Overview

Multi-approach transcriptomic analysis of chronic lymphocytic leukemia (CLL) 
response to dasatinib treatment, using data from 
[GSE151159](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151159) 
(Blatte et al., 2021).

We implement **three parallel pipelines** on the same dataset to compare 
data containers, normalization strategies, and statistical frameworks 
for differential expression analysis:

| Pipeline | Language | Data container | Normalization | DE method | Enrichment |
|----------|----------|----------------|---------------|-----------|------------|
| **1** | R | `SummarizedExperiment` → `DGEList` | TMM + voom | limma (eBayes) | GOstats / fgsea |
| **2** | R | `Seurat` object | LogNormalize / SCTransform | DESeq2 (Wald) | clusterProfiler |
| **3** | Python | `AnnData` (scanpy) | scanpy normalize | pyDESeq2 + diffxpy | GSEApy |

## Repository structure

```
CLL_RNAseq_MultiplePipeline/
├── data/                          # Raw data (git-ignored, see setup below)
│   ├── .gitignore
│   ├── GSE151159.rds              # SummarizedExperiment with counts + metadata
│   ├── counts_matrix.csv          # Exported by 00_export_shared_data.R
│   ├── sample_metadata.csv        # Exported by 00_export_shared_data.R
│   └── gene_metadata.csv          # Exported by 00_export_shared_data.R
├── scripts/
│   ├── 01_edgeR_limma.Rmd         # Pipeline 1: SE + edgeR/limma-voom
│   ├── 02_seurat_deseq2.Rmd       # Pipeline 2: Seurat + DESeq2
│   └── 03_scanpy_pydeseq2.Rmd     # Pipeline 3: AnnData + scanpy + pyDESeq2/diffxpy
├── results/
│   ├── plots/                     # Figures (PNG/PDF)
│   ├── tables/                    # DE results, gene lists (CSV)
│   ├── enrichment/                # GO/GSEA results (HTML/CSV)
│   └── reports/                   # Knitted HTML reports
├── 00_export_shared_data.R        # Run ONCE to extract CSVs from RDS
├── IEO_RNAseq_CLL.Rproj           # R project file
├── environment.yml                # Conda env for Pipeline 3
├── renv.lock                      
├── .gitignore
└── README.md
```

## Setup instructions

### 1. Clone the repository
```bash
git clone https://github.com/Sam-E18/CLL_RNAseq_MultiplePipeline.git
cd IEO_RNAseq_CLL_MultiPipeline
```

### 2. Get the data
Place `GSE151159.rds` in the `data/` directory. This file is available and can be downloaded from GEO accession GSE151159.

### 3. Set up the R environment (Pipelines 1 & 2)
```r
# In R, from the project root:
install.packages("renv")
renv::init()

# Core packages needed:
BiocManager::install(c(
  "SummarizedExperiment", "edgeR", "limma", "DESeq2",
  "sva", "GOstats", "org.Hs.eg.db", "fgsea", "clusterProfiler",
  "BiocStyle", "AnnotationDbi", "GenomicFeatures"
))
install.packages(c(
  "Seurat", "tidyverse", "pheatmap", "ggrepel",
  "knitr", "kableExtra", "rmarkdown", "here"
))

renv::snapshot()
```

### 4. Set up the Python environment (Pipeline 3)
```bash
conda env create -f environment.yml
conda activate ieo_rnaseq
```

### 5. Export shared data
```r
# In R, from the project root:
source("00_export_shared_data.R")
```
This creates `counts_matrix.csv`, `sample_metadata.csv`, and 
`gene_metadata.csv` in `data/` for use by all three pipelines.

### 6. Run the pipelines
```r
# Pipeline 1 (R):
rmarkdown::render("scripts/01_edgeR_limma.Rmd", output_dir = "results/reports/")

# Pipeline 2 (R):
rmarkdown::render("scripts/02_seurat_deseq2.Rmd", output_dir = "results/reports/")
```

```bash
# Pipeline 3 (Python, with conda active):
# Can be run as Jupyter notebook or as Rmd with reticulate
conda activate ieo_rnaseq
jupyter lab scripts/03_scanpy_pydeseq2.ipynb
```

## Dataset

- **GEO accession**: [GSE151159](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151159)
- **Organism**: *Homo sapiens* (CLL primary cells)
- **Design**: Responders vs. Non-responders to dasatinib (in vitro)
- **Samples**: 28 biological replicates (16 non-responders, 12 responders)
- **Sequencing**: Bulk RNA-seq

## References

- Blatte et al. (2021). Gene expression profiling predicts sensitivity of CLL cells to dasatinib.
- Love MI, Huber W, Anders S (2014). DESeq2. *Genome Biology* 15:550.
- Ritchie ME et al. (2015). limma powers DE analyses. *Nucleic Acids Research* 43(7):e47.
- Robinson MD, McCarthy DJ, Smyth GK (2010). edgeR. *Bioinformatics* 26(1):139-140.
- Stuart T et al. (2019). Comprehensive integration of single-cell data. *Cell* 177(7):1888-1902.
- Wolf FA, Angerer P, Theis FJ (2018). SCANPY. *Genome Biology* 19(1):15.

## Course resources

This project was developed as part of the IEO Transcriptomics course at UPF with some personal changes:
- [DESeq2 Vignette](https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [edgeR User Guide](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)
- [Seurat Tutorials](https://satijalab.org/seurat/articles/get_started.html)
- [Scanpy Tutorials](https://scanpy.readthedocs.io/)

## Authors

Samuel Escudero, Ivon Sanchez, Karim Hamed  
Universitat Pompeu Fabra
