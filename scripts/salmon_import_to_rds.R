# ============================================================
# salmon_import_to_rds.R
# ============================================================
# PURPOSE: Import Salmon quantification results into R using tximeta,
#          summarize to gene level, and save as an RDS file.
#
# PREREQUISITES: You must have already run Salmon on the cluster
#                (Steps 0-5 in the .Rmd guide).
#
# INPUT:  salmon_quant/<sample_name>/quant.sf (one per sample)
# OUTPUT: results/gene_counts_se.rds (SummarizedExperiment)
#         results/gene_counts_matrix.csv (count matrix)
#         results/sample_metadata.csv (sample info)
#
# USAGE:  Rscript salmon_import_to_rds.R
#         or source("salmon_import_to_rds.R") in RStudio
# ============================================================

# --- 1. Load libraries ---
cat("Loading libraries...\n")

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("tximeta")) BiocManager::install("tximeta")
if (!require("SummarizedExperiment")) BiocManager::install("SummarizedExperiment")

library(tximeta)
library(SummarizedExperiment)

# --- 2. Define paths ---
# CHANGE THIS to your actual project directory
project_dir <- "~/projects/rnaseq_project"
salmon_dir  <- file.path(project_dir, "salmon_quant")
results_dir <- file.path(project_dir, "results")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# --- 3. Build sample table ---
cat("\n--- Building sample table ---\n")

# Find all quant.sf files
quant_files <- list.files(
  salmon_dir,
  pattern = "quant.sf",
  recursive = TRUE,
  full.names = TRUE
)

# Extract sample names from directory names
sample_names <- basename(dirname(quant_files))

cat("Found", length(quant_files), "quantification files:\n")
cat(paste("  ", sample_names), sep = "\n")

# Build the coldata table
# IMPORTANT: Replace the condition assignments with YOUR experimental groups
# Option A: Manual assignment (if you know the order)
# condition <- factor(c(
#   "Control", "Control", "Control", "Control", "Control", "Control", "Control",
#   "Treatment", "Treatment", "Treatment", "Treatment", "Treatment", "Treatment", "Treatment"
# ))

# Option B: Read from a metadata file (recommended)
# If you have a CSV with sample names and conditions:
# meta <- read.csv(file.path(project_dir, "sample_metadata.csv"))
# condition <- factor(meta$condition)

# Option C: Extract from sample names (if condition is encoded in the name)
# For example: "control_rep1" and "treatment_rep1"
condition <- ifelse(grepl("control|ctrl", sample_names, ignore.case = TRUE),
                    "Control", "Treatment")
condition <- factor(condition)

coldata <- data.frame(
  names     = sample_names,
  files     = quant_files,
  condition = condition,
  stringsAsFactors = FALSE
)

cat("\nSample table:\n")
print(coldata[, c("names", "condition")])

# --- 4. Import with tximeta ---
cat("\n--- Importing Salmon quantifications with tximeta ---\n")

# tximeta automatically detects the reference transcriptome
# and links it to the correct gene annotations
se <- tryCatch({
  tximeta(coldata)
}, error = function(e) {
  cat("Auto-detection failed. Trying manual linked transcriptome...\n")
  cat("Error was:", e$message, "\n\n")

  # Manual linking: provide reference info explicitly
  # Adjust these paths to match YOUR reference files
  gtf_path   <- file.path(project_dir, "reference", "Homo_sapiens.GRCh38.109.gtf")
  fasta_path <- file.path(project_dir, "reference", "Homo_sapiens.GRCh38.cdna.all.fa.gz")
  index_dir  <- file.path(project_dir, "salmon_index")

  makeLinkedTxome(
    indexDir = index_dir,
    source   = "Ensembl",
    organism = "Homo sapiens",
    release  = "109",
    genome   = "GRCh38",
    fasta    = fasta_path,
    gtf      = gtf_path,
    write    = FALSE
  )

  tximeta(coldata)
})

cat("Transcript-level import complete:\n")
cat("  Transcripts:", nrow(se), "\n")
cat("  Samples:", ncol(se), "\n")
cat("  Assays:", paste(names(assays(se)), collapse = ", "), "\n")

# --- 5. Summarize to gene level ---
cat("\n--- Summarizing to gene level ---\n")

gse <- summarizeToGene(se)

cat("Gene-level SummarizedExperiment created:\n")
cat("  Genes:", nrow(gse), "\n")
cat("  Samples:", ncol(gse), "\n")
cat("  Assays:", paste(names(assays(gse)), collapse = ", "), "\n")
cat("    counts:    raw integer counts (for DESeq2, edgeR)\n")
cat("    abundance: TPM values (for visualization)\n")
cat("    length:    effective gene lengths\n")

# --- 6. Inspect the data ---
cat("\n--- Data inspection ---\n")

counts_matrix <- assay(gse, "counts")

# Library sizes
lib_sizes <- colSums(counts_matrix)
cat("\nLibrary sizes (millions of reads):\n")
print(round(lib_sizes / 1e6, 2))

# Detected genes
detected <- colSums(counts_matrix > 0)
cat("\nDetected genes per sample:\n")
print(detected)

# Summary statistics
cat("\nGenes with > 0 counts in all samples:",
    sum(rowSums(counts_matrix > 0) == ncol(gse)), "\n")
cat("Genes with > 10 total counts:",
    sum(rowSums(counts_matrix) > 10), "\n")
cat("Genes with > 100 total counts:",
    sum(rowSums(counts_matrix) > 100), "\n")

# Gene annotations
cat("\nGene annotation columns:", colnames(rowData(gse)), "\n")
cat("\nFirst 5 gene annotations:\n")
print(head(as.data.frame(rowData(gse)), 5))

# --- 7. Save as RDS ---
cat("\n--- Saving results ---\n")

# Save the full SummarizedExperiment (recommended)
rds_path <- file.path(results_dir, "gene_counts_se.rds")
saveRDS(gse, file = rds_path)
cat("RDS saved to:", rds_path, "\n")
cat("File size:", round(file.size(rds_path) / 1e6, 1), "MB\n")

# Save the count matrix as CSV (for Python/other tools)
csv_counts_path <- file.path(results_dir, "gene_counts_matrix.csv")
write.csv(as.data.frame(counts_matrix), csv_counts_path, row.names = TRUE)
cat("Count matrix CSV saved to:", csv_counts_path, "\n")

# Save sample metadata as CSV
csv_meta_path <- file.path(results_dir, "sample_metadata.csv")
write.csv(as.data.frame(colData(gse)), csv_meta_path, row.names = TRUE)
cat("Sample metadata CSV saved to:", csv_meta_path, "\n")

# Save gene annotations as CSV
csv_gene_path <- file.path(results_dir, "gene_metadata.csv")
write.csv(as.data.frame(rowData(gse)), csv_gene_path, row.names = TRUE)
cat("Gene metadata CSV saved to:", csv_gene_path, "\n")

# --- 8. How to reload later ---
cat("\n--- How to reload the data in future sessions ---\n")
cat("
# In any R session:
library(SummarizedExperiment)
gse <- readRDS('", rds_path, "')

# Access components:
counts_matrix <- assay(gse, 'counts')     # Raw counts
tpm_matrix    <- assay(gse, 'abundance')   # TPM values
sample_info   <- colData(gse)              # Sample metadata
gene_info     <- rowData(gse)              # Gene annotations

# Feed into DESeq2:
library(DESeq2)
dds <- DESeqDataSet(gse, design = ~ condition)
dds <- DESeq(dds)

# Feed into edgeR:
library(edgeR)
dge <- DGEList(counts = counts_matrix,
               samples = as.data.frame(colData(gse)),
               genes = as.data.frame(rowData(gse)))
", sep = "")

cat("\n============================================================\n")
cat("Pipeline complete! Your data is ready for downstream analysis.\n")
cat("============================================================\n")

# --- Session info ---
cat("\nSession info:\n")
sessionInfo()
