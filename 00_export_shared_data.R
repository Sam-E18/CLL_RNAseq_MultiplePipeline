# ============================================
# 00_export_shared_data.R
# ============================================
# PURPOSE: Extract count matrix and metadata from the GSE151159.rds
# (SummarizedExperiment) and save as CSV files that all three pipelines
# can consume. Run this ONCE before starting any pipeline.
#
# INPUT:  data/GSE151159.rds
# OUTPUT: data/counts_matrix.csv  (genes x samples, raw integer counts)
#         data/sample_metadata.csv (sample annotations with response2)
#         data/gene_metadata.csv   (gene annotations: ensembl_id, biotype, symbol)
# ============================================

library(SummarizedExperiment)

# --- Load the SummarizedExperiment object ---
# If you installed the IEOproject R package, you can also use:
#   cll <- readRDS(file.path(system.file("extdata", package="IEOproject"), "GSE151159.rds"))
# Otherwise, load directly from data/:
cat("Loading GSE151159.rds...\n")
cll <- readRDS("data/GSE151159.rds")
cat("Loaded:", class(cll), "with", nrow(cll), "genes and", ncol(cll), "samples\n\n")

# --- Create the response2 classification ---
desc2 <- colData(cll)$description.2
response_raw <- ifelse(
  grepl("^responders_final", desc2), "Responder",
  ifelse(grepl("^nonresponders_final", desc2), "NonResponder", NA)
)
colData(cll)$response2 <- factor(response_raw, levels = c("NonResponder", "Responder"))
cat("Response classification:\n")
print(table(colData(cll)$response2, useNA = "ifany"))

# --- Export raw count matrix ---
counts_mat <- as.data.frame(assays(cll)$counts)
write.csv(counts_mat, "data/counts_matrix.csv", row.names = TRUE)
cat("\nExported counts_matrix.csv:", nrow(counts_mat), "genes x", ncol(counts_mat), "samples\n")

# --- Export sample metadata ---
sample_meta <- as.data.frame(colData(cll))
write.csv(sample_meta, "data/sample_metadata.csv", row.names = TRUE)
cat("Exported sample_metadata.csv:", nrow(sample_meta), "samples\n")

# --- Export gene metadata ---
gene_meta <- as.data.frame(rowData(cll))
write.csv(gene_meta, "data/gene_metadata.csv", row.names = TRUE)
cat("Exported gene_metadata.csv:", nrow(gene_meta), "genes\n")

cat("\n✓ All shared data files exported to data/\n")
cat("  You can now run any of the three pipelines.\n")
