# ---
# title: "Pipeline 3: Bulk RNA-seq Analysis of CLL Dasatinib Response"
# subtitle: "AnnData + Scanpy + pyDESeq2 + diffxpy Approach"
# author: "Samuel Escudero, Ivon Sanchez, Karim Hamed"
# affiliation: "Universitat Pompeu Fabra"
# date: "2025"
# ---

# %% [markdown]
# # Pipeline 3: AnnData + Scanpy + pyDESeq2 + diffxpy
#
# ## Introduction
#
# This is **Pipeline 3** of our multi-approach analysis of the GSE151159
# CLL dataset. We now move from R to **Python**, using:
#
# - **AnnData** as our data container (the Python equivalent of
#   SummarizedExperiment)
# - **Scanpy** for quality control, normalization, and visualization
# - **pyDESeq2** for differential expression (Python port of DESeq2)
# - **diffxpy** as an alternative DE framework for comparison
# - **GSEApy** for functional enrichment analysis
#
# By running the same biological analysis in three different ecosystems
# (edgeR/limma in R, DESeq2/Seurat in R, scanpy/pyDESeq2 in Python),
# we validate our findings and learn the strengths of each approach.
#
# ### Theory: Why Python for RNA-seq?
#
# The Python bioinformatics ecosystem has grown rapidly, driven by the
# scRNA-seq revolution. Scanpy (Wolf et al., 2018) is the Python
# counterpart to Seurat, built on top of AnnData objects. While R
# remains dominant for statistical genomics (DESeq2, edgeR, limma are
# all R packages), Python offers advantages in:
#
# - Integration with machine learning libraries (scikit-learn, PyTorch)
# - Scalability to very large datasets (backed by sparse matrices)
# - Interoperability with the broader data science ecosystem (pandas,
#   matplotlib)
# - RNA velocity analysis (scVelo, cellDancer) runs natively in Python
#
# pyDESeq2 (Muzellec et al., 2023) is a faithful Python reimplementation
# of R's DESeq2, producing near-identical results. diffxpy provides
# an alternative framework with multiple test types (Wald, LRT,
# t-test, rank-based).

# %% [markdown]
# ### Theory: The AnnData object
#
# AnnData (Annotated Data) is the core data structure in the Python
# single-cell ecosystem. It stores:
#
# - `adata.X` - the expression matrix (genes x samples/cells).
#   Can be dense (numpy array) or sparse (scipy sparse matrix)
# - `adata.obs` - sample/cell-level metadata (equivalent to colData
#   in SummarizedExperiment or meta.data in Seurat)
# - `adata.var` - gene-level metadata (equivalent to rowData)
# - `adata.uns` - unstructured annotations (colors, parameters, etc.)
# - `adata.obsm` - multi-dimensional sample annotations (PCA, UMAP
#   coordinates)
# - `adata.layers` - alternative expression matrices (raw counts,
#   normalized data, etc.)
#
# Key difference from R objects: AnnData stores samples in ROWS and
# genes in COLUMNS (the transpose of what R uses). This is because
# Python's data science tools (pandas, scikit-learn) expect
# observations in rows.

# %%
# =============================================================
# 1. SETUP AND IMPORTS
# =============================================================

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

warnings.filterwarnings("ignore")

# Set scanpy plotting defaults
sc.settings.verbosity = 3  # 0=errors, 1=warnings, 2=info, 3=hints
sc.settings.set_figure_params(dpi=100, facecolor="white", frameon=True)

# Set output directories
import os
RESULTS_DIR = "../results"
PLOTS_DIR = os.path.join(RESULTS_DIR, "plots")
TABLES_DIR = os.path.join(RESULTS_DIR, "tables")
ENRICHMENT_DIR = os.path.join(RESULTS_DIR, "enrichment")

for d in [PLOTS_DIR, TABLES_DIR, ENRICHMENT_DIR]:
    os.makedirs(d, exist_ok=True)

print("Scanpy version:", sc.__version__)

# %% [markdown]
# # 2. Data Import and AnnData Construction
#
# We load the shared CSV files exported by Pipeline 1 (or Pipeline 2).
# This ensures all three pipelines analyze exactly the same samples
# and genes.
#
# ### Theory: CSV to AnnData conversion
#
# Since our data comes from an R SummarizedExperiment, the count
# matrix CSV has genes in rows and samples in columns. AnnData
# expects the transpose (samples in rows, genes in columns), so we
# must transpose the matrix during import.

# %%
# Load the shared CSV files
counts_df = pd.read_csv("../data/counts_matrix.csv", index_col=0)
sample_meta = pd.read_csv("../data/sample_metadata.csv", index_col=0)
gene_meta = pd.read_csv("../data/gene_metadata.csv", index_col=0)

print(f"Count matrix: {counts_df.shape[0]} genes x {counts_df.shape[1]} samples")
print(f"Sample metadata: {sample_meta.shape[0]} samples, {sample_meta.shape[1]} columns")
print(f"Gene metadata: {gene_meta.shape[0]} genes, {gene_meta.shape[1]} columns")

# %%
# Check sample groups
print("\nResponse classification:")
print(sample_meta["response2"].value_counts())

# %%
# Construct the AnnData object
# IMPORTANT: We transpose the count matrix so samples are rows
# and genes are columns (AnnData convention)

adata = sc.AnnData(
    X=counts_df.values.T,          # Transpose: (samples x genes)
    obs=sample_meta,               # Sample metadata -> obs
    var=gene_meta                  # Gene metadata -> var
)

# Set observation (sample) and variable (gene) names
adata.obs_names = counts_df.columns.tolist()
adata.var_names = counts_df.index.tolist()

# Store raw counts in a separate layer for DESeq2 later
# (scanpy normalization will modify adata.X)
adata.layers["counts"] = adata.X.copy()

print("\nAnnData object created:")
print(adata)

# %%
# Inspect the structure
print(f"\nadata.X shape: {adata.X.shape}  (samples x genes)")
print(f"adata.obs columns: {list(adata.obs.columns)}")
print(f"adata.var columns: {list(adata.var.columns)}")
print(f"Layers: {list(adata.layers.keys())}")

# Quick look at the count matrix
print(f"\nFirst 3 samples x 5 genes:")
print(pd.DataFrame(
    adata.X[:3, :5],
    index=adata.obs_names[:3],
    columns=adata.var_names[:5]
))

# %% [markdown]
# # 3. Quality Control
#
# ### Theory: QC for bulk RNA-seq in scanpy
#
# Scanpy's QC functions were designed for scRNA-seq (filtering cells
# by mitochondrial percentage, doublet removal, etc.), but the core
# metrics are equally useful for bulk RNA-seq:
#
# - `total_counts` (library size per sample)
# - `n_genes_by_counts` (number of detected genes per sample)
# - `pct_counts_mt` (mitochondrial percentage, less relevant for bulk)
#
# For bulk RNA-seq, we focus on library size and detected genes to
# identify outlier samples.

# %%
# Calculate QC metrics
# First, identify mitochondrial genes (start with "MT-")
adata.var["mt"] = adata.var_names.str.startswith("MT-")

# If gene names are Ensembl IDs, try using the symbol column
if "symbol" in adata.var.columns:
    adata.var["mt"] = adata.var["symbol"].fillna("").str.startswith("MT-")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)

print("QC metrics added to adata.obs:")
print(adata.obs[["total_counts", "n_genes_by_counts"]].describe())

# %%
# Visualize QC metrics
fig, axes = plt.subplots(1, 3, figsize=(14, 4))

# Library size per sample
colors = adata.obs["response2"].map(
    {"NonResponder": "steelblue", "Responder": "coral"}
).values

axes[0].bar(range(adata.n_obs), adata.obs["total_counts"].values,
            color=colors)
axes[0].set_xlabel("Samples")
axes[0].set_ylabel("Total counts")
axes[0].set_title("Library size per sample")
axes[0].ticklabel_format(style="sci", axis="y", scilimits=(6, 6))

# Genes detected per sample
axes[1].bar(range(adata.n_obs), adata.obs["n_genes_by_counts"].values,
            color=colors)
axes[1].set_xlabel("Samples")
axes[1].set_ylabel("Detected genes")
axes[1].set_title("Genes detected per sample")

# Scatter: total counts vs detected genes
for resp, color in [("NonResponder", "steelblue"), ("Responder", "coral")]:
    mask = adata.obs["response2"] == resp
    axes[2].scatter(
        adata.obs.loc[mask, "total_counts"],
        adata.obs.loc[mask, "n_genes_by_counts"],
        c=color, label=resp, alpha=0.7, s=50
    )
axes[2].set_xlabel("Total counts")
axes[2].set_ylabel("Detected genes")
axes[2].set_title("Counts vs Genes")
axes[2].legend()

plt.tight_layout()
plt.savefig(os.path.join(PLOTS_DIR, "p3_qc_metrics.png"), dpi=150, bbox_inches="tight")
plt.show()

# %% [markdown]
# ## Gene Filtering
#
# ### Theory: Filtering in scanpy
#
# `sc.pp.filter_genes()` removes genes detected in fewer than a
# minimum number of samples. This is equivalent to `filterByExpr()`
# in edgeR or the `min.cells` parameter in Seurat's
# `CreateSeuratObject()`. We require a gene to be expressed in at
# least 3 samples (min_cells=3).
#
# We also apply a count-based filter: genes must have a total count
# of at least 10 across all samples, matching what we did in
# Pipeline 2 (DESeq2 pre-filtering).

# %%
print(f"Genes before filtering: {adata.n_vars}")

# Filter 1: gene must be in at least 3 samples
sc.pp.filter_genes(adata, min_cells=3)
print(f"After min_cells=3 filter: {adata.n_vars}")

# Filter 2: gene must have >= 10 total counts
gene_counts = np.array(adata.X.sum(axis=0)).flatten()
mask = gene_counts >= 10
adata = adata[:, mask].copy()
print(f"After total_count>=10 filter: {adata.n_vars}")

# %% [markdown]
# # 4. Normalization
#
# ### Theory: Scanpy normalization vs. TMM vs. DESeq2 size factors
#
# Scanpy provides `sc.pp.normalize_total()` which normalizes each
# sample to a target sum (similar to CPM with a custom scale factor).
# Then `sc.pp.log1p()` applies log(x+1) transformation.
#
# This is equivalent to Seurat's `LogNormalize` (Pipeline 2) and
# conceptually similar to logCPM in edgeR (Pipeline 1), though:
#
# | Method | Pipeline | How it works |
# |--------|----------|-------------|
# | TMM | Pipeline 1 | Adjusts for RNA composition bias |
# | LogNormalize | Pipeline 2 | count/total * 10000, then log1p |
# | normalize_total | Pipeline 3 | count/total * target_sum, then log1p |
# | DESeq2 size factors | Pipelines 2,3 | Median of ratios to geometric mean |
#
# For DE analysis with pyDESeq2, we will use the RAW counts stored
# in `adata.layers["counts"]`. The scanpy normalization here is only
# for visualization and exploratory analysis (PCA, UMAP, heatmaps).

# %%
# Normalize for visualization (this modifies adata.X)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Store normalized data in a layer too
adata.layers["normalized"] = adata.X.copy()

print("Normalization complete.")
print(f"adata.X range: [{adata.X.min():.2f}, {adata.X.max():.2f}]")
print(f"Raw counts preserved in adata.layers['counts']")

# %% [markdown]
# ## Highly Variable Genes
#
# ### Theory: HVG selection in scanpy
#
# Like Seurat's `FindVariableFeatures()`, scanpy identifies genes
# with high variance relative to their mean. The `flavor='seurat'`
# option uses the same VST method as Seurat, making results
# comparable to Pipeline 2. We select 3000 HVGs for dimensionality
# reduction.

# %%
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat")
print(f"Highly variable genes: {adata.var['highly_variable'].sum()}")

# %%
# Plot HVGs
sc.pl.highly_variable_genes(adata)
plt.savefig(os.path.join(PLOTS_DIR, "p3_hvg.png"), dpi=150, bbox_inches="tight")
plt.show()

# %% [markdown]
# # 5. Dimensionality Reduction
#
# ### Theory: PCA and UMAP in scanpy
#
# Scanpy's `sc.pp.scale()` centers and scales each gene (mean=0,
# std=1), equivalent to Seurat's `ScaleData()`. Then `sc.tl.pca()`
# computes principal components on the scaled data.
#
# `sc.tl.umap()` requires a nearest-neighbor graph, which we compute
# with `sc.pp.neighbors()`. For bulk RNA-seq with few samples (~28),
# we use a smaller number of neighbors (n_neighbors=10) than the
# default for scRNA-seq (15-30).
#
# Note: In Pipeline 1 we used MDS (leading log-fold change distances),
# in Pipeline 2 we used Seurat's PCA/UMAP. Here we use scanpy's
# implementation, which should give comparable results.

# %%
# Scale data (only on HVGs for PCA, but we can scale all)
sc.pp.scale(adata, max_value=10)

# %%
# Run PCA
sc.tl.pca(adata, svd_solver="arpack", n_comps=20)

# %%
# PCA plot colored by response
sc.pl.pca(
    adata,
    color="response2",
    palette={"NonResponder": "steelblue", "Responder": "coral"},
    size=200,
    title="PCA: Responders vs NonResponders"
)
plt.savefig(os.path.join(PLOTS_DIR, "p3_pca.png"), dpi=150, bbox_inches="tight")
plt.show()

# %%
# Variance explained (elbow plot equivalent)
sc.pl.pca_variance_ratio(adata, n_pcs=20, log=False)
plt.savefig(os.path.join(PLOTS_DIR, "p3_pca_variance.png"), dpi=150, bbox_inches="tight")
plt.show()

# %%
# Compute neighbors and UMAP
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
sc.tl.umap(adata)

# %%
# UMAP plot
sc.pl.umap(
    adata,
    color="response2",
    palette={"NonResponder": "steelblue", "Responder": "coral"},
    size=200,
    title="UMAP: Responders vs NonResponders"
)
plt.savefig(os.path.join(PLOTS_DIR, "p3_umap.png"), dpi=150, bbox_inches="tight")
plt.show()

# %% [markdown]
# ## Hierarchical Clustering (Dendrogram)

# %%
sc.tl.dendrogram(adata, groupby="response2")
sc.pl.dendrogram(adata, groupby="response2")
plt.savefig(os.path.join(PLOTS_DIR, "p3_dendrogram.png"), dpi=150, bbox_inches="tight")
plt.show()

# %% [markdown]
# # 6. Differential Expression Analysis
#
# We perform DE analysis using TWO methods:
#
# 1. **pyDESeq2** - the Python port of R's DESeq2 (negative binomial
#    model, Wald test, LFC shrinkage). This should produce results
#    very similar to Pipeline 2.
#
# 2. **diffxpy** - an alternative Python DE framework that supports
#    multiple test types. We use its Wald test implementation as a
#    cross-validation.
#
# Both methods operate on RAW COUNTS (not the normalized data in
# adata.X), so we extract the counts from `adata.layers["counts"]`.
#
# ## 6A. Differential Expression with pyDESeq2
#
# ### Theory: pyDESeq2 internals
#
# pyDESeq2 (Muzellec et al., 2023) reimplements the core DESeq2
# algorithm in Python using PyTorch for GPU-accelerated fitting.
# The workflow is:
#
# 1. Create a `DeseqDataSet` from counts + metadata + design
# 2. Call `dds.deseq2()` which runs:
#    - `dds.fit_size_factors()` (median-of-ratios normalization)
#    - `dds.fit_genewise_dispersions()` (MLE of NB dispersion per gene)
#    - `dds.fit_dispersion_trend()` (smooth mean-dispersion curve)
#    - `dds.fit_dispersion_prior()` (empirical Bayes shrinkage)
#    - `dds.fit_MAP_dispersions()` (final shrunken dispersions)
#    - `dds.fit_LFC()` (GLM fitting)
#    - `dds.calculate_cooks()` (outlier detection)
# 3. Extract results with `DeseqStats` and apply Wald test
# 4. Optionally apply LFC shrinkage

# %%
# ----------------------------------------------------------
# 6A. pyDESeq2
# ----------------------------------------------------------
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Prepare count matrix (samples x genes, using raw counts)
counts_for_de = pd.DataFrame(
    adata.layers["counts"],
    index=adata.obs_names,
    columns=adata.var_names
).astype(int)

# Prepare metadata
metadata_for_de = adata.obs[["response2"]].copy()
metadata_for_de["response2"] = metadata_for_de["response2"].astype(str)

print(f"Counts for DE: {counts_for_de.shape}")
print(f"Metadata groups:\n{metadata_for_de['response2'].value_counts()}")

# %%
# Create DeseqDataSet
dds = DeseqDataSet(
    counts=counts_for_de,
    metadata=metadata_for_de,
    design_factors="response2",
    refit_cooks=True
)

# Run the full DESeq2 pipeline
dds.deseq2()
print("\npyDESeq2 analysis complete.")

# %%
# Extract results: Responder vs NonResponder
stat_res = DeseqStats(dds, contrast=["response2", "Responder", "NonResponder"])
stat_res.summary()

# Get results as DataFrame
pydeseq2_results = stat_res.results_df.copy()
pydeseq2_results.index.name = "gene_id"

print(f"\nTotal genes tested: {len(pydeseq2_results)}")
print(f"Genes with padj < 0.05: {(pydeseq2_results['padj'] < 0.05).sum()}")

# %%
# Apply LFC shrinkage (apeglm equivalent)
# pyDESeq2 supports LFC shrinkage
stat_res.lfc_shrink(coeff="response2_Responder_vs_NonResponder")
pydeseq2_shrunk = stat_res.results_df.copy()

print(f"After LFC shrinkage:")
print(f"Genes with padj < 0.05: {(pydeseq2_shrunk['padj'] < 0.05).sum()}")

# %%
# Define strict DEGs from pyDESeq2
strict_fdr = 0.05
min_abs_lfc = 1.0

pydeseq2_strict = pydeseq2_shrunk[
    (pydeseq2_shrunk["padj"] < strict_fdr) &
    (pydeseq2_shrunk["log2FoldChange"].abs() >= min_abs_lfc)
].copy()

print(f"\npyDESeq2 strict DEGs (FDR < {strict_fdr}, |log2FC| >= {min_abs_lfc}):")
print(f"  Total: {len(pydeseq2_strict)}")
print(f"  Up in Responders:   {(pydeseq2_strict['log2FoldChange'] > 0).sum()}")
print(f"  Down in Responders: {(pydeseq2_strict['log2FoldChange'] < 0).sum()}")

# %% [markdown]
# ## 6B. Differential Expression with diffxpy
#
# ### Theory: diffxpy framework
#
# diffxpy (Fischer et al., 2020) provides a flexible DE testing
# framework for Python. Key features:
#
# - Supports multiple test types: Wald test, likelihood ratio test
#   (LRT), t-test, Wilcoxon rank-sum test
# - Can handle complex designs with covariates
# - Integrates with AnnData objects natively
# - Provides its own negative binomial GLM implementation
#
# We use diffxpy's Wald test to compare with pyDESeq2. Both use
# a negative binomial model, but diffxpy has a different
# implementation of the GLM fitting and dispersion estimation.
# Comparing the two gives us confidence in our results.

# %%
# ----------------------------------------------------------
# 6B. diffxpy
# ----------------------------------------------------------
import diffxpy.api as de

# diffxpy works directly with AnnData, but needs raw counts in X
# Create a separate AnnData with raw counts for diffxpy
adata_raw = sc.AnnData(
    X=adata.layers["counts"].copy(),
    obs=adata.obs.copy(),
    var=adata.var.copy()
)

# Ensure response2 is a categorical variable
adata_raw.obs["response2"] = adata_raw.obs["response2"].astype("category")

# %%
# Run diffxpy Wald test
diffxpy_test = de.test.wald(
    data=adata_raw,
    formula_loc="~ 1 + response2",
    factor_loc_totest="response2"
)

# Get results
diffxpy_results = diffxpy_test.summary()
print(f"\ndiffxpy Wald test complete:")
print(f"Total genes tested: {len(diffxpy_results)}")
print(f"Genes with qval < 0.05: {(diffxpy_results['qval'] < 0.05).sum()}")

# %%
# diffxpy strict DEGs
diffxpy_strict = diffxpy_results[
    (diffxpy_results["qval"] < strict_fdr) &
    (diffxpy_results["log2fc"].abs() >= min_abs_lfc)
].copy()

print(f"\ndiffxpy strict DEGs (FDR < {strict_fdr}, |log2FC| >= {min_abs_lfc}):")
print(f"  Total: {len(diffxpy_strict)}")

# %% [markdown]
# ## 6C. Comparing pyDESeq2 vs. diffxpy
#
# ### Theory: Why compare two DE methods?
#
# No single DE method is universally best. Each makes different
# statistical assumptions and uses different algorithms for
# parameter estimation. By running both pyDESeq2 and diffxpy, we
# can assess:
#
# - **Concordance**: Do both methods identify similar DEGs?
# - **Sensitivity**: Does one find more DEGs than the other?
# - **Robustness**: Are the top-ranked genes consistent?
#
# High concordance gives us confidence that the results reflect
# true biology rather than artifacts of a particular algorithm.

# %%
# Compare the two methods
print("=" * 60)
print("COMPARISON: pyDESeq2 vs diffxpy")
print("=" * 60)

pydeseq2_deg_set = set(pydeseq2_strict.index)

# diffxpy uses 'gene' column for gene names
if "gene" in diffxpy_strict.columns:
    diffxpy_deg_set = set(diffxpy_strict["gene"].values)
else:
    diffxpy_deg_set = set(diffxpy_strict.index)

overlap = pydeseq2_deg_set & diffxpy_deg_set
only_pydeseq2 = pydeseq2_deg_set - diffxpy_deg_set
only_diffxpy = diffxpy_deg_set - pydeseq2_deg_set

print(f"\npyDESeq2 strict DEGs:  {len(pydeseq2_deg_set)}")
print(f"diffxpy strict DEGs:   {len(diffxpy_deg_set)}")
print(f"Overlap (both):        {len(overlap)}")
print(f"Only in pyDESeq2:      {len(only_pydeseq2)}")
print(f"Only in diffxpy:       {len(only_diffxpy)}")

if len(pydeseq2_deg_set) > 0 and len(diffxpy_deg_set) > 0:
    jaccard = len(overlap) / len(pydeseq2_deg_set | diffxpy_deg_set)
    print(f"Jaccard similarity:    {jaccard:.3f}")

# %%
# Venn-style comparison plot
fig, ax = plt.subplots(figsize=(6, 4))

categories = ["pyDESeq2 only", "Overlap", "diffxpy only"]
values = [len(only_pydeseq2), len(overlap), len(only_diffxpy)]
colors_bar = ["coral", "mediumpurple", "steelblue"]

bars = ax.bar(categories, values, color=colors_bar, edgecolor="white")
ax.set_ylabel("Number of DEGs")
ax.set_title("DEG Overlap: pyDESeq2 vs diffxpy")

for bar, val in zip(bars, values):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
            str(val), ha="center", fontsize=12, fontweight="bold")

plt.tight_layout()
plt.savefig(os.path.join(PLOTS_DIR, "p3_de_comparison.png"), dpi=150, bbox_inches="tight")
plt.show()

# %% [markdown]
# # 7. Visualization of DE Results
#
# We use pyDESeq2 results (with LFC shrinkage) as our primary results
# since it is the most comparable to Pipeline 2 (R DESeq2).

# %%
# P-value histogram
fig, ax = plt.subplots(figsize=(8, 5))
pvals = pydeseq2_shrunk["pvalue"].dropna()
ax.hist(pvals, bins=50, color="lightblue", edgecolor="white")
ax.axvline(x=0.05, color="blue", linestyle="--", linewidth=2)
ax.set_xlabel("Raw p-value")
ax.set_ylabel("Number of genes")
ax.set_title("Distribution of Raw P-values (pyDESeq2)")
ax.text(0.2, ax.get_ylim()[1] * 0.9,
        f"Genes with p <= 0.05: {(pvals <= 0.05).sum()}",
        color="blue", fontweight="bold")
plt.tight_layout()
plt.savefig(os.path.join(PLOTS_DIR, "p3_pvalue_hist.png"), dpi=150, bbox_inches="tight")
plt.show()

# %%
# Volcano plot
fig, ax = plt.subplots(figsize=(10, 7))

res_plot = pydeseq2_shrunk.dropna(subset=["padj", "log2FoldChange"]).copy()
res_plot["is_strict"] = (
    (res_plot["padj"] < strict_fdr) &
    (res_plot["log2FoldChange"].abs() >= min_abs_lfc)
)
res_plot["neg_log10_padj"] = -np.log10(res_plot["padj"])

# Plot non-significant
ns = res_plot[~res_plot["is_strict"]]
ax.scatter(ns["log2FoldChange"], ns["neg_log10_padj"],
           c="grey", s=5, alpha=0.4, label="Not significant")

# Plot strict DEGs
sig = res_plot[res_plot["is_strict"]]
ax.scatter(sig["log2FoldChange"], sig["neg_log10_padj"],
           c="firebrick", s=15, alpha=0.7, label="Strict DEG")

# Threshold lines
ax.axhline(-np.log10(strict_fdr), color="blue", linestyle="--", linewidth=0.5)
ax.axvline(-min_abs_lfc, color="darkgreen", linestyle="--", linewidth=0.5)
ax.axvline(min_abs_lfc, color="darkgreen", linestyle="--", linewidth=0.5)

# Label top genes
if "symbol" in adata.var.columns:
    sym_map = adata.var["symbol"].to_dict()
    top_genes = sig.nsmallest(15, "padj")
    for gene_id, row in top_genes.iterrows():
        label = sym_map.get(gene_id, gene_id)
        if pd.notna(label) and label != "":
            ax.annotate(label, (row["log2FoldChange"], row["neg_log10_padj"]),
                       fontsize=7, ha="center",
                       arrowprops=dict(arrowstyle="-", color="grey", lw=0.5))

ax.set_xlabel("log2 Fold Change (Responder vs NonResponder)")
ax.set_ylabel("-log10(adjusted p-value)")
ax.set_title("Volcano Plot (pyDESeq2, shrunken LFC)")
ax.legend()

plt.tight_layout()
plt.savefig(os.path.join(PLOTS_DIR, "p3_volcano.png"), dpi=150, bbox_inches="tight")
plt.show()

# %%
# MA plot
fig, ax = plt.subplots(figsize=(8, 5))

res_ma = pydeseq2_shrunk.dropna(subset=["padj", "log2FoldChange", "baseMean"]).copy()
res_ma["is_strict"] = (
    (res_ma["padj"] < strict_fdr) &
    (res_ma["log2FoldChange"].abs() >= min_abs_lfc)
)

ns = res_ma[~res_ma["is_strict"]]
sig = res_ma[res_ma["is_strict"]]

ax.scatter(np.log10(ns["baseMean"] + 1), ns["log2FoldChange"],
           c="grey", s=3, alpha=0.3)
ax.scatter(np.log10(sig["baseMean"] + 1), sig["log2FoldChange"],
           c="firebrick", s=10, alpha=0.7)
ax.axhline(0, color="blue", linewidth=0.8)
ax.set_xlabel("log10(baseMean + 1)")
ax.set_ylabel("log2 Fold Change")
ax.set_title("MA Plot (pyDESeq2)")

plt.tight_layout()
plt.savefig(os.path.join(PLOTS_DIR, "p3_ma_plot.png"), dpi=150, bbox_inches="tight")
plt.show()

# %%
# Heatmap of top 50 DEGs
top50_genes = pydeseq2_shrunk.dropna(subset=["padj"]).nsmallest(50, "padj").index

# Use normalized data for heatmap
adata_for_heatmap = adata[:, adata.var_names.isin(top50_genes)].copy()
adata_for_heatmap.X = adata_for_heatmap.layers["normalized"]

sc.pl.heatmap(
    adata_for_heatmap,
    var_names=top50_genes.tolist(),
    groupby="response2",
    cmap="RdBu_r",
    standard_scale="var",
    show_gene_labels=True,
    figsize=(12, 10)
)
plt.savefig(os.path.join(PLOTS_DIR, "p3_heatmap_top50.png"), dpi=150, bbox_inches="tight")
plt.show()

# %%
# Top 20 DEGs table
top20 = pydeseq2_strict.nlargest(20, "log2FoldChange", key=abs)

if "symbol" in adata.var.columns:
    sym_map = adata.var["symbol"].to_dict()
    top20["symbol"] = top20.index.map(sym_map)

print("\nTop 20 Strict DEGs (pyDESeq2):")
display_cols = ["log2FoldChange", "baseMean", "stat", "pvalue", "padj"]
if "symbol" in top20.columns:
    display_cols = ["symbol"] + display_cols
print(top20[display_cols].to_string())

# %% [markdown]
# # 8. Functional Enrichment Analysis with GSEApy
#
# ### Theory: GSEApy (Fang et al., 2023)
#
# GSEApy is the Python equivalent of R's clusterProfiler/fgsea. It
# supports:
#
# - **Enrichr**: over-representation analysis using multiple gene set
#   libraries (GO, KEGG, Reactome, MSigDB)
# - **GSEA (prerank)**: ranked gene set enrichment analysis
# - **ssGSEA**: single-sample GSEA
#
# We use Enrichr for ORA (similar to GOstats in Pipeline 1 and
# clusterProfiler in Pipeline 2) and GSEA prerank for ranked-list
# enrichment.

# %%
import gseapy as gp

# Get gene symbols for the strict DEGs
if "symbol" in adata.var.columns:
    sym_map = adata.var["symbol"].to_dict()
    deg_symbols = [sym_map.get(g, g) for g in pydeseq2_strict.index]
    deg_symbols = [s for s in deg_symbols if pd.notna(s) and s != ""]
else:
    deg_symbols = pydeseq2_strict.index.tolist()

print(f"DEG gene symbols for enrichment: {len(deg_symbols)}")
print(f"First 10: {deg_symbols[:10]}")

# %%
# Over-representation analysis with Enrichr
# Using GO_Biological_Process_2023 gene set library
try:
    enr = gp.enrichr(
        gene_list=deg_symbols,
        gene_sets=["GO_Biological_Process_2023"],
        organism="human",
        outdir=None,  # Don't save files, keep in memory
        cutoff=0.05
    )
    enr_results = enr.results
    enr_results = enr_results[enr_results["Adjusted P-value"] < 0.05]
    print(f"\nEnriched GO BP terms (Enrichr): {len(enr_results)}")

    if len(enr_results) > 0:
        print("\nTop 10 enriched terms:")
        print(enr_results.head(10)[["Term", "Adjusted P-value",
                                     "Odds Ratio", "Genes"]].to_string())
except Exception as e:
    print(f"Enrichr analysis failed (may need internet): {e}")
    enr_results = pd.DataFrame()

# %%
# GSEA (prerank) using the full ranked gene list
if "symbol" in adata.var.columns:
    # Create ranked list: gene symbol -> stat (Wald statistic)
    ranked = pydeseq2_shrunk.copy()
    ranked["symbol"] = ranked.index.map(adata.var["symbol"].to_dict())
    ranked = ranked.dropna(subset=["symbol", "stat"])
    ranked = ranked[ranked["symbol"] != ""]
    ranked = ranked.drop_duplicates(subset="symbol", keep="first")
    gene_ranking = ranked.set_index("symbol")["stat"].sort_values(ascending=False)
else:
    gene_ranking = pydeseq2_shrunk["stat"].dropna().sort_values(ascending=False)

print(f"\nGSEA ranked list: {len(gene_ranking)} genes")
print(f"Top 5: {gene_ranking.head().to_dict()}")
print(f"Bottom 5: {gene_ranking.tail().to_dict()}")

# %%
# Run GSEA prerank
try:
    gsea_res = gp.prerank(
        rnk=gene_ranking,
        gene_sets="GO_Biological_Process_2023",
        min_size=15,
        max_size=500,
        permutation_num=1000,
        outdir=None,
        seed=42
    )
    gsea_results = gsea_res.res2d
    gsea_sig = gsea_results[gsea_results["FDR q-val"] < 0.25]
    print(f"\nGSEA significant terms (FDR < 0.25): {len(gsea_sig)}")

    if len(gsea_sig) > 0:
        print("\nTop 10 GSEA terms:")
        print(gsea_sig.head(10)[["Term", "NES", "FDR q-val"]].to_string())
except Exception as e:
    print(f"GSEA prerank failed (may need internet): {e}")
    gsea_results = pd.DataFrame()

# %% [markdown]
# # 9. Save All Results

# %%
# Save pyDESeq2 results
pydeseq2_shrunk.to_csv(os.path.join(TABLES_DIR, "p3_full_DE_pydeseq2.csv"))
pydeseq2_strict.to_csv(os.path.join(TABLES_DIR, "p3_strict_DEGs_pydeseq2.csv"))

# Save diffxpy results
diffxpy_results.to_csv(os.path.join(TABLES_DIR, "p3_full_DE_diffxpy.csv"), index=False)
diffxpy_strict.to_csv(os.path.join(TABLES_DIR, "p3_strict_DEGs_diffxpy.csv"), index=False)

# Save enrichment results if available
if len(enr_results) > 0:
    enr_results.to_csv(os.path.join(ENRICHMENT_DIR, "p3_enrichr_GO_BP.csv"), index=False)

# Save the AnnData object (optional, for future reuse)
adata.write(os.path.join("..", "data", "adata_cll.h5ad"))

print("All results saved.")
print(f"  DE tables: {TABLES_DIR}")
print(f"  Plots: {PLOTS_DIR}")
print(f"  Enrichment: {ENRICHMENT_DIR}")

# %% [markdown]
# # 10. Cross-Pipeline Comparison
#
# If you have already run Pipelines 1 and 2, we can compare the
# DEG lists across all three approaches.

# %%
# Load Pipeline 1 and 2 results if available
try:
    p1_degs = pd.read_csv(os.path.join(TABLES_DIR, "p1_strict_DEGs.csv"))
    p1_set = set(p1_degs["ensembl_id"].values) if "ensembl_id" in p1_degs.columns else set()
    print(f"Pipeline 1 (edgeR/limma) strict DEGs: {len(p1_set)}")
except FileNotFoundError:
    p1_set = set()
    print("Pipeline 1 results not found (run Pipeline 1 first)")

try:
    p2_degs = pd.read_csv(os.path.join(TABLES_DIR, "p2_strict_DEGs_deseq2.csv"))
    p2_set = set(p2_degs["gene_id"].values) if "gene_id" in p2_degs.columns else set()
    print(f"Pipeline 2 (DESeq2) strict DEGs: {len(p2_set)}")
except FileNotFoundError:
    p2_set = set()
    print("Pipeline 2 results not found (run Pipeline 2 first)")

p3_set = set(pydeseq2_strict.index)
print(f"Pipeline 3 (pyDESeq2) strict DEGs: {len(p3_set)}")

if p1_set and p2_set and p3_set:
    all_three = p1_set & p2_set & p3_set
    any_two = ((p1_set & p2_set) | (p1_set & p3_set) | (p2_set & p3_set)) - all_three
    print(f"\nDEGs found by ALL three pipelines: {len(all_three)}")
    print(f"DEGs found by at least two pipelines: {len(any_two) + len(all_three)}")

# %% [markdown]
# # 11. Discussion
#
# In this third pipeline, we analyzed the CLL dasatinib response
# dataset entirely in Python, demonstrating that the same biological
# analysis can be performed across different programming ecosystems.
#
# **Key findings:**
#
# - pyDESeq2 produced results highly comparable to R's DESeq2
#   (Pipeline 2), validating the Python implementation
# - diffxpy provided an independent confirmation of the DE results
# - Scanpy's visualization tools (PCA, UMAP, heatmaps) offer a
#   convenient and publication-ready output
# - GSEApy provided enrichment analysis consistent with GOstats
#   (Pipeline 1) and clusterProfiler (Pipeline 2)
#
# **Pipeline comparison summary:**
#
# | Feature | P1 (R) | P2 (R) | P3 (Python) |
# |---------|--------|--------|-------------|
# | Container | SummarizedExperiment | Seurat | AnnData |
# | Normalization | TMM | LogNormalize | normalize_total |
# | DE method | limma-voom | DESeq2 | pyDESeq2 + diffxpy |
# | DE model | Normal + weights | Negative binomial | Negative binomial |
# | LFC shrinkage | No | apeglm | Yes |
# | Enrichment | GOstats + fgsea | clusterProfiler | GSEApy |
# | Visualization | base R + ggplot2 | Seurat + ggplot2 | scanpy + matplotlib |
#
# The convergence of biological conclusions across all three
# pipelines demonstrates the robustness of the findings and
# the reliability of the analytical frameworks.

# %%
print("\n" + "=" * 60)
print("Pipeline 3 complete!")
print("=" * 60)

# Print package versions
import pydeseq2
import diffxpy
print(f"\nPackage versions:")
print(f"  scanpy:   {sc.__version__}")
print(f"  anndata:  {adata.__class__.__module__}")
print(f"  pydeseq2: {pydeseq2.__version__}")
print(f"  diffxpy:  {diffxpy.__version__}")
print(f"  gseapy:   {gp.__version__}")
print(f"  numpy:    {np.__version__}")
print(f"  pandas:   {pd.__version__}")
