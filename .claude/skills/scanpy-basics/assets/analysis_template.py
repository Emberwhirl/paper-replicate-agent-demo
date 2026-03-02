#!/usr/bin/env python3
"""
Complete Single-Cell Analysis Template

This template provides a complete workflow for single-cell RNA-seq analysis
using scanpy, from data loading through clustering and cell type annotation.
Based on the latest scanpy best practices.

Customize the parameters and sections as needed for your specific dataset.
"""

# Core scverse libraries
import anndata as ad
import scanpy as sc
import numpy as np
import os

# ============================================================================
# CONFIGURATION
# ============================================================================

# File paths
INPUT_FILE = 'data/raw_counts.h5ad'  # Change to your input file
OUTPUT_DIR = 'results/'
FIGURES_DIR = 'figures/'

# QC parameters (permissive - revisit after clustering)
MIN_GENES = 100          # Minimum genes per cell
MIN_CELLS = 3            # Minimum cells per gene

# Analysis parameters
N_TOP_GENES = 2000       # Number of highly variable genes
BATCH_KEY = None         # Set to column name if you have multiple samples (e.g., "sample")

# Scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor='white')
sc.settings.figdir = FIGURES_DIR

# ============================================================================
# 1. LOAD DATA
# ============================================================================

print("=" * 80)
print("LOADING DATA")
print("=" * 80)

# Load data (adjust based on your file format)
adata = sc.read_h5ad(INPUT_FILE)
# adata = sc.read_10x_h5('data/filtered_feature_bc_matrix.h5')  # For 10X h5
# adata = sc.read_10x_mtx('data/filtered_gene_bc_matrices/')   # For 10X mtx

# For multiple samples, load and concatenate:
# samples = {
#     "sample1": "sample1_filtered_feature_bc_matrix.h5",
#     "sample2": "sample2_filtered_feature_bc_matrix.h5",
# }
# adatas = {}
# for sample_id, filename in samples.items():
#     sample_adata = sc.read_10x_h5(filename)
#     sample_adata.var_names_make_unique()
#     adatas[sample_id] = sample_adata
# adata = ad.concat(adatas, label="sample")
# adata.obs_names_make_unique()
# BATCH_KEY = "sample"

print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

# ============================================================================
# 2. QUALITY CONTROL
# ============================================================================

print("\n" + "=" * 80)
print("QUALITY CONTROL")
print("=" * 80)

# Normalize case to ensure cross-species compatibility and robustness

## mitochondrial genes
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
## ribosomal genes
adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
## hemoglobin genes
adata.var["hb"] = adata.var_names.str.upper().str.contains(r"^HB(?!P)")

# Calculate QC metrics with log1p transformation
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

# Visualize QC metrics before filtering
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
    save='_qc_before_filtering'
)

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt",
              save='_qc_scatter')

# Apply permissive filtering (revisit after clustering)
print(f"\nBefore filtering: {adata.n_obs} cells, {adata.n_vars} genes")

sc.pp.filter_cells(adata, min_genes=MIN_GENES)
sc.pp.filter_genes(adata, min_cells=MIN_CELLS)

print(f"After filtering: {adata.n_obs} cells, {adata.n_vars} genes")

# ============================================================================
# 3. DOUBLET DETECTION
# ============================================================================

print("\n" + "=" * 80)
print("DOUBLET DETECTION")
print("=" * 80)

# Run Scrublet for doublet detection
if BATCH_KEY:
    sc.pp.scrublet(adata, batch_key=BATCH_KEY)
else:
    sc.pp.scrublet(adata)

n_doublets = adata.obs['predicted_doublet'].sum()
print(f"Predicted doublets: {n_doublets} ({n_doublets/adata.n_obs*100:.1f}%)")

# ============================================================================
# 4. NORMALIZATION
# ============================================================================

print("\n" + "=" * 80)
print("NORMALIZATION")
print("=" * 80)

# Save raw count data in a layer
adata.layers["counts"] = adata.X.copy()

# Normalize to median total counts (log1PF normalization)
sc.pp.normalize_total(adata)

# Log-transform
sc.pp.log1p(adata)

# ============================================================================
# 5. FEATURE SELECTION
# ============================================================================

print("\n" + "=" * 80)
print("FEATURE SELECTION")
print("=" * 80)

# Identify highly variable genes
if BATCH_KEY:
    sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES, batch_key=BATCH_KEY)
else:
    sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES)

# Visualize
sc.pl.highly_variable_genes(adata, save='_hvg')

print(f"Selected {sum(adata.var.highly_variable)} highly variable genes")

# ============================================================================
# 6. DIMENSIONALITY REDUCTION
# ============================================================================

print("\n" + "=" * 80)
print("DIMENSIONALITY REDUCTION")
print("=" * 80)

# PCA
sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, save='_pca_variance')

# Check for batch effects or QC-driven variation
if BATCH_KEY:
    sc.pl.pca(
        adata,
        color=[BATCH_KEY, BATCH_KEY, "pct_counts_mt", "pct_counts_mt"],
        dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
        ncols=2,
        size=2,
        save='_pca_qc'
    )

# Compute neighborhood graph
sc.pp.neighbors(adata)

# UMAP
sc.tl.umap(adata)

# ============================================================================
# 7. CLUSTERING
# ============================================================================

print("\n" + "=" * 80)
print("CLUSTERING")
print("=" * 80)

# Leiden clustering with multiple resolutions
for res in [0.02, 0.5, 2.0]:
    sc.tl.leiden(adata, key_added=f"leiden_res_{res:.2f}", resolution=res, flavor="igraph", n_iterations=2)

# Visualize different resolutions
sc.pl.umap(
    adata,
    color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
    legend_loc="on data",
    save='_leiden_resolutions'
)

# Use resolution 0.5 as default
adata.obs['leiden'] = adata.obs['leiden_res_0.50']

print(f"Identified {len(adata.obs['leiden'].unique())} clusters (resolution 0.5)")

# ============================================================================
# 8. RE-ASSESS QUALITY CONTROL
# ============================================================================

print("\n" + "=" * 80)
print("RE-ASSESSING QUALITY CONTROL")
print("=" * 80)

# Visualize doublet scores
sc.pl.umap(
    adata,
    color=["leiden", "predicted_doublet", "doublet_score"],
    wspace=0.5,
    size=3,
    save='_qc_doublets'
)

# Visualize other QC metrics
sc.pl.umap(
    adata,
    color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
    wspace=0.5,
    ncols=2,
    save='_qc_metrics'
)

print("Review UMAP plots to identify clusters driven by poor quality or doublets")

# ============================================================================
# 9. MARKER GENE IDENTIFICATION
# ============================================================================

print("\n" + "=" * 80)
print("MARKER GENE IDENTIFICATION")
print("=" * 80)

# Find marker genes
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')

# Visualize top markers
sc.pl.rank_genes_groups_dotplot(adata, groupby='leiden', standard_scale='var', n_genes=5,
                                 save='_markers_dotplot')

# Get top markers for each cluster
for cluster in sorted(adata.obs['leiden'].unique()):
    print(f"\nCluster {cluster} top markers:")
    markers = sc.get.rank_genes_groups_df(adata, group=cluster).head(5)
    print(markers[['names', 'scores', 'pvals_adj']].to_string(index=False))

# ============================================================================
# 10. CELL TYPE ANNOTATION (CUSTOMIZE THIS SECTION)
# ============================================================================

print("\n" + "=" * 80)
print("CELL TYPE ANNOTATION")
print("=" * 80)

# Example marker genes for common cell types (customize for your data)
marker_genes = {
    "CD14+ Mono": ["FCN1", "CD14"],
    "CD16+ Mono": ["TCF7L2", "FCGR3A", "LYN"],
    "cDC2": ["CST3", "COTL1", "LYZ", "CLEC10A", "FCER1A"],
    "NK": ["GNLY", "NKG7", "CD247", "FCER1G", "TYROBP"],
    "B cells": ["MS4A1", "CD79A", "CD79B"],
    "Plasma cells": ["MZB1", "IGKC", "JCHAIN"],
    "CD4+ T": ["CD4", "IL7R", "TRBC2"],
    "CD8+ T": ["CD8A", "CD8B", "GZMK"],
    "T naive": ["LEF1", "CCR7", "TCF7"],
}

# Visualize marker genes at coarse resolution
sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.02", standard_scale="var",
              save='_markers_coarse')

# Hierarchical annotation - first annotate broad lineages
# Customize this mapping based on your marker analysis
adata.obs["cell_type_lvl1"] = adata.obs["leiden_res_0.02"].map(
    {
        "0": "Lymphocytes",
        "1": "Myeloid",
        "2": "Unknown",  # Update based on your data
    }
)
adata.obs["cell_type_lvl1"] = adata.obs["cell_type_lvl1"].fillna("Unknown")

# Visualize at finer resolution
sc.pl.dotplot(adata, marker_genes, groupby="leiden", standard_scale="var",
              save='_markers_fine')

# Visualize annotated cell types
sc.pl.umap(adata, color="cell_type_lvl1", legend_loc="on data", save='_celltypes')

# ============================================================================
# 11. ADDITIONAL ANALYSES (OPTIONAL)
# ============================================================================

print("\n" + "=" * 80)
print("ADDITIONAL ANALYSES")
print("=" * 80)

# PAGA trajectory analysis (optional)
# sc.tl.paga(adata, groups='leiden')
# sc.pl.paga(adata, color='leiden', save='_paga')

# Gene set scoring (optional)
# example_gene_set = ['CD3D', 'CD3E', 'CD3G']
# sc.tl.score_genes(adata, example_gene_set, score_name='T_cell_score')
# sc.pl.umap(adata, color='T_cell_score', save='_gene_set_score')

# ============================================================================
# 12. SAVE RESULTS
# ============================================================================

print("\n" + "=" * 80)
print("SAVING RESULTS")
print("=" * 80)

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Save processed AnnData object
adata.write(f'{OUTPUT_DIR}/processed_data.h5ad')
print(f"Saved processed data to {OUTPUT_DIR}/processed_data.h5ad")

# Export metadata
adata.obs.to_csv(f'{OUTPUT_DIR}/cell_metadata.csv')
adata.var.to_csv(f'{OUTPUT_DIR}/gene_metadata.csv')
print(f"Saved metadata to {OUTPUT_DIR}/")

# Export marker genes
for cluster in adata.obs['leiden'].unique():
    markers = sc.get.rank_genes_groups_df(adata, group=cluster)
    markers.to_csv(f'{OUTPUT_DIR}/markers_cluster_{cluster}.csv', index=False)
print(f"Saved marker genes to {OUTPUT_DIR}/")

# ============================================================================
# 13. SUMMARY
# ============================================================================

print("\n" + "=" * 80)
print("ANALYSIS SUMMARY")
print("=" * 80)

print(f"\nFinal dataset:")
print(f"  Cells: {adata.n_obs}")
print(f"  Genes: {adata.n_vars}")
print(f"  Clusters (res 0.5): {len(adata.obs['leiden'].unique())}")
print(f"  Predicted doublets: {adata.obs['predicted_doublet'].sum()}")

print(f"\nCell type distribution (level 1):")
print(adata.obs['cell_type_lvl1'].value_counts())

print("\n" + "=" * 80)
print("ANALYSIS COMPLETE")
print("=" * 80)
print("\nNext steps:")
print("1. Review QC UMAP plots to filter problematic clusters if needed")
print("2. Refine cell type annotations based on marker gene expression")
print("3. Consider batch integration if batch effects are observed")
