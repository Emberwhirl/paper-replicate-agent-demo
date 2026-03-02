# Standard Scanpy Workflow for Single-Cell Analysis

This document outlines the standard workflow for analyzing single-cell RNA-seq data using scanpy, based on the latest best practices.

## Complete Analysis Pipeline

### 1. Data Loading and Initial Setup

```python
# Core scverse libraries
import anndata as ad
import scanpy as sc

# Configure scanpy settings
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=100, facecolor='white')

# Load data (various formats)
adata = sc.read_10x_h5('path/to/data.h5')  # For 10X h5 format
# adata = sc.read_10x_mtx('path/to/data/')  # For 10X mtx directory
# adata = sc.read_h5ad('path/to/data.h5ad')  # For h5ad format

# For multiple samples, load and concatenate
samples = {
    "s1d1": "s1d1_filtered_feature_bc_matrix.h5",
    "s1d3": "s1d3_filtered_feature_bc_matrix.h5",
}
adatas = {}
for sample_id, filename in samples.items():
    sample_adata = sc.read_10x_h5(filename)
    sample_adata.var_names_make_unique()
    adatas[sample_id] = sample_adata

adata = ad.concat(adatas, label="sample")
adata.obs_names_make_unique()
```

### 2. Quality Control (QC)

Calculate QC metrics for mitochondrial, ribosomal, and hemoglobin genes:

```python
# Normalize case to ensure cross-species compatibility and robustness

## mitochondrial genes
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
## ribosomal genes
adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
## hemoglobin genes
adata.var["hb"] = adata.var_names.str.upper().str.contains(r"^HB(?!P)")

# Calculate QC metrics with log1p transformation
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

# Visualize QC metrics
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)

# Scatter plot colored by mitochondrial percentage
sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

# Apply permissive filtering - revisit thresholds after clustering
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)
```

### 3. Doublet Detection

```python
# Run Scrublet for doublet detection (use batch_key for multiple samples)
sc.pp.scrublet(adata, batch_key="sample")

# Results stored in adata.obs['doublet_score'] and adata.obs['predicted_doublet']
# Can filter doublets now or wait until after clustering to identify high-doublet clusters
```

### 4. Normalization

```python
# Save raw count data in a layer for later use
adata.layers["counts"] = adata.X.copy()

# Normalize to median total counts (log1PF normalization)
sc.pp.normalize_total(adata)

# Log-transform the data
sc.pp.log1p(adata)
```

### 5. Feature Selection

```python
# Identify highly variable genes (use batch_key for multiple samples)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")

# Visualize highly variable genes
sc.pl.highly_variable_genes(adata)
```

### 6. Dimensionality Reduction

```python
# Principal Component Analysis (PCA)
sc.tl.pca(adata)

# Visualize PCA variance ratio to assess number of PCs
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)

# Check for batch effects or QC-driven variation in PCs
sc.pl.pca(
    adata,
    color=["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
    dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
    ncols=2,
    size=2,
)

# Compute neighborhood graph
sc.pp.neighbors(adata)

# UMAP for visualization
sc.tl.umap(adata)

# Visualize by sample to check for batch effects
sc.pl.umap(adata, color="sample", size=2)
```

### 7. Clustering

```python
# Leiden clustering (using igraph flavor for better performance)
sc.tl.leiden(adata, flavor="igraph", n_iterations=2)

# Visualize clustering results
sc.pl.umap(adata, color=['leiden'])

# Try multiple resolutions to find optimal granularity
for res in [0.02, 0.5, 2.0]:
    sc.tl.leiden(adata, key_added=f"leiden_res_{res:.2f}", resolution=res, flavor="igraph")

sc.pl.umap(
    adata,
    color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
    legend_loc="on data",
)
```

### 8. Re-assess Quality Control

Visualize QC metrics on UMAP to identify clusters driven by poor quality or doublets:

```python
# Check doublet scores
sc.pl.umap(
    adata,
    color=["leiden", "predicted_doublet", "doublet_score"],
    wspace=0.5,
    size=3,
)

# Check other QC metrics
sc.pl.umap(
    adata,
    color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
    wspace=0.5,
    ncols=2,
)
```

### 9. Marker Gene Identification

```python
# Find marker genes for each cluster using Wilcoxon test
sc.tl.rank_genes_groups(adata, groupby="leiden_res_0.50", method="wilcoxon")

# Visualize top marker genes as dot plot
sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden_res_0.50", standard_scale="var", n_genes=5)

# Get marker gene dataframe for specific cluster
marker_genes_df = sc.get.rank_genes_groups_df(adata, group="0").head(5)

# Visualize specific markers on UMAP
cluster_genes = sc.get.rank_genes_groups_df(adata, group="7").head(5)["names"]
sc.pl.umap(
    adata,
    color=[*cluster_genes, "leiden_res_0.50"],
    legend_loc="on data",
    frameon=False,
    ncols=3,
)
```

### 10. Cell Type Annotation

Use hierarchical annotation with marker genes:

```python
# Define marker genes for expected cell types
marker_genes = {
    "CD14+ Mono": ["FCN1", "CD14"],
    "CD16+ Mono": ["TCF7L2", "FCGR3A", "LYN"],
    "cDC2": ["CST3", "COTL1", "LYZ", "DMXL2", "CLEC10A", "FCER1A"],
    "Erythroblast": ["MKI67", "HBA1", "HBB"],
    "Proerythroblast": ["CDK6", "SYNGR1", "HBM", "GYPA"],
    "NK": ["GNLY", "NKG7", "CD247", "FCER1G", "TYROBP", "KLRG1", "FCGR3A"],
    "ILC": ["ID2", "PLCG2", "GNLY", "SYNE1"],
    "Naive CD20+ B": ["MS4A1", "IL4R", "IGHD", "FCRL1", "IGHM"],
    "B cells": ["MS4A1", "ITGB1", "COL4A4", "PRDM1", "IRF4", "PAX5", "BCL11A", "BLK", "IGHD", "IGHM"],
    "Plasma cells": ["MZB1", "HSP90B1", "FNDC3B", "PRDM1", "IGKC", "JCHAIN"],
    "Plasmablast": ["XBP1", "PRDM1", "PAX5"],
    "CD4+ T": ["CD4", "IL7R", "TRBC2"],
    "CD8+ T": ["CD8A", "CD8B", "GZMK", "GZMA", "CCL5", "GZMB", "GZMH"],
    "T naive": ["LEF1", "CCR7", "TCF7"],
    "pDC": ["GZMB", "IL3RA", "COBLL1", "TCF4"],
}

# Visualize with dot plot at coarse resolution
sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.02", standard_scale="var")

# Annotate broad lineages first (hierarchical approach)
adata.obs["cell_type_lvl1"] = adata.obs["leiden_res_0.02"].map(
    {
        "0": "Lymphocytes",
        "1": "Monocytes",
        "2": "Erythroid",
        "3": "B Cells",
    }
)

# Then use finer resolution for detailed annotation
sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.50", standard_scale="var")

# Visualize annotated cell types
sc.pl.umap(adata, color="cell_type_lvl1", legend_loc="on data")
```

### 11. Saving Results

```python
# Save the processed AnnData object
adata.write('results/processed_data.h5ad')

# Export results to CSV
adata.obs.to_csv('results/cell_metadata.csv')
adata.var.to_csv('results/gene_metadata.csv')
```

## Additional Analysis Options

### Trajectory Inference

```python
# PAGA (Partition-based graph abstraction)
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, color=['leiden'])

# Diffusion pseudotime (DPT)
adata.uns['iroot'] = np.flatnonzero(adata.obs['leiden'] == '0')[0]
sc.tl.dpt(adata)
sc.pl.umap(adata, color=['dpt_pseudotime'])
```

### Differential Expression Between Conditions

For more robust differential expression, consider pseudo-bulking:

```python
# Simple comparison within a cell type
sc.tl.rank_genes_groups(adata, groupby='condition', groups=['treated'],
                         reference='control', method='wilcoxon')
sc.pl.rank_genes_groups(adata, groups=['treated'])

# For more conservative statistics, pseudo-bulk by sample
# bulked = sc.get.aggregate(adata, by=["sample", "cell_type"], func="sum", layer="counts")
# Then use pydeseq2 for differential expression
```

### Gene Set Scoring

```python
# Score cells for gene set expression
gene_set = ['CD3D', 'CD3E', 'CD3G']
sc.tl.score_genes(adata, gene_set, score_name='T_cell_score')
sc.pl.umap(adata, color='T_cell_score')
```

### Batch Integration

If batch effects are observed in UMAP, consider integration methods:

```python
# For batch integration, consider:
# - scanorama: https://github.com/brianhie/scanorama
# - scvi-tools: https://scvi-tools.org
```

## Key Parameters to Adjust

- **QC thresholds**: `min_genes`, `min_cells` - start permissive, revisit after clustering
- **Normalization**: Default uses median; can specify `target_sum` if needed
- **HVG parameters**: `n_top_genes`, `batch_key` for multi-sample datasets
- **PCA components**: Check variance ratio plot (can overestimate without significant downside)
- **Clustering resolution**: Higher values give more clusters (typically 0.02-2.0)
- **n_neighbors**: Affects granularity of UMAP and clustering (typically 10-30)

## Best Practices

1. Start with permissive QC filtering, revisit after clustering
2. Save raw counts in layers: `adata.layers["counts"] = adata.X.copy()`
3. Run doublet detection with batch_key for multi-sample data
4. Use batch_key in highly_variable_genes for multi-sample data
5. Use Leiden with igraph flavor for better performance
6. Try multiple clustering resolutions to find optimal granularity
7. Use hierarchical annotation (broad lineages first, then refine)
8. Visualize QC metrics on UMAP to identify problematic clusters
9. Consider pseudo-bulking for differential expression statistics
10. Save intermediate results at key steps
