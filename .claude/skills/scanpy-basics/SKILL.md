---
name: scanpy-basics
description: Standard single-cell RNA-seq analysis pipeline including QC, normalization, dimensionality reduction (PCA/UMAP/t-SNE), clustering, differential expression, and visualization. This is the default skill for most single-cell analysis tasks. Use specialized skills only when specifically requested: single-cell-rna-qc for MAD-based filtering, single-cell-integration-ingest-bbknn for batch correction, or scvi-tools for deep learning methods.
---

# Scanpy: Single-Cell Analysis

## Overview

Scanpy is a scalable Python toolkit for analyzing single-cell RNA-seq data, built on AnnData. Apply this skill for complete single-cell workflows including quality control, normalization, dimensionality reduction, clustering, marker gene identification, visualization, and trajectory analysis.

## When to Use This Skill

This is the **default skill for most single-cell analysis tasks**. Use this skill when:

- Analyzing single-cell RNA-seq data (.h5ad, 10X, CSV formats)
- Working with CellRanger output or any pre-filtered scRNA-seq data
- Performing quality control on scRNA-seq datasets
- Creating UMAP, t-SNE, or PCA visualizations
- Identifying cell clusters and finding marker genes
- Annotating cell types based on gene expression
- Conducting trajectory inference or pseudotime analysis
- Generating publication-quality single-cell plots

**Note:** Use `kallisto-bustools-import` first only if you have kallisto/bustools (kb-python) output with unfiltered empty droplets. For CellRanger output or pre-filtered data, start directly with this skill.

## Quick Start

### Basic Import and Setup

```python
# Core scverse libraries
import anndata as ad
import scanpy as sc

# Configure settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor='white')
sc.settings.figdir = './figures/'
```

### Loading Data

```python
# From 10X Genomics
adata = sc.read_10x_mtx('path/to/data/')
adata = sc.read_10x_h5('path/to/data.h5')

# From h5ad (AnnData format)
adata = sc.read_h5ad('path/to/data.h5ad')

# From CSV
adata = sc.read_csv('path/to/data.csv')

# Loading and combining multiple samples
samples = {
    "sample1": "sample1_filtered_feature_bc_matrix.h5",
    "sample2": "sample2_filtered_feature_bc_matrix.h5",
}
adatas = {}
for sample_id, filename in samples.items():
    sample_adata = sc.read_10x_h5(filename)
    sample_adata.var_names_make_unique()
    adatas[sample_id] = sample_adata

adata = ad.concat(adatas, label="sample")
adata.obs_names_make_unique()
```

### Understanding AnnData Structure

The AnnData object is the core data structure in scanpy:

```python
adata.X          # Expression matrix (cells × genes)
adata.obs        # Cell metadata (DataFrame)
adata.var        # Gene metadata (DataFrame)
adata.uns        # Unstructured annotations (dict)
adata.obsm       # Multi-dimensional cell data (PCA, UMAP)
adata.layers     # Additional data layers (e.g., raw counts)

# Access cell and gene names
adata.obs_names  # Cell barcodes
adata.var_names  # Gene names
```

## Standard Analysis Workflow

### 1. Quality Control

Identify and filter low-quality cells and genes. Calculate QC metrics for mitochondrial, ribosomal, and hemoglobin genes:

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

# Filter cells and genes with permissive thresholds
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)
```

### 2. Doublet Detection

Identify potential doublets using Scrublet:

```python
# Run scrublet for doublet detection (use batch_key for multiple samples)
sc.pp.scrublet(adata, batch_key="sample")

# Results are stored in adata.obs['doublet_score'] and adata.obs['predicted_doublet']
# Can filter doublets directly or wait until clustering to filter high-doublet clusters
```

**Use the QC script for automated analysis:**

```bash
python scripts/qc_analysis.py input_file.h5ad --output filtered.h5ad
```

### 3. Normalization

Apply count depth scaling with log1p transformation (log1PF normalization):

```python
# Save raw count data in a layer
adata.layers["counts"] = adata.X.copy()

# Normalize to median total counts
sc.pp.normalize_total(adata)

# Log-transform the data
sc.pp.log1p(adata)
```

### 4. Feature Selection

Identify highly variable genes for dimensionality reduction:

```python
# Identify highly variable genes (use batch_key for multiple samples)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")

# Visualize highly variable genes
sc.pl.highly_variable_genes(adata)
```

### 5. Dimensionality Reduction

```python
# PCA
sc.tl.pca(adata)

# Inspect variance contribution of PCs
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)

# Visualize PCs colored by sample and QC metrics
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
sc.pl.umap(adata, color="sample", size=2)
```

### 6. Clustering

```python
# Leiden clustering (using igraph for better performance)
sc.tl.leiden(adata, flavor="igraph", n_iterations=2)
sc.pl.umap(adata, color=["leiden"])

# Try multiple resolutions to find optimal granularity
for res in [0.02, 0.5, 2.0]:
    sc.tl.leiden(adata, key_added=f"leiden_res_{res:.2f}", resolution=res, flavor="igraph")

sc.pl.umap(
    adata,
    color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
    legend_loc="on data",
)
```

### 7. Re-assess Quality Control

Visualize QC metrics on UMAP to identify problematic clusters:

```python
sc.pl.umap(
    adata,
    color=["leiden", "predicted_doublet", "doublet_score"],
    wspace=0.5,
    size=3,
)

sc.pl.umap(
    adata,
    color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
    wspace=0.5,
    ncols=2,
)
```

### 8. Marker Gene Identification

```python
# Find marker genes for each cluster using Wilcoxon test
sc.tl.rank_genes_groups(adata, groupby="leiden_res_0.50", method="wilcoxon")

# Visualize top markers
sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden_res_0.50", standard_scale="var", n_genes=5)

# Get results as DataFrame
sc.get.rank_genes_groups_df(adata, group="0").head(5)
```

### 9. Cell Type Annotation

Use known marker genes and hierarchical annotation:

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

# Visualize markers with dot plot
sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.02", standard_scale="var")

# Hierarchical annotation - first annotate broad lineages
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

# Visualize annotated types
sc.pl.umap(adata, color="cell_type_lvl1", legend_loc="on data")
```

### 10. Save Results

```python
# Save processed data
adata.write('results/processed_data.h5ad')

# Export metadata
adata.obs.to_csv('results/cell_metadata.csv')
adata.var.to_csv('results/gene_metadata.csv')
```

## Common Tasks

### Creating Publication-Quality Plots

```python
# Set high-quality defaults
sc.settings.set_figure_params(dpi=300, frameon=False, figsize=(5, 5))
sc.settings.file_format_figs = 'pdf'

# UMAP with custom styling
sc.pl.umap(adata, color='cell_type',
           palette='Set2',
           legend_loc='on data',
           legend_fontsize=12,
           legend_fontoutline=2,
           frameon=False,
           save='_publication.pdf')

# Save plot by accessing returned figure
sc.pl.umap(adata, color=["leiden"], show=False).figure.savefig("output_path.pdf")

# Heatmap of marker genes
sc.pl.heatmap(adata, var_names=genes, groupby='cell_type',
              swap_axes=True, show_gene_labels=True,
              save='_markers.pdf')

# Dot plot
sc.pl.dotplot(adata, var_names=genes, groupby='cell_type',
              save='_dotplot.pdf')
```

Refer to `references/plotting_guide.md` for comprehensive visualization examples.

### Trajectory Inference

```python
# PAGA (Partition-based graph abstraction)
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, color='leiden')

# Diffusion pseudotime
adata.uns['iroot'] = np.flatnonzero(adata.obs['leiden'] == '0')[0]
sc.tl.dpt(adata)
sc.pl.umap(adata, color='dpt_pseudotime')
```

### Differential Expression Between Conditions

For more conservative statistical analysis, consider pseudo-bulking:

```python
# Compare treated vs control within cell types
adata_subset = adata[adata.obs['cell_type'] == 'T cells']
sc.tl.rank_genes_groups(adata_subset, groupby='condition',
                         groups=['treated'], reference='control')
sc.pl.rank_genes_groups(adata_subset, groups=['treated'])

# For more robust DE, pseudo-bulk by sample and use pydeseq2
# sc.get.aggregate(adata, by=["sample", "cell_type"], func="sum", layer="counts")
```

### Gene Set Scoring

```python
# Score cells for gene set expression
gene_set = ['CD3D', 'CD3E', 'CD3G']
sc.tl.score_genes(adata, gene_set, score_name='T_cell_score')
sc.pl.umap(adata, color='T_cell_score')
```

### Batch Correction

If batch effects are observed in UMAP, consider integration methods:

```python
# For batch integration when specifically requested:
# - single-cell-integration-ingest-bbknn skill: for ingest and BBKNN methods
# - scvi-tools skill: for deep learning integration (scVI, scANVI) when specifically requested
```

## Key Parameters to Adjust

### Quality Control

- `min_genes`: Minimum genes per cell (typically 100-200 for permissive filtering)
- `min_cells`: Minimum cells per gene (typically 3)
- `pct_counts_mt`: Mitochondrial threshold (consider revisiting after clustering)

### Normalization

- `target_sum`: Target counts per cell (default uses median, can specify e.g., 1e4)

### Feature Selection

- `n_top_genes`: Number of HVGs (typically 2000-3000)
- `batch_key`: Use for datasets with multiple samples/batches

### Dimensionality Reduction

- `n_pcs`: Number of principal components (check variance ratio plot, can overestimate)
- `n_neighbors`: Number of neighbors (typically 10-30)

### Clustering

- `resolution`: Clustering granularity (0.02-2.0, higher = more clusters)
- `flavor`: Use "igraph" for better performance
- `n_iterations`: Number of iterations (2 is often sufficient)

## Common Pitfalls and Best Practices

1. **Save raw counts in layers**: `adata.layers["counts"] = adata.X.copy()` before normalization
2. **Run doublet detection**: Use `sc.pp.scrublet()` with `batch_key` for multi-sample data
3. **Start with permissive QC**: Filter minimally, then revisit after clustering
4. **Use batch_key for multi-sample data**: In `highly_variable_genes()` and `scrublet()`
5. **Use Leiden with igraph**: `flavor="igraph", n_iterations=2` for better performance
6. **Try multiple clustering resolutions**: Find optimal granularity for your biology
7. **Use hierarchical annotation**: Annotate broad lineages first, then refine
8. **Visualize QC on UMAP**: Check for clusters driven by doublets or poor quality
9. **Check PCA variance ratio**: Determine optimal number of PCs (can overestimate)
10. **Save intermediate results**: Long workflows can fail partway through

## Bundled Resources

### scripts/qc_analysis.py

Automated quality control script that calculates metrics (including mitochondrial, ribosomal, and hemoglobin genes), generates plots, performs doublet detection, and filters data:

```bash
python scripts/qc_analysis.py input.h5ad --output filtered.h5ad \
    --min-genes 100 --min-cells 3
```

### references/standard_workflow.md

Complete step-by-step workflow with detailed explanations and code examples for:

- Data loading and combining multiple samples
- Quality control with visualization
- Doublet detection with Scrublet
- Normalization (log1PF)
- Feature selection with batch correction
- Dimensionality reduction (PCA, UMAP)
- Clustering with multiple resolutions (Leiden)
- QC re-assessment on embeddings
- Marker gene identification
- Hierarchical cell type annotation
- Differential expression with pseudo-bulking considerations

Read this reference when performing a complete analysis from scratch.

### references/api_reference.md

Quick reference guide for scanpy functions organized by module:

- Reading/writing data (`sc.read_*`, `adata.write_*`)
- Preprocessing (`sc.pp.*`) including doublet detection
- Tools (`sc.tl.*`)
- Plotting (`sc.pl.*`)
- AnnData structure, layers, and manipulation
- Settings and utilities

Use this for quick lookup of function signatures and common parameters.

### references/plotting_guide.md

Comprehensive visualization guide including:

- Quality control plots
- Doublet score visualization
- Dimensionality reduction visualizations
- Clustering visualizations
- Marker gene plots (heatmaps, dot plots, violin plots)
- Trajectory and pseudotime plots
- Publication-quality customization
- Multi-panel figures
- Color palettes and styling

Consult this when creating publication-ready figures.

### assets/analysis_template.py

Complete analysis template providing a full workflow from data loading through cell type annotation, following latest best practices. Copy and customize this template for new analyses:

```bash
cp assets/analysis_template.py my_analysis.py
# Edit parameters and run
python my_analysis.py
```

The template includes all standard steps with configurable parameters and helpful comments.

## Additional Resources

- **Official scanpy documentation**: https://scanpy.readthedocs.io/
- **Scanpy tutorials**: https://scanpy-tutorials.readthedocs.io/
- **scverse ecosystem**: https://scverse.org/ (related tools: squidpy, scvi-tools, cellrank)
- **Single Cell Best Practices**: https://www.sc-best-practices.org/

## Tips for Effective Analysis

1. **Start with the template**: Use `assets/analysis_template.py` as a starting point
2. **Run QC script first**: Use `scripts/qc_analysis.py` for initial quality assessment
3. **Consult references as needed**: Load workflow and API references into context
4. **Use batch_key for multi-sample data**: Apply consistently across QC and HVG selection
5. **Iterate on clustering**: Try multiple resolutions and visualization methods
6. **Validate biologically**: Check marker genes match expected cell types
7. **Document parameters**: Record QC thresholds and analysis settings
8. **Save checkpoints**: Write intermediate results at key steps
