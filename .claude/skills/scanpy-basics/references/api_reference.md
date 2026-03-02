# Scanpy API Quick Reference

Quick reference for commonly used scanpy functions organized by module.

## Import Convention

```python
import anndata as ad
import scanpy as sc
```

## Reading and Writing Data (sc.read_*)

### Reading Functions

```python
sc.read_10x_h5(filename)                    # Read 10X HDF5 file
sc.read_10x_mtx(path)                       # Read 10X mtx directory
sc.read_h5ad(filename)                      # Read h5ad (AnnData) file
sc.read_csv(filename)                       # Read CSV file
sc.read_excel(filename)                     # Read Excel file
sc.read_loom(filename)                      # Read loom file
sc.read_text(filename)                      # Read text file
sc.read_visium(path)                        # Read Visium spatial data
```

### Combining Multiple Samples

```python
# Use anndata.concat for combining multiple samples
adatas = {}
for sample_id, filename in samples.items():
    sample_adata = sc.read_10x_h5(filename)
    sample_adata.var_names_make_unique()
    adatas[sample_id] = sample_adata

adata = ad.concat(adatas, label="sample")
adata.obs_names_make_unique()
```

### Writing Functions

```python
adata.write_h5ad(filename)                  # Write to h5ad format
adata.write_csvs(dirname)                   # Write to CSV files
adata.write_loom(filename)                  # Write to loom format
adata.write_zarr(filename)                  # Write to zarr format
```

## Preprocessing (sc.pp.*)

### Quality Control

```python
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'], inplace=True, log1p=True)
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)
```

### Doublet Detection

```python
sc.pp.scrublet(adata)                            # Basic doublet detection
sc.pp.scrublet(adata, batch_key="sample")        # With batch correction
# Adds 'doublet_score' and 'predicted_doublet' to adata.obs
```

### Normalization and Transformation

```python
sc.pp.normalize_total(adata)                     # Normalize to median total counts
sc.pp.normalize_total(adata, target_sum=1e4)     # Normalize to specific target sum
sc.pp.log1p(adata)                               # Log(x + 1) transformation
sc.pp.sqrt(adata)                                # Square root transformation
```

### Feature Selection

```python
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=2000)
```

### Scaling and Regression (Optional)

```python
sc.pp.scale(adata, max_value=10)                      # Scale to unit variance
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])  # Regress out unwanted variation
```

### Dimensionality Reduction (Preprocessing)

```python
sc.pp.pca(adata, n_comps=50)                     # Principal component analysis
sc.pp.neighbors(adata)                            # Compute neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)  # With specific parameters
```

## Tools (sc.tl.*)

### Dimensionality Reduction

```python
sc.tl.pca(adata)                                  # PCA
sc.tl.pca(adata, svd_solver='arpack')             # PCA with specific solver
sc.tl.umap(adata)                                 # UMAP embedding
sc.tl.tsne(adata)                                 # t-SNE embedding
sc.tl.diffmap(adata)                              # Diffusion map
sc.tl.draw_graph(adata, layout='fa')              # Force-directed graph
```

### Clustering

```python
sc.tl.leiden(adata)                               # Leiden clustering (default)
sc.tl.leiden(adata, flavor="igraph", n_iterations=2)  # Faster with igraph
sc.tl.leiden(adata, resolution=0.5)               # With specific resolution
sc.tl.leiden(adata, key_added="leiden_res_0.50")  # Store with custom key
sc.tl.louvain(adata, resolution=0.5)              # Louvain clustering
```

### Marker Genes and Differential Expression

```python
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test')
sc.tl.rank_genes_groups(adata, groupby='leiden', method='logreg')

# Get results as dataframe
sc.get.rank_genes_groups_df(adata, group='0')
sc.get.rank_genes_groups_df(adata, group='0').head(5)
```

### Trajectory Inference

```python
sc.tl.paga(adata, groups='leiden')               # PAGA trajectory
sc.tl.dpt(adata)                                  # Diffusion pseudotime
```

### Gene Scoring

```python
sc.tl.score_genes(adata, gene_list, score_name='score')
sc.tl.score_genes_cell_cycle(adata, s_genes, g2m_genes)
```

### Aggregation

```python
# Pseudo-bulk aggregation for differential expression
sc.get.aggregate(adata, by=["sample", "cell_type"], func="sum", layer="counts")
```

## Plotting (sc.pl.*)

### Basic Embeddings

```python
sc.pl.umap(adata, color='leiden')                # UMAP plot
sc.pl.umap(adata, color='leiden', size=2)        # With point size
sc.pl.umap(adata, color=['leiden'], legend_loc='on data')
sc.pl.tsne(adata, color='gene_name')             # t-SNE plot
sc.pl.pca(adata, color='leiden')                 # PCA plot
sc.pl.pca(adata, color=['sample', 'pct_counts_mt'], dimensions=[(0, 1), (2, 3)])
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
sc.pl.diffmap(adata, color='leiden')             # Diffusion map plot
```

### Heatmaps and Dot Plots

```python
sc.pl.heatmap(adata, var_names=genes, groupby='leiden')
sc.pl.dotplot(adata, var_names=genes, groupby='leiden')
sc.pl.dotplot(adata, marker_genes, groupby='leiden', standard_scale='var')
sc.pl.matrixplot(adata, var_names=genes, groupby='leiden')
sc.pl.stacked_violin(adata, var_names=genes, groupby='leiden')
```

### Violin and Scatter Plots

```python
sc.pl.violin(adata, keys=['gene1', 'gene2'], groupby='leiden')
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)
sc.pl.scatter(adata, x='gene1', y='gene2', color='leiden')
sc.pl.scatter(adata, 'total_counts', 'n_genes_by_counts', color='pct_counts_mt')
```

### Marker Gene Visualization

```python
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
sc.pl.rank_genes_groups_violin(adata, groups='0')
sc.pl.rank_genes_groups_heatmap(adata, n_genes=10)
sc.pl.rank_genes_groups_dotplot(adata, groupby='leiden', n_genes=5, standard_scale='var')
```

### Highly Variable Genes

```python
sc.pl.highly_variable_genes(adata)
```

### Trajectory Visualization

```python
sc.pl.paga(adata, color='leiden')                # PAGA graph
sc.pl.dpt_timeseries(adata)                      # DPT timeseries
```

### QC Plots

```python
sc.pl.highest_expr_genes(adata, n_top=20)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'])
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
```

### Saving Plots

```python
# Using save parameter
sc.pl.umap(adata, color=['leiden'], save='_leiden.pdf')

# Accessing returned figure
sc.pl.umap(adata, color=["leiden"], show=False).figure.savefig("output_path.pdf")
```

## Common Parameters

### Color Parameters
- `color`: Variable(s) to color by (gene name, obs column)
- `palette`: Color palette to use
- `vmin`, `vmax`: Color scale limits
- `cmap`: Colormap for continuous variables

### Layout Parameters
- `basis`: Embedding basis ('umap', 'tsne', 'pca', etc.)
- `legend_loc`: Legend location ('on data', 'right margin', etc.)
- `size`: Point size
- `alpha`: Point transparency
- `wspace`: Horizontal space between panels
- `ncols`: Number of columns in multi-panel plots
- `dimensions`: PCA dimensions to plot, e.g., [(0, 1), (2, 3)]

### Saving Parameters
- `save`: Filename suffix to save plot
- `show`: Whether to show plot

## AnnData Structure

```python
adata.X                    # Expression matrix (cells × genes)
adata.obs                  # Cell annotations (DataFrame)
adata.var                  # Gene annotations (DataFrame)
adata.uns                  # Unstructured annotations (dict)
adata.obsm                 # Multi-dimensional cell annotations (e.g., PCA, UMAP)
adata.varm                 # Multi-dimensional gene annotations
adata.layers               # Additional data layers (e.g., raw counts)

# Using layers to store raw counts
adata.layers["counts"] = adata.X.copy()

# Access
adata.obs_names            # Cell barcodes
adata.var_names            # Gene names
adata.shape                # (n_cells, n_genes)

# Slicing
adata[cell_indices, gene_indices]
adata[:, adata.var_names.isin(gene_list)]
adata[adata.obs['leiden'] == '0', :]
```

## Settings

```python
sc.settings.verbosity = 3              # 0=error, 1=warning, 2=info, 3=hint
sc.settings.set_figure_params(dpi=100, facecolor='white')
sc.settings.autoshow = False           # Don't show plots automatically
sc.settings.autosave = True            # Autosave figures
sc.settings.figdir = './figures/'      # Figure directory
sc.settings.cachedir = './cache/'      # Cache directory
sc.settings.n_jobs = 8                 # Number of parallel jobs
sc.settings.file_format_figs = 'pdf'   # Default figure format
```

## Useful Utilities

```python
sc.logging.print_versions()            # Print version information
sc.logging.print_memory_usage()        # Print memory usage
adata.copy()                           # Create a copy of AnnData object
adata.var_names_make_unique()          # Make gene names unique
adata.obs_names_make_unique()          # Make cell names unique
```
