# Integration Methods: Ingest and BBKNN

This document provides detailed information about single-cell data integration methods using scanpy's ingest and BBKNN, including algorithm details, parameter tuning, and advanced scenarios.

## Algorithm Overview

### Ingest (sc.tl.ingest)

Ingest performs asymmetric, reference-based integration where query cells are mapped onto a pre-computed reference embedding.

**How it works:**

1. **Training phase**: Reference dataset has PCA, neighbors graph, and UMAP pre-computed
2. **Projection**: Query cells are projected into reference PCA space using trained loadings
3. **Neighbor assignment**: Each query cell finds its nearest neighbors among reference cells
4. **Label transfer**: Cell type labels are transferred via weighted k-NN voting
5. **Embedding**: Query cells are placed in reference UMAP using the same neighbor relationships

**Key characteristics:**

- Reference embedding is unchanged
- Query cells inherit reference structure
- Built-in label transfer mechanism
- Fast and deterministic

### BBKNN (Batch Balanced K-Nearest Neighbors)

BBKNN performs symmetric integration by modifying how neighbors are selected across batches.

**How it works:**

1. **PCA**: Standard dimensionality reduction on combined data
2. **Batch-aware neighbors**: Instead of finding k nearest neighbors globally, find k/n_batches neighbors from each batch
3. **Graph construction**: Build connectivity graph with balanced batch representation
4. **Clustering/UMAP**: Downstream methods operate on the corrected graph

**Key characteristics:**

- All batches treated equally
- Modifies neighbor graph, not expression data
- Strong batch mixing
- Computationally efficient

## When to Use Each Method

### Choose Ingest When:

- You have a well-annotated reference atlas
- You want to preserve reference structure (clusters, trajectories)
- Primary goal is label transfer
- Reference has comprehensive cell type coverage
- You're iteratively adding new samples

### Choose BBKNN When:

- No clear reference dataset exists
- All datasets should contribute equally
- You need maximum batch mixing
- Working with many batches (>5)
- Batch effects are strong but you trust the cell types

### Consider Other Methods When:

- Complex cross-technology batch effects → scVI, sysVI (scvi-tools skill)
- Need uncertainty quantification → scANVI (scvi-tools skill)
- Very large datasets (>1M cells) → scVI, Harmony
- Need to learn biological variation → VAE-based methods

## Parameter Tuning

### Ingest Parameters

**n_pcs (default: 50)**
- Number of principal components for training
- Higher values capture more variance but may include noise
- Recommendation: Use variance explained plot to choose; 30-50 usually sufficient

**n_neighbors (default: 15)**
- Number of neighbors for graph construction
- Affects cluster granularity and UMAP layout
- Lower values: Tighter clusters, more local structure
- Higher values: Smoother embedding, more global structure

**embedding_method (default: 'umap')**
- Method for embedding query cells: 'umap' or 'pca'
- 'umap': Query cells placed in reference UMAP space
- 'pca': Only PCA projection, UMAP computed fresh

### BBKNN Parameters

**Note:** BBKNN is a separate package (`pip install bbknn`). Import directly as `import bbknn` — the `sc.external.pp.bbknn` wrapper is deprecated in scanpy ≥1.10.

**neighbors_within_batch (default: 3)**
- Number of neighbors from each batch
- Total neighbors = neighbors_within_batch × n_batches
- Higher values: More batch mixing, potential over-integration
- Lower values: Less mixing, may retain some batch structure

```python
# Adjust for different scenarios
import bbknn
bbknn.bbknn(adata, batch_key='batch',
            neighbors_within_batch=5)  # Stronger mixing
bbknn.bbknn(adata, batch_key='batch',
            neighbors_within_batch=2)  # Gentler mixing
```

**n_pcs (default: 50)**
- Number of PCs to use for neighbor finding
- Same considerations as ingest

**trim (optional)**
- Trim graph edges based on connectivities
- Can help with noise and outliers

## Handling Complex Scenarios

### Multiple Query Batches

For integrating multiple query datasets onto a reference:

```python
from integration_core import intersect_genes, train_reference, ingest_query
import anndata as ad

# Load reference
adata_ref = ad.read_h5ad('reference.h5ad')
train_reference(adata_ref)

# Process each query
queries = ['query1.h5ad', 'query2.h5ad', 'query3.h5ad']
integrated = [adata_ref]

for path in queries:
    adata_query = ad.read_h5ad(path)
    _, adata_query = intersect_genes(adata_ref, adata_query)
    ingest_query(adata_query, adata_ref, label_key='celltype')
    integrated.append(adata_query)

# Combine
adata_combined = ad.concat(integrated, label='source')
```

### Mixed Integration Strategies

Use ingest for label transfer, then BBKNN for visualization:

```python
# First: Transfer labels using ingest
ingest_query(adata_query, adata_ref, label_key='celltype')

# Second: Combine and apply BBKNN for visualization
adata_combined = ad.concat([adata_ref, adata_query], label='batch')
sc.pp.pca(adata_combined)
import bbknn
bbknn.bbknn(adata_combined, batch_key='batch')
sc.tl.umap(adata_combined)
```

### Handling Gene Name Differences

When datasets use different gene naming conventions:

```python
import re

def harmonize_gene_names(adata, style='symbol'):
    """
    Harmonize gene names between datasets.
    
    style: 'symbol' for gene symbols, 'ensembl' for ENSEMBL IDs
    """
    if style == 'symbol':
        # Remove version numbers from ENSEMBL IDs
        adata.var_names = [re.sub(r'\.\d+$', '', g) for g in adata.var_names]
    
    # Convert to uppercase for case-insensitive matching
    adata.var_names = adata.var_names.str.upper()
    
    # Make unique
    adata.var_names_make_unique()
    
    return adata
```

### Novel Cell Types in Query

When query contains cell types not present in reference:

1. **Identify unmapped cells**: Look for query cells distant from all reference cells
2. **Cluster separately**: Run clustering on query subset
3. **Annotate manually**: Use marker genes to identify novel populations

```python
# After ingest, check embedding distances
distances = sc.pp.neighbors(adata_query, n_neighbors=10, key_added='ref_distances')
# High average distance → potentially novel cell type
```

## Evaluating Integration Quality

### Quantitative Metrics

**Label transfer accuracy** (for ingest with ground truth):

```python
from integration_core import evaluate_label_transfer

confusion = evaluate_label_transfer(adata_query, 'celltype', 'celltype_orig')
accuracy = confusion.values.diagonal().sum() / confusion.values.sum() * 100
print(f"Accuracy: {accuracy:.1f}%")
```

**Batch mixing metrics** (for BBKNN):

```python
# Average Silhouette Width (ASW) for batch
from sklearn.metrics import silhouette_score

asw_batch = silhouette_score(adata.obsm['X_pca'], adata.obs['batch'])
# Closer to 0 is better (batches well-mixed)
print(f"Batch ASW: {asw_batch:.3f}")

# ASW for cell types
asw_celltype = silhouette_score(adata.obsm['X_pca'], adata.obs['celltype'])
# Closer to 1 is better (cell types separated)
print(f"Cell type ASW: {asw_celltype:.3f}")
```

**kBET (k-Nearest Neighbor Batch Effect Test)**:
External tool for rigorous batch effect assessment. See: https://github.com/theislab/kBET

### Qualitative Assessment

**Visual inspection checklist:**

1. Do cell types form coherent clusters across batches?
2. Are batches well-mixed within cell type clusters?
3. Are rare populations preserved?
4. Do known biological differences remain visible?

**UMAP overlays:**

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Before integration (if available)
sc.pl.umap(adata_before, color='batch', ax=axes[0], title='Before', show=False)

# After integration
sc.pl.umap(adata_after, color='batch', ax=axes[1], title='After', show=False)

plt.tight_layout()
plt.savefig('integration_comparison.png', dpi=150)
```

## Troubleshooting

### Common Issues

**Issue: No common genes between datasets**
- Cause: Different gene naming (symbols vs ENSEMBL, species differences)
- Solution: Harmonize gene names before intersection

**Issue: Poor label transfer (low accuracy)**
- Cause: Query contains novel populations, or reference lacks coverage
- Solutions:
  - Expand reference with more cell types
  - Use BBKNN + de novo clustering instead
  - Check for strong batch effects requiring deep learning

**Issue: Over-integration (losing biological signal)**
- Cause: BBKNN mixing too aggressively
- Solutions:
  - Reduce `neighbors_within_batch`
  - Use ingest instead if reference is available
  - Check if datasets are actually comparable

**Issue: Batches still cluster separately after BBKNN**
- Cause: Very strong technical effects
- Solutions:
  - Increase `neighbors_within_batch`
  - Consider scVI/scANVI for stronger correction
  - Check for fundamental data quality issues

**Issue: Query cells don't map well to reference**
- Cause: Query has cell types absent from reference
- Solutions:
  - Identify and exclude novel populations
  - Expand reference dataset
  - Use BBKNN for those samples instead

## References

- Wolf et al. (2018): PAGA: Graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells
- Polański et al. (2020): BBKNN: fast batch alignment of single cell transcriptomes
- Luecken et al. (2022): Benchmarking atlas-level data integration in single-cell genomics

## See Also

- `scanpy-basics` skill: For standard single-cell analysis workflows
- `scvi-tools` skill: For deep learning integration methods (scVI, scANVI)
- `single-cell-rna-qc` skill: For quality control before integration

## Software dependencies & licensing

This skill *orchestrates* third-party open-source libraries that are **not** redistributed here — install them yourself. They remain under their own upstream licenses:

- **scanpy**, **anndata** — BSD-3-Clause
- **bbknn** — MIT

The skill itself (this documentation, the scripts, and the workflow) is MIT-licensed; see the `LICENSE` file in the skill root.
