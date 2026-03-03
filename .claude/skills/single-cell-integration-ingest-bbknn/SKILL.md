---
name: single-cell-integration-ingest-bbknn
description: Integrates multiple single-cell RNA-seq datasets to correct batch effects while preserving biological variation. Use when users need to combine datasets from different experiments, map query data onto a reference, transfer cell type labels, or visualize batch-corrected embeddings. Supports asymmetric integration with ingest (reference-based label transfer) and symmetric integration with BBKNN (batch-balanced neighbors). For deep learning integration methods, use scvi-tools instead.
argument-hint: "[reference.h5ad query.h5ad or combined.h5ad --method bbknn]"
allowed-tools: ["Read", "Write", "Edit", "Bash", "Grep", "Glob"]
---

# Single-Cell RNA-seq Data Integration

Automated workflows for integrating multiple single-cell datasets using scanpy's ingest and BBKNN methods.

## When to Use This Skill

Use when users:

- Want to combine/merge multiple scRNA-seq datasets
- Need to correct batch effects between experiments
- Want to map query cells onto an annotated reference
- Need to transfer cell type labels from reference to query data
- Ask about ingest, BBKNN, or batch correction in scanpy
- Want to compare cells across experimental conditions or studies
- Need to evaluate integration quality with confusion matrices

**Supported integration methods:**

- `sc.tl.ingest` - Asymmetric, reference-based integration with label transfer
- `bbknn.bbknn` - Symmetric batch-balanced k-nearest neighbors (via the `bbknn` package)

**When to use other skills instead:**

- For deep learning integration (scVI, scANVI, Harmony) → use `scvi-tools` skill
- For general scanpy analysis → use `scanpy-basics` skill
- For quality control before integration → use `single-cell-rna-qc` skill

**Default recommendation**: Use Approach 1 (ingest) when you have an annotated reference dataset. Use Approach 2 (BBKNN) when all batches should be treated symmetrically without a designated reference.

## Approach 1: Reference-Based Integration with Ingest (Recommended for Label Transfer)

Use `sc.tl.ingest` when you have an annotated reference dataset and want to:

- Map query cells onto the reference embedding
- Transfer cell type labels from reference to query
- Maintain the structure of a known biological trajectory

For standard reference-based integration, use the convenience script `scripts/integrate_ingest.py`:

```bash
python3 scripts/integrate_ingest.py reference.h5ad query.h5ad --output integrated.h5ad
```

**When to use this approach:**

- Well-annotated reference dataset exists
- Want to preserve reference structure and transfer labels
- Iterative/exploratory analysis with fast turnaround
- Want to maintain specific clusters or trajectories from reference

**Requirements:** anndata, scanpy

**Parameters:**

- `reference` - Reference h5ad file with annotations
- `query` - Query h5ad file to map onto reference
- `--output`, `-o` - Output file path (default: integrated.h5ad)
- `--label-key` - Observation column to transfer (default: celltype)
- `--output-dir` - Output directory for plots and results
- `--skip-plots` - Skip generating visualization plots

Use `--help` to see current default values.

**Outputs:**

- `integrated.h5ad` - Combined dataset with transferred labels and embeddings
- `umap_reference.png` - UMAP visualization of reference dataset
- `umap_query.png` - UMAP visualization of mapped query cells
- `umap_combined.png` - Joint visualization of reference and query

### Workflow Steps

The script performs the following steps:

1. **Load datasets** - Read reference and query h5ad files
2. **Intersect genes** - Find common gene space between datasets
3. **Train on reference** - Compute PCA, neighbors, and UMAP on reference
4. **Map query** - Project query cells onto reference embedding and transfer labels
5. **Concatenate** - Combine datasets for joint visualization
6. **Save results** - Export integrated data and visualizations

## Approach 2: Symmetric Integration with BBKNN

Use BBKNN when you want symmetric treatment of all batches without a designated reference:

```bash
python3 scripts/integrate_bbknn.py combined_data.h5ad --output integrated_bbknn.h5ad
```

**When to use this approach:**

- No clear reference dataset available
- Want maximum batch mixing
- Large number of batches to integrate
- All datasets should contribute equally to the embedding

**Requirements:** anndata, scanpy, bbknn (`pip install bbknn`)

**Parameters:**

- `input` - Combined h5ad file with batch column
- `--output`, `-o` - Output file path (default: integrated_bbknn.h5ad)
- `--batch-key` - Observation column indicating batch (default: batch)
- `--neighbors-within-batch` - Number of neighbors within each batch (default: 3)
- `--output-dir` - Output directory for plots and results
- `--skip-plots` - Skip generating visualization plots

**Outputs:**

- `integrated_bbknn.h5ad` - Dataset with batch-corrected neighbors graph
- `umap_bbknn.png` - UMAP visualization colored by batch and cell type

### Workflow Steps

The script performs the following steps:

1. **Load combined data** - Read h5ad with batch annotations
2. **Compute PCA** - Dimensionality reduction on combined data
3. **Apply BBKNN** - Replace standard neighbors with batch-balanced neighbors
4. **Compute UMAP** - Visualization on corrected graph
5. **Save results** - Export integrated data and visualizations

## Approach 3: Modular Building Blocks (For Custom Workflows)

For custom analysis or integration with other pipelines, use the modular utility functions from `scripts/integration_core.py`:

```python
# Run from scripts/ directory, or add scripts/ to sys.path if needed
import anndata as ad
from integration_core import (
    intersect_genes,
    train_reference,
    ingest_query,
    evaluate_label_transfer,
    combine_datasets
)

adata_ref = ad.read_h5ad('reference.h5ad')
adata_query = ad.read_h5ad('query.h5ad')

# Intersect to common genes
adata_ref, adata_query = intersect_genes(adata_ref, adata_query)

# Train on reference
train_reference(adata_ref)

# Map query
ingest_query(adata_query, adata_ref, label_key='celltype')

# Evaluate if ground truth available
confusion = evaluate_label_transfer(adata_query, 'celltype', 'celltype_orig')
```

**When to use this approach:**

- Custom workflow needed (skip steps, change order)
- Integration with other analysis pipelines
- Iterative integration of multiple batches
- Custom evaluation logic

**Available utility functions:**

From `integration_core.py`:

- `intersect_genes(adata_ref, adata_query)` - Find common gene space
- `train_reference(adata, n_pcs=50)` - Compute PCA, neighbors, UMAP on reference
- `ingest_query(adata_query, adata_ref, label_key)` - Map query and transfer labels
- `integrate_bbknn(adata, batch_key, neighbors_within_batch)` - Apply BBKNN
- `combine_datasets(adata_ref, adata_query, label_key)` - Concatenate with proper categories
- `evaluate_label_transfer(adata, transferred_key, original_key)` - Compute confusion matrix
- `print_integration_summary(adata, label)` - Print summary statistics

**Example custom workflows:**

**Example 1: Iterative integration of multiple batches**

```python
from integration_core import intersect_genes, train_reference, ingest_query
import anndata as ad

# Load all data with batch annotations
adata_all = ad.read_h5ad('all_samples.h5ad')

# Use batch 0 as reference
adata_ref = adata_all[adata_all.obs['batch'] == '0'].copy()
train_reference(adata_ref)

# Process each query batch
query_batches = ['1', '2', '3']
adatas = [adata_ref]

for batch_id in query_batches:
    adata = adata_all[adata_all.obs['batch'] == batch_id].copy()
    _, adata = intersect_genes(adata_ref, adata)
    adata.obs['celltype_orig'] = adata.obs['celltype']
    ingest_query(adata, adata_ref, label_key='celltype')
    adatas.append(adata)
    print(f"Integrated batch {batch_id}")

# Combine all batches
adata_integrated = ad.concat(adatas, label='batch', join='outer')
```

**Example 2: Evaluate integration quality with confusion matrix**

```python
from integration_core import evaluate_label_transfer
import pandas as pd

# After ingest, compare transferred vs original labels
confusion = evaluate_label_transfer(
    adata_query, 
    transferred_key='celltype',  # from ingest
    original_key='celltype_orig'  # ground truth
)
print(confusion)

# Calculate accuracy
accuracy = (confusion.values.diagonal().sum() / confusion.values.sum()) * 100
print(f"Label transfer accuracy: {accuracy:.1f}%")
```

**Example 3: Chain with QC pipeline**

```python
# After QC filtering (single-cell-rna-qc skill), integrate datasets
from integration_core import intersect_genes, train_reference, ingest_query

# Load QC-filtered datasets
adata_ref = ad.read_h5ad('reference_filtered.h5ad')
adata_query = ad.read_h5ad('query_filtered.h5ad')

# Proceed with integration
adata_ref, adata_query = intersect_genes(adata_ref, adata_query)
train_reference(adata_ref)
ingest_query(adata_query, adata_ref, label_key='celltype')
```

## Method Selection Guide

| Scenario | Recommended Method | Rationale |
|----------|-------------------|-----------|
| Well-annotated reference exists | **ingest** | Preserves reference structure, transfers labels |
| Want to maintain specific clusters/trajectories | **ingest** | Keeps reference embedding intact |
| No clear reference dataset | **BBKNN** | Symmetric treatment of all batches |
| Need maximum batch mixing | **BBKNN** | Better homogeneous mixing |
| Iterative/exploratory analysis | **ingest** | Fast, transparent workflow |
| Large number of batches | **BBKNN** | Handles many batches well |

## Best Practices

1. **Always intersect genes first** - Both methods require shared gene space between datasets
2. **QC before integration** - Run quality control on each dataset separately before combining
3. **Choose reference carefully** - For ingest, select reference with comprehensive cell type coverage
4. **Preserve original annotations** - Save original labels before transferring (`obs['celltype_orig'] = obs['celltype']`)
5. **Validate with confusion matrices** - Compare transferred labels against original annotations
6. **Consider cell type coverage** - Cells mapping to types absent in reference may be incorrectly labeled
7. **Check for over-correction** - Ensure biological differences aren't removed along with batch effects

## Common Pitfalls

| Issue | Cause | Solution |
|-------|-------|----------|
| KeyError on gene names | Gene naming differs between datasets | Harmonize gene names before intersection |
| Missing cell types in query | Reference lacks those populations | Use BBKNN or expand reference dataset |
| Poor cluster separation after BBKNN | Over-integration | Adjust `neighbors_within_batch` parameter |
| Cells from query cluster separately | Strong batch effect or novel populations | Increase reference size or use deep learning methods |

## Method Comparison Summary

| Feature | ingest | BBKNN |
|---------|--------|-------|
| Integration type | Asymmetric (reference → query) | Symmetric |
| Label transfer | ✓ Built-in | ✗ Requires separate step |
| Embedding preservation | ✓ Maintains reference structure | ✗ Recomputes embedding |
| Data matrix | Unchanged | Unchanged |
| Batch mixing | Moderate | Strong |
| Speed | Fast | Fast |
| Best for | Label transfer, trajectory preservation | Maximum mixing, no clear reference |

## Integration with Other Skills

### When to Use This Skill

Use this skill when **batch correction** is explicitly requested for integrating multiple single-cell datasets.

### Standard Analysis Workflow

For most analyses without batch correction needs:
- **CellRanger output or pre-filtered data** → Start directly with **scanpy-basics**
- **kallisto/bustools unfiltered data** → **kallisto-bustools-import** → **scanpy-basics**

### When Batch Correction is Requested

1. Load and preprocess data (use **kallisto-bustools-import** only if needed for kb-python unfiltered data, otherwise start with **scanpy-basics**)
2. **single-cell-integration-ingest-bbknn** (this skill) - Batch correction and integration
3. Continue with downstream analysis in **scanpy-basics**

### When to Use scvi-tools Instead

Use `scvi-tools` skill instead of this skill when:
- Deep learning-based integration is specifically requested (scVI, scANVI)
- Complex batch effects require probabilistic modeling
- Multi-modal data integration is needed

### Passing Data Between Skills

```python
# After QC, integrate
# python3 scripts/integrate_ingest.py ref_filtered.h5ad query_filtered.h5ad

# Continue with scanpy analysis
import scanpy as sc
adata = sc.read_h5ad('integrated.h5ad')
sc.tl.leiden(adata)
sc.pl.umap(adata, color=['leiden', 'celltype'])
```

## Reference Materials

For detailed methodology, parameter tuning, and advanced integration scenarios, see `references/integration_methods.md`. This reference provides:

- Detailed explanations of ingest and BBKNN algorithms
- Parameter tuning guidelines
- Handling complex integration scenarios
- Benchmarking integration quality
- When to escalate to deep learning methods

## Next Steps After Integration

- Clustering on integrated data (`sc.tl.leiden`)
- Differential expression between conditions
- Trajectory inference on combined dataset
- Cell type annotation refinement
- For complex batch effects, consider deep learning methods (scVI, scANVI) via `scvi-tools` skill
