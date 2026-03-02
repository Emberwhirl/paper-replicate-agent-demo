---
name: single-cell-rna-qc
description: Specialized quality control on single-cell RNA-seq data using MAD-based filtering and comprehensive visualizations. Use only when users specifically request detailed QC analysis, MAD-based outlier detection, or need specialized QC beyond what scanpy-basics provides. For standard QC as part of normal analysis, use scanpy-basics instead.
---

# Single-Cell RNA-seq Quality Control

Specialized QC workflow for single-cell RNA-seq data following scverse best practices with MAD-based filtering.

## When to Use This Skill

Use this skill **only** when users specifically request:

- MAD-based (Median Absolute Deviation) filtering or outlier detection
- Detailed QC metrics and comprehensive visualizations
- Specialized QC analysis beyond standard scanpy workflows
- Explicit scverse best practices for quality control

**For standard QC as part of a normal analysis workflow, use `scanpy-basics` instead** - it includes QC as part of the complete single-cell analysis pipeline.

**Supported input formats:**

- `.h5ad` files (AnnData format from scanpy/Python workflows)
- `.h5` files (10X Genomics Cell Ranger output)

**Default recommendation**: Use Approach 1 (complete pipeline) unless the user has specific custom requirements or explicitly requests non-standard filtering logic.

## Approach 1: Complete QC Pipeline (Recommended for Standard Workflows)

For standard QC following scverse best practices, use the convenience script `scripts/qc_analysis.py`:

```bash
python3 scripts/qc_analysis.py input.h5ad
# or for 10X Genomics .h5 files:
python3 scripts/qc_analysis.py raw_feature_bc_matrix.h5
```

The script automatically detects the file format and loads it appropriately.

**When to use this approach:**

- Standard QC workflow with adjustable thresholds (all cells filtered the same way)
- Batch processing multiple datasets
- Quick exploratory analysis
- User wants the "just works" solution

**Requirements:** anndata, scanpy, scipy, matplotlib, seaborn, numpy

**Parameters:**

Customize filtering thresholds using command-line parameters:

- `--output-dir` - Output directory
- `--mad-counts`, `--mad-genes`, `--mad-mt` - MAD thresholds for counts/genes/MT%
- `--mt-threshold` - Hard mitochondrial % cutoff
- `--min-cells` - Gene filtering threshold

Gene patterns use uppercase matching for cross-species compatibility (MT-, RPS/RPL, HB).

Use `--help` to see current default values.

**Outputs:**

All files are saved to `<input_basename>_qc_results/` directory by default (or to the directory specified by `--output-dir`):

- `qc_metrics_before_filtering.png` - Pre-filtering visualizations
- `qc_filtering_thresholds.png` - MAD-based threshold overlays
- `qc_metrics_after_filtering.png` - Post-filtering quality metrics
- `<input_basename>_filtered.h5ad` - Clean, filtered dataset ready for downstream analysis
- `<input_basename>_with_qc.h5ad` - Original data with QC annotations preserved

If copying outputs to `/mnt/user-data/outputs/` for user access, copy individual files (not the entire directory) so users can preview them directly as Claude.ai artifacts.

### Workflow Steps

The script performs the following steps:

1. **Calculate QC metrics** - Count depth, gene detection, mitochondrial/ribosomal/hemoglobin content
2. **Apply MAD-based filtering** - Permissive outlier detection using MAD thresholds for counts/genes/MT%
3. **Filter genes** - Remove genes detected in few cells
4. **Generate visualizations** - Comprehensive before/after plots with threshold overlays

## Approach 2: Modular Building Blocks (For Custom Workflows)

For custom analysis workflows or non-standard requirements, use the modular utility functions from `scripts/qc_core.py` and `scripts/qc_plotting.py`:

```python
# Run from scripts/ directory, or add scripts/ to sys.path if needed
import anndata as ad
from qc_core import calculate_qc_metrics, detect_outliers_mad, filter_cells
from qc_plotting import plot_qc_distributions  # Only if visualization needed

adata = ad.read_h5ad('input.h5ad')
calculate_qc_metrics(adata, inplace=True)
# ... custom analysis logic here
```

**When to use this approach:**

- Different workflow needed (skip steps, change order, apply different thresholds to subsets)
- Conditional logic (e.g., filter neurons differently than other cells)
- Partial execution (only metrics/visualization, no filtering)
- Integration with other analysis steps in a larger pipeline
- Custom filtering criteria beyond what command-line params support

**Available utility functions:**

From `qc_core.py` (core QC operations):

- `calculate_qc_metrics(adata, inplace=True)` - Calculate QC metrics using uppercase pattern matching for cross-species compatibility
- `detect_outliers_mad(adata, metric, n_mads, verbose=True)` - MAD-based outlier detection, returns boolean mask
- `apply_hard_threshold(adata, metric, threshold, operator='>', verbose=True)` - Apply hard cutoffs, returns boolean mask
- `filter_cells(adata, mask, inplace=False)` - Apply boolean mask to filter cells
- `filter_genes(adata, min_cells=20, min_counts=None, inplace=True)` - Filter genes by detection
- `print_qc_summary(adata, label='')` - Print summary statistics

From `qc_plotting.py` (visualization):

- `plot_qc_distributions(adata, output_path, title)` - Generate comprehensive QC plots
- `plot_filtering_thresholds(adata, outlier_masks, thresholds, output_path)` - Visualize filtering thresholds
- `plot_qc_after_filtering(adata, output_path)` - Generate post-filtering plots

**Example custom workflows:**

**Example 1: Only calculate metrics and visualize, don't filter yet**

```python
adata = ad.read_h5ad('input.h5ad')
calculate_qc_metrics(adata, inplace=True)
plot_qc_distributions(adata, 'qc_before.png', title='Initial QC')
print_qc_summary(adata, label='Before filtering')
```

**Example 2: Apply only MT% filtering, keep other metrics permissive**

```python
adata = ad.read_h5ad('input.h5ad')
calculate_qc_metrics(adata, inplace=True)

# Only filter high MT% cells
high_mt = apply_hard_threshold(adata, 'pct_counts_mt', 10, operator='>')
adata_filtered = filter_cells(adata, ~high_mt)
adata_filtered.write('filtered.h5ad')
```

**Example 3: Different thresholds for different subsets**

```python
adata = ad.read_h5ad('input.h5ad')
calculate_qc_metrics(adata, inplace=True)

# Apply type-specific QC (assumes cell_type metadata exists)
neurons = adata.obs['cell_type'] == 'neuron'
other_cells = ~neurons

# Neurons tolerate higher MT%, other cells use stricter threshold
neuron_qc = apply_hard_threshold(adata[neurons], 'pct_counts_mt', 15, operator='>')
other_qc = apply_hard_threshold(adata[other_cells], 'pct_counts_mt', 8, operator='>')
```

## Best Practices

1. **Be permissive with filtering** - Default thresholds intentionally retain most cells to avoid losing rare populations
2. **Inspect visualizations** - Always review before/after plots to ensure filtering makes biological sense
3. **Consider dataset-specific factors** - Some tissues naturally have higher mitochondrial content (e.g., neurons, cardiomyocytes)
4. **Check gene annotations** - Mitochondrial gene prefixes vary by species (mt- for mouse, MT- for human)
5. **Iterate if needed** - QC parameters may need adjustment based on the specific experiment or tissue type

## Reference Materials

For detailed QC methodology, parameter rationale, and troubleshooting guidance, see `references/scverse_qc_guidelines.md`. This reference provides:

- Detailed explanations of each QC metric and why it matters
- Rationale for MAD-based thresholds and why they're better than fixed cutoffs
- Guidelines for interpreting QC visualizations (histograms, violin plots, scatter plots)
- Species-specific considerations for gene annotations
- When and how to adjust filtering parameters
- Advanced QC considerations (ambient RNA correction, doublet detection)

Load this reference when users need deeper understanding of the methodology or when troubleshooting QC issues.

## Integration with Other Skills

### Workflow Positioning

This skill handles **specialized QC** (MAD-based filtering) as a dedicated step:

- **Before this skill**: Use `kallisto-bustools-import` if data is from kb-python with unfiltered empty droplets
- **After this skill**: Continue with `scanpy-basics` for normalization, clustering, and analysis

### Standard Analysis Workflows

**For kallisto/bustools data (unfiltered):**
1. **kallisto-bustools-import** - Load data, filter empty droplets
2. **single-cell-rna-qc** (this skill) - MAD-based QC (if specialized QC requested)
3. **scanpy-basics** - Normalization, clustering, annotation

**For CellRanger or pre-filtered data:**
1. **single-cell-rna-qc** (this skill) - MAD-based QC (if specialized QC requested)
2. **scanpy-basics** - Normalization, clustering, annotation

**For standard analysis without specialized QC requests:**
- Start directly with **scanpy-basics** - includes built-in QC

### Consistent QC Metrics

All skills use the same uppercase pattern matching for cross-species compatibility:

```python
# Mitochondrial genes
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
# Ribosomal genes
adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
# Hemoglobin genes (excluding HBP)
adata.var["hb"] = adata.var_names.str.upper().str.contains(r"^HB(?!P)")
```

This ensures consistent results across all skills in the repository.

## Next Steps After QC

After specialized QC, continue with standard analysis:

- Use `scanpy-basics` skill for normalization, dimensionality reduction, clustering, and annotation

**If additional specialized tasks are requested:**

- **Batch correction** for combining datasets → use `single-cell-integration-ingest-bbknn` skill
- **Deep learning integration** (scVI, scANVI) → use `scvi-tools` skill
