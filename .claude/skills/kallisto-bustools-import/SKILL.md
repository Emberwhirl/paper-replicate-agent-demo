---
name: kallisto-bustools-import
description: Imports kallisto/bustools (kb-python) workflow output and performs empty droplet filtering using the DropletUtils knee/inflection algorithm. Use when users have run kb count, kb-python, or kallisto/bustools and need to load the output, assess library saturation, generate knee plots, or filter empty droplets before downstream analysis.
argument-hint: "[path to h5ad or counts_unfiltered/ directory]"
allowed-tools: ["Read", "Write", "Edit", "Bash", "Grep", "Glob"]
---

# Kallisto/Bustools Import and Empty Droplet Filtering

Automated workflow for importing kb-python/kallisto/bustools output and filtering empty droplets using the DropletUtils barcodeRanks algorithm.

## When to Use This Skill

Use when users:

- Have run `kb count` or kallisto/bustools pipeline and need to import results
- Want to load `adata.h5ad` or MTX files from `counts_unfiltered/` directory
- Need to assess library saturation or sequencing depth
- Want to generate and interpret knee plots
- Need to filter empty droplets before downstream analysis
- Ask about UMI thresholds, knee/inflection points, or cell calling

**Supported input formats:**

- `.h5ad` files from kb-python output (typically in `counts_unfiltered/` or `counts_filtered/`)
- MTX format files (cells_x_genes.mtx, .barcodes.txt, .genes.txt)
- Standard AnnData format from kallisto/bustools workflows

**Default recommendation**: Use Approach 1 (complete pipeline) for standard empty droplet filtering. Use Approach 2 for custom thresholds, batch processing, or integration with other QC steps.

## Approach 1: Complete Import and Filtering Pipeline (Recommended)

For standard kb-python output processing, use the convenience script `scripts/kb_import.py`:

```bash
# From h5ad file
python3 scripts/kb_import.py counts_unfiltered/adata.h5ad

# From MTX files directly
python3 scripts/kb_import.py counts_unfiltered/ --from-mtx
```

**When to use this approach:**

- Standard kb-python output processing
- Automatic knee/inflection threshold detection using DropletUtils algorithm
- Quick exploratory analysis of kallisto/bustools results
- User wants the "just works" solution

**Requirements:** anndata, scanpy, scipy, matplotlib, numpy, pandas, scikit-learn

**Parameters:**

Customize filtering and visualization using command-line parameters:

- `--output-dir` - Output directory (default: `<input>_kb_import_results/`)
- `--from-mtx` - Read from MTX files instead of h5ad
- `--min-genes` - Minimum genes per cell for initial filter (default: None, leave for downstream QC)
- `--umi-threshold` - Manual UMI threshold override (auto-detected if not specified)
- `--expected-cells` - Expected number of cells (helps refine threshold detection)
- `--use-knee` - Use knee threshold instead of inflection (more conservative)
- `--lower` - Lower bound on UMI count for detection (default: 100)
- `--skip-saturation` - Skip library saturation plot
- `--skip-pca` - Skip initial PCA visualization

Use `--help` to see current default values.

**Outputs:**

All files are saved to `<input>_kb_import_results/` directory by default:

- `pca_embedding.png` - Initial 2D PCA projection of all barcodes
- `library_saturation.png` - Genes detected vs UMI counts (saturation assessment)
- `knee_plot.png` - Knee plot with both knee (blue) and inflection (green) thresholds
- `<input>_filtered.h5ad` - Filtered dataset ready for downstream QC and analysis
- `filtering_summary.txt` - Summary statistics including both thresholds

If copying outputs to `/mnt/user-data/outputs/` for user access, copy individual files (not the entire directory) so users can preview them directly.

### Workflow Steps

The script performs the following steps:

1. **Load kb-python output** - Read unfiltered count matrix (h5ad or MTX format)
2. **Initial visualization** - PCA embedding to visualize barcode distribution
3. **Library saturation check** - Assess if sequencing depth is sufficient
4. **Knee plot analysis** - Identify knee and inflection thresholds using DropletUtils algorithm
5. **Apply filtering** - Remove empty droplets based on inflection threshold (more permissive)
6. **Save results** - Export filtered data and visualizations

## Approach 2: Modular Building Blocks (For Custom Workflows)

For custom analysis or integration with other QC steps, use the modular utility functions from `scripts/kb_core.py` and `scripts/kb_plotting.py`:

```python
# Run from scripts/ directory, or add scripts/ to sys.path if needed
import anndata as ad
from kb_core import (
    load_kb_output,
    read_mtx_output,
    calculate_barcode_ranks,
    detect_knee_threshold,
    filter_empty_droplets
)
from kb_plotting import plot_knee, plot_saturation, plot_pca_embedding

# Load data
adata = load_kb_output('counts_unfiltered/adata.h5ad')
# Or from MTX files:
# adata = read_mtx_output('counts_unfiltered/')

# Calculate barcode ranks and thresholds using DropletUtils algorithm
df_ranks, knee, inflection = calculate_barcode_ranks(adata)
print(f"Knee: {knee:.0f} UMIs, Inflection: {inflection:.0f} UMIs")

# Filter using inflection (more permissive) or knee (more conservative)
adata_filtered = filter_empty_droplets(adata, umi_threshold=inflection)
```

**When to use this approach:**

- Custom threshold determination logic
- Integration with other QC pipelines (e.g., single-cell-rna-qc skill)
- Batch processing multiple samples
- Partial execution (only visualization, manual threshold)
- Combining with ambient RNA correction workflows

**Available utility functions:**

From `kb_core.py` (core operations):

- `load_kb_output(path)` - Load adata.h5ad from kb-python output
- `read_mtx_output(count_dir, name='cells_x_genes')` - Read MTX format files
- `calculate_barcode_ranks(adata, lower=100, ...)` - Calculate ranks and detect knee/inflection points (DropletUtils algorithm)
- `detect_knee_threshold(adata, expected_cells=None, use_inflection=True)` - Auto-detect threshold
- `filter_empty_droplets(adata, umi_threshold, min_genes=None)` - Apply threshold filtering
- `calculate_saturation_metrics(adata)` - Compute genes detected vs UMI counts
- `print_import_summary(adata, label='')` - Print summary statistics

From `kb_plotting.py` (visualization):

- `plot_pca_embedding(adata, output_path)` - 2D PCA projection of barcodes
- `plot_saturation(adata, output_path)` - Library saturation scatter plot
- `plot_knee(adata, knee_threshold, inflection_threshold, output_path)` - Knee plot with both thresholds
- `plot_knee_interactive(adata)` - Interactive threshold selection (Jupyter)

**Example custom workflows:**

**Example 1: Visualize before deciding on threshold**

```python
import anndata as ad
from kb_core import load_kb_output, calculate_barcode_ranks
from kb_plotting import plot_knee, plot_saturation

adata = load_kb_output('counts_unfiltered/adata.h5ad')
plot_saturation(adata, 'saturation.png')

# Calculate thresholds but don't filter yet
df_ranks, knee, inflection = calculate_barcode_ranks(adata)
plot_knee(adata, knee_threshold=knee, inflection_threshold=inflection,
          output_path='knee_explore.png')
# Inspect plots, then decide on threshold
```

**Example 2: Use expected cell count for threshold detection**

```python
from kb_core import load_kb_output, detect_knee_threshold, filter_empty_droplets

adata = load_kb_output('counts_unfiltered/adata.h5ad')
# Use expected cells from experimental design
threshold, _ = detect_knee_threshold(adata, expected_cells=10000)
adata_filtered = filter_empty_droplets(adata, threshold)
print(f"Retained {adata_filtered.n_obs} cells with threshold {threshold:.0f}")
```

**Example 3: Chain with downstream QC (scanpy-basics or single-cell-rna-qc)**

```python
import scanpy as sc
from kb_core import load_kb_output, calculate_barcode_ranks, filter_empty_droplets

# Step 1: Load and filter empty droplets
adata = load_kb_output('counts_unfiltered/adata.h5ad')
_, knee, inflection = calculate_barcode_ranks(adata)
adata = filter_empty_droplets(adata, umi_threshold=inflection)

# Step 2: Add QC metrics (consistent with scanpy-basics)
# Normalize gene names for cross-species compatibility
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.upper().str.contains(r"^HB(?!P)")

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

# Step 3: Save for downstream analysis
adata.write('adata_kb_filtered.h5ad')
```

**Example 4: Batch process multiple samples and concatenate**

```python
import re
from pathlib import Path
import anndata as ad
from kb_core import load_kb_output, calculate_barcode_ranks, filter_empty_droplets
from kb_plotting import plot_knee
import scanpy as sc

# Define sample directories
data_dir = Path('./data')
sample_dirs = sorted([d for d in data_dir.iterdir() if d.is_dir() and d.name.startswith("kb_")])
samples = [re.sub(r"^kb_standard_", "", d.name) for d in sample_dirs]

filtered_adatas = []

for sample_dir, sample_id in zip(sample_dirs, samples):
    print(f"Processing: {sample_id}")
    count_dir = sample_dir / "counts_unfiltered"
    
    # Load data (h5ad or MTX)
    h5ad_path = count_dir / "adata.h5ad"
    if h5ad_path.exists():
        adata = ad.read_h5ad(h5ad_path)
    else:
        from kb_core import read_mtx_output
        adata = read_mtx_output(count_dir)
    
    # Calculate thresholds and filter
    df_ranks, knee, inflection = calculate_barcode_ranks(adata, lower=10)
    adata_filt = filter_empty_droplets(adata, umi_threshold=inflection)
    
    # Remove genes with zero counts
    sc.pp.filter_genes(adata_filt, min_cells=1)
    
    # Add sample metadata
    adata_filt.obs['sample'] = sample_id
    
    # Save knee plot
    plot_knee(adata, knee_threshold=knee, inflection_threshold=inflection,
              output_path=f'knee_plots/{sample_id}_knee.png', title=sample_id)
    
    filtered_adatas.append(adata_filt)
    print(f"  {adata_filt.n_obs:,} cells retained")

# Concatenate all samples
adata_merged = ad.concat(filtered_adatas, label='orig.ident', keys=samples, index_unique='_')
print(f"Merged: {adata_merged.n_obs:,} cells x {adata_merged.n_vars:,} genes")

# Add QC metrics for downstream analysis
adata_merged.var["mt"] = adata_merged.var_names.str.upper().str.startswith("MT-")
adata_merged.var["ribo"] = adata_merged.var_names.str.upper().str.startswith(("RPS", "RPL"))
adata_merged.var["hb"] = adata_merged.var_names.str.upper().str.contains(r"^HB(?!P)")

adata_merged.write('adata_merged.h5ad')
```

## Understanding Knee vs Inflection Thresholds

The DropletUtils algorithm identifies two key points on the barcode rank plot:

### Knee Point (Blue, More Conservative)
- Represents the transition from real cells to empty droplets
- Found where the curve shows maximum curvature
- Retains fewer cells with higher UMI counts
- Use when you want to be more conservative

### Inflection Point (Green, More Permissive)
- Represents where the gradient is steepest
- Typically yields a lower UMI threshold
- Retains more cells, including those with lower UMI counts
- **Recommended default** - allows downstream QC to make final decisions

### When to Use Each

| Scenario | Recommended Threshold |
|----------|----------------------|
| Standard analysis | Inflection (default) |
| High ambient RNA contamination | Knee |
| Precious rare cell populations | Inflection |
| Conservative initial filtering | Knee |
| Will apply rigorous downstream QC | Inflection |

## Understanding Other Outputs

### PCA Embedding Plot

- Shows all barcodes in 2D PCA space
- Helps visualize the overall structure before filtering
- Dense clusters typically represent real cells

### Library Saturation Plot

- X-axis: Total UMI counts per barcode
- Y-axis: Number of genes detected
- **Saturated library**: Curve plateaus (more UMIs don't yield more genes)
- **Unsaturated library**: Linear relationship continues (may need deeper sequencing)

## Best Practices

1. **Use inflection threshold by default** - More permissive, lets downstream QC refine filtering
2. **Check library saturation first** - If unsaturated, consider deeper sequencing
3. **Inspect knee plot carefully** - Verify thresholds make biological sense
4. **Use expected cell count when known** - Experimental design predicts cell numbers
5. **Don't over-filter** - Let downstream QC handle additional filtering
6. **Chain with QC pipeline** - Empty droplet filtering is step 1; always follow with quality control
7. **Add consistent QC metrics** - Use uppercase gene name patterns for cross-species compatibility

## Integration with Other Skills

### When to Use This Skill

Use `kallisto-bustools-import` **only** when:
- Data was processed with kallisto/bustools (kb-python) pipeline
- Empty droplets have NOT been filtered yet (using `counts_unfiltered/` output)

### Standard Workflows

**For kallisto/bustools data with unfiltered empty droplets:**
1. **kallisto-bustools-import** (this skill) - Load data, filter empty droplets
2. **scanpy-basics** - Standard downstream analysis (QC, normalization, clustering, visualization)

**For CellRanger output or pre-filtered data:**
- Start directly with **scanpy-basics** - this skill is not needed

### When to Use Specialized Skills

- **single-cell-rna-qc** - Use when specific QC requests are made (MAD-based filtering, outlier detection)
- **single-cell-integration-ingest-bbknn** - Use when batch correction is explicitly requested
- **scvi-tools** - Use when deep learning-based integration is specifically requested

### Consistent QC Metrics Across Skills

All skills use the same gene annotation patterns for cross-species compatibility:

```python
# Standard pattern used in scanpy-basics and this skill
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.upper().str.contains(r"^HB(?!P)")
```

### Passing Data Between Skills

```python
# After kb import and filtering
adata.write('adata_kb_filtered.h5ad')

# Continue with scanpy-basics for standard analysis
# Or use specialized skills when specific requests are made
```

## Troubleshooting

### Common Issues

**No clear knee in plot:**

- Dataset may have low cell recovery
- Try specifying `--expected-cells` based on experimental design
- Consider using `--lower 10` to include more points in detection
- For difficult datasets, consider statistical methods (EmptyDrops)

**Too few cells after filtering:**

- Use `--umi-threshold` to manually set a lower threshold
- Check if sample had low cell viability
- Verify you're using `counts_unfiltered/` not `counts_filtered/`

**Memory errors with large datasets:**

- kb-python output is typically sparse; ensure sparse matrix operations
- Process samples individually before concatenating
- Use `--skip-pca` to reduce memory usage

**File not found errors:**

- Verify kb-python completed successfully
- Check for `counts_unfiltered/adata.h5ad` vs `counts_filtered/adata.h5ad`
- For MTX files, ensure all three files exist (.mtx, .barcodes.txt, .genes.txt)

**Threshold detection fails:**

- May have insufficient unique points above `lower` threshold
- Try lowering the `lower` parameter
- Manually specify threshold with `--umi-threshold`

## Reference Materials

For detailed methodology on knee plot analysis and empty droplet detection:

- Macosko et al., "Highly parallel genome-wide expression profiling of individual cells using nanoliter droplets" (2015) - Original Drop-seq paper introducing knee plots
- Lun et al., "EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data" (2019) - Statistical approach to empty droplet detection
- DropletUtils R package documentation - https://bioconductor.org/packages/DropletUtils

## Next Steps After Import

Standard downstream analysis (use `scanpy-basics` skill):

- Quality control, normalization, and log-transformation
- Doublet detection with Scrublet
- Feature selection and dimensionality reduction
- Clustering and cell type annotation

**Specialized requests:**

- **Specific QC requirements** (MAD-based filtering, detailed outlier detection) → use `single-cell-rna-qc` skill
- **Batch correction** for combining multiple datasets → use `single-cell-integration-ingest-bbknn` skill
- **Deep learning integration** (scVI, scANVI) → use `scvi-tools` skill
