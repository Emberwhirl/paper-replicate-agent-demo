#!/usr/bin/env python3
"""
Quality Control Analysis Script for Scanpy

Performs comprehensive quality control on single-cell RNA-seq data,
including calculating metrics, generating QC plots, doublet detection, and filtering cells.

Usage:
    python qc_analysis.py <input_file> [--output <output_file>]
"""

import argparse
import scanpy as sc
import matplotlib.pyplot as plt


def calculate_qc_metrics(adata):
    """
    Calculate QC metrics for mitochondrial, ribosomal, and hemoglobin genes.

    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix

    Returns:
    --------
    AnnData
        Annotated data matrix with QC metrics calculated
    """
    # Normalize case to ensure cross-species compatibility and robustness

    ## mitochondrial genes
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    ## ribosomal genes
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    ## hemoglobin genes
    adata.var["hb"] = adata.var_names.str.upper().str.contains(r"^HB(?!P)")

    # Calculate QC metrics with log1p transformation
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=['mt', 'ribo', 'hb'], inplace=True, log1p=True
    )

    print("\n=== QC Metrics Summary ===")
    print(f"Total cells: {adata.n_obs}")
    print(f"Total genes: {adata.n_vars}")
    print(f"Mean genes per cell: {adata.obs['n_genes_by_counts'].mean():.2f}")
    print(f"Mean counts per cell: {adata.obs['total_counts'].mean():.2f}")
    print(f"Mean mitochondrial %: {adata.obs['pct_counts_mt'].mean():.2f}")
    print(f"Mean ribosomal %: {adata.obs['pct_counts_ribo'].mean():.2f}")
    print(f"Mean hemoglobin %: {adata.obs['pct_counts_hb'].mean():.2f}")

    return adata


def generate_qc_plots(adata, output_prefix='qc'):
    """
    Generate comprehensive QC plots.

    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    output_prefix : str
        Prefix for saved figure files
    """
    # Create figure directory if it doesn't exist
    import os
    os.makedirs('figures', exist_ok=True)

    # Violin plots for QC metrics
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True, save=f'_{output_prefix}_violin.pdf')

    # Scatter plots
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt',
                  save=f'_{output_prefix}_mt_scatter.pdf')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',
                  color='pct_counts_mt', save=f'_{output_prefix}_genes_scatter.pdf')

    # Highest expressing genes
    sc.pl.highest_expr_genes(adata, n_top=20,
                              save=f'_{output_prefix}_highest_expr.pdf')

    print(f"\nQC plots saved to figures/ directory with prefix '{output_prefix}'")


def filter_data(adata, min_genes=100, max_genes=None,
                min_counts=None, max_counts=None, min_cells=3):
    """
    Filter cells and genes based on QC thresholds.
    Uses permissive filtering; QC thresholds should be revisited after clustering.

    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    min_genes : int
        Minimum number of genes per cell (default: 100)
    max_genes : int, optional
        Maximum number of genes per cell
    min_counts : int, optional
        Minimum number of counts per cell
    max_counts : int, optional
        Maximum number of counts per cell
    min_cells : int
        Minimum number of cells per gene (default: 3)

    Returns:
    --------
    AnnData
        Filtered annotated data matrix
    """
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars

    # Filter cells
    sc.pp.filter_cells(adata, min_genes=min_genes)
    if max_genes:
        adata = adata[adata.obs['n_genes_by_counts'] < max_genes, :]
    if min_counts:
        adata = adata[adata.obs['total_counts'] >= min_counts, :]
    if max_counts:
        adata = adata[adata.obs['total_counts'] < max_counts, :]

    # Filter genes
    sc.pp.filter_genes(adata, min_cells=min_cells)

    print(f"\n=== Filtering Results ===")
    print(f"Cells: {n_cells_before} -> {adata.n_obs} ({adata.n_obs/n_cells_before*100:.1f}% retained)")
    print(f"Genes: {n_genes_before} -> {adata.n_vars} ({adata.n_vars/n_genes_before*100:.1f}% retained)")

    return adata


def run_doublet_detection(adata, batch_key=None):
    """
    Run Scrublet doublet detection.

    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    batch_key : str, optional
        Column in adata.obs for batch/sample information

    Returns:
    --------
    AnnData
        Annotated data matrix with doublet scores
    """
    print("\n=== Doublet Detection ===")

    if batch_key and batch_key in adata.obs.columns:
        sc.pp.scrublet(adata, batch_key=batch_key)
        print(f"Ran Scrublet with batch_key='{batch_key}'")
    else:
        sc.pp.scrublet(adata)
        print("Ran Scrublet on full dataset")

    n_doublets = adata.obs['predicted_doublet'].sum()
    print(f"Predicted doublets: {n_doublets} ({n_doublets/adata.n_obs*100:.1f}%)")

    return adata


def main():
    parser = argparse.ArgumentParser(description='QC analysis for single-cell data')
    parser.add_argument('input', help='Input file (h5ad, 10X mtx, csv, etc.)')
    parser.add_argument('--output', default='qc_filtered.h5ad',
                        help='Output file name (default: qc_filtered.h5ad)')
    parser.add_argument('--min-genes', type=int, default=100,
                        help='Min genes per cell (default: 100)')
    parser.add_argument('--min-cells', type=int, default=3,
                        help='Min cells per gene (default: 3)')
    parser.add_argument('--batch-key', type=str, default=None,
                        help='Column name for batch/sample info (for doublet detection)')
    parser.add_argument('--skip-plots', action='store_true',
                        help='Skip generating QC plots')
    parser.add_argument('--skip-doublets', action='store_true',
                        help='Skip doublet detection')

    args = parser.parse_args()

    # Configure scanpy
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=100, facecolor='white')
    sc.settings.figdir = './figures/'

    print(f"Loading data from: {args.input}")

    # Load data based on file extension
    if args.input.endswith('.h5ad'):
        adata = sc.read_h5ad(args.input)
    elif args.input.endswith('.h5'):
        adata = sc.read_10x_h5(args.input)
    elif args.input.endswith('.csv'):
        adata = sc.read_csv(args.input)
    else:
        # Try reading as 10X mtx directory
        adata = sc.read_10x_mtx(args.input)

    print(f"Loaded data: {adata.n_obs} cells x {adata.n_vars} genes")

    # Calculate QC metrics
    adata = calculate_qc_metrics(adata)

    # Generate QC plots (before filtering)
    if not args.skip_plots:
        print("\nGenerating QC plots (before filtering)...")
        generate_qc_plots(adata, output_prefix='qc_before')

    # Filter data
    adata = filter_data(adata, min_genes=args.min_genes, min_cells=args.min_cells)

    # Run doublet detection
    if not args.skip_doublets:
        adata = run_doublet_detection(adata, batch_key=args.batch_key)

    # Generate QC plots (after filtering)
    if not args.skip_plots:
        print("\nGenerating QC plots (after filtering)...")
        generate_qc_plots(adata, output_prefix='qc_after')

    # Save filtered data
    print(f"\nSaving filtered data to: {args.output}")
    adata.write_h5ad(args.output)

    print("\n=== QC Analysis Complete ===")
    print("Note: Consider revisiting QC thresholds after clustering to identify")
    print("clusters driven by poor quality cells or doublets.")


if __name__ == "__main__":
    main()
