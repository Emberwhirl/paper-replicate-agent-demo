#!/usr/bin/env python3
"""
Symmetric integration using BBKNN (Batch Balanced K-Nearest Neighbors).

This is a convenience script that runs a complete BBKNN workflow using the
modular functions from integration_core.py.

Usage:
    python integrate_bbknn.py combined_data.h5ad --output integrated_bbknn.h5ad
"""

import anndata as ad
import scanpy as sc
import sys
import os
import argparse

# Import our modular utilities
from integration_core import (
    integrate_bbknn,
    print_integration_summary
)

print("=" * 80)
print("Single-Cell Integration: BBKNN (Batch Balanced K-Nearest Neighbors)")
print("=" * 80)

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description='Symmetric single-cell data integration using BBKNN',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
Examples:
  python integrate_bbknn.py combined.h5ad
  python integrate_bbknn.py combined.h5ad --output integrated_bbknn.h5ad
  python integrate_bbknn.py combined.h5ad --batch-key sample
  python integrate_bbknn.py combined.h5ad --neighbors-within-batch 5
    """
)

parser.add_argument('input', help='Combined h5ad file with batch column')
parser.add_argument('--output', '-o', default='integrated_bbknn.h5ad',
                    help='Output file path (default: integrated_bbknn.h5ad)')
parser.add_argument('--batch-key', default='batch',
                    help='Observation column indicating batch (default: batch)')
parser.add_argument('--neighbors-within-batch', type=int, default=3,
                    help='Number of neighbors from each batch (default: 3)')
parser.add_argument('--celltype-key', default=None,
                    help='Cell type column for visualization (optional)')
parser.add_argument('--output-dir', default=None,
                    help='Output directory for plots (default: same as output file)')
parser.add_argument('--skip-plots', action='store_true',
                    help='Skip generating visualization plots')

args = parser.parse_args()

# Verify input file exists
if not os.path.exists(args.input):
    print(f"\nError: Input file '{args.input}' not found!")
    sys.exit(1)

# Set up output directory
if args.output_dir:
    output_dir = args.output_dir
else:
    output_dir = os.path.dirname(args.output) or '.'

os.makedirs(output_dir, exist_ok=True)

# Display parameters
print(f"\nParameters:")
print(f"  Input: {args.input}")
print(f"  Batch key: {args.batch_key}")
print(f"  Neighbors within batch: {args.neighbors_within_batch}")
print(f"  Output: {args.output}")

# Load data
print("\n[1/3] Loading data...")
adata = ad.read_h5ad(args.input)
print(f"  Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

# Check for batch key
if args.batch_key not in adata.obs:
    print(f"\nError: Batch key '{args.batch_key}' not found in dataset.")
    print(f"Available columns: {list(adata.obs.columns)}")
    sys.exit(1)

n_batches = adata.obs[args.batch_key].nunique()
batch_counts = adata.obs[args.batch_key].value_counts()
print(f"  Found {n_batches} batches:")
for batch, count in batch_counts.items():
    print(f"    {batch}: {count} cells")

# Detect cell type key if not specified
celltype_key = args.celltype_key
if celltype_key is None:
    for key in ['celltype', 'cell_type', 'CellType', 'cluster', 'leiden']:
        if key in adata.obs:
            celltype_key = key
            break

if celltype_key:
    print(f"  Cell type column: {celltype_key} ({adata.obs[celltype_key].nunique()} types)")

# Apply BBKNN integration
print("\n[2/3] Applying BBKNN integration...")
integrate_bbknn(adata, batch_key=args.batch_key,
                neighbors_within_batch=args.neighbors_within_batch)

# Generate plots
if not args.skip_plots:
    print("\n[3/3] Generating visualizations...")

    sc.settings.figdir = output_dir

    # Prepare color columns
    color_cols = [args.batch_key]
    if celltype_key:
        color_cols.append(celltype_key)

    # UMAP colored by batch and cell type
    sc.pl.umap(adata, color=color_cols, wspace=0.5,
               save='_bbknn.png', show=False)
    print(f"  Saved: {os.path.join(output_dir, 'umap_bbknn.png')}")

    # Batch distribution plot
    if len(batch_counts) <= 6:
        # Side-by-side for small number of batches
        sc.pl.umap(adata, color=args.batch_key, groups=list(batch_counts.index),
                   ncols=3, save='_batches.png', show=False)
        print(f"  Saved: {os.path.join(output_dir, 'umap_batches.png')}")
else:
    print("\n[3/3] Skipping visualizations (--skip-plots)")

# Save results
print("\nSaving integrated data...")
adata.write(args.output)
print(f"  Saved: {args.output}")

# Print summary
print("\n" + "=" * 80)
print("Integration Summary")
print("=" * 80)

print_integration_summary(adata)

print("\n" + "=" * 80)
print("BBKNN Integration Complete!")
print("=" * 80)
print(f"\nOutput file: {args.output}")
print("\nNotes:")
print("  - BBKNN modifies the neighbors graph, not the expression matrix")
print("  - UMAP was computed on the batch-corrected graph")
print("  - Original data matrix is unchanged")
print("\nNext steps:")
print("  - Cluster with sc.tl.leiden(adata)")
print("  - Visualize with sc.pl.umap(adata, color=['batch', 'celltype'])")
print("  - Run differential expression between conditions")
