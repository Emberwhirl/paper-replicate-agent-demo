#!/usr/bin/env python3
"""
Reference-based integration using scanpy's ingest.

This is a convenience script that runs a complete ingest workflow using the
modular functions from integration_core.py.

Usage:
    python integrate_ingest.py reference.h5ad query.h5ad --output integrated.h5ad
"""

import anndata as ad
import scanpy as sc
import sys
import os
import argparse

# Import our modular utilities
from integration_core import (
    intersect_genes,
    train_reference,
    ingest_query,
    combine_datasets,
    print_integration_summary
)

print("=" * 80)
print("Single-Cell Integration: Ingest (Reference-Based)")
print("=" * 80)

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description='Reference-based single-cell data integration using ingest',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
Examples:
  python integrate_ingest.py reference.h5ad query.h5ad
  python integrate_ingest.py reference.h5ad query.h5ad --output integrated.h5ad
  python integrate_ingest.py reference.h5ad query.h5ad --label-key cell_type
  python integrate_ingest.py reference.h5ad query.h5ad --output-dir results/
    """
)

parser.add_argument('reference', help='Reference h5ad file with annotations')
parser.add_argument('query', help='Query h5ad file to map onto reference')
parser.add_argument('--output', '-o', default='integrated.h5ad',
                    help='Output file path (default: integrated.h5ad)')
parser.add_argument('--label-key', default='celltype',
                    help='Observation column to transfer (default: celltype)')
parser.add_argument('--output-dir', default=None,
                    help='Output directory for plots (default: same as output file)')
parser.add_argument('--skip-plots', action='store_true',
                    help='Skip generating visualization plots')
parser.add_argument('--n-pcs', type=int, default=50,
                    help='Number of principal components (default: 50)')
parser.add_argument('--n-neighbors', type=int, default=15,
                    help='Number of neighbors for graph (default: 15)')

args = parser.parse_args()

# Verify input files exist
for path, name in [(args.reference, 'Reference'), (args.query, 'Query')]:
    if not os.path.exists(path):
        print(f"\nError: {name} file '{path}' not found!")
        sys.exit(1)

# Set up output directory
if args.output_dir:
    output_dir = args.output_dir
else:
    output_dir = os.path.dirname(args.output) or '.'

os.makedirs(output_dir, exist_ok=True)

# Display parameters
print(f"\nParameters:")
print(f"  Reference: {args.reference}")
print(f"  Query: {args.query}")
print(f"  Label key: {args.label_key}")
print(f"  PCs: {args.n_pcs}, Neighbors: {args.n_neighbors}")
print(f"  Output: {args.output}")

# Load data
print("\n[1/5] Loading data...")
adata_ref = ad.read_h5ad(args.reference)
adata_query = ad.read_h5ad(args.query)

print(f"  Reference: {adata_ref.n_obs} cells × {adata_ref.n_vars} genes")
print(f"  Query: {adata_query.n_obs} cells × {adata_query.n_vars} genes")

# Check for label key in reference
if args.label_key not in adata_ref.obs:
    print(f"\nError: Label key '{args.label_key}' not found in reference dataset.")
    print(f"Available columns: {list(adata_ref.obs.columns)}")
    sys.exit(1)

n_celltypes = adata_ref.obs[args.label_key].nunique()
print(f"  Reference cell types: {n_celltypes}")

# Intersect genes
print("\n[2/5] Finding common gene space...")
adata_ref, adata_query = intersect_genes(adata_ref, adata_query)

# Train on reference
print("\n[3/5] Training reference model...")
train_reference(adata_ref, n_pcs=args.n_pcs, n_neighbors=args.n_neighbors)

# Ingest query
print("\n[4/5] Mapping query onto reference...")
ingest_query(adata_query, adata_ref, label_key=args.label_key)

# Combine datasets
print("\n[5/5] Combining datasets...")
adata_concat = combine_datasets(adata_ref, adata_query, label_key=args.label_key)

# Generate plots
if not args.skip_plots:
    print("\nGenerating visualizations...")

    sc.settings.figdir = output_dir

    # Reference UMAP
    sc.pl.umap(adata_ref, color=args.label_key, title='Reference Dataset',
               save='_reference.png', show=False)
    print(f"  Saved: {os.path.join(output_dir, 'umap_reference.png')}")

    # Query UMAP
    sc.pl.umap(adata_query, color=args.label_key, title='Query (Mapped)',
               save='_query.png', show=False)
    print(f"  Saved: {os.path.join(output_dir, 'umap_query.png')}")

    # Combined UMAP
    sc.pl.umap(adata_concat, color=['batch', args.label_key],
               title=['Batch', 'Cell Type'], wspace=0.5,
               save='_combined.png', show=False)
    print(f"  Saved: {os.path.join(output_dir, 'umap_combined.png')}")

# Save results
print("\nSaving integrated data...")
adata_concat.write(args.output)
print(f"  Saved: {args.output}")

# Print summary
print("\n" + "=" * 80)
print("Integration Summary")
print("=" * 80)

print_integration_summary(adata_concat)

# Show label transfer statistics
orig_key = f'{args.label_key}_orig'
if orig_key in adata_query.obs:
    from integration_core import evaluate_label_transfer
    print("\nLabel Transfer Evaluation (query cells only):")
    confusion = evaluate_label_transfer(adata_query, args.label_key, orig_key)

print("\n" + "=" * 80)
print("Integration Complete!")
print("=" * 80)
print(f"\nOutput file: {args.output}")
print("\nNext steps:")
print("  - Visualize with sc.pl.umap(adata, color=['batch', 'celltype'])")
print("  - Cluster with sc.tl.leiden(adata)")
print("  - Run differential expression between conditions")
