#!/usr/bin/env python3
"""
Complete kallisto/bustools import and empty droplet filtering pipeline.

Usage:
    python kb_import.py counts_unfiltered/adata.h5ad
    python kb_import.py input.h5ad --output-dir results/ --umi-threshold 500
    python kb_import.py input.h5ad --expected-cells 5000
    python kb_import.py counts_unfiltered/ --from-mtx  # Read from MTX files

This script:
1. Loads kb-python/kallisto/bustools output (h5ad or MTX format)
2. Generates PCA embedding visualization
3. Creates library saturation plot
4. Generates knee plot using DropletUtils algorithm and detects thresholds
5. Filters empty droplets using inflection threshold (more permissive)
6. Saves filtered data and summary statistics

Third-party dependencies & licensing
-------------------------------------
This skill *orchestrates* external open-source tools that are NOT redistributed
here -- install them yourself. They remain under their own upstream licenses:
  - kb-python (kallisto | bustools) -- BSD-2-Clause
  - DropletUtils                    -- GPL-3.0
  - anndata                         -- BSD-3-Clause
The skill code itself is MIT-licensed; see the LICENSE file in the skill root.
"""

import argparse
from pathlib import Path
import sys

# Add scripts directory to path for imports
script_dir = Path(__file__).parent
sys.path.insert(0, str(script_dir))

from kb_core import (
    load_kb_output,
    read_mtx_output,
    calculate_barcode_ranks,
    detect_knee_threshold,
    filter_empty_droplets,
    print_import_summary
)
from kb_plotting import (
    plot_pca_embedding,
    plot_saturation,
    plot_knee
)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Import kallisto/bustools output and filter empty droplets',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        'input',
        type=str,
        help='Path to adata.h5ad or directory containing MTX files'
    )
    parser.add_argument(
        '--from-mtx',
        action='store_true',
        help='Read from MTX files instead of h5ad'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default=None,
        help='Output directory (default: <input_stem>_kb_import_results/)'
    )
    parser.add_argument(
        '--min-genes',
        type=int,
        default=None,
        help='Minimum genes per cell for filtering (default: no filtering, '
             'leave for downstream QC)'
    )
    parser.add_argument(
        '--umi-threshold',
        type=float,
        default=None,
        help='Manual UMI threshold (auto-detected from knee if not specified)'
    )
    parser.add_argument(
        '--expected-cells',
        type=int,
        default=None,
        help='Expected number of cells (helps refine threshold detection)'
    )
    parser.add_argument(
        '--use-knee',
        action='store_true',
        help='Use knee threshold instead of inflection (more conservative)'
    )
    parser.add_argument(
        '--lower',
        type=float,
        default=100,
        help='Lower bound on UMI count for threshold detection'
    )
    parser.add_argument(
        '--skip-saturation',
        action='store_true',
        help='Skip library saturation plot'
    )
    parser.add_argument(
        '--skip-pca',
        action='store_true',
        help='Skip initial PCA visualization'
    )
    return parser.parse_args()


def main():
    args = parse_args()
    
    # Setup paths
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input path not found: {input_path}")
        sys.exit(1)
    
    # Determine output directory name
    if input_path.is_file():
        stem = input_path.stem.replace('.h5ad', '').replace(
            'adata', input_path.parent.name or 'data'
        )
    else:
        stem = input_path.name
    
    output_dir = (
        Path(args.output_dir) if args.output_dir
        else Path(f'{stem}_kb_import_results')
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n{'=' * 60}")
    print("Kallisto/Bustools Import Pipeline")
    print('=' * 60)
    print(f"Input:  {input_path}")
    print(f"Output: {output_dir}/")
    print('=' * 60 + '\n')
    
    # Step 1: Load data
    print("Step 1: Loading kb-python output...")
    if args.from_mtx or input_path.is_dir():
        adata = read_mtx_output(input_path)
    else:
        adata = load_kb_output(str(input_path))
    print_import_summary(adata, label='Before filtering')
    
    # Step 2: PCA visualization
    if not args.skip_pca:
        print("\nStep 2: Generating PCA embedding...")
        plot_pca_embedding(adata, str(output_dir / 'pca_embedding.png'))
    
    # Step 3: Library saturation
    if not args.skip_saturation:
        print("\nStep 3: Generating library saturation plot...")
        plot_saturation(adata, str(output_dir / 'library_saturation.png'))
    
    # Step 4: Knee plot and threshold detection using DropletUtils algorithm
    print("\nStep 4: Calculating barcode ranks and detecting thresholds...")
    
    if args.umi_threshold is not None:
        # Manual threshold
        threshold = args.umi_threshold
        knee_threshold = threshold
        inflection_threshold = threshold
        print(f"Using manual threshold: {threshold:.0f} UMIs")
    elif args.expected_cells is not None:
        # Use expected cells
        threshold, _ = detect_knee_threshold(
            adata,
            expected_cells=args.expected_cells,
            use_inflection=not args.use_knee
        )
        knee_threshold = threshold
        inflection_threshold = threshold
    else:
        # Auto-detect using DropletUtils algorithm
        df_ranks, knee_threshold, inflection_threshold = calculate_barcode_ranks(
            adata,
            lower=args.lower,
        )
        print(f"Knee threshold: {knee_threshold:.0f} UMIs")
        print(f"Inflection threshold: {inflection_threshold:.0f} UMIs")
        
        if args.use_knee:
            threshold = knee_threshold
            print(f"Using knee threshold (more conservative): {threshold:.0f}")
        else:
            threshold = inflection_threshold
            print(f"Using inflection threshold (more permissive): {threshold:.0f}")
    
    # Generate knee plot with both thresholds
    plot_knee(
        adata,
        knee_threshold=knee_threshold,
        inflection_threshold=inflection_threshold,
        output_path=str(output_dir / 'knee_plot.png'),
        title='Knee Plot with Thresholds'
    )
    
    # Step 5: Filter empty droplets
    print("\nStep 5: Filtering empty droplets...")
    adata_filtered = filter_empty_droplets(
        adata,
        umi_threshold=threshold,
        min_genes=args.min_genes,
        filter_genes_min_cells=1
    )
    print_import_summary(adata_filtered, label='After filtering')
    
    # Step 6: Save results
    print("\nStep 6: Saving results...")
    output_h5ad = output_dir / f'{stem}_filtered.h5ad'
    adata_filtered.write(str(output_h5ad))
    print(f"Saved filtered data to {output_h5ad}")
    
    # Save summary
    summary_path = output_dir / 'filtering_summary.txt'
    with open(summary_path, 'w') as f:
        f.write("Kallisto/Bustools Import Summary\n")
        f.write("=" * 40 + "\n\n")
        f.write(f"Input file: {input_path}\n")
        f.write(f"Knee threshold: {knee_threshold:.0f} UMIs\n")
        f.write(f"Inflection threshold: {inflection_threshold:.0f} UMIs\n")
        f.write(f"Applied threshold: {threshold:.0f} UMIs\n")
        if args.min_genes:
            f.write(f"Min genes filter: {args.min_genes}\n")
        f.write(f"\nBarcodes before: {adata.n_obs:,}\n")
        f.write(f"Barcodes after:  {adata_filtered.n_obs:,}\n")
        f.write(f"Barcodes removed: {adata.n_obs - adata_filtered.n_obs:,}\n")
        f.write(f"Retention rate: {100 * adata_filtered.n_obs / adata.n_obs:.1f}%\n")
        f.write(f"\nGenes before: {adata.n_vars:,}\n")
        f.write(f"Genes after:  {adata_filtered.n_vars:,}\n")
    print(f"Saved summary to {summary_path}")
    
    print(f"\n{'=' * 60}")
    print("Pipeline complete!")
    print(f"Results saved to: {output_dir}/")
    print('=' * 60 + '\n')


if __name__ == '__main__':
    main()