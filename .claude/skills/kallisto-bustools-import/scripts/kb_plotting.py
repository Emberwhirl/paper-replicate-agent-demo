"""
Visualization functions for kallisto/bustools import workflow.

Functions:
- plot_pca_embedding: 2D PCA projection of barcodes
- plot_saturation: Library saturation scatter plot
- plot_knee: Knee plot with knee and inflection threshold lines
- plot_knee_from_ranks: Knee plot using pre-computed barcode ranks DataFrame
"""

from pathlib import Path
from typing import Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse


def plot_pca_embedding(
    adata,
    output_path: Union[str, Path],
    title: str = 'PCA Embedding of Barcodes'
) -> None:
    """
    Generate 2D PCA projection of all barcodes.
    
    Parameters
    ----------
    adata : AnnData
        Count matrix from kb-python output
    output_path : str or Path
        Path to save the plot
    title : str
        Plot title
    """
    from sklearn.decomposition import TruncatedSVD
    
    # Perform SVD (works with sparse matrices)
    tsvd = TruncatedSVD(n_components=2, random_state=42)
    X_pca = tsvd.fit_transform(adata.X)
    
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.scatter(X_pca[:, 0], X_pca[:, 1], alpha=0.5, c='steelblue', s=1)
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_title(title)
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved PCA embedding to {output_path}")


def plot_saturation(
    adata,
    output_path: Union[str, Path],
    title: str = 'Library Saturation'
) -> None:
    """
    Generate library saturation plot (genes detected vs UMI counts).
    
    A saturated library shows a plateau where more UMIs don't yield more genes.
    An unsaturated library shows continued linear relationship.
    
    Parameters
    ----------
    adata : AnnData
        Count matrix from kb-python output
    output_path : str or Path
        Path to save the plot
    title : str
        Plot title
    """
    if sparse.issparse(adata.X):
        umi_counts = np.asarray(adata.X.sum(axis=1)).flatten()
        genes_detected = np.asarray((adata.X > 0).sum(axis=1)).flatten()
    else:
        umi_counts = np.asarray(adata.X.sum(axis=1)).flatten()
        genes_detected = np.asarray((adata.X > 0).sum(axis=1)).flatten()
    
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.scatter(umi_counts, genes_detected, color='steelblue', alpha=0.25, s=3)
    ax.set_xlabel('UMI Counts')
    ax.set_ylabel('Genes Detected')
    ax.set_title(title)
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    # Set reasonable axis limits
    ax.set_xlim((0.5, umi_counts.max() * 1.5))
    ax.set_ylim((0.5, genes_detected.max() * 1.5))
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved saturation plot to {output_path}")


def plot_knee(
    adata,
    knee_threshold: float = None,
    inflection_threshold: float = None,
    output_path: Union[str, Path] = None,
    title: str = 'Knee Plot',
    df_ranks: pd.DataFrame = None
) -> None:
    """
    Generate knee plot for empty droplet threshold visualization.
    
    Shows barcode rank vs total UMIs with knee and/or inflection thresholds.
    Uses the DropletUtils style with horizontal lines at rank cutoffs and
    vertical lines at UMI thresholds.
    
    Parameters
    ----------
    adata : AnnData
        Count matrix from kb-python output
    knee_threshold : float, optional
        Knee UMI threshold to display (blue line)
    inflection_threshold : float, optional
        Inflection UMI threshold to display (green line)
    output_path : str or Path, optional
        Path to save the plot. If None, displays interactively.
    title : str
        Plot title
    df_ranks : pd.DataFrame, optional
        Pre-computed barcode ranks DataFrame with 'rank' and 'total' columns.
        If None, computes from adata.
    """
    # Get rank and total data
    if df_ranks is not None:
        df_unique = df_ranks[["rank", "total"]].drop_duplicates()
        df_unique = df_unique[df_unique["total"] > 0].sort_values(
            "total", ascending=False
        )
        ranks = df_unique["rank"].values
        totals = df_unique["total"].values
        df_for_cutoff = df_ranks
    else:
        # Compute from adata
        if sparse.issparse(adata.X):
            umi_counts = np.asarray(adata.X.sum(axis=1)).flatten()
        else:
            umi_counts = np.asarray(adata.X.sum(axis=1)).flatten()
        
        sorted_counts = np.sort(umi_counts)[::-1]
        computed_ranks = np.arange(1, len(sorted_counts) + 1)
        
        # Filter for unique values
        unique_mask = np.concatenate(
            ([True], sorted_counts[1:] != sorted_counts[:-1])
        )
        totals = sorted_counts[unique_mask]
        ranks = computed_ranks[unique_mask]
        
        # Create DataFrame for cutoff calculation
        df_for_cutoff = pd.DataFrame({
            "rank": computed_ranks,
            "total": sorted_counts
        })
    
    fig, ax = plt.subplots(figsize=(9, 6))
    
    # Plot the curve (x=total, y=rank for log-log)
    ax.plot(totals, ranks, "k-", linewidth=1.5)
    
    # Add inflection threshold (green, more permissive)
    if inflection_threshold is not None:
        # Calculate rank cutoff
        inflection_rank_cutoff = df_for_cutoff.loc[
            df_for_cutoff["total"] > inflection_threshold, "rank"
        ].max()
        
        ax.axhline(
            y=inflection_rank_cutoff, color="forestgreen",
            linestyle="--", linewidth=1.5
        )
        ax.axvline(
            x=inflection_threshold, color="forestgreen",
            linestyle="--", linewidth=1.5,
            label=f"Inflection: {int(inflection_threshold):,} UMIs"
        )
    
    # Add knee threshold (blue, more conservative)
    if knee_threshold is not None:
        knee_rank_cutoff = df_for_cutoff.loc[
            df_for_cutoff["total"] > knee_threshold, "rank"
        ].max()
        
        ax.axhline(
            y=knee_rank_cutoff, color="dodgerblue",
            linestyle="--", linewidth=1.5
        )
        ax.axvline(
            x=knee_threshold, color="dodgerblue",
            linestyle="--", linewidth=1.5,
            label=f"Knee: {int(knee_threshold):,} UMIs"
        )
    
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Total UMIs")
    ax.set_ylabel("Barcode Rank")
    ax.set_title(title)
    
    if knee_threshold is not None or inflection_threshold is not None:
        ax.legend(loc="upper right")
    
    ax.minorticks_on()
    ax.grid(True, which='major', alpha=0.3)
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved knee plot to {output_path}")
    else:
        plt.show()


def plot_knee_interactive(adata, initial_threshold: int = 200):
    """
    Interactive knee plot for Jupyter notebooks with threshold slider.
    
    Parameters
    ----------
    adata : AnnData
        Count matrix from kb-python output
    initial_threshold : int
        Initial threshold value for slider
    
    Returns
    -------
    Widget
        Interactive widget (requires ipywidgets)
    """
    try:
        from ipywidgets import interact, IntSlider
    except ImportError:
        print("ipywidgets not available. Use plot_knee() for static plots.")
        return None
    
    if sparse.issparse(adata.X):
        umi_counts = np.asarray(adata.X.sum(axis=1)).flatten()
    else:
        umi_counts = np.asarray(adata.X.sum(axis=1)).flatten()
    
    sorted_counts = np.sort(umi_counts)[::-1]
    ranks = np.arange(1, len(sorted_counts) + 1)
    max_threshold = int(np.percentile(sorted_counts, 99))
    
    def update_plot(threshold):
        fig, ax = plt.subplots(figsize=(10, 7))
        ax.loglog(sorted_counts, ranks, linewidth=2, color='steelblue')
        ax.axvline(x=threshold, linewidth=2, color='forestgreen', linestyle='--')
        
        num_cells = np.sum(sorted_counts > threshold)
        ax.axhline(y=num_cells, linewidth=2, color='forestgreen', linestyle=':')
        
        ax.set_xlabel('Total UMIs')
        ax.set_ylabel('Barcode Rank')
        ax.set_title(
            f'Knee Plot | Threshold: {threshold:,} UMIs | Cells: {num_cells:,}'
        )
        ax.grid(True, which='both', alpha=0.3)
        plt.show()
    
    slider = IntSlider(
        value=initial_threshold,
        min=10,
        max=max_threshold,
        step=10,
        description='Threshold:',
        continuous_update=False
    )
    
    return interact(update_plot, threshold=slider)