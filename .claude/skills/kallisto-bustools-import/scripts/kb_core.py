"""
Core utility functions for kallisto/bustools import and empty droplet filtering.

This module implements the DropletUtils barcodeRanks algorithm for accurate
knee/inflection point detection, providing methodologically correct empty
droplet filtering for single-cell RNA-seq data.

Functions:
- load_kb_output: Load adata.h5ad from kb-python output
- read_mtx_output: Read kallisto/bustools MTX format output
- calculate_barcode_ranks: Calculate barcode ranks with knee/inflection detection
- detect_knee_threshold: Auto-detect UMI threshold using knee point
- filter_empty_droplets: Apply threshold filtering to remove empty droplets
- calculate_saturation_metrics: Compute genes detected vs UMI counts
- print_import_summary: Print summary statistics

Reference:
    Lun et al., "EmptyDrops: distinguishing cells from empty droplets in
    droplet-based single-cell RNA sequencing data" (2019)
"""

from pathlib import Path
from typing import Dict, Tuple, Union

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import io, sparse


# ==============================================================================
# Data Loading Functions
# ==============================================================================


def load_kb_output(path: Union[str, Path]) -> ad.AnnData:
    """
    Load AnnData from kb-python/kallisto/bustools output.
    
    Parameters
    ----------
    path : str or Path
        Path to adata.h5ad file (typically in counts_unfiltered/ directory)
    
    Returns
    -------
    AnnData
        Loaded count matrix with barcodes as observations and genes as variables
    """
    adata = ad.read_h5ad(path)
    print(f"Loaded {adata.n_obs:,} barcodes x {adata.n_vars:,} genes")
    return adata


def read_mtx_output(
    count_dir: Union[str, Path],
    name: str = "cells_x_genes"
) -> ad.AnnData:
    """
    Read kallisto/bustools MTX format output.
    
    The MTX file from kb-python is in cells (barcodes) x genes format.
    Returns AnnData with cells as observations and genes as variables.
    
    Parameters
    ----------
    count_dir : str or Path
        Path to directory containing MTX files (e.g., counts_unfiltered/)
    name : str
        Base name of output files (default: 'cells_x_genes')
    
    Returns
    -------
    AnnData
        Count matrix with cells as observations and genes as variables
    """
    count_dir = Path(count_dir)
    
    # Read MTX file (kb-python outputs cells x genes, no transpose needed)
    mtx_path = count_dir / f"{name}.mtx"
    mat = io.mmread(mtx_path).tocsr()
    
    # Read barcodes (cells) - these correspond to matrix rows
    barcodes_file = count_dir / f"{name}.barcodes.txt"
    barcodes = barcodes_file.read_text().strip().split('\n')
    
    # Read gene names (prefer .genes.names.txt if available for gene symbols)
    # These correspond to matrix columns
    genes_names_file = count_dir / f"{name}.genes.names.txt"
    genes_file = count_dir / f"{name}.genes.txt"
    
    if genes_names_file.exists():
        genes = genes_names_file.read_text().strip().split('\n')
    else:
        genes = genes_file.read_text().strip().split('\n')
    
    # Verify dimensions match
    if mat.shape[0] != len(barcodes):
        raise ValueError(f"MTX rows ({mat.shape[0]}) != barcodes ({len(barcodes)})")
    if mat.shape[1] != len(genes):
        raise ValueError(f"MTX cols ({mat.shape[1]}) != genes ({len(genes)})")
    
    # Create AnnData object (obs=cells/barcodes, var=genes)
    adata = ad.AnnData(
        X=mat,
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=genes)
    )
    
    print(f"Loaded {adata.n_obs:,} barcodes x {adata.n_vars:,} genes from MTX")
    return adata


# ==============================================================================
# DropletUtils barcodeRanks Algorithm Implementation
# ==============================================================================


def _interpolate_simple(
    prop: np.ndarray,
    left_val: np.ndarray,
    right_val: np.ndarray
) -> np.ndarray:
    """Linear interpolation between two values."""
    return left_val + prop * (right_val - left_val)


def _find_interval_left_open(target: np.ndarray, vec: np.ndarray) -> np.ndarray:
    """
    Port of R's findInterval(target, vec, left.open=TRUE).
    
    Returns the index of the interval for each target value.
    With left.open=TRUE, intervals are (vec[i], vec[i+1]].
    Returns 0 if target <= vec[0], and len(vec) if target > vec[-1].
    """
    idx = np.searchsorted(vec, target, side='left')
    return idx - 1


def _interpolate_on_curve(
    target_distance: np.ndarray,
    cumulative_distance: np.ndarray,
    stepwise_distance: np.ndarray,
    x: np.ndarray,
    y: np.ndarray
) -> Dict[str, np.ndarray]:
    """
    Port of DropletUtils:::.interpolate_on_curve
    
    Finds the (x, y) coordinates on the curve at the specified distances
    along the curve. Uses left-open intervals like R's findInterval.
    """
    target_distance = np.asarray(target_distance, dtype=float)
    cumulative_distance = np.asarray(cumulative_distance, dtype=float)
    stepwise_distance = np.asarray(stepwise_distance, dtype=float)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    # Find interval indices
    idx = _find_interval_left_open(target_distance, cumulative_distance)
    
    # Clamp to valid range
    idx = np.clip(idx, 0, len(cumulative_distance) - 2)
    idx_p1 = idx + 1
    
    # Calculate interpolation proportion
    prop = (target_distance - cumulative_distance[idx]) / stepwise_distance[idx_p1]

    return {
        "x": _interpolate_simple(prop, x[idx], x[idx_p1]),
        "y": _interpolate_simple(prop, y[idx], y[idx_p1]),
    }


def _reorder(vals: np.ndarray, lens: np.ndarray, o: np.ndarray) -> np.ndarray:
    """
    Port of DropletUtils:::.reorder
    
    Creates a repeated array, then permutes it to undo the sort order.
    """
    out = np.repeat(vals, lens)
    result = np.empty_like(out)
    result[o] = out
    return result


def calculate_barcode_ranks(
    adata: ad.AnnData,
    lower: float = 100,
    exclude_from: int = 50,
    window: float = 1.0,
    gradient_threshold: float = -1.0,
) -> Tuple[pd.DataFrame, float, float]:
    """
    Calculate barcode ranks and detect knee/inflection points.
    
    This is an exact port of the DropletUtils::barcodeRanks algorithm,
    which uses a sliding window approach along the log-log curve to
    identify the knee (transition point) and inflection (steepest descent).
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with counts in adata.X (cells/barcodes x genes)
    lower : float
        Lower bound on total UMI count for knee/inflection identification.
        Default: 100
    exclude_from : int
        Number of highest ranking barcodes to exclude from identification.
        These are assumed to be true cells. Default: 50
    window : float
        Length of the window (in log10 units) for knee/inflection point
        identification along the curve. Default: 1.0
    gradient_threshold : float
        Maximum threshold on gradient for elbow point identification.
        Default: -1.0
    
    Returns
    -------
    df_ranks : pd.DataFrame
        DataFrame with 'rank' and 'total' columns for each barcode,
        indexed by barcode names.
    knee : float
        Total UMI count at the knee point (recommended for filtering).
    inflection : float
        Total UMI count at the inflection point (steepest descent,
        more permissive threshold).
    
    Notes
    -----
    The knee point represents the transition from real cells to empty
    droplets. The inflection point is where the gradient is steepest.
    For most applications, filtering at the inflection threshold is
    recommended as it is more permissive and retains more cells.
    
    Reference
    ---------
    DropletUtils R package: https://bioconductor.org/packages/DropletUtils
    """
    # Total counts per barcode
    if sparse.issparse(adata.X):
        totals = np.asarray(adata.X.sum(axis=1)).flatten().astype(float)
    else:
        totals = np.asarray(adata.X.sum(axis=1)).flatten().astype(float)

    # Order by decreasing totals
    o = np.argsort(-totals)
    totals_sorted = totals[o]

    # Run-length encoding of tied totals
    change = np.concatenate(([True], totals_sorted[1:] != totals_sorted[:-1]))
    run_starts = np.where(change)[0]
    run_lengths = np.diff(np.append(run_starts, len(totals_sorted)))
    run_totals = totals_sorted[run_starts]
    
    # Get mid-rank of each run
    run_rank = np.cumsum(run_lengths) - (run_lengths - 1) / 2

    # Filter for knee/inflection detection
    # Use strict inequality (>) per DropletUtils convention - 'lower' is a lower bound
    keep = run_totals > lower
    keep[run_rank <= exclude_from] = False
    
    if np.sum(keep) < 2:
        raise ValueError(
            "Insufficient unique points for computing knee/inflection points. "
            "Try lowering the 'lower' threshold or check data quality."
        )

    # Log10 transform of filtered values
    y = np.log10(run_totals[keep])
    x = np.log10(run_rank[keep])

    # Calculate cumulative distance along the curve
    dist_along_curve = np.concatenate(
        ([0], np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2))
    )
    cumdist = np.cumsum(dist_along_curve)
    rhs_loc = cumdist + window
    to_scan = rhs_loc <= cumdist[-1]

    if not np.any(to_scan):
        # Edge case: just pick the last point
        knee = inflection = 10 ** y[-1]
    else:
        left_x = x[to_scan]
        left_y = y[to_scan]

        right_info = _interpolate_on_curve(
            rhs_loc[to_scan], cumdist, dist_along_curve, x, y
        )
        right_x = right_info["x"]
        right_y = right_info["y"]

        # Distance in 2D space between the ends of the window
        window_gap = np.sqrt((left_x - right_x) ** 2 + (left_y - right_y) ** 2)
        
        # Gradient and intercept of the line between window ends
        window_gradient = (right_y - left_y) / (right_x - left_x)
        window_intercept = right_y - window_gradient * right_x

        # Find midpoint of each window
        mid_info = _interpolate_on_curve(
            cumdist[to_scan] + window / 2, cumdist, dist_along_curve, x, y
        )
        mid_x = mid_info["x"]
        mid_y = mid_info["y"]

        # Check if midpoint is above the end-connecting line
        mid_above = mid_y > window_gradient * mid_x + window_intercept

        # Find elbow points: midpoint below line AND gradient below threshold
        has_elbow_index = np.where(
            (~mid_above) & (window_gradient < gradient_threshold)
        )[0]
        
        if len(has_elbow_index) == 0:
            infl_window = np.argmin(window_gradient)
            maybe_knee = np.where(mid_above)[0]
        else:
            first_elbow_window = has_elbow_index[0]
            # Only consider windows up to and including the first elbow
            maybe_knee = np.where(mid_above[: first_elbow_window + 1])[0]

            # For inflection: find windows starting before the first elbow midpoint
            before_first_elbow = np.searchsorted(
                left_x, mid_x[first_elbow_window], side='right'
            )
            # Ensure at least 1 window is considered to avoid empty slice
            before_first_elbow = max(before_first_elbow, 1)
            infl_window = np.argmin(window_gradient[:before_first_elbow])

        # Find knee window
        if len(maybe_knee) == 0:
            # Fallback if curve is unusual
            knee_window = infl_window
        else:
            knee_window = maybe_knee[np.argmin(window_gap[maybe_knee])]

        # Pick actual knee/inflection points based on midpoint of window
        knee = 10 ** mid_y[knee_window]
        inflection = 10 ** mid_y[infl_window]

    # Construct output DataFrame
    rank_out = _reorder(run_rank, run_lengths, o)
    total_out = _reorder(run_totals, run_lengths, o)

    df = pd.DataFrame(
        {"rank": rank_out, "total": total_out},
        index=adata.obs_names
    )

    return df, knee, inflection


# ==============================================================================
# Threshold Detection and Filtering
# ==============================================================================


def detect_knee_threshold(
    adata: ad.AnnData,
    expected_cells: int = None,
    use_inflection: bool = True,
    lower: float = 100,
    exclude_from: int = 50,
    window: float = 1.0,
) -> Tuple[float, float]:
    """
    Detect UMI thresholds from knee plot using DropletUtils algorithm.
    
    Parameters
    ----------
    adata : AnnData
        Count matrix from kb-python output
    expected_cells : int, optional
        Expected number of cells. If provided, uses this to estimate
        threshold instead of automatic detection.
    use_inflection : bool
        If True (default), returns inflection threshold as the primary
        threshold (more permissive). If False, returns knee threshold.
    lower : float
        Lower bound on total UMI count. Default: 100
    exclude_from : int
        Number of top barcodes to exclude. Default: 50
    window : float
        Window length in log10 units. Default: 1.0
    
    Returns
    -------
    threshold : float
        Recommended UMI threshold for filtering
    other_threshold : float
        The alternative threshold (inflection if use_inflection=False,
        knee if use_inflection=True)
    """
    if expected_cells is not None:
        # Use expected cells to determine threshold
        if sparse.issparse(adata.X):
            umi_counts = np.asarray(adata.X.sum(axis=1)).flatten()
        else:
            umi_counts = np.asarray(adata.X.sum(axis=1)).flatten()
        
        sorted_counts = np.sort(umi_counts)[::-1]
        
        if expected_cells < len(sorted_counts):
            threshold = float(sorted_counts[expected_cells - 1])
            print(f"Using expected cells ({expected_cells:,}) threshold: {threshold:.0f}")
            return max(threshold, 1.0), threshold
    
    # Auto-detect using DropletUtils algorithm
    _, knee, inflection = calculate_barcode_ranks(
        adata,
        lower=lower,
        exclude_from=exclude_from,
        window=window,
    )
    
    print(f"Knee threshold: {knee:.0f} UMIs")
    print(f"Inflection threshold: {inflection:.0f} UMIs")
    
    if use_inflection:
        return inflection, knee
    else:
        return knee, inflection


def filter_empty_droplets(
    adata: ad.AnnData,
    umi_threshold: float,
    min_genes: int = None,
    filter_genes_min_cells: int = 1,
    inplace: bool = False
) -> ad.AnnData:
    """
    Filter empty droplets based on UMI threshold.
    
    Parameters
    ----------
    adata : AnnData
        Count matrix from kb-python output
    umi_threshold : float
        Minimum UMI counts to retain a barcode (typically from
        detect_knee_threshold or calculate_barcode_ranks)
    min_genes : int, optional
        Minimum genes detected to retain a barcode. If None, no
        gene count filter is applied (recommended: let downstream
        QC handle this).
    filter_genes_min_cells : int
        Remove genes detected in fewer than this many cells. Default: 1
        (removes genes with zero counts in all retained cells).
    inplace : bool
        If True, modify adata in place (default: False)
    
    Returns
    -------
    AnnData
        Filtered count matrix with empty droplets removed
    """
    if not inplace:
        adata = adata.copy()
    
    n_before = adata.n_obs
    n_genes_before = adata.n_vars
    
    # Calculate UMI counts per barcode
    if sparse.issparse(adata.X):
        umi_counts = np.asarray(adata.X.sum(axis=1)).flatten()
    else:
        umi_counts = np.asarray(adata.X.sum(axis=1)).flatten()
    
    # Apply UMI threshold (use > for consistency with DropletUtils)
    adata = adata[umi_counts > umi_threshold].copy()
    
    # Apply minimum genes filter if specified
    if min_genes is not None:
        sc.pp.filter_cells(adata, min_genes=min_genes)
    
    # Filter genes with zero counts
    if filter_genes_min_cells > 0:
        sc.pp.filter_genes(adata, min_cells=filter_genes_min_cells)
    
    n_after = adata.n_obs
    n_genes_after = adata.n_vars
    
    print(f"Filtered {n_before:,} -> {n_after:,} barcodes "
          f"({100 * n_after / n_before:.1f}% retained)")
    print(f"Filtered {n_genes_before:,} -> {n_genes_after:,} genes")
    
    return adata


# ==============================================================================
# Utility Functions
# ==============================================================================


def calculate_saturation_metrics(
    adata: ad.AnnData
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate library saturation metrics (UMI counts vs genes detected).
    
    Parameters
    ----------
    adata : AnnData
        Count matrix from kb-python output
    
    Returns
    -------
    tuple
        (umi_counts, genes_detected) arrays for each barcode
    """
    if sparse.issparse(adata.X):
        umi_counts = np.asarray(adata.X.sum(axis=1)).flatten()
        genes_detected = np.asarray((adata.X > 0).sum(axis=1)).flatten()
    else:
        umi_counts = np.asarray(adata.X.sum(axis=1)).flatten()
        genes_detected = np.asarray((adata.X > 0).sum(axis=1)).flatten()
    
    return umi_counts, genes_detected


def get_knee_data(adata: ad.AnnData) -> Tuple[np.ndarray, np.ndarray]:
    """
    Get data for knee plot visualization.
    
    Parameters
    ----------
    adata : AnnData
        Count matrix from kb-python output
    
    Returns
    -------
    tuple
        (sorted_umi_counts, barcode_ranks) for plotting
    """
    if sparse.issparse(adata.X):
        umi_counts = np.asarray(adata.X.sum(axis=1)).flatten()
    else:
        umi_counts = np.asarray(adata.X.sum(axis=1)).flatten()
    
    sorted_counts = np.sort(umi_counts)[::-1]
    ranks = np.arange(1, len(sorted_counts) + 1)
    
    return sorted_counts, ranks


def print_import_summary(adata: ad.AnnData, label: str = '') -> None:
    """
    Print summary statistics for kb-python import.
    
    Parameters
    ----------
    adata : AnnData
        Count matrix
    label : str
        Optional label for the summary
    """
    if sparse.issparse(adata.X):
        umi_counts = np.asarray(adata.X.sum(axis=1)).flatten()
        genes_detected = np.asarray((adata.X > 0).sum(axis=1)).flatten()
        total_counts = adata.X.sum()
    else:
        umi_counts = np.asarray(adata.X.sum(axis=1)).flatten()
        genes_detected = np.asarray((adata.X > 0).sum(axis=1)).flatten()
        total_counts = adata.X.sum()
    
    header = f"Summary{f' ({label})' if label else ''}"
    print(f"\n{'=' * 50}")
    print(header)
    print('=' * 50)
    print(f"Barcodes:        {adata.n_obs:,}")
    print(f"Genes:           {adata.n_vars:,}")
    print(f"Total UMIs:      {total_counts:,.0f}")
    print(f"Median UMIs:     {np.median(umi_counts):,.0f}")
    print(f"Median genes:    {np.median(genes_detected):,.0f}")
    print(f"UMI range:       {umi_counts.min():,.0f} - {umi_counts.max():,.0f}")
    print(f"Genes range:     {genes_detected.min():,.0f} - {genes_detected.max():,.0f}")
    print('=' * 50 + '\n')