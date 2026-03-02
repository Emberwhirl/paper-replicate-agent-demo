#!/usr/bin/env python3
"""
Core utility functions for single-cell RNA-seq data integration.

This module provides building blocks for dataset integration using
scanpy's ingest and BBKNN methods.
"""

import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd


def intersect_genes(adata_ref, adata_query, verbose=True):
    """
    Intersect gene space between reference and query datasets.

    Parameters
    ----------
    adata_ref : AnnData
        Reference annotated data matrix
    adata_query : AnnData
        Query annotated data matrix
    verbose : bool
        Print intersection statistics (default: True)

    Returns
    -------
    tuple of AnnData
        (adata_ref, adata_query) subsetted to common genes
    """
    var_names = adata_ref.var_names.intersection(adata_query.var_names)

    if verbose:
        print(f"Reference genes: {adata_ref.n_vars}")
        print(f"Query genes: {adata_query.n_vars}")
        print(f"Common genes: {len(var_names)}")

    if len(var_names) == 0:
        raise ValueError("No common genes found between reference and query datasets. "
                        "Check if gene naming conventions differ (e.g., ENSEMBL vs symbols).")

    adata_ref_subset = adata_ref[:, var_names].copy()
    adata_query_subset = adata_query[:, var_names].copy()

    return adata_ref_subset, adata_query_subset


def train_reference(adata, n_pcs=50, n_neighbors=15, verbose=True):
    """
    Train embedding model on reference dataset.

    Computes PCA, neighbors graph, and UMAP on reference data.
    Modifies adata in place.

    Parameters
    ----------
    adata : AnnData
        Reference annotated data matrix (modified in place)
    n_pcs : int
        Number of principal components (default: 50)
    n_neighbors : int
        Number of neighbors for graph construction (default: 15)
    verbose : bool
        Print progress messages (default: True)

    Returns
    -------
    None
        Modifies adata in place
    """
    if verbose:
        print("Training reference model...")
        print(f"  Computing PCA (n_pcs={n_pcs})...")

    sc.pp.pca(adata, n_comps=n_pcs)

    if verbose:
        print(f"  Computing neighbors (n_neighbors={n_neighbors})...")

    sc.pp.neighbors(adata, n_neighbors=n_neighbors)

    if verbose:
        print("  Computing UMAP...")

    sc.tl.umap(adata)

    if verbose:
        print("  Reference model trained successfully")


def ingest_query(adata_query, adata_ref, label_key='celltype', embedding_method='umap',
                preserve_original=True, verbose=True):
    """
    Map query cells onto reference embedding and transfer labels.

    Parameters
    ----------
    adata_query : AnnData
        Query annotated data matrix (modified in place)
    adata_ref : AnnData
        Reference annotated data matrix (must have PCA and neighbors computed)
    label_key : str
        Observation column to transfer from reference (default: 'celltype')
    embedding_method : str
        Embedding method for ingest: 'umap' or 'pca' (default: 'umap')
    preserve_original : bool
        Save original labels as {label_key}_orig (default: True)
    verbose : bool
        Print progress messages (default: True)

    Returns
    -------
    None
        Modifies adata_query in place
    """
    if verbose:
        print(f"Ingesting query data (label_key='{label_key}')...")

    # Preserve original labels if requested and they exist
    if preserve_original and label_key in adata_query.obs:
        adata_query.obs[f'{label_key}_orig'] = adata_query.obs[label_key].copy()
        if verbose:
            print(f"  Preserved original labels as '{label_key}_orig'")

    # Run ingest
    sc.tl.ingest(adata_query, adata_ref, obs=label_key, embedding_method=embedding_method)

    # Copy color scheme if available
    color_key = f'{label_key}_colors'
    if color_key in adata_ref.uns:
        adata_query.uns[color_key] = adata_ref.uns[color_key]

    if verbose:
        print(f"  Transferred '{label_key}' labels to {adata_query.n_obs} query cells")


def integrate_bbknn(adata, batch_key='batch', neighbors_within_batch=3, verbose=True):
    """
    Apply BBKNN batch correction.

    Replaces standard neighbors graph with batch-balanced k-nearest neighbors.
    Modifies adata in place.

    Parameters
    ----------
    adata : AnnData
        Combined annotated data matrix with batch column (modified in place)
    batch_key : str
        Observation column indicating batch (default: 'batch')
    neighbors_within_batch : int
        Number of neighbors to use from each batch (default: 3)
    verbose : bool
        Print progress messages (default: True)

    Returns
    -------
    None
        Modifies adata in place

    Raises
    ------
    ImportError
        If bbknn is not installed
    ValueError
        If batch_key not found in adata.obs
    """
    try:
        import bbknn
    except ImportError:
        raise ImportError("BBKNN not installed. Install with: pip install bbknn")

    if batch_key not in adata.obs:
        raise ValueError(f"Batch key '{batch_key}' not found in adata.obs. "
                        f"Available columns: {list(adata.obs.columns)}")

    n_batches = adata.obs[batch_key].nunique()

    if verbose:
        print(f"Applying BBKNN integration...")
        print(f"  Batch key: '{batch_key}' ({n_batches} batches)")
        print(f"  Neighbors within batch: {neighbors_within_batch}")
        print("  Computing PCA...")

    sc.pp.pca(adata)

    if verbose:
        print("  Computing batch-balanced neighbors...")

    bbknn.bbknn(adata, batch_key=batch_key,
                neighbors_within_batch=neighbors_within_batch)

    if verbose:
        print("  Computing UMAP on corrected graph...")

    sc.tl.umap(adata)

    if verbose:
        print("  BBKNN integration complete")


def combine_datasets(adata_ref, adata_query, batch_key='batch',
                    ref_label='reference', query_label='query',
                    label_key='celltype', verbose=True):
    """
    Concatenate reference and query datasets for joint visualization.

    Parameters
    ----------
    adata_ref : AnnData
        Reference annotated data matrix
    adata_query : AnnData
        Query annotated data matrix (should have ingest applied)
    batch_key : str
        Name for batch observation column (default: 'batch')
    ref_label : str
        Label for reference batch (default: 'reference')
    query_label : str
        Label for query batch (default: 'query')
    label_key : str
        Cell type annotation column to harmonize (default: 'celltype')
    verbose : bool
        Print progress messages (default: True)

    Returns
    -------
    AnnData
        Combined dataset with harmonized categories
    """
    if verbose:
        print("Combining datasets...")
        print(f"  Reference: {adata_ref.n_obs} cells")
        print(f"  Query: {adata_query.n_obs} cells")

    adata_concat = ad.concat(
        [adata_ref, adata_query],
        label=batch_key,
        keys=[ref_label, query_label]
    )

    # Fix category ordering to match reference (with handling for novel categories in query)
    if label_key in adata_concat.obs and label_key in adata_ref.obs:
        ref_categories = list(adata_ref.obs[label_key].cat.categories)
        concat_categories = set(adata_concat.obs[label_key].cat.categories)
        
        # Find any new categories from query not in reference
        new_categories = concat_categories - set(ref_categories)
        
        # Combine: reference categories first, then any new ones
        all_categories = ref_categories + sorted(list(new_categories))
        
        try:
            adata_concat.obs[label_key] = (
                adata_concat.obs[label_key]
                .astype('category')
                .cat.reorder_categories(all_categories, ordered=False)
            )
        except ValueError:
            # If reordering fails, just ensure it's categorical
            adata_concat.obs[label_key] = adata_concat.obs[label_key].astype('category')

        # Copy colors if available
        color_key = f'{label_key}_colors'
        if color_key in adata_ref.uns:
            adata_concat.uns[color_key] = adata_ref.uns[color_key]

    if verbose:
        print(f"  Combined: {adata_concat.n_obs} total cells")

    return adata_concat


def evaluate_label_transfer(adata, transferred_key='celltype', original_key='celltype_orig',
                           conserved_only=False, verbose=True):
    """
    Evaluate label transfer quality using confusion matrix.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with transferred and original labels
    transferred_key : str
        Column name for transferred labels (default: 'celltype')
    original_key : str
        Column name for original labels (default: 'celltype_orig')
    conserved_only : bool
        Only evaluate cell types present in both datasets (default: False)
    verbose : bool
        Print evaluation summary (default: True)

    Returns
    -------
    pd.DataFrame
        Confusion matrix (rows: transferred, columns: original)
    """
    if transferred_key not in adata.obs:
        raise ValueError(f"Transferred label key '{transferred_key}' not found in adata.obs")
    if original_key not in adata.obs:
        raise ValueError(f"Original label key '{original_key}' not found in adata.obs. "
                        "Make sure to preserve original labels before label transfer.")

    obs = adata.obs[[transferred_key, original_key]].dropna()

    if len(obs) == 0:
        raise ValueError("No cells with both transferred and original labels found. "
                        "Check if labels are missing or all NaN.")

    if conserved_only:
        # Find cell types present in both
        transferred_types = set(obs[transferred_key].unique())
        original_types = set(obs[original_key].unique())
        conserved = transferred_types.intersection(original_types)

        if verbose:
            print(f"Conserved cell types: {len(conserved)}")

        if len(conserved) == 0:
            raise ValueError("No conserved cell types found between transferred and original labels.")

        obs = obs[
            obs[transferred_key].isin(conserved) &
            obs[original_key].isin(conserved)
        ]

        if len(obs) == 0:
            raise ValueError("No cells remaining after filtering to conserved cell types.")

    confusion = pd.crosstab(obs[transferred_key], obs[original_key])

    if verbose:
        # Calculate accuracy
        diagonal = np.diag(confusion.reindex(
            index=confusion.columns, fill_value=0
        ).values)
        total = confusion.values.sum()
        accuracy = (diagonal.sum() / total) * 100 if total > 0 else 0

        print(f"\nLabel Transfer Evaluation:")
        print(f"  Cells evaluated: {len(obs)}")
        print(f"  Accuracy: {accuracy:.1f}%")

    return confusion


def print_integration_summary(adata, label=''):
    """
    Print summary statistics for integrated dataset.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    label : str
        Label to prepend to output
    """
    if label:
        print(f"\n{label}:")
    print(f"  Cells: {adata.n_obs}")
    print(f"  Genes: {adata.n_vars}")

    # Check for batch information
    if 'batch' in adata.obs:
        batch_counts = adata.obs['batch'].value_counts()
        print(f"  Batches: {len(batch_counts)}")
        for batch, count in batch_counts.items():
            print(f"    {batch}: {count} cells")

    # Check for cell type information
    for key in ['celltype', 'cell_type', 'CellType']:
        if key in adata.obs:
            n_types = adata.obs[key].nunique()
            print(f"  Cell types ({key}): {n_types}")
            break

    # Check for embeddings
    embeddings = []
    if 'X_pca' in adata.obsm:
        embeddings.append(f"PCA ({adata.obsm['X_pca'].shape[1]})")
    if 'X_umap' in adata.obsm:
        embeddings.append("UMAP")
    if embeddings:
        print(f"  Embeddings: {', '.join(embeddings)}")
