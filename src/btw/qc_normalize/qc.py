"""
Quality Control (QC) metrics and filtering module for Bulk Transcriptomics (FR-2).
Calculates sample-level and gene-level summary statistics.
"""

from __future__ import annotations

from typing import Any, Dict, Optional, Tuple

import numpy as np
import pandas as pd
from btw import logger


def compute_sample_qc(
    counts: pd.DataFrame,
    min_count: int = 1,
) -> pd.DataFrame:
    """
    Compute sample-level Quality Control metrics.

    Parameters
    ----------
    counts : pd.DataFrame
        Raw count matrix of shape (n_genes, n_samples).
    min_count : int, default=1
        Threshold count to consider a gene detected/expressed.

    Returns
    -------
    pd.DataFrame
        Table indexed by sample ID with columns:
        - `library_size`: total counts per sample
        - `detected_genes`: count of genes with expression >= min_count
        - `detection_rate_pct`: percentage of detected genes
        - `mean_count`: mean expression across all genes
        - `median_count`: median expression
        - `q75_count`: 75th percentile
        - `max_count`: highest gene count in sample
        - `zero_genes_pct`: percentage of unexpressed genes (count < min_count)
    """
    total_genes = counts.shape[0]
    if total_genes == 0:
        raise ValueError("Count matrix has 0 genes.")

    lib_sizes = counts.sum(axis=0)
    detected = (counts >= min_count).sum(axis=0)
    detection_pct = (detected / total_genes) * 100.0
    zero_pct = 100.0 - detection_pct

    stats_df = pd.DataFrame(
        {
            "library_size": lib_sizes.astype(np.int64),
            "detected_genes": detected.astype(int),
            "detection_rate_pct": np.round(detection_pct, 2),
            "mean_count": np.round(counts.mean(axis=0), 2),
            "median_count": counts.median(axis=0),
            "q75_count": counts.quantile(0.75, axis=0),
            "max_count": counts.max(axis=0),
            "zero_genes_pct": np.round(zero_pct, 2),
        },
        index=counts.columns,
    )
    stats_df.index.name = "sample_id"

    logger.info(
        f"Calculated sample QC for {counts.shape[1]} samples. "
        f"Mean library size: {stats_df['library_size'].mean():,.0f}, "
        f"Mean detection rate: {stats_df['detection_rate_pct'].mean():.1f}%"
    )
    return stats_df


def compute_gene_qc(
    counts: pd.DataFrame,
    min_count: int = 1,
) -> pd.DataFrame:
    """
    Compute gene-level Quality Control metrics across all samples.

    Parameters
    ----------
    counts : pd.DataFrame
        Raw count matrix of shape (n_genes, n_samples).
    min_count : int, default=1
        Threshold to consider a gene expressed in a sample.

    Returns
    -------
    pd.DataFrame
        Table indexed by gene identifier with columns:
        - `total_counts`: sum across all samples
        - `mean_count`: average count per sample
        - `variance`: sample variance
        - `dispersion_ratio`: variance / mean ratio
        - `n_samples_expressing`: number of samples with count >= min_count
        - `fraction_samples_expressing`: proportion of samples expressing the gene
    """
    n_samples = counts.shape[1]
    tot = counts.sum(axis=1)
    mean = counts.mean(axis=1)
    var = counts.var(axis=1)
    # Avoid division by zero in dispersion
    disp = np.where(mean > 0, var / mean, 0.0)
    n_expressing = (counts >= min_count).sum(axis=1)
    frac_expressing = n_expressing / n_samples

    gene_df = pd.DataFrame(
        {
            "total_counts": tot,
            "mean_count": np.round(mean, 2),
            "variance": np.round(var, 2),
            "dispersion_ratio": np.round(disp, 2),
            "n_samples_expressing": n_expressing,
            "fraction_samples_expressing": np.round(frac_expressing, 3),
        },
        index=counts.index,
    )
    gene_df.index.name = "gene_id"
    return gene_df


def filter_low_expression_genes(
    counts: pd.DataFrame,
    min_counts: int = 10,
    min_samples: int = 3,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Filter out low-abundance genes that lack sufficient statistical power for DE.

    Parameters
    ----------
    counts : pd.DataFrame
        Raw count matrix of shape (n_genes, n_samples).
    min_counts : int, default=10
        Minimum count required per sample.
    min_samples : int, default=3
        Minimum number of samples that must reach `min_counts`.

    Returns
    -------
    filtered_counts : pd.DataFrame
        Subsetted count matrix with retained genes.
    filtering_summary : pd.DataFrame
        Summary table detailing before/after gene counts and retention rates.
    """
    mask = (counts >= min_counts).sum(axis=1) >= min_samples
    filtered_counts = counts.loc[mask]

    n_orig = counts.shape[0]
    n_retained = filtered_counts.shape[0]
    n_removed = n_orig - n_retained
    pct_retained = (n_retained / n_orig * 100.0) if n_orig > 0 else 0.0

    summary_df = pd.DataFrame(
        {
            "metric": [
                "Original genes",
                "Retained genes",
                "Removed low-count genes",
                "Retention rate (%)",
                "Filter threshold",
            ],
            "value": [
                str(n_orig),
                str(n_retained),
                str(n_removed),
                f"{pct_retained:.2f}%",
                f">= {min_counts} counts in >= {min_samples} samples",
            ],
        }
    )

    logger.info(
        f"Filtered low-expression genes: {n_retained}/{n_orig} retained ({pct_retained:.1f}%), "
        f"{n_removed} removed."
    )
    return filtered_counts, summary_df
