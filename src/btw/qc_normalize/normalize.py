"""
Normalization utilities for Bulk Transcriptomics (FR-2).
Thin helper wrapping PyDESeq2 median-of-ratios normalization, CPM, and log2 transforms.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple, Union

import numpy as np
import pandas as pd
from btw import logger

try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.preprocessing import deseq2_norm
    HAS_PYDESEQ2 = True
except ImportError:
    HAS_PYDESEQ2 = False
    DeseqDataSet = None
    deseq2_norm = None


@dataclass
class NormalizationResult:
    """Container holding normalized counts, size factors, and optional underlying DeseqDataSet."""
    normalized_counts: pd.DataFrame  # genes x samples
    size_factors: pd.Series          # sample_id -> size_factor
    method: str
    dds: Optional[Any] = None        # PyDESeq2 DeseqDataSet if created

    @property
    def log2_counts(self) -> pd.DataFrame:
        """Convenience property to obtain log2(normalized_counts + 1)."""
        return np.log2(self.normalized_counts + 1.0)


def compute_size_factors(
    counts: pd.DataFrame,
) -> pd.Series:
    """
    Calculate DESeq2 size factors using the Median-of-Ratios method.

    Parameters
    ----------
    counts : pd.DataFrame
        Gene count matrix of shape (n_genes, n_samples).

    Returns
    -------
    pd.Series
        Size factor for each sample, indexed by sample ID.
    """
    if HAS_PYDESEQ2 and deseq2_norm is not None:
        try:
            # PyDESeq2 expects samples x genes
            _, sfs = deseq2_norm(counts.T)
            return pd.Series(sfs, index=counts.columns, name="size_factors")
        except Exception as e:
            logger.debug(f"PyDESeq2 deseq2_norm returned error ({e}); falling back to native median-of-ratios.")

    # Native implementation of DESeq2 median-of-ratios
    # 1. Compute geometric mean for each gene across samples
    with np.errstate(divide="ignore", invalid="ignore"):
        log_counts = np.log(counts.values.astype(float))
        # Mask non-positive counts
        valid_mask = counts.values > 0
        log_counts_masked = np.where(valid_mask, log_counts, np.nan)
        # Geometric mean in log-space is mean of logs
        log_geom_means = np.nanmean(log_counts_masked, axis=1, keepdims=True)
        geom_means = np.exp(log_geom_means)

    # Genes with finite, non-zero geometric mean in all samples with data
    finite_genes = np.isfinite(geom_means.ravel()) & (geom_means.ravel() > 0)
    if not np.any(finite_genes):
        # Fallback to total count scaling if no gene is positive everywhere
        lib_sizes = counts.sum(axis=0)
        sfs = lib_sizes / lib_sizes.median()
        return pd.Series(sfs, index=counts.columns, name="size_factors")

    # 2. Ratio of counts to geometric mean
    ratios = counts.values[finite_genes, :] / geom_means[finite_genes, :]
    # 3. Median ratio per sample
    size_factors = np.nanmedian(ratios, axis=0)

    # Normalize so geometric mean of size factors is approximately 1
    geo_mean_sf = np.exp(np.mean(np.log(size_factors[size_factors > 0])))
    if geo_mean_sf > 0:
        size_factors = size_factors / geo_mean_sf

    return pd.Series(size_factors, index=counts.columns, name="size_factors")


def normalize_deseq2(
    counts: pd.DataFrame,
    metadata: Optional[pd.DataFrame] = None,
    design_factors: Optional[Union[str, list[str]]] = None,
    override_size_factors: Optional[pd.Series] = None,
    return_dds: bool = False,
    **dds_kwargs,
) -> NormalizationResult:
    """
    Thin wrapper around PyDESeq2 normalization with full parameter override capability.

    Parameters
    ----------
    counts : pd.DataFrame
        Gene count matrix of shape (n_genes, n_samples).
    metadata : pd.DataFrame, optional
        Sample metadata table.
    design_factors : str or list of str, optional
        Experimental design factor column(s) for DeseqDataSet if creating dds.
    override_size_factors : pd.Series, optional
        User-provided custom size factors overriding default calculation.
    return_dds : bool, default=False
        If True and metadata is provided, initializes and returns PyDESeq2 DeseqDataSet.
    **dds_kwargs
        Keyword arguments passed to DeseqDataSet constructor.

    Returns
    -------
    NormalizationResult
        Dataclass containing normalized counts, size factors, method name,
        and underlying DeseqDataSet (if requested).
    """
    if override_size_factors is not None:
        size_factors = override_size_factors.loc[counts.columns]
        logger.info("Using user-provided custom size factors for normalization.")
    else:
        size_factors = compute_size_factors(counts)
        logger.info(f"Calculated DESeq2 size factors (range: {size_factors.min():.3f} - {size_factors.max():.3f})")

    # Divide raw counts by size factors
    norm_counts = counts.div(size_factors, axis=1)

    dds = None
    if return_dds and metadata is not None:
        if not HAS_PYDESEQ2:
            logger.warning("pydeseq2 is not installed; cannot return DeseqDataSet object.")
        else:
            # Prepare PyDESeq2 DeseqDataSet (samples x genes)
            factor = design_factors if design_factors is not None else metadata.columns[0]
            dds = DeseqDataSet(
                counts=counts.T,
                metadata=metadata.loc[counts.columns],
                design_factors=factor,
                **dds_kwargs,
            )
            dds.fit_size_factors()
            logger.info("Instantiated PyDESeq2 DeseqDataSet with size factors fitted.")

    return NormalizationResult(
        normalized_counts=norm_counts,
        size_factors=size_factors,
        method="deseq2_median_of_ratios",
        dds=dds,
    )


def normalize_cpm(
    counts: pd.DataFrame,
    log2_transform: bool = False,
    prior_count: float = 1.0,
) -> pd.DataFrame:
    """
    Compute Counts Per Million (CPM) normalization.

    Parameters
    ----------
    counts : pd.DataFrame
        Raw count matrix of shape (n_genes, n_samples).
    log2_transform : bool, default=False
        Whether to return log2(CPM + prior_count).
    prior_count : float, default=1.0
        Pseudo-count added before log2 transformation.

    Returns
    -------
    pd.DataFrame
        CPM-normalized expression matrix.
    """
    lib_sizes = counts.sum(axis=0)
    cpm = (counts * 1e6).div(lib_sizes, axis=1)
    if log2_transform:
        cpm = np.log2(cpm + prior_count)
    return cpm
