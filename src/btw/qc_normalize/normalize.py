"""
Normalization utilities for Bulk Transcriptomics (FR-2).
Thin helper wrapping PyDESeq2 median-of-ratios normalization, CPM, and log2 transforms.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional, Union

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
    size_factors: pd.Series  # sample_id -> size_factor
    method: str = "deseq2_median_of_ratios"
    dds: Optional[Any] = None  # PyDESeq2 DeseqDataSet if created
    engine: str = "python"

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
            logger.debug(
                f"PyDESeq2 deseq2_norm returned error ({e}); falling back to native median-of-ratios."
            )

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
    engine: str = "r",
    fallback_to_python: bool = True,
    **dds_kwargs,
) -> NormalizationResult:
    """
    Normalize count matrix using DESeq2 Median-of-Ratios method.
    Supports both reference R DESeq2 and Python PyDESeq2 engines.

    Parameters
    ----------
    counts : pd.DataFrame
        Gene count matrix of shape (n_genes, n_samples).
    metadata : pd.DataFrame, optional
        Sample metadata table.
    design_factors : str or list of str, optional
        Experimental design factor column(s).
    override_size_factors : pd.Series, optional
        User-provided custom size factors overriding default calculation.
    return_dds : bool, default=False
        If True and metadata is provided, returns underlying DeseqDataSet.
    engine : str, default='r'
        Engine to use: 'r' (calls R DESeq2 via rpy2) or 'python' (PyDESeq2).
    fallback_to_python : bool, default=True
        If True and R engine is unavailable, automatically fall back to Python.
    **dds_kwargs
        Keyword arguments passed to DeseqDataSet constructor.

    Returns
    -------
    NormalizationResult
        Dataclass containing normalized counts, size factors, method name,
        and underlying DeseqDataSet (if requested).
    """
    used_engine = engine.lower()

    # Route R engine
    if used_engine == "r":
        try:
            import rpy2.robjects as ro

            from btw.r_interop.bridge import (
                check_r_package,
                pandas_to_r_df,
                pandas_to_r_matrix,
                require_r_package,
            )

            if check_r_package("DESeq2"):
                logger.info("Calculating size factors via R DESeq2 (estimateSizeFactors)...")
                deseq2_mod = ro.packages.importr("DESeq2")
                r_counts = pandas_to_r_matrix(counts, is_integer=True)

                if metadata is not None:
                    r_meta = pandas_to_r_df(metadata.loc[counts.columns])
                    factor = design_factors if design_factors is not None else metadata.columns[0]
                    design_formula = ro.r(f"as.formula('~{factor}')")
                    r_dds = deseq2_mod.DESeqDataSetFromMatrix(
                        countData=r_counts, colData=r_meta, design=design_formula
                    )
                else:
                    # Dummy metadata
                    dummy_meta = pd.DataFrame({"sample": counts.columns}, index=counts.columns)
                    r_meta = pandas_to_r_df(dummy_meta)
                    r_dds = deseq2_mod.DESeqDataSetFromMatrix(
                        countData=r_counts, colData=r_meta, design=ro.r("as.formula('~1')")
                    )

                r_dds = deseq2_mod.estimateSizeFactors(r_dds)
                sfs = list(ro.r("sizeFactors")(r_dds))
                size_factors = pd.Series(sfs, index=counts.columns, name="size_factors")
                logger.info(
                    f"Calculated R DESeq2 size factors (range: {size_factors.min():.3f} - {size_factors.max():.3f})"
                )

                norm_counts = counts.div(size_factors, axis=1)
                return NormalizationResult(
                    normalized_counts=norm_counts,
                    size_factors=size_factors,
                    method="deseq2_median_of_ratios",
                    dds=r_dds if return_dds else None,
                    engine="r",
                )
            else:
                if not fallback_to_python:
                    require_r_package("DESeq2", purpose="R DESeq2 normalization (FR-2)")
                logger.info(
                    "R package 'DESeq2' not found; falling back to Python normalization engine."
                )
                used_engine = "python"
        except Exception as e:
            if not fallback_to_python:
                raise e
            logger.info(
                f"R normalization encountered exception ({e}); falling back to Python normalization engine."
            )
            used_engine = "python"

    # Python engine execution
    if override_size_factors is not None:
        size_factors = override_size_factors.loc[counts.columns]
        logger.info("Using user-provided custom size factors for normalization.")
    else:
        size_factors = compute_size_factors(counts)
        logger.info(
            f"Calculated Python DESeq2 size factors (range: {size_factors.min():.3f} - {size_factors.max():.3f})"
        )

    norm_counts = counts.div(size_factors, axis=1)

    dds = None
    if return_dds and metadata is not None:
        if not HAS_PYDESEQ2:
            logger.warning("pydeseq2 is not installed; cannot return DeseqDataSet object.")
        else:
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
        engine="python",
    )


def vst_transform(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str = "condition",
    engine: str = "r",
    fallback_to_python: bool = True,
) -> pd.DataFrame:
    """
    Apply Variance Stabilizing Transformation (VST).
    Calls R DESeq2::vst when engine='r', with fallback to log2-normalized expression.
    """
    if engine.lower() == "r":
        try:
            from btw.r_interop.bridge import check_r_package, require_r_package
            from btw.r_interop.deseq2_r import run_r_vst

            if check_r_package("DESeq2"):
                return run_r_vst(counts, metadata, design=design)
            elif not fallback_to_python:
                require_r_package("DESeq2", purpose="VST transformation")
        except Exception as e:
            if not fallback_to_python:
                raise e

    # Fallback to log2(norm_counts + 1)
    norm_res = normalize_deseq2(counts, metadata=metadata, engine="python")
    return norm_res.log2_counts


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
