"""
Multiple testing correction utilities for Differential Expression Analysis (FR-3).
Integrates with statsmodels.stats.multitest for Benjamini-Hochberg, Bonferroni, and other FDR methods.
"""

from __future__ import annotations

from typing import List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from btw import logger

# Supported correction methods mapped to descriptions
SUPPORTED_METHODS = {
    "fdr_bh": "Benjamini/Hochberg (non-negative FDR)",
    "fdr_by": "Benjamini/Yekutieli (general FDR under dependence)",
    "bonferroni": "Bonferroni single-step FWER",
    "holm": "Holm-Bonferroni step-down FWER",
    "sidak": "Sidak single-step FWER",
    "holm-sidak": "Holm-Sidak step-down FWER",
    "fdr_tsbh": "Two-stage Benjamini/Hochberg FDR",
}


def adjust_pvalues(
    pvalues: Union[pd.Series, np.ndarray, List[float]],
    method: str = "fdr_bh",
    alpha: float = 0.05,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Apply multiple hypothesis testing correction to p-values, handling NaNs gracefully.

    Parameters
    ----------
    pvalues : Series, ndarray, or list of float
        Raw unadjusted p-values.
    method : str, default='fdr_bh'
        Multiple testing correction method. Options include 'fdr_bh', 'fdr_by',
        'bonferroni', 'holm', 'sidak', 'fdr_tsbh'.
    alpha : float, default=0.05
        Family-wise error rate or False Discovery Rate threshold.

    Returns
    -------
    rejected : np.ndarray of bool
        Boolean array indicating whether null hypothesis is rejected at level alpha.
    pvals_corrected : np.ndarray of float
        P-values adjusted for multiple testing. NaNs are preserved.
    """
    if method not in SUPPORTED_METHODS:
        valid = ", ".join(SUPPORTED_METHODS.keys())
        raise ValueError(f"Unsupported correction method '{method}'. Valid options: {valid}")

    pvals_arr = np.asarray(pvalues, dtype=float)
    if pvals_arr.size == 0:
        return np.array([], dtype=bool), np.array([], dtype=float)

    # Identify non-NaN and finite p-values
    valid_mask = np.isfinite(pvals_arr) & ~np.isnan(pvals_arr)

    rejected = np.zeros(pvals_arr.shape, dtype=bool)
    pvals_corrected = np.full(pvals_arr.shape, np.nan, dtype=float)

    if not np.any(valid_mask):
        logger.warning("All p-values are NaN or non-finite; returning unadjusted NaNs.")
        return rejected, pvals_corrected

    valid_pvals = pvals_arr[valid_mask]
    # Bound p-values within [0, 1]
    valid_pvals = np.clip(valid_pvals, 0.0, 1.0)

    reject_sub, padj_sub, _, _ = multipletests(valid_pvals, alpha=alpha, method=method)

    rejected[valid_mask] = reject_sub
    pvals_corrected[valid_mask] = padj_sub

    n_sig = np.sum(rejected)
    logger.debug(
        f"Multiple testing correction ({method}, alpha={alpha}): "
        f"{n_sig}/{len(valid_pvals)} tests rejected."
    )
    return rejected, pvals_corrected


def apply_multitest_correction(
    df: pd.DataFrame,
    pvalue_col: str = "pvalue",
    method: str = "fdr_bh",
    alpha: float = 0.05,
    output_padj_col: Optional[str] = None,
    output_sig_col: Optional[str] = None,
) -> pd.DataFrame:
    """
    Apply multiple testing correction to a DataFrame and add adjusted p-value and significance columns.

    Parameters
    ----------
    df : pd.DataFrame
        Differential expression table containing p-values.
    pvalue_col : str, default='pvalue'
        Name of column containing raw p-values.
    method : str, default='fdr_bh'
        Correction method from statsmodels.
    alpha : float, default=0.05
        Significance threshold.
    output_padj_col : str, optional
        Name of output column for adjusted p-values. Defaults to 'padj_{method}' or 'padj'.
    output_sig_col : str, optional
        Name of output column for boolean significance. Defaults to 'significant_{method}' or 'significant'.

    Returns
    -------
    pd.DataFrame
        Copy of input DataFrame with new adjusted p-value and significance columns.
    """
    if pvalue_col not in df.columns:
        raise KeyError(f"P-value column '{pvalue_col}' not found in DataFrame columns: {list(df.columns)}")

    padj_name = output_padj_col if output_padj_col is not None else "padj"
    sig_name = output_sig_col if output_sig_col is not None else "significant"

    rejected, padj = adjust_pvalues(df[pvalue_col].values, method=method, alpha=alpha)

    out_df = df.copy()
    out_df[padj_name] = padj
    out_df[sig_name] = rejected

    n_sig = rejected.sum()
    logger.info(
        f"Applied multiple testing correction ({method}, alpha={alpha}): "
        f"{n_sig} significant genes identified."
    )
    return out_df
