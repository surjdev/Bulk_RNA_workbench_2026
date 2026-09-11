"""
R DESeq2 wrapper through rpy2 for BTW (FR-3 & FR-10).
Serves as the reference implementation engine for Differential Expression analysis.
"""

from __future__ import annotations

from typing import Any, Optional, Tuple

import numpy as np
import pandas as pd

from btw import logger
from btw.r_interop.bridge import (
    check_r_package,
    get_r_package_version,
    is_r_available,
    pandas_to_r_df,
    pandas_to_r_matrix,
    r_to_pandas_df,
    require_r_package,
)

if is_r_available():
    import rpy2.robjects as ro

    _r = ro.r
else:
    _r = None


def run_r_deseq2(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    contrast: Tuple[str, str, str],
    design: Optional[str] = None,
    alpha: float = 0.05,
    lfc_threshold: float = 1.0,
    shrink_lfc: bool = True,
    shrink_type: str = "apeglm",
    fit_type: str = "parametric",
) -> Tuple[pd.DataFrame, Any]:
    """
    Run R Bioconductor DESeq2 differential expression analysis via rpy2.

    Parameters
    ----------
    counts : pd.DataFrame
        Raw count matrix of shape (genes x samples). Must contain non-negative integers.
    metadata : pd.DataFrame
        Sample metadata indexed by sample identifier.
    contrast : tuple of (factor, test_level, reference_level)
        Comparison specification (e.g. ('condition', 'treated', 'control')).
    design : str, optional
        R design formula without leading tilde (e.g. 'condition' or 'batch + condition').
        Defaults to contrast factor.
    alpha : float, default=0.05
        FDR significance threshold.
    lfc_threshold : float, default=1.0
        Log2 fold change threshold for significance labeling.
    shrink_lfc : bool, default=True
        Whether to perform empirical Bayes shrinkage on log2 fold changes.
    shrink_type : str, default='apeglm'
        Shrinkage method: 'apeglm', 'ashr', or 'normal'.
    fit_type : str, default='parametric'
        Dispersion fit type: 'parametric', 'local', or 'mean'.

    Returns
    -------
    results_df : pd.DataFrame
        Standardized DE result table with columns:
        ['baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'regulation'].
    r_dds : R object
        Underlying R DESeqDataSet handle for direct user access.
    """
    require_r_package(
        "DESeq2",
        purpose="Reference Differential Expression analysis (FR-3)",
        alternative="Set engine='python' in run_de() to use PyDESeq2 instead.",
    )

    factor, test_lvl, ref_lvl = contrast
    if design is None:
        design = factor

    # Align samples
    common_samples = [s for s in counts.columns if s in metadata.index]
    if len(common_samples) == 0:
        raise ValueError("No matching sample IDs between counts and metadata.")

    aligned_counts = counts[common_samples].copy()
    aligned_meta = metadata.loc[common_samples].copy()

    # Ensure factor column is categorical with ref_lvl as baseline
    aligned_meta[factor] = aligned_meta[factor].astype(str)
    unique_lvls = list(aligned_meta[factor].unique())
    if ref_lvl not in unique_lvls or test_lvl not in unique_lvls:
        raise ValueError(
            f"Contrast levels '{test_lvl}' or '{ref_lvl}' not found in metadata column '{factor}'."
        )

    # Order levels so ref_lvl is first
    ordered_lvls = [ref_lvl] + [lvl for lvl in unique_lvls if lvl != ref_lvl]

    logger.info(
        f"Executing R DESeq2 (v{get_r_package_version('DESeq2')}) on {len(aligned_counts)} genes x {len(aligned_meta)} samples..."
    )
    logger.info(f"Contrast: {factor} ({test_lvl} vs {ref_lvl}), Design: ~{design}")

    # Convert to R objects
    r_counts = pandas_to_r_matrix(aligned_counts, is_integer=True)
    r_meta = pandas_to_r_df(aligned_meta)

    # Convert factor column to R factor with reference level
    _r.assign(".tmp_meta", r_meta)
    _r(
        f".tmp_meta[['{factor}']] <- factor(.tmp_meta[['{factor}']], levels=c({', '.join([repr(lvl) for lvl in ordered_lvls])}))"
    )
    r_meta = _r(".tmp_meta")

    # Create DESeqDataSet in R
    deseq2_mod = ro.packages.importr("DESeq2")
    design_formula = _r(f"as.formula('~{design}')")
    dds = deseq2_mod.DESeqDataSetFromMatrix(
        countData=r_counts, colData=r_meta, design=design_formula
    )

    # Run DESeq
    dds = deseq2_mod.DESeq(dds, fitType=fit_type, quiet=True)

    # Contrast vector: c("factor", "test", "reference")
    contrast_vec = ro.StrVector([factor, test_lvl, ref_lvl])
    res = deseq2_mod.results(dds, contrast=contrast_vec, alpha=alpha)

    # LFC shrinkage if requested
    if shrink_lfc:
        coef_name = f"{factor}_{test_lvl}_vs_{ref_lvl}"
        res_names = list(_r("resultsNames")(dds))
        if coef_name in res_names and shrink_type == "apeglm" and check_r_package("apeglm"):
            try:
                res = deseq2_mod.lfcShrink(dds, coef=coef_name, type="apeglm", quiet=True)
                logger.debug("Applied apeglm LFC shrinkage in R DESeq2.")
            except Exception as e:
                logger.warning(f"apeglm shrinkage failed ({e}); retaining standard Wald results.")
        elif check_r_package("ashr") and shrink_type == "ashr":
            try:
                res = deseq2_mod.lfcShrink(dds, contrast=contrast_vec, type="ashr", quiet=True)
                logger.debug("Applied ashr LFC shrinkage in R DESeq2.")
            except Exception as e:
                logger.warning(f"ashr shrinkage failed ({e}); retaining standard Wald results.")
        else:
            try:
                res = deseq2_mod.lfcShrink(dds, contrast=contrast_vec, type="normal", quiet=True)
                logger.debug("Applied normal LFC shrinkage in R DESeq2.")
            except Exception:
                pass

    # Convert results to DataFrame
    res_df = r_to_pandas_df(res)
    res_df.index.name = "gene_id"

    # Map / Standardize column names
    col_mapping = {
        "baseMean": "baseMean",
        "log2FoldChange": "log2FoldChange",
        "lfcSE": "lfcSE",
        "stat": "stat",
        "pvalue": "pvalue",
        "padj": "padj",
    }
    for old, new in col_mapping.items():
        if old not in res_df.columns:
            res_df[new] = np.nan

    # Add regulation classification
    is_sig = (
        (res_df["padj"].notna())
        & (res_df["padj"] <= alpha)
        & (res_df["log2FoldChange"].abs() >= lfc_threshold)
    )
    res_df["regulation"] = "NS"
    res_df.loc[is_sig & (res_df["log2FoldChange"] > 0), "regulation"] = "UP"
    res_df.loc[is_sig & (res_df["log2FoldChange"] < 0), "regulation"] = "DOWN"

    return res_df, dds


def run_r_vst(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str = "condition",
    blind: bool = True,
) -> pd.DataFrame:
    """
    Apply Variance Stabilizing Transformation (VST) via R DESeq2::vst.

    Parameters
    ----------
    counts : pd.DataFrame
        Raw counts matrix (genes x samples).
    metadata : pd.DataFrame
        Sample metadata.
    design : str, default='condition'
        Experimental design formula variable.
    blind : bool, default=True
        Whether transformation is blind to design.

    Returns
    -------
    pd.DataFrame
        VST-transformed expression matrix.
    """
    require_r_package("DESeq2", purpose="Variance Stabilizing Transformation (FR-2)")

    r_counts = pandas_to_r_matrix(counts, is_integer=True)
    r_meta = pandas_to_r_df(metadata)

    deseq2_mod = ro.packages.importr("DESeq2")
    design_formula = _r(f"as.formula('~{design}')")
    dds = deseq2_mod.DESeqDataSetFromMatrix(
        countData=r_counts, colData=r_meta, design=design_formula
    )

    vsd = deseq2_mod.vst(dds, blind=blind)
    vst_mat = _r("assay")(vsd)

    vst_df = r_to_pandas_df(vst_mat)
    vst_df.index = counts.index
    vst_df.columns = counts.columns
    return vst_df
