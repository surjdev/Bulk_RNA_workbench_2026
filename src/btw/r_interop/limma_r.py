"""
R limma wrapper through rpy2 for BTW (FR-3 & FR-10).
Implements limma-voom and limma-trend for count-based and continuous transcriptomics.
"""

from __future__ import annotations

from typing import Any, Tuple

import numpy as np
import pandas as pd

from btw import logger
from btw.r_interop.bridge import (
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


def run_r_limma(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    contrast: Tuple[str, str, str],
    method: str = "voom",
    alpha: float = 0.05,
    lfc_threshold: float = 1.0,
) -> Tuple[pd.DataFrame, Any]:
    """
    Run R Bioconductor limma (voom or trend) differential expression analysis via rpy2.

    Parameters
    ----------
    counts : pd.DataFrame
        Gene expression counts (genes x samples).
    metadata : pd.DataFrame
        Sample annotations.
    contrast : tuple of (factor, test_level, reference_level)
        Comparison contrast.
    method : str, default='voom'
        'voom' (precision weights on count data) or 'trend' (logCPM with eBayes trend).
    alpha : float, default=0.05
        FDR significance threshold.
    lfc_threshold : float, default=1.0
        Log2 fold change threshold.

    Returns
    -------
    results_df : pd.DataFrame
        Standardized DE results table.
    fit : R object
        Underlying R MArrayLM object.
    """
    require_r_package(
        "limma",
        purpose=f"limma-{method} differential expression (FR-3)",
        alternative="Use method='deseq2' with R DESeq2 or PyDESeq2.",
    )
    if method == "voom":
        require_r_package("edgeR", purpose="edgeR DGEList for limma-voom normalization")

    factor, test_lvl, ref_lvl = contrast
    common_samples = [s for s in counts.columns if s in metadata.index]
    aligned_counts = counts[common_samples].copy()
    aligned_meta = metadata.loc[common_samples].copy()

    logger.info(
        f"Executing R limma-{method} (v{get_r_package_version('limma')}) for contrast: {factor} ({test_lvl} vs {ref_lvl})..."
    )

    r_counts = pandas_to_r_matrix(aligned_counts)
    r_meta = pandas_to_r_df(aligned_meta)

    _r.assign(".l_counts", r_counts)
    _r.assign(".l_meta", r_meta)

    # Design matrix: ~ 0 + factor
    _r(f".group <- factor(.l_meta[['{factor}']])")
    _r(".design <- model.matrix(~ 0 + .group)")
    _r("colnames(.design) <- levels(.group)")

    ro.packages.importr("limma")

    if method == "voom":
        ro.packages.importr("edgeR")
        _r(".dge <- edgeR::DGEList(counts=.l_counts)")
        _r(".dge <- edgeR::calcNormFactors(.dge)")
        _r(".v <- limma::voom(.dge, .design, plot=FALSE)")
        _r(".fit <- limma::lmFit(.v, .design)")
    else:
        # limma-trend
        _r(".logcpm <- edgeR::cpm(.l_counts, log=TRUE, prior.count=3)")
        _r(".fit <- limma::lmFit(.logcpm, .design)")

    # Build contrast
    contrast_str = f"{test_lvl} - {ref_lvl}"
    _r(f".cm <- limma::makeContrasts({contrast_str}, levels=.design)")
    _r(".fit2 <- limma::contrasts.fit(.fit, .cm)")

    if method == "trend":
        _r(".fit2 <- limma::eBayes(.fit2, trend=TRUE)")
    else:
        _r(".fit2 <- limma::eBayes(.fit2)")

    r_top = _r(f"limma::topTable(.fit2, coef='{contrast_str}', number=Inf, adjust.method='BH')")
    res_df = r_to_pandas_df(r_top)
    res_df.index.name = "gene_id"

    # Map column names to standard schema
    # limma topTable columns: logFC, AveExpr, t, P.Value, adj.P.Val, B
    rename_cols = {
        "logFC": "log2FoldChange",
        "AveExpr": "baseMean",
        "t": "stat",
        "P.Value": "pvalue",
        "adj.P.Val": "padj",
    }
    for old, new in rename_cols.items():
        if old in res_df.columns:
            res_df[new] = res_df[old]

    res_df["lfcSE"] = np.nan
    is_sig = (
        (res_df["padj"].notna())
        & (res_df["padj"] <= alpha)
        & (res_df["log2FoldChange"].abs() >= lfc_threshold)
    )
    res_df["regulation"] = "NS"
    res_df.loc[is_sig & (res_df["log2FoldChange"] > 0), "regulation"] = "UP"
    res_df.loc[is_sig & (res_df["log2FoldChange"] < 0), "regulation"] = "DOWN"

    return res_df, _r(".fit2")
