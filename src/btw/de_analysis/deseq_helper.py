"""
PyDESeq2 thin helper and differential expression analysis runner (FR-3).
Reduces boilerplate without obfuscating original PyDESeq2 objects (DeseqDataSet, DeseqStats).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from btw import logger
from btw.io.validator import validate_bulk_data

try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    HAS_PYDESEQ2 = True
except ImportError:
    HAS_PYDESEQ2 = False
    DeseqDataSet = None
    DeseqStats = None


@dataclass
class DEResult:
    """
    Standardized container holding Differential Expression analysis results
    while keeping the original PyDESeq2 objects directly accessible.
    """
    results_df: pd.DataFrame
    contrast: Tuple[str, str, str]  # (factor, test_level, ref_level)
    design_factor: str
    alpha: float
    lfc_threshold: float
    dds: Optional[Any] = None       # Original DeseqDataSet
    stat_res: Optional[Any] = None  # Original DeseqStats

    @property
    def genes(self) -> List[str]:
        """List of all tested gene identifiers."""
        return list(self.results_df.index)

    def get_degs(
        self,
        padj_cutoff: Optional[float] = None,
        lfc_cutoff: Optional[float] = None,
    ) -> pd.DataFrame:
        """
        Filter DE table for statistically significant differentially expressed genes.

        Parameters
        ----------
        padj_cutoff : float, optional
            Maximum adjusted p-value. Defaults to self.alpha.
        lfc_cutoff : float, optional
            Minimum absolute log2 fold-change. Defaults to self.lfc_threshold.

        Returns
        -------
        pd.DataFrame
            Subset of results_df containing significant genes.
        """
        alpha = padj_cutoff if padj_cutoff is not None else self.alpha
        lfc = lfc_cutoff if lfc_cutoff is not None else self.lfc_threshold

        df = self.results_df
        mask = (df["padj"].notna()) & (df["padj"] <= alpha) & (df["log2FoldChange"].abs() >= lfc)
        return df.loc[mask]

    @property
    def padj_threshold(self) -> float:
        """Significance adjusted p-value threshold."""
        return self.alpha

    def significant_genes(
        self,
        padj_cutoff: Optional[float] = None,
        lfc_cutoff: Optional[float] = None,
    ) -> List[str]:
        """List of gene identifiers for all significant DEGs."""
        return list(self.get_degs(padj_cutoff, lfc_cutoff).index)

    def up_regulated(
        self,
        padj_cutoff: Optional[float] = None,
        lfc_cutoff: Optional[float] = None,
    ) -> pd.DataFrame:
        """Subset of results_df containing significantly up-regulated genes."""
        degs = self.get_degs(padj_cutoff, lfc_cutoff)
        return degs.loc[degs["log2FoldChange"] > 0]

    def down_regulated(
        self,
        padj_cutoff: Optional[float] = None,
        lfc_cutoff: Optional[float] = None,
    ) -> pd.DataFrame:
        """Subset of results_df containing significantly down-regulated genes."""
        degs = self.get_degs(padj_cutoff, lfc_cutoff)
        return degs.loc[degs["log2FoldChange"] < 0]

    def get_up_genes(
        self,
        padj_cutoff: Optional[float] = None,
        lfc_cutoff: Optional[float] = None,
    ) -> List[str]:
        """Return list of significantly up-regulated gene IDs."""
        degs = self.get_degs(padj_cutoff, lfc_cutoff)
        return list(degs.loc[degs["log2FoldChange"] > 0].index)

    def get_down_genes(
        self,
        padj_cutoff: Optional[float] = None,
        lfc_cutoff: Optional[float] = None,
    ) -> List[str]:
        """Return list of significantly down-regulated gene IDs."""
        degs = self.get_degs(padj_cutoff, lfc_cutoff)
        return list(degs.loc[degs["log2FoldChange"] < 0].index)

    def summary_counts(
        self,
        padj_cutoff: Optional[float] = None,
        lfc_cutoff: Optional[float] = None,
    ) -> Dict[str, int]:
        """Summary count of UP, DOWN, and Non-Significant genes."""
        up = len(self.get_up_genes(padj_cutoff, lfc_cutoff))
        down = len(self.get_down_genes(padj_cutoff, lfc_cutoff))
        total = len(self.results_df)
        return {
            "total_tested": total,
            "significant_up": up,
            "significant_down": down,
            "total_degs": up + down,
            "not_significant": total - (up + down),
        }

    def summary(self) -> str:
        """Human-readable text summary of DE results."""
        counts = self.summary_counts()
        factor, test, ref = self.contrast
        return (
            f"--- DE Summary: {factor} ({test} vs {ref}) ---\n"
            f"Thresholds: padj <= {self.alpha}, |log2FC| >= {self.lfc_threshold}\n"
            f"Total genes analyzed: {counts['total_tested']}\n"
            f"  - Significant UP:   {counts['significant_up']}\n"
            f"  - Significant DOWN: {counts['significant_down']}\n"
            f"  - Total DEGs:       {counts['total_degs']}\n"
            f"  - Non-significant:  {counts['not_significant']}"
        )


def build_deseq_dataset(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design_factors: Union[str, List[str]] = "condition",
    ref_level: Optional[Dict[str, str]] = None,
    fit_model: bool = True,
    **dds_kwargs,
) -> DeseqDataSet:
    """
    Construct and optionally fit a PyDESeq2 DeseqDataSet from validated inputs.

    Parameters
    ----------
    counts : pd.DataFrame
        Gene count matrix of shape (n_genes, n_samples).
    metadata : pd.DataFrame
        Sample annotation table.
    design_factors : str or list of str, default='condition'
        Factor(s) to include in the negative binomial generalized linear model formula.
    ref_level : dict of {str: str}, optional
        Reference levels for categorical factors (e.g. {'condition': 'control'}).
    fit_model : bool, default=True
        Whether to immediately fit dispersions and LFC (dds.deseq2()).
    **dds_kwargs
        Additional parameters passed to DeseqDataSet (e.g. refit_cooks=True, n_cpus=1).

    Returns
    -------
    DeseqDataSet
        Fitted PyDESeq2 DeseqDataSet object.
    """
    if not HAS_PYDESEQ2:
        raise ImportError("pydeseq2 package is required to run differential expression analysis.")

    # Validate inputs
    val_report = validate_bulk_data(counts, metadata, auto_align=True, strict=True)
    aligned_counts = val_report.aligned_counts
    aligned_metadata = val_report.aligned_metadata.copy()

    # Apply reference levels if specified
    if ref_level is not None:
        for factor_col, base_lvl in ref_level.items():
            if factor_col in aligned_metadata.columns:
                unique_vals = list(aligned_metadata[factor_col].unique())
                if base_lvl in unique_vals:
                    # Move base_lvl to the first position for categorical ordering
                    unique_vals.remove(base_lvl)
                    new_order = [base_lvl] + unique_vals
                    aligned_metadata[factor_col] = pd.Categorical(
                        aligned_metadata[factor_col], categories=new_order, ordered=True
                    )

    logger.info(
        f"Initializing DeseqDataSet with {aligned_counts.shape[0]} genes and "
        f"{aligned_counts.shape[1]} samples. Design: {design_factors}"
    )

    # PyDESeq2 expects counts of shape (n_samples, n_genes)
    dds = DeseqDataSet(
        counts=aligned_counts.T,
        metadata=aligned_metadata,
        design_factors=design_factors,
        **dds_kwargs,
    )

    if fit_model:
        logger.info("Fitting negative binomial GLM, dispersions, and size factors via PyDESeq2...")
        dds.deseq2()
        logger.info("DeseqDataSet model fitting completed.")

    return dds


def run_deseq_stats(
    dds: DeseqDataSet,
    contrast: Tuple[str, str, str],
    alpha: float = 0.05,
    lfc_threshold: float = 1.0,
    **stats_kwargs,
) -> DEResult:
    """
    Execute DeseqStats for a specific contrast and return standardized DEResult.

    Parameters
    ----------
    dds : DeseqDataSet
        Fitted PyDESeq2 DeseqDataSet.
    contrast : tuple of (factor, test_level, ref_level)
        e.g. ('condition', 'treated', 'control')
    alpha : float, default=0.05
        Significance cutoff for padj.
    lfc_threshold : float, default=1.0
        Log2 fold change threshold (absolute) for designating significant DEGs.
    **stats_kwargs
        Additional arguments forwarded to DeseqStats (e.g. cooks_filter=True).

    Returns
    -------
    DEResult
        Standardized result container with formatted DataFrame and original objects.
    """
    if not HAS_PYDESEQ2:
        raise ImportError("pydeseq2 package is required.")

    if not isinstance(contrast, (tuple, list)) or len(contrast) != 3:
        raise ValueError("Contrast must be a 3-element tuple: (factor, test_level, ref_level)")

    factor, test_lvl, ref_lvl = contrast
    logger.info(f"Computing DeseqStats for contrast: {factor} ({test_lvl} vs {ref_lvl})")

    stat_res = DeseqStats(
        dds=dds,
        contrast=list(contrast),
        alpha=alpha,
        **stats_kwargs,
    )
    stat_res.summary()

    raw_df = stat_res.results_df.copy()

    # Standardize column naming and add regulation labels
    res_df = pd.DataFrame(
        {
            "baseMean": raw_df["baseMean"],
            "log2FoldChange": raw_df["log2FoldChange"],
            "lfcSE": raw_df["lfcSE"],
            "stat": raw_df["stat"],
            "pvalue": raw_df["pvalue"],
            "padj": raw_df["padj"],
        },
        index=raw_df.index,
    )
    res_df.index.name = "gene_id"

    # Compute significance and regulation direction
    is_sig = (
        (res_df["padj"].notna())
        & (res_df["padj"] <= alpha)
        & (res_df["log2FoldChange"].abs() >= lfc_threshold)
    )
    res_df["significant"] = is_sig

    conditions = [
        is_sig & (res_df["log2FoldChange"] > 0),
        is_sig & (res_df["log2FoldChange"] < 0),
    ]
    choices = ["UP", "DOWN"]
    res_df["regulation"] = np.select(conditions, choices, default="NS")

    res = DEResult(
        results_df=res_df,
        contrast=(str(factor), str(test_lvl), str(ref_lvl)),
        design_factor=str(factor),
        alpha=alpha,
        lfc_threshold=lfc_threshold,
        dds=dds,
        stat_res=stat_res,
    )

    logger.info(res.summary())
    return res


def run_de(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    contrast: Tuple[str, str, str],
    alpha: float = 0.05,
    lfc_threshold: float = 1.0,
    design_factors: Optional[Union[str, List[str]]] = None,
    **kwargs,
) -> DEResult:
    """
    Convenience wrapper to run complete Differential Expression pipeline in one call:
    Validates input -> Builds & fits DeseqDataSet -> Runs DeseqStats -> Formats DEResult.

    Parameters
    ----------
    counts : pd.DataFrame
        Gene count matrix (n_genes, n_samples).
    metadata : pd.DataFrame
        Sample metadata.
    contrast : tuple of (factor, test_level, ref_level)
        Experimental comparison.
    alpha : float, default=0.05
        FDR significance threshold.
    lfc_threshold : float, default=1.0
        Log2 fold change cutoff.
    design_factors : str or list of str, optional
        Factors for GLM model. Defaults to contrast[0].
    **kwargs
        Additional arguments passed to DeseqDataSet or DeseqStats.

    Returns
    -------
    DEResult
        Complete analysis result.
    """
    factor = design_factors if design_factors is not None else contrast[0]
    ref_level = {contrast[0]: contrast[2]}

    dds = build_deseq_dataset(
        counts=counts,
        metadata=metadata,
        design_factors=factor,
        ref_level=ref_level,
        fit_model=True,
    )

    return run_deseq_stats(
        dds=dds,
        contrast=contrast,
        alpha=alpha,
        lfc_threshold=lfc_threshold,
        **kwargs,
    )
