"""
Standardized schema and data models for Functional & Pathway Enrichment (FR-5).
Ensures outputs across ORA, GSEA, and decoupler activity inference share a unified tabular structure.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd
from btw import logger


@dataclass
class EnrichmentResult:
    """
    Standardized container holding functional pathway enrichment results
    across ORA (Enrichr, goatools), GSEA (prerank), and Activity Inference (decoupler).
    """
    results_df: pd.DataFrame
    source_method: str               # 'gseapy_enrichr', 'goatools_ora', 'gseapy_prerank', 'decoupler_ulm', etc.
    gene_set_database: str           # 'KEGG_2021_Human', 'MSigDB_Hallmark', 'PROGENy', etc.
    alpha: float = 0.05
    raw_output: Optional[Any] = None # Underlying raw object (e.g. gseapy.Enrichr or decoupler tuple)
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def terms(self) -> List[str]:
        """List of all tested terms/pathways."""
        return list(self.results_df["term"]) if "term" in self.results_df.columns else []

    def get_significant(self, padj_cutoff: Optional[float] = None) -> pd.DataFrame:
        """
        Filter enrichment table for significant terms based on adjusted p-value.

        Parameters
        ----------
        padj_cutoff : float, optional
            P-adjusted threshold. Defaults to self.alpha.

        Returns
        -------
        pd.DataFrame
            Subset of significant enrichment terms.
        """
        cutoff = padj_cutoff if padj_cutoff is not None else self.alpha
        df = self.results_df
        if "padj" not in df.columns:
            return df
        return df.loc[df["padj"].notna() & (df["padj"] <= cutoff)]

    def summary(self) -> str:
        """Summary text of enrichment results."""
        sig_df = self.get_significant()
        total = len(self.results_df)
        sig_cnt = len(sig_df)
        top_terms = list(sig_df["term"][:5]) if not sig_df.empty else ["None"]
        return (
            f"--- Enrichment Summary: {self.source_method} ({self.gene_set_database}) ---\n"
            f"Threshold: padj <= {self.alpha}\n"
            f"Total pathways/terms tested: {total}\n"
            f"Significant terms: {sig_cnt}\n"
            f"Top terms: {', '.join(top_terms)}"
        )


def standardize_enrichment_table(
    df: pd.DataFrame,
    source_method: str,
    term_col: str = "term",
    score_col: str = "score",
    pvalue_col: str = "pvalue",
    padj_col: str = "padj",
    genes_col: Optional[str] = "genes",
    count_col: Optional[str] = "gene_count",
) -> pd.DataFrame:
    """
    Standardize varied third-party enrichment outputs into BTW's uniform schema:
    ['term', 'score', 'pvalue', 'padj', 'gene_count', 'genes', 'source_method']

    Parameters
    ----------
    df : pd.DataFrame
        Input enrichment table from gseapy, goatools, or decoupler.
    source_method : str
        Method identifier string.
    term_col : str, default='term'
        Column for pathway name.
    score_col : str, default='score'
        Column for enrichment score / NES / combined score / activity.
    pvalue_col : str, default='pvalue'
        Column for nominal p-value.
    padj_col : str, default='padj'
        Column for FDR / adjusted p-value.
    genes_col : str, optional
        Column for overlapping / leading-edge genes.
    count_col : str, optional
        Column for gene count.

    Returns
    -------
    pd.DataFrame
        Standardized DataFrame.
    """
    if df.empty:
        return pd.DataFrame(
            columns=["term", "score", "pvalue", "padj", "gene_count", "genes", "source_method"]
        )

    out = pd.DataFrame()
    out["term"] = df[term_col].astype(str) if term_col in df.columns else df.index.astype(str)
    out["score"] = df[score_col].astype(float) if score_col in df.columns else 0.0
    out["pvalue"] = df[pvalue_col].astype(float) if pvalue_col in df.columns else np.nan
    out["padj"] = df[padj_col].astype(float) if padj_col in df.columns else np.nan

    def _parse_count(val):
        if pd.isna(val):
            return 0
        if isinstance(val, (int, np.integer)):
            return int(val)
        if isinstance(val, (float, np.floating)):
            return int(val)
        s = str(val).strip()
        if "/" in s:
            s = s.split("/")[0].strip()
        try:
            return int(float(s))
        except Exception:
            return 0

    if count_col is not None and count_col in df.columns:
        out["gene_count"] = df[count_col].apply(_parse_count)
    else:
        out["gene_count"] = 0

    if genes_col is not None and genes_col in df.columns:
        out["genes"] = df[genes_col]
        # If count was not explicitly provided or resulted in all 0s, count from genes string/list
        if count_col is None or count_col not in df.columns:
            out["gene_count"] = out["genes"].apply(
                lambda x: len(x) if isinstance(x, (list, set)) else len(str(x).split(";")) if str(x) != "nan" and str(x) != "" else 0
            )
    else:
        out["genes"] = ""

    out["source_method"] = source_method

    # Sort by significance (padj ascending, then absolute score descending)
    out = out.sort_values(by=["padj", "pvalue", "score"], ascending=[True, True, False]).reset_index(drop=True)
    return out
