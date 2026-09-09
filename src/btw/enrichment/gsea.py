"""
Gene Set Enrichment Analysis (GSEA) module (FR-5.2).
Prepares ranked gene lists from DE metrics and wraps gseapy.prerank with MSigDB integration.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd
from btw import logger
from btw.de_analysis.deseq_helper import DEResult
from btw.enrichment.schema import EnrichmentResult, standardize_enrichment_table

try:
    import gseapy as gp
    HAS_GSEAPY = True
except ImportError:
    HAS_GSEAPY = False


def prepare_ranked_gene_list(
    de_data: Union[DEResult, pd.DataFrame],
    rank_by: str = "stat",
    drop_na: bool = True,
) -> pd.Series:
    """
    Transform Differential Expression results into a sorted ranked gene list for GSEA.

    Parameters
    ----------
    de_data : DEResult or pd.DataFrame
        Differential expression table containing effect sizes and statistics.
    rank_by : {'stat', 'log2FoldChange', 'signed_pvalue'}, default='stat'
        Ranking metric:
        - 'stat': Wald test statistic from PyDESeq2.
        - 'log2FoldChange': Log2 fold change.
        - 'signed_pvalue': sign(log2FC) * -log10(pvalue).
    drop_na : bool, default=True
        Whether to drop genes with NaN values in the selected ranking metric.

    Returns
    -------
    pd.Series
        Ranked gene series sorted in descending order, indexed by gene symbol/ID.
    """
    if isinstance(de_data, DEResult):
        df = de_data.results_df.copy()
    else:
        df = de_data.copy()

    rank_by_lower = rank_by.lower()

    if rank_by_lower in ["stat", "wald"]:
        if "stat" not in df.columns:
            raise KeyError("Column 'stat' missing from DE data. Use rank_by='log2FoldChange' instead.")
        scores = df["stat"].astype(float)
    elif rank_by_lower in ["log2foldchange", "lfc", "log2fc"]:
        if "log2FoldChange" not in df.columns:
            raise KeyError("Column 'log2FoldChange' missing from DE data.")
        scores = df["log2FoldChange"].astype(float)
    elif rank_by_lower in ["signed_pvalue", "signal"]:
        if "pvalue" not in df.columns or "log2FoldChange" not in df.columns:
            raise KeyError("Columns 'pvalue' and 'log2FoldChange' required for signed_pvalue ranking.")
        pvals = np.clip(df["pvalue"].values.astype(float), 1e-300, 1.0)
        neg_log10 = -np.log10(pvals)
        sign = np.sign(df["log2FoldChange"].values.astype(float))
        scores = pd.Series(sign * neg_log10, index=df.index)
    else:
        raise ValueError(f"Unknown ranking metric '{rank_by}'. Choose from 'stat', 'log2FoldChange', 'signed_pvalue'.")

    scores.name = "ranking_metric"
    scores.index = scores.index.astype(str)

    if drop_na:
        scores = scores.dropna()

    # Deduplicate genes by keeping the maximum absolute score
    if scores.index.duplicated().any():
        scores = scores.groupby(scores.index).apply(lambda s: s.iloc[np.argmax(np.abs(s.values))])

    # Sort descending
    ranked = scores.sort_values(ascending=False)
    logger.info(
        f"Prepared ranked gene list of {len(ranked)} genes (rank_by='{rank_by}', "
        f"range: [{ranked.min():.2f}, {ranked.max():.2f}])."
    )
    return ranked


def run_prerank(
    ranked_genes: Union[pd.Series, DEResult, pd.DataFrame],
    gene_sets: Union[str, List[str], Dict[str, List[str]]] = "MSigDB_Hallmark_2020",
    rank_by: str = "stat",
    min_size: int = 5,
    max_size: int = 500,
    permutation_num: int = 100,
    seed: int = 42,
    alpha: float = 0.05,
    outdir: Optional[str] = None,
    **kwargs,
) -> EnrichmentResult:
    """
    Execute Gene Set Enrichment Analysis (GSEA) using gseapy.prerank.

    Parameters
    ----------
    ranked_genes : Series, DEResult, or pd.DataFrame
        Ranked gene series or DE result to rank automatically.
    gene_sets : str, list of str, or dict
        MSigDB collection (e.g. 'MSigDB_Hallmark_2020', 'KEGG_2021_Human')
        or custom gene set dictionary.
    rank_by : str, default='stat'
        Metric for ranking if DEResult/DataFrame is provided.
    min_size : int, default=5
        Minimum gene set size.
    max_size : int, default=500
        Maximum gene set size.
    permutation_num : int, default=100
        Number of permutations for empirical p-value calculation.
    seed : int, default=42
        Random seed for reproducibility.
    alpha : float, default=0.05
        FDR cutoff.
    outdir : str, optional
        Output folder for GSEA plots.
    **kwargs
        Additional arguments forwarded to gseapy.prerank.

    Returns
    -------
    EnrichmentResult
        Standardized GSEA enrichment results.
    """
    if not HAS_GSEAPY:
        raise ImportError("gseapy is required for GSEA prerank. Run pip install gseapy.")

    # 1. Prepare ranked list if needed
    if isinstance(ranked_genes, pd.Series):
        rnk = ranked_genes.copy()
    else:
        rnk = prepare_ranked_gene_list(ranked_genes, rank_by=rank_by)

    logger.info(f"Running gseapy.prerank with {len(rnk)} ranked genes against {gene_sets}...")

    prerank_res = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        min_size=min_size,
        max_size=max_size,
        permutation_num=permutation_num,
        seed=seed,
        outdir=outdir,
        no_plot=True,
        verbose=False,
        **kwargs,
    )

    raw_df = prerank_res.res2d.copy()
    db_name = str(gene_sets) if isinstance(gene_sets, (str, dict)) else ";".join(gene_sets)

    # Standardize output
    std_df = standardize_enrichment_table(
        raw_df,
        source_method="gseapy_prerank",
        term_col="Term",
        score_col="NES",
        pvalue_col="NOM p-val",
        padj_col="FDR q-val",
        genes_col="Lead_genes",
        count_col="Tag %",
    )

    # Map Tag % or Matched size to gene_count if available
    if "Matched size" in raw_df.columns:
        std_df["gene_count"] = raw_df["Matched size"].astype(int).values

    return EnrichmentResult(
        results_df=std_df,
        source_method="gseapy_prerank",
        gene_set_database=db_name,
        alpha=alpha,
        raw_output=prerank_res,
        metadata={"rank_by": rank_by, "permutation_num": permutation_num},
    )
