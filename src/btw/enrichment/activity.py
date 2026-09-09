"""
Activity and Pathway Inference module via decoupler (FR-5.3).
Computes pathway activities (PROGENy) and transcription factor activities (CollecTRI/DoRothEA) from DE statistics.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from btw import logger
from btw.de_analysis.correction import adjust_pvalues
from btw.de_analysis.deseq_helper import DEResult
from btw.enrichment.schema import EnrichmentResult, standardize_enrichment_table

try:
    import decoupler as dc
    HAS_DECOUPLER = True
except ImportError:
    HAS_DECOUPLER = False


def prepare_decoupler_input(
    de_data: Union[DEResult, pd.DataFrame],
    metric: str = "stat",
    contrast_name: Optional[str] = None,
) -> pd.DataFrame:
    """
    Format differential expression output into the matrix structure required by decoupler
    (shape: 1 row x n_genes).

    Parameters
    ----------
    de_data : DEResult or pd.DataFrame
        Differential expression table containing gene statistics.
    metric : {'stat', 'log2FoldChange'}, default='stat'
        Metric to infer pathway activities from.
    contrast_name : str, optional
        Row label representing the contrast.

    Returns
    -------
    pd.DataFrame
        DataFrame of shape (1, n_genes) with genes as columns and contrast as index.
    """
    if isinstance(de_data, DEResult):
        df = de_data.results_df.copy()
        row_id = contrast_name if contrast_name is not None else "_vs_".join(de_data.contrast[1:])
    else:
        df = de_data.copy()
        row_id = contrast_name if contrast_name is not None else "contrast"

    metric_col = "stat" if metric == "stat" else "log2FoldChange"
    if metric_col not in df.columns:
        metric_col = "log2FoldChange" if "log2FoldChange" in df.columns else df.columns[0]

    series = df[metric_col].dropna()
    mat = pd.DataFrame([series.values], index=[row_id], columns=series.index)
    mat.columns = mat.columns.astype(str)

    logger.info(f"Prepared decoupler input matrix: {mat.shape[0]} contrast x {mat.shape[1]} genes (metric='{metric_col}').")
    return mat


def run_decoupler_activity(
    mat: pd.DataFrame,
    net: pd.DataFrame,
    method: str = "ulm",
    min_n: int = 3,
    alpha: float = 0.05,
    network_name: str = "Prior_Network",
    **method_kwargs,
) -> EnrichmentResult:
    """
    Run biological activity inference on a contrast statistic matrix against a prior network.

    Parameters
    ----------
    mat : pd.DataFrame
        Matrix of shape (n_contrasts, n_genes).
    net : pd.DataFrame
        Network table containing ['source', 'target', 'weight'].
    method : {'ulm', 'mlm', 'ora'}, default='ulm'
        Statistical inference method from decoupler.
    min_n : int, default=3
        Minimum targets per regulon/pathway.
    alpha : float, default=0.05
        FDR cutoff.
    network_name : str, default='Prior_Network'
        Name of network database.
    **method_kwargs
        Additional arguments passed to decoupler method.

    Returns
    -------
    EnrichmentResult
        Standardized enrichment table containing activity scores and p-values.
    """
    if not HAS_DECOUPLER:
        raise ImportError("decoupler package is required for activity inference. Run pip install decoupler.")

    # Validate network columns
    required_cols = {"source", "target", "weight"}
    if not required_cols.issubset(net.columns):
        raise ValueError(f"Network DataFrame must contain columns: {required_cols}. Found: {list(net.columns)}")

    method_lower = method.lower()
    logger.info(f"Running decoupler method='{method_lower}' on {mat.shape[1]} genes against {len(net)} network edges...")

    if hasattr(dc, "mt"):
        # decoupler >= 2.0
        method_func = getattr(dc.mt, method_lower, None)
        if method_func is None:
            raise ValueError(f"Unsupported decoupler method '{method}'. Valid: 'ulm', 'mlm', 'ora'.")
        acts, pvals = method_func(mat, net, tmin=min_n, verbose=False, **method_kwargs)
    elif hasattr(dc, f"run_{method_lower}"):
        # decoupler < 2.0
        run_func = getattr(dc, f"run_{method_lower}")
        acts, pvals = run_func(mat=mat, net=net, min_n=min_n, verbose=False, **method_kwargs)
    else:
        raise ValueError(f"Unsupported decoupler method '{method}'. Valid: 'ulm', 'mlm', 'ora'.")

    # Format first row (single contrast) into summary table
    row_idx = acts.index[0]
    scores_s = acts.loc[row_idx]
    pvals_s = pvals.loc[row_idx]

    res_df = pd.DataFrame(
        {
            "term": scores_s.index,
            "score": scores_s.values,
            "pvalue": pvals_s.values,
        }
    )

    # Multiple testing correction
    _, padj = adjust_pvalues(res_df["pvalue"].values, method="fdr_bh", alpha=alpha)
    res_df["padj"] = padj

    std_df = standardize_enrichment_table(
        res_df,
        source_method=f"decoupler_{method_lower}",
        term_col="term",
        score_col="score",
        pvalue_col="pvalue",
        padj_col="padj",
    )

    return EnrichmentResult(
        results_df=std_df,
        source_method=f"decoupler_{method_lower}",
        gene_set_database=network_name,
        alpha=alpha,
        raw_output=(acts, pvals),
        metadata={"method": method_lower, "contrast": str(row_idx)},
    )


def run_progeny_activity(
    de_data: Union[DEResult, pd.DataFrame],
    net: Optional[pd.DataFrame] = None,
    organism: str = "human",
    top: int = 100,
    metric: str = "stat",
    method: str = "ulm",
    **kwargs,
) -> EnrichmentResult:
    """
    Infer signaling pathway activity scores (e.g. MAPK, PI3K, TNFa) using PROGENy network weights.

    Parameters
    ----------
    de_data : DEResult or pd.DataFrame
        Differential expression data.
    net : pd.DataFrame, optional
        Pre-loaded PROGENy network. If None, retrieved via decoupler.
    organism : str, default='human'
        Target species.
    top : int, default=100
        Number of top footprint genes per pathway.
    metric : str, default='stat'
        DE metric to infer from.
    method : str, default='ulm'
        Decoupler inference algorithm.

    Returns
    -------
    EnrichmentResult
        Standardized pathway activity scores.
    """
    if net is None:
        if not HAS_DECOUPLER:
            raise ImportError("decoupler is required.")
        logger.info(f"Retrieving PROGENy network for {organism} (top={top})...")
        if hasattr(dc, "op") and hasattr(dc.op, "progeny"):
            net = dc.op.progeny(organism=organism, top=top)
        elif hasattr(dc, "get_progeny"):
            net = dc.get_progeny(organism=organism, top=top)
        else:
            raise AttributeError("decoupler has no progeny network retrieval function.")

    mat = prepare_decoupler_input(de_data, metric=metric)
    return run_decoupler_activity(
        mat=mat,
        net=net,
        method=method,
        network_name="PROGENy_Pathways",
        **kwargs,
    )


def run_tf_activity(
    de_data: Union[DEResult, pd.DataFrame],
    net: Optional[pd.DataFrame] = None,
    organism: str = "human",
    metric: str = "stat",
    method: str = "ulm",
    **kwargs,
) -> EnrichmentResult:
    """
    Infer Transcription Factor (TF) regulon activities using CollecTRI / DoRothEA network weights.

    Parameters
    ----------
    de_data : DEResult or pd.DataFrame
        Differential expression statistics.
    net : pd.DataFrame, optional
        Pre-loaded regulon network. If None, retrieved via decoupler.
    organism : str, default='human'
        Organism name.
    metric : str, default='stat'
        Input ranking metric.
    method : str, default='ulm'
        Inference algorithm.

    Returns
    -------
    EnrichmentResult
        Standardized TF activity results.
    """
    if net is None:
        if not HAS_DECOUPLER:
            raise ImportError("decoupler is required.")
        logger.info(f"Retrieving CollecTRI gene regulatory network for {organism}...")
        if hasattr(dc, "op") and hasattr(dc.op, "collectri"):
            net = dc.op.collectri(organism=organism)
        elif hasattr(dc, "get_collectri"):
            net = dc.get_collectri(organism=organism)
        else:
            raise AttributeError("decoupler has no collectri network retrieval function.")

    mat = prepare_decoupler_input(de_data, metric=metric)
    return run_decoupler_activity(
        mat=mat,
        net=net,
        method=method,
        network_name="CollecTRI_TFs",
        **kwargs,
    )
