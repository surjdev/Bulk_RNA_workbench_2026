"""
Over-Representation Analysis (ORA) module (FR-5.1).
Integrates gseapy.enrichr and goatools with custom background gene sets and offline statistical tests.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Union

import pandas as pd
from scipy.stats import fisher_exact

from btw import logger
from btw.de_analysis.correction import adjust_pvalues
from btw.de_analysis.deseq_helper import DEResult
from btw.enrichment.schema import EnrichmentResult, standardize_enrichment_table

try:
    import gseapy as gp

    HAS_GSEAPY = True
except ImportError:
    HAS_GSEAPY = False

try:
    import goatools  # noqa: F401

    HAS_GOATOOLS = True
except ImportError:
    HAS_GOATOOLS = False


def extract_significant_genes(
    de_data: Union[DEResult, pd.DataFrame],
    padj_cutoff: float = 0.05,
    lfc_cutoff: float = 1.0,
    direction: str = "both",
) -> List[str]:
    """
    Extract significant gene identifiers from DE analysis results.

    Parameters
    ----------
    de_data : DEResult or pd.DataFrame
        Differential expression results.
    padj_cutoff : float, default=0.05
        Maximum adjusted p-value.
    lfc_cutoff : float, default=1.0
        Minimum absolute log2 fold change.
    direction : {'both', 'up', 'down'}, default='both'
        Filter by regulation direction.

    Returns
    -------
    list of str
        Identifiers of significant genes.
    """
    if isinstance(de_data, DEResult):
        if direction == "up":
            return de_data.get_up_genes(padj_cutoff, lfc_cutoff)
        elif direction == "down":
            return de_data.get_down_genes(padj_cutoff, lfc_cutoff)
        else:
            return list(de_data.get_degs(padj_cutoff, lfc_cutoff).index)

    df = de_data
    is_sig = (
        (df["padj"].notna())
        & (df["padj"] <= padj_cutoff)
        & (df["log2FoldChange"].abs() >= lfc_cutoff)
    )
    if direction == "up":
        mask = is_sig & (df["log2FoldChange"] > 0)
    elif direction == "down":
        mask = is_sig & (df["log2FoldChange"] < 0)
    else:
        mask = is_sig

    return list(df.loc[mask].index)


def run_custom_ora(
    gene_list: List[str],
    gene_sets: Dict[str, List[str]],
    background: Optional[List[str]] = None,
    alpha: float = 0.05,
) -> EnrichmentResult:
    """
    Run exact Fisher's exact test for Over-Representation Analysis using custom gene sets.
    Operates completely offline without requiring network access.

    Parameters
    ----------
    gene_list : list of str
        Significant query genes.
    gene_sets : dict of {str: list of str}
        Mapping of pathway/term names to member genes.
    background : list of str, optional
        Universe of tested genes. If None, union of gene_sets + gene_list is used.
    alpha : float, default=0.05
        FDR significance cutoff.

    Returns
    -------
    EnrichmentResult
        Standardized enrichment results.
    """
    query_set = set(gene_list)

    if background is None:
        bg_set = set().union(*gene_sets.values()).union(query_set)
    else:
        bg_set = set(background)

    # Filter query to background
    query_in_bg = query_set & bg_set
    n_query = len(query_in_bg)
    n_bg = len(bg_set)

    rows = []
    for term, term_genes in gene_sets.items():
        term_set = set(term_genes) & bg_set
        n_term = len(term_set)

        if n_term == 0:
            continue

        overlap = query_in_bg & term_set
        k = len(overlap)

        # 2x2 Contingency table:
        # [[overlap, query_not_in_term],
        #  [term_not_in_query, neither]]
        a = k
        b = n_query - k
        c = n_term - k
        d = (n_bg - n_query) - (n_term - k)

        if d < 0:
            d = 0

        table = [[a, b], [c, d]]
        odds_ratio, pval = fisher_exact(table, alternative="greater")

        rows.append(
            {
                "term": term,
                "score": odds_ratio,
                "pvalue": pval,
                "gene_count": k,
                "genes": ";".join(sorted(list(overlap))),
            }
        )

    res_df = pd.DataFrame(rows)
    if not res_df.empty:
        _, padj = adjust_pvalues(res_df["pvalue"].values, method="fdr_bh", alpha=alpha)
        res_df["padj"] = padj
    else:
        res_df["padj"] = []

    std_df = standardize_enrichment_table(
        res_df,
        source_method="custom_ora_fisher",
        term_col="term",
        score_col="score",
        pvalue_col="pvalue",
        padj_col="padj",
        genes_col="genes",
        count_col="gene_count",
    )

    return EnrichmentResult(
        results_df=std_df,
        source_method="custom_ora_fisher",
        gene_set_database="custom",
        alpha=alpha,
        metadata={"n_query": n_query, "n_background": n_bg},
    )


def run_enrichr(
    gene_list: Union[List[str], DEResult, pd.DataFrame],
    gene_sets: Union[str, List[str], Dict[str, List[str]]] = "KEGG_2021_Human",
    background: Optional[Union[List[str], int]] = None,
    organism: str = "human",
    padj_cutoff: float = 0.05,
    lfc_cutoff: float = 1.0,
    direction: str = "both",
    outdir: Optional[str] = None,
    **kwargs,
) -> EnrichmentResult:
    """
    Execute Over-Representation Analysis (ORA) via Enrichr API or custom dictionary.

    Parameters
    ----------
    gene_list : list of str, DEResult, or pd.DataFrame
        Query gene list or DE result (genes will be extracted automatically).
    gene_sets : str, list of str, or dict
        Enrichr library name (e.g. 'KEGG_2021_Human', 'MSigDB_Hallmark_2020')
        or dictionary of custom gene sets.
    background : list of str or int, optional
        Custom background gene universe.
    organism : str, default='human'
        Target organism.
    padj_cutoff : float, default=0.05
        Cutoff when extracting genes from DE result.
    lfc_cutoff : float, default=1.0
        Log2FC cutoff when extracting genes from DE result.
    direction : {'both', 'up', 'down'}, default='both'
        Regulation filter when extracting genes.
    outdir : str, optional
        Output directory for Enrichr artifacts.
    **kwargs
        Additional arguments passed to gseapy.enrichr.

    Returns
    -------
    EnrichmentResult
        Standardized enrichment result.
    """
    # 1. Resolve query gene list
    if isinstance(gene_list, (DEResult, pd.DataFrame)):
        sig_genes = extract_significant_genes(
            gene_list, padj_cutoff=padj_cutoff, lfc_cutoff=lfc_cutoff, direction=direction
        )
        logger.info(f"Extracted {len(sig_genes)} significant genes ({direction}) from DE results.")
    else:
        sig_genes = list(gene_list)

    if not sig_genes:
        logger.warning("Empty gene list provided for ORA enrichment analysis.")
        empty_std = standardize_enrichment_table(pd.DataFrame(), source_method="gseapy_enrichr")
        return EnrichmentResult(
            results_df=empty_std,
            source_method="gseapy_enrichr",
            gene_set_database=str(gene_sets),
            alpha=padj_cutoff,
        )

    # 2. If gene_sets is a dictionary, use custom Fisher exact test
    if isinstance(gene_sets, dict):
        bg_list = list(background) if isinstance(background, (list, set)) else None
        return run_custom_ora(sig_genes, gene_sets, background=bg_list, alpha=padj_cutoff)

    # 3. Use gseapy.enrichr
    if not HAS_GSEAPY:
        raise ImportError("gseapy is required for Enrichr integration. Run pip install gseapy.")

    logger.info(f"Running gseapy.enrichr with {len(sig_genes)} genes against {gene_sets}...")
    try:
        enr = gp.enrichr(
            gene_list=sig_genes,
            gene_sets=gene_sets,
            organism=organism,
            background=background,
            outdir=outdir,
            no_plot=True,
            **kwargs,
        )
        raw_df = enr.results.copy()
        db_name = str(gene_sets) if isinstance(gene_sets, str) else ";".join(gene_sets)

        std_df = standardize_enrichment_table(
            raw_df,
            source_method="gseapy_enrichr",
            term_col="Term",
            score_col="Combined Score",
            pvalue_col="P-value",
            padj_col="Adjusted P-value",
            genes_col="Genes",
            count_col="Overlap",
        )

        return EnrichmentResult(
            results_df=std_df,
            source_method="gseapy_enrichr",
            gene_set_database=db_name,
            alpha=padj_cutoff,
            raw_output=enr,
        )
    except Exception as e:
        logger.warning(
            f"gseapy.enrichr encountered error ({e}); check network connectivity or parameters."
        )
        raise


def run_clusterprofiler(
    gene_list: Union[List[str], DEResult, pd.DataFrame],
    org_db: str = "org.Hs.eg.db",
    key_type: str = "SYMBOL",
    ont: str = "BP",
    pvalue_cutoff: float = 0.05,
    qvalue_cutoff: float = 0.2,
    background: Optional[List[str]] = None,
    padj_cutoff: float = 0.05,
    lfc_cutoff: float = 1.0,
    direction: str = "both",
    fallback_to_python: bool = True,
    fallback_gene_sets: Union[str, Dict[str, List[str]]] = "KEGG_2021_Human",
    **kwargs,
) -> EnrichmentResult:
    """
    Execute Over-Representation Analysis (ORA) using R Bioconductor clusterProfiler::enrichGO (FR-5 & FR-10).
    If clusterProfiler or org_db is unavailable and fallback_to_python=True, gracefully falls back to Python ORA.

    Parameters
    ----------
    gene_list : list of str, DEResult, or pd.DataFrame
        Significant gene identifiers or DE analysis result object.
    org_db : str, default='org.Hs.eg.db'
        Bioconductor organism annotation database package name.
    key_type : str, default='SYMBOL'
        Identifier type ('SYMBOL', 'ENSEMBL', 'ENTREZID').
    ont : str, default='BP'
        Gene ontology sub-ontology: 'BP' (Biological Process), 'MF' (Molecular Function), 'CC' (Cellular Component).
    pvalue_cutoff : float, default=0.05
        P-value threshold.
    qvalue_cutoff : float, default=0.2
        Q-value threshold.
    background : list of str, optional
        Custom universe of background genes.
    padj_cutoff : float, default=0.05
        Cutoff when extracting genes from DE result.
    lfc_cutoff : float, default=1.0
        Log2FC cutoff when extracting genes from DE result.
    direction : {'both', 'up', 'down'}, default='both'
        Regulation filter when extracting genes.
    fallback_to_python : bool, default=True
        Whether to fall back to Python Enrichr / Fisher exact test if clusterProfiler or R is unavailable.
    fallback_gene_sets : str or dict, default='KEGG_2021_Human'
        Gene sets to use if falling back to Python.
    **kwargs
        Additional arguments forwarded to fallback function or clusterProfiler.

    Returns
    -------
    EnrichmentResult
        Standardized enrichment container.
    """
    if isinstance(gene_list, (DEResult, pd.DataFrame)):
        sig_genes = extract_significant_genes(
            gene_list, padj_cutoff=padj_cutoff, lfc_cutoff=lfc_cutoff, direction=direction
        )
        logger.info(
            f"Extracted {len(sig_genes)} significant genes ({direction}) from DE results for clusterProfiler."
        )
    else:
        sig_genes = list(gene_list)

    if not sig_genes:
        logger.warning("Empty gene list provided for clusterProfiler ORA analysis.")
        empty_std = standardize_enrichment_table(pd.DataFrame(), source_method="R_clusterProfiler")
        return EnrichmentResult(
            results_df=empty_std,
            source_method="R_clusterProfiler",
            gene_set_database=f"{org_db}_{ont}",
            alpha=pvalue_cutoff,
        )

    try:
        from btw.r_interop.bridge import check_r_package, is_r_available, require_r_package

        if is_r_available() and check_r_package("clusterProfiler") and check_r_package(org_db):
            from btw.r_interop.clusterprofiler_r import run_r_clusterprofiler

            return run_r_clusterprofiler(
                gene_list=sig_genes,
                org_db=org_db,
                key_type=key_type,
                ont=ont,
                pvalue_cutoff=pvalue_cutoff,
                qvalue_cutoff=qvalue_cutoff,
                background=background,
            )
        else:
            if not fallback_to_python:
                require_r_package(
                    "clusterProfiler", purpose="Bioconductor GO/KEGG Enrichment Analysis (FR-5)"
                )
            logger.info(
                f"R package 'clusterProfiler' or '{org_db}' not available; falling back to Python ORA."
            )
    except Exception as e:
        if not fallback_to_python:
            raise e
        logger.info(f"clusterProfiler execution failed ({e}); falling back to Python ORA.")

    # Fallback to Python
    if isinstance(fallback_gene_sets, dict):
        bg_list = list(background) if isinstance(background, (list, set)) else None
        return run_custom_ora(
            sig_genes, fallback_gene_sets, background=bg_list, alpha=pvalue_cutoff
        )
    else:
        return run_enrichr(
            sig_genes,
            gene_sets=fallback_gene_sets,
            background=background,
            padj_cutoff=padj_cutoff,
            lfc_cutoff=lfc_cutoff,
            direction=direction,
            **kwargs,
        )
