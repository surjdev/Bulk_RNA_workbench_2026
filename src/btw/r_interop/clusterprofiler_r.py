"""
R clusterProfiler wrapper through rpy2 for BTW (FR-5 & FR-10).
Serves as the optional Bioconductor reference implementation for Gene Ontology / KEGG over-representation analysis.
"""

from __future__ import annotations

from typing import List, Optional

import numpy as np
import pandas as pd

from btw import logger
from btw.enrichment.schema import EnrichmentResult, standardize_enrichment_table
from btw.r_interop.bridge import (
    get_r_package_version,
    is_r_available,
    r_to_pandas_df,
    require_r_package,
)

if is_r_available():
    import rpy2.robjects as ro

    _r = ro.r
else:
    _r = None


def run_r_clusterprofiler(
    gene_list: List[str],
    org_db: str = "org.Hs.eg.db",
    key_type: str = "SYMBOL",
    ont: str = "BP",
    pvalue_cutoff: float = 0.05,
    qvalue_cutoff: float = 0.2,
    background: Optional[List[str]] = None,
) -> EnrichmentResult:
    """
    Execute R Bioconductor clusterProfiler::enrichGO over-representation analysis via rpy2.

    Parameters
    ----------
    gene_list : list of str
        List of significant query gene identifiers.
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

    Returns
    -------
    EnrichmentResult
        Standardized enrichment container adhering to the BTW unified schema.
    """
    require_r_package(
        "clusterProfiler",
        purpose="Bioconductor GO/KEGG Enrichment Analysis (FR-5)",
        alternative="Use run_enrichr() or run_custom_ora() in Python.",
    )
    require_r_package(
        org_db,
        purpose=f"Organism annotation package {org_db}",
        alternative="Install via BiocManager::install('{org_db}')",
    )

    logger.info(
        f"Executing R clusterProfiler::enrichGO (v{get_r_package_version('clusterProfiler')}, ont={ont}) on {len(gene_list)} genes..."
    )

    r_genes = ro.StrVector([str(g) for g in gene_list])
    _r.assign(".cp_genes", r_genes)

    universe_cmd = ""
    if background:
        _r.assign(".cp_universe", ro.StrVector([str(g) for g in background]))
        universe_cmd = ", universe=.cp_universe"

    cmd = (
        f"clusterProfiler::enrichGO("
        f"  gene=.cp_genes, "
        f"  OrgDb='{org_db}', "
        f"  keyType='{key_type}', "
        f"  ont='{ont}', "
        f"  pvalueCutoff={pvalue_cutoff}, "
        f"  qvalueCutoff={qvalue_cutoff}"
        f"  {universe_cmd}"
        f")"
    )
    r_res = _r(cmd)
    r_df = _r("as.data.frame")(r_res)
    df = r_to_pandas_df(r_df)

    if df.empty:
        std_df = pd.DataFrame(
            columns=["term", "score", "pvalue", "padj", "gene_count", "genes", "source_method"]
        )
    else:
        # clusterProfiler columns: ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count
        # Map: Description -> term, p.adjust -> padj, Count -> gene_count, geneID -> genes, -log10(p.adjust) -> score
        df["score"] = -np.log10(df["p.adjust"].clip(lower=1e-300))
        std_df = standardize_enrichment_table(
            df=df,
            source_method=f"clusterProfiler_{ont}",
            term_col="Description",
            score_col="score",
            pvalue_col="pvalue",
            padj_col="p.adjust",
            genes_col="geneID",
            count_col="Count",
        )

    return EnrichmentResult(
        results_df=std_df,
        source_method=f"R_clusterProfiler_{ont}",
        gene_set_database=f"{org_db}_{ont}",
        alpha=pvalue_cutoff,
        raw_output=r_res,
        metadata={"org_db": org_db, "key_type": key_type, "ont": ont},
    )
