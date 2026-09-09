"""
Functional and Pathway Enrichment subpackage for BTW (FR-5).
Unifies ORA (gseapy/goatools), GSEA (prerank), and Activity Inference (decoupler).
"""

from btw.enrichment.activity import (
    prepare_decoupler_input,
    run_decoupler_activity,
    run_progeny_activity,
    run_tf_activity,
)
from btw.enrichment.cache import get_memory_cache
from btw.enrichment.gsea import prepare_ranked_gene_list, run_prerank
from btw.enrichment.ora import (
    extract_significant_genes,
    run_custom_ora,
    run_enrichr,
)
from btw.enrichment.schema import (
    EnrichmentResult,
    standardize_enrichment_table,
)
from btw.enrichment.viz import (
    plot_enrichment_barplot,
    plot_enrichment_dotplot,
)

__all__ = [
    "EnrichmentResult",
    "standardize_enrichment_table",
    "extract_significant_genes",
    "run_enrichr",
    "run_custom_ora",
    "prepare_ranked_gene_list",
    "run_prerank",
    "prepare_decoupler_input",
    "run_decoupler_activity",
    "run_progeny_activity",
    "run_tf_activity",
    "plot_enrichment_dotplot",
    "plot_enrichment_barplot",
    "get_memory_cache",
]
