"""
Gene Annotation and Reference Mapping module (FR-6).
"""

from btw.annotation.gtf import create_gene_map_from_gtf, parse_gtf_file
from btw.annotation.mapper import annotate_de_results, map_gene_ids

__all__ = [
    "map_gene_ids",
    "annotate_de_results",
    "parse_gtf_file",
    "create_gene_map_from_gtf",
]
