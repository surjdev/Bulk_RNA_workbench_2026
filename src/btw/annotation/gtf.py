"""
GTF parsing and reference annotation utilities for BTW (FR-6).
Provides fast loading of gene models and ID mappings from custom GTF/GFF files.
"""

from __future__ import annotations

import io
from pathlib import Path
from typing import Dict, List, Optional, Set, Union

import pandas as pd
from btw import logger

try:
    import gtfparse

    HAS_GTFPARSE = True
except ImportError:
    HAS_GTFPARSE = False


def parse_gtf_file(
    gtf_path_or_buffer: Union[str, Path, io.StringIO],
    features: Optional[Union[List[str], Set[str]]] = None,
    result_type: str = "pandas",
    expand_attributes: bool = True,
    **kwargs,
) -> pd.DataFrame:
    """
    Parse GTF reference file into a standardized pandas DataFrame.

    Parameters
    ----------
    gtf_path_or_buffer : str, Path, or StringIO
        Path to GTF/GFF file (can be gzip compressed) or text buffer.
    features : list or set of str, optional
        Filter features to include (e.g. {'gene', 'transcript', 'exon'}).
        If None, all features are loaded.
    result_type : {'pandas', 'polars'}, default='pandas'
        Target DataFrame type.
    expand_attributes : bool, default=True
        Whether to expand semi-colon separated GTF attributes into columns.
    **kwargs
        Additional arguments passed to gtfparse.read_gtf.

    Returns
    -------
    pd.DataFrame
        Parsed annotations table.
    """
    if not HAS_GTFPARSE:
        raise ImportError("gtfparse is required for GTF parsing. Run pip install gtfparse.")

    feat_set = set(features) if features is not None else None
    logger.info(f"Parsing GTF reference (features={feat_set})...")

    df = gtfparse.read_gtf(
        gtf_path_or_buffer,
        expand_attribute_column=expand_attributes,
        features=feat_set,
        result_type="pandas" if result_type == "pandas" else result_type,
        **kwargs,
    )

    if not isinstance(df, pd.DataFrame) and hasattr(df, "to_pandas"):
        df = df.to_pandas()

    logger.info(f"Successfully loaded {len(df)} genomic records from GTF.")
    return df


def create_gene_map_from_gtf(
    gtf_path_or_buffer: Union[str, Path, io.StringIO],
    from_attr: str = "gene_id",
    to_attr: str = "gene_name",
    feature: str = "gene",
    drop_duplicates: bool = True,
) -> pd.Series:
    """
    Extract a 1-to-1 gene ID to symbol mapping dictionary/Series from a GTF file.

    Parameters
    ----------
    gtf_path_or_buffer : str, Path, or StringIO
        Path to GTF file.
    from_attr : str, default='gene_id'
        Source identifier column.
    to_attr : str, default='gene_name'
        Target symbol column.
    feature : str, default='gene'
        Feature row type to extract from.
    drop_duplicates : bool, default=True
        Whether to drop duplicate source IDs.

    Returns
    -------
    pd.Series
        Mapping with index=from_attr and values=to_attr.
    """
    df = parse_gtf_file(gtf_path_or_buffer, features=[feature])

    if from_attr not in df.columns or to_attr not in df.columns:
        raise ValueError(
            f"Requested attributes ('{from_attr}', '{to_attr}') not found in GTF columns: {list(df.columns)}"
        )

    sub = df[[from_attr, to_attr]].dropna()
    if drop_duplicates:
        sub = sub.drop_duplicates(subset=[from_attr])

    mapping = sub.set_index(from_attr)[to_attr]
    logger.info(f"Constructed gene ID mapping for {len(mapping)} genes ({from_attr} -> {to_attr}).")
    return mapping
