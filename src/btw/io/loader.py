"""
Data loader module for Bulk Transcriptomics Workbench.
Supports CSV, TSV, Excel, and Parquet file ingestion with format auto-detection.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union

import pandas as pd

from btw import logger
from btw.io.validator import ValidationReport, validate_bulk_data


@dataclass
class BulkDataset:
    """Container holding verified gene count matrix and sample metadata."""

    counts: pd.DataFrame
    metadata: pd.DataFrame
    validation_report: Optional[ValidationReport] = None

    @property
    def n_genes(self) -> int:
        return self.counts.shape[0]

    @property
    def n_samples(self) -> int:
        return self.counts.shape[1]

    @property
    def samples(self) -> list[str]:
        return list(self.counts.columns)

    @property
    def genes(self) -> list[str]:
        return list(self.counts.index)

    def __iter__(self):
        """Allow tuple unpacking: counts, metadata = dataset"""
        return iter((self.counts, self.metadata))


def _detect_separator(filepath: Path) -> str:
    """Heuristically infer delimiter from file extension or head snippet."""
    ext = filepath.suffix.lower()
    if ext in [".tsv", ".tab"]:
        return "\t"
    elif ext == ".csv":
        return ","
    # Peek first line
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            first_line = f.readline()
            if "\t" in first_line:
                return "\t"
            if "," in first_line:
                return ","
            if ";" in first_line:
                return ";"
    except Exception:
        pass
    return ","


def load_counts(
    file_path: Union[str, Path],
    gene_col: Union[str, int] = 0,
    sep: Optional[str] = None,
    **kwargs,
) -> pd.DataFrame:
    """
    Load raw count matrix from CSV, TSV, Excel, or Parquet file.

    Parameters
    ----------
    file_path : str or Path
        Path to the count matrix file.
    gene_col : str or int, default=0
        Column name or 0-based column index to use as gene identifier index.
    sep : str, optional
        Delimiter for text files. If None, inferred automatically.
    **kwargs
        Additional keyword arguments forwarded to pandas reader.

    Returns
    -------
    pd.DataFrame
        Gene expression matrix indexed by gene identifiers.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Count matrix file not found: {path}")

    ext = path.suffix.lower()
    logger.info(f"Loading count matrix from {path} (type: {ext})")

    if ext in [".xlsx", ".xls"]:
        df = pd.read_excel(path, **kwargs)
        if isinstance(gene_col, int):
            col_name = df.columns[gene_col]
            df = df.set_index(col_name)
        else:
            df = df.set_index(gene_col)
    elif ext == ".parquet":
        df = pd.read_parquet(path, **kwargs)
        if isinstance(gene_col, str) and gene_col in df.columns:
            df = df.set_index(gene_col)
    else:
        actual_sep = sep if sep is not None else _detect_separator(path)
        index_col = gene_col if isinstance(gene_col, int) else None
        df = pd.read_csv(path, sep=actual_sep, index_col=index_col, **kwargs)
        if isinstance(gene_col, str) and gene_col in df.columns:
            df = df.set_index(gene_col)

    # Clean index name and ensure string identifiers
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)

    logger.info(
        f"Loaded count matrix with shape {df.shape} ({df.shape[0]} genes, {df.shape[1]} samples)"
    )
    return df


def load_metadata(
    file_path: Union[str, Path],
    sample_col: Optional[Union[str, int]] = 0,
    sep: Optional[str] = None,
    **kwargs,
) -> pd.DataFrame:
    """
    Load sample metadata annotations from CSV, TSV, Excel, or Parquet file.

    Parameters
    ----------
    file_path : str or Path
        Path to the metadata file.
    sample_col : str or int, optional, default=0
        Column name or index to use as sample identifier index.
    sep : str, optional
        Delimiter for text files. If None, inferred automatically.
    **kwargs
        Additional arguments passed to pandas reader.

    Returns
    -------
    pd.DataFrame
        Metadata table indexed by sample identifiers.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Metadata file not found: {path}")

    ext = path.suffix.lower()
    logger.info(f"Loading metadata from {path} (type: {ext})")

    if ext in [".xlsx", ".xls"]:
        df = pd.read_excel(path, **kwargs)
        if sample_col is not None:
            col_name = df.columns[sample_col] if isinstance(sample_col, int) else sample_col
            df = df.set_index(col_name)
    elif ext == ".parquet":
        df = pd.read_parquet(path, **kwargs)
        if sample_col is not None and isinstance(sample_col, str) and sample_col in df.columns:
            df = df.set_index(sample_col)
    else:
        actual_sep = sep if sep is not None else _detect_separator(path)
        index_col = sample_col if isinstance(sample_col, int) else None
        df = pd.read_csv(path, sep=actual_sep, index_col=index_col, **kwargs)
        if sample_col is not None and isinstance(sample_col, str) and sample_col in df.columns:
            df = df.set_index(sample_col)

    df.index = df.index.astype(str)
    logger.info(f"Loaded metadata for {df.shape[0]} samples with columns: {list(df.columns)}")
    return df


def load_dataset(
    counts_path: Union[str, Path],
    metadata_path: Union[str, Path],
    gene_col: Union[str, int] = 0,
    sample_col: Optional[Union[str, int]] = 0,
    validate: bool = True,
    strict: bool = True,
    **kwargs,
) -> BulkDataset:
    """
    Convenience function to load both count matrix and metadata, validate schema,
    and align samples.

    Parameters
    ----------
    counts_path : str or Path
        Path to count matrix.
    metadata_path : str or Path
        Path to sample metadata.
    gene_col : str or int, default=0
        Gene ID column identifier.
    sample_col : str or int, default=0
        Sample ID column identifier.
    validate : bool, default=True
        Whether to run validation checks.
    strict : bool, default=True
        Whether to raise ValidationError if schema checks fail.

    Returns
    -------
    BulkDataset
        Container with aligned counts, metadata, and optional ValidationReport.
    """
    counts = load_counts(counts_path, gene_col=gene_col, **kwargs)
    metadata = load_metadata(metadata_path, sample_col=sample_col, **kwargs)

    report = None
    if validate:
        report = validate_bulk_data(counts, metadata, strict=strict, auto_align=True)
        counts = report.aligned_counts if report.aligned_counts is not None else counts
        metadata = report.aligned_metadata if report.aligned_metadata is not None else metadata

    return BulkDataset(counts=counts, metadata=metadata, validation_report=report)
