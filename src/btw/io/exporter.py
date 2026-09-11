"""
Data export module for Bulk Transcriptomics Workbench.
Supports exporting tabular results, DE tables, and matrices to CSV, TSV, Excel, and Parquet.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Union

import pandas as pd

from btw import logger


def export_table(
    df: pd.DataFrame,
    output_path: Union[str, Path],
    index: bool = True,
    sep: Optional[str] = None,
    **kwargs,
) -> Path:
    """
    Export DataFrame to disk in CSV, TSV, Excel, or Parquet format based on extension.

    Parameters
    ----------
    df : pd.DataFrame
        Data table to export.
    output_path : str or Path
        Destination file path.
    index : bool, default=True
        Whether to write row index.
    sep : str, optional
        Delimiter for text formats.
    **kwargs
        Additional arguments passed to specific pandas exporter.

    Returns
    -------
    Path
        Absolute path to the exported file.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ext = path.suffix.lower()

    if ext in [".xlsx", ".xls"]:
        df.to_excel(path, index=index, engine="openpyxl", **kwargs)
    elif ext == ".parquet":
        df.to_parquet(path, index=index, **kwargs)
    elif ext in [".tsv", ".tab"]:
        delimiter = sep if sep is not None else "\t"
        df.to_csv(path, sep=delimiter, index=index, **kwargs)
    else:  # Default to CSV
        delimiter = sep if sep is not None else ","
        df.to_csv(path, sep=delimiter, index=index, **kwargs)

    logger.info(f"Exported table ({df.shape[0]} rows, {df.shape[1]} cols) to {path.resolve()}")
    return path.resolve()


def export_excel_multisheet(
    sheets: Dict[str, pd.DataFrame],
    output_path: Union[str, Path],
    index: bool = True,
    **kwargs,
) -> Path:
    """
    Export multiple DataFrames into separate worksheets in a single Excel workbook.
    Useful for multi-contrast DE summaries and combined QC/DE reports.

    Parameters
    ----------
    sheets : dict of {str: pd.DataFrame}
        Mapping of sheet name to DataFrame.
    output_path : str or Path
        Destination Excel workbook path (.xlsx).
    index : bool, default=True
        Whether to write row names as index.
    **kwargs
        Additional options passed to to_excel.

    Returns
    -------
    Path
        Absolute path to created workbook.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(path, engine="openpyxl", **kwargs) as writer:
        for sheet_name, data in sheets.items():
            # Truncate sheet name to Excel's 31-character limit
            safe_name = sheet_name[:31]
            data.to_excel(writer, sheet_name=safe_name, index=index)

    logger.info(f"Exported {len(sheets)} sheets to Excel workbook: {path.resolve()}")
    return path.resolve()
