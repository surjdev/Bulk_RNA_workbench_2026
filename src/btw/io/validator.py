"""
Data validation module for Bulk Transcriptomics Workbench.
Validates count matrices and sample metadata schemas according to FR-1.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np
import pandas as pd

from btw import logger


class ValidationError(ValueError):
    """Raised when bulk RNA-seq input data violates required schema or integrity constraints."""

    pass


@dataclass
class ValidationReport:
    """Summary report of data integrity and schema validation checks."""

    is_valid: bool
    n_genes: int
    n_samples: int
    sample_ids: List[str]
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    aligned_counts: Optional[pd.DataFrame] = None
    aligned_metadata: Optional[pd.DataFrame] = None

    def summary(self) -> str:
        status = "PASSED" if self.is_valid else "FAILED"
        msg = [
            f"--- Validation Report [{status}] ---",
            f"Genes: {self.n_genes}, Samples: {self.n_samples}",
        ]
        if self.errors:
            msg.append(f"Errors ({len(self.errors)}):")
            for err in self.errors:
                msg.append(f"  - [ERROR] {err}")
        if self.warnings:
            msg.append(f"Warnings ({len(self.warnings)}):")
            for warn in self.warnings:
                msg.append(f"  - [WARN]  {warn}")
        return "\n".join(msg)


def validate_bulk_data(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    sample_id_col: Optional[str] = None,
    allow_float: bool = False,
    auto_align: bool = True,
    strict: bool = False,
) -> ValidationReport:
    """
    Validate count matrix and sample metadata integrity and consistency.

    Parameters
    ----------
    counts : pd.DataFrame
        Gene expression matrix of shape (n_genes, n_samples), where rows are genes
        and columns are sample identifiers.
    metadata : pd.DataFrame
        Sample annotation table. Sample identifiers can either be the DataFrame index
        or located in a designated column specified by `sample_id_col`.
    sample_id_col : str, optional
        Name of column containing sample IDs in metadata if not the index.
    allow_float : bool, default=False
        If False, raises warning/error if count values are non-integer floats.
    auto_align : bool, default=True
        If True, reorders metadata rows to match count matrix columns exactly.
    strict : bool, default=False
        If True, raises ValidationError on any error; otherwise returns ValidationReport.

    Returns
    -------
    ValidationReport
        Comprehensive report containing validation status, errors, warnings,
        and aligned DataFrames.
    """
    errors: List[str] = []
    warnings: List[str] = []

    # 1. Type and dimension checks
    if not isinstance(counts, pd.DataFrame):
        errors.append(f"`counts` must be a pandas DataFrame, got {type(counts).__name__}")
        return ValidationReport(
            is_valid=False, n_genes=0, n_samples=0, sample_ids=[], errors=errors
        )

    if not isinstance(metadata, pd.DataFrame):
        errors.append(f"`metadata` must be a pandas DataFrame, got {type(metadata).__name__}")
        return ValidationReport(
            is_valid=False,
            n_genes=counts.shape[0],
            n_samples=counts.shape[1],
            sample_ids=[],
            errors=errors,
        )

    if counts.empty:
        errors.append("Count matrix is empty.")
    if metadata.empty:
        errors.append("Metadata table is empty.")

    if errors:
        if strict:
            raise ValidationError("\n".join(errors))
        return ValidationReport(
            is_valid=False,
            n_genes=counts.shape[0],
            n_samples=counts.shape[1],
            sample_ids=list(counts.columns),
            errors=errors,
        )

    # 2. Duplicate gene identifiers check
    if counts.index.duplicated().any():
        dup_count = counts.index.duplicated().sum()
        sample_dups = list(counts.index[counts.index.duplicated()][:5])
        errors.append(
            f"Duplicate gene identifiers found in counts index ({dup_count} duplicates, e.g. {sample_dups}). "
            "Gene identifiers must be unique."
        )

    # 3. Resolve metadata sample IDs
    meta_df = metadata.copy()
    if sample_id_col is not None:
        if sample_id_col not in meta_df.columns:
            errors.append(
                f"Specified `sample_id_col` '{sample_id_col}' not found in metadata columns."
            )
        else:
            meta_df = meta_df.set_index(sample_id_col)

    if meta_df.index.duplicated().any():
        dup_samples = list(meta_df.index[meta_df.index.duplicated()][:5])
        errors.append(f"Duplicate sample IDs found in metadata index (e.g. {dup_samples}).")

    # Normalize column/index names to strings for reliable matching
    counts_norm = counts.copy()
    counts_norm.columns = [str(c) for c in counts_norm.columns]
    meta_df.index = [str(i) for i in meta_df.index]

    # 4. Sample correspondence check
    count_samples_set = set(counts_norm.columns)
    meta_samples_set = set(meta_df.index)

    missing_in_meta = count_samples_set - meta_samples_set
    missing_in_counts = meta_samples_set - count_samples_set

    if missing_in_meta:
        errors.append(
            f"Samples present in counts but missing in metadata ({len(missing_in_meta)}): "
            f"{sorted(list(missing_in_meta))[:5]}"
        )
    if missing_in_counts:
        errors.append(
            f"Samples present in metadata but missing in counts ({len(missing_in_counts)}): "
            f"{sorted(list(missing_in_counts))[:5]}"
        )

    # 5. Missing / NaN values check in counts
    nan_count = counts_norm.isna().sum().sum()
    if nan_count > 0:
        errors.append(
            f"Count matrix contains {nan_count} NaN/null values. Impute or filter before analysis."
        )

    # 6. Negative values check
    numeric_counts = counts_norm.select_dtypes(include=[np.number])
    if numeric_counts.shape != counts_norm.shape:
        errors.append("Count matrix contains non-numeric columns.")
    else:
        if (numeric_counts < 0).any().any():
            neg_count = (numeric_counts < 0).sum().sum()
            errors.append(
                f"Count matrix contains {neg_count} negative values. Raw counts must be >= 0."
            )

        # 7. Non-integer check (for raw counts)
        if not allow_float:
            is_non_integer = np.any(~np.isclose(numeric_counts.values % 1, 0, atol=1e-5))
            if is_non_integer:
                warnings.append(
                    "Count matrix contains floating point values with fractional parts. "
                    "PyDESeq2 expects raw unnormalized integer counts."
                )

    # 8. All-zero genes warning
    all_zero_genes = (
        (numeric_counts == 0).all(axis=1).sum() if numeric_counts.shape == counts_norm.shape else 0
    )
    if all_zero_genes > 0:
        warnings.append(
            f"{all_zero_genes} genes have zero counts across all samples ({all_zero_genes / len(counts_norm) * 100:.1f}%)."
        )

    is_valid = len(errors) == 0

    aligned_counts: Optional[pd.DataFrame] = None
    aligned_metadata: Optional[pd.DataFrame] = None

    if is_valid:
        common_samples = [s for s in counts_norm.columns if s in meta_df.index]
        if auto_align:
            aligned_counts = counts_norm[common_samples]
            aligned_metadata = meta_df.loc[common_samples]
            logger.info(
                f"Validation passed: {aligned_counts.shape[0]} genes across {len(common_samples)} samples aligned."
            )
        else:
            aligned_counts = counts_norm
            aligned_metadata = meta_df
    else:
        logger.error(f"Validation failed with {len(errors)} error(s).")

    report = ValidationReport(
        is_valid=is_valid,
        n_genes=counts.shape[0],
        n_samples=counts.shape[1],
        sample_ids=list(counts.columns),
        errors=errors,
        warnings=warnings,
        aligned_counts=aligned_counts,
        aligned_metadata=aligned_metadata,
    )

    if strict and not is_valid:
        raise ValidationError(report.summary())

    return report
