"""
Unit tests for I/O and Validation module (FR-1).
"""

import numpy as np
import pandas as pd
import pytest

from btw.io import (
    BulkDataset,
    ValidationError,
    export_excel_multisheet,
    export_table,
    load_counts,
    load_dataset,
    load_metadata,
    validate_bulk_data,
)


def test_validation_valid_data(synthetic_data):
    """Test validation on pristine synthetic data passes without error."""
    counts, metadata = synthetic_data
    report = validate_bulk_data(counts, metadata)
    assert report.is_valid
    assert len(report.errors) == 0
    assert report.n_genes == 100
    assert report.n_samples == 6
    assert report.aligned_counts is not None
    assert report.aligned_metadata is not None


def test_validation_sample_mismatch(synthetic_data):
    """Test validation correctly flags mismatched samples."""
    counts, metadata = synthetic_data
    # Rename one sample in metadata
    bad_meta = metadata.copy()
    bad_meta.index = ["ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "WRONG_SAMPLE"]

    report = validate_bulk_data(counts, bad_meta)
    assert not report.is_valid
    assert any("missing in metadata" in err for err in report.errors)
    assert any("missing in counts" in err for err in report.errors)


def test_validation_nan_values(synthetic_data):
    """Test validation catches NaN values in count matrix."""
    counts, metadata = synthetic_data
    bad_counts = counts.copy().astype(float)
    bad_counts.iloc[0, 0] = np.nan

    report = validate_bulk_data(bad_counts, metadata)
    assert not report.is_valid
    assert any("NaN/null" in err for err in report.errors)


def test_validation_negative_values(synthetic_data):
    """Test validation catches negative values in count matrix."""
    counts, metadata = synthetic_data
    bad_counts = counts.copy()
    bad_counts.iloc[5, 2] = -10

    report = validate_bulk_data(bad_counts, metadata)
    assert not report.is_valid
    assert any("negative values" in err for err in report.errors)


def test_validation_duplicate_genes(synthetic_data):
    """Test validation catches duplicate gene IDs."""
    counts, metadata = synthetic_data
    bad_index = list(counts.index)
    bad_index[1] = bad_index[0]  # Duplicate first gene
    bad_counts = counts.copy()
    bad_counts.index = bad_index

    report = validate_bulk_data(bad_counts, metadata)
    assert not report.is_valid
    assert any("Duplicate gene identifiers" in err for err in report.errors)


def test_validation_strict_mode(synthetic_data):
    """Test strict=True raises ValidationError on failure."""
    counts, metadata = synthetic_data
    bad_counts = counts.copy()
    bad_counts.iloc[0, 0] = -5

    with pytest.raises(ValidationError):
        validate_bulk_data(bad_counts, metadata, strict=True)


def test_io_load_and_export_table(synthetic_data, temp_dir):
    """Test round-trip export and load for CSV, TSV, and Parquet."""
    counts, metadata = synthetic_data

    # 1. CSV test
    csv_path = temp_dir / "counts.csv"
    export_table(counts, csv_path)
    loaded_csv = load_counts(csv_path)
    assert loaded_csv.shape == counts.shape
    pd.testing.assert_frame_equal(loaded_csv, counts)

    # 2. TSV test
    tsv_path = temp_dir / "metadata.tsv"
    export_table(metadata, tsv_path)
    loaded_tsv = load_metadata(tsv_path)
    assert loaded_tsv.shape == metadata.shape
    pd.testing.assert_frame_equal(loaded_tsv, metadata)

    # 3. Parquet test (if pyarrow/fastparquet is installed)
    try:
        import pyarrow  # noqa: F401

        parquet_path = temp_dir / "counts.parquet"
        export_table(counts, parquet_path)
        loaded_parquet = load_counts(parquet_path)
        assert loaded_parquet.shape == counts.shape
    except ImportError:
        pass  # pyarrow not in test environment yet


def test_io_excel_multisheet_export(synthetic_data, temp_dir):
    """Test exporting multiple tables to a multi-sheet Excel workbook."""
    counts, metadata = synthetic_data
    excel_path = temp_dir / "summary_workbook.xlsx"

    sheets = {
        "Raw_Counts": counts,
        "Sample_Metadata": metadata,
    }
    out_path = export_excel_multisheet(sheets, excel_path)
    assert out_path.exists()

    # Read back sheets
    with pd.ExcelFile(excel_path) as xls:
        assert "Raw_Counts" in xls.sheet_names
        assert "Sample_Metadata" in xls.sheet_names
        df1 = pd.read_excel(xls, "Raw_Counts", index_col=0)
        assert df1.shape == counts.shape


def test_load_dataset_convenience(synthetic_data, temp_dir):
    """Test load_dataset loads, validates, and aligns dataset automatically."""
    counts, metadata = synthetic_data
    c_path = temp_dir / "raw_counts.csv"
    m_path = temp_dir / "samples.csv"

    export_table(counts, c_path)
    export_table(metadata, m_path)

    dataset = load_dataset(c_path, m_path, validate=True)
    assert isinstance(dataset, BulkDataset)
    assert dataset.n_genes == 100
    assert dataset.n_samples == 6
    assert len(dataset.genes) == 100
    assert len(dataset.samples) == 6
    assert dataset.validation_report.is_valid
