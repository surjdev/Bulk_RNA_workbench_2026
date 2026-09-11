"""
Unit tests for QC and Normalization module (FR-2).
"""

import numpy as np
import pandas as pd

from btw.qc_normalize import (
    NormalizationResult,
    compute_gene_qc,
    compute_sample_qc,
    compute_size_factors,
    filter_low_expression_genes,
    normalize_cpm,
    normalize_deseq2,
)


def test_sample_qc_metrics(synthetic_data):
    """Test sample QC statistics calculation."""
    counts, _ = synthetic_data
    qc_df = compute_sample_qc(counts, min_count=1)

    assert isinstance(qc_df, pd.DataFrame)
    assert qc_df.shape == (6, 8)
    assert list(qc_df.index) == list(counts.columns)

    # Check metrics existence
    for col in ["library_size", "detected_genes", "detection_rate_pct", "mean_count"]:
        assert col in qc_df.columns

    # Verify library size is non-negative and matches sum
    assert (qc_df["library_size"] == counts.sum(axis=0)).all()
    assert (qc_df["detection_rate_pct"] <= 100.0).all()
    assert (qc_df["detection_rate_pct"] >= 0.0).all()


def test_gene_qc_metrics(synthetic_data):
    """Test gene QC statistics calculation."""
    counts, _ = synthetic_data
    gene_df = compute_gene_qc(counts, min_count=5)

    assert isinstance(gene_df, pd.DataFrame)
    assert gene_df.shape == (100, 6)
    assert list(gene_df.index) == list(counts.index)

    for col in [
        "total_counts",
        "mean_count",
        "variance",
        "dispersion_ratio",
        "n_samples_expressing",
    ]:
        assert col in gene_df.columns

    assert (gene_df["n_samples_expressing"] <= 6).all()
    assert (gene_df["n_samples_expressing"] >= 0).all()


def test_filter_low_expression(synthetic_data):
    """Test filtering of low-count genes."""
    counts, _ = synthetic_data
    filtered, summary = filter_low_expression_genes(counts, min_counts=10, min_samples=3)

    assert isinstance(filtered, pd.DataFrame)
    assert isinstance(summary, pd.DataFrame)
    assert filtered.shape[0] <= counts.shape[0]
    assert filtered.shape[1] == counts.shape[1]

    # Every retained gene must have at least 10 counts in at least 3 samples
    valid_retained = ((filtered >= 10).sum(axis=1) >= 3).all()
    assert valid_retained


def test_size_factors_and_deseq2_normalization(synthetic_data):
    """Test size factors calculation and DESeq2 median-of-ratios normalization."""
    counts, metadata = synthetic_data
    sfs = compute_size_factors(counts)

    assert isinstance(sfs, pd.Series)
    assert len(sfs) == 6
    assert (sfs > 0).all()

    # Run thin normalization helper
    result = normalize_deseq2(counts, metadata=metadata, design_factors="condition")
    assert isinstance(result, NormalizationResult)
    assert result.normalized_counts.shape == counts.shape
    assert result.method == "deseq2_median_of_ratios"

    # Test log2 property
    log2_df = result.log2_counts
    assert log2_df.shape == counts.shape
    assert (log2_df >= 0).all().all()

    # Test custom override of size factors
    custom_sfs = pd.Series([1.0, 1.0, 1.0, 1.0, 1.0, 1.0], index=counts.columns)
    override_result = normalize_deseq2(counts, override_size_factors=custom_sfs)
    pd.testing.assert_frame_equal(override_result.normalized_counts, counts.astype(float))


def test_cpm_normalization(synthetic_data):
    """Test CPM and log2(CPM+1) normalization."""
    counts, _ = synthetic_data
    cpm = normalize_cpm(counts, log2_transform=False)

    assert isinstance(cpm, pd.DataFrame)
    assert cpm.shape == counts.shape

    # Each sample column should sum to approximately 1,000,000
    column_sums = cpm.sum(axis=0)
    for s in column_sums:
        assert np.isclose(s, 1e6, rtol=1e-3)

    # Log2 transform check
    log_cpm = normalize_cpm(counts, log2_transform=True, prior_count=1.0)
    assert (log_cpm >= 0).all().all()
