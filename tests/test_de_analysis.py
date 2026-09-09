"""
Unit tests for Differential Expression Integration module (FR-3).
"""

import numpy as np
import pandas as pd
import pytest
from btw.de_analysis import (
    DEResult,
    MultiContrastResult,
    adjust_pvalues,
    apply_multitest_correction,
    build_deseq_dataset,
    run_de,
    run_deseq_stats,
    run_multiple_contrasts,
)


def test_adjust_pvalues_methods():
    """Test multiple testing correction methods with valid and NaN inputs."""
    pvals = [0.001, 0.01, 0.04, 0.06, 0.5, np.nan, 0.9]

    # Benjamini-Hochberg
    rejected_bh, padj_bh = adjust_pvalues(pvals, method="fdr_bh", alpha=0.05)
    assert len(padj_bh) == len(pvals)
    assert np.isnan(padj_bh[5])  # NaN preserved
    assert not rejected_bh[5]
    assert padj_bh[0] < 0.05
    assert rejected_bh[0]

    # Bonferroni (stricter)
    rejected_bonf, padj_bonf = adjust_pvalues(pvals, method="bonferroni", alpha=0.05)
    assert padj_bonf[0] >= padj_bh[0]
    assert np.isnan(padj_bonf[5])

    # Unsupported method raises ValueError
    with pytest.raises(ValueError):
        adjust_pvalues(pvals, method="invalid_method_name")


def test_apply_multitest_correction_df():
    """Test applying multiple testing correction onto a pandas DataFrame."""
    df = pd.DataFrame(
        {
            "gene": ["g1", "g2", "g3", "g4"],
            "pvalue": [0.0001, 0.02, 0.04, 0.8],
        }
    )

    res_df = apply_multitest_correction(df, pvalue_col="pvalue", method="fdr_bh", alpha=0.05)
    assert "padj" in res_df.columns
    assert "significant" in res_df.columns
    assert res_df.loc[0, "significant"]
    assert not res_df.loc[3, "significant"]


def test_build_deseq_dataset(synthetic_data):
    """Test constructing and fitting DeseqDataSet from synthetic data."""
    counts, metadata = synthetic_data
    dds = build_deseq_dataset(
        counts=counts,
        metadata=metadata,
        design_factors="condition",
        fit_model=True,
    )
    assert dds is not None
    # PyDESeq2 stores samples x genes
    assert dds.n_obs == counts.shape[1]
    assert dds.n_vars == counts.shape[0]


def test_run_deseq_stats_and_deresult(synthetic_data):
    """Test running DeseqStats and validating the DEResult container."""
    counts, metadata = synthetic_data
    dds = build_deseq_dataset(
        counts=counts,
        metadata=metadata,
        design_factors="condition",
        fit_model=True,
    )

    contrast = ("condition", "treated", "control")
    de_res = run_deseq_stats(dds, contrast=contrast, alpha=0.05, lfc_threshold=1.0)

    assert isinstance(de_res, DEResult)
    assert de_res.contrast == ("condition", "treated", "control")
    assert de_res.alpha == 0.05
    assert de_res.dds is dds
    assert de_res.stat_res is not None

    df = de_res.results_df
    assert df.shape[0] == 100
    for col in ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "significant", "regulation"]:
        assert col in df.columns

    # Test filtering methods
    degs = de_res.get_degs(padj_cutoff=0.05, lfc_cutoff=1.0)
    assert isinstance(degs, pd.DataFrame)

    up_genes = de_res.get_up_genes()
    down_genes = de_res.get_down_genes()
    assert isinstance(up_genes, list)
    assert isinstance(down_genes, list)

    summary_cnts = de_res.summary_counts()
    assert summary_cnts["total_tested"] == 100
    assert summary_cnts["significant_up"] == len(up_genes)
    assert summary_cnts["significant_down"] == len(down_genes)

    summary_text = de_res.summary()
    assert "DE Summary" in summary_text


def test_run_de_all_in_one(synthetic_data):
    """Test one-step run_de convenience function."""
    counts, metadata = synthetic_data
    contrast = ("condition", "treated", "control")

    res = run_de(counts, metadata, contrast=contrast, alpha=0.05, lfc_threshold=1.0)
    assert isinstance(res, DEResult)
    assert len(res.genes) == 100
    assert res.results_df.shape[0] == 100


def test_run_multiple_contrasts(synthetic_data, temp_dir):
    """Test batch multiple contrasts execution and MultiContrastResult features."""
    counts, metadata = synthetic_data
    # Two contrasts: condition (treated vs control) and batch (batch2 vs batch1)
    contrasts = [
        ("condition", "treated", "control"),
        ("batch", "batch2", "batch1"),
    ]

    multi_res = run_multiple_contrasts(
        counts=counts,
        metadata=metadata,
        contrasts=contrasts,
        alpha=0.05,
        lfc_threshold=1.0,
    )

    assert isinstance(multi_res, MultiContrastResult)
    assert len(multi_res.contrast_names) == 2
    assert "treated_vs_control" in multi_res.contrast_names
    assert "batch2_vs_batch1" in multi_res.contrast_names

    # Check retrieval of single contrast
    c1 = multi_res.get_contrast("treated_vs_control")
    assert isinstance(c1, DEResult)

    # Check summary table
    summary_tbl = multi_res.summary_table()
    assert isinstance(summary_tbl, pd.DataFrame)
    assert summary_tbl.shape[0] == 2

    # Check DEG sets
    deg_sets = multi_res.get_deg_sets()
    assert "treated_vs_control" in deg_sets
    assert isinstance(deg_sets["treated_vs_control"], set)

    # Check master table
    assert multi_res.master_table is not None
    assert "log2FoldChange_treated_vs_control" in multi_res.master_table.columns
    assert "log2FoldChange_batch2_vs_batch1" in multi_res.master_table.columns

    # Check Excel multi-sheet export
    excel_path = temp_dir / "multi_contrast_report.xlsx"
    out_path = multi_res.export_excel(excel_path)
    assert out_path.exists()
