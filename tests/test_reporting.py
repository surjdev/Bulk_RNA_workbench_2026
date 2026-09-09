"""
Unit tests for Reporting Module (FR-9).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pytest
from btw.reporting import (
    export_publication_bundle,
    generate_html_report,
    generate_markdown_report,
)
from btw.viz.style import save_figure


def test_markdown_and_html_reports(de_result_fixture, tmp_path):
    """Test generating Markdown and standalone HTML reports."""
    # Create a dummy plot to embed
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.plot([1, 2], [3, 4])
    plot_path = tmp_path / "dummy_plot.png"
    save_figure(fig, plot_path)
    plt.close(fig)

    fig_dict = {"Sample Trend": plot_path}
    qc_data = {"Total Samples": 6, "Sequencing Depth": "30M reads"}

    # 1. Test Markdown Report
    md_out = tmp_path / "report.md"
    md_text = generate_markdown_report(
        title="Test Transcriptomics Report",
        de_result=de_result_fixture,
        qc_metrics=qc_data,
        figure_paths=fig_dict,
        output_path=md_out,
    )
    assert md_out.exists()
    assert "Test Transcriptomics Report" in md_text
    assert "Executive Summary" in md_text
    assert "Significant UP" in md_text

    # 2. Test HTML Report
    html_out = tmp_path / "report.html"
    html_text = generate_html_report(
        title="Test Transcriptomics HTML",
        de_result=de_result_fixture,
        qc_metrics=qc_data,
        figure_paths=fig_dict,
        embed_images=True,
        output_path=html_out,
    )
    assert html_out.exists()
    assert "<!DOCTYPE html>" in html_text
    assert "data:image/png;base64," in html_text
    assert "Differential Expression Summary" in html_text


def test_export_publication_bundle(de_result_fixture, tmp_path):
    """Test complete publication bundle export with Excel, figures, and HTML."""
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.scatter([1, 2], [3, 4])
    plot_path = tmp_path / "pca_figure.png"
    save_figure(fig, plot_path)
    plt.close(fig)

    bundle_dir = tmp_path / "bundle_output"
    exported = export_publication_bundle(
        de_result=de_result_fixture,
        output_dir=bundle_dir,
        figure_paths={"PCA Clustering": plot_path},
        bundle_name="Study_Report",
    )

    assert exported.exists()
    assert (bundle_dir / "tables" / "differential_expression_master.csv").exists()
    assert (bundle_dir / "tables" / "complete_analysis_workbook.xlsx").exists()
    assert (bundle_dir / "figures" / "pca_figure.png").exists()
    assert (bundle_dir / "Study_Report.html").exists()
    assert (bundle_dir / "Study_Report.md").exists()
