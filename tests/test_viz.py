"""
Unit tests for Publication-Grade Visualization module (FR-4).
"""

import matplotlib

matplotlib.use("Agg")  # Non-interactive backend for testing
import matplotlib.pyplot as plt
import pytest

from btw.de_analysis import run_de, run_multiple_contrasts
from btw.viz import (
    compute_pca,
    plot_dispersion,
    plot_heatmap,
    plot_ma,
    plot_pca,
    plot_upset,
    plot_venn,
    plot_volcano,
    save_figure,
    set_publication_style,
)


@pytest.fixture
def de_result_fixture(synthetic_data):
    """Run DE analysis on synthetic data for visualization tests."""
    counts, metadata = synthetic_data
    return run_de(counts, metadata, contrast=("condition", "treated", "control"))


def test_style_and_save_figure(temp_dir):
    """Test publication style application and saving figures."""
    set_publication_style(palette="nature", dpi=150)
    fig, ax = plt.subplots()
    ax.plot([1, 2, 3], [4, 5, 6])

    png_path = temp_dir / "test_fig.png"
    out = save_figure(fig, png_path, dpi=150)
    assert out.exists()
    plt.close(fig)


def test_volcano_plot_static_and_interactive(de_result_fixture):
    """Test volcano plot generation in static and interactive modes."""
    # 1. Static mode
    fig, ax = plot_volcano(
        de_result_fixture,
        padj_cutoff=0.05,
        lfc_cutoff=1.0,
        top_n_labels=5,
        interactive=False,
    )
    assert isinstance(fig, plt.Figure)
    assert isinstance(ax, plt.Axes)
    assert len(ax.collections) > 0  # Scatter points exist
    plt.close(fig)

    # 2. Interactive mode (Plotly)
    plotly_fig = plot_volcano(
        de_result_fixture,
        padj_cutoff=0.05,
        lfc_cutoff=1.0,
        interactive=True,
    )
    assert plotly_fig is not None
    assert len(plotly_fig.data) == 3  # UP, DOWN, NS traces


def test_pca_computation_and_plots(synthetic_data):
    """Test PCA calculation and static/interactive plotting."""
    counts, metadata = synthetic_data

    # Test computation
    pca, pca_df = compute_pca(counts, metadata, n_components=2, top_n_variable_genes=50)
    assert pca_df.shape == (6, 4)  # PC1, PC2 + condition, batch
    assert len(pca.explained_variance_ratio_) == 2

    # Test static plot
    fig, ax, df_res = plot_pca(
        counts,
        metadata,
        color_by="condition",
        shape_by="batch",
        top_n_variable_genes=50,
        interactive=False,
    )
    assert isinstance(fig, plt.Figure)
    assert df_res.shape[0] == 6
    plt.close(fig)

    # Test interactive plot
    plotly_fig, _ = plot_pca(
        counts,
        metadata,
        color_by="condition",
        interactive=True,
    )
    assert plotly_fig is not None


def test_ma_plot_static_and_interactive(de_result_fixture):
    """Test MA plot in both static and interactive modes."""
    fig, ax = plot_ma(de_result_fixture, interactive=False)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)

    plotly_fig = plot_ma(de_result_fixture, interactive=True)
    assert plotly_fig is not None


def test_dispersion_plot(de_result_fixture):
    """Test dispersion plot using fitted PyDESeq2 DeseqDataSet."""
    fig, ax = plot_dispersion(de_result_fixture)
    assert isinstance(fig, plt.Figure)
    assert ax.get_xscale() == "log"
    assert ax.get_yscale() == "log"
    plt.close(fig)


def test_heatmap_plot(synthetic_data, de_result_fixture):
    """Test hierarchical clustered heatmap with metadata tracks."""
    counts, metadata = synthetic_data

    # Top DEGs heatmap
    g = plot_heatmap(
        counts,
        metadata=metadata,
        de_result=de_result_fixture,
        top_n_degs=20,
        annotation_cols=["condition", "batch"],
        z_score=True,
    )
    assert g is not None
    assert hasattr(g, "ax_heatmap")
    plt.close(g.fig)


def test_overlaps_upset_and_venn(synthetic_data):
    """Test cross-contrast DEG overlap plots (UpSet and Venn)."""
    counts, metadata = synthetic_data
    contrasts = [
        ("condition", "treated", "control"),
        ("batch", "batch2", "batch1"),
    ]
    multi_res = run_multiple_contrasts(counts, metadata, contrasts=contrasts)

    # Test UpSet plot
    axes_dict = plot_upset(multi_res, padj_cutoff=0.05, lfc_cutoff=1.0)
    assert axes_dict is not None
    plt.close("all")

    # Test Venn diagram (or fallback)
    fig, ax = plot_venn(multi_res)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)
