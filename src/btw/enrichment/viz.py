"""
Centralized enrichment visualization module (FR-5).
Provides unified Dot Plots and Bar Plots applicable across ORA, GSEA, and decoupler activity results.
"""

from __future__ import annotations

from typing import Any, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from btw import logger
from btw.enrichment.schema import EnrichmentResult
from btw.viz.style import set_publication_style

try:
    import plotly.express as px
    import plotly.graph_objects as go
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False


def _resolve_enrichment_df(
    enrichment: Union[EnrichmentResult, pd.DataFrame],
    padj_cutoff: Optional[float] = None,
    top_n: int = 15,
) -> pd.DataFrame:
    """Extract, filter by padj_cutoff, and sort top N terms from enrichment input."""
    if isinstance(enrichment, EnrichmentResult):
        df = enrichment.results_df.copy()
        alpha = padj_cutoff if padj_cutoff is not None else enrichment.alpha
    else:
        df = enrichment.copy()
        alpha = padj_cutoff if padj_cutoff is not None else 0.05

    if df.empty:
        return df

    # Filter by significance if padj available
    if "padj" in df.columns:
        valid = df.loc[df["padj"].notna() & (df["padj"] <= alpha)]
        if not valid.empty:
            df = valid
        else:
            logger.info("No terms below significance threshold; plotting top terms by nominal p-value.")

    # Sort by padj ascending (most significant first)
    if "padj" in df.columns and "pvalue" in df.columns:
        df = df.sort_values(by=["padj", "pvalue"], ascending=[True, True])
    elif "pvalue" in df.columns:
        df = df.sort_values(by="pvalue", ascending=True)
    elif "score" in df.columns:
        df = df.sort_values(by="score", key=abs, ascending=False)

    return df.head(top_n).copy()


def plot_enrichment_dotplot(
    enrichment: Union[EnrichmentResult, pd.DataFrame],
    top_n: int = 15,
    size_by: str = "gene_count",
    color_by: str = "padj",
    padj_cutoff: Optional[float] = None,
    interactive: bool = False,
    figsize: Tuple[float, float] = (8.0, 6.0),
    title: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
) -> Any:
    """
    Generate unified Dot Plot showing enriched pathways/terms ranked by significance.
    Dot size reflects gene count; dot color reflects -log10(padj).

    Parameters
    ----------
    enrichment : EnrichmentResult or pd.DataFrame
        Standardized enrichment table.
    top_n : int, default=15
        Number of top terms to display.
    size_by : str, default='gene_count'
        Column mapped to dot diameter.
    color_by : str, default='padj'
        Column mapped to color gradient.
    padj_cutoff : float, optional
        Significance cutoff.
    interactive : bool, default=False
        Whether to generate interactive Plotly chart.
    figsize : tuple, default=(8, 6)
        Matplotlib figure size.
    title : str, optional
        Plot title.
    ax : plt.Axes, optional
        Target matplotlib axes.

    Returns
    -------
    tuple of (plt.Figure, plt.Axes) or plotly.graph_objects.Figure
    """
    plot_df = _resolve_enrichment_df(enrichment, padj_cutoff=padj_cutoff, top_n=top_n)

    if plot_df.empty:
        raise ValueError("Enrichment table has 0 rows to plot.")

    # Compute -log10(padj) for color scale
    pvals = np.clip(plot_df[color_by].values.astype(float), 1e-100, 1.0)
    plot_df["neg_log10_padj"] = -np.log10(pvals)

    # Invert order so highest ranked term is at top of y-axis
    plot_df = plot_df.iloc[::-1].reset_index(drop=True)

    default_title = title if title is not None else "Pathway Enrichment Dot Plot"

    # -------------------------------------------------------------
    # Interactive Plotly Mode
    # -------------------------------------------------------------
    if interactive:
        if not HAS_PLOTLY:
            raise ImportError("plotly is required for interactive dot plots.")

        size_col = size_by if size_by in plot_df.columns else None
        fig = px.scatter(
            plot_df,
            x="score",
            y="term",
            size=size_col,
            color="neg_log10_padj",
            color_continuous_scale="Viridis",
            labels={"score": "Enrichment Score / Activity", "neg_log10_padj": "-log₁₀(padj)", "term": "Pathway / Term"},
            title=default_title,
            hover_data=["pvalue", "padj", "gene_count"],
            template="plotly_white",
        )
        return fig

    # -------------------------------------------------------------
    # Static Matplotlib Mode
    # -------------------------------------------------------------
    set_publication_style()

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # Determine point sizes
    if size_by in plot_df.columns:
        counts = pd.to_numeric(plot_df[size_by], errors="coerce").fillna(0).values.astype(float)
        counts = np.maximum(counts, 0.0)
        min_c, max_c = float(counts.min()), float(counts.max())
        if max_c <= min_c:
            sizes = np.full_like(counts, 120.0, dtype=float)
        else:
            sizes = 40.0 + (counts - min_c) / (max_c - min_c) * 220.0
            sizes = np.clip(sizes, 30.0, 300.0)
    else:
        sizes = np.full(len(plot_df), 120.0, dtype=float)

    scatter = ax.scatter(
        plot_df["score"],
        np.arange(len(plot_df)),
        s=sizes,
        c=plot_df["neg_log10_padj"],
        cmap="viridis",
        edgecolors="#222222",
        linewidths=0.7,
        alpha=0.9,
    )

    ax.set_yticks(np.arange(len(plot_df)))
    ax.set_yticklabels(plot_df["term"])
    ax.set_xlabel("Enrichment Score / Activity")
    ax.set_title(default_title)
    ax.grid(axis="x", linestyle=":", alpha=0.6)

    # Colorbar
    cbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label("-log₁₀ adjusted p-value", rotation=270, labelpad=15)

    return fig, ax


def plot_enrichment_barplot(
    enrichment: Union[EnrichmentResult, pd.DataFrame],
    top_n: int = 15,
    metric: str = "score",
    color_by: str = "padj",
    padj_cutoff: Optional[float] = None,
    interactive: bool = False,
    figsize: Tuple[float, float] = (8.0, 5.0),
    title: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
) -> Any:
    """
    Generate unified Bar Plot for enrichment results ranked by effect size or significance.

    Parameters
    ----------
    enrichment : EnrichmentResult or pd.DataFrame
        Standardized enrichment table.
    top_n : int, default=15
        Number of top terms.
    metric : str, default='score'
        Bar length value ('score', 'gene_count', '-log10_padj').
    color_by : str, default='padj'
        Color gradient mapping.
    padj_cutoff : float, optional
        Significance threshold.
    interactive : bool, default=False
        Whether to generate interactive Plotly chart.
    figsize : tuple, default=(8, 5)
        Matplotlib figure dimensions.
    title : str, optional
        Plot title.
    ax : plt.Axes, optional
        Pre-existing axes.

    Returns
    -------
    tuple of (plt.Figure, plt.Axes) or plotly.graph_objects.Figure
    """
    plot_df = _resolve_enrichment_df(enrichment, padj_cutoff=padj_cutoff, top_n=top_n)

    if plot_df.empty:
        raise ValueError("Enrichment table has 0 rows to plot.")

    pvals = np.clip(plot_df[color_by].values.astype(float), 1e-100, 1.0)
    plot_df["neg_log10_padj"] = -np.log10(pvals)

    # Invert order for top-to-bottom bar alignment
    plot_df = plot_df.iloc[::-1].reset_index(drop=True)

    default_title = title if title is not None else "Enriched Pathways Significance"

    # -------------------------------------------------------------
    # Interactive Plotly Mode
    # -------------------------------------------------------------
    if interactive:
        if not HAS_PLOTLY:
            raise ImportError("plotly is required.")

        fig = px.bar(
            plot_df,
            x=metric,
            y="term",
            color="neg_log10_padj",
            color_continuous_scale="Tealrose",
            orientation="h",
            title=default_title,
            labels={metric: metric.capitalize(), "neg_log10_padj": "-log₁₀(padj)", "term": "Pathway"},
            template="plotly_white",
        )
        return fig

    # -------------------------------------------------------------
    # Static Matplotlib Mode
    # -------------------------------------------------------------
    set_publication_style()

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    bars = ax.barh(
        np.arange(len(plot_df)),
        plot_df[metric],
        color="#3C5488",
        edgecolor="#222222",
        linewidth=0.6,
        alpha=0.85,
    )

    ax.set_yticks(np.arange(len(plot_df)))
    ax.set_yticklabels(plot_df["term"])
    ax.set_xlabel(metric.capitalize())
    ax.set_title(default_title)
    ax.grid(axis="x", linestyle=":", alpha=0.5)

    return fig, ax
