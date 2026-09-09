"""
Volcano plot visualization module for Differential Expression (FR-4).
Supports static publication-grade matplotlib plots with adjustText labels and interactive Plotly charts.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from btw import logger
from btw.de_analysis.deseq_helper import DEResult
from btw.viz.style import REGULATION_COLORS, set_publication_style

try:
    from adjustText import adjust_text
    HAS_ADJUST_TEXT = True
except ImportError:
    HAS_ADJUST_TEXT = False

try:
    import plotly.graph_objects as go
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False


def _extract_de_dataframe(de_data: Union[DEResult, pd.DataFrame]) -> pd.DataFrame:
    """Extract and validate standardized DataFrame from DEResult or raw DataFrame."""
    if isinstance(de_data, DEResult):
        return de_data.results_df.copy()
    elif isinstance(de_data, pd.DataFrame):
        return de_data.copy()
    else:
        raise TypeError(f"Expected DEResult or pd.DataFrame, got {type(de_data).__name__}")


def plot_volcano(
    de_data: Union[DEResult, pd.DataFrame],
    padj_col: str = "padj",
    lfc_col: str = "log2FoldChange",
    padj_cutoff: float = 0.05,
    lfc_cutoff: float = 1.0,
    top_n_labels: int = 15,
    highlight_genes: Optional[List[str]] = None,
    interactive: bool = False,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (7.0, 6.0),
    ax: Optional[plt.Axes] = None,
    colors: Optional[Dict[str, str]] = None,
) -> Any:
    """
    Generate Volcano plot highlighting significantly up- and down-regulated genes.

    Parameters
    ----------
    de_data : DEResult or pd.DataFrame
        Differential expression results table.
    padj_col : str, default='padj'
        Column name for adjusted p-value.
    lfc_col : str, default='log2FoldChange'
        Column name for log2 fold-change.
    padj_cutoff : float, default=0.05
        P-adjusted cutoff line.
    lfc_cutoff : float, default=1.0
        Log2FC cutoff lines (absolute value).
    top_n_labels : int, default=15
        Number of top significant genes to annotate with labels.
    highlight_genes : list of str, optional
        Specific gene symbols/IDs to force-label.
    interactive : bool, default=False
        If True, returns an interactive Plotly Figure; otherwise returns (fig, ax).
    title : str, optional
        Plot title.
    figsize : tuple of (width, height), default=(7, 6)
        Matplotlib figure size.
    ax : matplotlib.axes.Axes, optional
        Pre-existing axes to draw on.
    colors : dict, optional
        Custom color mapping for {'UP': ..., 'DOWN': ..., 'NS': ...}.

    Returns
    -------
    tuple of (plt.Figure, plt.Axes) or plotly.graph_objects.Figure
        Figure objects according to interactive parameter.
    """
    df = _extract_de_dataframe(de_data)

    if padj_col not in df.columns or lfc_col not in df.columns:
        raise KeyError(f"Columns '{padj_col}' and '{lfc_col}' must be present in DE data.")

    palette = colors if colors is not None else REGULATION_COLORS

    # Filter out NaNs
    valid = df[padj_col].notna() & df[lfc_col].notna() & (df[padj_col] > 0)
    plot_df = df.loc[valid].copy()

    # Compute -log10(padj)
    plot_df["neg_log10_padj"] = -np.log10(plot_df[padj_col])

    # Assign regulation status
    is_sig = (plot_df[padj_col] <= padj_cutoff) & (plot_df[lfc_col].abs() >= lfc_cutoff)
    conds = [
        is_sig & (plot_df[lfc_col] > 0),
        is_sig & (plot_df[lfc_col] < 0),
    ]
    plot_df["status"] = np.select(conds, ["UP", "DOWN"], default="NS")

    n_up = (plot_df["status"] == "UP").sum()
    n_down = (plot_df["status"] == "DOWN").sum()
    n_ns = (plot_df["status"] == "NS").sum()

    neg_log10_cutoff = -np.log10(padj_cutoff)
    default_title = title if title is not None else f"Volcano Plot (UP: {n_up}, DOWN: {n_down})"

    # -------------------------------------------------------------
    # Interactive Plotly Mode
    # -------------------------------------------------------------
    if interactive:
        if not HAS_PLOTLY:
            raise ImportError("plotly is required for interactive volcano plots. Run pip install plotly.")

        fig = go.Figure()

        for status, color in [("NS", palette["NS"]), ("DOWN", palette["DOWN"]), ("UP", palette["UP"])]:
            sub = plot_df[plot_df["status"] == status]
            count = len(sub)
            fig.add_trace(
                go.Scatter(
                    x=sub[lfc_col],
                    y=sub["neg_log10_padj"],
                    mode="markers",
                    name=f"{status} ({count})",
                    marker=dict(color=color, size=5, opacity=0.7),
                    text=sub.index,
                    hovertemplate=(
                        "<b>%{text}</b><br>"
                        "log2FC: %{x:.2f}<br>"
                        "-log10(padj): %{y:.2f}<br>"
                        "<extra></extra>"
                    ),
                )
            )

        # Add dashed threshold lines
        fig.add_hline(y=neg_log10_cutoff, line_dash="dash", line_color="#555555", line_width=1)
        fig.add_vline(x=lfc_cutoff, line_dash="dash", line_color="#555555", line_width=1)
        fig.add_vline(x=-lfc_cutoff, line_dash="dash", line_color="#555555", line_width=1)

        fig.update_layout(
            title=default_title,
            xaxis_title="log₂ Fold Change",
            yaxis_title="-log₁₀ (adjusted p-value)",
            template="plotly_white",
            hovermode="closest",
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

    # Plot Non-significant first as background
    ns_mask = plot_df["status"] == "NS"
    ax.scatter(
        plot_df.loc[ns_mask, lfc_col],
        plot_df.loc[ns_mask, "neg_log10_padj"],
        color=palette["NS"],
        s=12,
        alpha=0.4,
        label=f"NS ({n_ns})",
        edgecolors="none",
    )

    # Plot Down-regulated
    down_mask = plot_df["status"] == "DOWN"
    ax.scatter(
        plot_df.loc[down_mask, lfc_col],
        plot_df.loc[down_mask, "neg_log10_padj"],
        color=palette["DOWN"],
        s=20,
        alpha=0.8,
        label=f"DOWN ({n_down})",
        edgecolors="none",
    )

    # Plot Up-regulated
    up_mask = plot_df["status"] == "UP"
    ax.scatter(
        plot_df.loc[up_mask, lfc_col],
        plot_df.loc[up_mask, "neg_log10_padj"],
        color=palette["UP"],
        s=20,
        alpha=0.8,
        label=f"UP ({n_up})",
        edgecolors="none",
    )

    # Add threshold guidelines
    ax.axhline(neg_log10_cutoff, linestyle="--", color="#666666", linewidth=0.8, alpha=0.7)
    ax.axvline(lfc_cutoff, linestyle="--", color="#666666", linewidth=0.8, alpha=0.7)
    ax.axvline(-lfc_cutoff, linestyle="--", color="#666666", linewidth=0.8, alpha=0.7)

    # Candidate genes to label
    genes_to_label: List[str] = []
    if highlight_genes is not None:
        for g in highlight_genes:
            if g in plot_df.index:
                genes_to_label.append(g)

    if top_n_labels > 0:
        sig_candidates = plot_df.loc[is_sig].sort_values(
            by=["neg_log10_padj", lfc_col], ascending=[False, False]
        )
        for g in sig_candidates.index[:top_n_labels]:
            if g not in genes_to_label:
                genes_to_label.append(g)

    # Adjust text labels to prevent collisions
    texts = []
    for g in genes_to_label:
        row = plot_df.loc[g]
        texts.append(
            ax.text(
                row[lfc_col],
                row["neg_log10_padj"],
                str(g),
                fontsize=8,
                fontweight="medium",
            )
        )

    if HAS_ADJUST_TEXT and texts:
        adjust_text(
            texts,
            ax=ax,
            arrowprops=dict(arrowstyle="-", color="#444444", lw=0.6, alpha=0.8),
        )

    ax.set_xlabel("log₂ Fold Change")
    ax.set_ylabel("-log₁₀ adjusted p-value")
    ax.set_title(default_title)
    ax.legend(loc="upper right", frameon=False)

    return fig, ax
