"""
MA plot visualization module for Differential Expression (FR-4).
Plots log2 fold change versus mean expression across all samples.
"""

from __future__ import annotations

from typing import Any, Dict, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from btw.de_analysis.deseq_helper import DEResult
from btw.viz.style import REGULATION_COLORS, set_publication_style

try:
    import plotly.graph_objects as go

    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False


def plot_ma(
    de_data: Union[DEResult, pd.DataFrame],
    base_mean_col: str = "baseMean",
    lfc_col: str = "log2FoldChange",
    padj_col: str = "padj",
    padj_cutoff: float = 0.05,
    lfc_cutoff: float = 1.0,
    interactive: bool = False,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (7.0, 5.5),
    ax: Optional[plt.Axes] = None,
    colors: Optional[Dict[str, str]] = None,
) -> Any:
    """
    Generate an MA plot showing log2 fold changes over mean expression levels.

    Parameters
    ----------
    de_data : DEResult or pd.DataFrame
        Differential expression results table.
    base_mean_col : str, default='baseMean'
        Mean expression column.
    lfc_col : str, default='log2FoldChange'
        Log2 fold change column.
    padj_col : str, default='padj'
        Adjusted p-value column.
    padj_cutoff : float, default=0.05
        Significance cutoff.
    lfc_cutoff : float, default=1.0
        Log2FC threshold line.
    interactive : bool, default=False
        Whether to generate interactive Plotly chart.
    title : str, optional
        Plot title.
    figsize : tuple, default=(7, 5.5)
        Matplotlib figure dimensions.
    ax : plt.Axes, optional
        Target matplotlib axes.
    colors : dict, optional
        Color map.

    Returns
    -------
    tuple of (plt.Figure, plt.Axes) or plotly.graph_objects.Figure
    """
    if isinstance(de_data, DEResult):
        df = de_data.results_df.copy()
    else:
        df = de_data.copy()

    for c in [base_mean_col, lfc_col, padj_col]:
        if c not in df.columns:
            raise KeyError(f"Required column '{c}' missing from DE data.")

    palette = colors if colors is not None else REGULATION_COLORS

    valid = df[base_mean_col].notna() & df[lfc_col].notna() & (df[base_mean_col] > 0)
    plot_df = df.loc[valid].copy()

    is_sig = (
        (plot_df[padj_col].notna())
        & (plot_df[padj_col] <= padj_cutoff)
        & (plot_df[lfc_col].abs() >= lfc_cutoff)
    )
    conds = [
        is_sig & (plot_df[lfc_col] > 0),
        is_sig & (plot_df[lfc_col] < 0),
    ]
    plot_df["status"] = np.select(conds, ["UP", "DOWN"], default="NS")

    n_up = (plot_df["status"] == "UP").sum()
    n_down = (plot_df["status"] == "DOWN").sum()
    default_title = title if title is not None else f"MA Plot (UP: {n_up}, DOWN: {n_down})"

    # -------------------------------------------------------------
    # Interactive Plotly Mode
    # -------------------------------------------------------------
    if interactive:
        if not HAS_PLOTLY:
            raise ImportError("plotly is required for interactive MA plots.")

        fig = go.Figure()
        for status, color in [
            ("NS", palette["NS"]),
            ("DOWN", palette["DOWN"]),
            ("UP", palette["UP"]),
        ]:
            sub = plot_df[plot_df["status"] == status]
            fig.add_trace(
                go.Scatter(
                    x=sub[base_mean_col],
                    y=sub[lfc_col],
                    mode="markers",
                    name=f"{status} ({len(sub)})",
                    marker=dict(color=color, size=5, opacity=0.7),
                    text=sub.index,
                    hovertemplate=(
                        "<b>%{text}</b><br>"
                        "Base Mean: %{x:.1f}<br>"
                        "log2FC: %{y:.2f}<br>"
                        "<extra></extra>"
                    ),
                )
            )

        fig.add_hline(y=0, line_color="#333333", line_width=1)
        fig.add_hline(y=lfc_cutoff, line_dash="dash", line_color="#666666", line_width=1)
        fig.add_hline(y=-lfc_cutoff, line_dash="dash", line_color="#666666", line_width=1)

        fig.update_layout(
            title=default_title,
            xaxis_title="Mean of Normalized Counts",
            xaxis_type="log",
            yaxis_title="log₂ Fold Change",
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

    # Plot non-significant points
    ns_sub = plot_df[plot_df["status"] == "NS"]
    ax.scatter(
        ns_sub[base_mean_col],
        ns_sub[lfc_col],
        color=palette["NS"],
        s=10,
        alpha=0.35,
        label=f"NS ({len(ns_sub)})",
    )

    # Plot down-regulated points
    down_sub = plot_df[plot_df["status"] == "DOWN"]
    ax.scatter(
        down_sub[base_mean_col],
        down_sub[lfc_col],
        color=palette["DOWN"],
        s=18,
        alpha=0.8,
        label=f"DOWN ({len(down_sub)})",
    )

    # Plot up-regulated points
    up_sub = plot_df[plot_df["status"] == "UP"]
    ax.scatter(
        up_sub[base_mean_col],
        up_sub[lfc_col],
        color=palette["UP"],
        s=18,
        alpha=0.8,
        label=f"UP ({len(up_sub)})",
    )

    ax.axhline(0, color="#333333", linewidth=0.8)
    ax.axhline(lfc_cutoff, linestyle="--", color="#666666", linewidth=0.8, alpha=0.7)
    ax.axhline(-lfc_cutoff, linestyle="--", color="#666666", linewidth=0.8, alpha=0.7)

    ax.set_xscale("log")
    ax.set_xlabel("Mean of Normalized Counts (log scale)")
    ax.set_ylabel("log₂ Fold Change")
    ax.set_title(default_title)
    ax.legend(loc="upper right", frameon=False)

    return fig, ax
