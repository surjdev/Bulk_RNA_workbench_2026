"""
Diagnostic evaluation and visual QC for batch effect correction (FR-7).
Provides side-by-side PCA comparisons and quantitative batch silhouette metrics.
"""

from __future__ import annotations

from typing import Any, Dict, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from btw import logger
from btw.viz.style import PALETTES, set_publication_style
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score


def evaluate_batch_effect(
    data: pd.DataFrame,
    metadata: pd.DataFrame,
    batch_col: str,
    condition_col: Optional[str] = None,
    n_components: int = 2,
) -> Dict[str, float]:
    """
    Compute quantitative metrics assessing batch effect and biological clustering.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix (genes x samples or samples x genes).
    metadata : pd.DataFrame
        Sample metadata.
    batch_col : str
        Metadata column containing batch assignments.
    condition_col : str, optional
        Metadata column containing biological condition.
    n_components : int, default=2
        Number of principal components.

    Returns
    -------
    dict
        Evaluation metrics: 'batch_silhouette', 'condition_silhouette', 'pc1_variance', 'pc2_variance'.
    """
    # Orient to (samples x genes)
    if set(metadata.index).issubset(data.columns):
        X = data.T.loc[metadata.index].values
    else:
        X = data.loc[metadata.index].values

    # Handle NaNs or zeros if any
    X = np.nan_to_num(X, nan=0.0)

    pca = PCA(n_components=n_components)
    coords = pca.fit_transform(X)
    var_exp = pca.explained_variance_ratio_

    metrics: Dict[str, float] = {
        "pc1_variance": float(var_exp[0]),
        "pc2_variance": float(var_exp[1]) if len(var_exp) > 1 else 0.0,
    }

    # Batch silhouette score (lower is better after correction, indicating less batch separation)
    batches = metadata[batch_col].astype(str).values
    if len(np.unique(batches)) > 1 and len(batches) > len(np.unique(batches)):
        try:
            metrics["batch_silhouette"] = float(silhouette_score(coords, batches))
        except Exception:
            metrics["batch_silhouette"] = 0.0
    else:
        metrics["batch_silhouette"] = 0.0

    # Biological silhouette score (higher is better, indicating strong biological grouping)
    if condition_col is not None and condition_col in metadata.columns:
        conds = metadata[condition_col].astype(str).values
        if len(np.unique(conds)) > 1 and len(conds) > len(np.unique(conds)):
            try:
                metrics["condition_silhouette"] = float(silhouette_score(coords, conds))
            except Exception:
                metrics["condition_silhouette"] = 0.0

    return metrics


def compare_pca_batch(
    data_before: pd.DataFrame,
    data_after: pd.DataFrame,
    metadata: pd.DataFrame,
    batch_col: str,
    condition_col: Optional[str] = None,
    figsize: Tuple[float, float] = (13.0, 5.5),
    title_prefix: str = "Batch Effect Comparison",
    palette: str = "nature",
    ax: Optional[Tuple[plt.Axes, plt.Axes]] = None,
) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes], Dict[str, Any]]:
    """
    Generate side-by-side publication-grade PCA plots before and after batch correction.

    Parameters
    ----------
    data_before : pd.DataFrame
        Uncorrected expression matrix.
    data_after : pd.DataFrame
        Batch-corrected expression matrix.
    metadata : pd.DataFrame
        Sample annotations.
    batch_col : str
        Column denoting technical batch.
    condition_col : str, optional
        Column denoting biological condition.
    figsize : tuple, default=(13, 5.5)
        Figure width and height.
    title_prefix : str
        Figure main title prefix.
    palette : str, default='nature'
        Color scheme name.
    ax : tuple of (ax1, ax2), optional
        Pre-existing subplot axes.

    Returns
    -------
    tuple of (plt.Figure, (ax_before, ax_after), metrics_dict)
    """
    set_publication_style()

    if ax is None:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    else:
        ax1, ax2 = ax
        fig = ax1.figure

    # Sample alignment
    if set(metadata.index).issubset(data_before.columns):
        X_before = data_before.T.loc[metadata.index].values
        X_after = data_after.T.loc[metadata.index].values
    else:
        X_before = data_before.loc[metadata.index].values
        X_after = data_after.loc[metadata.index].values

    # Run PCA
    pca_b = PCA(n_components=2).fit(X_before)
    coords_b = pca_b.transform(X_before)
    var_b = pca_b.explained_variance_ratio_

    pca_a = PCA(n_components=2).fit(X_after)
    coords_a = pca_a.transform(X_after)
    var_a = pca_a.explained_variance_ratio_

    # Metrics
    metrics_before = evaluate_batch_effect(data_before, metadata, batch_col, condition_col)
    metrics_after = evaluate_batch_effect(data_after, metadata, batch_col, condition_col)

    # Color mapping for batches
    unique_batches = metadata[batch_col].unique()
    colors = PALETTES.get(palette, PALETTES["nature"])
    batch_color_map = {b: colors[i % len(colors)] for i, b in enumerate(unique_batches)}

    # Marker mapping for conditions
    markers = ["o", "s", "^", "D", "v", "P"]
    unique_conds = metadata[condition_col].unique() if condition_col else [None]
    cond_marker_map = {c: markers[i % len(markers)] for i, c in enumerate(unique_conds)}

    # Plot helper
    for target_ax, coords, var_exp, title, metrics in [
        (ax1, coords_b, var_b, "Before Correction", metrics_before),
        (ax2, coords_a, var_a, "After ComBat Correction", metrics_after),
    ]:
        for i, sample in enumerate(metadata.index):
            b_val = metadata.loc[sample, batch_col]
            c_val = metadata.loc[sample, condition_col] if condition_col else None
            c_color = batch_color_map[b_val]
            m_style = cond_marker_map[c_val]

            target_ax.scatter(
                coords[i, 0],
                coords[i, 1],
                color=c_color,
                marker=m_style,
                s=110,
                alpha=0.85,
                edgecolors="#222222",
                linewidth=0.8,
            )

        target_ax.set_xlabel(f"PC1 ({var_exp[0] * 100:.1f}% variance)")
        target_ax.set_ylabel(f"PC2 ({var_exp[1] * 100:.1f}% variance)")
        target_ax.set_title(f"{title}\n(Batch Silhouette: {metrics['batch_silhouette']:.3f})")
        target_ax.grid(True, linestyle=":", alpha=0.5)

    # Legends
    # Batch legend (color)
    batch_handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=batch_color_map[b], markersize=9, label=f"Batch: {b}")
        for b in unique_batches
    ]
    cond_handles = []
    if condition_col:
        cond_handles = [
            plt.Line2D([0], [0], marker=cond_marker_map[c], color="w", markerfacecolor="#555555", markersize=9, label=f"{c}")
            for c in unique_conds
        ]

    ax2.legend(handles=batch_handles + cond_handles, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=True)
    fig.suptitle(title_prefix, fontsize=13, fontweight="bold", y=0.98)
    fig.tight_layout()

    combined_metrics = {
        "before": metrics_before,
        "after": metrics_after,
        "batch_silhouette_reduction": metrics_before["batch_silhouette"] - metrics_after["batch_silhouette"],
    }
    logger.info(
        f"Batch QC: Silhouette reduced from {metrics_before['batch_silhouette']:.3f} to {metrics_after['batch_silhouette']:.3f}."
    )

    return fig, (ax1, ax2), combined_metrics
