"""
Dispersion plot visualization module for PyDESeq2 outputs (FR-4).
Visualizes genewise dispersions, fitted trend, and MAP shrinkage estimates.
"""

from __future__ import annotations

from typing import Any, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from btw import logger
from btw.de_analysis.deseq_helper import DEResult
from btw.viz.style import set_publication_style


def plot_dispersion(
    source: Union[Any, DEResult],
    figsize: Tuple[float, float] = (7.0, 5.5),
    ax: Optional[plt.Axes] = None,
    title: Optional[str] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot dispersion estimates versus mean of normalized counts from PyDESeq2.

    Parameters
    ----------
    source : DeseqDataSet or DEResult
        PyDESeq2 object containing dispersion fits.
    figsize : tuple, default=(7, 5.5)
        Figure size.
    ax : plt.Axes, optional
        Pre-existing axes.
    title : str, optional
        Plot title.

    Returns
    -------
    tuple of (plt.Figure, plt.Axes)
    """
    dds = source.dds if isinstance(source, DEResult) else source

    if dds is None or not hasattr(dds, "varm"):
        raise ValueError("Provided object does not contain fitted PyDESeq2 dispersion tables.")

    # Check for required dispersion attributes in var (PyDESeq2 0.5+) or varm
    var_df = getattr(dds, "var", pd.DataFrame())
    varm = getattr(dds, "varm", {})

    genewise_disp = None
    final_disp = None
    fitted_disp = None
    mean_counts = None

    if "genewise_dispersions" in var_df.columns:
        genewise_disp = var_df["genewise_dispersions"].values
        final_disp = var_df["dispersions"].values
        fitted_disp = var_df["fitted_dispersions"].values if "fitted_dispersions" in var_df.columns else None
        mean_counts = var_df["_normed_means"].values if "_normed_means" in var_df.columns else None
    elif "genewise_dispersions" in varm:
        genewise_disp = np.asarray(varm["genewise_dispersions"])
        final_disp = np.asarray(varm["dispersions"])
        fitted_disp = np.asarray(varm["fitted_dispersions"]) if "fitted_dispersions" in varm else None

    if genewise_disp is None or final_disp is None:
        raise KeyError("DeseqDataSet has not had dispersions fitted. Run dds.deseq2() first.")

    # Mean normalized counts per gene
    if mean_counts is None:
        if hasattr(dds, "layers") and "normed_counts" in dds.layers:
            mean_counts = np.mean(dds.layers["normed_counts"], axis=0)
        else:
            mean_counts = np.mean(dds.X, axis=0)

    set_publication_style()

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # Filter non-positive values for log scale
    valid = (mean_counts > 0) & (genewise_disp > 0) & (final_disp > 0)
    means_clean = mean_counts[valid]
    gw_clean = genewise_disp[valid]
    fn_clean = final_disp[valid]

    # Genewise estimates (black circles)
    ax.scatter(
        means_clean,
        gw_clean,
        color="#222222",
        s=12,
        alpha=0.45,
        label="Gene-wise",
        edgecolors="none",
    )

    # Fitted trend (red line / dots)
    if fitted_disp is not None:
        fit_clean = fitted_disp[valid]
        # Sort by mean for smooth curve
        sort_idx = np.argsort(means_clean)
        ax.plot(
            means_clean[sort_idx],
            fit_clean[sort_idx],
            color="#E64B35",
            linewidth=1.8,
            label="Fitted trend",
        )

    # Final MAP estimates (blue circles)
    ax.scatter(
        means_clean,
        fn_clean,
        facecolors="none",
        edgecolors="#3C5488",
        s=18,
        alpha=0.6,
        linewidths=0.7,
        label="Final (MAP)",
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Mean of Normalized Counts")
    ax.set_ylabel("Dispersion")
    ax.set_title(title if title is not None else "PyDESeq2 Dispersion Estimates")
    ax.legend(loc="upper right", frameon=False)

    return fig, ax
