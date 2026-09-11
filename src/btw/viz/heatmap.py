"""
Hierarchical clustered heatmap visualization for Bulk Transcriptomics (FR-4).
Clusters top DEGs or variable genes across biological samples with metadata annotation bars.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from btw import logger
from btw.de_analysis.deseq_helper import DEResult
from btw.qc_normalize.normalize import NormalizationResult
from btw.viz.style import NATURE_PALETTE, set_publication_style


def _build_col_colors(
    metadata: pd.DataFrame,
    samples: List[str],
    annotation_cols: List[str],
) -> Tuple[pd.DataFrame, Dict[str, Dict[str, str]]]:
    """Build column color annotations DataFrame and legend mapping."""
    sub_meta = metadata.loc[samples, annotation_cols].copy()
    col_colors = pd.DataFrame(index=samples)
    legend_dict: Dict[str, Dict[str, str]] = {}

    for i, col in enumerate(annotation_cols):
        unique_vals = list(sub_meta[col].dropna().unique())
        colors = NATURE_PALETTE[i * 4 : i * 4 + len(unique_vals)]
        if len(colors) < len(unique_vals):
            # Fallback to general seaborn palette
            colors = sns.color_palette("tab10", len(unique_vals)).as_hex()

        mapping = {val: colors[j] for j, val in enumerate(unique_vals)}
        col_colors[col] = sub_meta[col].map(mapping)
        legend_dict[col] = mapping

    return col_colors, legend_dict


def plot_heatmap(
    data: Union[pd.DataFrame, NormalizationResult],
    metadata: Optional[pd.DataFrame] = None,
    de_result: Optional[DEResult] = None,
    genes: Optional[List[str]] = None,
    top_n_degs: Optional[int] = None,
    top_n_variable: int = 50,
    annotation_cols: Optional[List[str]] = None,
    z_score: bool = True,
    cluster_rows: bool = True,
    cluster_cols: bool = True,
    cmap: str = "vlag",
    figsize: Tuple[float, float] = (8.5, 7.5),
    show_gene_labels: bool = True,
    title: Optional[str] = None,
    **clustermap_kwargs,
) -> sns.matrix.ClusterGrid:
    """
    Generate hierarchical clustered heatmap of top DEGs or most variable genes.

    Parameters
    ----------
    data : pd.DataFrame or NormalizationResult
        Gene expression matrix (genes x samples, ideally log2 or VST counts).
    metadata : pd.DataFrame, optional
        Sample annotations.
    de_result : DEResult, optional
        If provided, automatically selects top DEGs from results_df.
    genes : list of str, optional
        Explicit gene IDs to plot.
    top_n_degs : int, optional
        Number of top DEGs to select when de_result is provided.
    top_n_variable : int, default=50
        Fallback number of top variable genes if no gene list is given.
    annotation_cols : list of str, optional
        Metadata columns to include as colored track headers above samples.
    z_score : bool, default=True
        Whether to standardize gene expression to z-scores across samples (row-wise).
    cluster_rows : bool, default=True
        Hierarchically cluster genes.
    cluster_cols : bool, default=True
        Hierarchically cluster samples.
    cmap : str, default='vlag'
        Diverging colormap (e.g. 'vlag', 'coolwarm', 'RdBu_r').
    figsize : tuple, default=(8.5, 7.5)
        Figure size.
    show_gene_labels : bool, default=True
        Whether to display gene names on row labels.
    title : str, optional
        Title on figure.

    Returns
    -------
    seaborn.matrix.ClusterGrid
        Clustermap object with dendrograms and heatmap.
    """
    if isinstance(data, NormalizationResult):
        mat = data.log2_counts
    else:
        mat = data.copy()

    # Determine gene subset
    if genes is not None and len(genes) > 0:
        chosen_genes = [g for g in genes if g in mat.index]
        if not chosen_genes:
            raise ValueError("None of the specified genes were found in the expression matrix.")
        sub_mat = mat.loc[chosen_genes]
    elif de_result is not None:
        degs = de_result.get_degs()
        if degs.empty:
            logger.warning(
                "No significant DEGs found in DEResult; falling back to top variable genes."
            )
            variances = mat.var(axis=1)
            n_sel = min(top_n_variable, mat.shape[0])
            chosen_genes = list(variances.nlargest(n_sel).index)
        else:
            n_top = top_n_degs if top_n_degs is not None else min(50, len(degs))
            sorted_degs = degs.sort_values(by="padj")
            chosen_genes = list(sorted_degs.index[:n_top])
        sub_mat = mat.loc[chosen_genes]
    else:
        variances = mat.var(axis=1)
        n_sel = min(top_n_variable, mat.shape[0])
        chosen_genes = list(variances.nlargest(n_sel).index)
        sub_mat = mat.loc[chosen_genes]

    # Row-wise Z-score standardization: (x - mean) / std
    if z_score:
        means = sub_mat.mean(axis=1)
        stds = sub_mat.std(axis=1)
        stds = np.where(stds == 0, 1.0, stds)
        plot_mat = sub_mat.sub(means, axis=0).div(stds, axis=0)
        cbar_label = "Z-score"
    else:
        plot_mat = sub_mat
        cbar_label = "Expression"

    # Build column color tracks if metadata provided
    col_colors = None
    legend_dict = None
    if metadata is not None and annotation_cols:
        valid_cols = [c for c in annotation_cols if c in metadata.columns]
        if valid_cols:
            col_colors, legend_dict = _build_col_colors(
                metadata=metadata,
                samples=list(plot_mat.columns),
                annotation_cols=valid_cols,
            )

    set_publication_style()

    # Build ClusterMap
    g = sns.clustermap(
        plot_mat,
        col_colors=col_colors,
        row_cluster=cluster_rows,
        col_cluster=cluster_cols,
        cmap=cmap,
        figsize=figsize,
        yticklabels=show_gene_labels,
        xticklabels=True,
        cbar_kws={"label": cbar_label, "orientation": "vertical"},
        **clustermap_kwargs,
    )

    if title is not None:
        g.fig.suptitle(title, y=1.02, fontsize=12, fontweight="bold")

    # Add legends for annotation tracks if present
    if legend_dict is not None and len(legend_dict) > 0:
        for label, mapping in legend_dict.items():
            handles = [
                plt.Line2D(
                    [0], [0], marker="s", color="w", markerfacecolor=col, markersize=8, label=val
                )
                for val, col in mapping.items()
            ]
            g.ax_heatmap.legend(
                handles=handles,
                title=label,
                bbox_to_anchor=(1.25, 1.0),
                loc="upper left",
                frameon=False,
            )

    return g
