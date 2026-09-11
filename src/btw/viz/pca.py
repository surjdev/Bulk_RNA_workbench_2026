"""
Principal Component Analysis (PCA) visualization module for Bulk Transcriptomics (FR-4).
Uses scikit-learn PCA to compute variance-explained embeddings with static and interactive views.
"""

from __future__ import annotations

from typing import Any, Optional, Tuple, Union

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA

from btw import logger
from btw.qc_normalize.normalize import NormalizationResult
from btw.viz.style import NATURE_PALETTE, set_publication_style

try:
    import plotly.express as px
    import plotly.graph_objects as go  # noqa: F401

    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False


def compute_pca(
    data: Union[pd.DataFrame, NormalizationResult],
    metadata: pd.DataFrame,
    n_components: int = 2,
    top_n_variable_genes: int = 500,
) -> Tuple[PCA, pd.DataFrame]:
    """
    Perform PCA on the most variable genes across samples.

    Parameters
    ----------
    data : pd.DataFrame or NormalizationResult
        Gene expression matrix of shape (n_genes, n_samples).
    metadata : pd.DataFrame
        Sample annotations table.
    n_components : int, default=2
        Number of principal components to calculate.
    top_n_variable_genes : int, default=500
        Number of top variable genes by variance to include.

    Returns
    -------
    pca : PCA
        Fitted scikit-learn PCA estimator.
    pca_df : pd.DataFrame
        Sample coordinate table with PC scores and metadata columns.
    """
    if isinstance(data, NormalizationResult):
        mat = data.log2_counts
    else:
        mat = data.copy()

    # Align samples with metadata
    common_samples = [s for s in mat.columns if s in metadata.index]
    if len(common_samples) < 3:
        raise ValueError(f"Need at least 3 matching samples for PCA, found {len(common_samples)}")

    mat = mat[common_samples]
    meta_aligned = metadata.loc[common_samples]

    # Select top variable genes
    variances = mat.var(axis=1)
    n_top = min(top_n_variable_genes, mat.shape[0])
    top_genes = variances.nlargest(n_top).index
    mat_subset = mat.loc[top_genes]

    logger.info(
        f"Computing PCA on {n_top} most variable genes across {len(common_samples)} samples..."
    )

    # Samples as rows, genes as columns
    x_matrix = mat_subset.T.values

    pca = PCA(n_components=n_components)
    coords = pca.fit_transform(x_matrix)

    col_names = [f"PC{i + 1}" for i in range(n_components)]
    pca_df = pd.DataFrame(coords, index=common_samples, columns=col_names)

    # Attach metadata columns
    for col in meta_aligned.columns:
        pca_df[col] = meta_aligned[col]

    return pca, pca_df


def plot_pca(
    data: Union[pd.DataFrame, NormalizationResult],
    metadata: pd.DataFrame,
    color_by: str,
    shape_by: Optional[str] = None,
    n_components: int = 2,
    top_n_variable_genes: int = 500,
    interactive: bool = False,
    show_sample_labels: bool = True,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (7.0, 5.5),
    ax: Optional[plt.Axes] = None,
    **kwargs,
) -> Any:
    """
    Generate static or interactive PCA plot with variance explained per component.

    Parameters
    ----------
    data : pd.DataFrame or NormalizationResult
        Gene expression counts (e.g. log2-normalized).
    metadata : pd.DataFrame
        Sample annotation table.
    color_by : str
        Metadata column to color sample points.
    shape_by : str, optional
        Metadata column for marker shapes.
    n_components : int, default=2
        Number of principal components (2 or 3).
    top_n_variable_genes : int, default=500
        Number of variable genes to select.
    interactive : bool, default=False
        If True, returns a Plotly Figure.
    show_sample_labels : bool, default=True
        Whether to print text sample IDs next to dots.
    title : str, optional
        Plot title.
    figsize : tuple, default=(7, 5.5)
        Matplotlib figure size.
    ax : plt.Axes, optional
        Pre-existing axes.

    Returns
    -------
    tuple of (plt.Figure, plt.Axes, pd.DataFrame) or (plotly.graph_objects.Figure, pd.DataFrame)
    """
    pca, pca_df = compute_pca(
        data=data,
        metadata=metadata,
        n_components=n_components,
        top_n_variable_genes=top_n_variable_genes,
    )

    var_ratios = pca.explained_variance_ratio_ * 100
    pc1_lbl = f"PC1 ({var_ratios[0]:.1f}%)"
    pc2_lbl = f"PC2 ({var_ratios[1]:.1f}%)"

    default_title = (
        title
        if title is not None
        else f"PCA Plot (Top {min(top_n_variable_genes, len(data))} Genes)"
    )

    # -------------------------------------------------------------
    # Interactive Plotly Mode
    # -------------------------------------------------------------
    if interactive:
        if not HAS_PLOTLY:
            raise ImportError(
                "plotly is required for interactive PCA plots. Run pip install plotly."
            )

        if n_components >= 3:
            pc3_lbl = f"PC3 ({var_ratios[2]:.1f}%)"
            fig = px.scatter_3d(
                pca_df,
                x="PC1",
                y="PC2",
                z="PC3",
                color=color_by,
                symbol=shape_by,
                hover_name=pca_df.index,
                title=default_title,
                labels={"PC1": pc1_lbl, "PC2": pc2_lbl, "PC3": pc3_lbl},
                template="plotly_white",
            )
        else:
            fig = px.scatter(
                pca_df,
                x="PC1",
                y="PC2",
                color=color_by,
                symbol=shape_by,
                hover_name=pca_df.index,
                text=pca_df.index if show_sample_labels else None,
                title=default_title,
                labels={"PC1": pc1_lbl, "PC2": pc2_lbl},
                template="plotly_white",
            )
            if show_sample_labels:
                fig.update_traces(textposition="top center")

        return fig, pca_df

    # -------------------------------------------------------------
    # Static Matplotlib Mode
    # -------------------------------------------------------------
    set_publication_style()

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    sns.scatterplot(
        data=pca_df,
        x="PC1",
        y="PC2",
        hue=color_by,
        style=shape_by,
        s=90,
        palette=NATURE_PALETTE[: len(pca_df[color_by].unique())],
        ax=ax,
        edgecolor="#222222",
        linewidth=0.8,
    )

    if show_sample_labels:
        for sample_id, row in pca_df.iterrows():
            ax.annotate(
                str(sample_id),
                (row["PC1"], row["PC2"]),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=8,
                alpha=0.85,
            )

    ax.set_xlabel(pc1_lbl)
    ax.set_ylabel(pc2_lbl)
    ax.set_title(default_title)
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False)

    return fig, ax, pca_df
