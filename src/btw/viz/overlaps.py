"""
Cross-contrast DEG overlap visualization module (FR-4).
Provides UpSet plots and Venn diagrams for multi-contrast comparisons using upsetplot.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Set, Tuple

import matplotlib.pyplot as plt
import pandas as pd
from btw import logger
from btw.de_analysis.contrasts import MultiContrastResult
from btw.viz.style import set_publication_style

try:
    from upsetplot import UpSet, from_contents
    HAS_UPSETPLOT = True
except ImportError:
    HAS_UPSETPLOT = False

try:
    from matplotlib_venn import venn2, venn3
    HAS_VENN = True
except ImportError:
    HAS_VENN = False


def _extract_deg_sets(
    source: Any,
    padj_cutoff: float = 0.05,
    lfc_cutoff: float = 1.0,
    direction: str = "both",
) -> Dict[str, Set[str]]:
    """Extract dictionary of gene sets from MultiContrastResult or raw dictionary."""
    if isinstance(source, MultiContrastResult):
        return source.get_deg_sets(padj_cutoff=padj_cutoff, lfc_cutoff=lfc_cutoff, direction=direction)
    elif isinstance(source, dict):
        return {str(k): set(v) for k, v in source.items()}
    else:
        raise TypeError(f"Expected MultiContrastResult or dict of sets, got {type(source).__name__}")


def plot_upset(
    source: Any,
    padj_cutoff: float = 0.05,
    lfc_cutoff: float = 1.0,
    direction: str = "both",
    min_subset_size: int = 0,
    sort_by: str = "cardinality",
    show_counts: bool = True,
    figsize: Tuple[float, float] = (8.0, 5.0),
    title: Optional[str] = None,
) -> Any:
    """
    Generate an UpSet plot visualizing intersections among multiple DEG sets.

    Parameters
    ----------
    source : MultiContrastResult or dict of {str: set of str}
        Source containing DEG sets per contrast.
    padj_cutoff : float, default=0.05
        Significance cutoff.
    lfc_cutoff : float, default=1.0
        Log2FC threshold.
    direction : {'both', 'up', 'down'}, default='both'
        Regulation direction.
    min_subset_size : int, default=0
        Minimum size of intersection subset to display.
    sort_by : {'cardinality', 'degree'}, default='cardinality'
        Ordering of intersection subsets.
    show_counts : bool, default=True
        Whether to show numbers atop subset bars.
    figsize : tuple, default=(8, 5)
        Plot dimensions.
    title : str, optional
        Plot title.

    Returns
    -------
    dict
        UpSet plot axes dictionary.
    """
    if not HAS_UPSETPLOT:
        raise ImportError("upsetplot is required for UpSet plots. Run pip install upsetplot.")

    deg_sets = _extract_deg_sets(source, padj_cutoff=padj_cutoff, lfc_cutoff=lfc_cutoff, direction=direction)

    if not deg_sets:
        raise ValueError("No gene sets provided to plot_upset.")

    # Check that at least one set is non-empty
    total_genes = sum(len(s) for s in deg_sets.values())
    if total_genes == 0:
        logger.warning("All DEG sets are empty. UpSet plot cannot be constructed.")

    # Convert sets into upsetplot Series format
    upset_data = from_contents(deg_sets)

    set_publication_style()
    fig = plt.figure(figsize=figsize)

    upset = UpSet(
        upset_data,
        min_subset_size=min_subset_size,
        sort_by=sort_by,
        show_counts=show_counts,
        facecolor="#3C5488",
    )
    axes_dict = upset.plot(fig=fig)

    if title is not None:
        fig.suptitle(title, fontsize=12, fontweight="bold", y=0.98)

    return axes_dict


def plot_venn(
    source: Any,
    padj_cutoff: float = 0.05,
    lfc_cutoff: float = 1.0,
    direction: str = "both",
    figsize: Tuple[float, float] = (6.0, 5.0),
    title: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Generate a 2-way or 3-way Venn diagram for DEG comparisons.

    Parameters
    ----------
    source : MultiContrastResult or dict of {str: set of str}
        Source containing DEG sets.
    padj_cutoff : float, default=0.05
        Cutoff for padj.
    lfc_cutoff : float, default=1.0
        Cutoff for log2FC.
    direction : {'both', 'up', 'down'}, default='both'
        Direction filter.
    figsize : tuple, default=(6, 5)
        Figure size.
    title : str, optional
        Title.
    ax : plt.Axes, optional
        Pre-existing axes.

    Returns
    -------
    tuple of (plt.Figure, plt.Axes)
    """
    deg_sets = _extract_deg_sets(source, padj_cutoff=padj_cutoff, lfc_cutoff=lfc_cutoff, direction=direction)
    names = list(deg_sets.keys())
    n_sets = len(names)

    set_publication_style()
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    if HAS_VENN:
        if n_sets == 2:
            venn2(
                subsets=[deg_sets[names[0]], deg_sets[names[1]]],
                set_labels=(names[0], names[1]),
                ax=ax,
                set_colors=("#E64B35", "#4DBBD5"),
            )
        elif n_sets == 3:
            venn3(
                subsets=[deg_sets[names[0]], deg_sets[names[1]], deg_sets[names[2]]],
                set_labels=(names[0], names[1], names[2]),
                ax=ax,
                set_colors=("#E64B35", "#4DBBD5", "#00A087"),
            )
        else:
            raise ValueError(f"Venn diagram supports 2 or 3 sets, but got {n_sets}. Use plot_upset() instead.")
    else:
        # Fallback table visualization if matplotlib-venn not installed
        rows = []
        for i, name_i in enumerate(names):
            for j, name_j in enumerate(names):
                if i <= j:
                    overlap = len(deg_sets[name_i] & deg_sets[name_j])
                    rows.append({"Set A": name_i, "Set B": name_j, "Shared Genes": overlap})
        overlap_df = pd.DataFrame(rows)
        ax.axis("off")
        ax.table(
            cellText=overlap_df.values,
            colLabels=overlap_df.columns,
            loc="center",
            cellLoc="center",
        )
        ax.set_title("DEG Overlap Summary Table")

    if title is not None:
        ax.set_title(title, fontsize=12, fontweight="bold")

    return fig, ax
