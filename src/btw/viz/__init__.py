"""
Publication-grade biological visualization subpackage for BTW (FR-4).
"""

from btw.viz.dispersion import plot_dispersion
from btw.viz.heatmap import plot_heatmap
from btw.viz.ma import plot_ma
from btw.viz.overlaps import plot_upset, plot_venn
from btw.viz.pca import compute_pca, plot_pca
from btw.viz.style import (
    CELL_PALETTE,
    NATURE_PALETTE,
    REGULATION_COLORS,
    SCIENCE_PALETTE,
    save_figure,
    set_publication_style,
)
from btw.viz.volcano import plot_volcano

__all__ = [
    "set_publication_style",
    "save_figure",
    "NATURE_PALETTE",
    "CELL_PALETTE",
    "SCIENCE_PALETTE",
    "REGULATION_COLORS",
    "plot_volcano",
    "compute_pca",
    "plot_pca",
    "plot_ma",
    "plot_dispersion",
    "plot_heatmap",
    "plot_upset",
    "plot_venn",
]
