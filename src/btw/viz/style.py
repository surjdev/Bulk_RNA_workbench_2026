"""
Publication-grade visual styling and curated color palettes for Bulk Transcriptomics (FR-4).
Inspired by leading journals (Nature, Cell, Science).
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Union

import matplotlib.pyplot as plt
import seaborn as sns
from btw import logger

# Curated publication palettes
NATURE_PALETTE: List[str] = [
    "#E64B35",  # Coral red
    "#4DBBD5",  # Cyan
    "#00A087",  # Teal green
    "#3C5488",  # Navy
    "#F39B7F",  # Peach
    "#8491B4",  # Slate
    "#91D1C2",  # Mint
    "#DC0000",  # Crimson
    "#7E6148",  # Warm brown
    "#B09C85",  # Sand
]

CELL_PALETTE: List[str] = [
    "#1F77B4",  # Blue
    "#FF7F0E",  # Orange
    "#2CA02C",  # Green
    "#D62728",  # Red
    "#9467BD",  # Purple
    "#8C564B",  # Brown
    "#E377C2",  # Pink
    "#7F7F7F",  # Gray
    "#BCBD22",  # Olive
    "#17BECF",  # Sky blue
]

SCIENCE_PALETTE: List[str] = [
    "#0C5DA5",  # Deep blue
    "#00B945",  # Bright green
    "#FF9500",  # Bright orange
    "#FF2C00",  # Red
    "#845B97",  # Violet
    "#474747",  # Charcoal
    "#9E9E9E",  # Light gray
]

REGULATION_COLORS: Dict[str, str] = {
    "UP": "#E64B35",    # Vivid red
    "DOWN": "#4DBBD5",  # Soft blue
    "NS": "#B0B0B0",    # Neutral gray
}

PALETTES: Dict[str, List[str]] = {
    "nature": NATURE_PALETTE,
    "cell": CELL_PALETTE,
    "science": SCIENCE_PALETTE,
}


def set_publication_style(
    palette: str = "nature",
    dpi: int = 300,
    font_family: str = "sans-serif",
) -> None:
    """
    Configure matplotlib and seaborn global parameters for publication-quality figures.

    Parameters
    ----------
    palette : {'nature', 'cell', 'science'}, default='nature'
        Color palette theme.
    dpi : int, default=300
        Resolution for raster outputs.
    font_family : str, default='sans-serif'
        Font family name.
    """
    sns.set_theme(style="ticks")

    palettes = {
        "nature": NATURE_PALETTE,
        "cell": CELL_PALETTE,
        "science": SCIENCE_PALETTE,
    }
    chosen_palette = palettes.get(palette.lower(), NATURE_PALETTE)
    sns.set_palette(chosen_palette)

    plt.rcParams.update(
        {
            "figure.dpi": dpi,
            "savefig.dpi": dpi,
            "font.family": font_family,
            "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.titleweight": "bold",
            "axes.labelsize": 10,
            "axes.labelweight": "medium",
            "axes.edgecolor": "#222222",
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.frameon": False,
            "legend.fontsize": 9,
        }
    )
    logger.debug(f"Applied publication styling: palette={palette}, dpi={dpi}")


def save_figure(
    fig: plt.Figure,
    output_path: Union[str, Path],
    dpi: int = 300,
    transparent: bool = False,
    bbox_inches: str = "tight",
    **kwargs,
) -> Path:
    """
    Save figure with publication-grade resolution (PNG, PDF, SVG).

    Parameters
    ----------
    fig : matplotlib Figure
        Figure object to save.
    output_path : str or Path
        Destination path.
    dpi : int, default=300
        Image resolution.
    transparent : bool, default=False
        Whether background is transparent.
    bbox_inches : str, default='tight'
        Bounding box cropping.

    Returns
    -------
    Path
        Absolute path to saved figure.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        path,
        dpi=dpi,
        transparent=transparent,
        bbox_inches=bbox_inches,
        **kwargs,
    )
    logger.info(f"Saved publication figure ({path.suffix}) to {path.resolve()}")
    return path.resolve()
