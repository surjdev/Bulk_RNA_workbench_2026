"""
Downstream Bulk Transcriptomics Workbench (`btw`)
A modular Python utility layer for reproducible downstream bulk RNA-seq analysis.
"""

from __future__ import annotations

import os
import random
import sys
from typing import Optional

import numpy as np
__version__ = "0.1.0"

try:
    from loguru import logger
    # Configure default logger formatting
    logger.remove()
    logger.add(
        sys.stderr,
        format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
        level="INFO",
    )
except ImportError:
    import logging
    logger = logging.getLogger("btw")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)-8s | %(name)s - %(message)s",
    )


def set_seed(seed: int = 42) -> int:
    """
    Set random seeds across Python, NumPy, and optional frameworks from a single point
    to ensure reproducible results across the analysis workflow.

    Parameters
    ----------
    seed : int, default=42
        Integer seed value.

    Returns
    -------
    int
        The seed that was set.
    """
    random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    np.random.seed(seed)

    # If PyTorch is installed (e.g. used by decoupler), set seed as well
    if "torch" in sys.modules:
        import torch  # type: ignore

        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)

    logger.info(f"Central random seed set to {seed}")
    return seed


__all__ = ["__version__", "logger", "set_seed"]
