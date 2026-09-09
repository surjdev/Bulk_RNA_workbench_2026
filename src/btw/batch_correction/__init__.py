"""
Batch Effect Correction module for BTW (FR-7).
"""

from btw.batch_correction.combat import run_combat
from btw.batch_correction.diagnostic import compare_pca_batch, evaluate_batch_effect

__all__ = [
    "run_combat",
    "compare_pca_batch",
    "evaluate_batch_effect",
]
