"""
Differential Expression (DE) Analysis subpackage for BTW (FR-3).
"""

from btw.de_analysis.contrasts import MultiContrastResult, run_multiple_contrasts
from btw.de_analysis.correction import (
    SUPPORTED_METHODS,
    adjust_pvalues,
    apply_multitest_correction,
)
from btw.de_analysis.deseq_helper import (
    DEResult,
    build_deseq_dataset,
    run_de,
    run_deseq_stats,
)

__all__ = [
    "build_deseq_dataset",
    "run_deseq_stats",
    "run_de",
    "DEResult",
    "run_multiple_contrasts",
    "MultiContrastResult",
    "adjust_pvalues",
    "apply_multitest_correction",
    "SUPPORTED_METHODS",
]
