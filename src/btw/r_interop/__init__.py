"""
R-Python Interoperability subpackage for BTW (FR-10).
Connects Python analysis workflows to reference R/Bioconductor engines.
"""

from btw.r_interop.bridge import (
    MissingRPackageError,
    RInteropError,
    check_r_package,
    get_r_package_version,
    get_r_version,
    is_r_available,
    pandas_to_r_df,
    pandas_to_r_matrix,
    r_to_pandas_df,
    require_r_package,
    set_r_seed,
)
from btw.r_interop.clusterprofiler_r import run_r_clusterprofiler
from btw.r_interop.combat_r import run_r_combat
from btw.r_interop.deseq2_r import run_r_deseq2, run_r_vst
from btw.r_interop.limma_r import run_r_limma
from btw.r_interop.wgcna_r import run_r_wgcna

__all__ = [
    "is_r_available",
    "check_r_package",
    "require_r_package",
    "get_r_version",
    "get_r_package_version",
    "set_r_seed",
    "pandas_to_r_df",
    "pandas_to_r_matrix",
    "r_to_pandas_df",
    "RInteropError",
    "MissingRPackageError",
    "run_r_deseq2",
    "run_r_vst",
    "run_r_limma",
    "run_r_combat",
    "run_r_wgcna",
    "run_r_clusterprofiler",
]
