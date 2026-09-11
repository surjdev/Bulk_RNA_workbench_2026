"""
Central R-Python Interoperability Bridge for BTW (FR-10).
Provides safe rpy2 execution, bidirectional pandas <-> R conversions,
package availability checks, and informative error diagnostics.
"""

from __future__ import annotations

from typing import Any, Optional

import pandas as pd

from btw import logger

# Lazy / Safe rpy2 import
HAS_RPY2 = False
_R_AVAILABLE = False
_r = None
_robjects = None
_pandas2ri = None

try:
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    HAS_RPY2 = True
    _robjects = ro
    _r = ro.r
    _pandas2ri = pandas2ri
    # Verify R session is responsive
    _ = ro.r("R.version.string")[0]
    _R_AVAILABLE = True
except (ImportError, Exception):
    HAS_RPY2 = False
    _R_AVAILABLE = False


class RInteropError(RuntimeError):
    """Base exception for R interoperability issues."""

    pass


class MissingRPackageError(RInteropError):
    """
    Raised when an R/Bioconductor package required for an analysis engine is not installed.
    Provides clear, human-readable installation instructions.
    """

    def __init__(self, package_name: str, purpose: str = "", alternative: str = ""):
        self.package_name = package_name
        self.purpose = purpose
        self.alternative = alternative

        msg = (
            f"\n{'=' * 70}\n"
            f" [R Dependency Error] R package '{package_name}' is not installed.\n"
            f"{'=' * 70}\n"
        )
        if purpose:
            msg += f"Purpose: {purpose}\n\n"
        msg += (
            f"To install this package in your R environment, run:\n"
            f"  R -e \"if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager'); "
            f"BiocManager::install('{package_name}')\"\n"
        )
        if alternative:
            msg += f"\nAlternative: {alternative}\n"
        msg += f"{'=' * 70}\n"
        super().__init__(msg)


def is_r_available() -> bool:
    """Return True if both rpy2 and an active R installation are accessible."""
    return _R_AVAILABLE


def get_r_version() -> Optional[str]:
    """Retrieve installed R version string, or None if R is unavailable."""
    if not is_r_available():
        return None
    try:
        return str(_r("R.version.string")[0])
    except Exception:
        return None


def check_r_package(package_name: str) -> bool:
    """
    Check whether a specific R package is installed in the current R library path.

    Parameters
    ----------
    package_name : str
        Name of R package (e.g. 'DESeq2', 'limma', 'WGCNA', 'sva', 'clusterProfiler').

    Returns
    -------
    bool
        True if package is available and loadable.
    """
    if not is_r_available():
        return False
    try:
        res = _r(f"nzchar(system.file(package='{package_name}'))")[0]
        return bool(res)
    except Exception:
        return False


def get_r_package_version(package_name: str) -> Optional[str]:
    """Retrieve version string of an installed R package, or None."""
    if not check_r_package(package_name):
        return None
    try:
        ver = _r(f"as.character(packageVersion('{package_name}'))")[0]
        return str(ver)
    except Exception:
        return None


def require_r_package(package_name: str, purpose: str = "", alternative: str = "") -> None:
    """
    Assert that an R package is installed; raises MissingRPackageError if absent.

    Parameters
    ----------
    package_name : str
        R package name.
    purpose : str, optional
        Explanation of why the package is needed.
    alternative : str, optional
        Alternative configuration advice (e.g. 'set engine=\"python\"').
    """
    if not is_r_available():
        raise RInteropError(
            "R or rpy2 is not accessible. Please ensure R 4.x and rpy2 are installed in your environment."
        )
    if not check_r_package(package_name):
        raise MissingRPackageError(package_name, purpose=purpose, alternative=alternative)


def set_r_seed(seed: int) -> None:
    """Set the random seed in the active R session."""
    if is_r_available():
        try:
            _r(f"set.seed({int(seed)})")
            logger.debug(f"R random seed synchronized to {seed}")
        except Exception as e:
            logger.warning(f"Failed to synchronize R seed: {e}")


def pandas_to_r_df(df: pd.DataFrame) -> Any:
    """
    Convert a pandas DataFrame to an R data.frame with row names preserved.

    Parameters
    ----------
    df : pd.DataFrame
        Input table.

    Returns
    -------
    rpy2.robjects.vectors.DataFrame
        Equivalent R data.frame.
    """
    if not is_r_available():
        raise RInteropError("rpy2 / R session is not available.")

    # Convert columns to standard types
    clean_df = df.copy()
    row_names = [str(x) for x in clean_df.index]

    with (ro.default_converter + pandas2ri.converter).context():
        r_df = ro.conversion.get_conversion().py2rpy(clean_df)

    # Explicitly assign row names in R
    r_df = _r["rownames<-"](r_df, ro.StrVector(row_names))
    return r_df


def pandas_to_r_matrix(df: pd.DataFrame, is_integer: bool = False) -> Any:
    """
    Convert a numeric pandas DataFrame into an R matrix with rownames and colnames.

    Parameters
    ----------
    df : pd.DataFrame
        Numeric expression table (e.g. genes x samples).
    is_integer : bool, default=False
        Whether matrix should be typed as integer.

    Returns
    -------
    rpy2.robjects.vectors.Matrix
        R numeric/integer matrix.
    """
    if not is_r_available():
        raise RInteropError("rpy2 / R session is not available.")

    clean_df = df.copy()
    row_names = ro.StrVector([str(x) for x in clean_df.index])
    col_names = ro.StrVector([str(x) for x in clean_df.columns])

    values = clean_df.values.flatten(order="F")  # Fortran order for R matrix
    if is_integer:
        r_vec = ro.IntVector([int(v) for v in values])
    else:
        r_vec = ro.FloatVector([float(v) for v in values])

    dimnames = _r["list"](row_names, col_names)
    r_mat = _r["matrix"](r_vec, nrow=clean_df.shape[0], ncol=clean_df.shape[1], dimnames=dimnames)
    return r_mat


def r_to_pandas_df(r_obj: Any) -> pd.DataFrame:
    """
    Convert an R data.frame or DataFrame to a pandas DataFrame with index intact.

    Parameters
    ----------
    r_obj : R object
        R data.frame, DataFrame, or MArrayLM output.

    Returns
    -------
    pd.DataFrame
        Equivalent pandas DataFrame.
    """
    if not is_r_available():
        raise RInteropError("rpy2 / R session is not available.")

    # Check if object has rownames
    has_rownames = False
    row_names_list = []
    try:
        r_rownames = _r("rownames")(r_obj)
        if r_rownames is not ro.NULL and len(r_rownames) > 0:
            has_rownames = True
            row_names_list = list(r_rownames)
    except Exception:
        pass

    with (ro.default_converter + pandas2ri.converter).context():
        py_df = ro.conversion.get_conversion().rpy2py(r_obj)

    if not isinstance(py_df, pd.DataFrame):
        py_df = pd.DataFrame(py_df)

    if has_rownames and len(row_names_list) == len(py_df):
        py_df.index = row_names_list

    return py_df
