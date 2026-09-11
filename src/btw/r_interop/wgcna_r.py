"""
R WGCNA wrapper through rpy2 for BTW (FR-8 & FR-10).
Serves as the reference implementation engine for co-expression network analysis.
"""

from __future__ import annotations

from typing import Any, Tuple

import pandas as pd

from btw import logger
from btw.r_interop.bridge import (
    get_r_package_version,
    is_r_available,
    pandas_to_r_matrix,
    r_to_pandas_df,
    require_r_package,
)

if is_r_available():
    import rpy2.robjects as ro

    _r = ro.r
else:
    _r = None


def run_r_wgcna(
    data: pd.DataFrame,
    power: int = 6,
    min_module_size: int = 10,
    network_type: str = "unsigned",
    merge_cut_height: float = 0.25,
) -> Tuple[pd.Series, pd.DataFrame, Any]:
    """
    Execute original R WGCNA blockwiseModules algorithm via rpy2.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix of shape (samples x genes) or (genes x samples).
    power : int, default=6
        Soft-thresholding power beta.
    min_module_size : int, default=10
        Minimum cluster size for module detection.
    network_type : str, default='unsigned'
        'unsigned', 'signed', or 'signed hybrid'.
    merge_cut_height : float, default=0.25
        Dendrogram cut height for module merging.

    Returns
    -------
    module_labels : pd.Series
        Mapping of gene identifiers to module color strings.
    module_eigengenes : pd.DataFrame
        Module eigengenes across samples.
    r_net : R object
        Raw R blockwiseModules result list.
    """
    require_r_package(
        "WGCNA",
        purpose="Reference Co-expression Network analysis (FR-8)",
        alternative="Use engine='python' for Python co-expression helper.",
    )

    # Orientation: WGCNA in R strictly expects (samples x genes)
    if data.shape[0] > data.shape[1] and data.shape[1] <= 100:
        # Typical bulk RNA: 10,000+ genes x 6 samples -> transpose
        expr = data.T.copy()
    else:
        expr = data.copy()

    logger.info(
        f"Executing R WGCNA (v{get_r_package_version('WGCNA')}) on {expr.shape[1]} genes across "
        f"{expr.shape[0]} samples (power={power}, minModuleSize={min_module_size})..."
    )

    r_expr = pandas_to_r_matrix(expr)
    _r.assign(".wgcna_expr", r_expr)

    # Disable WGCNA multi-threading in R session to avoid deadlocks
    _r("options(stringsAsFactors = FALSE)")
    _r("WGCNA::disableWGCNAThreads()")

    cmd = (
        f"WGCNA::blockwiseModules("
        f"  .wgcna_expr, "
        f"  power={power}, "
        f"  networkType='{network_type}', "
        f"  minModuleSize={min_module_size}, "
        f"  mergeCutHeight={merge_cut_height}, "
        f"  numericLabels=FALSE, "
        f"  pamRespectsDendro=FALSE, "
        f"  saveTOMs=FALSE, "
        f"  verbose=0"
        f")"
    )
    r_net = _r(cmd)

    # Extract colors (genes -> module color)
    r_colors = _r("as.character(.wgcna_net$colors)")
    module_labels = pd.Series(list(r_colors), index=expr.columns, name="module")

    # Extract MEs (samples x MEs)
    r_mes = _r(".wgcna_net$MEs")
    mes_df = r_to_pandas_df(r_mes)
    mes_df.index = expr.index

    return module_labels, mes_df, r_net
