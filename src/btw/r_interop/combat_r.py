"""
R sva::ComBat wrapper through rpy2 for BTW (FR-7 & FR-10).
Serves as the reference implementation engine for batch effect correction.
"""

from __future__ import annotations

from typing import Any, List, Optional, Union

import numpy as np
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


def run_r_combat(
    data: pd.DataFrame,
    batch: Union[str, List[Any], pd.Series],
    metadata: Optional[pd.DataFrame] = None,
    biological_factor: Optional[Union[str, List[Any], pd.Series]] = None,
    is_count: Optional[bool] = None,
) -> pd.DataFrame:
    """
    Execute original R Bioconductor sva::ComBat or sva::ComBat_seq via rpy2.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix (genes x samples).
    batch : str, list, or Series
        Batch indicator per sample.
    metadata : pd.DataFrame, optional
        Sample metadata table.
    biological_factor : str, list, or Series, optional
        Biological condition to protect.
    is_count : bool, optional
        If True, calls sva::ComBat_seq; if False, calls sva::ComBat.
        If None, inferred from data dtype.

    Returns
    -------
    pd.DataFrame
        Batch-adjusted expression matrix.
    """
    require_r_package(
        "sva",
        purpose="Reference ComBat batch correction (FR-7)",
        alternative="Use engine='python' for inmoose.pycombat.",
    )

    df = data.copy()

    # Resolve batch series
    if isinstance(batch, str):
        if metadata is None or batch not in metadata.columns:
            raise ValueError(f"Batch column '{batch}' not found in metadata.")
        batch_series = metadata.loc[df.columns, batch]
    elif isinstance(batch, pd.Series):
        batch_series = batch.loc[df.columns]
    else:
        batch_series = pd.Series(batch, index=df.columns)

    # Incur is_count if not specified
    if is_count is None:
        is_count = bool(np.all(np.equal(np.mod(df.values, 1), 0)) and (df.values.min() >= 0))

    logger.info(
        f"Executing R sva::ComBat (v{get_r_package_version('sva')}, is_count={is_count}) on "
        f"{df.shape[0]} genes across {df.shape[1]} samples in {batch_series.nunique()} batches..."
    )

    r_dat = pandas_to_r_matrix(df, is_integer=is_count)
    r_batch = ro.StrVector([str(b) for b in batch_series])

    _r.assign(".combat_dat", r_dat)
    _r.assign(".combat_batch", r_batch)

    # Biological covariate model matrix
    if biological_factor is not None:
        if isinstance(biological_factor, str):
            if metadata is None or biological_factor not in metadata.columns:
                raise ValueError(f"Biological column '{biological_factor}' not found in metadata.")
            bio_series = metadata.loc[df.columns, biological_factor]
        elif isinstance(biological_factor, pd.Series):
            bio_series = biological_factor.loc[df.columns]
        else:
            bio_series = pd.Series(biological_factor, index=df.columns)

        _r.assign(".combat_bio", ro.StrVector([str(b) for b in bio_series]))
        _r(".combat_mod <- model.matrix(~ as.factor(.combat_bio))")
        mod_arg = ".combat_mod"
    else:
        mod_arg = "NULL"

    ro.packages.importr("sva")

    if is_count:
        if biological_factor is not None:
            r_res = _r(
                "sva::ComBat_seq(counts=.combat_dat, batch=.combat_batch, group=as.factor(.combat_bio))"
            )
        else:
            r_res = _r("sva::ComBat_seq(counts=.combat_dat, batch=.combat_batch)")
    else:
        r_res = _r(f"sva::ComBat(dat=.combat_dat, batch=.combat_batch, mod={mod_arg})")

    adj_df = r_to_pandas_df(r_res)
    adj_df.index = df.index
    adj_df.columns = df.columns
    return adj_df
