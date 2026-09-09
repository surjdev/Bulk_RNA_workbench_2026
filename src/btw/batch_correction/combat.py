"""
ComBat and ComBat-Seq Batch Effect Correction utilities for BTW (FR-7).
Integrates inmoose.pycombat for removing technical batch variations while preserving biology.
"""

from __future__ import annotations

from typing import List, Optional, Union

import numpy as np
import pandas as pd
from btw import logger

try:
    import inmoose.pycombat as pc

    HAS_INMOOSE = True
except ImportError:
    HAS_INMOOSE = False


def run_combat(
    data: pd.DataFrame,
    batch: Union[str, List[str], pd.Series],
    metadata: Optional[pd.DataFrame] = None,
    biological_factor: Optional[Union[str, List[str], pd.Series]] = None,
    is_count: Optional[bool] = None,
    **combat_kwargs,
) -> pd.DataFrame:
    """
    Apply ComBat (for continuous/normalized data) or ComBat-Seq (for integer RNA-Seq count data)
    to adjust for technical batch effects.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix. Supports (genes x samples) or (samples x genes).
    batch : str, list, or Series
        Batch indicator per sample. If string and metadata is provided, extracted from metadata[batch].
    metadata : pd.DataFrame, optional
        Sample metadata containing batch and biological covariates.
    biological_factor : str, list, or Series, optional
        Biological condition to protect during adjustment (e.g. treatment vs control).
    is_count : bool, optional
        If True, runs ComBat-Seq (negative binomial GLM). If False, runs parametric ComBat.
        If None, automatically detected based on whether data contains integers.
    **combat_kwargs
        Additional arguments passed to inmoose.pycombat_seq or inmoose.pycombat_norm.

    Returns
    -------
    pd.DataFrame
        Batch-corrected expression matrix matching the input shape and index orientation.
    """
    if not HAS_INMOOSE:
        raise ImportError("inmoose package is required for ComBat batch correction. Run pip install inmoose.")

    df = data.copy()

    # Determine orientation: ComBat expects (genes x samples)
    transposed = False
    if metadata is not None:
        if set(metadata.index).issubset(df.index) and not set(metadata.index).issubset(df.columns):
            # Input is (samples x genes)
            df = df.T
            transposed = True
            logger.info("Transposed input matrix to (genes x samples) for ComBat processing.")

    samples = list(df.columns)
    n_genes, n_samples = df.shape

    # Resolve batch vector
    if isinstance(batch, str):
        if metadata is None or batch not in metadata.columns:
            raise ValueError(f"Batch column '{batch}' not found in metadata.")
        batch_series = metadata.loc[samples, batch]
        batch_vec = batch_series.tolist()
    elif isinstance(batch, pd.Series):
        batch_vec = batch.loc[samples].tolist() if set(samples).issubset(batch.index) else batch.tolist()
    else:
        batch_vec = list(batch)

    if len(batch_vec) != n_samples:
        raise ValueError(f"Batch vector length ({len(batch_vec)}) does not match number of samples ({n_samples}).")

    # Resolve biological covariate model
    covar_mod = None
    if biological_factor is not None:
        if isinstance(biological_factor, str):
            if metadata is None or biological_factor not in metadata.columns:
                raise ValueError(f"Biological factor '{biological_factor}' not found in metadata.")
            covar_series = metadata.loc[samples, biological_factor]
            covar_mod = pd.DataFrame({biological_factor: covar_series.values})
        elif isinstance(biological_factor, (pd.Series, pd.DataFrame)):
            covar_mod = biological_factor
        elif isinstance(biological_factor, list):
            covar_mod = pd.DataFrame({"covariate": biological_factor})

    # Auto-detect count data if not specified
    if is_count is None:
        # Sample non-null values to test if they are integer counts
        sample_vals = df.iloc[: min(20, n_genes), : min(10, n_samples)].values
        is_count = bool(np.all(np.equal(np.mod(sample_vals, 1), 0)) and np.all(sample_vals >= 0))
        logger.info(f"Auto-detected data type: {'integer counts (ComBat-Seq)' if is_count else 'normalized/continuous (ComBat)'}")

    if is_count:
        logger.info(f"Executing ComBat-Seq on {n_genes} genes across {n_samples} samples...")
        corrected = pc.pycombat_seq(
            counts=df,
            batch=batch_vec,
            covar_mod=covar_mod,
            **combat_kwargs,
        )
    else:
        logger.info(f"Executing ComBat on {n_genes} genes across {n_samples} samples...")
        corrected = pc.pycombat_norm(
            counts=df,
            batch=batch_vec,
            covar_mod=covar_mod,
            **combat_kwargs,
        )

    if not isinstance(corrected, pd.DataFrame):
        corrected = pd.DataFrame(corrected, index=df.index, columns=df.columns)
    else:
        corrected.index = df.index
        corrected.columns = df.columns

    if transposed:
        corrected = corrected.T

    logger.info("ComBat batch correction successfully completed.")
    return corrected
