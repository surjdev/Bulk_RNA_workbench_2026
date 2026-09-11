"""
Caching infrastructure for computationally heavy enrichment and annotation tasks (NFR-Performance).
Uses joblib.Memory to avoid redundant API requests and repeated permutation tests.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

from joblib import Memory

DEFAULT_CACHE_DIR = Path(".btw_cache")


def get_memory_cache(
    cache_dir: Optional[Union[str, Path]] = None,
    enabled: bool = True,
    verbose: int = 0,
) -> Memory:
    """
    Instantiate or retrieve a joblib.Memory cache manager.

    Parameters
    ----------
    cache_dir : str or Path, optional
        Target directory for persistent cache files. Defaults to '.btw_cache'.
    enabled : bool, default=True
        If False, creates an inactive cache (location=None).
    verbose : int, default=0
        Verbosity level for joblib cache hits/misses.

    Returns
    -------
    Memory
        joblib Memory instance.
    """
    if not enabled:
        return Memory(location=None, verbose=verbose)

    target_dir = Path(cache_dir) if cache_dir is not None else DEFAULT_CACHE_DIR
    target_dir.mkdir(parents=True, exist_ok=True)
    return Memory(location=str(target_dir), verbose=verbose)
