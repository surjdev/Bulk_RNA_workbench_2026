"""
Gene ID and symbol mapping utilities for BTW (FR-6).
Provides integration with MyGene.info with automatic local caching.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Union

import pandas as pd

from btw import logger
from btw.de_analysis.contrasts import DEResult

try:
    import mygene

    HAS_MYGENE = True
except ImportError:
    HAS_MYGENE = False


def map_gene_ids(
    gene_list: List[str],
    from_id: str = "ensembl.gene",
    to_id: str = "symbol",
    species: str = "human",
    cache_dir: Optional[Union[str, Path]] = None,
    use_cache: bool = True,
    mapping_dict: Optional[Dict[str, str]] = None,
    **query_kwargs,
) -> pd.DataFrame:
    """
    Map a list of gene identifiers (e.g. Ensembl IDs) to symbols, Entrez, or names using MyGene.info.

    Parameters
    ----------
    gene_list : list of str
        List of input gene IDs.
    from_id : str, default='ensembl.gene'
        Source identifier namespace ('ensembl.gene', 'symbol', 'entrezgene').
    to_id : str, default='symbol'
        Target identifier field ('symbol', 'name', 'entrezgene').
    species : str, default='human'
        Target organism ('human', 'mouse', 'rat', 9606, 10090).
    cache_dir : str or Path, optional
        Directory for local cache. If None, default memory/sqlite cache is used.
    use_cache : bool, default=True
        Whether to cache query results locally to reduce network calls.
    mapping_dict : dict, optional
        Pre-existing offline dictionary {from_id: to_id} to bypass external queries.
    **query_kwargs
        Additional arguments passed to mygene.MyGeneInfo.querymany.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ['query', 'symbol', 'name', ...] indexed by query ID.
    """
    genes = list(dict.fromkeys(gene_list))  # Deduplicate preserving order

    # 1. Offline custom dictionary override
    if mapping_dict is not None:
        logger.info(f"Using provided offline mapping dictionary for {len(genes)} genes.")
        mapped_rows = []
        for g in genes:
            mapped_rows.append(
                {
                    "query": g,
                    to_id: mapping_dict.get(g, g),
                }
            )
        res_df = pd.DataFrame(mapped_rows).set_index("query")
        return res_df

    if not HAS_MYGENE:
        raise ImportError(
            "mygene package is required for online gene mapping. Run pip install mygene."
        )

    logger.info(
        f"Querying MyGene.info for {len(genes)} genes ({from_id} -> {to_id}, species={species})..."
    )

    mg = mygene.MyGeneInfo()

    if use_cache:
        c_dir = Path(cache_dir) if cache_dir is not None else Path(".cache_mygene")
        c_dir.mkdir(parents=True, exist_ok=True)
        cache_file = str(c_dir / "mygene_cache.sqlite")
        try:
            mg.set_caching(cache_db=cache_file)
        except Exception as e:
            logger.warning(f"Could not enable sqlite caching on MyGeneInfo: {e}")

    # Fields to query
    fields_to_request = list(set(["symbol", "name", "entrezgene", to_id]))

    try:
        raw_res = mg.querymany(
            genes,
            scopes=from_id,
            fields=",".join(fields_to_request),
            species=species,
            as_dataframe=True,
            returnall=False,
            **query_kwargs,
        )

        if isinstance(raw_res, pd.DataFrame):
            res_df = raw_res.copy()
            if not res_df.index.name:
                res_df.index.name = "query"
            # Deduplicate multiple hits by keeping the first hit
            if not res_df.index.is_unique:
                res_df = res_df[~res_df.index.duplicated(keep="first")]
        else:
            res_df = pd.DataFrame(raw_res)
            if "query" in res_df.columns:
                res_df = res_df.drop_duplicates(subset=["query"]).set_index("query")

        logger.info(f"Retrieved annotations for {len(res_df)} queried genes.")
        return res_df

    except Exception as exc:
        logger.warning(
            f"MyGene.info online query failed ({exc}). Falling back to identity mapping."
        )
        return pd.DataFrame({"query": genes, to_id: genes}).set_index("query")


def annotate_de_results(
    de_data: Union[DEResult, pd.DataFrame],
    id_col: Optional[str] = None,
    from_id: str = "ensembl.gene",
    to_id: str = "symbol",
    species: str = "human",
    mapping_dict: Optional[Dict[str, str]] = None,
    cache_dir: Optional[Union[str, Path]] = None,
    **kwargs,
) -> pd.DataFrame:
    """
    Annotate differential expression results with gene symbols and gene metadata.

    Parameters
    ----------
    de_data : DEResult or pd.DataFrame
        Differential expression analysis result table.
    id_col : str, optional
        Column containing gene IDs. If None, DataFrame index is used.
    from_id : str, default='ensembl.gene'
        Source identifier scope.
    to_id : str, default='symbol'
        Target gene symbol column name.
    species : str, default='human'
        Target organism.
    mapping_dict : dict, optional
        Pre-built dictionary {gene_id: symbol}.
    cache_dir : str or Path, optional
        Local cache folder.
    **kwargs
        Additional arguments passed to map_gene_ids.

    Returns
    -------
    pd.DataFrame
        Annotated DataFrame with mapped symbol inserted as first data column.
    """
    if isinstance(de_data, DEResult):
        df = de_data.results_df.copy()
    else:
        df = de_data.copy()

    if id_col is not None and id_col in df.columns:
        gene_ids = df[id_col].astype(str).tolist()
    else:
        gene_ids = df.index.astype(str).tolist()

    mapping_df = map_gene_ids(
        gene_list=gene_ids,
        from_id=from_id,
        to_id=to_id,
        species=species,
        mapping_dict=mapping_dict,
        cache_dir=cache_dir,
        **kwargs,
    )

    # Merge annotations
    out_df = df.copy()
    if to_id in mapping_df.columns:
        symbol_map = mapping_df[to_id].dropna().to_dict()
        if id_col is not None and id_col in out_df.columns:
            out_df[to_id] = [symbol_map.get(str(g), str(g)) for g in out_df[id_col]]
        else:
            out_df[to_id] = [symbol_map.get(str(g), str(g)) for g in out_df.index]

        # Put symbol as first column for readability
        cols = [to_id] + [c for c in out_df.columns if c != to_id]
        out_df = out_df[cols]

    if "name" in mapping_df.columns and "name" not in out_df.columns:
        name_map = mapping_df["name"].dropna().to_dict()
        if id_col is not None and id_col in out_df.columns:
            out_df["gene_name"] = [name_map.get(str(g), "") for g in out_df[id_col]]
        else:
            out_df["gene_name"] = [name_map.get(str(g), "") for g in out_df.index]

    logger.info(f"Annotated DE table ({len(out_df)} rows) with gene symbols.")
    return out_df
