"""
Network export and graph conversion utilities for BTW (FR-8).
Converts co-expression modules into NetworkX graphs and exports to Cytoscape SIF / edge lists.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import networkx as nx
import numpy as np
import pandas as pd

from btw import logger
from btw.network.wgcna_helper import WGCNAClusterResult


def module_to_networkx(
    network_data: Union[WGCNAClusterResult, pd.DataFrame],
    module: Optional[str] = None,
    threshold: float = 0.1,
    top_n_edges: Optional[int] = None,
) -> nx.Graph:
    """
    Convert a co-expression module or TOM matrix into a NetworkX Graph.

    Parameters
    ----------
    network_data : WGCNAClusterResult or pd.DataFrame
        WGCNA result object or TOM similarity matrix.
    module : str, optional
        Target module name (e.g. 'turquoise', 'blue').
        If None and WGCNAClusterResult is given, uses all genes across all non-grey modules.
    threshold : float, default=0.1
        Minimum edge weight (TOM) to include in graph.
    top_n_edges : int, optional
        If set, retains only the top N highest weighted edges.

    Returns
    -------
    nx.Graph
        NetworkX undirected graph with edge weights and node attributes (module, degree).
    """
    if isinstance(network_data, WGCNAClusterResult):
        if module is not None:
            genes = network_data.get_module_genes(module)
            if not genes:
                raise ValueError(
                    f"No genes found in module '{module}'. Available: {network_data.modules}"
                )
            tom_sub = network_data.tom_matrix.loc[genes, genes]
            module_map = {g: module for g in genes}
        else:
            genes = [g for g, m in network_data.module_labels.items() if m != "grey"]
            tom_sub = network_data.tom_matrix.loc[genes, genes]
            module_map = network_data.module_labels.to_dict()
    else:
        tom_sub = network_data.copy()
        module_map = {g: "module_1" for g in tom_sub.index}

    gene_list = list(tom_sub.index)
    tom_vals = tom_sub.values

    # Extract upper triangle edges
    i_upper, j_upper = np.triu_indices(len(gene_list), k=1)
    weights = tom_vals[i_upper, j_upper]

    # Filter by threshold
    valid_mask = weights >= threshold
    sources = [gene_list[i] for i in i_upper[valid_mask]]
    targets = [gene_list[j] for j in j_upper[valid_mask]]
    edge_weights = weights[valid_mask]

    edge_df = pd.DataFrame(
        {
            "source": sources,
            "target": targets,
            "weight": edge_weights,
        }
    )

    if top_n_edges is not None and len(edge_df) > top_n_edges:
        edge_df = (
            edge_df.sort_values(by="weight", ascending=False)
            .head(top_n_edges)
            .reset_index(drop=True)
        )

    # Build NetworkX graph
    G = nx.Graph()

    # Add nodes with attributes
    for g in gene_list:
        if g in edge_df["source"].values or g in edge_df["target"].values:
            G.add_node(g, module=module_map.get(g, "unknown"))

    # Add edges
    for _, row in edge_df.iterrows():
        G.add_edge(row["source"], row["target"], weight=float(row["weight"]))

    # Annotate degree
    for node in G.nodes():
        G.nodes[node]["degree"] = G.degree[node]

    logger.info(
        f"Constructed NetworkX graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges (module={module})."
    )
    return G


def export_cytoscape_sif(
    graph_or_edges: Union[nx.Graph, pd.DataFrame],
    output_path: Union[str, Path],
    interaction_type: str = "coexpression",
) -> Path:
    """
    Export network to Cytoscape Simple Interaction Format (.sif).

    Format per line:
    NodeA <tab> interaction_type <tab> NodeB

    Parameters
    ----------
    graph_or_edges : nx.Graph or DataFrame
        Input network graph or DataFrame with ['source', 'target'] columns.
    output_path : str or Path
        Destination .sif file path.
    interaction_type : str, default='coexpression'
        Label for edge interaction type.

    Returns
    -------
    Path
        Path to exported .sif file.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    lines = []
    if isinstance(graph_or_edges, nx.Graph):
        for u, v in graph_or_edges.edges():
            lines.append(f"{u}\t{interaction_type}\t{v}")
    elif isinstance(graph_or_edges, pd.DataFrame):
        if not {"source", "target"}.issubset(graph_or_edges.columns):
            raise ValueError("DataFrame must contain ['source', 'target'] columns.")
        for _, row in graph_or_edges.iterrows():
            lines.append(f"{row['source']}\t{interaction_type}\t{row['target']}")
    else:
        raise TypeError("Input must be a networkx.Graph or pandas DataFrame.")

    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + ("\n" if lines else ""))

    logger.info(f"Exported {len(lines)} interactions to Cytoscape SIF: {path}")
    return path


def export_edge_list(
    graph_or_edges: Union[nx.Graph, pd.DataFrame],
    output_path: Union[str, Path],
    sep: str = "\t",
) -> Path:
    """
    Export network edges with weights to a tabular file (TSV or CSV).

    Parameters
    ----------
    graph_or_edges : nx.Graph or DataFrame
        Input network graph or DataFrame.
    output_path : str or Path
        Destination file path.
    sep : str, default='\t'
        Delimiter.

    Returns
    -------
    Path
        Path to exported edge list file.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if isinstance(graph_or_edges, nx.Graph):
        edge_data = []
        for u, v, data in graph_or_edges.edges(data=True):
            edge_data.append(
                {
                    "source": u,
                    "target": v,
                    "weight": data.get("weight", 1.0),
                    "source_module": graph_or_edges.nodes[u].get("module", ""),
                    "target_module": graph_or_edges.nodes[v].get("module", ""),
                }
            )
        df = pd.DataFrame(edge_data)
    elif isinstance(graph_or_edges, pd.DataFrame):
        df = graph_or_edges.copy()
    else:
        raise TypeError("Input must be a networkx.Graph or pandas DataFrame.")

    df.to_csv(path, sep=sep, index=False)
    logger.info(f"Exported {len(df)} edges to {path}")
    return path
