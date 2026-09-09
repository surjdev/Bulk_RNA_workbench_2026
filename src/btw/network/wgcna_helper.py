"""
Weighted Gene Co-expression Network Analysis (WGCNA) helper module for BTW (FR-8).
Provides standard WGCNA topological overlap matrix (TOM), module detection,
module eigengene extraction, and integration with PyWGCNA.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from btw import logger
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr
from sklearn.decomposition import PCA

try:
    import PyWGCNA

    HAS_PYWGCNA = True
except ImportError:
    HAS_PYWGCNA = False


# Standard WGCNA module color palette
WGCNA_COLORS = [
    "turquoise",
    "blue",
    "brown",
    "yellow",
    "green",
    "red",
    "black",
    "pink",
    "magenta",
    "purple",
    "greenyellow",
    "tan",
    "salmon",
    "cyan",
    "midnightblue",
    "lightcyan",
    "darkred",
    "darkgreen",
    "lightgreen",
    "orange",
]


@dataclass
class WGCNAClusterResult:
    """
    Standardized container holding WGCNA co-expression network analysis results.
    """
    module_labels: pd.Series            # Gene ID -> Module Color string
    module_eigengenes: pd.DataFrame     # Samples x Modules (ME<color>)
    tom_matrix: pd.DataFrame            # Genes x Genes topological overlap
    adjacency_matrix: pd.DataFrame      # Genes x Genes soft-thresholded adjacency
    power: int = 6
    network_type: str = "unsigned"
    module_trait_cor: Optional[pd.DataFrame] = None
    module_trait_pval: Optional[pd.DataFrame] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def modules(self) -> List[str]:
        """Unique module names/colors excluding grey unassigned genes."""
        mods = [m for m in self.module_labels.unique() if m != "grey"]
        return sorted(mods)

    @property
    def all_modules(self) -> List[str]:
        """All unique module names including grey."""
        return sorted(self.module_labels.unique().tolist())

    def get_module_genes(self, module: str) -> List[str]:
        """Retrieve list of gene IDs assigned to a specific module."""
        return self.module_labels[self.module_labels == module].index.tolist()

    def module_sizes(self) -> pd.Series:
        """Count of genes in each co-expression module."""
        return self.module_labels.value_counts()

    def summary(self) -> str:
        """Text summary of detected co-expression network modules."""
        sizes = self.module_sizes()
        n_total = len(self.module_labels)
        n_mods = len([m for m in sizes.index if m != "grey"])
        n_grey = sizes.get("grey", 0)

        lines = [
            f"--- WGCNA Co-expression Summary (Power={self.power}, Type={self.network_type}) ---",
            f"Total genes analyzed: {n_total}",
            f"Identified co-expression modules: {n_mods} (unassigned grey genes: {n_grey})",
            "Module Sizes:",
        ]
        for mod, count in sizes.items():
            pct = (count / n_total) * 100
            lines.append(f"  - Module {mod:12s}: {count:4d} genes ({pct:.1f}%)")
        return "\n".join(lines)


def compute_adjacency(
    data: pd.DataFrame,
    power: int = 6,
    network_type: str = "unsigned",
    cor_method: str = "pearson",
) -> pd.DataFrame:
    """
    Calculate soft-thresholded adjacency matrix from gene expression data.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix of shape (n_samples, n_genes) or (n_genes, n_samples).
    power : int, default=6
        Soft-thresholding exponent beta.
    network_type : {'unsigned', 'signed'}, default='unsigned'
        Network topology type.
    cor_method : {'pearson', 'spearman'}, default='pearson'
        Correlation method.

    Returns
    -------
    pd.DataFrame
        Symmetric Adjacency matrix of shape (n_genes, n_genes).
    """
    # Orient to (n_samples, n_genes)
    if data.shape[0] > data.shape[1]:
        # Typically n_genes > n_samples, so if rows > cols, transpose
        df = data.T
    else:
        df = data.copy()

    # Calculate correlation
    cor_df = df.corr(method=cor_method).fillna(0.0)

    if network_type == "unsigned":
        adj = np.abs(cor_df.values) ** power
    elif network_type == "signed":
        adj = (0.5 * (1.0 + cor_df.values)) ** power
    else:
        raise ValueError(f"Invalid network_type '{network_type}'. Choose 'unsigned' or 'signed'.")

    np.fill_diagonal(adj, 1.0)
    return pd.DataFrame(adj, index=cor_df.index, columns=cor_df.columns)


def compute_tom(adjacency: Union[pd.DataFrame, np.ndarray]) -> np.ndarray:
    """
    Calculate Topological Overlap Matrix (TOM) from soft-thresholded adjacency.

    Parameters
    ----------
    adjacency : pd.DataFrame or np.ndarray
        Symmetric adjacency matrix with values in [0, 1].

    Returns
    -------
    np.ndarray
        TOM similarity matrix of shape (n_genes, n_genes).
    """
    vals = adjacency.values if isinstance(adjacency, pd.DataFrame) else adjacency
    A = np.array(vals, copy=True, dtype=float)
    np.fill_diagonal(A, 1.0)

    # Numerator: L_ij + A_ij = (A @ A) - A
    L_plus_A = np.dot(A, A) - A

    # Denominator: min(k_i, k_j) + 1 - A_ij
    k = A.sum(axis=1) - 1.0
    min_k = np.minimum.outer(k, k)
    denom = min_k + 1.0 - A

    # Avoid zero division
    denom = np.maximum(denom, 1e-12)

    tom = L_plus_A / denom
    np.fill_diagonal(tom, 1.0)
    return np.clip(tom, 0.0, 1.0)


def detect_coexpression_modules(
    data: pd.DataFrame,
    power: int = 6,
    network_type: str = "unsigned",
    min_module_size: int = 10,
    max_modules: int = 15,
    sample_metadata: Optional[pd.DataFrame] = None,
    traits: Optional[Union[List[str], pd.DataFrame]] = None,
) -> WGCNAClusterResult:
    """
    Perform complete WGCNA co-expression network analysis:
    Adjacency -> TOM -> Hierarchical Clustering -> Module Eigengenes -> Trait Correlation.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix. If rows > columns, treated as (genes x samples) and transposed.
    power : int, default=6
        Soft thresholding power.
    network_type : {'unsigned', 'signed'}, default='unsigned'
        Co-expression topology type.
    min_module_size : int, default=10
        Minimum number of genes to form a distinct module.
    max_modules : int, default=15
        Upper bound on the number of modules formed.
    sample_metadata : pd.DataFrame, optional
        Sample metadata for module-trait correlation analysis.
    traits : list of str or DataFrame, optional
        Traits to correlate against module eigengenes.

    Returns
    -------
    WGCNAClusterResult
        Structured container with module assignments, TOM, and eigengenes.
    """
    # 1. Orientation: ensure samples as rows, genes as columns
    if data.shape[0] > data.shape[1]:
        expr_df = data.T.copy()
        logger.info(f"Transposed expression matrix to ({expr_df.shape[0]} samples x {expr_df.shape[1]} genes).")
    else:
        expr_df = data.copy()

    gene_names = list(expr_df.columns)
    sample_names = list(expr_df.index)
    n_genes = len(gene_names)
    n_samples = len(sample_names)

    logger.info(f"Running WGCNA on {n_genes} genes across {n_samples} samples (power={power}, type={network_type})...")

    # 2. Adjacency and TOM
    adj_df = compute_adjacency(expr_df, power=power, network_type=network_type)
    tom = compute_tom(adj_df)
    tom_df = pd.DataFrame(tom, index=gene_names, columns=gene_names)

    # 3. Dissimilarity and Hierarchical Clustering
    diss_tom = np.clip(1.0 - tom, 0.0, 1.0)
    np.fill_diagonal(diss_tom, 0.0)

    # Convert condensed distance matrix
    condensed_dist = squareform(diss_tom, checks=False)
    z_linkage = linkage(condensed_dist, method="average")

    # 4. Dynamic/Hierarchical cut to find modules
    # Determine number of clusters based on gene count and min_module_size
    target_k = min(max_modules, max(2, n_genes // max(min_module_size, 5)))
    clusters = fcluster(z_linkage, t=target_k, criterion="maxclust")

    # 5. Map cluster IDs to WGCNA color names
    unique_clusters, counts = np.unique(clusters, return_counts=True)
    # Sort clusters by size descending
    sorted_clusters = unique_clusters[np.argsort(-counts)]

    cluster_to_color: Dict[int, str] = {}
    color_idx = 0
    for cl_id in sorted_clusters:
        cl_size = np.sum(clusters == cl_id)
        if cl_size >= min_module_size and color_idx < len(WGCNA_COLORS):
            cluster_to_color[cl_id] = WGCNA_COLORS[color_idx]
            color_idx += 1
        else:
            cluster_to_color[cl_id] = "grey"

    module_labels = pd.Series(
        [cluster_to_color[cl] for cl in clusters],
        index=gene_names,
        name="module",
    )

    # 6. Compute Module Eigengenes (ME) via PCA
    detected_modules = [m for m in WGCNA_COLORS if m in module_labels.values]
    me_dict: Dict[str, np.ndarray] = {}

    for mod in detected_modules:
        mod_genes = module_labels[module_labels == mod].index
        if len(mod_genes) < 2:
            continue
        sub_expr = expr_df[mod_genes].values
        # Standardize genes (Z-score)
        std_vals = np.std(sub_expr, axis=0, ddof=1)
        std_vals[std_vals == 0] = 1.0
        sub_z = (sub_expr - np.mean(sub_expr, axis=0)) / std_vals

        pca = PCA(n_components=1)
        eigengene = pca.fit_transform(sub_z).flatten()

        # Ensure positive average correlation with module genes
        corrs = [np.corrcoef(eigengene, sub_expr[:, j])[0, 1] for j in range(sub_expr.shape[1])]
        if np.nanmean(corrs) < 0:
            eigengene = -eigengene

        me_dict[f"ME{mod}"] = eigengene

    me_df = pd.DataFrame(me_dict, index=sample_names)

    # 7. Module-Trait Correlation if metadata / traits provided
    trait_cor_df = None
    trait_pval_df = None
    if sample_metadata is not None and not me_df.empty:
        common_samples = [s for s in sample_names if s in sample_metadata.index]
        if common_samples:
            sub_me = me_df.loc[common_samples]
            sub_meta = sample_metadata.loc[common_samples]

            if traits is not None:
                trait_cols = [c for c in (traits if isinstance(traits, list) else traits.columns) if c in sub_meta.columns]
            else:
                # Include numeric columns or convert binary categories to 0/1
                trait_cols = []
                for col in sub_meta.columns:
                    if pd.api.types.is_numeric_dtype(sub_meta[col]):
                        trait_cols.append(col)
                    elif len(sub_meta[col].unique()) == 2:
                        trait_cols.append(col)

            if trait_cols:
                cor_mat = np.zeros((len(sub_me.columns), len(trait_cols)))
                pval_mat = np.zeros((len(sub_me.columns), len(trait_cols)))

                for i, me_name in enumerate(sub_me.columns):
                    for j, t_name in enumerate(trait_cols):
                        t_vals = pd.to_numeric(sub_meta[t_name], errors="coerce").fillna(0).values
                        r, p = pearsonr(sub_me[me_name].values, t_vals)
                        cor_mat[i, j] = r
                        pval_mat[i, j] = p

                trait_cor_df = pd.DataFrame(cor_mat, index=sub_me.columns, columns=trait_cols)
                trait_pval_df = pd.DataFrame(pval_mat, index=sub_me.columns, columns=trait_cols)
                logger.info(f"Computed module-trait correlations across {len(trait_cols)} traits.")

    result = WGCNAClusterResult(
        module_labels=module_labels,
        module_eigengenes=me_df,
        tom_matrix=tom_df,
        adjacency_matrix=adj_df,
        power=power,
        network_type=network_type,
        module_trait_cor=trait_cor_df,
        module_trait_pval=trait_pval_df,
    )
    logger.info(f"Detected {len(result.modules)} co-expression modules: {', '.join(result.modules)}")
    return result


def run_pywgcna(
    data: pd.DataFrame,
    metadata: pd.DataFrame,
    species: str = "human",
    min_module_size: int = 15,
    **kwargs,
) -> Any:
    """
    Execute PyWGCNA full pipeline wrapper.

    Parameters
    ----------
    data : pd.DataFrame
        Expression matrix.
    metadata : pd.DataFrame
        Sample annotations.
    species : str, default='human'
        Target species ('human', 'mouse').
    min_module_size : int, default=15
        Minimum cluster size.
    **kwargs
        Additional arguments passed to PyWGCNA.WGCNA.

    Returns
    -------
    PyWGCNA.WGCNA
        Fitted PyWGCNA analysis object.
    """
    if not HAS_PYWGCNA:
        raise ImportError("PyWGCNA package is required. Run pip install pywgcna.")

    logger.info(f"Initializing PyWGCNA pipeline for species='{species}'...")
    # Orient genes as columns, samples as rows
    if data.shape[0] > data.shape[1]:
        expr = data.T.copy()
    else:
        expr = data.copy()

    pyw = PyWGCNA.WGCNA(
        name="BTW_WGCNA",
        species=species,
        minModuleSize=min_module_size,
        **kwargs,
    )
    pyw.geneExp = expr
    pyw.sampleInfo = metadata
    return pyw
