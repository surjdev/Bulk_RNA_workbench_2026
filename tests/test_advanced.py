"""
Comprehensive Unit Tests for Part 5: Advanced Downstream Modules.
Covers:
- FR-6: Gene Annotation & GTF parsing (src/btw/annotation)
- FR-7: ComBat Batch Effect Correction & PCA QC (src/btw/batch_correction)
- FR-8: WGCNA Co-expression Networks & Cytoscape SIF Export (src/btw/network)
"""

import io
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import pytest
from btw.annotation import (
    annotate_de_results,
    create_gene_map_from_gtf,
    map_gene_ids,
    parse_gtf_file,
)
from btw.batch_correction import (
    compare_pca_batch,
    evaluate_batch_effect,
    run_combat,
)
from btw.de_analysis.contrasts import DEResult
from btw.network import (
    WGCNAClusterResult,
    compute_adjacency,
    compute_tom,
    detect_coexpression_modules,
    export_cytoscape_sif,
    export_edge_list,
    module_to_networkx,
)


# =====================================================================
# FR-6: Annotation & GTF Parsing Tests
# =====================================================================

GTF_MOCK_DATA = """chr1\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tgene_id "ENSG000001"; gene_name "TP53"; gene_biotype "protein_coding";
chr1\tHAVANA\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "ENSG000001"; transcript_id "ENST000001"; gene_name "TP53";
chr1\tHAVANA\tgene\t3000\t4000\t.\t-\t.\tgene_id "ENSG000002"; gene_name "BRCA1"; gene_biotype "protein_coding";
chr2\tHAVANA\tgene\t5000\t6000\t.\t+\t.\tgene_id "ENSG000003"; gene_name "MYC"; gene_biotype "protein_coding";
"""


def test_gtf_parsing_and_mapping():
    """Test parsing GTF format text into DataFrame and extracting gene maps."""
    buf = io.StringIO(GTF_MOCK_DATA)
    df = parse_gtf_file(buf, features=["gene"])

    assert isinstance(df, pd.DataFrame)
    assert len(df) == 3
    assert "gene_id" in df.columns
    assert "gene_name" in df.columns
    assert "TP53" in df["gene_name"].values

    # Test create_gene_map_from_gtf
    buf2 = io.StringIO(GTF_MOCK_DATA)
    gene_map = create_gene_map_from_gtf(buf2, from_attr="gene_id", to_attr="gene_name")
    assert isinstance(gene_map, pd.Series)
    assert gene_map["ENSG000001"] == "TP53"
    assert gene_map["ENSG000002"] == "BRCA1"
    assert gene_map["ENSG000003"] == "MYC"


def test_gene_mapping_offline(de_result_fixture):
    """Test mapping gene IDs to symbols with offline mapping dictionary."""
    mapping_dict = {
        "GENE_001": "CDK1",
        "GENE_002": "CCNB1",
        "GENE_003": "AURKA",
        "GENE_004": "PLK1",
    }

    genes = ["GENE_001", "GENE_002", "GENE_003", "GENE_004", "GENE_099"]
    res_df = map_gene_ids(genes, mapping_dict=mapping_dict)

    assert isinstance(res_df, pd.DataFrame)
    assert res_df.loc["GENE_001", "symbol"] == "CDK1"
    assert res_df.loc["GENE_002", "symbol"] == "CCNB1"
    # Unmapped gene falls back to query ID
    assert res_df.loc["GENE_099", "symbol"] == "GENE_099"

    # Test annotate_de_results
    annotated = annotate_de_results(de_result_fixture, mapping_dict=mapping_dict)
    assert isinstance(annotated, pd.DataFrame)
    assert "symbol" in annotated.columns
    assert annotated.columns[0] == "symbol"  # Placed as first column
    assert annotated.loc["GENE_001", "symbol"] == "CDK1"


# =====================================================================
# FR-7: ComBat Batch Effect Correction Tests
# =====================================================================


@pytest.fixture
def batch_dataset_fixture():
    """Create synthetic count dataset with clear technical batch effect."""
    np.random.seed(42)
    n_genes, n_samples = 60, 6
    genes = [f"G_{i:02d}" for i in range(n_genes)]
    samples = [f"S_{i:02d}" for i in range(n_samples)]

    counts = np.random.negative_binomial(n=25, p=0.08, size=(n_genes, n_samples))
    # Inject large technical offset into batch b2 (samples 3, 4, 5)
    counts[:, 3:] += 150

    counts_df = pd.DataFrame(counts, index=genes, columns=samples)
    metadata = pd.DataFrame(
        {
            "batch": ["b1", "b1", "b1", "b2", "b2", "b2"],
            "condition": ["ctrl", "treat", "ctrl", "ctrl", "treat", "treat"],
        },
        index=samples,
    )
    return counts_df, metadata


def test_run_combat_seq_counts(batch_dataset_fixture):
    """Test ComBat-Seq on integer count data."""
    counts_df, metadata = batch_dataset_fixture

    corrected = run_combat(
        data=counts_df,
        batch="batch",
        metadata=metadata,
        biological_factor="condition",
        is_count=True,
    )

    assert isinstance(corrected, pd.DataFrame)
    assert corrected.shape == counts_df.shape
    assert list(corrected.index) == list(counts_df.index)
    assert list(corrected.columns) == list(counts_df.columns)
    # Check that batch differences in mean expression are significantly reduced
    mean_diff_before = abs(counts_df.iloc[:, 3:].values.mean() - counts_df.iloc[:, :3].values.mean())
    mean_diff_after = abs(corrected.iloc[:, 3:].values.mean() - corrected.iloc[:, :3].values.mean())
    assert mean_diff_after < mean_diff_before


def test_run_combat_norm(batch_dataset_fixture):
    """Test parametric ComBat on continuous normalized expression data."""
    counts_df, metadata = batch_dataset_fixture
    log_norm = np.log2(counts_df + 1.0)

    corrected = run_combat(
        data=log_norm,
        batch=metadata["batch"],
        is_count=False,
    )

    assert isinstance(corrected, pd.DataFrame)
    assert corrected.shape == counts_df.shape


def test_batch_diagnostic_pca(batch_dataset_fixture):
    """Test PCA comparison plots and batch silhouette metrics."""
    counts_df, metadata = batch_dataset_fixture
    log_norm = np.log2(counts_df + 1.0)

    corrected = run_combat(
        data=log_norm,
        batch="batch",
        metadata=metadata,
        is_count=False,
    )

    # Evaluate metrics
    m_before = evaluate_batch_effect(log_norm, metadata, batch_col="batch", condition_col="condition")
    m_after = evaluate_batch_effect(corrected, metadata, batch_col="batch", condition_col="condition")

    assert "batch_silhouette" in m_before
    assert "batch_silhouette" in m_after
    # Batch silhouette should decrease after correction
    assert m_after["batch_silhouette"] <= m_before["batch_silhouette"]

    # Side-by-side plot
    fig, (ax1, ax2), metrics = compare_pca_batch(
        data_before=log_norm,
        data_after=corrected,
        metadata=metadata,
        batch_col="batch",
        condition_col="condition",
    )
    assert isinstance(fig, plt.Figure)
    assert isinstance(ax1, plt.Axes)
    assert isinstance(ax2, plt.Axes)
    assert "batch_silhouette_reduction" in metrics
    plt.close(fig)


# =====================================================================
# FR-8: WGCNA Co-expression Network Tests
# =====================================================================


@pytest.fixture
def coexpression_dataset_fixture():
    """Create expression data with distinct co-expressed gene modules."""
    np.random.seed(42)
    n_samples = 12
    n_genes = 40
    samples = [f"Sample_{i:02d}" for i in range(n_samples)]
    genes = [f"Gene_{i:03d}" for i in range(n_genes)]

    # Latent module signals
    signal_m1 = np.sin(np.linspace(0, 3 * np.pi, n_samples))
    signal_m2 = np.cos(np.linspace(0, 3 * np.pi, n_samples))

    expr = np.random.normal(0, 0.2, size=(n_samples, n_genes))
    # Module 1 (first 15 genes)
    for j in range(15):
        expr[:, j] += signal_m1 + np.random.normal(0, 0.1, n_samples)
    # Module 2 (next 15 genes)
    for j in range(15, 30):
        expr[:, j] += signal_m2 + np.random.normal(0, 0.1, n_samples)

    expr_df = pd.DataFrame(expr, index=samples, columns=genes)
    metadata = pd.DataFrame(
        {
            "phenotype_score": signal_m1 * 2.0,
        },
        index=samples,
    )
    return expr_df, metadata


def test_wgcna_adjacency_and_tom(coexpression_dataset_fixture):
    """Test soft-thresholded adjacency and Topological Overlap Matrix."""
    expr_df, _ = coexpression_dataset_fixture

    adj = compute_adjacency(expr_df, power=6, network_type="unsigned")
    assert adj.shape == (40, 40)
    assert np.all(adj.values >= 0.0)
    assert np.all(adj.values <= 1.0)
    assert np.allclose(np.diag(adj.values), 1.0)

    tom = compute_tom(adj)
    assert tom.shape == (40, 40)
    assert np.all(tom >= 0.0)
    assert np.all(tom <= 1.0)
    assert np.allclose(np.diag(tom), 1.0)


def test_wgcna_module_detection_and_eigengenes(coexpression_dataset_fixture):
    """Test dynamic module detection, module sizes, and eigengene correlation."""
    expr_df, metadata = coexpression_dataset_fixture

    res = detect_coexpression_modules(
        data=expr_df,
        power=6,
        min_module_size=8,
        sample_metadata=metadata,
        traits=["phenotype_score"],
    )

    assert isinstance(res, WGCNAClusterResult)
    assert len(res.module_labels) == 40
    assert len(res.modules) >= 1
    assert "grey" in res.all_modules or len(res.modules) > 0

    # Module eigengenes
    assert isinstance(res.module_eigengenes, pd.DataFrame)
    assert len(res.module_eigengenes) == 12

    # Module sizes & summary
    sizes = res.module_sizes()
    assert sizes.sum() == 40
    summary_text = res.summary()
    assert "WGCNA Co-expression Summary" in summary_text

    # Module-trait correlation
    assert res.module_trait_cor is not None
    assert "phenotype_score" in res.module_trait_cor.columns


def test_networkx_and_cytoscape_export(coexpression_dataset_fixture, tmp_path):
    """Test graph conversion with NetworkX and export to Cytoscape SIF and edge list."""
    expr_df, _ = coexpression_dataset_fixture

    res = detect_coexpression_modules(expr_df, power=6, min_module_size=8)
    top_mod = res.modules[0]

    # Convert to NetworkX
    G = module_to_networkx(res, module=top_mod, threshold=0.05)
    assert isinstance(G, nx.Graph)
    assert G.number_of_nodes() > 0
    assert G.number_of_edges() > 0

    # Node attributes
    sample_node = list(G.nodes())[0]
    assert "module" in G.nodes[sample_node]
    assert "degree" in G.nodes[sample_node]

    # Export Cytoscape SIF
    sif_file = tmp_path / "network.sif"
    exported_sif = export_cytoscape_sif(G, sif_file, interaction_type="coexpressed")
    assert exported_sif.exists()
    sif_content = exported_sif.read_text(encoding="utf-8")
    assert "coexpressed" in sif_content

    # Export Edge list
    edge_file = tmp_path / "edges.tsv"
    exported_edges = export_edge_list(G, edge_file)
    assert exported_edges.exists()
    edge_df = pd.read_csv(exported_edges, sep="\t")
    assert {"source", "target", "weight"}.issubset(edge_df.columns)


def test_combat_transposed_and_auto_count(batch_dataset_fixture):
    """Test ComBat auto-detecting count data and handling transposed samples x genes orientation."""
    counts_df, metadata = batch_dataset_fixture
    # Transpose to (samples x genes)
    counts_transposed = counts_df.T

    corrected = run_combat(
        data=counts_transposed,
        batch="batch",
        metadata=metadata,
        is_count=None,  # Auto-detection
    )
    assert corrected.shape == counts_transposed.shape
    assert list(corrected.index) == list(counts_transposed.index)


def test_annotate_de_results_with_id_col(de_result_fixture):
    """Test annotating DE results when gene ID is in a column instead of index."""
    df = de_result_fixture.results_df.reset_index()
    # Now 'gene_id' is a column
    mapping_dict = {"GENE_001": "CDK1", "GENE_002": "CCNB1"}
    annotated = annotate_de_results(df, id_col="gene_id", mapping_dict=mapping_dict)
    assert "symbol" in annotated.columns
    assert annotated.loc[0, "symbol"] == "CDK1"


def test_export_dataframe_sif_and_edge_list(tmp_path):
    """Test exporting network directly from DataFrame with source/target columns."""
    df = pd.DataFrame(
        {
            "source": ["GeneA", "GeneB"],
            "target": ["GeneB", "GeneC"],
            "weight": [0.85, 0.72],
        }
    )
    sif_path = export_cytoscape_sif(df, tmp_path / "from_df.sif", interaction_type="pp")
    assert sif_path.exists()
    assert "GeneA\tpp\tGeneB" in sif_path.read_text(encoding="utf-8")

    edge_path = export_edge_list(df, tmp_path / "from_df.tsv")
    assert edge_path.exists()

