"""
Unit tests for Functional & Pathway Enrichment module (FR-5).
Tests ORA, GSEA, decoupler activity inference, standardized schema, and unified plots.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import pytest

from btw.enrichment import (
    EnrichmentResult,
    extract_significant_genes,
    get_memory_cache,
    plot_enrichment_barplot,
    plot_enrichment_dotplot,
    prepare_decoupler_input,
    prepare_ranked_gene_list,
    run_custom_ora,
    run_decoupler_activity,
    run_enrichr,
    run_prerank,
    standardize_enrichment_table,
)


@pytest.fixture
def custom_gene_sets():
    """Custom biological pathways for deterministic testing."""
    return {
        "Cell_Proliferation": ["GENE_001", "GENE_002", "GENE_003", "GENE_004", "GENE_005"],
        "DNA_Repair": ["GENE_006", "GENE_007", "GENE_008", "GENE_009", "GENE_010"],
        "Metabolism": ["GENE_080", "GENE_081", "GENE_082", "GENE_083", "GENE_084"],
    }


def test_schema_standardization():
    """Test standardizing third-party tables into uniform schema."""
    raw_df = pd.DataFrame(
        {
            "Term": ["Pathway_A", "Pathway_B"],
            "NES": [2.4, -1.8],
            "NOM p-val": [0.001, 0.04],
            "FDR q-val": [0.005, 0.08],
            "Lead_genes": ["g1;g2", "g3"],
            "Overlap": [2, 1],
        }
    )

    std_df = standardize_enrichment_table(
        raw_df,
        source_method="test_source",
        term_col="Term",
        score_col="NES",
        pvalue_col="NOM p-val",
        padj_col="FDR q-val",
        genes_col="Lead_genes",
    )

    for col in ["term", "score", "pvalue", "padj", "gene_count", "genes", "source_method"]:
        assert col in std_df.columns

    res = EnrichmentResult(
        results_df=std_df,
        source_method="test_source",
        gene_set_database="test_db",
        alpha=0.05,
    )
    sig_df = res.get_significant()
    assert len(sig_df) == 1
    assert "Pathway_A" in list(sig_df["term"])
    assert "Enrichment Summary" in res.summary()


def test_memory_cache(temp_dir):
    """Test joblib.Memory caching wrapper."""
    cache = get_memory_cache(cache_dir=temp_dir, enabled=True)

    @cache.cache
    def heavy_function(x):
        return x * 2

    val1 = heavy_function(10)
    val2 = heavy_function(10)
    assert val1 == val2 == 20


def test_ora_extraction_and_custom_fisher(de_result_fixture, custom_gene_sets):
    """Test ORA gene extraction and Fisher exact test."""
    # 1. Extraction
    up_genes = extract_significant_genes(de_result_fixture, direction="up")
    assert len(up_genes) == 15
    assert "GENE_001" in up_genes

    # 2. Custom ORA
    bg_genes = [f"GENE_{i:03d}" for i in range(1, 101)]
    res = run_custom_ora(up_genes, custom_gene_sets, background=bg_genes, alpha=0.05)

    assert isinstance(res, EnrichmentResult)
    assert len(res.results_df) == 3
    # Cell Proliferation has all 5 genes in up_genes, should be highly enriched
    top_term = res.results_df.iloc[0]["term"]
    assert top_term in ["Cell_Proliferation", "DNA_Repair"]
    assert res.results_df.iloc[0]["padj"] < 0.05

    # 3. Test run_enrichr dispatcher with dict
    res_enr = run_enrichr(de_result_fixture, gene_sets=custom_gene_sets, direction="up")
    assert isinstance(res_enr, EnrichmentResult)
    assert len(res_enr.results_df) == 3


def test_gsea_ranked_list_and_prerank(de_result_fixture, custom_gene_sets):
    """Test GSEA ranked gene list generation and prerank execution."""
    # Test rank_by='stat'
    ranked_stat = prepare_ranked_gene_list(de_result_fixture, rank_by="stat")
    assert isinstance(ranked_stat, pd.Series)
    assert len(ranked_stat) == 100
    # Top ranked gene should have high positive stat
    assert ranked_stat.iloc[0] > ranked_stat.iloc[-1]

    # Test rank_by='log2FoldChange'
    ranked_lfc = prepare_ranked_gene_list(de_result_fixture, rank_by="log2FoldChange")
    assert len(ranked_lfc) == 100

    # Test prerank with custom dictionary
    gsea_res = run_prerank(
        ranked_genes=ranked_stat,
        gene_sets=custom_gene_sets,
        min_size=2,
        permutation_num=20,
        seed=42,
    )
    assert isinstance(gsea_res, EnrichmentResult)
    assert gsea_res.source_method == "gseapy_prerank"
    assert len(gsea_res.results_df) > 0


def test_decoupler_input_and_activity(de_result_fixture):
    """Test formatting DE stats for decoupler and running activity inference."""
    # 1. Format input matrix
    mat = prepare_decoupler_input(de_result_fixture, metric="stat")
    assert mat.shape == (1, 100)

    # 2. Build mock network
    net = pd.DataFrame(
        {
            "source": ["TF_1", "TF_1", "TF_1", "TF_2", "TF_2", "TF_2"],
            "target": ["GENE_001", "GENE_002", "GENE_003", "GENE_016", "GENE_017", "GENE_018"],
            "weight": [1.0, 1.2, 0.9, 1.0, 1.1, 0.8],
        }
    )

    # 3. Run decoupler ULM
    act_res = run_decoupler_activity(mat, net, method="ulm", min_n=2)
    assert isinstance(act_res, EnrichmentResult)
    assert act_res.source_method == "decoupler_ulm"
    assert len(act_res.results_df) == 2
    # TF_1 targets up-regulated genes (GENE_001-003) -> positive activity score
    tf1_row = act_res.results_df.loc[act_res.results_df["term"] == "TF_1"].iloc[0]
    assert tf1_row["score"] > 0
    # TF_2 targets down-regulated genes (GENE_016-018) -> negative activity score
    tf2_row = act_res.results_df.loc[act_res.results_df["term"] == "TF_2"].iloc[0]
    assert tf2_row["score"] < 0


def test_enrichment_visualizations(de_result_fixture, custom_gene_sets):
    """Test unified Dot Plot and Bar Plot across enrichment results."""
    up_genes = extract_significant_genes(de_result_fixture, direction="up")
    res = run_custom_ora(up_genes, custom_gene_sets)

    # 1. Dot Plot (static)
    fig_dot, ax_dot = plot_enrichment_dotplot(res, top_n=5, interactive=False)
    assert isinstance(fig_dot, plt.Figure)
    assert isinstance(ax_dot, plt.Axes)
    plt.close(fig_dot)

    # 2. Dot Plot (interactive)
    plotly_dot = plot_enrichment_dotplot(res, top_n=5, interactive=True)
    assert plotly_dot is not None

    # 3. Bar Plot (static)
    fig_bar, ax_bar = plot_enrichment_barplot(res, top_n=5, interactive=False)
    assert isinstance(fig_bar, plt.Figure)
    plt.close(fig_bar)

    # 4. Bar Plot (interactive)
    plotly_bar = plot_enrichment_barplot(res, top_n=5, interactive=True)
    assert plotly_bar is not None
