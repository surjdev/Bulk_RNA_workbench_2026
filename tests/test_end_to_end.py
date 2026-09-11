"""
End-to-End Workbench Integration Tests for BTW (FR-1 through FR-9 & Acceptance Criteria).
Verifies the complete pipeline from raw counts to publication bundle and confirms
that underlying original library objects remain fully accessible.
"""

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sklearn.decomposition import PCA

from btw import set_seed
from btw.batch_correction import run_combat
from btw.de_analysis import run_de
from btw.enrichment.ora import run_custom_ora
from btw.enrichment.schema import EnrichmentResult
from btw.io import export_table, load_dataset, validate_bulk_data
from btw.qc_normalize import compute_sample_qc, filter_low_expression_genes, normalize_deseq2
from btw.reporting import export_publication_bundle
from btw.viz import compute_pca, plot_heatmap, plot_pca, plot_volcano, save_figure


def test_complete_end_to_end_pipeline(synthetic_data, tmp_path):
    """
    Test end-to-end execution of Bulk Transcriptomics downstream analysis:
    Raw counts -> QC & Filtering -> Normalization -> DE -> Viz -> Enrichment -> Batch -> Reporting.
    """
    set_seed(42)
    raw_counts, metadata = synthetic_data

    # --- Step 1: Data I/O & Schema Validation (FR-1) ---
    counts_file = tmp_path / "raw_counts.csv"
    meta_file = tmp_path / "sample_metadata.csv"
    export_table(raw_counts, counts_file, index=True)
    export_table(metadata, meta_file, index=True)

    dataset = load_dataset(counts_file, meta_file)
    loaded_counts, loaded_meta = dataset.counts, dataset.metadata
    val_report = dataset.validation_report or validate_bulk_data(loaded_counts, loaded_meta)
    assert val_report.is_valid
    assert val_report.n_genes == 100
    assert val_report.n_samples == 6

    # --- Step 2: QC & Normalization (FR-2) ---
    qc_stats = compute_sample_qc(loaded_counts)
    assert len(qc_stats) == 6
    assert "library_size" in qc_stats.columns

    filtered_counts, filter_summary = filter_low_expression_genes(
        loaded_counts, min_counts=5, min_samples=2
    )
    assert len(filtered_counts) <= 100
    assert not filter_summary.empty

    norm_res = normalize_deseq2(filtered_counts)
    norm_counts = norm_res.normalized_counts
    assert norm_counts.shape == filtered_counts.shape
    assert len(norm_res.size_factors) == 6

    # --- Step 3: Differential Expression Integration (FR-3) ---
    de_result = run_de(
        counts=filtered_counts,
        metadata=loaded_meta,
        contrast=("condition", "treated", "control"),
        alpha=0.05,
        lfc_threshold=1.0,
    )
    assert de_result is not None
    assert len(de_result.results_df) > 0
    up_genes = de_result.get_up_genes()
    down_genes = de_result.get_down_genes()
    assert len(up_genes) > 0
    assert len(down_genes) > 0

    # Acceptance Criteria #2: Verify underlying original library objects are accessible
    assert isinstance(de_result.dds, DeseqDataSet)
    assert isinstance(de_result.stat_res, DeseqStats)
    assert hasattr(de_result.dds, "varm")

    # --- Step 4: Publication Visualization (FR-4) ---
    fig_dir = tmp_path / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # 4a. Volcano Plot
    volcano_fig, _ = plot_volcano(
        de_result,
        interactive=False,
        top_n_labels=5,
        title="E2E Volcano Plot",
    )
    volcano_path = fig_dir / "volcano_plot.png"
    save_figure(volcano_fig, volcano_path)
    assert volcano_path.exists()

    # 4b. PCA Plot & Underlying PCA Model Check
    pca_model, pca_df = compute_pca(norm_counts, loaded_meta, n_components=2)
    assert isinstance(pca_model, PCA)
    assert len(pca_model.explained_variance_ratio_) == 2

    pca_fig, _, _ = plot_pca(
        norm_counts,
        metadata=loaded_meta,
        color_by="condition",
        title="E2E PCA Plot",
        interactive=False,
    )
    pca_path = fig_dir / "pca_plot.png"
    save_figure(pca_fig, pca_path)
    assert pca_path.exists()

    # 4c. Clustered Heatmap
    heatmap_grid = plot_heatmap(
        norm_counts,
        metadata=loaded_meta,
        de_result=de_result,
        top_n_degs=20,
        annotation_cols=["condition"],
        title="E2E Top DEGs Heatmap",
    )
    heatmap_path = fig_dir / "heatmap.png"
    save_figure(heatmap_grid.fig, heatmap_path)
    assert heatmap_path.exists()

    # --- Step 5: Functional / Pathway Enrichment (FR-5) ---
    # Construct a local mock gene set for deterministic offline test
    custom_gene_sets = {
        "PATHWAY_INFLAMMATION": up_genes[:8] + ["GENE_099"],
        "PATHWAY_METABOLISM": down_genes[:8] + ["GENE_098"],
        "PATHWAY_RANDOM": [f"GENE_{i:03d}" for i in range(50, 65)],
    }
    enr_result = run_custom_ora(
        gene_list=up_genes,
        gene_sets=custom_gene_sets,
        background=list(filtered_counts.index),
    )
    assert isinstance(enr_result, EnrichmentResult)
    assert len(enr_result.results_df) > 0
    assert "PATHWAY_INFLAMMATION" in enr_result.results_df["term"].values

    # --- Step 6: Batch Effect Correction (FR-7) ---
    combat_counts = run_combat(
        data=norm_counts,
        batch="batch",
        metadata=loaded_meta,
        biological_factor="condition",
        is_count=False,
    )
    assert combat_counts.shape == norm_counts.shape

    # --- Step 7: Reporting & Export Publication Bundle (FR-9) ---
    bundle_dir = tmp_path / "final_publication_bundle"
    figure_mapping = {
        "Volcano Plot": volcano_path,
        "PCA Analysis": pca_path,
        "Top DEGs Heatmap": heatmap_path,
    }

    out_bundle = export_publication_bundle(
        de_result=de_result,
        output_dir=bundle_dir,
        enrichment_results={"Custom_Signatures": enr_result},
        figure_paths=figure_mapping,
        bundle_name="End_to_End_Bulk_RNA_Report",
    )

    assert out_bundle.exists()
    assert (bundle_dir / "tables" / "differential_expression_master.csv").exists()
    assert (bundle_dir / "tables" / "complete_analysis_workbook.xlsx").exists()
    assert (bundle_dir / "figures" / "volcano_plot.png").exists()
    assert (bundle_dir / "figures" / "pca_plot.png").exists()
    assert (bundle_dir / "figures" / "heatmap.png").exists()
    assert (bundle_dir / "End_to_End_Bulk_RNA_Report.html").exists()
    assert (bundle_dir / "End_to_End_Bulk_RNA_Report.md").exists()

    # Verify content in generated reports
    html_content = (bundle_dir / "End_to_End_Bulk_RNA_Report.html").read_text(encoding="utf-8")
    md_content = (bundle_dir / "End_to_End_Bulk_RNA_Report.md").read_text(encoding="utf-8")

    assert "End to End Bulk RNA Report" in html_content
    assert "data:image/png;base64," in html_content
    assert "Executive Summary" in md_content
    assert "Engine:" in html_content
    assert "Execution Engine:" in md_content


def test_pipeline_reproducibility(synthetic_data):
    """
    Acceptance Criteria #3: Verify that repeating the pipeline under the same
    random seed produces strictly identical numerical DE and PCA results.
    """
    raw_counts, metadata = synthetic_data

    # Run 1
    set_seed(42)
    de1 = run_de(raw_counts, metadata, contrast=("condition", "treated", "control"))
    pca1, _ = compute_pca(raw_counts, metadata, n_components=2)

    # Run 2
    set_seed(42)
    de2 = run_de(raw_counts, metadata, contrast=("condition", "treated", "control"))
    pca2, _ = compute_pca(raw_counts, metadata, n_components=2)

    # Assert exact numerical equality across runs
    pd.testing.assert_frame_equal(de1.results_df, de2.results_df)
    np.testing.assert_allclose(pca1.explained_variance_ratio_, pca2.explained_variance_ratio_)
