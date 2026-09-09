"""
Publication and Presentation Export utilities for BTW (FR-9).
Compiles tables, high-resolution figures, and reports into a unified deliverable directory.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Dict, Optional, Union

import pandas as pd
from btw import logger
from btw.de_analysis.contrasts import DEResult
from btw.enrichment.schema import EnrichmentResult
from btw.io.exporter import export_table
from btw.reporting.report import generate_html_report, generate_markdown_report


def export_publication_bundle(
    de_result: DEResult,
    output_dir: Union[str, Path] = "results/publication_bundle",
    enrichment_results: Optional[Dict[str, EnrichmentResult]] = None,
    figure_paths: Optional[Dict[str, Union[str, Path]]] = None,
    bundle_name: str = "Bulk_RNA_Analysis_Report",
) -> Path:
    """
    Compile a complete publication bundle containing:
    1. Multi-sheet Excel workbook with all DE results and enrichment tables.
    2. CSV / TSV copies for programmatic consumption.
    3. High-resolution figure directory (copying PNG/SVG/PDF files).
    4. Standalone Executive HTML Report with embedded interactive cards.
    5. Executive Markdown Report.

    Parameters
    ----------
    de_result : DEResult
        Differential expression result object.
    output_dir : str or Path, default='results/publication_bundle'
        Root directory for the bundle.
    enrichment_results : dict of {str: EnrichmentResult}, optional
        Collection of pathway enrichment results.
    figure_paths : dict of {str: str/Path}, optional
        Figure captions mapped to image file paths.
    bundle_name : str, default='Bulk_RNA_Analysis_Report'
        Base filename for generated documents.

    Returns
    -------
    Path
        Path to the generated bundle directory.
    """
    bundle_path = Path(output_dir)
    tables_dir = bundle_path / "tables"
    figs_dir = bundle_path / "figures"

    bundle_path.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    figs_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Assembling publication bundle at: {bundle_path.resolve()}...")

    # 1. Export DE Master Table
    de_csv = tables_dir / "differential_expression_master.csv"
    export_table(de_result.results_df, de_csv, index=True)

    # 2. Multi-sheet Excel
    excel_path = tables_dir / "complete_analysis_workbook.xlsx"
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        de_result.results_df.to_excel(writer, sheet_name="All_Genes")
        de_result.up_regulated().to_excel(writer, sheet_name="Significant_UP")
        de_result.down_regulated().to_excel(writer, sheet_name="Significant_DOWN")

        if enrichment_results:
            for name, enr in enrichment_results.items():
                sheet_safe = name[:30].replace("/", "_")
                enr.results_df.to_excel(writer, sheet_name=sheet_safe, index=False)

    logger.info(f"Saved complete multi-sheet Excel to: {excel_path}")

    # 3. Copy figures to bundle
    copied_figures = {}
    if figure_paths:
        for caption, src_path in figure_paths.items():
            src = Path(src_path)
            if src.exists():
                dst = figs_dir / src.name
                shutil.copy2(src, dst)
                copied_figures[caption] = dst
            else:
                logger.warning(f"Figure file not found: {src}")

    # 4. Generate Markdown & HTML Reports
    first_enr = list(enrichment_results.values())[0] if enrichment_results else None

    md_path = bundle_path / f"{bundle_name}.md"
    generate_markdown_report(
        title=bundle_name.replace("_", " "),
        de_result=de_result,
        enrichment_result=first_enr,
        figure_paths=copied_figures,
        output_path=md_path,
    )

    html_path = bundle_path / f"{bundle_name}.html"
    generate_html_report(
        title=bundle_name.replace("_", " "),
        de_result=de_result,
        enrichment_result=first_enr,
        figure_paths=copied_figures,
        embed_images=True,
        output_path=html_path,
    )

    logger.info(f"Publication bundle successfully assembled with {len(copied_figures)} figures and reports.")
    return bundle_path
