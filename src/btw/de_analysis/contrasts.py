"""
Multi-contrast batch execution and aggregation module (FR-3).
Enables running multiple biological comparisons in a loop and consolidating results into a single workbench.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union

import pandas as pd

from btw import logger
from btw.de_analysis.deseq_helper import DEResult, build_deseq_dataset, run_de, run_deseq_stats
from btw.io.exporter import export_excel_multisheet


def _format_contrast_name(contrast: Tuple[str, str, str]) -> str:
    """Generate a clean label for a contrast tuple, e.g. 'treated_vs_control'."""
    factor, test, ref = contrast
    return f"{test}_vs_{ref}"


@dataclass
class MultiContrastResult:
    """
    Container consolidating multiple Differential Expression comparisons.
    Provides cross-contrast DEG queries, summary tables, and unified export.
    """

    contrasts: Dict[str, DEResult] = field(default_factory=dict)
    master_table: Optional[pd.DataFrame] = None

    @property
    def contrast_names(self) -> List[str]:
        """List of all contrast labels in this collection."""
        return list(self.contrasts.keys())

    def get_contrast(self, name: str) -> DEResult:
        """Retrieve the DEResult for a specific contrast."""
        if name not in self.contrasts:
            raise KeyError(f"Contrast '{name}' not found. Available: {self.contrast_names}")
        return self.contrasts[name]

    def summary_table(self) -> pd.DataFrame:
        """
        Produce a summary table comparing DEG numbers across all contrasts.

        Returns
        -------
        pd.DataFrame
            Overview table containing contrast details and DEG counts.
        """
        rows = []
        for name, de_res in self.contrasts.items():
            factor, test, ref = de_res.contrast
            counts = de_res.summary_counts()
            rows.append(
                {
                    "contrast_name": name,
                    "factor": factor,
                    "test_level": test,
                    "ref_level": ref,
                    "alpha": de_res.alpha,
                    "lfc_threshold": de_res.lfc_threshold,
                    "total_tested": counts["total_tested"],
                    "sig_up": counts["significant_up"],
                    "sig_down": counts["significant_down"],
                    "total_degs": counts["total_degs"],
                }
            )
        return pd.DataFrame(rows)

    def get_deg_sets(
        self,
        padj_cutoff: Optional[float] = None,
        lfc_cutoff: Optional[float] = None,
        direction: str = "both",
    ) -> Dict[str, Set[str]]:
        """
        Extract DEG gene sets per contrast for Venn and UpSet plot comparisons.

        Parameters
        ----------
        padj_cutoff : float, optional
            P-adjusted cutoff.
        lfc_cutoff : float, optional
            Log2FC cutoff.
        direction : {'both', 'up', 'down'}, default='both'
            Filter by regulation direction.

        Returns
        -------
        dict of {str: set of str}
            Mapping of contrast name to set of significant gene IDs.
        """
        deg_sets = {}
        for name, de_res in self.contrasts.items():
            if direction == "up":
                genes = set(de_res.get_up_genes(padj_cutoff, lfc_cutoff))
            elif direction == "down":
                genes = set(de_res.get_down_genes(padj_cutoff, lfc_cutoff))
            else:
                genes = set(de_res.get_degs(padj_cutoff, lfc_cutoff).index)
            deg_sets[name] = genes
        return deg_sets

    def export_excel(self, output_path: Union[str, Path]) -> Path:
        """
        Export all contrast DE tables and summary into a single multi-sheet Excel file.

        Parameters
        ----------
        output_path : str or Path
            Destination .xlsx file path.

        Returns
        -------
        Path
            Path to the saved workbook.
        """
        sheets = {"Summary": self.summary_table()}
        if self.master_table is not None:
            sheets["Master_Table"] = self.master_table

        for name, de_res in self.contrasts.items():
            sheets[f"DE_{name}"[:31]] = de_res.results_df

        return export_excel_multisheet(sheets, output_path)


def run_multiple_contrasts(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    contrasts: List[Tuple[str, str, str]],
    design_factors: Optional[Union[str, List[str]]] = None,
    alpha: float = 0.05,
    lfc_threshold: float = 1.0,
    engine: str = "r",
    method: str = "deseq2",
    fallback_to_python: bool = True,
    **kwargs,
) -> MultiContrastResult:
    """
    Execute multiple differential expression contrasts in batch and assemble a master table.
    Supports both reference R engines (DESeq2, limma) and Python PyDESeq2 engine.

    Parameters
    ----------
    counts : pd.DataFrame
        Gene count matrix (n_genes, n_samples).
    metadata : pd.DataFrame
        Sample metadata.
    contrasts : list of tuple of (factor, test_level, ref_level)
        List of experimental contrasts to test.
    design_factors : str or list of str, optional
        GLM design factors. If None, inferred from factors present in contrasts.
    alpha : float, default=0.05
        Significance threshold.
    lfc_threshold : float, default=1.0
        Log2 fold change threshold.
    engine : str, default='r'
        'r' (calls R DESeq2/limma via rpy2) or 'python' (PyDESeq2).
    method : str, default='deseq2'
        'deseq2', 'limma_voom', or 'limma_trend'.
    fallback_to_python : bool, default=True
        Whether to fall back to PyDESeq2 if R packages are missing.
    **kwargs
        Additional arguments passed to engine.

    Returns
    -------
    MultiContrastResult
        Aggregated results across all contrasts.
    """
    if not contrasts:
        raise ValueError("Must provide at least one contrast tuple.")

    # Deduplicate design factors needed
    if design_factors is None:
        needed_factors = list(dict.fromkeys([c[0] for c in contrasts]))
        design = needed_factors[0] if len(needed_factors) == 1 else needed_factors
    else:
        design = design_factors

    results_dict: Dict[str, DEResult] = {}
    master_cols = []

    # If pure Python engine requested explicitly
    if engine.lower() == "python":
        logger.info(f"Preparing shared PyDESeq2 DeseqDataSet for {len(contrasts)} contrast(s)...")
        dds = build_deseq_dataset(
            counts=counts,
            metadata=metadata,
            design_factors=design,
            fit_model=True,
            **kwargs,
        )
        for c in contrasts:
            name = _format_contrast_name(c)
            logger.info(f"Running contrast: {name} {c}")
            de_res = run_deseq_stats(
                dds=dds,
                contrast=c,
                alpha=alpha,
                lfc_threshold=lfc_threshold,
                **kwargs,
            )
            results_dict[name] = de_res
            c_df = de_res.results_df[["log2FoldChange", "padj", "regulation"]].copy()
            c_df.columns = [f"{col}_{name}" for col in c_df.columns]
            master_cols.append(c_df)
    else:
        # R engine per contrast (or fallback handled inside run_de)
        for c in contrasts:
            name = _format_contrast_name(c)
            logger.info(f"Running contrast via run_de (engine={engine}): {name} {c}")
            de_res = run_de(
                counts=counts,
                metadata=metadata,
                contrast=c,
                design_factors=design,
                alpha=alpha,
                lfc_threshold=lfc_threshold,
                engine=engine,
                method=method,
                fallback_to_python=fallback_to_python,
                **kwargs,
            )
            results_dict[name] = de_res
            c_df = de_res.results_df[["log2FoldChange", "padj", "regulation"]].copy()
            c_df.columns = [f"{col}_{name}" for col in c_df.columns]
            master_cols.append(c_df)

    # Combine into master table
    master_df = pd.concat(master_cols, axis=1)

    multi_res = MultiContrastResult(
        contrasts=results_dict,
        master_table=master_df,
    )

    logger.info(
        f"Completed multi-contrast analysis: {len(results_dict)} contrast(s) evaluated. "
        f"Master table has {master_df.shape[1]} metrics across {master_df.shape[0]} genes."
    )
    return multi_res
