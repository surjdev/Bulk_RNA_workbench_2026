"""
Configuration management module for BTW.
Provides helper functions and dataclasses to load, validate, and access configuration settings.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional, Union

import yaml
from btw import logger


@dataclass
class QCConfig:
    min_library_size: int = 100000
    min_detected_genes: int = 500
    min_count_per_gene: int = 10
    min_samples_expressing: int = 3


@dataclass
class NormalizationConfig:
    method: str = "deseq2"
    fit_type: str = "parametric"


@dataclass
class DEConfig:
    design_factor: str = "condition"
    reference_level: str = "control"
    alpha: float = 0.05
    lfc_threshold: float = 1.0
    padj_method: str = "fdr_bh"


@dataclass
class VizConfig:
    style: str = "nature"
    dpi: int = 300
    font_family: str = "sans-serif"
    volcano_top_n_genes: int = 15
    heatmap_top_n: int = 50


@dataclass
class AppConfig:
    project_name: str = "Bulk Transcriptomics Workbench"
    version: str = "0.1.0"
    random_seed: int = 42
    qc: QCConfig = field(default_factory=QCConfig)
    normalization: NormalizationConfig = field(default_factory=NormalizationConfig)
    de_analysis: DEConfig = field(default_factory=DEConfig)
    visualization: VizConfig = field(default_factory=VizConfig)
    caching_enabled: bool = True
    cache_dir: str = ".btw_cache"


def load_config(config_path: Optional[Union[str, Path]] = None) -> AppConfig:
    """
    Load configuration from a YAML file, falling back to default values.

    Parameters
    ----------
    config_path : str or Path, optional
        Path to the YAML config file. If None, default values are returned.

    Returns
    -------
    AppConfig
        Populated configuration dataclass.
    """
    if config_path is None:
        logger.info("No config path provided; using default AppConfig.")
        return AppConfig()

    path = Path(config_path)
    if not path.exists():
        logger.warning(f"Config file not found at {path}; using default AppConfig.")
        return AppConfig()

    with open(path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}

    proj = data.get("project", {})
    qc = data.get("qc", {})
    norm = data.get("normalization", {})
    de = data.get("de_analysis", {})
    viz = data.get("visualization", {})
    cache = data.get("caching", {})

    app_cfg = AppConfig(
        project_name=proj.get("name", "Bulk Transcriptomics Workbench"),
        version=proj.get("version", "0.1.0"),
        random_seed=proj.get("random_seed", 42),
        qc=QCConfig(
            min_library_size=qc.get("min_library_size", 100000),
            min_detected_genes=qc.get("min_detected_genes", 500),
            min_count_per_gene=qc.get("min_count_per_gene", 10),
            min_samples_expressing=qc.get("min_samples_expressing", 3),
        ),
        normalization=NormalizationConfig(
            method=norm.get("method", "deseq2"),
            fit_type=norm.get("fit_type", "parametric"),
        ),
        de_analysis=DEConfig(
            design_factor=de.get("design_factor", "condition"),
            reference_level=de.get("reference_level", "control"),
            alpha=de.get("alpha", 0.05),
            lfc_threshold=de.get("lfc_threshold", 1.0),
            padj_method=de.get("padj_method", "fdr_bh"),
        ),
        visualization=VizConfig(
            style=viz.get("style", "nature"),
            dpi=viz.get("dpi", 300),
            font_family=viz.get("font_family", "sans-serif"),
            volcano_top_n_genes=viz.get("volcano", {}).get("top_n_genes", 15),
            heatmap_top_n=viz.get("heatmap", {}).get("top_n_variable_genes", 50),
        ),
        caching_enabled=cache.get("enabled", True),
        cache_dir=cache.get("cache_dir", ".btw_cache"),
    )

    logger.info(f"Loaded configuration from {path}")
    return app_cfg
