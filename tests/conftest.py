"""
Global Pytest fixtures for Bulk Transcriptomics Workbench (`btw`).
Provides synthetic count matrices, metadata, and temporary test assets.
"""

import tempfile
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
import pytest

from btw import set_seed
from btw.config import AppConfig, load_config


@pytest.fixture(autouse=True)
def reset_random_seed():
    """Reset random seed before every test for absolute reproducibility."""
    set_seed(42)


@pytest.fixture
def synthetic_data() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Generate synthetic bulk RNA-seq count matrix and metadata.
    Counts: 100 genes x 6 samples (3 control, 3 treated).
    Injects 15 up-regulated genes and 15 down-regulated genes in treated.
    """
    np.random.seed(42)
    genes = [f"GENE_{i:03d}" for i in range(1, 101)]
    samples = ["ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3"]

    # Baseline counts drawn from negative binomial distribution
    base_counts = np.random.negative_binomial(5, 0.01, size=(100, 6))

    # Inject differential expression
    # First 15 genes up-regulated in treated (treat_1, treat_2, treat_3)
    base_counts[:15, 3:] = (base_counts[:15, 3:] * 4.5).astype(int)
    # Next 15 genes down-regulated in treated
    base_counts[15:30, 3:] = (base_counts[15:30, 3:] * 0.2).astype(int)

    # Ensure non-negative integers
    counts = np.clip(base_counts, a_min=0, a_max=None).astype(int)

    count_df = pd.DataFrame(counts, index=genes, columns=samples)

    metadata_df = pd.DataFrame(
        {
            "sample_id": samples,
            "condition": ["control", "control", "control", "treated", "treated", "treated"],
            "batch": ["batch1", "batch2", "batch1", "batch2", "batch1", "batch2"],
        }
    ).set_index("sample_id")

    return count_df, metadata_df


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test file I/O."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_config(temp_dir: Path) -> AppConfig:
    """Create a sample YAML configuration file and return AppConfig."""
    cfg_yaml = temp_dir / "test_config.yaml"
    cfg_yaml.write_text(
        """
project:
  name: "Test Workbench"
  version: "0.1.0"
  random_seed: 123
qc:
  min_library_size: 50000
  min_detected_genes: 200
de_analysis:
  design_factor: "condition"
  reference_level: "control"
  alpha: 0.01
"""
    )
    return load_config(cfg_yaml)


@pytest.fixture
def de_result_fixture(synthetic_data):
    """Generate fitted DE result for visualization, enrichment, and annotation testing."""
    from btw.de_analysis import run_de

    counts, metadata = synthetic_data
    return run_de(counts, metadata, contrast=("condition", "treated", "control"))
