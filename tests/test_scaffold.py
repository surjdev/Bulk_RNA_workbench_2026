"""
Smoke tests for BTW package scaffold, configuration, and reproducibility fixtures.
"""

import numpy as np
import pandas as pd
import pytest
from btw import __version__, logger, set_seed
from btw.config import AppConfig, load_config


def test_package_metadata():
    """Verify version string is present and valid."""
    assert __version__ == "0.1.0"
    assert logger is not None


def test_set_seed_reproducibility():
    """Verify that set_seed ensures reproducible random numbers."""
    set_seed(42)
    val1 = np.random.rand()
    set_seed(42)
    val2 = np.random.rand()
    assert val1 == val2

    # Different seed yields different results
    set_seed(99)
    val3 = np.random.rand()
    assert val1 != val3


def test_config_loader(temp_dir, sample_config):
    """Test loading configuration defaults and from YAML file."""
    # Test default config
    default_cfg = load_config(None)
    assert isinstance(default_cfg, AppConfig)
    assert default_cfg.project_name == "Bulk Transcriptomics Workbench"
    assert default_cfg.qc.min_library_size == 100000

    # Test overridden config from fixture
    assert sample_config.project_name == "Test Workbench"
    assert sample_config.random_seed == 123
    assert sample_config.qc.min_library_size == 50000
    assert sample_config.de_analysis.alpha == 0.01

    # Test non-existent path fallback
    non_existent = load_config(temp_dir / "does_not_exist.yaml")
    assert isinstance(non_existent, AppConfig)


def test_synthetic_data_fixture(synthetic_data):
    """Verify synthetic counts and metadata fixtures meet requirements."""
    counts, metadata = synthetic_data
    assert isinstance(counts, pd.DataFrame)
    assert isinstance(metadata, pd.DataFrame)

    assert counts.shape == (100, 6)
    assert metadata.shape == (6, 2)

    # Check sample IDs alignment
    assert list(counts.columns) == list(metadata.index)

    # Check counts are non-negative integers
    assert (counts >= 0).all().all()
    assert np.issubdtype(counts.values.dtype, np.integer)
