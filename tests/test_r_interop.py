"""
Unit tests for BTW R-Python Interoperability (FR-10).
Validates bridge detection, type conversions, error diagnostics, and graceful fallback to Python engines.
"""

import numpy as np
import pandas as pd
import pytest

from btw.batch_correction.combat import run_combat
from btw.de_analysis.deseq_helper import run_de
from btw.enrichment.ora import run_clusterprofiler
from btw.network.wgcna_helper import run_wgcna
from btw.qc_normalize.normalize import normalize_deseq2
from btw.r_interop.bridge import (
    MissingRPackageError,
    check_r_package,
    get_r_version,
    is_r_available,
    pandas_to_r_df,
    pandas_to_r_matrix,
    r_to_pandas_df,
    require_r_package,
    set_r_seed,
)


def test_r_availability():
    """Verify that R runtime is detected when installed."""
    assert is_r_available() is True
    version = get_r_version()
    assert version is not None
    assert len(version.split(".")) >= 2


def test_r_package_checks():
    """Verify package checking and require_r_package behavior."""
    assert check_r_package("base") is True
    assert check_r_package("stats") is True
    assert check_r_package("__definitely_not_a_real_r_package__") is False

    with pytest.raises(MissingRPackageError) as exc_info:
        require_r_package("__definitely_not_a_real_r_package__", purpose="Testing diagnostics")
    assert "__definitely_not_a_real_r_package__" in str(exc_info.value)
    assert "Testing diagnostics" in str(exc_info.value)
    assert "install.packages" in str(exc_info.value) or "BiocManager" in str(exc_info.value)


def test_r_conversions():
    """Verify roundtrip conversion between pandas DataFrame and R data.frame."""
    df = pd.DataFrame(
        {
            "gene_a": [10.5, 20.2, 30.1],
            "gene_b": [5.0, 15.0, 25.0],
        },
        index=["sample_1", "sample_2", "sample_3"],
    )
    df.index.name = "sample_id"

    r_df = pandas_to_r_df(df)
    back_df = r_to_pandas_df(r_df)

    assert back_df.shape == df.shape
    assert list(back_df.columns) == list(df.columns)
    assert list(back_df.index) == list(df.index)
    np.testing.assert_allclose(back_df.values, df.values, rtol=1e-5)


def test_r_matrix_conversion():
    """Verify pandas DataFrame to R matrix conversion."""
    import rpy2.robjects as ro

    df = pd.DataFrame(
        {
            "s1": [100, 200],
            "s2": [300, 400],
        },
        index=["g1", "g2"],
    )
    r_mat = pandas_to_r_matrix(df)
    dims = tuple(ro.r("dim")(r_mat))
    assert dims == (2, 2)
    colnames = list(ro.r("colnames")(r_mat))
    assert colnames == ["s1", "s2"]
    rownames = list(ro.r("rownames")(r_mat))
    assert rownames == ["g1", "g2"]


def test_set_r_seed():
    """Verify that setting R seed executes without error and affects R randomness."""
    import rpy2.robjects as ro

    set_r_seed(123)
    val1 = float(ro.r("runif(1)")[0])

    set_r_seed(123)
    val2 = float(ro.r("runif(1)")[0])

    assert val1 == val2


def test_normalize_deseq2_fallback(synthetic_data):
    """Verify normalize_deseq2 engine routing and fallback."""
    counts, metadata = synthetic_data
    # With fallback_to_python=True, succeeds even if R DESeq2 is absent
    norm_res = normalize_deseq2(counts, metadata, engine="r", fallback_to_python=True)
    assert norm_res.normalized_counts.shape == counts.shape
    assert (norm_res.normalized_counts >= 0).all().all()

    # With fallback_to_python=False and missing DESeq2, raises MissingRPackageError
    with pytest.raises(MissingRPackageError):
        normalize_deseq2(counts, metadata, engine="r", fallback_to_python=False)


def test_run_de_engine_routing(synthetic_data):
    """Verify run_de engine switching and fallback behavior."""
    counts, metadata = synthetic_data
    # 1. Explicit Python engine
    res_py = run_de(counts, metadata, contrast=("condition", "treated", "control"), engine="python")
    assert res_py.engine == "python"
    assert res_py.method == "pydeseq2"
    assert "pvalue" in res_py.results_df.columns
    assert "padj" in res_py.results_df.columns

    # 2. R engine with fallback
    res_r_fallback = run_de(
        counts,
        metadata,
        contrast=("condition", "treated", "control"),
        engine="r",
        fallback_to_python=True,
    )
    assert res_r_fallback is not None
    assert "padj" in res_r_fallback.results_df.columns

    # 3. R engine without fallback raises MissingRPackageError
    with pytest.raises(MissingRPackageError):
        run_de(
            counts,
            metadata,
            contrast=("condition", "treated", "control"),
            engine="r",
            fallback_to_python=False,
        )


def test_combat_engine_fallback(synthetic_data):
    """Verify ComBat engine routing and fallback."""
    counts, metadata = synthetic_data
    # R engine with fallback succeeds
    corrected = run_combat(
        data=counts,
        batch=metadata["batch"],
        metadata=metadata,
        biological_factor="condition",
        engine="r",
        fallback_to_python=True,
    )
    assert corrected.shape == counts.shape

    # Without fallback, raises MissingRPackageError
    with pytest.raises(MissingRPackageError):
        run_combat(
            data=counts,
            batch=metadata["batch"],
            metadata=metadata,
            biological_factor="condition",
            engine="r",
            fallback_to_python=False,
        )


def test_wgcna_engine_fallback(synthetic_data):
    """Verify WGCNA engine routing and fallback."""
    counts, _ = synthetic_data
    # Subset to 20 genes for quick test
    small_counts = counts.iloc[:20]

    res = run_wgcna(small_counts, engine="r", fallback_to_python=True, min_module_size=3)
    assert res is not None
    assert hasattr(res, "module_labels")
    assert len(res.module_labels) == 20

    with pytest.raises(MissingRPackageError):
        run_wgcna(small_counts, engine="r", fallback_to_python=False)


def test_clusterprofiler_fallback():
    """Verify clusterProfiler engine routing and fallback to Python ORA."""
    query_genes = ["GENE_001", "GENE_002", "GENE_003"]
    custom_sets = {
        "PATH_A": ["GENE_001", "GENE_002", "GENE_004"],
        "PATH_B": ["GENE_005", "GENE_006"],
    }

    # Fallback to Python custom ORA
    res = run_clusterprofiler(
        gene_list=query_genes,
        fallback_to_python=True,
        fallback_gene_sets=custom_sets,
    )
    assert res is not None
    assert not res.results_df.empty
    assert "PATH_A" in res.results_df["term"].values

    # Without fallback, raises MissingRPackageError
    with pytest.raises(MissingRPackageError):
        run_clusterprofiler(gene_list=query_genes, fallback_to_python=False)
