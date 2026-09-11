# Downstream Bulk Transcriptomics Workbench (`btw`)

[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.14-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.3%2B-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)

A publication-grade downstream **Bulk Transcriptomics Workbench (`btw`)** built for Jupyter, supporting a **Hybrid R/Python Architecture** in compliance with [`SRS_Downstream_Bulk_Transcriptomics_Workbench_v2.md`](SRS_Downstream_Bulk_Transcriptomics_Workbench_v2.md).

---

## 🎯 Architecture & Design Philosophy (SRS v2)

1. **Hybrid R/Python on Jupyter (Language Matrix):**
   - If an R/Bioconductor package is the gold-standard reference implementation (`DESeq2`, `limma-voom`/`limma-trend`, `sva::ComBat`, `WGCNA`, `clusterProfiler`), `btw` invokes it via `rpy2` as default (`engine="r"`).
   - Pure-Python implementations (`pydeseq2`, `inmoose`, `PyWGCNA`, `gseapy`, `decoupler`) serve as fully supported alternatives (`engine="python"`).
   - **Seamless Fallback:** If R or specific Bioconductor packages are missing in the user's environment and `fallback_to_python=True` (default in config), `btw` automatically logs and falls back to the Python alternative without crashing (NFR-10, Acceptance Criteria #6).
2. **No Object Concealment (Acceptance Criteria #2):**
   - `btw` never conceals native R objects (`DESeqDataSet`, `MArrayLM`) or Python objects (`DeseqDataSet`, `DeseqStats`, scikit-learn `PCA`, `networkx.Graph`). All handles remain directly accessible.
3. **Dataframe In / Dataframe Out:**
   - All R functions accept standard `pandas.DataFrame` and return clean `pandas.DataFrame` results without requiring researchers to write R scripts.
4. **Engine Provenance Tracking (FR-9):**
   - Every analysis artifact, HTML report, and Markdown summary includes clear badges showing the execution engine (e.g. `Engine: R (DESeq2 v1.42.0)` or `Engine: PYTHON (pydeseq2 v0.4.11)`).

---

## 📊 Analytical Method Matrix

| Task | R Reference Engine (`engine="r"`) | Python Alternative (`engine="python"`) | Fallback Mechanism |
|:---|:---|:---|:---|
| **Differential Expression** | `DESeq2` (Wald test + apeglm/ashr), `limma-voom`, `limma-trend` | `pydeseq2` | Automatic fallback to `pydeseq2` |
| **Size-Factor Normalization** | `DESeq2::estimateSizeFactors` | Median-of-Ratios (`compute_size_factors`) | Automatic fallback to Python |
| **Variance Stabilization** | `DESeq2::vst` | Log2(norm + 1) or Python VST | Automatic fallback to log2 |
| **Batch Correction** | `sva::ComBat` & `sva::ComBat_seq` | `inmoose.pycombat` | Automatic fallback to `inmoose` |
| **Co-expression Network** | `WGCNA::blockwiseModules` | Python co-expression helper & `PyWGCNA` | Automatic fallback to Python helper |
| **Pathway Enrichment** | `clusterProfiler::enrichGO` / `enrichKEGG` | `gseapy.enrichr`, custom Fisher ORA, `decoupler` | Automatic fallback to Python ORA |

---

## 📂 Project Structure

```text
bulk-transcriptomics-workbench/
├── environment.yml                  # Conda/Mamba dependencies (Python + R + Bioconductor)
├── pyproject.toml                   # Package metadata & optional extras [r-stats]
├── configs/
│   └── default_config.yaml          # YAML configuration with engine & fallback settings
├── src/
│   └── btw/                         # Core Python workbench package
│       ├── __init__.py              # Central logger & synchronized set_seed (Python + R)
│       ├── config.py                # Config schemas (AppConfig, DEConfig, BatchConfig, NetworkConfig)
│       ├── r_interop/               # FR-10: Isolated R/Bioconductor bridges
│       │   ├── bridge.py            # Package checks, type conversions, MissingRPackageError
│       │   ├── deseq2_r.py          # R DESeq2 Wald test & LFC shrinkage
│       │   ├── limma_r.py           # R limma-voom & limma-trend
│       │   ├── combat_r.py          # R sva::ComBat & ComBat_seq
│       │   ├── wgcna_r.py           # R WGCNA::blockwiseModules
│       │   └── clusterprofiler_r.py # R clusterProfiler::enrichGO
│       ├── io/                      # FR-1: Loader, Validator, Exporter
│       ├── qc_normalize/            # FR-2: QC summaries & normalization (R / Python)
│       ├── de_analysis/             # FR-3: Wald test, limma, multi-contrasts
│       ├── viz/                     # FR-4: Publication-grade visualizations
│       ├── enrichment/              # FR-5: ORA (clusterProfiler / Enrichr), GSEA, decoupler
│       ├── annotation/              # FR-6: Gene ID mapper & GTF parser
│       ├── batch_correction/        # FR-7: ComBat (sva / inmoose) & diagnostic PCA
│       ├── network/                 # FR-8: WGCNA (R / Python) & Cytoscape export
│       └── reporting/               # FR-9: Automated summaries & publication bundle
├── notebooks/
│   └── examples/
│       └── 00_end_to_end_pipeline.ipynb # Complete runnable hybrid pipeline
└── tests/                           # Pytest test suite (63 unit & integration tests)
```

---

## 🚀 Installation & Quick Start

### Option 1: Full Hybrid Environment (Conda / Mamba - Recommended)

```bash
# Clone repository
git clone https://github.com/surjdev/Bulk_RNA_workbench_2026.git
cd Bulk_RNA_workbench_2026

# Create environment with Python 3.11+, R 4.3+, and Bioconductor packages
conda env create -f environment.yml
conda activate btw_env

# Install BTW in editable mode with R interop
pip install -e ".[r-stats]"
```

### Option 2: Python-Only Environment (Minimal)

```bash
pip install -e .
```
> When running in Python-only mode without R, all analysis steps automatically execute using high-performance Python engines (`pydeseq2`, `inmoose`, `gseapy`, `decoupler`).

### Installing Bioconductor Packages in Existing R

If you already have R 4.3+ installed, install the reference Bioconductor libraries:
```R
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "limma", "edgeR", "sva", "WGCNA", "clusterProfiler", "apeglm"))
```

---

## 🧪 Quick Usage Example

```python
import pandas as pd
from btw import set_seed
from btw.de_analysis import run_de

# 1. Synchronize random seed across Python and R
set_seed(42)

# 2. Run Differential Expression using R reference engine (with automatic Python fallback)
de_result = run_de(
    counts=counts_df,
    metadata=metadata_df,
    contrast=("condition", "treated", "control"),
    engine="r",                # "r" or "python"
    method="deseq2",           # "deseq2", "limma_voom", "limma_trend"
    fallback_to_python=True,   # Gracefully falls back if R package missing
)

print(de_result.summary())

# Direct access to native reference objects
if de_result.engine == "r":
    print("R DESeqDataSet handle:", de_result.dds)
else:
    print("PyDESeq2 DeseqDataSet:", de_result.dds)
```

---

## 🔬 Testing & Quality Assurance

```bash
# Run all tests with coverage report
pytest tests/ -v --cov=src/btw --cov-report=term-missing

# Run R-Python interoperability tests specifically
pytest tests/test_r_interop.py -v
```

---

## 📄 License

MIT License. Designed for reproducible, rigorous transcriptomics research.
