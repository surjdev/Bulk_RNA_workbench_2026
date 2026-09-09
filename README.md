# Downstream Bulk Transcriptomics Workbench (`btw`)

[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)

ชุดเครื่องมือ Python (`btw`) สำหรับงานวิเคราะห์ downstream ของ **Bulk RNA-seq / Transcriptomics** บน Jupyter Notebook ออกแบบตามข้อกำหนดใน [`SRS_Downstream_Bulk_Transcriptomics_Workbench.md`](SRS_Downstream_Bulk_Transcriptomics_Workbench.md)

---

## 🎯 ปรัชญาการออกแบบ (Design Philosophy)
- **Thin Utility Layer (ไม่ห่อหุ้มจนปิดกั้น):** ลดงานซ้ำ (boilerplate) ในการเตรียมข้อมูล, วิเคราะห์, พล็อตภาพ และทำรายงาน แต่ผู้ใช้ยังคงเข้าถึงและควบคุม object ดั้งเดิม (`pydeseq2`, `scikit-learn`, `scipy`, `networkx`) ได้เสมอ
- **Schema & Validation First:** ตรวจสอบความถูกต้องของ count matrix และ sample metadata ก่อนประมวลผลเสมอ ป้องกัน silent errors
- **Publication-Ready:** กราฟและตารางผลลัพธ์จัดแต่งตามมาตรฐานวารสารวิชาการชั้นนำ (Nature/Cell style palettes, 300+ DPI, non-overlapping labels)
- **Multi-Method Pathway Enrichment:** รองรับทั้ง ORA (Enrichr, goatools), GSEA (prerank), และ Activity inference (decoupler, PROGENy, DoRothEA) ภายใต้ unified output format
- **High Reproducibility:** รองรับการจัดการ environment ด้วยทั้ง Conda (`environment.yml`) และ Pixi (`pixi.toml`) พร้อมตั้งค่า random seed จากส่วนกลาง

---

## 📂 โครงสร้างโปรเจกต์ (Project Structure)

```text
bulk-transcriptomics-workbench/
├── environment.yml                  # Conda/Mamba dependencies (Python 3.11)
├── pyproject.toml                   # Package metadata & tool configurations
├── pixi.toml                        # Fast environment & task runner
├── HANDOFF.md                       # AI Handoff protocol & milestone tracking
├── SRS_Downstream_Bulk_Transcriptomics_Workbench.md
├── src/
│   └── btw/                         # แพ็กเกจหลัก
│       ├── __init__.py              # Central logger & set_seed
│       ├── config.py                # Config management
│       ├── io/                      # FR-1: Loader, Validator, Exporter
│       ├── qc_normalize/            # FR-2: QC summaries & normalization helpers
│       ├── de_analysis/             # FR-3: Differential expression & multi-contrasts
│       ├── viz/                     # FR-4: Publication-grade biological visualizations
│       ├── enrichment/              # FR-5: ORA, GSEA, decoupler & unified dotplots
│       ├── annotation/              # FR-6: Gene ID mapping & GTF parsing
│       ├── batch_correction/        # FR-7: ComBat batch correction & diagnostic PCA
│       ├── network/                 # FR-8: WGCNA & Cytoscape network export
│       └── reporting/               # FR-9: Automated summaries & multi-format export
├── notebooks/
│   └── examples/                    # Example notebooks แยกตามแต่ละ module
├── configs/                         # YAML configuration defaults
└── tests/                           # Unit & integration test suite
```

---

## 🚀 การติดตั้งและเริ่มต้นใช้งาน (Installation)

### ตัวเลือกที่ 1: ติดตั้งผ่าน Pixi (แนะนำ)
```bash
# ติดตั้ง dependencies และเริ่ม JupyterLab
pixi run lab

# รันชุดการทดสอบ
pixi run test
```

### ตัวเลือกที่ 2: ติดตั้งผ่าน Conda / Mamba
```bash
conda env create -f environment.yml
conda activate btw_env
pip install -e .
```

---

## 🧪 การทดสอบ (Testing & Code Quality)

```bash
# รัน Unit Tests พร้อมวัด Coverage
pytest tests/ -v --cov=src/btw

# ตรวจสอบ Code Quality
ruff check src/ tests/
```
