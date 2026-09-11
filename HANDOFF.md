# 🤝 AI Handoff & Milestone Verification Protocol
## Project: Bulk Transcriptomics Workbench (`btw`)
**Repository:** `/home/surj/Workspace/Bulk_RNA_workbench_2026`  
**Specification:** [`SRS_Downstream_Bulk_Transcriptomics_Workbench.md`](SRS_Downstream_Bulk_Transcriptomics_Workbench.md)  
**Implementation Plan:** [`implementation_plan.md`](file:///home/surj/.gemini/antigravity-ide/brain/14e29b1d-042d-4e62-ba2b-e7f922feeba7/implementation_plan.md)  
**Last Updated:** 2026-09-09  

---

## 1. Milestone Tracking Matrix

| Part | Milestone Description | Scope (FR/NFR) | Status | Test Coverage | Sign-off Status |
|:---:|:---|:---|:---:|:---:|:---:|
| **Part 0** | Scaffold, Environment & Foundation | NFR-1, 3, 5, 6, 8 | ✅ Completed | 100% (Scaffold) | Verified & Passed |
| **Part 1** | Data I/O, Schema Validation & QC | FR-1, FR-2 | ✅ Completed | 82% (Module) | Verified & Passed |
| **Part 2** | Differential Expression Integration | FR-3 | ✅ Completed | 93% (Module) | Verified & Passed |
| **Part 3** | Publication-Grade Visualization | FR-4 | ✅ Completed | 85% (Module), 84% (Total) | Verified & Passed |
| **Part 4** | Pathway & Functional Enrichment | FR-5 (5.1, 5.2, 5.3) | ✅ Completed | 74% (Module), 81% (Total) | Verified & Passed |
| **Part 5** | Annotation, Batch & Co-expression | FR-6, FR-7, FR-8 | ✅ Completed | 82% (Module), 80% (Total) | Verified & Passed |
| **Part 6** | Reporting & End-to-End Workbench | FR-9, Acceptance | ✅ Completed | 96% (Module), 82% (Total) | Verified & Passed |
| **Part 7** | SRS v2 Rebuild (Hybrid R/Python Architecture) | FR-10, NFR-10, All FRs | ✅ Completed | 74% (Full Project, 2,689 stmts) | Verified & Passed |

*Status options: `⏳ Pending` \| `🔄 In Progress` \| `✅ Completed` \| `⚠️ Blocked`*

---

## 2. AI Pre-Flight Inspection Protocol (สิ่งที่ AI ต้องตรวจก่อนเริ่มงานทุกครั้ง)

เมื่อ AI เข้ามารับงานต่อในแต่ละ Part หรือเริ่มต้นเซสชันใหม่ **ต้องทำตามขั้นตอนตรวจสอบต่อไปนี้ก่อนเริ่มลงมือเขียนโค้ด:**

1. **ตรวจสอบ Git Working Tree:**
   ```bash
   git status
   ```
   *ตรวจดูว่ามี uncommitted changes หรือ conflict ตกค้างหรือไม่*

2. **ตรวจสอบ Python Environment & Dependencies:**
   ```bash
   python3 --version
   pixi --version || conda --version
   ```
   *ตรวจสอบว่ากำลังทำงานบน Python 3.10+ (แนะนำ 3.11) ที่ติดตั้ง dependencies ครบถ้วน*

3. **รัน Test Suite เดิมเพื่อเช็ค Regression:**
   ```bash
   pytest tests/ -v
   ```
   *ต้องมั่นใจว่าการทดสอบของ Part ก่อนหน้าผ่าน 100% (ปัจจุบัน 31/31 passed) ก่อนเริ่ม Part ถัดไปเสมอ*

4. **ตรวจสอบสถานะในตาราง Milestone:**
   *ดูว่า Part ปัจจุบันคือ Part ใด และผู้ใช้ส่งคำสั่งให้ทำ Part ไหน*

---

## 3. Part-by-Part Verification & Acceptance Gates

### 🎯 Part 0: Project Scaffold & Environment Foundation
- **Deliverables:**
  - `pyproject.toml`, `environment.yml`, `pixi.toml`
  - `.pre-commit-config.yaml`, `.gitignore`
  - `README.md` (ปรับปรุงชื่อและเนื้อหาโครงการ)
  - `src/btw/__init__.py` (version, `set_seed()`, `logger` ด้วย `loguru`/`logging`)
  - `src/btw/config.py`, `configs/default_config.yaml`
  - `tests/conftest.py`, `tests/test_scaffold.py`
- **คำสั่งตรวจสอบ (Verification Command):**
  ```bash
  pytest tests/test_scaffold.py -v
  ```
- **Acceptance Gate Checklist:**
  - [x] Package ติดตั้งแบบ editable ได้ (`pip install -e .` หรือ `pixi run`)
  - [x] `import btw` สำเร็จ และ `btw.set_seed(42)` กำหนด seed ได้จากส่วนกลาง
  - [x] Synthetic bulk RNA count matrix และ metadata fixture ใน `conftest.py` พร้อมใช้งาน

---

### 🎯 Part 1: Data I/O & Preprocessing (FR-1, FR-2)
- **Deliverables:**
  - `src/btw/io/{loader,validator,exporter}.py`
  - `src/btw/qc_normalize/{qc,normalize}.py`
  - `tests/test_io.py`, `tests/test_qc_normalize.py`
  - `notebooks/examples/01_io_and_qc.ipynb`
- **คำสั่งตรวจสอบ (Verification Command):**
  ```bash
  pytest tests/test_io.py tests/test_qc_normalize.py -v --cov=src/btw/io --cov=src/btw/qc_normalize
  ```
- **Acceptance Gate Checklist:**
  - [x] ตรวจจับ sample mismatch ระหว่าง count matrix และ metadata ได้ถูกต้อง
  - [x] ตรวจจับ NaN, negative values, duplicate gene IDs ได้แม่นยำ
  - [x] PyDESeq2 size-factor normalization wrapper ทำงานได้และ override พารามิเตอร์ได้
  - [x] Example notebook `notebooks/examples/01_io_and_qc.ipynb` รันผ่านสมบูรณ์

---

### 🎯 Part 2: Differential Expression Integration (FR-3)
- **Deliverables:**
  - `src/btw/de_analysis/{deseq_helper,contrasts,correction}.py`
  - `tests/test_de_analysis.py`
  - `notebooks/examples/02_differential_expression.ipynb`
- **คำสั่งตรวจสอบ (Verification Command):**
  ```bash
  pytest tests/test_de_analysis.py -v --cov=src/btw/de_analysis
  ```
- **Acceptance Gate Checklist:**
  - [x] Wrapper คืนทั้ง `DeseqDataSet`, `DeseqStats` ดั้งเดิม และ standard summary DataFrame (`DEResult`)
  - [x] Multi-contrast batch execution รันหลาย comparison ได้ใน loop เดียว พร้อม Master Table
  - [x] Multiple testing correction ผ่าน `statsmodels` (เช่น Benjamini-Hochberg, Bonferroni) ทำงานถูกต้อง
  - [x] Example notebook `notebooks/examples/02_differential_expression.ipynb` รันผ่านสมบูรณ์

---

### 🎯 Part 3: Publication-Grade Visualization (FR-4)
- **Deliverables:**
  - `src/btw/viz/{style,volcano,pca,ma,dispersion,heatmap,overlaps}.py`
  - `tests/test_viz.py`
  - `notebooks/examples/03_visualization.ipynb`
- **คำสั่งตรวจสอบ (Verification Command):**
  ```bash
  pytest tests/test_viz.py -v --cov=src/btw/viz
  ```
- **Acceptance Gate Checklist:**
  - [x] Volcano plot ใส่ label ยีนสำคัญได้โดยข้อความไม่ทับซ้อน (`adjustText`)
  - [x] PCA plot แสดง variance explained และระบุกลุ่มตัวอย่างตาม metadata ได้ (2D/3D)
  - [x] รองรับทั้ง static (`matplotlib`/`seaborn`) และ interactive (`plotly`)
  - [x] Heatmap ทำ hierarchical clustering แถว/คอลัมน์ได้อย่างถูกต้อง พร้อม metadata color tracks
  - [x] Example notebook `notebooks/examples/03_visualization.ipynb` รันผ่านสมบูรณ์

---

### 🎯 Part 4: Pathway & Functional Enrichment (FR-5)
- **Deliverables:**
  - `src/btw/enrichment/{schema,ora,gsea,activity,viz,cache}.py`
  - `tests/test_enrichment.py`
  - `notebooks/examples/04_pathway_enrichment.ipynb`
- **คำสั่งตรวจสอบ (Verification Command):**
  ```bash
  pytest tests/test_enrichment.py -v --cov=src/btw/enrichment
  ```
- **Acceptance Gate Checklist:**
  - [x] ครบทั้ง 3 แนวทาง: ORA (`enrichr`/`goatools`), GSEA (`prerank`), Activity (`decoupler`)
  - [x] ทั้ง 3 วิธี output ผลลัพธ์ใน schema ตารางมาตรฐานเดียวกัน (`EnrichmentResult`)
  - [x] Unified Dot plot / Bar plot สามารถพล็อตเปรียบเทียบผลจากทั้ง 3 วิธีได้
  - [x] ระบบ Caching ผ่าน `joblib.Memory` ป้องกันการ query ซ้ำ
  - [x] Example notebook `notebooks/examples/04_pathway_enrichment.ipynb` รันผ่านสมบูรณ์

---

### 🎯 Part 5: Advanced Downstream (FR-6, FR-7, FR-8)
- **Deliverables:**
  - `src/btw/annotation/{mapper,gtf}.py`
  - `src/btw/batch_correction/{combat,diagnostic}.py`
  - `src/btw/network/{wgcna_helper,export}.py`
  - `tests/test_advanced.py`
  - `notebooks/examples/05_advanced_downstream.ipynb`
- **คำสั่งตรวจสอบ (Verification Command):**
  ```bash
  pytest tests/test_advanced.py -v
  ```
- **Acceptance Gate Checklist:**
  - [x] Gene ID mapping ↔ symbol ทำงานได้พร้อม local cache
  - [x] ComBat batch correction แสดงผล PCA ก่อนและหลัง correction ได้
  - [x] PyWGCNA module detection แปลงเข้า `networkx` และ export Cytoscape `.sif` ได้
  - [x] Example notebook `notebooks/examples/05_advanced_downstream.ipynb` รันผ่านสมบูรณ์

---

### 🎯 Part 6: Reporting & Acceptance Signoff (FR-9 & Acceptance)
- **Deliverables:**
  - `src/btw/reporting/{report,export}.py`
  - `notebooks/examples/00_end_to_end_pipeline.ipynb`
  - `tests/test_end_to_end.py`
- **คำสั่งตรวจสอบ (Verification Command):**
  ```bash
  pytest tests/ -v --cov=src/btw --cov-report=term-missing
  ```
- **Acceptance Gate Checklist:**
  - [x] Coverage รวมของ `src/btw` ≥ 70% (ทำได้ 82% ทั่วทั้ง 2,131 statements)
  - [x] End-to-End Notebook รันผ่านตั้งแต่ raw counts จนถึงออก report สมบูรณ์ (`notebooks/examples/00_end_to_end_pipeline.ipynb`)
  - [x] Object ดั้งเดิมของ library ต่างๆ ยังคงเข้าถึงได้ตาม Acceptance Criteria (`DeseqDataSet`, `DeseqStats`, `PCA`, `NetworkX`)
  - [x] ผลการวิเคราะห์ reproducible เมื่อรันซ้ำด้วย seed เดียวกัน (ยืนยันผ่าน `test_pipeline_reproducibility`)

---

## 4. AI Handoff Log (บันทึกการส่งต่องานระหว่างเซสชัน)

### 📝 Handoff History

#### [2026-09-09] - Part 7: SRS v2 Rebuild (Hybrid R/Python Architecture) Completed & Verified
- **Milestones Worked On:**
  - **Part 7:** SRS v2 Rebuild ("Hybrid R/Python on Jupyter", FR-10, NFR-10, Acceptance Criteria #6)
- **Completed Deliverables:**
  - `src/btw/r_interop/`: Isolated subpackage connecting BTW to reference R/Bioconductor engines:
    - `bridge.py`: Detection (`is_r_available`, `get_r_version`), package checking (`check_r_package`, `require_r_package`), vectorized conversions (`pandas_to_r_df`, `pandas_to_r_matrix`, `r_to_pandas_df`), synchronized seed (`set_r_seed`), and actionable `MissingRPackageError`.
    - `deseq2_r.py`: Reference R `DESeq2` Wald test and `apeglm`/`ashr` LFC shrinkage + VST transform (`run_r_deseq2`, `run_r_vst`).
    - `limma_r.py`: Reference R `limma-voom` and `limma-trend` via `edgeR::DGEList` and `limma::eBayes` (`run_r_limma`).
    - `combat_r.py`: Reference R `sva::ComBat` and `sva::ComBat_seq` batch correction (`run_r_combat`).
    - `wgcna_r.py`: Reference R `WGCNA::blockwiseModules` co-expression network analysis (`run_r_wgcna`).
    - `clusterprofiler_r.py`: Reference R Bioconductor `clusterProfiler::enrichGO` (`run_r_clusterprofiler`).
    - `__init__.py`: Public export of all R interop functions.
  - **Core Modules Engine Routing:**
    - `src/btw/qc_normalize/normalize.py`: Added `engine="r"|"python"` with automatic Python fallback to `normalize_deseq2` and added `vst_transform()`.
    - `src/btw/de_analysis/deseq_helper.py`: Added `engine`, `method`, `metadata` to `DEResult` container, added engine provenance to `summary()`, and integrated `engine="r"|"python"` dispatch in `run_de()`.
    - `src/btw/de_analysis/contrasts.py`: Multi-contrast batch execution supporting `engine="r"` and `engine="python"`.
    - `src/btw/batch_correction/combat.py`: `run_combat` defaults to R `sva::ComBat` with fallback to `inmoose`.
    - `src/btw/network/wgcna_helper.py`: `detect_coexpression_modules` / `run_wgcna` defaults to R `WGCNA` with fallback to Python co-expression helper.
    - `src/btw/enrichment/ora.py`: Added `run_clusterprofiler` with fallback to Python Enrichr / custom ORA.
    - `src/btw/reporting/report.py`: Injected engine provenance badges into both Markdown reports (`Execution Engine: ...`) and standalone HTML reports (`.engine-badge`).
    - `src/btw/config.py` & `configs/default_config.yaml`: Added `engine`, `method`, and `fallback_to_python` configurations.
    - `src/btw/__init__.py`: Synchronized `set_seed()` across Python and R sessions.
  - **Testing & Quality:**
    - `tests/test_r_interop.py`: 10 comprehensive unit tests for bridge, roundtrip matrix/df conversions, `MissingRPackageError`, and all engine fallbacks.
    - Updated `tests/test_de_analysis.py` and `tests/test_end_to_end.py`.
    - Entire project test suite: **63 passed, 0 failed in 33s**, with **74% total test coverage** across 2,689 statements.
  - **Notebook & Documentation:**
    - Updated `notebooks/examples/00_end_to_end_pipeline.ipynb`: All 21 cells executed cleanly without errors.
    - Updated `README.md` and `pyproject.toml` (`r-stats = ["rpy2>=3.5.0", "tzlocal>=5.0"]`).
- **All Acceptance Criteria Met:**
  1. [x] Reference R packages invoked directly via rpy2 as default engines (`engine="r"`).
  2. [x] Native R objects (`DESeqDataSet`, `MArrayLM`) and Python objects (`DeseqDataSet`, `DeseqStats`, `PCA`, `networkx.Graph`) remain directly accessible without restrictive wrapping.
  3. [x] Pure-Python fallback operates seamlessly when `fallback_to_python=True` without crashes.
  4. [x] Informative, actionable `MissingRPackageError` provided when packages are absent and fallback is disabled.
  5. [x] Engine provenance tracked and reported across all outputs and reports.
  6. [x] 100% tests passing with coverage > 70% (74%).

#### [2026-09-09] - Part 6 Completed and Full Workbench Verified & Signed Off
- **Milestones Worked On:** 
  - **Part 6:** Reporting & Acceptance Signoff (FR-9 & Acceptance Gates)
- **Completed Deliverables:**
  - `src/btw/reporting/report.py`: แก้ไข `sig_df = de_result.get_degs()`, เพิ่ม fallback pure-python `_df_to_markdown` ไม่พึ่งพา `tabulate` ภายนอก
  - `src/btw/reporting/export.py`: ยืนยันการทำงานของ `export_publication_bundle` (Multi-sheet Excel, CSVs, high-res figures, Standalone HTML report with base64 embedded images, Markdown report)
  - `src/btw/io/loader.py`: เพิ่ม `__iter__` method ให้ `BulkDataset` รองรับ tuple unpacking `counts, meta = load_dataset(...)`
  - `src/btw/enrichment/schema.py`: เพิ่ม `significant_terms` helper method ให้กับ `EnrichmentResult`
  - `src/btw/viz/style.py`: ปรับปรุง `save_figure` ให้ unwrap matplotlib tuple `(fig, ax)` และ seaborn `ClusterGrid` ได้อัตโนมัติ
  - `tests/test_reporting.py`: ยูนิตเทสต์รายงานทั้ง Markdown และ HTML และการส่งออก publication bundle ผ่าน 100%
  - `tests/test_end_to_end.py`: ทดสอบบูรณาการ end-to-end ทั้งระบบ (I/O, QC, Normalization, PyDESeq2 DE, Publication Viz, ORA Enrichment, ComBat Batch Correction, Cytoscape Network, Publication Bundle) และทดสอบความเที่ยงตรงซ้ำ (Reproducibility under seed=42)
  - `notebooks/examples/00_end_to_end_pipeline.ipynb`: สมุดงานตัวอย่าง end-to-end ฉบับสมบูรณ์ (17 cells) สาธิตครบทุก FR ตั้งแต่ต้นจนจบ รันผ่าน 100% ปราศจาก cell output noise ใน git
- **Test & Verification Results:**
  - Command: `pytest tests/ -v --cov=src/btw --cov-report=term-missing`
  - Results: **52 passed in 13.66s**
  - Code Coverage: **82%** (2,131 statements, 392 missed) เหนือเกณฑ์ขั้นต่ำ (≥70%)
- **All Acceptance Criteria Met:**
  1. [x] Pipeline รันได้ตั้งแต่ต้นจนจบใน notebook เดียว (`00_end_to_end_pipeline.ipynb`)
  2. [x] Native objects (`DeseqDataSet`, `DeseqStats`, `PCA`, `NetworkX`) ยังคงเปิดให้เข้าถึงได้ตรง
  3. [x] ผลลัพธ์ reproducible 100% (ผ่าน `tests/test_end_to_end.py::test_pipeline_reproducibility`)
  4. [x] Test suite ผ่าน 100% และ coverage ≥ 70% (ได้ 82%)
  5. [x] มี example notebooks ครบทุกโมดูล (00 ถึง 05) พร้อมคำอธิบาย
- **Status:** All parts (Part 0 through Part 6) are 100% complete and fully verified.

#### [2026-09-09 05:18] - Part 3 Completed and Verified
- **Agent/Session ID:** 14e29b1d-042d-4e62-ba2b-e7f922feeba7
- **Milestones Worked On:** 
  - **Part 3:** Publication-Grade Visualization Module (FR-4)
- **Completed Deliverables:**
  - `src/btw/viz/style.py`: Palettes (Nature, Cell, Science), `set_publication_style`, `save_figure` (PNG 300 DPI, SVG, PDF)
  - `src/btw/viz/volcano.py`: `plot_volcano` (adjustText static labels & interactive Plotly mode)
  - `src/btw/viz/pca.py`: `compute_pca` & `plot_pca` (variance explained %, 2D/3D scatter, static & interactive)
  - `src/btw/viz/ma.py`: `plot_ma` (mean expression vs log2FC, static & interactive)
  - `src/btw/viz/dispersion.py`: `plot_dispersion` (genewise, fitted trend, MAP shrinkage)
  - `src/btw/viz/heatmap.py`: `plot_heatmap` (hierarchical clustering, Z-score, metadata color tracks)
  - `src/btw/viz/overlaps.py`: `plot_upset` & `plot_venn` (cross-contrast DEG intersections)
  - `src/btw/viz/__init__.py`: Public visualization API
  - `tests/test_viz.py`: 7 unit tests covering all visualization functions
  - `notebooks/examples/03_visualization.ipynb`: Complete runnable demonstration notebook
- **Test & Verification Results:**
  - Command: `pytest tests/ -v --cov=src/btw --cov-report=term-missing`
  - Results: **31 passed in 14.71s**
  - Code Coverage: **84%** across 1,075 statements
  - Notebook smoke test: All cells executed without error.
- **Known Issues / Blockers:** None.
- **Next Recommended Action:**
  - เริ่มพัฒนา **Part 4: Functional & Pathway Enrichment Module (FR-5)** ซึ่งจะครอบคลุม ORA (`enrichr`/`goatools`), GSEA (`prerank`), และ Activity inference (`decoupler`) พร้อม unified schema และ dotplot/barplot

---

#### [2026-09-09 05:12] - Part 2 Completed and Verified
- **Agent/Session ID:** 14e29b1d-042d-4e62-ba2b-e7f922feeba7
- **Milestones Worked On:** Part 2: Differential Expression (DE) Integration (FR-3)
- **Completed Deliverables:** `deseq_helper.py`, `contrasts.py`, `correction.py`, `tests/test_de_analysis.py`, `02_differential_expression.ipynb`
- **Test & Verification Results:** 24 passed in 9.30s, Coverage 84%.

---

#### [2026-09-09 05:07] - Part 0 & Part 1 Completed and Verified
- **Agent/Session ID:** 14e29b1d-042d-4e62-ba2b-e7f922feeba7
- **Milestones Worked On:** Part 0 (Scaffold) & Part 1 (I/O, Validation, QC)
- **Completed Deliverables:** Configs, `src/btw/{__init__,config,io,qc_normalize}`, `tests/`, `notebooks/examples/01_io_and_qc.ipynb`
- **Test & Verification Results:** 18 passed in 0.34s, Coverage 82%.
