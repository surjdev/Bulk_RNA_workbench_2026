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
| **Part 6** | Reporting & End-to-End Workbench | FR-9, Acceptance | ⏳ Pending | - | Ready to Start |

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
  - [ ] Coverage รวมของ `src/btw` ≥ 70%
  - [ ] End-to-End Notebook รันผ่านตั้งแต่ raw counts จนถึงออก report สมบูรณ์
  - [ ] Object ดั้งเดิมของ library ต่างๆ ยังคงเข้าถึงได้ตาม Acceptance Criteria

---

## 4. AI Handoff Log (บันทึกการส่งต่องานระหว่างเซสชัน)

### 📝 Handoff History

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
