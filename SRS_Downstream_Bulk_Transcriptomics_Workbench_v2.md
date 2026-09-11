# Software Requirements Specification (SRS)
## Downstream Bulk Transcriptomics Workbench (Hybrid R/Python on Jupyter)

**เวอร์ชันเอกสาร:** 2.0 — ปรับเป็น Hybrid R/Python บน Jupyter Notebook เดียว
**สถานะ:** Draft สำหรับทีม Software Engineer
**ขอบเขต:** เครื่องมือเสริม (utility layer) สำหรับงานวิเคราะห์ downstream ของ Bulk RNA-seq/Transcriptomics บน Jupyter Notebook

---

## 0. หลักการเลือกภาษา (Language Selection Principle)

**กฎหลัก:** ถ้ามี R/Bioconductor package ที่เป็น reference implementation ของ method นั้น ต้องเรียก R ตัวนั้นผ่าน `rpy2` **ห้าม reimplement หรือใช้ Python port แทน** ส่วนงานที่ไม่มี R reference ชัดเจน หรือเป็น infra/glue/visualization ให้ทำใน Python

| หมวดงาน | ภาษา | เหตุผล |
|---|---|---|
| Differential Expression core (DESeq2) | **R** | `DESeq2` คือต้นฉบับ; `PyDESeq2` เป็น port — ให้เปลี่ยนมาเรียก R โดยตรง |
| limma-trend / limma-voom | **R** | ไม่มี Python เทียบเท่า |
| Batch correction (ComBat) | **R** | `sva::ComBat` คือต้นฉบับ; `inmoose`/`pycombat` เป็น port |
| Co-expression network (WGCNA) | **R** | `WGCNA` คือต้นฉบับ; `PyWGCNA` เป็น port |
| GO enrichment แบบละเอียด (topology-aware) | **R** | `clusterProfiler`/`topGO` เป็น field-standard ที่ community อ้างอิงมากกว่า goatools ในหลายกรณี — ให้เป็น option เสริมควบคู่กับ goatools เดิม |
| ORA/GSEA ผ่าน Enrichr/MSigDB | **Python** | `gseapy` เป็น thin wrapper รอบ REST API อยู่แล้ว ไม่มีประเด็น "ต้นฉบับ" ต่างภาษา |
| Topology-based activity inference (decoupler) | **Python** | `decoupler` (Python, Saezlab) ถูก maintain คู่ขนานกับ `decoupleR` (R) อย่างสม่ำเสมอ ไม่มีเหตุผลต้องย้าย |
| I/O, validation, QC summary | **Python** | glue logic ธรรมดา |
| Visualization | **Python** | matplotlib/seaborn/plotly เป็น presentation layer |
| Gene annotation / ID mapping | **Python** | mygene/pybiomart ใช้งานสะดวกจาก Python อยู่แล้ว |
| Reporting/export | **Python** | nbconvert/papermill |

โครงสร้างยังคงเป็น **Python package (`btw`) บน Jupyter** เหมือนเดิม ส่วนที่เป็น R จะถูกเรียกผ่าน submodule ที่ isolate ชัดเจน (ดูข้อ 6)

---

## 1. Background / วัตถุประสงค์

ทีมมีการวิเคราะห์ downstream ของ Bulk Transcriptomics โดยใช้ PyDESeq2 เป็นเครื่องมือหลักสำหรับ Differential Expression (DE) analysis อยู่แล้ว แต่ยังขาดชุดเครื่องมือ/สาธารณูปโภค (utility) ที่ช่วยให้ workflow ตั้งแต่ normalization, visualization, enrichment analysis, ไปจนถึง reporting ทำงานร่วมกันได้สะดวกและ reproducible มากขึ้น

**การเปลี่ยนแปลงสำคัญในเวอร์ชันนี้:** แทนที่จะยึด PyDESeq2 เป็น engine เดียวสำหรับ DE, workbench จะเรียก **R `DESeq2`/`limma`/`edgeR` ต้นฉบับผ่าน `rpy2`** เป็นค่าเริ่มต้น เพื่อความถูกต้องทางสถิติที่ตรงกับ literature ที่ทีมอ้างอิง ส่วน PyDESeq2 ยังคงเปิดให้เลือกใช้ได้แบบ explicit สำหรับกรณีที่ไม่ต้องการพึ่ง R เลย

**สิ่งที่ต้องการ:** ชุด Python package/utility ภายใน repo เดียว ที่นักวิเคราะห์เรียกใช้ผ่าน Jupyter Notebook ได้โดยตรง โดยส่วนที่เรียก R จะถูกซ่อนอยู่หลัง thin wrapper — ผู้ใช้ไม่ต้องเขียน R เอง

**สิ่งที่ไม่ต้องการ:** GUI, web application, หรือ wrapper ที่ครอบปิดการเข้าถึง library ต้นทาง (DESeq2, PyDESeq2, scikit-learn ฯลฯ) — ผู้ใช้ต้องยังคงเรียกใช้ library เหล่านั้นได้โดยตรงเสมอ

---

## 2. Scope & Non-Goals

### In-Scope
- Utility functions/module สำหรับแต่ละขั้นตอนของ downstream analysis (ดูหัวข้อ 3)
- R interop layer สำหรับ module ที่ reference implementation เป็น R (DE, batch correction, co-expression network)
- Infrastructure สนับสนุนการทำงานบน Jupyter (reproducibility, caching, config, versioning) — รวม environment เดียวที่ pin ทั้ง Python และ R/Bioconductor
- Documentation และ example notebook สำหรับแต่ละ module รวมถึงตัวอย่างการใช้ R-backed function

### Non-Goals (Out of Scope)
- ไม่สร้าง GUI/Dashboard/Web app
- ไม่ครอบ (wrap) DESeq2/PyDESeq2/scikit-learn จนผู้ใช้เรียก API เดิมไม่ได้
- ไม่รวม survival/clinical analysis (lifelines) — ไม่มี clinical data ในสโคปปัจจุบัน
- ไม่ทำ pipeline orchestration แบบเต็มรูปแบบ (Snakemake/Nextflow) ในเฟสแรก — ใช้ papermill สำหรับ parameterized notebook ก่อน
- ไม่ reimplement สถิติของ R package ใด ๆ ใน Python

---

## 3. Functional Requirements (FR)

### FR-1: Data I/O & Validation Module — **ภาษา: Python**
- โหลด count matrix (CSV/TSV/Excel) และ sample metadata พร้อม schema validation (เช่น ตรวจสอบว่า sample ID ตรงกันระหว่าง count matrix และ metadata)
- แจ้ง warning/error ที่ชัดเจนเมื่อข้อมูลไม่ตรงตามที่คาดหวัง (missing sample, NaN, duplicate gene ID)
- รองรับการ export ผลลัพธ์กลับเป็น CSV/Excel/Parquet

### FR-2: QC & Normalization Utilities — **ภาษา: Python (glue) + R (engine)**
- ฟังก์ชันสรุป QC เบื้องต้น (library size, gene detection rate, distribution plot) → Python
- Normalization: เรียก `DESeq2::estimateSizeFactors`/`vst` ผ่าน rpy2 เป็นค่าเริ่มต้น พร้อม default parameter ที่ทีมตกลงกัน แต่ต้องเปิดให้ override ได้ทุก parameter
- ต้องรองรับ fallback ไป PyDESeq2 normalization ได้ผ่าน config flag เดียว (ไม่ต้องแก้โค้ด notebook)

### FR-3: Differential Expression Integration — **ภาษา: R (ค่าเริ่มต้น) / Python (ทางเลือก)**
- Helper function สำหรับสร้าง R `DESeqDataSet` และรัน `DESeq`/`results` ผ่าน `rpy2` จาก metadata ที่ผ่านการ validate แล้ว (ลด boilerplate ไม่ใช่ครอบ API) — **นี่คือ engine เริ่มต้น**
- ฟังก์ชันช่วยจัดการ multiple contrasts (loop รันหลาย comparison แล้วรวมผลเป็นตารางเดียว) — glue logic ทำใน Python, เรียก R ต่อ contrast
- ทางเลือกเสริม: `PyDESeq2` (Python-native) สำหรับทีมที่ต้องการหลีกเลี่ยง R dependency ทั้งหมด — ต้อง log ให้ชัดว่า run ไหนใช้ engine ใด เพื่อไม่ให้ผลจากสอง engine ถูกเทียบกันโดยไม่รู้ตัว
- limma-trend / limma-voom → R `limma` ผ่าน rpy2 (ไม่มี Python เทียบเท่า)
- Multiple testing correction utilities (เชื่อมกับ statsmodels สำหรับกรณีนอกเหนือจาก DESeq2/limma) → Python

### FR-4: Visualization Module — **ภาษา: Python**
- ฟังก์ชันสำเร็จรูป (แต่ปรับแต่งได้) สำหรับ:
  - Volcano plot (พร้อม label gene สำคัญด้วย adjustText)
  - MA plot
  - PCA plot (จาก scikit-learn PCA, ทำงานบน VST-transformed data จาก DESeq2/limma)
  - Heatmap (top DEGs, clustering ด้วย scipy/seaborn)
  - Dispersion plot (จาก DESeq2 output ที่แปลงกลับมาเป็น DataFrame)
- รองรับทั้ง static (matplotlib/seaborn) และ interactive (plotly) mode
- ฟังก์ชันเปรียบเทียบ DEG set หลาย contrast (upsetplot/venn)

### FR-5: Functional / Pathway Enrichment Module

**FR-5.1: Over-Representation Analysis (ORA) — ภาษา: Python**
- Wrapper บาง ๆ รอบ gseapy (`enrichr`) สำหรับเชื่อม Enrichr API (KEGG, Reactome, GO, MSigDB Hallmark ฯลฯ)
- Integration กับ goatools สำหรับ GO enrichment แบบละเอียด
- **เสริมใหม่:** option เรียก R `clusterProfiler` ผ่าน rpy2 สำหรับกรณีที่ต้องการผลตรงกับ paper ที่ใช้ clusterProfiler อ้างอิง (field-standard สำหรับ GO/KEGG ในหลาย journal)
- ฟังก์ชันช่วยสร้าง significant gene list จากผล DE โดยอัตโนมัติ พร้อมให้ปรับ threshold ได้ (padj, log2FC cutoff) และให้ผู้ใช้กำหนด background gene set เองได้

**FR-5.2: Gene Set Enrichment Analysis (GSEA) — ภาษา: Python**
- Wrapper บาง ๆ รอบ gseapy (`gsea`/`prerank`)
- ฟังก์ชันแปลงผล DE เป็น ranked gene list โดยอัตโนมัติ (rank ตาม stat จาก DESeq2/limma output)
- รองรับการดึง gene set collection จาก MSigDB ผ่าน gseapy โดยตรง

**FR-5.3: Topology-based / Activity Inference — ภาษา: Python**
- Integration กับ `decoupler` รองรับหลาย method (ULM, MLM, ORA ภายใน decoupler เอง)
- เชื่อมกับ PROGENy, DoRothEA/CollecTRI ผ่าน decoupler โดยตรง
- ฟังก์ชันแปลงผล DE (stat/log2FC matrix) ให้อยู่ในรูปแบบ input ที่ decoupler ต้องการโดยอัตโนมัติ

**Cross-cutting requirement:**
- ทั้ง 3 sub-module ต้อง output ผลลัพธ์ในรูปแบบตารางมาตรฐานเดียวกัน (pathway/term name, score, p-value, padj, source method) เพื่อให้เปรียบเทียบผลข้าม method ได้ง่าย รวมถึงผลจาก R clusterProfiler ต้อง normalize เข้า schema เดียวกันนี้ด้วย
- มีฟังก์ชัน visualization กลาง (Python) สำหรับแสดงผล enrichment ที่ใช้ได้กับผลลัพธ์จากทุกแนวทาง

### FR-6: Gene Annotation / ID Mapping Module — **ภาษา: Python**
- ฟังก์ชัน map gene ID ↔ symbol ↔ annotation ผ่าน mygene หรือ pybiomart
- Caching ผล annotation ไว้ local (ลดการเรียก API ซ้ำ)
- รองรับ parse GTF file (gtfparse) สำหรับกรณีใช้ custom reference

### FR-7: Batch Effect Correction Module — **ภาษา: R**
- เรียก R `sva::ComBat` ผ่าน rpy2 เป็น engine หลัก (ต้นฉบับของ method)
- ฟังก์ชันเปรียบเทียบ PCA ก่อน/หลัง correction เพื่อ visual QC → Python (plot จาก DataFrame ที่แปลงกลับมาแล้ว)

### FR-8: Network / Co-expression Module — **ภาษา: R (engine) + Python (glue/export)**
- เรียก R `WGCNA` ผ่าน rpy2 สำหรับ module detection (soft-thresholding, TOM, module-trait correlation) — ต้นฉบับของ method
- ฟังก์ชันช่วยแปลงผล WGCNA module (แปลงกลับเป็น DataFrame แล้ว) เป็น input สำหรับ `networkx` (Python) เพื่อ visualize/analyze เพิ่มเติม
- ฟังก์ชัน export network เป็นไฟล์ที่ใช้กับ Cytoscape ได้ (เช่น .sif หรือ edge list) → Python

### FR-9: Reporting / Export Module — **ภาษา: Python**
- รวมผลลัพธ์ (DE table, enrichment result, plot) เป็นรายงานสรุป (HTML ผ่าน nbconvert หรือ markdown)
- Export ตาราง/รูปในรูปแบบพร้อมใช้งานสำหรับ presentation (PNG/SVG resolution สูง, Excel)
- ระบุใน report ทุกครั้งว่าแต่ละผลลัพธ์มาจาก engine ใด (R หรือ Python) เพื่อความโปร่งใสเวลาตีพิมพ์/ตรวจสอบย้อนหลัง

### FR-10: R↔Python Interop Layer (ใหม่)
- ทุกฟังก์ชันที่เรียก R ต้อง: (1) รับ input เป็น `pandas.DataFrame`, (2) แปลงเป็น R object ภายในฟังก์ชันผ่าน `rpy2`, (3) รัน R function, (4) แปลงผลกลับเป็น `pandas.DataFrame` ก่อน return — ผู้ใช้ไม่ต้องยุ่งกับ R syntax เลย
- ต้อง handle กรณี R package ไม่ได้ติดตั้งด้วย error message ที่ชัดเจน พร้อมคำแนะนำวิธีติดตั้ง (ไม่ใช่ traceback ดิบจาก rpy2)
- R dependency ต้อง isolate เป็น optional extra (`pip install btw[r-stats]`) ไม่ให้ทั้ง package ล้มถ้าไม่มี R environment
- โมดูล R ทั้งหมดอยู่ใต้ `btw/r_interop/` แยกจาก pure-Python modules อย่างชัดเจน เพื่อให้ตรวจสอบ dependency boundary ได้ง่าย

---

## 4. Non-Functional Requirements (NFR)

| หมวด | ข้อกำหนด |
|---|---|
| **Reproducibility** | ใช้ conda/mamba (`environment.yml`) เดียว pin เวอร์ชันทุก dependency ทั้ง Python และ R/Bioconductor (`DESeq2`, `limma`, `sva`, `WGCNA`, `clusterProfiler`); random seed ต้อง set ได้จากส่วนกลางทั้งสองฝั่งภาษา |
| **Performance** | ฟังก์ชันที่คำนวณหนัก (WGCNA, PCA, enrichment) ต้องรองรับ caching ผ่าน `joblib.Memory`; การเรียก rpy2 ต้อง batch ข้อมูลแทนการเรียกทีละแถวเพื่อลด overhead ของการแปลง object ข้ามภาษา |
| **Compatibility** | รองรับ Python 3.10+ และ R 4.x, รันได้ปกติใน JupyterLab (ไม่ใช้ library ที่ conflict กับ ipykernel) |
| **Extensibility** | ทุก module ต้องออกแบบให้ override/substitute เครื่องมือได้ (เช่น สลับจาก DESeq2(R) เป็น PyDESeq2 ได้ผ่าน config โดยไม่กระทบ pipeline อื่น) |
| **Testability** | มี unit test (pytest) ครอบคลุม utility function ที่เขียนเอง + regression test เทียบผล rpy2 wrapper กับรัน R script ตรง ๆ coverage ขั้นต่ำ 70% สำหรับโค้ดใน `src/` |
| **Code Quality** | ใช้ black, ruff (หรือ flake8), isort ผ่าน pre-commit hook (ฝั่ง Python); R glue script ที่มีอยู่ (ถ้ามี) ควรสั้นและเป็น function call ตรงไปตรงมา ไม่มี control flow ซับซ้อน |
| **Documentation** | ทุก module มี docstring (Google/NumPy style) และมี example notebook อย่างน้อย 1 ไฟล์ต่อ module รวมถึงระบุชัดว่า module ไหนเรียก R engine |
| **Logging** | ใช้ loguru สำหรับ log ขั้นตอนสำคัญ (เช่น จำนวน gene ที่ผ่าน filter, parameter ที่ใช้รัน DE, **engine ที่ใช้ (R/Python) ต่อการรันแต่ละครั้ง**) |
| **R dependency isolation** | R/Bioconductor ต้องเป็น optional extra แยกจาก core install เสมอ (ดู FR-10) |

---

## 5. Technical Stack / Dependencies

### 5.1 Python (core)
`pandas`, `numpy`, `scipy`, `statsmodels`, `scikit-learn`, `rpy2`

### 5.2 Python (optional DE, explicit alternative to R engine)
`pydeseq2`

### 5.3 Visualization
`matplotlib`, `seaborn`, `plotly`, `adjustText`, `upsetplot`

### 5.4 Enrichment / Pathway (Python)
`gseapy`, `goatools`, `decoupler`

### 5.5 Annotation (Python)
`mygene`, `pybiomart`, `gtfparse`

### 5.6 R / Bioconductor (ผ่าน rpy2, optional extra `[r-stats]`)
```
DESeq2         # DE analysis — engine เริ่มต้น
limma          # limma-trend / limma-voom
edgeR          # ทางเลือกเสริมสำหรับ count-based DE
sva            # ComBat batch correction
WGCNA          # co-expression network
clusterProfiler # GO/KEGG enrichment (ทางเลือกเสริมของ FR-5.1)
```

### 5.7 Workbench Infrastructure
`jupytext`, `papermill`, `ipywidgets`, `joblib`, `dvc`, `pytest`, `black`, `ruff`, `isort`, `pre-commit`, `loguru`

> **หมายเหตุ:** ตัดกลุ่ม survival analysis (`lifelines`) ออกจากสโคปนี้ เนื่องจากปัจจุบันไม่มี clinical/survival data ประกอบการวิเคราะห์ หากมีในอนาคตให้เพิ่มเป็น module แยก (FR-11 เดิม, เลื่อนเป็น FR-12 ในเวอร์ชันนี้)

---

## 6. Proposed Project Structure

```
bulk-transcriptomics-workbench/
├── environment.yml                # pin ทั้ง Python + R/Bioconductor
├── pyproject.toml
├── README.md
├── src/
│   └── btw/
│       ├── io/                    # FR-1 (Python)
│       ├── qc_normalize/          # FR-2 (Python glue + R engine)
│       ├── de_analysis/           # FR-3 (R engine เริ่มต้น, PyDESeq2 ทางเลือก)
│       ├── viz/                   # FR-4 (Python)
│       ├── enrichment/            # FR-5 (Python, + R clusterProfiler ทางเลือก)
│       ├── annotation/            # FR-6 (Python)
│       ├── batch_correction/      # FR-7 (R engine)
│       ├── network/               # FR-8 (R engine + Python glue)
│       ├── reporting/             # FR-9 (Python)
│       └── r_interop/             # FR-10: rpy2 wrapper กลาง, ใช้ร่วมกันทุก module ข้างบน
├── notebooks/
│   └── examples/                  # example notebook ต่อ module รวมตัวอย่าง R-backed path
├── tests/
│   ├── unit/                      # Python unit test
│   └── r_regression/              # เทียบผล rpy2 wrapper กับรัน R script ตรง ๆ
├── configs/                       # yaml config ต่อ project/analysis (รวม engine selection flag)
└── .pre-commit-config.yaml
```

---

## 7. Acceptance Criteria

1. นักวิเคราะห์สามารถรัน pipeline ตั้งแต่ raw count matrix → DE result (ผ่าน R DESeq2) → enrichment result → report ได้ภายใน notebook เดียว โดยเรียก utility function จาก `btw` package โดยไม่ต้องเขียน R เอง
2. ทุก utility function ยังคงเปิดให้เข้าถึง object ดั้งเดิมของ DESeq2 (ผ่าน rpy2 handle)/PyDESeq2/scikit-learn ได้ (ไม่ถูกซ่อนไว้)
3. ผลลัพธ์ (DE table, plot, enrichment) reproducible เมื่อรันซ้ำด้วย config เดียวกัน ไม่ว่าจะเลือก engine ใด
4. Unit test + R regression test ผ่านทั้งหมดใน CI และ coverage ≥ 70%
5. มี example notebook ครบทุก module พร้อมคำอธิบาย รวมถึงตัวอย่างสลับ engine ระหว่าง R DESeq2 กับ PyDESeq2
6. Core Python-only workflow (I/O, viz, gseapy-based enrichment) ยังติดตั้งและรันได้แม้ไม่มี R environment (R extra optional จริงตาม FR-10)

---

## 8. Out of Scope (เพื่อความชัดเจน)

- Survival/clinical analysis (lifelines) — ไม่อยู่ในสโคปเนื่องจากไม่มี clinical data
- Web dashboard / GUI
- Full workflow orchestration (Snakemake/Nextflow) ในเฟสแรก
- Single-cell analysis (scanpy) — นอกสโคป bulk transcriptomics (ดู scRNA-seq workbench แยกต่างหาก)
- Reimplement สถิติของ R package ใด ๆ ใน Python
