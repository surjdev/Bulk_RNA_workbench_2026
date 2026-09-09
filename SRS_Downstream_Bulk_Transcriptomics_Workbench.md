# Software Requirements Specification (SRS)
## Downstream Bulk Transcriptomics Workbench (Jupyter + Python)

**เวอร์ชันเอกสาร:** 1.0
**สถานะ:** Draft สำหรับทีม Software Engineer
**ขอบเขต:** เครื่องมือเสริม (utility layer) สำหรับงานวิเคราะห์ downstream ของ Bulk RNA-seq/Transcriptomics บน Jupyter Notebook

---

## 1. Background / วัตถุประสงค์

ทีมมีการวิเคราะห์ downstream ของ Bulk Transcriptomics โดยใช้ PyDESeq2 เป็นเครื่องมือหลักสำหรับ Differential Expression (DE) analysis อยู่แล้ว แต่ยังขาดชุดเครื่องมือ/สาธารณูปโภค (utility) ที่ช่วยให้ workflow ตั้งแต่ normalization, visualization, enrichment analysis, ไปจนถึง reporting ทำงานร่วมกับ numpy/pandas/scikit-learn ได้สะดวกและ reproducible มากขึ้น

**สิ่งที่ต้องการ:** ชุด Python package/utility ภายใน repo เดียว ที่นักวิเคราะห์เรียกใช้ผ่าน Jupyter Notebook ได้โดยตรง

**สิ่งที่ไม่ต้องการ:** GUI, web application, หรือ wrapper ที่ครอบปิดการเข้าถึง library ต้นทาง (PyDESeq2, scikit-learn ฯลฯ) — ผู้ใช้ต้องยังคงเรียกใช้ library เหล่านั้นได้โดยตรงเสมอ utility ที่สร้างขึ้นทำหน้าที่ "ลดงานซ้ำ" (boilerplate) เท่านั้น ไม่ใช่ "แทนที่" การควบคุมของผู้ใช้

---

## 2. Scope & Non-Goals

### In-Scope
- Utility functions/module สำหรับแต่ละขั้นตอนของ downstream analysis (ดูหัวข้อ 3)
- Infrastructure สนับสนุนการทำงานบน Jupyter (reproducibility, caching, config, versioning)
- Documentation และ example notebook สำหรับแต่ละ module

### Non-Goals (Out of Scope)
- ไม่สร้าง GUI/Dashboard/Web app
- ไม่ครอบ (wrap) PyDESeq2 จนผู้ใช้เรียก API เดิมไม่ได้
- ไม่รวม survival/clinical analysis (lifelines) — **ไม่มี clinical data ในสโคปปัจจุบัน**
- ไม่ทำ pipeline orchestration แบบเต็มรูปแบบ (Snakemake/Nextflow) ในเฟสแรก — ใช้ papermill สำหรับ parameterized notebook ก่อน

---

## 3. Functional Requirements (FR)

### FR-1: Data I/O & Validation Module
- โหลด count matrix (CSV/TSV/Excel) และ sample metadata พร้อม schema validation (เช่น ตรวจสอบว่า sample ID ตรงกันระหว่าง count matrix และ metadata)
- แจ้ง warning/error ที่ชัดเจนเมื่อข้อมูลไม่ตรงตามที่คาดหวัง (missing sample, NaN, duplicate gene ID)
- รองรับการ export ผลลัพธ์กลับเป็น CSV/Excel/Parquet

### FR-2: QC & Normalization Utilities
- ฟังก์ชันสรุป QC เบื้องต้น (library size, gene detection rate, distribution plot)
- Wrapper บาง ๆ (thin helper) สำหรับเรียก normalization ของ PyDESeq2 พร้อม default parameter ที่ทีมตกลงกัน แต่ต้องเปิดให้ override ได้ทุก parameter

### FR-3: Differential Expression Integration
- Helper function สำหรับสร้าง `DeseqDataSet` และรัน `DeseqStats` จาก metadata ที่ผ่านการ validate แล้ว (ลด boilerplate ไม่ใช่ครอบ API)
- ฟังก์ชันช่วยจัดการ multiple contrasts (loop รันหลาย comparison แล้วรวมผลเป็นตารางเดียว)
- Multiple testing correction utilities (เชื่อมกับ statsmodels สำหรับกรณีนอกเหนือจาก PyDESeq2)

### FR-4: Visualization Module
- ฟังก์ชันสำเร็จรูป (แต่ปรับแต่งได้) สำหรับ:
  - Volcano plot (พร้อม label gene สำคัญด้วย adjustText)
  - MA plot
  - PCA plot (จาก scikit-learn PCA)
  - Heatmap (top DEGs, clustering ด้วย scipy/seaborn)
  - Dispersion plot (จาก PyDESeq2 output)
- รองรับทั้ง static (matplotlib/seaborn) และ interactive (plotly) mode
- ฟังก์ชันเปรียบเทียบ DEG set หลาย contrast (upsetplot/venn)

### FR-5: Functional / Pathway Enrichment Module

ครอบคลุม 3 แนวทางของ pathway analysis เพื่อรองรับคำถามงานวิจัยที่ต่างกัน:

**FR-5.1: Over-Representation Analysis (ORA)**
- Wrapper บาง ๆ รอบ gseapy (`enrichr` function) สำหรับเชื่อม Enrichr API (KEGG, Reactome, GO, MSigDB Hallmark ฯลฯ)
- Integration กับ goatools สำหรับ GO enrichment แบบละเอียด (control parameter เช่น propagate_counts, background gene set)
- ฟังก์ชันช่วยสร้าง significant gene list จากผล DE โดยอัตโนมัติ พร้อมให้ปรับ threshold ได้ (padj, log2FC cutoff) และให้ผู้ใช้กำหนด background gene set เองได้ (ไม่ default เป็น whole genome เสมอไป)

**FR-5.2: Gene Set Enrichment Analysis (GSEA)**
- Wrapper บาง ๆ รอบ gseapy (`gsea`/`prerank` function)
- ฟังก์ชันแปลงผล DE เป็น ranked gene list โดยอัตโนมัติ (rank ตาม stat หรือ log2FC จาก PyDESeq2 output) โดยไม่ต้องตัด threshold
- รองรับการดึง gene set collection จาก MSigDB (Hallmark, C2, C5 ฯลฯ) ผ่าน gseapy โดยตรง

**FR-5.3: Topology-based / Activity Inference**
- Integration กับ decoupler รองรับหลาย method (ULM, MLM, ORA ภายใน decoupler เอง)
- เชื่อมกับ PROGENy สำหรับ pathway activity score (เช่น MAPK, PI3K, TNFa)
- เชื่อมกับ DoRothEA/CollecTRI สำหรับ Transcription Factor activity inference
- ฟังก์ชันแปลงผล DE (stat/log2FC matrix) ให้อยู่ในรูปแบบ input ที่ decoupler ต้องการโดยอัตโนมัติ

**Cross-cutting requirement:**
- ทั้ง 3 sub-module ต้อง output ผลลัพธ์ในรูปแบบตารางมาตรฐานเดียวกัน (pathway/term name, score, p-value, padj, source method) เพื่อให้เปรียบเทียบผลข้าม method ได้ง่าย
- มีฟังก์ชัน visualization กลางสำหรับแสดงผล enrichment (dot plot / bar plot ranked by significance) ที่ใช้ได้กับผลลัพธ์จากทั้ง 3 แนวทาง

### FR-6: Gene Annotation / ID Mapping Module
- ฟังก์ชัน map gene ID ↔ symbol ↔ annotation ผ่าน mygene หรือ pybiomart
- Caching ผล annotation ไว้ local (ลดการเรียก API ซ้ำ)
- รองรับ parse GTF file (gtfparse) สำหรับกรณีใช้ custom reference

### FR-7: Batch Effect Correction Module
- Integration กับ inmoose/pycombat สำหรับ ComBat batch correction
- ฟังก์ชันเปรียบเทียบ PCA ก่อน/หลัง correction เพื่อ visual QC

### FR-8: Network / Co-expression Module
- Integration กับ PyWGCNA สำหรับ module detection (soft-thresholding, TOM, module-trait correlation)
- ฟังก์ชันช่วยแปลงผล WGCNA module เป็น input สำหรับ networkx เพื่อ visualize/analyze เพิ่มเติม
- ฟังก์ชัน export network เป็นไฟล์ที่ใช้กับ Cytoscape ได้ (เช่น .sif หรือ edge list)

### FR-9: Reporting / Export Module
- รวมผลลัพธ์ (DE table, enrichment result, plot) เป็นรายงานสรุป (HTML ผ่าน nbconvert หรือ markdown)
- Export ตาราง/รูปในรูปแบบพร้อมใช้งานสำหรับ presentation (PNG/SVG resolution สูง, Excel)

---

## 4. Non-Functional Requirements (NFR)

| หมวด | ข้อกำหนด |
|---|---|
| **Reproducibility** | ใช้ conda/mamba (`environment.yml`) หรือ poetry (`pyproject.toml`) pin เวอร์ชันทุก dependency; random seed ต้อง set ได้จากส่วนกลาง |
| **Performance** | ฟังก์ชันที่คำนวณหนัก (WGCNA, PCA, enrichment) ต้องรองรับ caching ผ่าน `joblib.Memory` เพื่อไม่ให้รันซ้ำโดยไม่จำเป็น |
| **Compatibility** | รองรับ Python 3.10+ และรันได้ปกติใน JupyterLab (ไม่ใช้ library ที่ conflict กับ ipykernel) |
| **Extensibility** | ทุก module ต้องออกแบบให้ override/substitute เครื่องมือได้ (เช่น สลับจาก gseapy เป็น goatools ได้โดยไม่กระทบ pipeline อื่น) |
| **Testability** | มี unit test (pytest) ครอบคลุม utility function ที่เขียนเอง (ไม่ต้อง test ตัว library ภายนอกซ้ำ) coverage ขั้นต่ำ 70% สำหรับโค้ดใน `src/` |
| **Code Quality** | ใช้ black, ruff (หรือ flake8), isort ผ่าน pre-commit hook |
| **Documentation** | ทุก module มี docstring (Google/NumPy style) และมี example notebook อย่างน้อย 1 ไฟล์ต่อ module |
| **Logging** | ใช้ loguru สำหรับ log ขั้นตอนสำคัญ (เช่น จำนวน gene ที่ผ่าน filter, parameter ที่ใช้รัน DE) |

---

## 5. Technical Stack / Dependencies

### Core Analysis
`pandas`, `numpy`, `scipy`, `pydeseq2`, `statsmodels`, `scikit-learn`

### Visualization
`matplotlib`, `seaborn`, `plotly`, `adjustText`, `upsetplot`

### Enrichment / Pathway
`gseapy`, `goatools`, `decoupler`

### Annotation
`mygene`, `pybiomart`, `gtfparse`

### Batch Correction
`inmoose` (หรือ `pycombat`)

### Network / Co-expression
`PyWGCNA`, `networkx`

### Workbench Infrastructure
`jupytext`, `papermill`, `ipywidgets`, `joblib`, `dvc`, `pytest`, `black`, `ruff`, `isort`, `pre-commit`, `loguru`

> **หมายเหตุ:** ตัดกลุ่ม survival analysis (`lifelines`) ออกจากสโคปนี้ เนื่องจากปัจจุบันไม่มี clinical/survival data ประกอบการวิเคราะห์ หากมีในอนาคตให้เพิ่มเป็น module แยก (FR-10)

---

## 6. Proposed Project Structure

```
bulk-transcriptomics-workbench/
├── environment.yml
├── pyproject.toml
├── README.md
├── src/
│   └── btw/                      # package หลัก
│       ├── io/                   # FR-1
│       ├── qc_normalize/         # FR-2
│       ├── de_analysis/          # FR-3
│       ├── viz/                  # FR-4
│       ├── enrichment/           # FR-5
│       ├── annotation/           # FR-6
│       ├── batch_correction/     # FR-7
│       ├── network/              # FR-8
│       └── reporting/            # FR-9
├── notebooks/
│   └── examples/                 # example notebook ต่อ module
├── tests/
├── configs/                      # yaml config ต่อ project/analysis
└── .pre-commit-config.yaml
```

---

## 7. Acceptance Criteria

1. นักวิเคราะห์สามารถรัน pipeline ตั้งแต่ raw count matrix → DE result → enrichment result → report ได้ภายใน notebook เดียว โดยเรียก utility function จาก `btw` package
2. ทุก utility function ยังคงเปิดให้เข้าถึง object ดั้งเดิมของ PyDESeq2/scikit-learn ได้ (ไม่ถูกซ่อนไว้)
3. ผลลัพธ์ (DE table, plot, enrichment) reproducible เมื่อรันซ้ำด้วย config เดียวกัน
4. Unit test ผ่านทั้งหมดใน CI และ coverage ≥ 70%
5. มี example notebook ครบทุก module พร้อมคำอธิบาย

---

## 8. Out of Scope (เพื่อความชัดเจน)

- Survival/clinical analysis (lifelines) — ไม่อยู่ในสโคปเนื่องจากไม่มี clinical data
- Web dashboard / GUI
- Full workflow orchestration (Snakemake/Nextflow) ในเฟสแรก
- Single-cell analysis (scanpy) — นอกสโคป bulk transcriptomics
