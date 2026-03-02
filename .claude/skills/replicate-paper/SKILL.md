# Skill: /replicate-paper

**Trigger:** `/replicate-paper [paper.md] [data_dir] [code_dir]` or "replicate this paper"

**Purpose:** Full 6-phase autonomous replication of a bulk or single-cell RNA-seq paper. Produces Python and R scripts plus a polished validation report. When the paper introduces a novel computational method or workflow, also produces a detailed plain-language explanation of the method.

---

## Invocation

```
/replicate-paper papers/AuthorYear/original_paper.md data/AuthorYear/
```

With original code:
```
/replicate-paper papers/AuthorYear/original_paper.md data/AuthorYear/ papers/AuthorYear/code/
```

Or with just: "replicate this paper" (Claude will ask for paths if not provided).

**Input format:** Papers must be provided as Markdown (`.md`) files. The user is responsible for converting PDFs to Markdown before invoking this skill (e.g., using docling or other conversion tools). Figures referenced by the paper should be placed in `papers/AuthorYear/original_paper_attachment/`.

---

## The 6-Phase Pipeline

### Phase 1: Intake

**Goal:** Understand exactly what needs to be replicated.

1. Read the paper Markdown file (all sections: Abstract, Methods, Results, Supplementary)
2. Identify **every table and figure** that presents empirical results
3. For each: record the gold standard values, statistics (fold changes, adjusted p-values, cluster counts, trajectory labels, etc.), sample sizes, and source location
4. Save targets to `quality_reports/[paper_name]_replication_targets.md`
5. Summarize: original software, data source (GEO/SRA accession), sample N (cells or samples), key methods, genome/annotation version, any replication package available
6. **If the paper introduces a novel computational method:** note this explicitly and flag Phase 1b below

**Output:** `quality_reports/[paper_name]_replication_targets.md`

---

### Phase 1b: Methods Explanation (for novel computational methods)

**Trigger:** Run this phase if the paper's focus is a novel algorithm, pipeline, or computational workflow.

**Goal:** Produce a clear, detailed, self-contained explanation of the method — understandable to a computational biologist who did not read the paper.

Cover all of the following:

1. **Motivation** — What problem does the method solve? What was missing or insufficient in prior approaches?
2. **High-level intuition** — What is the method doing conceptually, in plain language (one paragraph, no jargon)?
3. **Step-by-step algorithm** — Walk through each computational step in order:
   - Inputs (data format, required preprocessing)
   - Each transformation, model, or inference step
   - Key parameters and their biological meaning
   - Outputs (what is returned and how to interpret it)
4. **Key assumptions** — What does the method assume about the data? When would it fail or give misleading results?
5. **Comparison to alternatives** — How does it differ from the closest competing methods (e.g., Seurat vs. scran normalization; DESeq2 vs. edgeR; UMAP vs. t-SNE; Monocle vs. Slingshot)?
6. **Worked example** — Annotate the paper's Figure 1 / schematic diagram in plain language if one is present
7. **Practical guidance** — What are the most important parameters to tune? What are common failure modes?

Save to: `reports/[paper_name]_methods_explanation.md`

**Output:** `reports/[paper_name]_methods_explanation.md`

---

### Phase 2: Data Audit

**Goal:** Confirm what we can and cannot replicate given the available data.

1. Locate or download the dataset (GEO/SRA accession from paper); document source in script header
2. Compare to paper's described sample:
   - Total N (samples for bulk; cells and samples for single-cell)
   - Key variable distributions (library size, mitochondrial fraction, gene detection rate)
   - Preprocessing steps applied (trimming, alignment tool, quantification method)
3. Apply inclusion/exclusion criteria as stated in Methods (cell QC thresholds, sample filters); document each step's effect on N
4. If variables are missing or differ: document the gap; flag as a known discrepancy
5. Save audit summary to `quality_reports/[paper_name]_data_audit.md`

**Output:** `quality_reports/[paper_name]_data_audit.md`

---

### Phase 3: Code Analysis

**Goal:** Map the paper's methods to our dataset before writing a single line of code.

1. Read original R/Python code from `[code_dir]` if provided (may be individual scripts or a full package structure)
2. Map each variable name in original code → corresponding variable in our dataset
3. Identify methodological steps: preprocessing, normalization, feature selection, model fitting, clustering/DE testing, visualization
4. Flag any steps where original code differs from Methods text (use the paper, not the code, as ground truth)
5. Document the mapping in `quality_reports/[paper_name]_variable_map.md`

**Output:** `quality_reports/[paper_name]_variable_map.md`

---

### Phase 4: Translation

**Goal:** Produce clean, reproducible Python and R scripts that implement the paper's analysis.

**Rules:**
- Line-by-line translation first — **no improvements during replication**
- Follow `python-code-conventions.md` and `r-code-conventions.md` exactly
- Set seed: `random.seed(YYYYMMDD)` + `numpy.random.seed(YYYYMMDD)` (Python); `set.seed(YYYYMMDD)` (R)
- Use `pathlib.Path` for all Python paths; `here::here()` for all R paths
- Comment every non-obvious translation decision
- Refer to `replication-protocol.md` pitfall tables

**Python script:** `replications/[paper_name]/python/replicate.py`

Structure:
```python
# Replication: [Paper Author (Year)]
# Date: YYYY-MM-DD
# Data: [GEO/SRA accession]
# Genome: [assembly + GTF version]
# Python version: X.Y.Z
# Key packages: scanpy X.X, anndata X.X, pydeseq2 X.X

from pathlib import Path
import random
import numpy as np
import pandas as pd

random.seed(YYYYMMDD)
np.random.seed(YYYYMMDD)

DATA_DIR = Path(__file__).parents[3] / "data"
RESULTS_DIR = Path(__file__).parent / "results"
FIGURES_DIR = Path(__file__).parent / "figures"
RESULTS_DIR.mkdir(exist_ok=True)
FIGURES_DIR.mkdir(exist_ok=True)

# --- 1. Load Data ---
# --- 2. Quality Control ---
# --- 3. Normalization / Preprocessing ---
# --- 4. Analysis (DE / Clustering / Trajectory) ---
# --- 5. Save Results ---
```

**R script:** `replications/[paper_name]/R/replicate.R`

Structure:
```r
# Replication: [Paper Author (Year)]
# Date: YYYY-MM-DD
# Data: [GEO/SRA accession]
# Genome: [assembly + GTF version]
# R version: X.Y.Z
# Key packages: DESeq2 X.X, Seurat X.X, edgeR X.X

library(here)
library(tidyverse)
library(DESeq2)   # or Seurat, edgeR, etc.

set.seed(YYYYMMDD)

data_dir <- here("data")
results_dir <- here("replications", "[paper_name]", "R", "results")
figures_dir <- here("replications", "[paper_name]", "R", "figures")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load Data ---
# --- 2. Quality Control ---
# --- 3. Normalization / Preprocessing ---
# --- 4. Analysis (DE / Clustering / Trajectory) ---
# --- 5. Save Results ---
```

**Outputs:**
- `replications/[paper_name]/python/replicate.py`
- `replications/[paper_name]/R/replicate.R`
- `replications/[paper_name]/python/results/` (parquet/h5ad files)
- `replications/[paper_name]/R/results/` (rds files)

---

### Phase 5: Validation

**Goal:** Run both scripts and compare results to gold standard targets.

1. Execute Python script: `python replications/[paper_name]/python/replicate.py`
2. Execute R script: `Rscript replications/[paper_name]/R/replicate.R`
3. Load results; compare to targets using tolerance thresholds from `replication-protocol.md`:
   - Integers (N cells/samples, gene counts): exact
   - Fold changes / log fold changes: ±0.05
   - Adjusted p-values: same significance bracket (< 0.05, < 0.01, < 0.001)
   - Cluster counts / trajectory node counts: exact
   - Cluster proportions / percentages: ±1pp
4. For each mismatch: investigate root cause before proceeding
5. Save `replications/[paper_name]/validation_report.md`

**Output:** `replications/[paper_name]/validation_report.md`

---

### Phase 6: Report

**Goal:** Produce a polished, self-contained replication report.

Report structure:
```markdown
# Replication Report: [Paper Author (Year)]
**Date:** [YYYY-MM-DD]
**Replicator:** Claude (domain-reviewer verified)

## Paper Summary
[1 paragraph: research question, organism/tissue, data type (bulk/scRNA-seq), key finding]

## Methods Summary
[Bullet list: sample, QC thresholds, normalization, key analysis steps, software and versions]

## Data
[Bullet list: accession, N samples/cells after QC, any discrepancies vs. paper sample]

## Methods Explanation
[If paper presents a novel method: link to `reports/[paper_name]_methods_explanation.md`]

## Results Comparison

| Target | Table/Fig | Paper Value | Our Value (Python) | Our Value (R) | Diff | Status |
|--------|-----------|-------------|-------------------|---------------|------|--------|

## Discrepancies
[Each discrepancy: what, investigated how, resolved or not]

## Corrective Steps Taken
[Any adjustments made during validation and why]

## Verdict
**[REPLICATED / PARTIAL / FAILED]**
- Targets matched: N / Total
- Remaining discrepancies: [list or "none"]

## Reproducibility
- Python: X.Y.Z | scanpy X.X | anndata X.X | pydeseq2 X.X
- R: X.Y.Z | DESeq2 X.X | Seurat X.X
- Data: [accession, genome assembly, GTF version]
- Seed: YYYYMMDD
```

Save to: `reports/[paper_name]_replication_report.md`

After saving: run domain-reviewer agent on the report.

---

## Quality Gate

After Phase 6, score the output. Minimum 80/100 to commit.

**Auto-commit if score >= 80:**
```
git add replications/[paper_name]/ reports/[paper_name]_replication_report.md reports/[paper_name]_methods_explanation.md quality_reports/[paper_name]_*.md
git commit -m "Replicate [Paper Author (Year)] -- [VERDICT]: N/Total targets matched"
```

---

## Failure Modes & Recovery

| Failure | Recovery |
|---------|---------|
| Script syntax error | Fix before proceeding |
| N mismatch > 5% | Stop, audit QC thresholds and filtering steps |
| All fold changes off by same factor | Check normalization method (CPM vs. TPM vs. VST vs. log-normalized) |
| DE gene lists differ substantially | Check filtering (min counts, min cells), model formula, reference level |
| Cluster assignments differ | Check resolution parameter, number of PCs, random seed |
| Cannot install Bioconductor package | Document, note in report, use closest alternative |
| Data accession unavailable | Document gap; attempt proxy dataset; flag as ASSUMED in report |
