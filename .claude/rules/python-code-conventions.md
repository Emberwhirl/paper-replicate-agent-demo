---
paths:
  - "replications/**/*.py"
  - "scripts/**/*.py"
---

# Python Scientific Coding Conventions

**Scope:** All Python scripts in `replications/` and `scripts/`.

---

## Paths

- Always use `pathlib.Path` — never string concatenation for paths
- Always use paths relative to the script file or project root:

```python
from pathlib import Path

# Relative to script
DATA_DIR = Path(__file__).parents[3] / "data"

# Or relative to project root (if script is run from root)
DATA_DIR = Path("data")
```

- Never hardcode absolute paths (e.g., `/Users/zhuch/...`) — this breaks on other machines

---

## Reproducibility

At the top of every script that has stochastic elements:

```python
import random
import numpy as np

random.seed(YYYYMMDD)      # Use today's date as integer, e.g., 20260220
np.random.seed(YYYYMMDD)
```

If using other RNG-dependent libraries, seed them explicitly:
```python
import torch
torch.manual_seed(YYYYMMDD)
```

---

## Imports

All imports at the top of the file. Never inside functions. Order:
1. Standard library (`os`, `random`, `pathlib`)
2. Third-party (`numpy`, `pandas`, `scanpy`, `anndata`, `pydeseq2`)
3. Local modules

```python
# Standard library
from pathlib import Path
import random

# Third party
import numpy as np
import pandas as pd
import scanpy as sc        # scRNA-seq
# import pydeseq2          # bulk RNA-seq

# Local
from scripts.utils import load_counts
```

---

## Pandas

- Always specify `dtype` on `pd.read_csv()` for columns that will be used as keys or binary indicators:

```python
df = pd.read_csv(DATA_DIR / "ukb.csv", dtype={"eid": str, "event": int})
```

- Never use `inplace=True` — it causes silent failures and unclear code:

```python
# Bad
df.dropna(inplace=True)

# Good
df = df.dropna()
```

- Use `.copy()` when slicing a DataFrame you will modify:

```python
subset = df[df["age"] >= 40].copy()
```

---

## Modeling

Use established scientific Python libraries, not ad-hoc implementations:

| Task | Library |
|------|---------|
| Single-cell analysis (Python) | `scanpy`, `anndata` |
| Bulk RNA-seq DE (Python) | `pydeseq2` |
| Pseudotime / trajectory | `scFates`, `cellrank` |
| Batch integration | `harmonypy`, `scvi-tools` |
| Doublet detection | `scrublet`, `scDblFinder` (via rpy2 or standalone) |
| Gene set enrichment | `gseapy` |
| Visualization | `matplotlib`, `seaborn`, `scanpy.pl` |

Save all model results as structured files, not just printed output:

```python
# Good — AnnData object (single-cell)
adata.write_h5ad(RESULTS_DIR / "adata_processed.h5ad")

# Good — DE results table
results_df.to_parquet(RESULTS_DIR / "de_results_group1_vs_group2.parquet")
```

---

## Figures

```python
import matplotlib.pyplot as plt
import matplotlib as mpl

# Okabe-Ito colorblind-safe palette
OKABE_ITO = ["#E69F00", "#56B4E9", "#009E73", "#F0E442",
             "#0072B2", "#D55E00", "#CC79A7", "#000000"]

fig, ax = plt.subplots(figsize=(8, 6))
# ... plotting code ...

fig.savefig(FIGURES_DIR / "figure1.png", dpi=300, bbox_inches="tight",
            facecolor="white")
plt.close(fig)
```

- 300 DPI minimum
- White or transparent background (`facecolor="white"` or `"none"`)
- `bbox_inches="tight"` to avoid clipping
- Always close figure after saving to free memory
- Never use `plt.show()` in scripts (breaks headless execution)

---

## Comments

Comment WHY, not WHAT. The code shows what; comments explain non-obvious decisions.

```python
# Bad
# Drop rows where count is zero
df = df[df["count"] > 0]

# Good
# Paper applies a minimum count filter of 10 reads in at least 3 samples (Methods, p. 4)
# This removes lowly-expressed genes before normalization
keep = (counts >= 10).sum(axis=1) >= 3
counts = counts[keep]
```

Always comment translation decisions or parameter choices:

```python
# REPLICATION NOTE: Paper uses VST normalization (DESeq2::varianceStabilizingTransformation)
# not log-CPM. pydeseq2 does not expose VST; using R DESeq2 via subprocess for this step.
```

---

## Script Structure Template

```python
#!/usr/bin/env python3
"""
Replication: [Paper Author (Year)]
Date: YYYY-MM-DD
Data: [GEO/SRA accession]
Genome: [assembly + GTF version]
Python version: 3.X.Y
Key packages: see requirements.txt

Replicates: [Table X, Figure Y]
"""

# ── Imports ────────────────────────────────────────────────────────────────
from pathlib import Path
import random
import numpy as np
import pandas as pd
import scanpy as sc        # for scRNA-seq; replace with pydeseq2 for bulk

# ── Reproducibility ────────────────────────────────────────────────────────
random.seed(20260220)
np.random.seed(20260220)

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).parents[3]
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = Path(__file__).parent / "results"
FIGURES_DIR = Path(__file__).parent / "figures"
RESULTS_DIR.mkdir(exist_ok=True)
FIGURES_DIR.mkdir(exist_ok=True)

# ── 1. Load Data ───────────────────────────────────────────────────────────

# ── 2. Quality Control ─────────────────────────────────────────────────────

# ── 3. Normalization / Preprocessing ──────────────────────────────────────

# ── 4. Analysis (DE / Clustering / Trajectory) ────────────────────────────

# ── 5. Save Results ────────────────────────────────────────────────────────
```
