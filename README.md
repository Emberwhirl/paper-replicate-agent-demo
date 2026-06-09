# RNA-seq Empirical Replication Agent

A structured **Claude Code** workflow for empirically replicating published **bulk and single-cell RNA-seq** studies.

You describe a paper. Claude reads it, identifies every table and figure to reproduce, plans the analysis, writes the R/Python scripts, runs them, checks the outputs against the published numbers within explicit tolerances, documents every discrepancy, and hands back a validation report — like a research contractor who owns the whole pipeline. When a paper introduces a novel computational method, Claude also writes a plain-language explanation of how it works.

The goal is simple: take the friction out of reproducing a paper, so your time goes to the science instead of the plumbing.

> Forked from and built on the excellent UK Biobank replication workflow by [**maxwell2732 (朱晨 Chen Zhu | 遗传社科研究)**](https://github.com/maxwell2732/), then re-targeted for RNA-seq and extended with a suite of single-cell skills built on the [scverse](https://scverse.org/) ecosystem. Heartfelt thanks to the original author for sharing such a thoughtful framework.

---

## What you get

- **A repeatable replication pipeline** — 6 phases from paper intake to a polished report, with quality gates and tolerance checks built in.
- **Bulk RNA-seq support** — DESeq2 / edgeR / limma in R, or pydeseq2 in Python.
- **Single-cell RNA-seq support** — a dedicated set of [scverse](https://scverse.org/)-based skills covering QC, normalization, clustering, batch integration, and deep-learning models (scVI/scANVI/totalVI/…).
- **Guardrails** — plan-first workflow, automated code/manuscript review agents, and an 80/90/95 quality scoring gate so nothing ships half-finished.

---

## Quick Start

**1. Clone and open**

```bash
git clone https://github.com/Emberwhirl/paper-replicate-agent-demo.git
cd paper-replicate-agent-demo
claude
```

**2. Drop your paper in** `papers/[PaperName]/original_paper.md` and your data in `data/`.

**3. Describe the task.** For example:

> I want to replicate **Author (Year)**. The paper is in `papers/Author2024/original_paper.md` and the RNA-seq data is in `data/Author2024/`. Enter plan mode, read the paper, identify all empirical targets, and plan the replication.

Claude reads the paper and config, inventories the data, lists every target table/figure, enters **plan mode**, drafts a step-by-step plan, and waits for your approval. After you approve, it implements: runs scripts, checks outputs against tolerance thresholds, scores the result, and writes a validation report.

**Single-cell?** You can also invoke a specific skill directly, e.g.:

> `/scanpy-basics data/pbmc/` — run the standard scRNA-seq pipeline (QC → clustering → DE).
>
> `/scvi-tools data/integrated.h5ad --model scvi` — deep-learning batch integration.

---

## How It Works

### Contractor mode

You describe a task; Claude does the work and reports back. For complex or ambiguous requests it first writes a short **requirements spec** (MUST / SHOULD / MAY) for you to approve. Then it plans → implements → verifies → reviews → scores → summarizes. Say **"just do it"** and it auto-commits once the work clears the quality gate.

### The replication pipeline (6 phases)

| Phase | What happens |
|-------|-------------|
| **0. Paper intake** | Read the paper; record every empirical target in `quality_reports/[paper]_replication_targets.md` |
| **1. Inventory & data audit** | Confirm the data matches the paper's described sample (N, groups, key variables) |
| **2. Translate & execute** | Write R/Python that matches the original method specification exactly |
| **3. Verify match** | Compare outputs to targets within tolerance thresholds |
| **4. Document discrepancies** | Investigate every near-miss and explain it |
| **5. Report** | Save `replications/[paper]/validation_report.md` and a polished `reports/[paper]_replication_report.md` |

### Quality gates

Every script and report is scored 0–100. Scores below threshold block the action:

| Score | Gate |
|-------|------|
| **80** | Commit |
| **90** | Pull request |
| **95** | Excellence (aspirational) |

### Tolerance thresholds

| Quantity | Tolerance |
|----------|-----------|
| Sample sizes (N, events) | Exact match |
| Point estimates (HR, OR, β, log2FC) | ±0.01 |
| Standard errors | ±0.05 |
| P-values | Same significance bracket |
| Percentages | ±0.1 pp |

---

## Skills

The single-cell skills are the main RNA-seq addition to the original framework; the rest are the general replication and review tooling.

| Skill | What it does |
|-------|-------------|
| `/replicate-paper` | Full 6-phase replication pipeline (Markdown input; explains novel methods) |
| `/data-analysis` | End-to-end R analysis with publication-ready output |
| `/review-r` | R code quality review |
| `/review-paper` | Manuscript review: structure, methods, referee objections |
| `/proofread` | Grammar / typo / consistency review of reports |
| `/devils-advocate` | Challenge design decisions before committing |
| `/commit` | Stage, commit, open a PR, and merge to main |
| `/scanpy-basics` | **Default scRNA-seq pipeline**: QC, normalization, PCA/UMAP/t-SNE, clustering, DE, visualization |
| `/single-cell-rna-qc` | Specialized MAD-based QC with detailed diagnostics |
| `/single-cell-integration-ingest-bbknn` | Batch integration via `ingest` (label transfer) or BBKNN (symmetric) |
| `/scvi-tools` | Deep-learning models: scVI, scANVI, totalVI, PeakVI, MultiVI, DestVI, veloVI |
| `/kallisto-bustools-import` | Import kb-python output and filter empty droplets (knee/inflection) |

<details>
<summary><strong>Also included: 4 agents, 15 rules, 6 hooks, 7 templates</strong> (click to expand)</summary>

### Agents (`.claude/agents/`)

| Agent | What it does |
|-------|-------------|
| `domain-reviewer` | Senior bioinformatics referee (Nature Methods / Genome Biology standard): normalization, batch effects, statistical models, code–method alignment |
| `r-reviewer` | R code quality, reproducibility, RNA-seq correctness |
| `proofreader` | Grammar, typos, consistency in reports |
| `verifier` | End-to-end task completion verification |

### Rules (`.claude/rules/`)

**Always-on** (load every session):

| Rule | What it enforces |
|------|-----------------|
| `plan-first-workflow` | Plan mode for non-trivial tasks + context preservation |
| `orchestrator-protocol` | Contractor mode: implement → verify → review → fix → score |
| `session-logging` | Three logging triggers: post-plan, incremental, end-of-session |
| `meta-governance` | Template vs working-project decisions (generic vs specific) |

**Path-scoped** (load only when matching files are touched):

| Rule | Triggers on | What it enforces |
|------|------------|-----------------|
| `replication-protocol` | `replications/**`, `scripts/**` | 6-phase replication + RNA-seq pitfall tables |
| `quality-gates` | `*.R`, `*.py`, `reports/**` | 80/90/95 scoring + tolerance thresholds |
| `r-code-conventions` | `*.R` | R standards, reproducibility, RNA-seq pitfalls |
| `python-code-conventions` | `*.py` | Python scientific coding standards |
| `orchestrator-research` | `*.R`, `explorations/**` | Simplified orchestrator for research (no multi-round reviews) |
| `verification-protocol` | `replications/**`, `reports/**` | Replication task completion checklist |
| `markdown-processing` | `master_supporting_docs/`, `papers/` | Safe large-Markdown paper handling |
| `proofreading-protocol` | `*.md`, `quality_reports/**` | Propose-first, then apply with approval |
| `knowledge-base-template` | `*.R`, `*.py`, `replications/**` | Gene-ID registry, normalization registry, pitfalls |
| `exploration-folder-protocol` | `explorations/` | Structured sandbox for experimental work |
| `exploration-fast-track` | `explorations/` | Lightweight exploration workflow (60/100 threshold) |

### Hooks (`.claude/hooks/`)

`log-reminder.py`, `notify.sh`, `post-merge.sh`, `pre-compact.sh`, `protect-files.sh`, `run-python.sh` — session-log reminders, notifications, context-survival prompts, file protection, and a portable Python launcher.

### Templates (`templates/`)

`session-log.md`, `quality-report.md`, `exploration-readme.md`, `archive-readme.md`, `requirements-spec.md`, `constitutional-governance.md`, `skill-template.md`.

</details>

---

## Prerequisites

| Tool | Required for | Install |
|------|-------------|---------|
| [Claude Code](https://docs.anthropic.com/en/docs/claude-code/overview) | Everything | `npm install -g @anthropic-ai/claude-code` |
| R (≥ 4.5) | R replication scripts | [r-project.org](https://www.r-project.org/) |
| Python (≥ 3.11) | Python replication scripts | [python.org](https://www.python.org/) |
| [gh CLI](https://cli.github.com/) | PR workflow | [cli.github.com](https://cli.github.com/) |
| Bioconductor | Bulk RNA-seq: `DESeq2`, `edgeR`, `limma` | `BiocManager::install(c("DESeq2","edgeR","limma"))` |
| Seurat / Bioconductor | scRNA-seq in R: `Seurat`, `SingleCellExperiment`, `scran` | `install.packages("Seurat")` |
| scanpy / anndata | Python scRNA-seq | `pip install scanpy anndata` |
| scvi-tools | Deep-learning single-cell | `pip install scvi-tools` |
| pydeseq2 | Python bulk RNA-seq | `pip install pydeseq2` |

> **Python environments:** this project standardizes on `conda` — invoke Python via `conda run --no-capture-output -n <env> python ...` rather than bare `python`. See `CLAUDE.md`.

---

## Data Setup

Create a `data/` directory at the project root and place datasets there. It is **gitignored** (data can be large or controlled-access).

1. **Public data:** download from GEO (e.g. `GEOquery::getGEO("GSExxxxxx")`) or SRA; record accession IDs in your scripts.
2. **Counts / matrices:** put them in `data/[PaperName]/` with a short README describing the source.
3. **Genome annotation:** record the genome assembly and GTF version (e.g. Ensembl 110, GRCh38) in the script header.
4. **Large files:** consider symlinking `data/` to external storage if files exceed GitHub limits.

---

## Folder Structure

```
paper-replicate-agent-demo/
    ├── papers/                  # Markdown papers + original replication packages
    │   └── [PaperName]/
    │       ├── original_paper.md
    │       ├── original_paper_attachment/   # Figures referenced by the paper
    │       └── code/                        # Original analysis code (if provided)
    ├── data/                    # RNA-seq data (gitignored)
    ├── replications/            # Our R/Python replication scripts + outputs
    ├── reports/                 # Polished final reports
    ├── quality_reports/         # Plans, specs, session logs, replication targets
    ├── explorations/            # Sandbox for experimental analyses
    ├── master_supporting_docs/  # Reference papers and methods docs
    ├── scripts/                 # Utility scripts and shared R functions
    └── .claude/                 # Skills, agents, rules, hooks
```

---

## License

This workflow is released under the **MIT License**; the bundled skills are licensed individually, according to where each one came from:

- **Framework skills** (`replicate-paper`, `data-analysis`, `review-r`, `review-paper`, `proofread`, `devils-advocate`, `commit`) — **MIT**, co-crediting **maxwell2732 (朱晨 Chen Zhu | 遗传社科研究)**, who created the UK Biobank workflow this project was forked from, and **Emberwhirl**, who re-targeted it for RNA-seq.
- **Single-cell skills built here** (`scanpy-basics`, `single-cell-integration-ingest-bbknn`, `kallisto-bustools-import`) — **MIT © Emberwhirl**. They build on the open [scverse](https://scverse.org/) ecosystem; `scanpy-basics`, for instance, is adapted from the official scanpy tutorials (BSD-3-Clause).
- **Skills adapted from Anthropic's life-sciences collection** (`scvi-tools`, `single-cell-rna-qc`) — distributed under their upstream **Apache License 2.0**, each with a `NOTICE` file recording the source ([anthropics/life-sciences](https://github.com/anthropics/life-sciences)) and the local changes.

The single-cell skills *orchestrate* third-party open-source software — scanpy, scvi-tools, bbknn, kb-python, DropletUtils, and friends — which is installed separately rather than bundled here. Those tools keep their own licenses (mostly permissive, a few under other terms); each skill notes its relevant dependencies in its own files.

Use it freely for your research — and a nod to the upstream authors (Chen Zhu, the scverse community, and Anthropic) is always appreciated. 🙏
