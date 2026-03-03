# CLAUDE.MD -- RNA-seq Empirical Replication Agent

**Project:** RNA-seq Empirical Replication Agent
**Institution:** [YOUR INSTITUTION]
**Branch:** main

---

## Core Principles

- **Plan first** -- enter plan mode before non-trivial tasks; save plans to `quality_reports/plans/`
- **Verify after** -- run scripts and confirm outputs match targets at the end of every task
- **Replicate before extending** -- match published results exactly before any modifications
- **Quality gates** -- nothing ships below 80/100
- **[LEARN] tags** -- when corrected, save `[LEARN:category] wrong → right` to MEMORY.md

---

## Folder Structure

```
my-rnaseq-agent/
├── CLAUDE.md                    # This file
├── .claude/                     # Rules, skills, agents, hooks
├── papers/                      # Source Markdown papers and original replication packages
│   └── [PaperName]/
│       ├── original_paper.md           # Paper in Markdown format
│       ├── supplementary.md            # Supplementary material in Markdown format
│       ├── original_paper_attachment/  # Figures/images referenced by the Markdown files
│       ├── code/                       # Original analysis code (if provided — may be a package)
│       │   ├── *.R / *.py              # Individual scripts, or full package structure
│       │   └── README.md
│       └── README.md
├── data/                        # Datasets (gitignored — can be large/controlled-access)
├── replications/                # Our replication scripts and outputs
│   └── [PaperName]/
│       ├── R/replicate.R
│       ├── R/figures/
│       ├── R/results/
│       ├── python/replicate.py  # (if Python replication needed)
│       └── validation_report.md
├── reports/                     # Polished final replication reports
├── scripts/                     # Utility scripts (quality_score.py, helpers)
│   └── R/                       # Shared R utility functions
├── quality_reports/             # Plans, session logs, replication targets
│   ├── plans/
│   ├── specs/
│   ├── session_logs/
│   ├── merges/
│   └── [Paper]_replication_targets.md
├── explorations/                # Exploratory analysis sandbox
│   └── ARCHIVE/
├── master_supporting_docs/      # Methodology reference papers and slides
│   ├── supporting_papers/
│   └── supporting_slides/
└── templates/                   # Session log, quality report templates
```

---

## Commands

```bash
# Run R replication script
Rscript replications/[PaperName]/R/replicate.R

# Run Python replication script
python replications/[PaperName]/python/replicate.py

# Quality score
python scripts/quality_score.py replications/[PaperName]/R/replicate.R
```

---

## Quality Thresholds

| Score | Gate | Meaning |
|-------|------|---------|
| 80 | Commit | Good enough to save |
| 90 | PR | Ready for deployment |
| 95 | Excellence | Aspirational |

---

## Skills Quick Reference

| Command | What It Does |
|---------|-------------|
| `/replicate-paper [paper.md] [data] [code]` | Full 6-phase replication pipeline |
| `/data-analysis [dataset]` | End-to-end R analysis |
| `/review-r [file]` | R code quality review |
| `/review-paper [file]` | Manuscript review |
| `/commit [msg]` | Stage, commit, PR, merge |
| `/proofread [file]` | Grammar/typo review of reports |
| `/devils-advocate [topic]` | Challenge design decisions |

---

## Active Replications

| Paper | Status | Targets | Pass | Fail | Notes |
|-------|--------|---------|------|------|-------|
| [AuthorYear] | — | — | — | — | bulk/scRNA-seq |

---

## RNA-seq Data Notes

- **Data location:** `data/` (gitignored — can be large or controlled-access)
- **Public data:** Download from GEO (accession GSExxxxxx) or SRA; document accession IDs in replication scripts
- **Genome/annotation:** Record genome assembly (e.g., GRCh38) and GTF version (e.g., Ensembl 110) used for alignment/quantification
- **Gene IDs:** Use Ensembl gene IDs as primary keys; map to gene symbols only for display; verify ID versions match between paper and your annotation
- **Normalization:** Match paper's normalization method exactly (raw counts, CPM, TPM, VST, log-normalized); document deviations
- **Batch correction:** Apply batch correction (e.g., ComBat, Harmony) only if the paper does; document batch variables used
- **Filtering:** Match paper's minimum count/cell thresholds before running differential expression or clustering
