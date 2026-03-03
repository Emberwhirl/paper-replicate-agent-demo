# Session Log: Workflow Adaptation for RNA-seq Replication
**Date:** 2026-02-20
**Goal:** Adapt lecture-slide repo for empirical replication of computational biology papers using bulk and single-cell RNA-seq data.
**Plan:** `quality_reports/plans/` (approved in preceding plan session)

---

## Summary

Implemented full workflow adaptation. 7 files modified, 2 files created, 4 directories added.

## Changes Made

### Files Modified
1. **CLAUDE.md** — Rewrote for RNA-seq replication context: new folder structure, Python/R commands, replication-focused skills table, "Active Replications" tracker
2. **.claude/WORKFLOW_QUICK_REF.md** — Filled all placeholders: pathlib/here conventions, seeds, 300 DPI figures, Okabe-Ito palette, tolerance thresholds
3. **.claude/agents/domain-reviewer.md** — Customized for senior computational biology referee (Cell/Nature/Science standard): Lens 1 = data processing/QC, Lens 2 = normalization/statistical model, Lens 4 = code-method alignment pitfalls
4. **.claude/rules/quality-gates.md** — Added Python scoring table, updated paths to include replications/, filled tolerance thresholds
5. **.claude/rules/replication-protocol.md** — Added Phase 0 (Paper Intake), expanded RNA-seq pitfall tables for bulk and single-cell, updated report template to include Python script path
6. **.gitignore** — Added `data/` (datasets may be large or controlled-access), `replications/**/*.parquet`, `.pkl`, `.rds`

### Files Created
7. **.claude/skills/replicate-paper/SKILL.md** — 6-phase replication pipeline: Intake → Data Audit → Code Analysis → Translation → Validation → Report
8. **.claude/rules/python-code-conventions.md** — Python scientific coding standards: pathlib, seeds, imports, pandas discipline, statsmodels/lifelines/pymc, 300 DPI figures, Okabe-Ito palette, WHY-not-WHAT comments, script template

### Directories Created
- `papers/` — Source PDFs and replication packages
- `data/` — Datasets (gitignored)
- `replications/` — Replication scripts and outputs
- `reports/` — Polished final reports

## Design Decisions

- **Left lecture-slide skills in place** (compile-latex, deploy, etc.) — dead infrastructure but harmless; avoids breakage risk from deletion
- **data/ gitignored** — Datasets may be large or controlled-access; must never be committed
- **Replication intermediate outputs gitignored** (parquet, pkl, rds in replications/) — large binary files; reproducible from script
- **settings.json pre-modified** by user before session: replaced LaTeX/Quarto permissions with Python tooling (python, pip, pytest, ruff, black, mypy)

## Verification

- No stale Beamer/LaTeX/Quarto references in CLAUDE.md: ✓
- replicate-paper skill has 6 phases: ✓
- data/ in .gitignore: ✓
- RNA-seq pitfall tables in replication-protocol: ✓
- git diff --stat: 7 files changed + 5 untracked new files: ✓

## Status

Implementation complete. Ready for first `/replicate-paper` invocation.

---
**Context compaction () at 21:20**
Check git log and quality_reports/plans/ for current state.
