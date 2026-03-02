---
paths:
  - "replications/**/*.R"
  - "replications/**/*.py"
  - "scripts/**/*.R"
  - "scripts/**/*.py"
---

# Replication-First Protocol

**Core principle:** Replicate original results to the dot BEFORE extending.

---

## Phase 0: Paper Intake

Before any coding:

- [ ] Read the full paper Markdown file; identify **every table and figure** that presents empirical results
- [ ] Record gold standard values in `quality_reports/[paper_name]_replication_targets.md`:

```markdown
## Replication Targets: [Paper Author (Year)]

| Target | Table/Figure | Value | Stat | N | Notes |
|--------|-------------|-------|------|---|-------|
| Top DE gene LFC | Table 2, Col 1 | 3.42 | padj=0.0001 | 12,845 cells | Primary comparison |
```

- [ ] Note the **Methods section** in full: organism, tissue, data type (bulk/scRNA-seq), QC thresholds, normalization method, model type, statistical test, software versions
- [ ] Identify the **original code language** (R / Python) and whether a replication package / GEO supplementary code is provided
- [ ] If the paper introduces a novel computational method, flag for Phase 1b (methods explanation) in the replicate-paper skill

---

## Phase 1: Inventory & Data Audit

- [ ] Read the paper's replication README (if provided)
- [ ] Inventory replication package: language, data files, scripts, outputs
- [ ] Locate or download data (GEO accession, SRA run IDs); document in script header
- [ ] Load provided dataset; compare to paper's described sample:
  - N (total cells/samples, per group, after QC)
  - Library size distribution, mitochondrial fraction, gene detection rate
  - QC thresholds — apply them in the paper's stated order
- [ ] Document any **discrepancies between available data and paper description** before coding

---

## Phase 2: Translate & Execute

- [ ] Follow `r-code-conventions.md` and `python-code-conventions.md` for all coding standards
- [ ] Translate line-by-line initially — **do NOT improve during replication**
- [ ] Match original specification exactly: normalization method, feature selection, model formula, test type
- [ ] Save all intermediate datasets as `.rds` (R) or `.h5ad` / `.parquet` (Python)

### Bulk RNA-seq Pitfalls

| Step | Common Trap | Correct Approach |
|------|-------------|-----------------|
| Normalization | Using CPM when paper uses VST or TMM | Match paper's normalization exactly; document choice |
| Filtering | Applying different min-count threshold | Use paper's exact `filterByExpr` or count cutoff |
| Model formula | Omitting batch or covariates | Match paper's design matrix exactly |
| DE test | Using Wald test when paper uses LRT | Specify `test="LRT"` if paper used likelihood ratio |
| Log fold change | Using unshrunken LFC instead of `lfcShrink` | Apply shrinkage only if paper does; use same shrinkage method |
| Reference level | Wrong reference level in factor | Explicitly `relevel()` to match paper's reference group |
| Multiple testing | Using BH when paper uses Bonferroni | Match paper's `p.adjust.method` exactly |

### Single-Cell RNA-seq Pitfalls

| Step | Common Trap | Correct Approach |
|------|-------------|-----------------|
| Doublet removal | Skipping doublet detection if paper applies it | Apply same doublet tool (DoubletFinder, scDblFinder) with same parameters |
| Normalization | Using `NormalizeData` (log-normalize) when paper uses scran | Match normalization method; pooling-based scran ≠ per-cell log-normalize |
| Highly variable genes | Different n_top_genes / flavor | Match paper's HVG selection method and number exactly |
| PCA dims | Using different number of PCs for downstream | Match paper's `dims` parameter for neighbor graph |
| Clustering resolution | Using default resolution instead of paper's value | Set `resolution` to paper's value; results are resolution-sensitive |
| UMAP seed | No fixed seed for UMAP | Set `seed.use` (Seurat) or `random_state` (scanpy) to match paper |
| Cluster labels | Misassigning cell type labels | Verify marker genes used for annotation match paper's Supplementary Table |
| Trajectory | Wrong root cell or start cluster | Set root as specified in paper; trajectory topology is root-sensitive |

### Python-Specific Pitfalls

| Step | Common Trap | Correct Approach |
|------|-------------|-----------------|
| AnnData layers | Using `.X` when paper uses raw counts in `layers["counts"]` | Check which layer the paper normalizes from |
| scanpy neighbors | Different `n_neighbors` or metric | Match `sc.pp.neighbors` parameters exactly |
| pydeseq2 | Different convergence tolerance | Note package version; results may differ slightly from R DESeq2 |

---

## Phase 3: Verify Match

### Tolerance Thresholds

| Type | Tolerance | Rationale |
|------|-----------|-----------|
| Integers (N cells, samples, genes) | Exact match | No reason for any difference |
| Log fold changes (LFC) | ±0.05 | Rounding in paper display + shrinkage estimator variation |
| Adjusted p-values | Same significance bracket (< 0.05, < 0.01, < 0.001) | Exact p may differ across software versions |
| Number of DE genes | ±5% | Minor filtering differences |
| Cluster count | Exact match | Resolution is deterministic given seed |
| Cluster proportions | ±1pp | Rounding in paper display |

### If Mismatch

**Do NOT proceed to extensions.** Isolate which step introduces the difference:
1. Sample size mismatch → check QC thresholds, filtering order, minimum count/cell parameters
2. LFC mismatch → check normalization method, reference level, LFC shrinkage
3. Adjusted p-value bracket mismatch → check multiple testing method, pre-filtering of genes
4. Cluster count mismatch → check resolution, number of PCs, neighbor graph parameters, random seed

Document all investigations even if unresolved.

### RNA-seq-Specific Considerations

**Bulk RNA-seq:**
- **Genome/annotation version:** Verify alignment and GTF versions match paper; gene-level counts depend on annotation
- **Strand specificity:** Confirm strandedness setting (unstranded, forward, reverse) matches paper's library prep protocol
- **Pseudocount:** Match paper's pseudocount for log-transformation (typically +1 but confirm)
- **Batch variable:** Include batch covariate only if paper includes it; batch-corrected vs. uncorrected gives different results
- **Replicate structure:** Verify sample groupings and contrast direction match paper exactly

**Single-cell RNA-seq:**
- **Cell ranger / STARsolo version:** Counts can differ between versions; note version used
- **Ambient RNA removal:** Apply CellBender or SoupX only if paper does; document parameters
- **Doublet method:** Match tool (DoubletFinder, scDblFinder) and `pN` / `pK` parameters
- **Cell cycle regression:** Apply only if paper applies it; affects clustering
- **Integration method:** Match paper's batch integration tool (Harmony, Seurat CCA, scVI); parameters matter
- **Cell type annotation:** Use same marker gene lists; note if labels differ — do not relabel silently

### Replication Report

Save to `replications/[paper_name]/validation_report.md` AND final polished version to `reports/[paper_name]_replication_report.md`:

```markdown
# Replication Report: [Paper Author (Year)]
**Date:** [YYYY-MM-DD]
**Original language:** [R / Python]
**R replication:** [replications/[paper]/R/replicate.R]
**Python replication:** [replications/[paper]/python/replicate.py]

## Summary
- **Targets checked / Passed / Failed:** N / M / K
- **Overall:** [REPLICATED / PARTIAL / FAILED]

## Results Comparison

| Target | Table/Figure | Paper | Ours | Diff | Status |
|--------|-------------|-------|------|------|--------|

## Discrepancies (if any)
- **Target:** X | **Investigation:** ... | **Resolution:** ...

## Environment
- R version, Python version, key packages (with versions), data accession, genome assembly, GTF version, seed
```

---

## Phase 4: Only Then Extend

After replication is verified (all targets PASS or discrepancies documented with explanation):

- [ ] Commit replication scripts: `"Replicate [Paper] Table X -- all targets match"`
- [ ] Now extend with additional analyses (alternative normalizations, additional cell types, gene set enrichment)
- [ ] Each extension builds on the verified baseline
