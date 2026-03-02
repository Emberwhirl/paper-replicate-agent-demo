---
paths:
  - "Slides/**/*.tex"
  - "Quarto/**/*.qmd"
  - "scripts/**/*.R"
  - "replications/**/*.R"
  - "replications/**/*.py"
---

# Project Knowledge Base: RNA-seq Replication

<!-- Fill in the tables below with YOUR project-specific content.
     Claude reads this before creating/modifying any replication content. -->

## Notation Registry

| Rule | Convention | Example | Anti-Pattern |
|------|-----------|---------|-------------|
| Gene identifiers | Ensembl IDs as primary keys | ENSG00000141510 | Gene symbols (not unique) |
| Log fold change | log2 FC unless stated otherwise | LFC = 1.5 | Mixing log2 and ln |
| Adjusted p-values | BH-adjusted unless paper specifies | padj < 0.05 | Unadjusted p-values for DE |

## Symbol / Abbreviation Reference

| Symbol | Meaning | Context |
|--------|---------|---------|
| LFC | Log2 fold change | DE analysis |
| padj | BH-adjusted p-value | DE analysis |
| VST | Variance-stabilizing transformation (DESeq2) | Normalization |
| CPM | Counts per million | Normalization |
| TPM | Transcripts per million | Normalization |
| HVG | Highly variable genes | scRNA-seq feature selection |
| UMAP | Uniform Manifold Approximation and Projection | Dimensionality reduction |
| PCA | Principal component analysis | Dimensionality reduction |

## Active Replication Projects

| Paper | Data Type | Accession | Status | Notes |
|-------|-----------|-----------|--------|-------|
| | bulk RNA-seq | | | |
| | scRNA-seq | | | |

## Normalization Methods Registry

| Method | Package | When to Use | Key Parameter |
|--------|---------|-------------|--------------|
| DESeq2 VST | DESeq2 | Bulk RNA-seq visualization + DE | `blind=FALSE` for DE-aware |
| DESeq2 rlog | DESeq2 | Small sample bulk RNA-seq | Slow on large datasets |
| TMM | edgeR | Bulk RNA-seq with edgeR | Default in `calcNormFactors` |
| log-normalize | Seurat / scanpy | scRNA-seq | `scale.factor=10000` |
| scran pooling | scran | scRNA-seq (reduces zero-inflation bias) | `min.mean=0.1` |

## Pitfalls Registry

| Pitfall | Impact | Prevention |
|---------|--------|------------|
| Wrong reference level | Reversed LFC direction | Always `relevel()` to match paper |
| No seed for clustering/UMAP | Non-reproducible results | `set.seed()` + `seed.use` in Seurat |
| Mixing gene ID versions | Gene ID mismatches | Lock annotation version in script header |
| Ambient RNA not removed | Inflated gene expression in scRNA-seq | Apply CellBender/SoupX only if paper does |
| Doublets not filtered | Artifactual cell clusters | Match paper's doublet tool and thresholds |

## Anti-Patterns (Don't Do This)

| Anti-Pattern | What Happened | Correction |
|-------------|---------------|-----------|
| Using gene symbols as join keys | Many-to-one / missing matches | Use Ensembl IDs as primary keys |
| Running DE without pre-filtering | Inflated FDR, slow computation | Apply `filterByExpr` or equivalent first |
| Hardcoding Bioconductor versions | Breaks on updates | Document versions; use `BiocManager::install(..., version=)` |

## R Code Pitfalls

| Bug | Impact | Fix |
|-----|--------|-----|
| `DESeqDataSetFromMatrix` without `~1` design | Error if no covariates | Use `design = ~1` for intercept-only |
| `FindAllMarkers` without `only.pos=TRUE` | Slow and bloated output | Set `only.pos=TRUE` or `min.pct=0.25` |
| `RunUMAP` without `seed.use` | Non-reproducible UMAP | Always set `seed.use = 42` or paper's seed |
