---
paths:
  - "Quarto/**/*.qmd"
  - "replications/**"
  - "reports/**"
  - "docs/**"
---

# Task Completion Verification Protocol

**At the end of EVERY task, Claude MUST verify the output works correctly.** This is non-negotiable.

## For Quarto/HTML Slides:
1. Run `./scripts/sync_to_docs.sh` (or `./scripts/sync_to_docs.sh LectureN`) to render and deploy
2. Open the HTML in browser: `open docs/slides/LectureX.html`
3. Verify images display by reading 2-3 image files to confirm valid content
4. Check HTML source for correct image paths
5. Check for overflow by scanning dense slides
6. Report verification results

## For R Replication Scripts:
1. Run `Rscript replications/[paper]/R/replicate.R`
2. Verify output files (RDS, CSV, PNG/PDF figures) were created with non-zero size
3. Spot-check key results (DE gene counts, LFCs, cluster counts) for reasonable values
4. Confirm figures are saved at the expected paths

## For Python Replication Scripts:
1. Run `python replications/[paper]/python/replicate.py`
2. Verify output files (h5ad, parquet, PNG figures) were created with non-zero size
3. Spot-check key results for reasonable values

## For Replication Reports (.md):
1. Open the report and confirm all sections are present
2. Verify the results comparison table is populated (not empty)
3. Confirm the verdict (REPLICATED / PARTIAL / FAILED) is stated
4. Confirm the environment section lists package versions and data accession

## Common Pitfalls:
- **Assuming success**: Always verify output files exist AND contain correct content
- **Missing seed**: Stochastic steps (UMAP, clustering) produce different results without fixed seed
- **Wrong normalization**: Check that the normalization step matches the paper before running DE

## Verification Checklist:
```
[ ] Output file created successfully
[ ] No runtime errors
[ ] Figures/results saved at expected paths
[ ] Key statistics within tolerance of paper targets
[ ] Reported results to user
```
