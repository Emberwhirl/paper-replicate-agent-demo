---
paths:
  - "master_supporting_docs/**"
  - "papers/**"
---

# Markdown Paper Processing

## The Standard Workflow

**Step 1: Confirm the file is Markdown**
```bash
ls papers/AuthorYear/original_paper.md
```
If only a PDF is available, convert it to Markdown before proceeding:
```bash
# Using pandoc (recommended)
pandoc input.pdf -o output.md

# Or use a PDF-to-Markdown tool of your choice
# Verify headings, tables, and math render correctly after conversion
```

**Step 2: Read the Markdown file**
- Use the Read tool directly on `.md` files — no splitting required
- For very long papers (>500 KB), read in sections:
  - Abstract + Introduction
  - Methods (most critical for replication)
  - Results
  - Supplementary Materials

**Step 3: Extract key information**
- Paper sections: Abstract, Introduction, Methods, Results, Discussion, Supplementary
- For each results section: identify every table and figure with numerical targets
- Record software versions, package names, and data accession numbers from Methods

**Step 4: Save targets**
- Record all replication targets in `quality_reports/[paper_name]_replication_targets.md`

## Handling Large Supplementary Materials

Some bioinformatics papers have extensive supplementary files:
- **Supplementary Tables:** Usually provided as CSV/Excel; download and inspect directly
- **Supplementary Methods:** Often embedded in the main `.md` or as a separate `.md` file
- **Supplementary Figures:** Note figure numbers but don't need to read image files

## Error Handling

**If Markdown conversion from PDF is imperfect:**
1. Note sections where conversion may have garbled math or tables
2. Manually check Methods and Results sections for correctness
3. Document any reading gaps in `quality_reports/[paper_name]_data_audit.md`

**If the paper is only available as PDF:**
1. Convert to Markdown using pandoc: `pandoc paper.pdf -o paper.md`
2. Verify the conversion quality, especially for equations and tables
3. Save the `.md` file to `papers/[PaperName]/original_paper.md`

## Important: No PDF Required

This workflow uses Markdown as the primary input format. Do NOT attempt to read PDF files directly. Always convert to Markdown first.

