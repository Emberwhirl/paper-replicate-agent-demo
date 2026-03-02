---
paths:
  - "master_supporting_docs/**"
  - "papers/**"
---

# Markdown Paper Processing

## Expected Paper Structure

Papers must be provided as Markdown (`.md`) files before this workflow is invoked. The user is responsible for converting PDFs to Markdown using their preferred tool (e.g., docling) before starting a replication session.

Expected folder layout for each paper:

```
papers/[PaperName]/
├── original_paper.md           # Full paper in Markdown format
├── supplementary.md            # Supplementary methods/results (if provided separately)
├── original_paper_attachment/  # Figures and images referenced by the Markdown files
│   ├── figure1.png
│   └── ...
├── code/                       # Original analysis code (if provided by authors)
│   ├── *.R / *.py              # Individual scripts, or a full R/Python package
│   └── README.md
└── README.md                   # Notes on data source, accession, conversion details
```

## The Standard Workflow

**Step 1: Confirm the Markdown file exists**
```bash
ls papers/AuthorYear/original_paper.md
```
If the file is missing, ask the user to provide it. Do NOT attempt to read or convert PDF files.

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
- If figures are referenced, look for them in `papers/[PaperName]/original_paper_attachment/`

**Step 4: Save targets**
- Record all replication targets in `quality_reports/[paper_name]_replication_targets.md`

## Handling Large Supplementary Materials

Some bioinformatics papers have extensive supplementary files:
- **Supplementary Tables:** Usually provided as CSV/Excel; download and inspect directly
- **Supplementary Methods:** Often embedded in the main `.md` or as a separate `.md` file
- **Supplementary Figures:** Note figure numbers; images are in `original_paper_attachment/` if provided

## Important: Markdown Only

This workflow accepts Markdown as the sole paper input format. If a user provides a PDF path, ask them to convert it to Markdown first using a tool like docling, then re-invoke the skill with the `.md` path.

