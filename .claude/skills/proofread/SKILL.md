---
name: proofread
description: Run the proofreading protocol on Markdown reports and replication write-ups. Checks grammar, typos, consistency, and academic writing quality. Produces a report without editing files.
disable-model-invocation: true
argument-hint: "[filename or path to .md report]"
allowed-tools: ["Read", "Grep", "Glob", "Write", "Task"]
---

# Proofread Markdown Reports

Run the mandatory proofreading protocol on Markdown reports and replication write-ups. This produces a report of all issues found WITHOUT editing any source files.

## Steps

1. **Identify files to review:**
   - If `$ARGUMENTS` is a specific filename or path: review that file only
   - If `$ARGUMENTS` is "all": review all `.md` files in `reports/` and `replications/`

2. **For each file, launch the proofreader agent** that checks for:

   **GRAMMAR:** Subject-verb agreement, articles (a/an/the), prepositions, tense consistency
   **TYPOS:** Misspellings (including tool/gene names), duplicated words, punctuation errors
   **CONSISTENCY:** Gene ID format, terminology, notation, table/figure cross-references
   **ACADEMIC QUALITY:** Informal language, incomplete sentences, unsupported claims
   **REPLICATION REPORT:** Accession numbers, genome version, seed, verdict consistency

3. **Produce a detailed report** for each file listing every finding with:
   - Location (section heading or line number)
   - Current text (what's wrong)
   - Proposed fix (what it should be)
   - Category and severity

4. **Save each report** to `quality_reports/`:
   - `quality_reports/[FILENAME_WITHOUT_EXT]_proofread.md`

5. **IMPORTANT: Do NOT edit any source files.**
   Only produce the report. Fixes are applied separately after user review.

6. **Present summary** to the user:
   - Total issues found per file
   - Breakdown by category
   - Most critical issues highlighted
