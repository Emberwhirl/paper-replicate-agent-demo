---
name: review-paper
description: Comprehensive manuscript review covering argument structure, computational methods, statistical correctness, citation completeness, and potential referee objections. Optimized for bioinformatics and computational biology papers.
disable-model-invocation: true
argument-hint: "[paper filename in papers/ or master_supporting_docs/ — must be .md format]"
allowed-tools: ["Read", "Grep", "Glob", "Write", "Task"]
---

# Manuscript Review

Produce a thorough, constructive review of an academic manuscript — the kind of report a top computational biology journal referee would write.

**Input:** `$ARGUMENTS` — path to a paper (`.md`), or a filename in `papers/` or `master_supporting_docs/`. Papers must be in Markdown format.

---

## Steps

1. **Locate and read the manuscript.** Check:
   - Direct path from `$ARGUMENTS`
   - `papers/$ARGUMENTS/original_paper.md`
   - `master_supporting_docs/supporting_papers/$ARGUMENTS`
   - Glob for partial matches

2. **Read the full paper** end-to-end (Markdown sections: Abstract → Introduction → Methods → Results → Discussion → Supplementary).

3. **Evaluate across 6 dimensions** (see below).

4. **Generate 3-5 "referee objections"** — the tough questions a top referee would ask.

5. **Produce the review report.**

6. **Save to** `quality_reports/paper_review_[sanitized_name].md`

---

## Review Dimensions

### 1. Argument Structure
- Is the research question clearly stated?
- Does the introduction motivate the question effectively?
- Is the logical flow sound (question → method → results → conclusion)?
- Are the conclusions supported by the evidence?
- Are limitations acknowledged?

### 2. Computational Methods
- Is the novel method (if any) clearly described step-by-step?
- Are key algorithmic choices (distance metrics, model assumptions, optimization strategy) justified?
- Are comparison to existing methods (DESeq2 vs. edgeR; Seurat vs. scran) fair and comprehensive?
- Is the method's time and memory complexity discussed?
- Is software availability and documentation addressed (GitHub link, version)?

### 3. Statistical Correctness
- Is the normalization method appropriate for the data type and downstream analysis?
- Are multiple testing corrections applied (BH, Bonferroni, IHW)?
- Are batch effects accounted for if present?
- Is the sample size sufficient for the claims made?
- Are robustness checks or sensitivity analyses included?

### 4. Literature Positioning
- Are the key papers cited (DESeq2, edgeR, Seurat, scanpy, relevant method papers)?
- Is prior work characterized accurately?
- Is the contribution clearly differentiated from existing work?
- Any missing citations that a referee would flag?

### 5. Writing Quality
- Clarity and concision
- Academic tone
- Consistent notation throughout
- Abstract effectively summarizes the paper
- Tables and figures are self-contained (clear labels, notes, accession numbers)

### 6. Reproducibility and Data Availability
- Is data deposited in a public repository (GEO, SRA, Zenodo)?
- Is accession number provided?
- Is code available (GitHub, Zenodo)?
- Is the software environment (R/Python versions, package versions) documented?
- Would an independent group be able to reproduce the main results?

---

## Output Format

```markdown
# Manuscript Review: [Paper Title]

**Date:** [YYYY-MM-DD]
**Reviewer:** review-paper skill
**File:** [path to manuscript]

## Summary Assessment

**Overall recommendation:** [Strong Accept / Accept / Revise & Resubmit / Reject]

[2-3 paragraph summary: main contribution, strengths, and key concerns]

## Strengths

1. [Strength 1]
2. [Strength 2]
3. [Strength 3]

## Major Concerns

### MC1: [Title]
- **Dimension:** [Methods / Statistics / Argument / Literature / Writing / Reproducibility]
- **Issue:** [Specific description]
- **Suggestion:** [How to address it]
- **Location:** [Section/page/table if applicable]

[Repeat for each major concern]

## Minor Concerns

### mc1: [Title]
- **Issue:** [Description]
- **Suggestion:** [Fix]

[Repeat]

## Referee Objections

These are the tough questions a top referee would likely raise:

### RO1: [Question]
**Why it matters:** [Why this could be fatal]
**How to address it:** [Suggested response or additional analysis]

[Repeat for 3-5 objections]

## Specific Comments

[Line-by-line or section-by-section comments, if any]

## Summary Statistics

| Dimension | Rating (1-5) |
|-----------|-------------|
| Argument Structure | [N] |
| Computational Methods | [N] |
| Statistical Correctness | [N] |
| Literature | [N] |
| Writing | [N] |
| Reproducibility | [N] |
| **Overall** | **[N]** |
```

---

## Principles

- **Be constructive.** Every criticism should come with a suggestion.
- **Be specific.** Reference exact sections, equations, tables.
- **Think like a referee at Nature Methods or Genome Biology.** What would make them reject?
- **Distinguish fatal flaws from minor issues.** Not everything is equally important.
- **Acknowledge what's done well.** Good research deserves recognition.
- **Do NOT fabricate details.** If you can't read a section clearly, say so.
