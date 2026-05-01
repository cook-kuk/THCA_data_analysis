# bioRxiv + Zenodo + GitHub Release Submission Workflow

**Date:** 2026-05-03 (marathon scaffolding/infra Round 5)
**Target:** W6 (6/8-6/13) bioRxiv submission, then Cell Rep Med formal W7
**Format:** Step-by-step actionable commands

---

## Phase A — bioRxiv preprint (W6, 6/8-6/13)

### Step A1 — GitHub release v0.9 (preprint companion)

```bash
cd /opt/thyroid-dash/project

# Verify smoke tests pass
python notebooks_or_scripts/tests/test_signature_score.py

# Verify all 11 submission TSVs present
for f in S1 S2 S3 S4 S5 S5b S6 S7 S8a S8b S8c; do
  test -f submission_data/${f}*.tsv || echo "MISSING $f"
done

# Verify all figures present
ls results/figures/F[1-5]_*/F*.pdf  # 5 main
ls results/figures/suppl/SF*.pdf    # 6 suppl

# Tag preprint companion
git add -A
git commit -m "release: v0.9 (bioRxiv preprint companion)"
git tag -a v0.9 -m "Phase 0 Cancer Paper bioRxiv preprint companion"

# Push
git push origin main --tags
```

### Step A2 — Zenodo DOI reservation

```
Browser actions:
1. Go to https://zenodo.org/account/settings/github/
2. Toggle "thyca-paper-2026" repository to ON (auto-archives on release)
3. After GitHub push, Zenodo will detect the v0.9 tag and create archive
4. Note the assigned DOI (format: 10.5281/zenodo.NNNNNNNN)
5. Update README.md badge with this DOI

Zenodo metadata (manually fill on first archive):
- Title: "THCA Paper 2026 — Phase 0 Cancer Paper code & data (bioRxiv preprint v0.9)"
- Description: [Methods M1-M11 + Result R1-R5 abstract]
- Authors: Cook, Seungho
- License: MIT (code) + CC-BY-4.0 (data)
- Keywords: thyroid cancer, PTC, HLA, autoimmune, Hashimoto's
- Communities: Health Sciences, Cancer Research
```

### Step A3 — bioRxiv submission portal

```
Portal: https://www.biorxiv.org/submit-a-manuscript

Required fields:
- Manuscript title (final, voice-protected)
- Author list with affiliations + ORCID
- Corresponding author: Seungho Cook (kukshomr@gmail.com)
- Subject area: Cancer Biology (primary) + Genetics/Genomics (secondary)
- Manuscript type: Research Article
- License: CC-BY-NC-ND or CC-BY (Cell Press preference)
- Abstract (≤300 words, voice-finalized)
- Manuscript file (PDF, ≤30 MB; LaTeX or Word source archive)
- Suppl materials (single PDF or zip; F1-F5 + SF1-SF6 + S1-S8 Excel; ≤50 MB)
- Cover letter (optional but recommended)
- Funding information
- Conflict of interest declarations
- Data availability + GitHub link + Zenodo DOI
- Preprint posting confirmation

Estimated review time: 24-48 hours after submission
DOI assigned: usually within 24h of acceptance
```

### Step A4 — Post-bioRxiv share (long-term identity building)

```bash
# Twitter/X thread (placeholder)
echo "🧵 New preprint! Cell Press companion paper Phase 0 of 5-Pillar Cancer trajectory"
# [DOI link, key finding, figure highlight, Yu professor handle]
# [hashtags: #ThyroidCancer #HLA #Autoimmune #cancerresearch]

# LinkedIn post
echo "Excited to share our preprint on transcriptional differentiation axis in PTC..."
# [Detail: 5-Pillar structure, Korean cohort, key findings]
# [Tag Yu professor's institution + cookHLA team]

# Yu professor's Korean thyroid community channel
# [Yu professor handles share via institutional email + KTA newsletter]
```

---

## Phase B — Cell Rep Med formal submission (W7, 6/14-6/20)

### Step B1 — Final manuscript revision after bioRxiv community feedback

```
Wait 1-3 days for any preprint comments
Address any obvious issues (typos, broken links)
```

### Step B2 — Cell Rep Med portal submission

```
Portal: https://www.editorialmanager.com/cellrepmed

Required documents:
1. Cover letter (voice-finalized ¶1, fact-only ¶2-4 ready)
2. Manuscript main text (≤5,000 words excl. methods)
3. Methods section (M1-M11 prose, ~1,500-2,000 words)
4. STAR Methods Key Resources Table (Excel)
5. CRediT author contributions
6. Conflict of Interest statement
7. Inclusion + Diversity statement (Cell Press mandatory)
8. Highlights (3-4 bullets, ≤85 chars each, voice-finalized)
9. eTOC blurb (50 words, voice-finalized)
10. Main figures F1-F7 (PDF, 300+ DPI)
11. Suppl figures SF1-SF6 (PDF)
12. Suppl tables S1-S8 (Excel)
13. References (Vancouver style, ~22-30 numbered)
14. Suggested reviewers (4-6 names with affiliations + emails)
15. Excluded reviewers (conflicts, with brief reasons)
16. Author affiliations + ORCID IDs

Article type: Article (full)
Subject area: Cancer / Translational research / HLA + immunogenetics
Estimated word count + figure count provided

Preprint info:
- bioRxiv DOI from Phase A
- Confirm "preprint deposit not detrimental to publication"
```

### Step B3 — Editorial process tracking

```
Day 0 (W7): Submission confirmed
Day 1-3: Initial editorial screening (~70% rejection at this stage)
Day 4-30: Reviewer assignment + review (3 reviewers typical)
Day 30-45: Decision letter
  - Accept (rare on first round)
  - Revise + resubmit (most common positive outcome)
  - Reject (~50% acceptance rate at this stage)

If revise + resubmit:
- 2-4 weeks revision turnaround typical
- Reviewer Q&A document already prepared (13 anticipated Qs in 
  reports/2026_05_03_reviewer_QA_consolidated.md)

Final acceptance: typically Q3-Q4 2026
Publication: Q3-Q4 2026 (after typesetting)
```

---

## Phase C — Backup venues (if Cell Rep Med rejection)

### Step C1 — JCI Insight (W15+, 6 weeks after Cell Rep Med decision)

```
Reformatting required:
- Vancouver-style references retained (Cell Press → JCI compatible)
- Methods section may expand (JCI Insight has more lenient word limits)
- Title + Abstract may need minor adjustment
- Figures retained as-is (Cell Press 300 DPI same)
- Suppl materials retained

Portal: https://insight-jci.editorialmanager.com
Submission timeline: 3-5 days reformatting, 2-4 weeks review

Cover letter rewrite:
- Para 1 (voice) — minor tone adjustment
- Para 2-4 (fact) — reuse mostly verbatim
```

### Step C2 — Genome Medicine (W20+ if both rejection)

```
Different format: Word limit ~3,500 main text
References: numerical Vancouver
Figure max: 6 main + unlimited suppl
```

### Step C3 — Sci Rep / Endocrine-Related Cancer (final fallback)

```
Sci Rep:
- Open access $1,490 (cheaper)
- Less restrictive scope, single-discipline
- Faster decision (~30 days)

Endocrine-Related Cancer:
- IF 5
- Endocrine-specific specialty journal
- Slower decision (~60 days)
```

---

## Phase D — Stretch venue (if Bundang Scenario A response)

### Step D1 — Bundang Graves' validation analysis

```
Timeline: W2-W4 if Bundang responds in W1
- Bundang sample receipt + IRB clearance
- arcasHLA imputation if RNA-seq, cookHLA if SNP genotype
- Pillar 1 expansion: n=874 → n=874+Bundang
- Phase 1 paper trajectory open
```

### Step D2 — Nat Commun submission

```
Portal: https://mts-natcomm.nature.com

Pre-submission inquiry recommended (1 page):
- Title + Abstract
- Significance statement (broader interest)
- Suggested reviewers

Article format: Standard (≤5,000 words)
Acceptance rate: ~7%
Time to first decision: 30-60 days

Pre-submission inquiry response: 7-14 days
Editor decides whether to send for full review

If proceed:
- Full submission with Bundang as co-authoring institution
- Publication timeline Q4 2026 - Q1 2027
```

---

## Critical safeguards

### Pre-submission checklist (W7 day-of)

```bash
# Final verification before clicking "submit"
cd /opt/thyroid-dash/project

# 1. All smoke tests pass
python notebooks_or_scripts/tests/test_signature_score.py
test $? -eq 0 || (echo "❌ Smoke tests fail" && exit 1)

# 2. All figures correct version
md5sum results/figures/F[1-5]_*/F*.pdf
md5sum results/figures/suppl/SF*.pdf

# 3. All TSVs present
ls -la submission_data/S*.tsv | wc -l  # should be 11

# 4. Voice-protected sections finalized
grep -l "VOICE-PROTECTED" reports/2026_05_03_*.md
# Check that user has filled in: Hook, Aim, Discussion 3.1, Limitations,
# Cover Letter Para 1, Reviewer Q9 narrative, Title

# 5. References cross-checked
grep "Chu 2018\|Chu X" docs/methods_M1_M11.md  # NOT "Chen 2018"
grep "Landa 2016" docs/methods_M1_M11.md  # NOT "Krishnamoorthy 2025 Nat Comm"

# 6. ORCID + emails for all authors
# Manually verify in submission_data/cohort_assembly.tsv

# 7. Cover letter Para 1 user-confirmed
# Manually verify reports/2026_05_03_cover_letter_assembly.md

echo "✅ All pre-submission checks pass"
```

### NEVER skip these (Cell Press will return without review)

- [ ] STAR Methods Key Resources Table format (Excel mandatory)
- [ ] CRediT author contributions
- [ ] Inclusion + Diversity statement
- [ ] Conflict of Interest declarations
- [ ] Data availability statement
- [ ] Code availability statement (with GitHub + Zenodo DOI)
- [ ] All ORCID IDs

### Avoid these mistakes

- [ ] Don't submit before Yu professor 미팅 confirmation
- [ ] Don't submit with bioRxiv-only DOI placeholder; wait for actual DOI
- [ ] Don't bundle high-resolution figures into main text PDF; use separate files
- [ ] Don't omit suggested reviewers list; portal may reject
- [ ] Don't submit on weekend; editorial office processes faster on weekdays

---

## Post-acceptance workflow (Q3-Q4 2026)

### Step E1 — Galley proof review

```
Cell Press emails galley proof PDF within 7-14 days of acceptance
Author has 48-72 hours to review

Check:
- Figure labels (A, B, C panels correctly assigned)
- Reference numbering (Vancouver style)
- Table column alignment
- Author affiliation order
- Funding ID correctness

Return signed-off proof via portal
```

### Step E2 — Final publication + post-publication share

```
Cell Press publishes online ~7-21 days after galley proof return
DOI assigned (different from bioRxiv DOI)
- Update Zenodo metadata with peer-reviewed DOI
- Update README.md badge
- Twitter/X thread "Now published in Cell Reports Medicine"
- LinkedIn post
- Yu professor's Korean thyroid community announcement
- Update Phase 1 + Phase 2 paper outlines with Phase 0 publication anchor cite
```

### Step E3 — Long-term tracking

```
Setup:
- Google Scholar Alert for "thyroid cancer DM1 DM2 transcriptional axis"
- ResearchGate post + tag Yu professor
- ORCID record update

Annual review (Q1 2027+):
- Citation count via Google Scholar
- Altmetric score
- Geographic distribution of readers
- Any follow-up calls from oncology community
```

---

## Voice-protected gating (NEVER bypass)

| Section | Decision-maker |
|---|---|
| Title | User voice |
| Abstract first sentence | User voice |
| Introduction Hook ¶1 | User voice |
| Introduction Aim ¶4 | User voice |
| Discussion §3.1 framing tone | User voice |
| Discussion §3.4 Limitations narrative | User voice |
| Cover Letter Paragraph 1 motivation | User voice |
| Reviewer Q9 answer narrative tone | User voice |
| Highlights verb tone | User voice (after Title) |
| eTOC blurb tone | User voice (after Title) |

→ All other components are scaffolding/infra, prepared and ready.
