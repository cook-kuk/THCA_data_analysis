---
title: "Decision lock — Lee 2014 fill (Pillar I v2 G3)"
date: 2026-05-04
status: DECISION ONLY — execution not started. Documents the user's (b) approval; subsequent actual baseline lookup + integration awaits G2 + G6 + advisor sign-off.
related_spec: PAPER2_PILLAR1_V2_FOREST_META_SPEC_2026_05_04.md §4 (Alleles), §8 (Decision gates G2-G6)
---

# Decision lock — Lee 2014 fill (Pillar I v2 Decision Gate G3)

## Decision (locked 2026-05-04)

**G3 = YES — include Lee 2014 fill → 6-allele forest path approved.**

Forest meta will report 6 alleles (not 4):

1. DPB1\*05:01
2. A\*02:07
3. B\*46:01
4. DRB1\*07:01
5. **C\*01:02** ← added via Lee 2014 fill
6. **DQB1\*02:01** ← added via Lee 2014 fill

## What this decision means (and does not mean)

**This decision IS:**
- An approval to use Lee 2014 (Lee MN et al. *Tissue Antigens* 2014) Korean baseline frequencies for C\*01:02 and DQB1\*02:01 as part of the AFND South Korea baseline pool when Pillar I v2 forest meta is eventually executed.
- A commitment to a 6-allele forest plot (Figure 1 / Pillar I figure for Paper 2).

**This decision IS NOT:**
- An execution authorization. No baseline lookup, no Korean PTC pool allele frequency computation, no forest plot, no manuscript paragraph fill performed yet.
- A bypass of remaining decision gates G1 (advisor approval), G2 (baseline pool composition), G4 (Harbin exclusion), G5 (Figure 1 re-render), G6 (manuscript paragraph framing).
- A commitment to specific Lee 2014 frequency values; those will be extracted at execution time and verified against the published table.

## Lee 2014 reference (verification required at execution)

**Citation (to verify at execution)**:
- Lee MN et al. (or Lee KU / Lee HJ — first author exact verification needed)
- *Tissue Antigens* 2014
- KOTRY (Korean Organ Transplantation Registry) donor cohort, Korean general population reference
- Reports HLA-A, B, C, DRB1, DQB1, DPB1 4-digit allele frequencies
- N typically reported per allele class (~5,000+ for HLA-A/B/C; ~2,000-5,000 for HLA-DRB1/DQB1/DPB1)

**Verification gate**:
At execution, the citation will be confirmed against PubMed (no current PMID hardcoded — avoid wrong-cite repeat of Krishnamoorthy 2025 / Chen 2018 misattributions).

**Expected baseline frequency ranges (per general published Korean HLA literature, NOT yet verified against Lee 2014 specifically)**:
- C\*01:02: ~12-18% Korean carrier rate (TBD verify)
- DQB1\*02:01: ~10-15% Korean carrier rate (TBD verify)

These ranges are place-holders — actual values will be extracted from Lee 2014 Table N (TBD) at execution.

## Cross-reference with G2 (AFND South Korea baseline pool composition)

G3 YES decision implies that **Lee 2014** is included as ONE of the AFND South Korea baseline pool sources (alongside KOTRY / KBP / cookHLA — composition decision still open at G2).

- If G2 selects "Lee 2014 only" → Lee 2014 frequencies become the entire baseline reference (single-source)
- If G2 selects "AFND pooled (KOTRY + KBP + cookHLA + Lee 2014)" → Lee 2014 frequencies are merged via inverse-variance weighting per allele
- If G2 selects "exclude Lee 2014" → contradicts G3, requires reopening

G2 decision still pending advisor input.

## Status summary

| Gate | Status |
|---|---|
| G1 advisor approval | ⏸ pending |
| G2 baseline pool composition | ⏸ pending |
| **G3 Lee 2014 fill** | ✅ **YES (locked 2026-05-04)** |
| G4 Harbin Korean fallback exclusion | ⏸ pending |
| G5 Figure 1 re-render | ⏸ pending |
| G6 manuscript paragraph framing | ⏸ pending |

**No execution started. Spec frozen. Awaiting remaining gates.**
