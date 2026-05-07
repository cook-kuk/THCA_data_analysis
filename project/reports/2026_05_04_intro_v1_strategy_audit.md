---
title: "Introduction v1 audit against final Paper 1 strategy memo"
date: 2026-05-04
target: project/manuscript_v8/03_introduction.md (clean draft, status updated 2026-05-04)
audit_type: section-by-section revision plan ONLY (no prose rewrite, no new intro generation, no manuscript modification)
strategy_memo: 8 rules (driver-orthogonal axis lead / 8-gene = readout not discovery / BRAF-RAS-neg = motivation not headline / H&E-DM1 dropped / fusion claims not headline / methylation not demonstrated → "consistent with DNMT-associated repression" only / TROP2 tumor-level / GSE76039 external advanced-disease replication not continuum proof)
status: AUDIT ONLY — no edits applied
---

# Introduction v1 — Strategy compliance audit

---

## Strategy reference (8 rules)

| # | Rule |
|---|---|
| 1 | Lead with **driver-orthogonal transcriptional differentiation axis** |
| 2 | 8-gene = **compact RAI-lineage readout**, not discovery premise |
| 3 | BRAF/RAS-negative dark matter = **motivation/supplement**, not headline |
| 4 | **H&E-DM1 dropped** |
| 5 | **Fusion-map claims NOT headline** (fusion calls incomplete) |
| 6 | **Methylation NOT demonstrated** — use "consistent with DNMT-associated repression" only |
| 7 | TROP2 = **tumor-level** vulnerability, not spot-level |
| 8 | GSE76039 = **external advanced-disease replication**, NOT monotonic PTC→PDTC→ATC proof |

---

## 1.1 Clinical context — verdict: MODIFY

### Risky phrases

| Phrase (exact) | Why risky | Strategy rule |
|---|---|---|
| "...for a substantial **BRAF/RAS-negative compartment without mechanistic sub-stratification**." (Hook) | Headlines BRAF/RAS-negative compartment as the primary unmet need; strategy positions it as motivation/supplement only | Rule 3 |
| "**5-20% develop recurrent or persistent disease and approximately 10% develop distant metastases**" | Numerical claim attributed solely to Haugen 2016 — verify whether these specific ranges are given in ATA 2015 or sourced elsewhere | Cite verification |
| "**Bethesda III/IV indeterminate cytology affects 15-30% of fine-needle aspiration biopsies**" | Specific 15-30% range needs Cibas 2017 verification (currently `[Cibas 2017 — verify]` placeholder) | Cite verification |
| "[**SEER cite — verify**]" | Placeholder — incidence claim without resolved citation | Cite verification |

### Safe replacement concept (no prose)

- **Hook reframe**: lead with the **mechanistic-stratification gap on the differentiation axis itself**, not the BRAF/RAS-negative substrate. Position BRAF/RAS-negative compartment as one of multiple sites where this gap is most clinically consequential.
- Retain ATA 2015 + ATA 2025 citation as the molecular-modifier gap evidence (BRAF V600E sole molecular modifier), but anchor the gap on "axis missing", not "compartment missing".
- Retain Bethesda III/IV mention as a parallel diagnostic gap, but optional — not central to the differentiation-axis hook.

### Citation checks

- ✅ Haugen 2016 ATA 2015 (PMID 26462967, PMC4739132) — confirmed
- ⚠ Ringel 2025 ATA 2025 (PMID 40844370) — full author list verify (BibTeX note)
- ⚠ SEER incidence claim — placeholder, alternate cites: Davies & Welch 2014 *JAMA*, Lim et al. 2017 *JAMA Oncol*
- ⚠ Cibas 2017 Bethesda — verify exact wording for 15-30% (PMID 29091573)
- ⚠ "5-20% recurrence / ~10% distant met" — verify specific Haugen 2016 source location

### Paper 2/3/4 overlap

- ❌ No HT/GD/ICI/Pan-Asian/HLA-II language present — clean of Paper 2/3/4 contamination
- "BRAF/RAS-negative" framing used throughout — does not overlap Paper 2/3/4 (Paper 1 territory)

---

## 1.2 Existing molecular framework — verdict: MODIFY

### Risky phrases

| Phrase (exact) | Why risky | Strategy rule |
|---|---|---|
| "...a 16-gene thyroid differentiation core (TDS-core) capturing canonical RAI uptake biology — **including SLC5A5/NIS, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, and DIO1**" | Enumerates the 8 panel genes inside 1.2; previews the panel before motivation/aim is established. Frames panel as discovery premise rather than readout instrument | Rule 2 |
| "The **BRAF-RAS Score (BRS), originally derived from 273 transcripts capturing this axis (Cancer Genome Atlas Research Network, 2014; Yoo et al., 2016)**" | Attribution conflates TCGA 2014 origin and Yoo 2016 elaboration — verify whether 273-gene BRS is from TCGA 2014 alone or jointly | Cite verification |
| "**approximately 23% of TCGA-THCA tumors and up to 37% of Korean cohorts harbor neither BRAF V600E nor RAS hotspot mutations**" | Centers BRAF/RAS-negative scope as the framework gap; per strategy this is supplement, not headline | Rule 3 (mild) |

### Safe replacement concept (no prose)

- **Drop the 8-gene enumeration** from 1.2 — defer to 1.4 (Aim) where panel is applied, or to Methods. Keep TDS-core 16-gene reference (Yoo 2016) without listing panel members.
- Frame TCGA 2014 + Yoo 2016 BRS framework as the **driver-axis reference** that this paper extends with an **orthogonal differentiation axis** (sets up Rule 1 lead).
- BRAF/RAS-negative gap mention OK as one example of where the existing framework leaves clinical questions open; do not anchor section on it.

### Citation checks

- ✅ TCGA 2014 *Cell* (PMID 25417114) — confirmed
- ✅ Yoo 2016 *PLoS Genet* (PMID 27494611) — confirmed
- ⚠ Liu et al. 2017 — full citation TBD (BibTeX note)
- ⚠ BRS-273 attribution — verify TCGA 2014 vs Yoo 2016 origin

### Paper 2/3/4 overlap

- ❌ No Paper 2/3/4 contamination
- Korean cohort mention here is generic (37% BRAF/RAS-neg) — neutral, no HT-specific language

---

## 1.3 Dark matter / advanced-disease continuum — verdict: MAJOR REFRAME

### Risky phrases

| Phrase (exact) | Why risky | Strategy rule |
|---|---|---|
| "**Whether the dedifferentiation phenotype Landa described in advanced disease has an upstream signature within the primary BRAF/RAS-negative PTC compartment has not been systematically tested.**" | Frames Paper 1 as **proving a PTC→PDTC→ATC continuum** (asking whether ATC dedifferentiation has a PTC-stage upstream signature). Strategy explicitly says GSE76039 is replication of advanced disease, NOT monotonic continuum proof. Landa 2016 cannot bear continuum-proof framing in this paper. | Rule 8 |
| "Such an upstream marker, if it existed, would provide both **a mechanistic axis for the dark matter**" | "Mechanistic axis" overstates — strategy positions panel as **compact readout**, not mechanistic axis. "Mechanism" is what the 8-gene reads off, not what it defines. | Rule 2 |
| "...and an early-stage candidate biomarker for **fusion-targeted therapy**..." | Forward-promises fusion-targeted therapy in Intro; strategy says fusion claims NOT headline (fusion calls incomplete). Forward-promise in Intro = headline elevation. | Rule 5 |
| "...and **epigenetic-targeted RAI re-induction strategies**." | Forward-promises epigenetic-targeted RAI re-induction; strategy says methylation NOT demonstrated. "Epigenetic-targeted" implies demonstrated methylation mechanism. | Rule 6 |
| "thyroid differentiation transcripts (TG, TSHR, TPO, PAX8, SLC26A4, DIO1, DUOX2) are profoundly suppressed in ATC" | Landa 2016 quote — verify exact wording matches Landa 2016 (audit memory `v17_landa2016_cite_save.md` confirms 7-gene list, but verify quotation accuracy) | Cite verification |
| "**The dark matter compartment is enriched in East Asian populations: in TCGA-THCA (predominantly North American/European), 28.4% of primary tumors are BRAF/RAS-negative, rising to 37.8% in Korean cohorts**" | Centers BRAF/RAS-negative scope and Korean enrichment — Korean cohort framing here belongs to Paper 2 main thesis (Hashimoto-overlap context), not Paper 1 lead | Rule 3 + Paper 2 boundary |

### Safe replacement concept (no prose)

- **MAJOR REFRAME**: Pivot the section's central question from **"upstream of dedifferentiation continuum"** to **"compact readout for the differentiation axis itself within differentiated thyroid carcinoma"**.
- Position Landa 2016 as **advanced-disease reference cohort** providing comparator for differentiation-transcript suppression — NOT as the trajectory's downstream that this paper extends upstream.
- **Drop forward-looking promise language** ("fusion-targeted therapy", "epigenetic-targeted RAI re-induction") — relocate to Discussion as candidate forward implications with rule-compliant framing ("consistent with DNMT-associated repression"; fusion enrichment as supportive observation, not headline claim).
- **GSE76039 introduction**: if used in Results, frame here as additional advanced-disease cohort (not as part of monotonic continuum). Currently 1.3 doesn't mention GSE76039 — decide whether to mention here or defer to Methods/Results.
- Korean enrichment 37.8% mention: tone down to single-line generalization observation, NOT thesis-defining; keeps this Paper 1 territory clean from Paper 2 substrate.

### Citation checks

- ⚠ Xing dark matter cite — BibTeX has `Xing2014` keyed to Xing M. *Nat Rev Cancer.* 2013;13(3):184-199. Strategy memo treats Xing as 2014 dark matter origin; verify whether dark matter concept is from this *Nat Rev Cancer* 2013, Xing 2014 *NEJM*, or Xing 2015 *JCO*
- ✅ Landa 2016 *JCI* (PMID 26878173) — confirmed; verify quotation of differentiation transcript list
- ⚠ Wang et al. 2024 — Wang YL *Endocr Connect* 13(11):e240301 (PMID 39235852) — verify (BibTeX note flags D3-P4 audit context)

### Paper 2/3/4 overlap

- ⚠ "37.8% in Korean cohorts" — Korean cohort substrate enrichment risks overlap with Paper 2 (Hashimoto-overlap PTC HLA-II axis). Tone down to neutral statement, do NOT use as Paper 1 thesis anchor.
- ⚠ "epigenetic-targeted RAI re-induction" — Paper 1 territory by strategy, but per Rule 6 tone down to "consistent with DNMT-associated repression"; do NOT promise as Intro headline.
- ❌ No HT/GD/ICI/Pan-Asian/HLA-II language present — clean

---

## 1.4 Aim and preview — verdict: REWRITE REQUIRED

This section has **multiple strategy violations** in headline claims. Re-architecture needed before line-level edit.

### Risky phrases

| Phrase (exact) | Why risky | Strategy rule |
|---|---|---|
| "**Here, we apply an 8-gene RAI-responsiveness panel** — independently selected from canonical thyroid differentiation biology (Yoo et al., 2016) before access to the Landa 2016 ATC-silenced gene list" | "We apply" framing is OK in principle (matches Rule 2 readout positioning), BUT panel is introduced **before** the differentiation axis is articulated — sequence inverts Rule 1 (axis lead, panel as how-we-measure-it). The "before access to Landa" reverse-causality clause is critical and must remain (defends 5/8 overlap) but should appear after axis lead. | Rule 1 sequence |
| "...across **TCGA-THCA (n=504), MSK-IMPACT thyroid (n=117), Korean cohorts (n=874), and external single-cell datasets**." | (a) **MSK-IMPACT thyroid (n=117)** ambiguously labels Landa 2016 advanced-disease cohort (84 PDTC + 33 ATC). Without "advanced disease" qualifier, reader may assume primary tumors. (b) **"external single-cell datasets"** is vague — name them or omit cohort enumeration in Aim (defer to Methods). (c) **"Korean cohorts (n=874)"** preserves substrate scope but does not signal Paper 2 boundary; OK at this level of generality. | Rule 8 (MSK = advanced); cohort precision |
| "We show that this panel **resolves the BRAF/RAS-negative compartment** into a DM1/DM2 axis" | Centers BRAF/RAS-negative as primary scope of resolution. Strategy positions BRAF/RAS-negative as motivation, not headline — DM1/DM2 axis should be presented as resolving the differentiation continuum across all PTC, with BRAF/RAS-negative compartment as one application. | Rule 3 |
| "...**orthogonal to canonical driver classification, stable to candidate-pool restriction, and neutral to driver-transcript abundance.**" | This is the differentiation-axis defense set — strong, ON-strategy. Keep this language. | ✅ Rule 1 compliant |
| "Within this compartment, **DM1 captures most tyrosine-kinase-fusion-positive tumors, including 81.8% of TCGA RET-fusion-positive cases**" | ★ Headline claim of fusion enrichment with specific 81.8% RET-fusion number. Strategy Rule 5 explicitly says fusion claims NOT headline because fusion calls incomplete. The 81.8% number specifically is one of the listed special-attention phrases. | Rule 5 (★) |
| "...while also **harboring fusion-independent promoter hypermethylation of thyroid differentiation genes**." | ★ Headline claim of "promoter hypermethylation" as a demonstrated mechanism. Strategy Rule 6 says methylation NOT demonstrated — only "consistent with DNMT-associated repression" framing allowed. "Promoter hypermethylation" is a specific molecular claim that exceeds the evidence the strategy permits. | Rule 6 (★) |
| "These findings define DM1 as a **fusion-driven, epigenetically silenced dark-matter subtype**" | ★★★ Combines fusion-driven (Rule 5 violation) + epigenetically silenced (Rule 6 violation) + dark-matter headline (Rule 3 violation) into a single subtype-defining phrase. This is the headline subtype definition that strategy memo most directly contradicts. | Rules 3 + 5 + 6 (★★★) |
| "...and support a clinically interpretable framework for **reflex fusion testing and prospective evaluation of epigenetic-targeted RAI re-induction**." | Forward-loading two clinical actionability claims (fusion reflex + epigenetic re-induction) in Aim. Both are downstream applications, not Aim-level claims. Reflex fusion testing per Rule 5 is supportive, not headline. Epigenetic-targeted re-induction per Rule 6 lacks demonstrated mechanism support. | Rules 5 + 6 |

### Safe replacement concept (no prose) — re-architecture spec

Section 1.4 needs **structural rewrite**, not line-level edit. Suggested new shape:

1. **Sentence 1 (axis lead)**: Articulate the question Paper 1 answers — does a compact RAI-lineage transcriptional readout reveal a clinically meaningful sub-stratification of differentiated thyroid carcinoma orthogonal to canonical BRAF/RAS driver classification? (Rule 1 lead)
2. **Sentence 2 (panel as instrument)**: Introduce 8-gene panel as the compact readout instrument, derived independently from canonical RAI biology (Yoo 2016) prior to Landa 2016 exposure — preserves reverse-causality defense. (Rule 2)
3. **Sentence 3 (cohorts)**: List discovery + validation cohorts with explicit advanced-disease labeling for Landa 2016 (n=117) and external sc replication (Pu 2021 named, Lu 2023 named). Include GSE76039 if used in Results — labeled as **external advanced-disease replication**. Korean cohort mention general, no HT specifics. (Rule 8)
4. **Sentence 4 (axis findings)**: Report the differentiation-axis defense set — orthogonal to BRAF/RAS, stable to panel choice, driver-mRNA-neutral. (Rule 1, ON-strategy claim)
5. **Sentence 5 (sub-stratification observations)**: Report DM1 properties as **observed enrichments / consistencies** — DM1 is enriched for actionable molecular alterations consistent with broader dedifferentiation programs. **Do not headline 81.8% RET capture, "fusion-driven", or "promoter hypermethylation"**.
6. **Sentence 6 (forward implication, restrained)**: Frame applications as **forward implications evaluated in Discussion** — clinical sub-stratification framework; supportive observations consistent with DNMT-associated repression of differentiation transcripts. **Do not promise** "epigenetic-targeted RAI re-induction" at Aim level.

### Special phrase resolution

| Special phrase | Status |
|---|---|
| "Here, we apply an 8-gene..." | KEEP framing ("apply" matches readout role) but reorder so axis articulation precedes panel introduction |
| "BRAF/RAS-negative compartment" | TONE DOWN to motivation/supplement throughout 1.4 (currently 2 occurrences in 1.4 + 4+ across other sections) |
| "MSK-IMPACT n=117" | RELABEL with explicit "advanced-disease" qualifier (Landa 2016 cohort is 84 PDTC + 33 ATC, not primary tumors) |
| "external single-cell datasets" | NAME explicitly (Pu 2021 + Lu 2023) or DEFER to Methods — vague phrasing weak |
| "81.8% RET-fusion" | REMOVE from Aim headline; relocate to Results 2.3 |
| "promoter hypermethylation" | REPLACE with "transcriptional repression consistent with DNMT-associated mechanisms" (Rule 6) |
| "fusion-driven, epigenetically silenced dark-matter subtype" | REWRITE — the subtype definition triggers 3 strategy rules. Suggested concept: "a clinically aggressive sub-stratum of differentiated thyroid carcinoma defined by a compact RAI-lineage readout and enriched for actionable molecular alterations" |
| "ATA 2025" | (1.1) KEEP citation (Ringel 2025 PMID 40844370), verify full author list |
| "Bethesda III/IV" | (1.1) KEEP as parallel diagnostic gap mention; verify Cibas 2017 source |

### Citation checks

- ✅ Yoo 2016 — confirmed
- ✅ Landa 2016 — confirmed (verify advanced-disease labeling consistency)
- ✅ Pu 2021 *Nat Commun* (PMID 34663816) — confirmed (verify author list spelling)
- ⚠ Lu 2023 — full citation TBD (BibTeX note: GSE193581 first author + journal)
- ⚠ GSE76039 — if introduced in 1.4 or later, citation TBD

### Paper 2/3/4 overlap

- ❌ Aim itself does not mention HT / GD / ICI / Pan-Asian / HLA-II — clean of cross-paper substrate language
- ⚠ "fusion-driven" + "epigenetically silenced" subtype framing risks overstating Paper 1 scope; tone-down per Rules 5 + 6 also tightens cross-paper hygiene (Paper 2 and Paper 4 will have their own mechanism claims)
- ⚠ "Korean cohorts (n=874)" — neutral substrate mention OK at this generality; do NOT attach HT-specific framing in Aim (Paper 2 territory)

---

## Cross-section summary

| Section | Verdict | Headline issue |
|---|---|---|
| 1.1 Clinical context | MODIFY | Hook centers BRAF/RAS-negative compartment (Rule 3); 4 cite verifications needed |
| 1.2 Existing molecular framework | MODIFY | 8-gene enumeration premature (Rule 2); BRS-273 attribution verify |
| 1.3 Dark matter / advanced disease | **MAJOR REFRAME** | Frames paper as PTC→ATC continuum proof (Rule 8); forward-promises fusion + epigenetic therapies (Rules 5 + 6); centers Korean BRAF/RAS-neg enrichment as thesis |
| 1.4 Aim and preview | **REWRITE REQUIRED** | Headlines 81.8% RET (Rule 5) + promoter hypermethylation (Rule 6) + "fusion-driven, epigenetically silenced dark-matter subtype" (Rules 3 + 5 + 6); MSK-IMPACT advanced-disease ambiguity (Rule 8) |

## Citation verification queue (consolidated)

| # | Citation | Status |
|---|---|---|
| 1 | SEER incidence (1.1 placeholder) | Resolve to Davies & Welch 2014 *JAMA* OR Lim et al. 2017 *JAMA Oncol* OR retain SEER URL |
| 2 | Cibas 2017 Bethesda 15-30% (1.1 placeholder) | PMID 29091573 verify exact wording |
| 3 | Haugen 2016 "5-20% recurrence / ~10% distant met" (1.1) | Verify Haugen 2016 source location for these specific figures |
| 4 | Ringel 2025 ATA 2025 (1.1) | PMID 40844370, full author list verify |
| 5 | TCGA 2014 vs Yoo 2016 BRS-273 attribution (1.2) | Verify which paper defines 273-gene BRS |
| 6 | Liu et al. 2017 (1.2 + 1.3) | Full citation TBD |
| 7 | Xing dark matter cite (1.3) | Disambiguate Xing 2013 *Nat Rev Cancer* vs Xing 2014 *NEJM* vs Xing 2015 *JCO* |
| 8 | Wang et al. 2024 (1.3) | Wang YL *Endocr Connect* 13(11):e240301 (PMID 39235852) verify |
| 9 | Landa 2016 differentiation transcript quotation (1.3) | Verify exact wording matches Landa 2016 |
| 10 | Pu 2021 author list (1.4) | PMID 34663816 verify spelling |
| 11 | Lu 2023 (1.4) | Full citation TBD; GSE193581 first author + journal verify |
| 12 | GSE76039 (if introduced in 1.x) | Citation TBD |

## Cross-paper boundary verification

- **Paper 2 (HT-overlap PTC) overlap risk in Intro**: Korean enrichment 37.8% framing in 1.3 is the highest-risk substrate-scope overlap. Tone down to neutral observation, do NOT anchor as Paper 1 thesis. No HT-specific terms in current draft (clean).
- **Paper 3 (ICI vulnerability) overlap risk**: None found in current Intro draft.
- **Paper 4 (Korean GD HLA) overlap risk**: None found. No GD/Graves/TSAb language present.
- **Forbidden words scan**: Cross-reference with `_terminology_scan_hits.tsv` — Intro file (`03_introduction.md`) currently has 0 actionable hits for forbidden Paper 2/3/4 terms after user's 2026-05-04 cleanup.

## Marathon mode compliance (audit-only deliverable)

- ✅ No prose rewrite
- ✅ No new introduction generated
- ✅ Manuscript not modified
- ✅ Section-by-section revision plan only
- ✅ Citation verification deferred to user's prep day reading work (Landa, Yoo, ATA already on user's 5/4 read list per `_prep_04_reading_urls.md`)

# END AUDIT — revision plan only. Manuscript edits await user keyboard (voice-protected) + cite verification (5/4 prep).
