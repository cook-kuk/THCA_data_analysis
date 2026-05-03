# Paper 1 — BibTeX + Citation Cross-Check Audit (T2)

**Date:** 2026-05-04 (marathon scaffolding/infra; read-only audit)
**Scope:** `2026_05_03_references.bib` (29 entries) ↔ all manuscript prose files (Outline / Methods / Results / Discussion / Captions / QA / STAR).
**Mode:** Audit only. **No file modified.**
**Voice-risk:** ZERO.

---

## 1. Bib summary

- **Total `@article` entries:** 29 (verified via `grep -c "^@" references.bib`)
- **Header claim (line 371):** "~22 papers + tools" — ⚠ inconsistent with actual count.
- **Cell Press style target:** numerical references.

---

## 2. Bib entry inventory

| # | Citation key | Author | Year | Journal | Status |
|---|---|---|---|---|---|
| 1 | `Cancer2014TCGA` | TCGA Research Network | 2014 | Cell | ✓ used in Outline ¶2 |
| 2 | `Yoo2016SNUGMI` | Yoo SK et al. | 2016 | PLOS Genet | ✓ used in M2, Q2, Outline |
| 3 | `Lee2024GSE213647` | Lee Y, others | 2024 | GEO | ⚠ author placeholder "Y and others"; needs full author list |
| 4 | `Lim2025GSE286332` | Lim DW, Kim SM | 2025 | GEO; PMID 41113708 | ✓ used widely |
| 5 | `Chu2018JMG` | Chu X et al. | 2018 | J Med Genet | ✓ used in M9, Q13, R1, F1 |
| 6 | `Cook2021cookHLA` | Cook S et al. | 2021 | Nat Commun | ✓ used in R1 (Cook et al. 2021) |
| 7 | `Orenbuch2020arcasHLA` | Orenbuch R et al. | 2020 | Bioinformatics | ✓ used in M9, STAR |
| 8 | `Robinson2020IPDIMGT` | Robinson J et al. | 2020 | NAR | ✓ STAR only |
| 9 | `Muzellec2023PyDESeq2` | Muzellec B et al. | 2023 | Bioinformatics | ✓ STAR only; M3 cites tool version not paper |
| 10 | `Fang2023gseapy` | Fang Z et al. | 2023 | Bioinformatics | ✓ STAR only; M3 cites tool version |
| 11 | `Liberzon2015Hallmark` | Liberzon A et al. | 2015 | Cell Syst | ✓ STAR only |
| 12 | `Frankish2021GENCODE` | Frankish A et al. | 2021 | NAR | ✓ STAR only |
| 13 | `Landa2016JCI` | Landa I et al. | 2016 | JCI | ✓ Outline ¶1 + Discussion §3.1 (voice-protected) |
| 14 | `Pan2025NatComm` | Pan Z et al. + Ge M | 2025 | Nat Commun | ✓ Outline ¶2 ("Pan/Ge 2025") |
| 15 | `Riesco2014EJE` | Riesco-Eizaguirre & Santisteban | 2014 | EJE | ✓ M2, Q1 |
| 16 | `Cabrita2020TLS` | Cabrita R et al. | 2020 | Nature | ✓ M7, R5c, Caption F5C |
| 17 | `Pfister2013Bimodality` | Pfister R et al. | 2013 | Front Psychol | ✓ M5 |
| 18 | `Haugen2016ATA` | Haugen BR et al. | 2016 | Thyroid | ✓ Outline ¶1 (ATA 2015 cheatsheet) |
| 19 | `Tuttle2019JCEM` | Tuttle RM, Alzahrani AS | 2019 | JCEM | ⚠ NOT YET CITED in any prose file |
| 20 | `Yi2016KTA` | Yi KH et al. | 2016 | Endocrinol Metab | ⚠ NOT YET CITED |
| 21 | `Chen2024EndocrConnect` | Chen XF et al. + Wang Y | 2024 | Endocr Connect | ⚠ NOT YET CITED in main prose (only in `D3P4_wang2024_audit.md` decision record) |
| 22 | `Pedregosa2011sklearn` | Pedregosa F et al. | 2011 | JMLR | ✓ STAR only |
| 23 | `Virtanen2020scipy` | Virtanen P et al. | 2020 | Nat Methods | ✓ STAR only |
| 24 | `Harris2020numpy` | Harris CR et al. | 2020 | Nature | ✓ STAR only |
| 25 | `Bray2016kallisto` | Bray NL et al. | 2016 | Nat Biotechnol | ✓ STAR only |
| 26 | `Langmead2012bowtie2` | Langmead B, Salzberg SL | 2012 | Nat Methods | ✓ STAR only |
| 27 | `Li2009samtools` | Li H et al. | 2009 | Bioinformatics | ✓ STAR only |
| 28 | `Lee2014TissueAntigens` | Lee KW et al. | 2014 | Tissue Antigens | ✓ Q13 (Korean baseline reference) |
| 29 | `Kim2014KoreanGraves` | Kim TH et al. | 2014 | Korean Endocrinol Soc | ⚠ NOT YET CITED in main prose; memory mentions 56% Korean GD DPB1*05:01 replication |

---

## 3. Cited but not in bib — STATUS: ✓ NONE FOUND

All "Author Year" patterns observed in the prose map to a bib entry:
- "TCGA Cancer Network 2014" → `Cancer2014TCGA`
- "Yoo et al. 2016" / "Yoo 2016" → `Yoo2016SNUGMI`
- "Lee et al. 2024" / "Lee Y" → `Lee2024GSE213647`
- "Lim DW et al. 2025" → `Lim2025GSE286332`
- "Chu et al. 2018" / "Chu 2018" → `Chu2018JMG`
- "Cook et al. 2021" → `Cook2021cookHLA`
- "Orenbuch et al. 2020" → `Orenbuch2020arcasHLA`
- "Cabrita et al. 2020" / "Cabrita 2020" → `Cabrita2020TLS`
- "Pfister et al. 2013" → `Pfister2013Bimodality`
- "Riesco-Eizaguirre & Santisteban 2014" → `Riesco2014EJE`
- "Landa 2016" / "Krishnamoorthy/Landa 2016" → `Landa2016JCI`
- "Pan/Ge 2025" / "Pan 2025" → `Pan2025NatComm`
- "Lee et al. 2014 Tissue Antigens" → `Lee2014TissueAntigens`

**Verdict:** No orphan citation found. Bib coverage is complete for currently-drafted prose.

---

## 4. ⚠ In bib but not yet cited — 4 entries

| Bib key | Reason for inclusion | Likely use |
|---|---|---|
| `Tuttle2019JCEM` | ATA 2015 + risk-stratification context | Discussion §3.1 / §3.2 (voice-protected) — likely user voice will cite |
| `Yi2016KTA` | Korean Thyroid Association 2016 guidelines | Discussion §3.2 (ATA-vs-KTA comparison) — likely user voice will cite |
| `Chen2024EndocrConnect` | Wang/Chen 2024 Shanghai n=2,844 mutation comparator | Discussion §3.3 (Asian mutation context) — D3P4 audit memo decided to use as comparator only |
| `Kim2014KoreanGraves` | Korean GD DPB1*05:01 ~36% baseline | Q13 + Discussion §3.3 (Pan-Asian susceptibility) — could replace or augment Lee2014TissueAntigens reference |

**Reconciliation candidate (no decision):**

- **(a)** Keep all 4. They are domain-aligned and likely to be cited in voice-protected Discussion. Acceptable to keep in bib as reserve.
- **(b)** Drop unused entries before submission to keep bib lean. Move to `references_reserve.bib`.
- **(c)** Decide per-entry during voice-protected Discussion drafting (user keyboard).

Recommend **(c)** — defer to user voice draft. Track in checklist.

---

## 5. ⚠ Author placeholder in `Lee2024GSE213647`

```
@article{Lee2024GSE213647,
  author = {Lee, Y and others},   ⚠ placeholder
  ...
}
```

**Action required:** Pull full author list from GSE213647 GEO record (Lee Y et al. 2024) or PubMed, before bioRxiv 6/8 freeze. Same may apply to `Kim2014KoreanGraves` and `Lee2014TissueAntigens` (both have minimal "and others" patterns).

---

## 6. ⚠ Title typo / formatting risk in bib

Spot-check identified:

- `Cancer2014TCGA` — `author = {{Cancer Genome Atlas Research Network}}` — double braces correct for institutional author. ✓
- `Lim2025GSE286332` — title is "Dysregulation of Vitamin D and Its Signaling in Hashimoto's Thyroiditis in Korean Population" — confirm with PubMed PMID 41113708. ⚠ Vitamin-D framing may not match our usage of GSE286332 as "PTC vs PTC+HT" — verify the actual published title before submission.
- `Pan2025NatComm` — `pages = {s41467-025-58910-3}` — pages format unusual; should be article number e.g., `1234` or `e58910`. Verify against Nat Comm convention.
- `Chu2018JMG` — `note` contains "★ Citation correction: previously misattributed as 'Chen 2018'". The ★ Unicode character may break some bibtex pipelines. ⚠ Recommend ASCII alternative ("(NOTE)").

---

## 7. Cross-paper citation discipline (Paper 2 / Paper 3 / Paper 4 boundary)

Paper 1 prose **must not** cite:
- Paper 2 (HT-overlap PTC pathology) — Paper 1 may use Paper 2's GSE286332 *data* but not its analytic claims.
- Paper 3 (ICI dark matter) — frozen, not citable from Paper 1.
- Paper 4 (GD HLA backlog) — gated, not yet citable.

**Audit verdict:** ✓ no cross-paper citation found in current Paper 1 prose. Paper 1 stays self-contained.

---

## 8. Voice-protected note

Discussion §3.1 (Krishnamoorthy/Landa 2016 framing) is voice-protected. When user drafts §3.1, they may need to add `@Krishnamoorthy2025NatComm` if Krishnamoorthy 2025 is being cited (currently absent from bib but flagged in memory `v17_landa2016_cite_save.md` as a misattribution-corrected reference). **Possible bib gap.**

**Action:** Before §3.1 drafting, decide:
- **(a)** Cite Landa 2016 only (current bib state), per `v17_landa2016_cite_save.md` correction, or
- **(b)** Cite both Landa 2016 + Krishnamoorthy 2025 — then add Krishnamoorthy2025 entry to bib.

---

## 9. Summary of audit findings

| # | Severity | Issue | Action |
|---|---|---|---|
| 1 | **MEDIUM** | `Lee2024GSE213647` author placeholder | Pull full author list from GEO/PubMed |
| 2 | **MEDIUM** | `Lim2025GSE286332` title vs paper usage mismatch (Vitamin D framing) | Verify against PMID 41113708 |
| 3 | LOW | Bib header "~22" vs actual 29 | Update header (cosmetic) |
| 4 | LOW | 4 bib entries (Tuttle2019, Yi2016, Chen2024, Kim2014) not yet cited | Defer to Discussion drafting |
| 5 | LOW | `Pan2025NatComm` pages format unusual | Verify Nat Comm convention |
| 6 | LOW | ★ Unicode in `Chu2018JMG` note | Replace with ASCII |
| 7 | LOW | Possible Krishnamoorthy 2025 bib gap (depends on §3.1 voice) | Decide at §3.1 drafting time |

**No orphan citation, no missing-bib citation found.** Bib is in good shape for current prose state.

---

## 10. Recommended next action

1. **Pre-bioRxiv (W6 6/8):** Fix MEDIUM #1 + #2 (author list + title verification).
2. **Cosmetic pass:** Fix LOW #3 + #6 in same diff.
3. **Discussion drafting time:** Decide #4 + #7 with user voice.
4. **Cell Rep Med formal submission (W7 6/14):** Re-audit with this same checklist.

---

Track A completed. No marathon violation. Bib + citation audit is read-only — no manuscript file modified.
