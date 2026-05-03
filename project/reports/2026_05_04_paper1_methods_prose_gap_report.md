# Paper 1 — Methods Prose vs Scaffold Gap Report (T4)

**Date:** 2026-05-04 (marathon scaffolding/infra; read-only audit)
**Scope:** Compare `2026_05_03_methods_M1_M11_scaffold.md` (124 lines) vs `2026_05_03_methods_M1_M11_prose.md` (72 lines) for section coverage gaps.
**Mode:** Audit only. **No file modified.** Originally proposed as "gap-fill edit"; audit found no real coverage gap, so converted to report-only.
**Voice-risk:** ZERO (no prose generated).

---

## 1. Why this was converted from edit to audit

The candidate task T4 originally proposed: "Methods prose gap-fill — generate fact-only prose for missing scaffold sections." After comparing both files section-by-section, **no missing section identified**. Both files cover M1, M2, M3, M4, M5, M5b, M6, M7, M8, M9, M10, M11 with substantively equivalent content. The line-count difference (124 vs 72) is driven by:

- Scaffold has **header metadata + voice-protected reminder section + quality-check at bottom** (~30 lines of non-Methods scaffolding).
- Scaffold uses denser bullet formatting + extra blank lines.
- Prose uses tighter paragraph form.

Generating "missing prose" would have produced redundant content, with non-zero verb-tone leak risk. **Edit cancelled. Report-only output below.**

---

## 2. Section-by-section coverage matrix

| Section | Scaffold | Prose | Coverage verdict |
|---|---|---|---|
| M1 Cohort Assembly | line 9–11 | line 9–10 | ✓ equivalent — both list 5 cohorts with n + accession |
| M2 TIERA67 candidate pool | line 15–17 | line 12–14 | ✓ equivalent — both list 7 categories + 8-gene panel members |
| M3 DEG + GSEA | line 21–25 | line 16–20 | ✓ equivalent — PyDESeq2 v0.5.4 + gseapy 1.1.13 + parameters |
| M4 DM1/DM2 + classifier | line 29–31 | line 22–24 | ✓ equivalent — KMeans + LogReg + AUC=0.962 |
| M5 Hashimoto signature transfer | line 35–37 | line 26–28 | ✓ equivalent — 150 up + 50 dn signature, 4 thresholds |
| M5b Confounder residualization | line 39–43 | line 30–32 | ✓ equivalent — Stromal + immune-proxy residualization |
| M6 Mediation Baron-Kenny | line 47–49 | line 34–36 | ✓ equivalent — 3-equation framework + 5,000 bootstrap |
| M7 BCR + TLS | line 53–55 | line 38–40 | ✓ equivalent — IGH/IGK/IGL extraction + Cabrita TLS |
| M8 DM1 sub-cluster | line 59–65 | line 42–48 | ✓ equivalent — KMeans k=2 + Welch t + sub-A/B labels |
| M9 Pan-Asian forest meta | line 69–77 | line 50–58 | ✓ equivalent — arcasHLA + DerSimonian-Laird + 4-scenario sensitivity |
| M10 Pan-genome MAD | line 81–85 | line 60–64 | ✓ equivalent — top-N MAD + KMeans + ARI/NMI |
| M11 Stats + reproducibility | line 89–95 | line 66–72 | ✓ equivalent — Cohen d + Wilson CI + BH-FDR + Python 3.12 environment |
| Voice-protected reminder | line 99–111 | (absent) | scaffold-only meta |
| Quality check ✅ | line 115–125 | (absent) | scaffold-only meta |

**Gaps:** none.

---

## 3. Substantive content drift between scaffold and prose — STATUS: minor

Spot-check identified 3 sentences appearing in scaffold but slightly trimmed in prose. None are factually consequential.

| Topic | Scaffold | Prose | Drift |
|---|---|---|---|
| M5b OR=0.289 statement | "(OR = 0.289, p = 8×10⁻⁹), confirming that the autoimmune-PTC axis is not a generic immune-infiltration artifact" | (line 32 of prose) "yielded OR = 0.289 (Fisher p = 8×10⁻⁹), confirming that the autoimmune-PTC axis is not a generic immune-infiltration artifact" | ✓ same content |
| M6 single-predictor dominance | "with HLA-II as the dominant single-predictor" | "with HLA-II as dominant single-predictor (Suppl Table S5/S5b)" | ✓ prose adds Suppl Table cite — improvement |
| M9 Suppl Figure S6 cite | "Sensitivity analyses recomputed the meta-analysis under four scenarios" | "(Suppl Figure S6)" added | ✓ prose adds Suppl Figure cite — improvement |
| M11 GitHub Actions CI | (absent) | "Reproducibility smoke tests covering all five pillars are provided in `tests/test_signature_score.py` and verified continuously via GitHub Actions CI" | ✓ prose adds CI mention — improvement |

**Verdict:** Prose is **slightly more polished** than scaffold (cross-references inserted, CI mention added). No regression.

---

## 4. ⚠ MINOR FINDING — Methods prose missing an explicit cite to references.bib

**Observation:** Methods prose mentions tool names + versions (PyDESeq2 v0.5.4, gseapy v1.1.13, arcasHLA v0.6.0, etc.) but does not insert author-year citations inline. STAR Methods Key Resources Table covers tool citations.

**Cell Press convention:** numerical references in main text body. Methods may rely on STAR Methods table for tool citations. Acceptable.

**Reconciliation candidates (no decision):**
- **(a)** Keep current state — STAR table carries citations.
- **(b)** Add inline "(Muzellec et al. 2023)", "(Fang et al. 2023)", "(Orenbuch et al. 2020)" citations in Methods M3 and M9. Lengthens prose by ~5 sentences.

Recommend **(a)** — STAR Methods table is the conventional location.

---

## 5. ⚠ MINOR FINDING — Verb-tone consistency

Spot-check of Results section verbs (R1–R5):
- "We assembled" / "We analyzed" / "We tested" — used uniformly ✓
- "Demonstrating" / "indicating" / "consistent with" / "confirming" / "extending" — used in conclusion clauses, varied
- Borderline subjective: "★ exceptional in this cohort" (Caption F2D line 66) — could read as voice-touched. ⚠

**Reconciliation:** ★ markers in captions are convention markers, not subjective claims. Borderline. User may want to remove or replace with neutral language during voice-protected polish.

**This audit does not modify any caption.** Flag only.

---

## 6. Summary of gap report

| # | Severity | Finding |
|---|---|---|
| 1 | **NONE** | Methods prose covers all M1–M11 + M5b sections — no gap |
| 2 | LOW | Prose is *slightly more polished* than scaffold (added Suppl Table/Figure cites + CI mention) |
| 3 | LOW | Methods relies on STAR Methods table for tool citations — acceptable Cell Press convention |
| 4 | LOW | "★ exceptional" in F2D caption is borderline subjective — flag for user voice polish |

**No prose generation needed. T4 closed as audit.**

---

## 7. Recommended next action

1. **Do not modify Methods prose** — it is complete and already slightly improved over scaffold.
2. **Optionally remove "★ exceptional"** from F2D caption during voice-protected polish (single-word edit).
3. **Optionally add inline tool citations** in Methods M3/M9 if Cell Press reviewer requests — currently not needed.

---

Track A completed. No marathon violation. Methods prose audit is read-only — no manuscript file modified.
