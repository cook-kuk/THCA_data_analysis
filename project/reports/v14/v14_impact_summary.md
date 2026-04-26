# v14 — Impact summary (paper-quality)

_Prepared 2026-04-25. One page. For internal use + to drive paper edits._

## 1. Three core findings

**1. TROP2 high-expression CCLE thyroid lines are all BRAF V600E-mutant.**
The three highest-TACSTD2 lines in the CCLE thyroid panel — B-CPAP (z=+2.42), BHT-101 (z=+1.94), 8505C (z=+0.36) — all carry BRAF V600E; the seven lowest-TACSTD2 lines are wild-type, KRAS, or NRAS. This is the first cell-line-level corroboration of the v13 hypothesis that TROP2 up-regulation in PTC is BRAF-axis-driven.

**2. SN-38 class preferentially kills BRAF-like thyroid lines in PRISM.**
In the PRISM Repurposing 19Q4 primary screen, topotecan and irinotecan (the SN-38 payload class that is delivered by sacituzumab govitecan) produce ~1.0 log-fold-change greater killing in BRAF-like CCLE thyroid lines than in RAS-like lines (topotecan ΔLFC = −1.10; irinotecan ΔLFC = −0.96). This is a novel, orthogonal validation of the v13 TROP2-ADC repurposing rationale and is independent of TROP2-directed delivery.

**3. Three of the five top BRAF-selective compounds are v13 mechanism anchors.**
Of the five compounds with BRAF-selective killing (ΔLFC < −0.5), three — topotecan (TACSTD2 payload), simvastatin (LDLR axis), atorvastatin (LDLR axis) — are direct v13 mechanism anchors, and a fourth (kaempferol) is a v13 CYP1B1 flavonoid. The v13 prioritisation is therefore concordant with an external, unsupervised drug-screen read-out.

## 2. How to integrate into the Bioinformatics paper — Results §5.Y

- Use `reports/v14/v14_ccle_section.md` verbatim as the first draft of Results §5.Y, immediately after v13 (§5.X).
- Open with the **honest caveat** (BRS52 66.7 % on cell lines; preserves v8.1's negative verdict) so the two positive read-outs are framed as orthogonal, not as signature replication.
- Add one figure: CCLE TACSTD2 z-score vs mutation label (bar or dot plot, 13 lines) — data already in `ccle_8target_wide.tsv`.
- Add one table: BRAF-selective compounds with ΔLFC (5 rows) — data already in `ccle_drug_sensitivity_by_subtype.tsv`.
- In the discussion, connect §5.Y.4 topotecan/irinotecan finding back to the §5.X TROP2 story with one sentence: *"the PRISM read-out cannot distinguish TROP2-directed delivery from tumour-intrinsic SN-38 sensitivity in BRAF-like cells — this is an explicit open question for xenograft follow-up."*

## 3. How to update the TROP2 standalone draft

`reports/v13/v13_trop2_standalone_draft.md` currently argues the BRAF-like × TACSTD2-high two-gate biomarker on the basis of transcriptomic stratification alone. v14 adds cell-line-level corroboration:

- **Insert into §2 (TROP2 biology):** one sentence noting that CCLE thyroid lines recapitulate the BRAF-axis TROP2 up-regulation (top three lines all V600E), despite CCLE BRS52 classification being imperfect.
- **Insert into §5 (Proposed correlative design):** add a sentence that the BRAF-like / TROP2-high signature is *doubly motivated* — by tumour-level RNA stratification AND by PRISM-level SN-38 sensitivity — so the correlative overlay tests two non-redundant mechanistic hypotheses.
- **Add as §6.5 (new subsection "Payload vs delivery question"):** flag that the PRISM signal is payload-class-level, and a xenograft comparison of sacituzumab govitecan (TROP2 + SN-38) vs. free topotecan in BRAF-like vs. RAS-like PTC models would directly quantify the targeted-delivery contribution. This is a wet-lab v15 follow-up, not a computational result.

## 4. Open questions for v15 / wet-lab

- **Payload vs delivery.** Does TROP2-directed ADC add value over free SN-38 in BRAF-like PTC, or is the BRAF-like sensitivity a cell-intrinsic SN-38 vulnerability that any Topo-I inhibitor would exploit? Direct xenograft comparison needed (SG vs free topotecan vs vehicle, in BCPAP / 8505C / FTC-133 xenografts).
- **LDLR mechanism.** Simvastatin and atorvastatin are strongly BRAF-selective in PRISM but the driver is unclear — statins up-regulate LDLR via SREBP (their intended mechanism) or hit a BRAF-pathway-linked metabolic vulnerability? MEK-i combination test and cholesterol-rescue experiments in paired BRAF-like vs RAS-like lines would resolve.
- **CBD / dronabinol gap.** PRISM 19Q4 does not cover cannabidiol (our flagship v7 CYP1B1 hit) or dronabinol (v13 Tier-A). A scheduled one-shot agent (2026-05-23) will re-check newer PRISM releases; if still absent, a targeted viability assay (BCPAP / 8505C / CAL-62, 3 doses × 2 compounds) closes the CYP1B1 mechanism arm for the price of one experiment.
- **BRS52 cell-line failure.** BRS52 fails at 66.7 % on CCLE — is this a generic cell-line drift artefact, or does it specifically implicate the BRS52 genes that are not expressed in established PTC lines (e.g., TPO, DIO1, TG)? An ablation with a cell-line-compatible classifier (drop thyroid-differentiation-specific genes) could rescue accuracy and inform what the signature is *actually* reading in primary tumour.
