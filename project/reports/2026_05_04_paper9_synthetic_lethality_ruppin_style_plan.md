# Paper 9 — Thyroid Synthetic-Lethal Vulnerability Map (Ruppin-style, conceptual)

Date: 2026-05-04
Author: Seungho Cook
Status: STRATEGIC DESIGN ONLY. No data download. No analysis. No GPU. No manuscript text.
Scope: post-marathon (after 2026-06-13). Anchor on Paper 1 once it is in print.
Forbidden in this document: copying ISLE/SELECT prose, claims, figures, or signatures. Concepts are referenced; sentences are not.

---

## 1. Executive verdict

GO, conditional on Paper 1 in print.

Paper 9 is a **forward-looking translational paper** that takes the DM1 / RAI-lineage-collapsed state defined in Paper 1 as a **state-of-vulnerability** and asks: "given this state, what genes/drugs kill the tumor selectively, and which compensatory partners must we co-target?" The design is Ruppin-style in spirit (transcriptome-guided synthetic lethality / synthetic rescue, population-level co-occurrence patterns, cross-cohort transferability) but specifically thyroid-tuned and DM1-anchored.

Why this is a real paper and not an extension memo:
- Paper 1 defines a state but does not prosecute drug vulnerability.
- Existing pan-cancer SL inference engines (ISLE, SELECT, etc.) systematically under-sample thyroid; their inferred pairs are not directly transferable.
- DM1-like dedifferentiation (FOXE1/NKX2-1/PAX8/HHEX collapse + STAT3/AP-1/DNMT activation) is an unusual transcriptional state — lineage collapse plus compensatory rewiring — and is precisely the regime where SL theory has the most leverage (lost redundancy → narrow survival corridor).

Risks: lack of thyroid-specific CRISPR depth in DepMap; PRISM/GDSC drug coverage is sparse for thyroid; SL inference is statistically delicate. These are *manageable* but require honest framing (Sec 12).

---

## 2. Why Paper 9 exists after Paper 1, not in parallel

| Paper | Asks | Establishes |
|---|---|---|
| Paper 1 | What is DM1? | A reproducible transcriptional state + 8-gene readout |
| Paper 2 | How does Hashimoto-overlap PTC project onto DM1 morphologically/biologically? | DM1 visibility on H&E + immune layer |
| Paper 3 | Why is the dark thyroid cancer ICI-vulnerable? | Immunological vulnerability axis |
| Paper 4 (backlog) | Korean GD HLA Pan-Asian | independent immunogenetic line |
| **Paper 9** | **Given DM1, what kills it?** | **Druggable / wet-lab-testable vulnerability map** |

Sequencing rationale:
- Paper 9 needs the DM1 state to be *citable*, not "in prep". Otherwise reviewers fight the state definition instead of the SL claim.
- Paper 1 in print = anchor citation. Paper 9 inherits authority.
- Running Paper 9 in parallel would dilute Paper 1's identity ("8-gene readout becomes the basis for a drug map") — better as the *next* paper's headline.

---

## 3. Ruppin-style conceptual template (used as scaffolding, not as text)

We borrow concepts and method-class, not language. The conceptual scaffolding is publicly described in the field and may be re-implemented as long as we do not lift sentences, claim independent priority over their algorithms, or imply equivalence to ISLE/SELECT outputs.

Concepts we adopt:
1. **Synthetic lethality (SL) as a population-level statistical signal.** Pairs (A, B) where co-loss is depleted in tumors that survive (because co-loss kills the cell) — operationalized via co-occurrence / co-expression / survival statistics.
2. **Synthetic rescue (SR).** Partner whose up-regulation rescues fitness when the target is lost; predicts resistance.
3. **Transcriptome as therapy-selection layer.** Move beyond actionable mutations: use gene-expression state to predict drug response, especially in mutation-poor tumors (which DM1 partially is).
4. **Cross-data triangulation.** SL/SR candidates must show a coherent signal across (i) cell-line dependency, (ii) drug response, (iii) patient outcome — not just one.

Concepts we explicitly *do not* adopt or claim:
- We do not claim ISLE-equivalent or SELECT-equivalent inference.
- We do not lift their precomputed SL pair lists as ground truth.
- We do not benchmark our pipeline against theirs as a horse-race; we treat them as inspiration, not adversary.

If Ruppin-tools are run *for sanity check* (recommended but optional), they are cited as third-party reference, not as our pipeline.

---

## 4. Thyroid-specific adaptation

Why thyroid breaks pan-cancer SL out of the box:
- Thyroid lines in DepMap are few (≈12 named, fewer with full CRISPR depth); pan-cancer SL pairs are dominated by signal from over-represented lineages (breast, colon, lung).
- Thyroid mutation density is low; co-mutation-based SL inference is anemic.
- Thyroid lineage TFs (PAX8, FOXE1, NKX2-1) are not standard SL nodes in pan-cancer maps.

Adaptations:
1. **State-conditioned SL.** Instead of asking "is gene X SL with Y pan-cancer", ask "is X SL with Y *given a DM1-like transcriptional context*". Define DM1-like state at the cell-line level using the 8-gene compact readout + TF-collapse score, then condition SL inference on that state.
2. **Lineage-imputed pseudobulk.** Augment thin thyroid cell-line panel with thyroid PDX/organoid expression where available, treated as held-out validation, not training.
3. **Cross-tissue borrowing for power.** Where the DM1-like transcriptional state appears in non-thyroid lineages (e.g., dedifferentiated neuroendocrine, lineage-plastic prostate), borrow strength for SL inference but always validate the hit *back* in thyroid context.
4. **TF-collapse as the SL anchor, not driver mutation.** This is the single biggest conceptual departure from canonical Ruppin work, which is mutation-anchored. We are state-anchored.

---

## 5. Main hypothesis (single sentence + decomposition)

**H0:** DM1-like thyroid tumors carry a non-redundant set of synthetic-lethal vulnerabilities arising from lineage collapse (FOXE1/NKX2-1/PAX8/HHEX loss) and compensatory rewiring (STAT3/AP-1/DNMT activation), and these vulnerabilities are predictable from the 8-gene compact readout without requiring driver-mutation profiling.

Decomposition:
- **H1 (state):** DM1-like state is reproducibly callable in CCLE.
- **H2 (dependency):** A subset of CRISPR dependencies in DepMap is significantly enriched in DM1-like lines.
- **H3 (compensation):** A subset of these dependencies has a synthetic-rescue partner whose up-regulation predicts resistance.
- **H4 (drug):** Drug response in PRISM/GDSC/CTRP recapitulates H2 for druggable targets.
- **H5 (outcome):** In TCGA-THCA and GSE76039, high target-pathway expression in DM1-classified tumors stratifies survival — motivating but not proving therapeutic benefit.

---

## 6. Candidate target classes

For each: rationale (DM1-link), prior thyroid evidence, drug class, what would falsify, SL-vs-non-SL classification.

### 6.1 NAMPT / NAD salvage
- **Rationale:** STAT3-high + high transcriptional output → high NAD turnover. NAMPT is the rate-limiter of NAD salvage. Cancers with metabolic plasticity often become NAD-addicted.
- **Drugs:** FK866 (APO866), KPT-9274, OT-82.
- **SL?** Yes, candidate SL with PARP-pathway / OXPHOS partners; SR with NAPRT (NAD de novo arm) — NAPRT-low + NAMPT-high is the classic vulnerability window.
- **Falsifier:** if NAPRT is uniformly high in DM1 lines, vulnerability vanishes.

### 6.2 STAT3 / JAK
- **Rationale:** Direct readout of the DM1 inflammatory-rewiring signature.
- **Drugs:** ruxolitinib (JAK1/2), TYK2-selective (deucravacitinib-class), STAT3-selective (TTI-101, napabucasin caveat).
- **SL?** Often *not* SL — STAT3 inhibition alone gives partial response. Real SL targets are likely **upstream cytokine receptors (IL6R, OSMR)** or **parallel TFs (NFκB)** in DM1 context.
- **Falsifier:** STAT3 single-agent inhibitor sensitivity does not separate DM1-like vs. non-DM1 lines.

### 6.3 DNMT / epigenetic repression
- **Rationale:** DNMT activation in DM1 → hypermethylation of lineage TFs (FOXE1 promoter is a known target). DNMT inhibitors may *reactivate* lineage program.
- **Drugs:** decitabine, guadecitabine; combination with HDACi (vorinostat) or LSD1i.
- **SL?** Conceptually different: this is **lineage-rewiring**, not pure SL. We classify it as a **redifferentiation vulnerability** and frame separately.
- **Falsifier:** DNMTi does not restore FOXE1/PAX8/NKX2-1/HHEX expression in DM1-like cell lines.

### 6.4 LYN / SRC-family kinases
- **Rationale:** Compensatory tyrosine-kinase signaling in mutation-poor dedifferentiated thyroid; SFK activation is reported in advanced/anaplastic thyroid cancer.
- **Drugs:** dasatinib (broad SFK), bosutinib, saracatinib.
- **SL?** Plausible SL partner with PI3K/AKT or with STAT3.
- **Falsifier:** SFK dependency in DepMap is not enriched in DM1-like lines.

### 6.5 KCNN4 (KCa3.1)
- **Rationale:** Calcium-activated K+ channel; up-regulated in some dedifferentiated solid tumors; reported in immune-tumor crosstalk.
- **Drugs:** senicapoc, TRAM-34 (tool compound).
- **SL?** More likely a *modifier* than a hard SL; included because of CRISPR-screen reproducibility in adjacent contexts.
- **Falsifier:** No DepMap dependency signal at all in DM1-like context.

### 6.6 DNA damage / PARP
- **Rationale:** Lineage plasticity often correlates with replication stress and BRCA-ness signatures.
- **Drugs:** olaparib, talazoparib, niraparib; ATR inhibitors (ceralasertib).
- **SL?** Yes if replication-stress signature is elevated in DM1; SR partner is BRCA1/2 / HR pathway intact tumors that resist.
- **Falsifier:** No HRD signature / no replication-stress markers (CHEK1, ATR, RPA pathway) elevated in DM1-like tumors.

### 6.7 Redox / OXPHOS-glycolysis compensation
- **Rationale:** Lineage collapse often forces metabolic reprogramming. DM1 may rely on either OXPHOS (mitochondrial) or compensatory glutamine/glycolysis.
- **Drugs:** IACS-010759 (OXPHOS), CB-839 (glutaminase), 2-DG (research tool).
- **SL?** Complementary metabolic SL pair: GLS1 (glutaminase) ↔ MYC, or OXPHOS ↔ glycolytic capacity (LDHA/HK2).
- **Falsifier:** Metabolic flux signature in DM1 is no different from generic high-proliferation thyroid lines.

### 6.8 TROP2 — explicit non-SL ADC vulnerability (carved out)
- **Rationale:** TROP2 is a tumor-cell-surface marker, not a synthetic-lethal partner.
- **Drugs:** sacituzumab govitecan (ADC, payload SN-38).
- **SL?** **NO.** This is *expression-level vulnerability* for ADC delivery. We include it in the paper as a **parallel translational opportunity**, explicitly distinguished from the SL-map. Mixing the two categories has been a recurring failure mode in similar papers and we will not.
- **What we say:** "TROP2-high DM1 tumors are candidates for ADC therapy, by an orthogonal mechanism to the synthetic-lethal targets above."

---

## 7. Data needed later (no download now)

Acquisition is **post-marathon**. List for planning only.

| Source | Use | Acquisition |
|---|---|---|
| DepMap CRISPR (DependencyMap, Broad) | Per-gene dependency in cell lines, conditioned on DM1-like state | Public; figshare / depmap.org |
| CCLE expression | Define DM1-like cell lines | Public |
| PRISM Repurposing | Drug-response across ~5,000 compounds, broad cell-line panel | Public |
| GDSC1/GDSC2 | Drug-response, ~700 compounds, deeper IC50 fits | Public |
| CTRP v2 | Redundant drug-response cross-check | Public |
| TCGA-THCA RNA + mutation + clinical | Patient-cohort validation; outcome stratification | dbGaP/UCSC Xena (already locally available) |
| GSE76039 | Advanced-disease replication; we already use this in Paper 1 | Already on disk |
| Thyroid PDX expression panels (where published) | Held-out validation, not training | Selective public sets |
| Organoid datasets (where available) | Wet-lab grounding for top hits | Selective public |

What we will NOT pull until Paper 1 is in print:
- DepMap, PRISM, GDSC, CTRP, organoid sets — all deferred. This document reserves the slot, it does not start the work.

---

## 8. Analysis plan (when work begins)

Sequenced, not started. Each step has an explicit GO/NO-GO criterion.

**Step A — Define DM1-like state at cell-line level.**
- Score CCLE lines on the 8-gene compact readout + TF-collapse composite (FOXE1/NKX2-1/PAX8/HHEX).
- Calibrate against thyroid lines we already trust as DM1-like (anaplastic, poorly differentiated).
- GO if ≥6 thyroid lines + ≥30 cross-tissue lines hit DM1-like threshold.
- NO-GO recovery: relax cross-tissue threshold; re-anchor on TF-collapse only.

**Step B — Compute target dependency.**
- For each gene, contrast CRISPR dependency score in DM1-like vs. non-DM1-like lines.
- Multiple-testing correction across the genome.
- Stratify by lineage to detect thyroid-specific vs. shared.
- GO if any of {NAMPT, NAPRT, STAT3, JAK1/2, DNMT1, DNMT3B, LYN, FYN, SRC, GLS1, ATR} shows DM1-conditioned enrichment.

**Step C — Infer SL / SR pairs.**
- For top dependencies, identify partners whose co-loss / co-up-regulation modulates dependency.
- Cross-validate against Ruppin pan-cancer SL pair lists *as a sanity probe only* — not as ground truth.
- GO if at least 2 partner-gated SL pairs survive.

**Step D — Rank druggable vulnerabilities.**
- Cross-reference Step B/C hits with PRISM/GDSC drug-target maps.
- Score by: (i) effect size in DM1-like lines, (ii) target druggability, (iii) clinical-stage drug availability, (iv) thyroid-cell-line confirmation.
- Output: ranked top-10 vulnerability list.

**Step E — Validate in patient cohorts.**
- In TCGA-THCA: stratify by DM1 classifier; test whether high target-pathway score associates with survival, recurrence, RAI-refractoriness flags.
- Replicate in GSE76039 (advanced-disease, smaller N).
- GO if at least 1 target shows concordant direction across both.

**Step F — Stage results into figure plan F1–F6.**

---

## 9. Figure plan F1–F6

- **F1 — Concept.** Schematic: DM1 state → lost redundancy → narrowed compensatory pathways → druggable nodes. Adapted from Ruppin-style SL/SR concept; visually distinct from any published ISLE/SELECT figure.
- **F2 — DM1-like state in CCLE.** Heatmap of 8-gene + TF-collapse score across lines; thyroid lines flagged; cross-tissue DM1-like lines marked.
- **F3 — Genome-wide CRISPR dependency, DM1 vs. non-DM1.** Volcano. Highlight top 20 hits; color-code by target class (NAD, JAK/STAT, epigenetic, SFK, metabolic, DNA damage).
- **F4 — Top SL/SR pairs.** Ten-row mini-grid: target gene × partner gene, with effect sizes for SL and SR axes.
- **F5 — Drug response (PRISM / GDSC).** IC50 distributions for top druggable targets, DM1-like vs. non-DM1-like. Add NAMPT/NAPRT, JAK, DNMT, SFK panels.
- **F6 — Patient validation (TCGA-THCA + GSE76039).** Survival stratification by DM1 + target-pathway score. Forest plot for hazard ratios.

Supplementary: per-target sensitivity panels; cell-line list with DM1 score; SL/SR statistical tables; sanity-check overlap with Ruppin pan-cancer SL lists (as reference, not benchmark).

---

## 10. Wet-lab validation options (not committed; describes the landscape)

Tier 1 — minimum to defend the paper:
- siRNA / CRISPR knock-down in 2–3 thyroid cell lines (one DM1-like, one non-DM1-like control) for the top 2 SL targets.
- Drug-response curves (IC50 + viability) for matching small-molecule inhibitors.

Tier 2 — strong but optional:
- Thyroid organoid panel (if accessible via collaborator) — drug response on top 3 hits.
- Combination experiments (target + SR partner inhibitor) to demonstrate the SR rescue logic.

Tier 3 — beyond first paper, for follow-up:
- PDX in vivo for top 1 hit.
- Patient-derived organoid bank cross-screen.

Decision: paper is *defendable* with Tier 1; Tier 2 is a stretch for first round; Tier 3 belongs to a follow-up.

Cost / collaborator question: thyroid wet-lab access. Open. If no wet-lab partner, paper still publishable as **computational + patient-cohort** vulnerability map, framed as hypothesis-generating, in a translational journal (JCI Insight, Mol Cancer Therap, NPJ Precision Oncology) rather than a flagship outlet.

---

## 11. Business / IP angle

- **Composition-of-matter:** not patentable for repurposed approved drugs (NAMPT inhibitors in trials, ruxolitinib, decitabine, dasatinib, olaparib, IACS-010759, CB-839 — all third-party).
- **Method-of-use / method-of-treatment:** potentially patentable: "method of selecting thyroid cancer patients for [target-class] therapy by measuring the 8-gene DM1 readout + TF-collapse score." This is the actionable IP slot.
- **Companion diagnostic:** the 8-gene readout itself is the most defensible asset; CDx claims tied to specific drug classes are filable separately.
- **Trial design:** umbrella / basket trial design where the DM1 score routes patients to one of 3–4 target arms (NAD, JAK/STAT, epigenetic, ADC). Pitchable to industry partners after Paper 9 publishes.
- **Defensive citation:** cite Ruppin-group methods as the conceptual prior art *generously* to avoid any infringement appearance and to position our work as thyroid-specific, state-anchored, clinically translatable.

---

## 12. What not to claim

Hard rules. Violating any of these wrecks reviewer trust.

1. Do **not** claim direct equivalence with ISLE/SELECT outputs. We are conceptually inspired, not replicating.
2. Do **not** claim clinical efficacy. Patient-cohort survival stratification is hypothesis-generating, not proof of therapeutic benefit.
3. Do **not** claim DM1-readout replaces driver-mutation profiling. Frame as **complementary**, especially in mutation-poor / RAI-refractory settings.
4. Do **not** mix TROP2-ADC vulnerability with synthetic lethality. They are different mechanistic categories. Carved out in Sec 6.8.
5. Do **not** over-extend to Hashimoto-overlap or Graves' subtypes. That is Paper 2 / Paper 4 territory; Paper 9 is specifically the **DM1 dark-matter axis**.
6. Do **not** use the H&E-DM1 angle. NO-GO per existing standing decision. Stay in molecular / drug-response / outcome space.
7. Do **not** claim wet-lab validation we have not done. If Tier 1 is missing, label the paper "computational + cohort-validated, wet-lab pending" explicitly in the abstract.
8. Do **not** lift sentences, signature lists, or precomputed pair tables from Ruppin-group papers. Concept yes, text/data no.
9. Do **not** start any of this before Paper 1 is in print and the marathon is closed (2026-06-13 floor).

---

## 13. Why this is Paper 9, not Paper 1 / 2 / 3

| Question | Answer |
|---|---|
| Could this be folded into Paper 1? | No. Paper 1 establishes biology; folding in SL inflates scope and invites reviewer mutiny on the state definition. |
| Could this be Paper 2? | No. Paper 2 is the Hashimoto-overlap PTC paper (HT-isolated per 2026-05-04 Yu decision). Different cohort, different question. |
| Could this be Paper 3? | No. Paper 3 is the ICI-vulnerability dark thyroid cancer paper (FROZEN bundle). Synthetic lethality is a different therapeutic axis; mixing them dilutes both. |
| Why number 9 specifically? | Reserves a slot in the long-range pipeline. Papers 5–8 are reserved for [paper-numbering doc owner-defined slots]; Paper 9 is the first translational follow-up to the DM1 axis. The number is conventional, not literal — reassign if numbering changes. |
| Why not just publish target list as a follow-up letter? | Letter-format does not support the figure-and-cohort weight of F1–F6 + supplements. Full article is the right vehicle. |

---

## 14. Twelve-week future execution plan (post-marathon)

Floor: cannot start before manuscript marathon closes (2026-06-13). Earliest practical start: 2026-06-16. Twelve-week sprint plan, each week is a checkpoint, not a deadline.

| Wk | Phase | Output |
|---|---|---|
| 1 | DepMap + CCLE pull, integrity checks | Local mirror, version-pinned |
| 2 | Define DM1-like state in CCLE; score thyroid lines | DM1 score table, cell-line ranking |
| 3 | Genome-wide DM1-conditioned dependency analysis | Volcano (F3 draft) |
| 4 | SL / SR pair inference | Pair table, sanity cross-check vs. Ruppin lists |
| 5 | PRISM / GDSC pull and integration | Drug-response matrix |
| 6 | Drug-response stratification by DM1; F5 draft | Per-target IC50 panels |
| 7 | TCGA-THCA + GSE76039 patient validation | Survival forest plots (F6 draft) |
| 8 | Cross-checks: lineage, batch, confounders | Sensitivity tables |
| 9 | Figure assembly F1–F6, supplements | Figure pack v1 |
| 10 | Draft writing (intro, results, discussion) | Manuscript v0 |
| 11 | Internal review with Yu / collaborator(s); revise | Manuscript v1 |
| 12 | Submission package: cover letter, COI, key resources | Submitted |

Wet-lab Tier 1 is a parallel track running weeks 4–12 if a collaborator is available; otherwise the paper is submitted as computational + cohort-validated and Tier 1 is the response-to-reviewer plan.

Stop conditions:
- If Step A (DM1-like state in CCLE) fails: pivot to bulk pseudobulk validation only and reframe scope. Do not fabricate state membership.
- If Step C (SL pairs) yields zero hits: paper still publishable as a *negative* result + ranked dependency list, but in a lower-tier venue. Decide explicitly at week 4 checkpoint.
- If wet-lab Tier 1 fails for top hits: replace top hits with rank-2/3 vulnerabilities, do not bury the failure.

---

## End of strategic plan

This document is design only. No data has been pulled, no analysis run, no GPU launched, no manuscript text written. Paper 1, 2, 3, 4 not touched. Voice-protected sections not authored.

Next action gates (in order):
1. Paper 1 in print.
2. Marathon closed (≥ 2026-06-13).
3. Yu sign-off on Paper 9 slot.
4. Then and only then: Step A.
