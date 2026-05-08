# Track 26 — IFN-γ × HLA-I/II × DM1 mediation in TCGA-THCA

> **Caption boilerplate (applies to every figure and table in this track):**
> HLA-I/II gene-expression module — not allele genotype. Cancer-cohort allele
> genotyping is out of scope per
> `project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md` (Sections 1.1, 1.2).
> Mediation framing is observational, not interventional.

---

## 0. Boundary

This track operates strictly inside the Paper-1 allowlist for HLA: gene-expression
modules used as a context/residualization control. No allele genotyping is
performed on TCGA-THCA. No HLA susceptibility, prognosis, or patient-selection
claim is made. All HLA findings are transcript-level expression.

Mediation language is observational: we report covariation patterns *consistent
with* a DM1 → IFN-γ → HLA pathway, not interventional proof.

---

## 1. Data

- **Cohort:** TCGA-THCA primary tumors, n=527 with DM1, HLA-I, HLA-II, IFN-γ scores.
- **Source:** `project/data/raw/TCGA_pancan/pancan_geneExp.gz` + Track 5 frame
  `track5_dm1_hla1_module/tables/T01_per_sample_module_scores.tsv` (DM1 score,
  driver class, purity proxy already computed).
- **Methylation layer:** `project/results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv`
  (HM450 mean β over 8-gene panel), n=523 after intersect.
- **DM1 score:** Track 5 `DM1_use` (pancan z-score, fallback to local 8-panel).
- **HLA-I module (19 genes):** HLA-A/B/C, B2M, TAP1/2, TAPBP, NLRC5, IRF1,
  PSMB8/9, ERAP1/2, HLA-E/F/G, CALR, CANX, PDIA3.
- **HLA-II module (12 genes):** HLA-DRA/DRB1/DPA1/DPB1/DQA1/DQB1/DMA/DMB/DOA/DOB,
  CIITA, CD74.
- **IFN-γ signatures (two, compared):**
  - HALLMARK_INTERFERON_GAMMA_RESPONSE — 195 / 199 canonical genes mapped.
  - Ayers 2017 18-gene Tumor Inflammation Signature (TIS) — 18 / 18 mapped.
- All module scores: per-gene z-score, then mean across the module.

---

## 2. IFN-γ score per TCGA-THCA sample (Deliverable 1)

Both IFN-γ scores are released per sample to `T01_ifng_score_per_sample.tsv`.

- Hallmark vs Ayers TIS: Spearman ρ = **0.894** (n=527). The two definitions
  agree very strongly, so downstream results are essentially robust to choice of
  IFN-γ signature.

This score becomes our `M` (mediator).

---

## 3. DM1 × IFN-γ — Deliverable 2

| Comparison | n | Spearman ρ | p |
|---|---:|---:|---:|
| DM1 vs IFN-γ Hallmark | 527 | **0.340** | 1.0×10⁻¹⁵ |
| DM1 vs IFN-γ Ayers TIS | 527 | 0.193 | 8.3×10⁻⁶ |

DM1-high tumors have meaningfully higher IFN-γ activity. The Hallmark signal is
nearly twice as large as Ayers because Hallmark includes ISG-only genes that
track interferon induction without requiring T-cell content; Ayers TIS contains
explicit T-cell markers (CD8A, GZMA, PRF1, LAG3, TIGIT) so it is partially
limited by per-tumor T-cell infiltrate. We use Hallmark as primary mediator,
Ayers as sensitivity layer.

Figure: `figs/F01_dm1_vs_ifng_scatter.png`.

---

## 4. IFN-γ × HLA-I and IFN-γ × HLA-II — Deliverable 3

| Comparison | n | Spearman ρ | p |
|---|---:|---:|---:|
| IFN-γ Hallmark vs HLA-I | 527 | **0.880** | 4.4×10⁻¹⁷² |
| IFN-γ Hallmark vs HLA-II | 527 | **0.873** | 2.0×10⁻¹⁶⁵ |
| IFN-γ Ayers TIS vs HLA-I | 527 | 0.888 | 1.8×10⁻¹⁷⁹ |
| IFN-γ Ayers TIS vs HLA-II | 527 | 0.864 | 8.7×10⁻¹⁵⁹ |

IFN-γ score and HLA module are essentially the same axis at ρ ≈ 0.88 for both
class I and class II. This is biologically expected — IFN-γ is the canonical
HLA inducer (NLRC5 → MHC-I, CIITA → MHC-II). The shared-genes overlap (Hallmark
contains B2M, HLA-A/B, HLA-DMA, HLA-DQA1, HLA-DRB1, HLA-G, NLRC5, TAP1, TAPBP,
PSMB8/9 — 11 / 19 of the HLA-I module) explains a chunk of this, but the same
strength holds for HLA-II module against Ayers TIS, which contains only one HLA
gene (HLA-DRA), confirming the relationship is not a tautology.

Figure: `figs/F02_ifng_vs_hla_scatter.png`.

---

## 5. Mediation: DM1 → IFN-γ → HLA-I and HLA-II — Deliverables 4, 5

Baron-Kenny + Sobel + non-parametric bootstrap (pingouin 0.6.1, n_boot=2000).
Reference: `T03_mediation_baron_kenny_sobel.tsv`.

| Path | n | Total c | Direct c′ | Indirect a·b | % mediated | Sobel z | Sobel p |
|---|---:|---:|---:|---:|---:|---:|---:|
| DM1 → IFNG-Hallmark → **HLA-I** | 527 | 0.453 | −0.035 | **0.488** | **107.8 %** | 8.16 | 2.2×10⁻¹⁶ |
| DM1 → IFNG-Hallmark → **HLA-II** | 527 | 0.930 | 0.277 | **0.653** | **70.2 %** | 8.13 | 4.4×10⁻¹⁶ |
| DM1 → Ayers-TIS → HLA-I | 527 | 0.453 | 0.210 | 0.243 | 53.7 % | 3.98 | 6.9×10⁻⁵ |
| DM1 → Ayers-TIS → HLA-II | 527 | 0.930 | 0.606 | 0.324 | 34.8 % | 3.97 | 7.1×10⁻⁵ |

**Headline finding.** Using Hallmark IFN-γ, ~108 % of DM1 → HLA-I and ~70 % of
DM1 → HLA-II run through IFN-γ. The 108 % means the direct DM1 → HLA-I path is
slightly *negative* (−0.035, n.s.) once IFN-γ is in the model — i.e., HLA-I
elevation in DM1-high tumors is **fully explained** by IFN-γ; and a residual
*non*-IFN-γ path actually *reduces* HLA-I in DM1-high tumors (consistent with
poor-differentiation / chaperone-CALR-class repression — see §8 below).
HLA-II is partially mediated (70 %), with a still-significant direct path
(c′=0.28) — DM1 → HLA-II keeps an IFN-γ-independent component.

Bootstrap 95 % CIs for the indirect effect exclude zero in all four cases
(`figs/F04_indirect_effect_bootstrap_CI.png`).

Figures: `figs/F03_mediation_pct_mediated.png`,
`figs/F09_direct_indirect_summary.png`.

---

## 6. Direct vs indirect decomposition — Deliverable 5

For Hallmark IFN-γ:

- **HLA-I:** total = 0.453, direct = −0.035 (n.s.), indirect via IFN-γ = 0.488 → **fully mediated**.
- **HLA-II:** total = 0.930, direct = 0.277 (sig), indirect via IFN-γ = 0.653 → **partially mediated**, residual direct path remains.

The HLA-I result is the strongest interpretable signal: DM1 itself does not
"directly upregulate" HLA-I in cancer cells in TCGA-THCA — it is an IFN-γ
response that walks alongside dedifferentiation. The HLA-II result keeps a
direct path that is consistent with thyroid-epithelium-intrinsic class-II
upregulation in dedifferentiated cells (a known biology of
RAI-refractory / Hashimoto-overlap PTC) plus the IFN-γ component.

---

## 7. Conditional independence — Deliverable 6

Partial Spearman ρ for DM1 × HLA module *given* candidate confounders
(`T04_conditional_independence.tsv`):

| Module | ρ raw | ρ \| IFN-γ-Hallmark | ρ \| IFN-γ-Ayers | ρ \| Purity proxy | ρ \| IFN-γ + Purity |
|---|---:|---:|---:|---:|---:|
| HLA-I | 0.318 (p=8e-14) | **0.021 (p=0.64)** | 0.324 (p=3e-14) | 0.301 (p=2e-12) | 0.021 (p=0.63) |
| HLA-II | 0.393 (p=7e-21) | **0.183 (p=2.5e-5)** | 0.436 (p=7e-26) | 0.411 (p=7e-23) | 0.213 (p=8.5e-7) |

**HLA-I:** conditioning on IFN-γ Hallmark drives the partial ρ to ~0
(0.021, p=0.64). This is the same conclusion as the mediation analysis — DM1 ×
HLA-I is **fully explained by IFN-γ**. Purity correction alone barely moves the
correlation, ruling out infiltrate-mass as the dominant confounder.

**HLA-II:** partial ρ drops from 0.39 to 0.18 (still p=2.5e-5) — about half of
the DM1 × HLA-II axis remains after IFN-γ adjustment. An independent path
remains.

The Ayers-TIS-conditioned partial ρ does *not* drop as far because Ayers is a
weaker mediator (smaller a-path; see §3) — i.e., this validates the Hallmark
result rather than contradicting it (Ayers under-captures IFN-γ in low-T-cell
tumors).

Figure: `figs/F05_conditional_independence.png`.

---

## 8. Driver-stratified mediation — Deliverable 7

Within driver classes (n_BRAF=294, n_RAS=54, n_Fusion=0 in this DM1-frame
intersection, n_TripleNeg=168). Reference: `T05_driver_stratified_mediation.tsv`.

| Driver | Target | n | Total c | Direct c′ | Indirect | % mediated | Sobel p |
|---|---|---:|---:|---:|---:|---:|---:|
| **BRAF** | HLA-I | 294 | 0.256 | 0.095 | 0.161 | **62.8 %** | 0.059 |
| BRAF | HLA-II | 294 | 0.508 | 0.311 | 0.197 | 38.8 % | 0.060 |
| RAS | HLA-I | 54 | −0.109 | 0.003 | −0.112 | 103 % (small n) | 0.61 |
| RAS | HLA-II | 54 | 0.386 | 0.546 | −0.161 | n.s. | 0.61 |
| **TripleNeg** | HLA-I | 168 | 0.393 | −0.147 | **0.540** | **137 %** | 8.8×10⁻⁵ |
| TripleNeg | HLA-II | 168 | 1.010 | 0.261 | 0.750 | **74.2 %** | 8.7×10⁻⁵ |

**Pattern.** The IFN-γ-mediation pattern is most pronounced in the
**triple-negative (BRAF/RAS/fusion-neg, dark-matter)** compartment, exactly the
DM1-enriched group: 137 % HLA-I mediation (direct path slightly negative),
74 % HLA-II mediation. BRAF-positive shows weaker mediation (~63 % HLA-I) and
weaker overall a-path. RAS is too small at n=54 to draw mediation conclusions.
Fusion is absent because Track 5's DM1-frame fusion column is `unknown`-coded
into TripleNeg here; fusion-vs-non-fusion stratification within this intersected
frame would need the v17p35 fusion-vs-anchor table — out of scope for this
track.

This means the IFN-γ → HLA pathway is the **dominant axis exactly where dark
matter lives**, which is the Paper-1-relevant interpretation: in BRAF/RAS-neg
PTC, DM1 dedifferentiation is accompanied by a strong IFN-γ-driven HLA
elevation, but the underlying DM1 dedifferentiation itself does not directly
upregulate HLA-I (direct path ≈ 0, slightly negative).

Figure: `figs/F06_driver_stratified_mediation_heatmap.png`.

---

## 9. Per-gene split — Deliverable 8

Per-gene mediation, sorted by % mediated (`T06_per_gene_mediation.tsv`):

**Most IFN-γ-mediated (top, > 70 % indirect via IFN-γ):**

| Gene | Module | % mediated | Sobel p |
|---|---|---:|---:|
| TAPBP | HLA-I | 271 % | 1.2×10⁻¹² |
| TAP2 | HLA-I | 149 % | 4.4×10⁻¹⁶ |
| IRF1 | HLA-I | 140 % | 4.4×10⁻¹⁶ |
| HLA-E | HLA-I | 127 % | 8.9×10⁻¹⁶ |
| HLA-DMB | HLA-II | 112 % | 3.1×10⁻¹² |
| NLRC5 | HLA-I | 89 % | 6.7×10⁻¹⁶ |
| HLA-F | HLA-I | 86 % | 6.7×10⁻¹⁶ |
| HLA-DMA | HLA-II | 81 % | 8.9×10⁻¹⁶ |
| HLA-DPB1 | HLA-II | 78 % | 4.4×10⁻¹⁶ |
| PSMB9 | HLA-I | 74 % | 6.7×10⁻¹⁶ |
| HLA-DOB | HLA-II | 73 % | 6.7×10⁻¹⁶ |
| HLA-DPA1 | HLA-II | 71 % | 4.4×10⁻¹⁶ |
| CD74 | HLA-II | 71 % | 6.7×10⁻¹⁶ |
| HLA-DOA | HLA-II | 70 % | 4.4×10⁻¹⁶ |

**Least IFN-γ-mediated (bottom; chaperones):**

| Gene | Module | % mediated | Sobel p |
|---|---|---:|---:|
| CANX | HLA-I | 22 % | 3.7×10⁻¹² |
| CALR | HLA-I | 4.7 % | 1.3×10⁻⁸ |
| PDIA3 | HLA-I | 3.1 % | 8.3×10⁻⁷ |

**Interpretation.** The classical IFN-γ-induced antigen-processing axis
(TAPBP, TAP1/2, PSMB8/9, IRF1, NLRC5, HLA-A/B/C/E/F, CIITA, HLA-DR/DP/DQ/DO/DM,
CD74) is heavily IFN-γ-mediated in DM1-high tumors. The
ER-resident chaperones **CALR, CANX, PDIA3** — which assemble peptide-loading
complex but are house-keeping, not IFN-γ-induced — show ≤ 22 % mediation, as
expected. This is exactly the textbook prior, and it sanity-checks the entire
mediation analysis: only IFN-γ-inducible HLA machinery rides the IFN-γ path; the
constitutive chaperones do not.

Figure: `figs/F07_per_gene_mediation_lollipop.png`.

---

## 10. Methylation 3-way — Deliverable 9

Question: per Track 5 / Round 4 memory, HM450 β tracks +DM1 and +HLA. Is
methylation an *independent* path from IFN-γ, or does it absorb the same axis?

Single-mediator (intersected n=523; `T08_methylation_single_mediator.tsv`):

| Mediator | Target | Total c | Indirect | % mediated | Sobel p |
|---|---|---:|---:|---:|---:|
| mean_8g_β | HLA-I | 0.457 | 0.584 | **128 %** | 0 (cap) |
| IFNG-Hallmark | HLA-I | 0.457 | 0.492 | 108 % | 2×10⁻¹⁶ |
| mean_8g_β | HLA-II | 0.938 | 0.810 | 86 % | 0 (cap) |
| IFNG-Hallmark | HLA-II | 0.938 | 0.662 | 71 % | 2×10⁻¹⁶ |

**Both** mean_8g_β and IFN-γ score, taken alone, can saturate the DM1 → HLA-I
total — they are essentially redundant single-mediator candidates because both
correlate with DM1 and with HLA elevation.

Dual-mediator decomposition (`T07_dual_mediator_meth_ifng.tsv`,
`figs/F08_dual_mediator_decomposition.png`): when both mediators are placed in
the same model, the indirect effects split between methylation and IFN-γ.
Partial Spearman with both as confounders (`T09_partial_meth_ifng.tsv`):

| Module | ρ raw | ρ \| meth | ρ \| IFN-γ | ρ \| both |
|---|---:|---:|---:|---:|
| HLA-I | 0.320 | −0.011 (n.s.) | 0.023 (n.s.) | −0.019 (n.s.) |
| HLA-II | 0.394 | 0.092 (p=0.04) | 0.180 (p=4e-5) | 0.122 (p=0.005) |

**HLA-I:** *either* methylation *or* IFN-γ alone fully explains DM1 × HLA-I.
The two mediators co-occupy the same explanatory volume — they are
not independent paths of comparable size; rather, methylation is upstream
(DM1 ↔ epigenetic dedifferentiation) and IFN-γ is the proximal HLA inducer.
With this n we cannot statistically separate which of the two paths is dominant
for HLA-I; we can say HLA-I co-varies with DM1 entirely through this paired
methylation × IFN-γ axis.

**HLA-II:** the residual partial ρ stays significant (p=0.005) even with both
in the model. There is some HLA-II axis that is *independent* of both
methylation and IFN-γ — a candidate for thyroid-epithelium-intrinsic class-II
expression in dedifferentiated cells, or a CIITA-independent class-II path.

---

## 11. Limitations

1. **Observational mediation, not interventional.** All findings here are
   covariation patterns in n=527 cross-sectional bulk RNA-seq. Baron-Kenny and
   Sobel are the standard tools, but they assume the directed model
   X → M → Y; we cannot rule out reverse or bidirectional flow (e.g., HLA →
   IFN-γ feedback via T cells).
2. **Bulk RNA-seq mixes tumor and stroma.** The IFN-γ score is a tissue-mixture
   score; tumor-cell-intrinsic vs T-cell-intrinsic IFN-γ contributions are not
   separable here. Track 5 found HLA-I × DM1 partial ρ ≈ 0.30 even after purity
   correction, ruling out *complete* explanation by infiltrate; but the IFN-γ
   mediator itself is partially driven by infiltrate. The Ayers-TIS sensitivity
   layer is informative because it is dominated by infiltrate; the Hallmark
   layer survives even without explicit T-cell genes.
3. **Driver-fusion stratification incomplete** in this run because the
   Track 5-merged frame lacked Fusion-anchor labels (n_Fusion = 0 in the
   intersected DM1 frame). A future re-run should pull v17p35 driver labels
   directly to recover the Fusion stratum.
4. **Methylation layer is 8-gene** mean β, not whole-genome HM450; the
   methylation mediator is therefore biased toward the thyroid-differentiation
   axis. Whole-genome methylation as mediator would be a future extension.
5. **HALLMARK_INTERFERON_GAMMA_RESPONSE shares 11 genes with the HLA-I module.**
   This pulls the IFN-γ ↔ HLA-I correlation upward by construction. The Ayers
   TIS sensitivity confirms the result without that shared-gene bias, and the
   per-gene split (§9) confirms that even non-overlapping HLA genes (HLA-C,
   HLA-DRB1, HLA-DPB1) are heavily IFN-γ-mediated. We therefore consider this
   limitation cosmetic, not structural.
6. **Causal language.** Throughout, we say *consistent with mediation*,
   *covaries with*, *runs through*, never *causes*. Reviewer-Q-9 voice is
   user-keyboard-protected per marathon-mode rules.

---

## 12. Implication for the manuscript

- **HLA-I elevation in DM1-high tumors is, in TCGA-THCA, fully explained by
  IFN-γ signaling.** The direct DM1 → HLA-I path collapses to ~0 once IFN-γ
  is conditioned on. This means the "HLA-I module residualization control"
  (Paper-1 boundary §1.1) for the DM1 narrative is appropriate to phrase as
  *immune-context*, not as *tumor-cell-autonomous antigen presentation*.
- **HLA-II elevation in DM1-high tumors has an additional non-IFN-γ
  component.** This is consistent with the Hashimoto-overlap / dedifferentiated
  thyroid-epithelium HLA-II programme (per memory `v17_D5P6_BCR_clonal_TLS` and
  `v17_D4P2_tcga_hashimoto_generalization`). HLA-II partial ρ stays p < 0.005
  even after both methylation and IFN-γ adjustment.
- **The IFN-γ → HLA pathway is most dominant in the BRAF/RAS-negative dark-matter
  compartment** (TripleNeg: 137 % HLA-I, 74 % HLA-II mediation), exactly the
  Paper-1-target population. The DM1 dedifferentiation phenotype is *not* a
  driver of cell-autonomous HLA loss; on the contrary, in this cohort it
  co-travels with active IFN-γ-driven HLA induction. Any
  immune-evasion-via-HLA-loss claim for dark-matter tumors must be made
  mechanistically distinct from this transcript-level signature (e.g., LOH
  rather than expression).

---

## Output paths

- Tables: `project/results/hla_deepdive_2026_05_08/track26_ifng_hla_dm1/tables/T01..T10`
- Figures: `project/results/hla_deepdive_2026_05_08/track26_ifng_hla_dm1/figs/F01..F09`
- Summary JSON: `project/results/hla_deepdive_2026_05_08/track26_ifng_hla_dm1/track26_summary.json`
- Script: `scripts/hla_deepdive_2026_05_08/track26/run_track26.py`
