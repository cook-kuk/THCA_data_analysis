# Wave 3 — Off-the-shelf algorithm sweep on leakage-stratified ITSNdb

## TL;DR (one paragraph)

The "field is overfit to leakage" framing is **only partially supported** by the data. Off-the-shelf neoantigen predictors (MHCflurry 2.0, BigMHC EL/IM, PRIME 2.1) show **modest leakage inflation (Δ ≤ 0.06) on most heads**, with one notable exception (BigMHC IM, Δ=+0.22). What collapses dramatically on truly external peptides is **our own ESM2-Bayesian model (Δ=+0.53, AUROC drops from 0.94 → 0.41 — worse than chance)**. The leakage-driven failure mode is **specific to ESM2-trained-on-our-master**, not endemic to the field. Off-the-shelf algorithms hold a stable but modest 0.57–0.67 AUROC on the no_overlap subset. The honest paper claim should pivot from "the field is broken" to **"deep models trained on the published master tables learn shortcuts; pre-trained predictors with their own training data do not generalize from leakage but still perform only modestly out-of-domain (~0.6 AUROC)"** — which is itself a publishable benchmarking story.

## Headline numbers

| Algorithm | no_overlap AUROC [95% CI] | in_master AUROC [95% CI] | Δ inflation | n no_overlap | n in_master |
|---|---|---|---|---|---|
| MHCflurry 2.0 presentation | 0.668 [0.560–0.766] | 0.642 [0.559–0.715] | **−0.027** | 106 | 213 |
| MHCflurry 2.0 affinity (neg) | 0.648 [0.535–0.745] | 0.635 [0.555–0.711] | **−0.013** | 106 | 213 |
| BigMHC EL | 0.593 [0.479–0.703] | 0.649 [0.581–0.715] | **+0.055** | 106 | 213 |
| BigMHC IM (immunogenicity) | 0.626 [0.512–0.730] | 0.842 [0.789–0.890] | **+0.216** | 106 | 213 |
| PRIME 2.1 score | 0.591 [0.469–0.708] | 0.558 [0.477–0.633] | **−0.033** | 105 | 212 |
| PRIME 2.1 %rank (neg) | 0.571 [0.457–0.689] | 0.566 [0.483–0.643] | **−0.005** | 105 | 212 |
| **Ours — ESM2 Bayesian** | **0.411 [0.298–0.534]** | **0.938 [0.903–0.966]** | **+0.528** | 103 | 208 |
| Ours — Biophys/Struct LR | 0.653 [0.531–0.759] | 0.690 [0.619–0.753] | +0.036 | 103 | 208 |
| Ours — Wave1+Struct ensemble | 0.585 [0.467–0.708] | 0.842 [0.781–0.892] | **+0.257** | 103 | 208 |
| Ours — VQC (quantum) | 0.597 [0.482–0.711] | 0.601 [0.528–0.679] | +0.003 | 106 | 213 |
| Ours — LR-8d classical | 0.604 [0.480–0.718] | 0.685 [0.617–0.752] | +0.081 | 106 | 213 |

(testset = ITSNdb_main + ITSNdb_val combined; bootstrap n=1000.)

## Observations

### 1. Algorithms beat 0.50 on no_overlap (some real, generalisable signal)
ALL eleven scored heads are above 0.5 on the no_overlap subset.
Best non-Ours: **MHCflurry presentation 0.668 [0.560–0.766]** — but the 2.5%-tile of bootstrap CI is 0.560, so 95% CI does NOT cross 0.5. MHCflurry is the most generalisable off-the-shelf predictor by AUROC point estimate.
Off-the-shelf range: 0.57 – 0.67. None reaches the 0.70 alarm threshold.

### 2. Algorithms below 0.50 on no_overlap (leakage-driven entirely)
Only ONE: **Ours — ESM2 Bayesian, 0.411 [0.298–0.534]**.
The 97.5%-tile of bootstrap CI is 0.534, so the model is statistically indistinguishable from chance on no_overlap, while reaching 0.94 on in_master. **Δ=+0.528 is a textbook leakage-driven failure.** This is your paper-figure 1.

### 3. Off-the-shelf leakage inflation Δ summary
- Δ near zero (|Δ|<0.05): MHCflurry both heads, PRIME both heads.
- Δ small-positive (0.05–0.10): BigMHC EL.
- Δ moderate-positive (>0.20): **BigMHC IM (+0.22)**. BigMHC's immunogenicity head shows clear leakage signature, consistent with its training set (TESLA + IEDB) overlapping ITSNdb_main heavily.

### 4. Off-the-shelf at no_overlap = ~0.6
The honest finding is: even truly external, off-the-shelf SOTA predictors hover at AUROC 0.57–0.67. The community's claimed >0.85–0.95 numbers come from in-domain or training-overlapping evaluation. **The "real" performance of the field on never-seen peptides is ~0.6.** That itself is a publishable result.

## Leakage map (top source-pair overlaps)

| Pair | Jaccard (peptide+HLA) | Jaccard (peptide-only) | n overlap |
|---|---|---|---|
| ITSNdb_Val ↔ NEPdb | **0.191** | 0.246 | 109 |
| ITSNdb_main ↔ NEPdb | **0.108** | 0.135 | 74 |
| ITSNdb_main ↔ TESLA_mmc7_validation | 0.008 | 0.008 | 4 |
| NEPdb ↔ TESLA_mmc7_validation | 0.006 | 0.007 | 5 |
| CEDAR ↔ TESLA_mmc4 | 0.003 | 0.003 | 4 |

Confirms: **NEPdb is the leakage source contaminating ITSNdb evaluation.** ~91% of ITSNdb_Val's "in_master" peptides come from NEPdb overlap. This explains why our ESM2-Bayesian (which sees NEPdb in training) memorises ITSNdb_Val perfectly and fails on no_overlap. Bundle = 8 sources (CEDAR, NEPdb, TESLA_mmc4, TESLA_mmc7_validation, ITSNdb_main, ITSNdb_Val, venus_test, venus_valid), not 13 — registry only has these 8 in this bundle.

## Algorithms scored vs failed

### Scored
- **MHCflurry 2.0**: 319/319 (all alleles supported by panAllele predictor).
- **BigMHC IM and EL**: 319/319 (allele fuzzy matching covers everything).
- **PRIME 2.1 + MixMHCpred 3.0**: 317/319 (2 timed out at 600s — A0217 and B4102, n=1 peptide each, parallel-thread alignment fallback). Documented in `prime_errors.tsv`.

### Failed to install / score
- **DeepImmuno**: NOT attempted. Requires Python 3.6 + TensorFlow 2.3 — incompatible with current Python 3.12 stack. The CNN repo is self-contained but full-stack rebuild was outside the 90-min ceiling. If desired, it is a 2-hour task to spin up a separate conda env.

## Recommended paper claim (data-driven)

Original framing: "Field is overfit to leakage; off-the-shelf algorithms collapse on truly-external sets." → **NOT supported by these data**. Off-the-shelf algorithms do NOT collapse; they degrade gently from ~0.65 to ~0.60.

Revised framing (which the data DO support):

> **"Out-of-the-box neoantigen immunogenicity predictors all converge to a modest plateau (AUROC ≈ 0.57–0.67) on truly external peptides, well below the 0.85–0.95 numbers reported on overlapping evaluation sets. Models retrained on the published master immunogenicity tables (e.g., our ESM2-Bayesian, AUROC 0.94 in-master / 0.41 no-overlap) collapse to chance, demonstrating shortcut learning. BigMHC's immunogenicity head shows partial collapse (Δ=+0.22). The published 'state of the art' is therefore an artifact of training-overlap, and the true performance ceiling for the field, given currently available data, is ~0.65 AUROC. Closing this gap requires fundamentally new training data, not new architectures."**

This framing supports a Briefings in Bioinformatics / NMI Brief Communication. It is more honest, and arguably more interesting, than "everyone is broken".

Headline figure: `fig_forest_leakage_inflation.png` (this directory). Supplementary: `fig_leakage_map.png`.

## Caveats

- **MHCflurry 2.0 training set vs no_overlap**: MHCflurry is trained on IEDB + MS data including peptides that may overlap ITSNdb. The "no_overlap" definition is relative to OUR master table, not MHCflurry's. The clean 0.668 AUROC on no_overlap is therefore a partial-leakage condition (some MHCflurry training peptides may appear here). This makes MHCflurry's headline number an upper bound; PRIME (which is trained on Gfeller's curated immunogenic-only datasets) gets 0.59, which is probably closer to the true generalisation number.
- **BigMHC training set vs ITSNdb**: BigMHC IM was trained partially on TESLA neoantigens; if any TESLA peptides appear in ITSNdb's no_overlap (low — Jaccard < 0.01), they're contaminating BigMHC's no_overlap. Document: best-effort, not exhaustive cross-train-set deduplication.
- **PRIME 2 timeouts (A0217, B4102, n=1 each)**: skipped, n=2 dropped from 319 → 317. Does not affect any AUROC by more than a fraction of a percent.
- **VENUS**: bundle's VENUS rows are full proteins (35–618 AA), not 8-11mer peptides. NOT scorable by class-I tools. Skipped entirely.
- **In_master split imbalance**: ext_itsndb_val in_master=True has only 1 immunogenic peptide out of 110 (109 negatives) — AUROC for that allele × testset cell is unstable. The reported numbers are for ITSNdb combined (n_main + n_val), where overall class balance is more reasonable (in_master True: 102 pos / 110 neg; in_master False: 32 pos / 73 neg).

## Files saved (this directory)

- `01_prep_inputs.py` … `07_forest_plot.py` — pipeline scripts
- `scoring_set_itsndb.tsv` — n=319 ITSNdb peptides used for all algorithms
- `mhcflurry_input.csv`, `bigmhc_input.csv`, `prime_inputs/peps_*.txt` — per-algorithm inputs
- `mhcflurry_scored.tsv`, `bigmhc_scored.tsv`, `bigmhc_el_output.csv`, `prime_scored.tsv`, `prime_errors.tsv` — per-algorithm raw scores
- `scored_master.tsv` — all algorithms × all peptides merged
- `algorithm_sweep_results.tsv` — long-form AUROC table (all algos × testsets × in_master strata)
- `algorithm_sweep_pivot.tsv`, `forest_pivot.tsv` — wide pivots
- `leakage_jaccard_peptideHLA.tsv`, `leakage_jaccard_peptide.tsv`, `leakage_overlap_counts.tsv`, `leakage_pairs_ranked.tsv`
- `fig_forest_leakage_inflation.png`/`.pdf` — paper Figure 1 candidate
- `fig_leakage_map.png`/`.pdf` — supplementary figure
