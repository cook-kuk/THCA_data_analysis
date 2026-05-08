# MHCflurry 2.0 — leakage-stratified ITSNdb scoring

## Version + install

- `mhcflurry==2.2.1` (pip; latest 2.x line — install spec `pip install mhcflurry`)
- Models: `models_class1_presentation` (pre-2.0 release, 2020-06-11), already cached
  at `~/.local/share/mhcflurry/4/2.2.0/`. No fresh download needed.
- Predictor: `Class1PresentationPredictor.load()` (combined affinity + processing →
  `presentation_score` ∈ [0, 1]).

## Bundle filter

Source: `project/results/p_neo_bayesian_2026_05_09/bundle.tsv`

Splits scored:

| split             | n   |
|-------------------|-----|
| ext_itsndb_main   | 199 |
| ext_itsndb_val    | 120 |
| ext_venus_test    |  78 |
| ext_venus_valid   |  78 |

ITSNdb total = **319** rows after split filter (note: spec said 311, actual bundle
is 319 — used the 319 actually present).

8–11-mer + standard AA filter:
- 156 rows dropped for `pep_len ∉ [8, 11]` — **all 156 were VenusVaccine** (lengths
  31–618 aa; these are full antigens, not MHC-I epitopes, and are out of scope for
  Class I presentation prediction).
- 0 rows dropped for non-standard amino acids.

→ Final scored set = **319 ITSNdb rows, 0 VenusVaccine rows.**

## Skipped HLAs

`predictor.supported_alleles` covers 14,883 alleles; **0 ITSNdb rows skipped** —
all HLA-A/B/C alleles in the leakage-stratified ITSNdb evaluation pool are
supported by the pan-allele model.

## Stratification

`in_master` flag = peptide overlaps the master training pool used by the wave-1/2
in-house bayesian model. Used here to split MHCflurry's external-set AUROC into:

| stratum             | n   | n_pos | label note                                  |
|---------------------|-----|-------|---------------------------------------------|
| ITSNdb_combined     | 319 | 136   | full external set                           |
| ITSNdb_no_overlap   | 106 |  33   | `in_master == False` (truly external)       |
| ITSNdb_in_master    | 213 | 103   | `in_master == True` (training-overlapping)  |

## Results (presentation_score → AUROC, 1000-bootstrap 95% CI)

| testset             | n   | n_pos | AUROC | 95% CI lo | 95% CI hi |
|---------------------|-----|-------|-------|-----------|-----------|
| ITSNdb_combined     | 319 | 136   | 0.637 | 0.573     | 0.700     |
| **ITSNdb_no_overlap** | 106 |  33   | **0.668** | 0.536     | 0.787     |
| ITSNdb_in_master    | 213 | 103   | 0.642 | 0.566     | 0.717     |

Δ (in_master − no_overlap) = **−0.026**; CIs heavily overlap. MHCflurry shows
*no* leakage signature on this stratification — its AUROC is essentially the
same on overlap vs no-overlap subsets (slightly higher on no-overlap, in fact).
This is consistent with the fact that "in_master" reflects overlap with the
*in-house* bayesian training set, not MHCflurry's own training distribution.

## VenusVaccine

Not scored: 156/156 peptides are full antigens (length 31–618), incompatible with
Class I 8–11-mer windowing. Would need explicit sliding-window epitope generation
+ a defined evaluation protocol to produce a meaningful AUROC, which is outside
the brief.

## Runtime

- Predictor load: ~9 s
- `predict()` on 319 rows (paired sample-name mode): **14.8 s**
- Total wall: **28.7 s** on CPU (no GPU used).

## Scripts + outputs

- `run_score.py` — entry point (paired sample-name mode for one-allele-per-row
  inference; sort by `peptide_num` to align scores back to bundle rows).
- `predictions.tsv` — 319 rows (peptide, hla, label, in_master, score_mhcflurry,
  source, split).
- `auroc_summary.tsv` — 3 ITSNdb rows (combined / no_overlap / in_master).
- `_run_meta.json` — runtime + skip stats.

## Caveats

- ITSNdb_no_overlap n_pos = 33 ⇒ wide CI (0.536–0.787). Δ direction is small and
  not significant under the bootstrap.
- `presentation_score` combines affinity × processing; for pure
  binding-affinity-only comparison use `affinity` (lower = better — sign flip).
  Brief asked for `presentation_score`, used as-is.
