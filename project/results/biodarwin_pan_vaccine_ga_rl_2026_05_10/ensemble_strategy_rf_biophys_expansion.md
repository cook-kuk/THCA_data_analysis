# BioDarwin ensemble strategy expansion: RF_biophys as the AUPRC anchor

Generated: 2026-05-11

## Why this branch exists

`RF_biophys` is not the final model. It is the strongest **ranking anchor** we currently have in the in-house stack when the goal is to pull positives toward the top of the list.

In the current wave11 comparisons it shows:

| dataset | AUPRC | AUROC | training overlap |
|---|---:|---:|---:|
| cedar_partial | 0.9967 | 0.9555 | 1.0000 |
| nepdb | 0.9300 | 0.9731 | 1.0000 |
| itsndb | 0.6501 | 0.7531 | 0.5768 |
| tesla_mmc4 | 0.3641 | 0.8839 | 1.0000 |

Those numbers tell us two different things at once:

1. It has real ranking power, especially in AUPRC.
2. It is not a clean OOD truth source, because the strongest rows are overlap-heavy.

So the right move is not to discard it and not to let it dominate. The right move is to **extract its ranking behavior and route it through a leakage-aware ensemble**.

## Core idea

Use `RF_biophys` as the **high-recall ranker** in a routed ensemble:

1. `RF_biophys` proposes candidate ordering.
2. Public-safe comparators and overlap-aware penalties regularize that ordering.
3. GA/RL learns which mode to trust for each data regime.

This keeps the RF strength where it is useful:

- early enrichment
- AUPRC gain on imbalanced sets
- robust ordering within a narrow source family

and prevents it from becoming the only decision-maker:

- overlap-heavy source memorization
- false confidence on pseudo-held-out rows
- brittle transfer to genuinely new industrial data

## What RF_biophys contributes

The current implementation is simple by design:

- 24 peptide biophysical features
- HLA one-hot encoding
- Random Forest classifier

That simplicity is useful because it gives us a strong, interpretable baseline that is cheap to score and easy to gate.

The useful signal from RF is not just its absolute score. It is:

- `rf_rank`
- `rf_margin`
- agreement with `BigMHC_IM`
- disagreement against `PRIME` / `MHCflurry`
- source-family stability

That means RF should be promoted from “model” to “feature family”.

## Proposed architecture

### 1. Discovery mode

Use for academic-like or internal candidate generation.

Recommended components:

- `RF_biophys`
- `BigMHC_IM`
- `BioDarwin_PanVax`
- `CROSS_claimsafe`

Suggested logic:

```text
discovery_score =
  a1 * rf_biophys_norm
  + a2 * bigmhc_im_norm
  + a3 * cross_claimsafe_norm
  + a4 * biodarwin_full_mode_norm
  + a5 * disagreement_bonus
  - a6 * redundancy_penalty
```

What this mode is for:

- pull top candidates fast
- maximize AUPRC-like ranking
- build candidate pools for wetlab

### 2. Public-anchor mode

Use for industrial/public locked fallback.

Recommended components:

- `BigMHC_IM`
- `RF_biophys`
- `BioDarwin_public_anchor_v3`

Suggested logic:

```text
public_anchor_score =
  0.50 * bigmhc_im_norm
  + 0.25 * rf_biophys_norm
  + 0.25 * public_anchor_v3_norm
  - overlap_risk_penalty
```

Why RF stays here:

- it gives rank compression that often helps early retrieval
- it can recover signal when BigMHC is too conservative
- it provides a cheap orthogonal baseline

Why it is not higher than that:

- overlap-sensitive datasets can inflate it
- it should help, not overrule, the public anchor

### 3. Disagreement mode

Use for 96-well selection.

The wetlab set should not be only top-score items. It should mix:

- high consensus
- high RF-only uplift
- high BigMHC-only uplift
- high model disagreement
- low-overlap novelty

Suggested selection score:

```text
plate_score =
  top_score
  + disagreement_weight * |rf_biophys_norm - bigmhc_im_norm|
  + novelty_weight * novelty_proxy
  + diversity_weight * source_family_diversity
  - overlap_weight * training_overlap_fraction
```

This is the part that can turn RF from a “good leaderboard model” into a useful experimental design component.

## GA/RL objective

The GA/RL controller should not optimize one raw metric.

It should optimize a vector objective:

- `AUPRC` on pre-freeze academic validation
- `Top10` / `Top24` on industrial locked mini-sets
- calibration stability
- overlap penalty
- disagreement coverage
- class balance robustness

In practical terms:

```text
fitness =
  w1 * validation_AUPRC
  + w2 * industrial_top10
  + w3 * consensus_stability
  + w4 * disagreement_gain
  - w5 * overlap_sensitivity
  - w6 * redundancy
```

RF_biophys should receive a **positive prior** in the search, but not a hard-coding.

## What to evolve

Instead of evolving a single score formula, evolve these knobs:

1. RF weight by regime
2. BigMHC weight by regime
3. public-anchor weight by regime
4. disagreement threshold
5. overlap penalty strength
6. source-family diversity bonus
7. top-k plate composition

This is the part where GA/RL actually helps.

It lets the system learn:

- when RF is a booster
- when RF is just noise
- when RF should only define candidate ordering

## Ablation plan

The clean comparison table should include:

| mode | RF included | overlap penalty | public anchor | expected behavior |
|---|---|---|---|---|
| base | no | no | no | weakest baseline |
| RF-anchor | yes | no | no | best raw AUPRC, but unsafe |
| RF + penalty | yes | yes | no | safer internal use |
| public-anchor | partial | yes | yes | industrial fallback |
| full mode-bank | yes | yes | yes | best overall candidate system |

The most important ablation is:

- `RF-only`
- `RF + overlap penalty`
- `RF + BigMHC`
- `RF + BigMHC + public anchor`

That tells us whether RF is carrying real signal or just overlap.

## Practical interpretation

The strongest use of RF_biophys is not “model winner”.

It is:

- a strong ranker for high-prevalence or overlap-heavy families
- a teacher signal for the GA/RL controller
- a top-k booster for discovery mode
- a disagreement partner for wetlab selection

So the expansion is:

**RF_biophys becomes the AUPRC backbone, but the system decides when to trust it.**

## Next implementation step

Implement a `rf_anchor` lane with three outputs:

1. `rf_biophys_norm`
2. `rf_vs_bigmhc_disagreement`
3. `rf_overlap_penalty`

Then let the GA/RL controller learn a routed score:

```text
final_score = mode_router(source, overlap, class, length) -> weighted ensemble
```

That is the most direct way to keep the RF advantage while reducing the leakage risk.
