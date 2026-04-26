# v16 Xenograft Design for Sacituzumab Govitecan Follow-up

## Objective

Translate the in silico `v16_sg_simulation` findings into a first-pass in vivo
efficacy experiment testing whether **BRAF-like / TACSTD2-high thyroid tumors**
show greater response to **sacituzumab govitecan (SG)** than RAS-like /
TACSTD2-lower comparators.

This document is an experimental design brief, not an approved protocol.

## In Silico Rationale

- `TACSTD2` is strongly higher in the BRAF-like contrast (`log2FC = +4.87` in the
  virtual-treatment summary), supporting greater ADC internalization potential.
- The MAPK-output proxy (`DUSP4/DUSP5/DUSP6`) is also higher in the BRAF-like
  state (`mean log2FC = +1.64`), consistent with a more proliferative /
  signaling-active background that could accentuate SN-38 vulnerability.
- The heuristic composite from
  [virtual_treatment_summary.json](/home/seungho/personal/THCA_data_analysis/project/results/v16_sg_simulation/virtual_treatment_summary.json)
  yields a **51.51x predicted BRAF-like fold sensitivity** versus RAS-like.
- Structural and sequence work products supporting target plausibility were
  generated under
  [project/results/v16_sg_simulation](/home/seungho/personal/THCA_data_analysis/project/results/v16_sg_simulation).

## Primary Hypothesis

SG will produce greater tumor growth inhibition in **TROP2/TACSTD2-high,
BRAF-like thyroid cancer xenografts** than in matched RAS-like or
TROP2-lower comparators.

## Model Selection

### Arm A: BRAF-like / SG-sensitive candidate

- Choose 1 thyroid cancer model with:
- `BRAF V600E` genotype.
- High `TACSTD2` / TROP2 mRNA and confirmable surface protein expression.
- Adequate in vitro growth and xenograft take rate.

Suggested screening rule before implantation:

- RT-qPCR for `TACSTD2`.
- Flow cytometry or IHC for TROP2 membrane signal.
- Optional confirmation of `DUSP4`, `DUSP5`, `DUSP6` as MAPK-output markers.

### Arm B: Comparator / lower-sensitivity candidate

- Choose 1 thyroid cancer model with:
- RAS-like genotype or phenotype.
- Lower TROP2/TACSTD2 abundance than Arm A.
- Similar engraftment robustness if possible.

### Optional Arm C: Target-dependence control

- BRAF-like model with `TACSTD2` knockdown or CRISPRi.
- This is the cleanest mechanistic control if resources permit.

## Minimal Study Design

### Recommended first pass

- Two models: one BRAF-like, one comparator.
- Four treatment groups per model:
- `vehicle`
- `naked anti-TROP2 antibody control` if available
- `irinotecan or SN-38-linked comparator` if feasible
- `sacituzumab govitecan`

If budget is tight, the absolute minimum is:

- `vehicle`
- `SG`

within each model.

### Animal allocation

- Start with `n = 8-10` mice per group for the pilot efficacy readout.
- Randomize once tumors reach approximately `100-150 mm^3`.
- Balance groups by baseline tumor volume.

## Dosing Framework

Use the vendor / institutional precedent dose as the starting anchor rather than
the heuristic model. The in silico work does **not** justify inventing a novel
dose.

Operational recommendations:

- Follow published SG mouse schedules from solid-tumor xenograft literature.
- Record body weight at each dosing event.
- Prespecify dose-hold rules for weight loss and clinical toxicity.

## Endpoints

### Primary endpoint

- Tumor growth inhibition versus vehicle.

### Key secondary endpoints

- Time to volume doubling.
- Objective regression fraction.
- Body-weight change / tolerability.
- Terminal TROP2 IHC and `TACSTD2` qPCR.
- DNA-damage pharmacodynamic markers:
- `gamma-H2AX`
- cleaved caspase-3
- Ki-67 reduction

### Mechanistic secondary endpoints

- Compare pre/post-treatment `TACSTD2`, `DUSP4`, `DUSP5`, `DUSP6`.
- Assess whether SG response tracks baseline TROP2 abundance more strongly than
  genotype alone.

## Statistical Plan

- Analyze each model separately first.
- Then test a model-by-treatment interaction across pooled data.
- Prespecify:
- exclusion criteria for failed engraftment
- censoring rules
- humane endpoint criteria

Recommended analyses:

- Mixed-effects model or repeated-measures ANOVA for longitudinal tumor volume.
- AUC or day-21 tumor-volume comparison as a simple confirmatory endpoint.
- Effect sizes with confidence intervals, not p-values alone.

## Go / No-Go Criteria

Advance to expanded validation if all of the following are observed:

- SG beats vehicle in the BRAF-like model.
- The BRAF-like model shows stronger response than the comparator model.
- Baseline TROP2/TACSTD2 level tracks with efficacy directionally.
- Tolerability is acceptable under the selected schedule.

Stop or redesign if:

- Both models respond equally despite clear TROP2 separation.
- Neither model responds despite confirmed target expression.
- Toxicity prevents interpretable exposure.

## Main Caveats

- The `51.51x` selectivity estimate is a **heuristic ranking score**, not a real
  efficacy multiplier.
- The virtual model lacks explicit repair-capacity, payload PK, and linker
  stability terms.
- The current work does not include true perturbation modeling (`GEARS`) or
  antibody-antigen structural validation with AF3 / wet-lab binding assays.

## Recommended Pre-Xenograft Checklist

- Confirm TROP2 surface expression in candidate lines.
- Confirm `BRAF-like` versus comparator transcriptional state.
- Verify SG access, formulation, and institution-approved handling workflow.
- Lock a simple statistical analysis plan before dosing.
- Bank baseline RNA / FFPE material for post hoc biomarker work.

## Related Outputs

- [v16_index.json](/home/seungho/personal/THCA_data_analysis/project/results/v16_sg_simulation/v16_index.json)
- [virtual_treatment_summary.json](/home/seungho/personal/THCA_data_analysis/project/results/v16_sg_simulation/virtual_treatment_summary.json)
- [braf_tacstd2_pathway.json](/home/seungho/personal/THCA_data_analysis/project/results/v16_sg_simulation/braf_tacstd2_pathway.json)
- [ab_trop2_summary.json](/home/seungho/personal/THCA_data_analysis/project/results/v16_sg_simulation/ab_trop2/ab_trop2_summary.json)
- [sn38_topoi_summary.json](/home/seungho/personal/THCA_data_analysis/project/results/v16_sg_simulation/sn38_topoi_docking/sn38_topoi_summary.json)
