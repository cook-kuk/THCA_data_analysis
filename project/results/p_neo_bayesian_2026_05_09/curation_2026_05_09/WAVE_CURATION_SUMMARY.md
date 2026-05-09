# Neo Wave Curation Summary - 2026-05-09

## Process cleanup

Stopped stale `p_neo_bayesian_2026_05_09` wave jobs:

- `wave3_algorithm_sweep` PRIME tail jobs for A0217 / B4102, already documented as timeout/skipped in `WAVE3_REPORT.md`.
- `wave6/track6b_simple_reps.py` and associated Claude log watchers.
- `wave5a/train_student_distill.py`, incomplete distillation run and watchers.
- Remote `neo_wave4a` log tail watcher.

After cleanup, no `p_neo_bayesian_2026_05_09` process remained active.

## Fast curation outputs

- `quick_no_overlap_ranking.tsv`: primary fast ranking on ITSNdb no-overlap.
- `quick_method_metrics.tsv`: combined / in-master / no-overlap metrics, no bootstrap.
- `quick_wave8_artifact_stress.tsv`: Wave 8 exact-overlap stress audit.
- `quick_wave_inventory.tsv`: wave directory inventory.
- `curated_predictions_itsndb_long.tsv`: row-level prediction cache.
- `curated_predictions_itsndb_wide.tsv`: wide row-level prediction cache.

The full bootstrap script was patched to be lighter but still ran too slowly for cleanup mode; quick outputs are the working source of truth for this pass.

## Primary no-overlap ranking

Top methods on `ITSNdb_no_overlap` from row-level predictions:

| Rank | Method | AUROC | n | Read |
|---:|---|---:|---:|---|
| 1 | Wave8_TCR_SelfSim_full | 0.735 | 106 | Candidate only; exact-reference sensitive |
| 2 | Wave8_TCR_SelfSim_no_exact | 0.732 | 106 | Candidate only; still TCR-reference sensitive |
| 3 | Wave8_TCR_motif_only | 0.722 | 106 | Driven by TCR reference hits |
| 4 | MHCflurry | 0.668 | 106 | Public pretrained comparator; training-overlap unresolved |
| 5 | Structure_LR | 0.653 | 103 | Best custom robust baseline with known local training split |
| 6 | Wave8_TCR_only | 0.644 | 106 | Drops after exact-hit removal |
| 7 | Stack_LR_inmaster_E3a | 0.631 | 106 | Not fair primary, trained on in-master |
| 8 | BigMHC_IM | 0.626 | 106 | Off-the-shelf, leakage-inflated on in-master |
| 9 | Wave8_SelfSim_full | 0.625 | 106 | Modest, label-artifact risk |
| 10 | Stack_mean_E1 | 0.623 | 106 | Does not beat best single tool |

## Wave 8 artifact audit

Wave 8 looks best only before removing exact TCR-reference hits.

| Subset | n | n_pos | TCR exact hits | self exact hits | Wave8 full | Wave8 no-exact | MHCflurry | Structure_LR |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| no_overlap_all | 106 | 33 | 13 | 1 | 0.735 | 0.732 | 0.668 | 0.653 |
| no_overlap_no_same_hla_tcr_exact | 93 | 22 | 0 | 1 | 0.633 | 0.634 | 0.618 | 0.678 |
| no_overlap_no_tcr_or_self_exact | 92 | 22 | 0 | 0 | 0.642 | 0.637 | 0.621 | 0.679 |

Conclusion: Wave 8 is not paper-ready as a new SOTA method. Its headline gain is heavily dependent on 13 exact TCR-reference hits. After strict removal, `Structure_LR` is the cleaner local-training leader.

## MHCflurry training-overlap caveat

`MHCflurry` must not be called a "clean" external baseline here. Local model provenance shows:

- Installed package: `mhcflurry 2.2.1`.
- Downloaded presentation model URL: `https://github.com/openvax/mhcflurry/releases/download/pre-2.0/models_class1_presentation.20200611.tar.bz2`.
- Model `info.txt`: trained on 2020-06-11 with package `mhcflurry 1.7.0`.
- Official MHCflurry docs describe pretrained models as fit to mass-spec-identified MHC-I ligands and peptide/MHC affinity measurements deposited in IEDB plus other sources; the presentation score combines binding-affinity and antigen-processing predictors.

Therefore the `MHCflurry` no-overlap AUROC = 0.668 is only no-overlap relative to **our master table**, not relative to MHCflurry's hidden/public training corpus. It should be reported as a public-pretrained upper-bound comparator, not as definitive external generalization.

## Wave dispositions

| Wave | Disposition | Keep for paper? | Reason |
|---|---|---|---|
| wave1 / root | archive negative | yes, as failure baseline | ESM2-Bayesian no-overlap ~0.41 but in-master ~0.94; clean leakage-failure example |
| wave2 | keep partial | yes | Structure_LR and calibration/OOD panels are useful; quantum is negative |
| wave25 | feature source | supplement only | PWM features feed stacking/ablation |
| wave3_algorithm_sweep + wave3_* | keep | yes, with training-corpus caveat | Core off-the-shelf comparison; MHCflurry/BigMHC/others may have unobserved public-training overlap |
| wave4a | archive negative | supplement only | DG variants inflate in-master and fail no-overlap |
| wave4a_lora_fix | no claim | no | Script only |
| wave4b | keep negative/control | supplement | kNN/TENT/conformal do not solve no-overlap; conformal undercoverage is useful |
| wave4c | archive negative | supplement | Stacking does not beat MHCflurry |
| wave5a | stopped incomplete | no | No completed result table |
| wave5b | keep ablation | supplement | Multitask improves slightly but remains below chance/no-overlap |
| wave5c | no claim | no | Scripts only |
| wave6 | stopped partial | notes only | Partial negative encoder-hierarchy work |
| wave7 | archive negative | supplement | Synthesis classifier fails to beat MHCflurry/Structure_LR |
| wave8 | candidate retest | not yet | Strong apparent signal, exact-overlap sensitive |
| wave9 | keep | yes, with training-corpus caveat | Adds T-SCAPE/MHCnuggets/stabpan; public-training overlap unresolved |
| wave10 | canonical | yes | Best existing mega comparison and figures |
| wave11 | keep supplement | yes, with low-power flags | Expanded test sets; many all-positive or low-power bundles |
| wave12 / wave15 / wave16 | absent on disk | no | Claude UI labels only; no local output found |

## Paper-grade read

The honest paper axis is not "new architecture beats the field." The clean result is:

1. Methods trained on the in-house/public master table can show huge in-master AUROC and collapse on no-overlap.
2. Some public pretrained tools (`MHCflurry`, `PRIME`, `NetMHCpan`, `TransPHLA`) have smaller inflation gaps against **our** overlap flag, but their own training corpora are not fully audited.
3. `Structure_LR` is the best custom robust baseline with known local training split; `MHCflurry` is a public-pretrained upper-bound comparator, not a clean external baseline.
4. Wave 8 TCR/self-similarity is a retest candidate, not a claim, until exact-reference leakage is controlled on independent external sets.
