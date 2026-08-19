# CROSS-Neo v1 Practical Handoff

Date: 2026-05-10

## Verdict

Decision remains HOLD for manuscript-level promotion because public comparator
training overlap is still unresolved and source-heldout behavior is uneven.

The strongest practical v1 signal is now:

- locked HLA-stratified main candidate:
  `v1_gated_moe / prespecified_equal_weight_C_QK_no_anchor`
- bounded hard-decoy repair:
  `hard_decoy_rule_aux_C_QK_no_anchor` and
  `hard_decoy_nested_meta_C_QK_aux`

No external validation or quantum advantage should be claimed.

## Hard-Decoy Repair

The first decoy/focal branch was too standalone and collapsed on locked splits.
The repaired branch generates train-fold hard decoys and uses them only as a
bounded auxiliary signal on top of C+QK.

Key locked results:

- HLA-stratified: AUPRC 0.521, top10 0.600
- near-peptide holdout: AUPRC 0.539, top10 0.600 with nested meta
- exact peptide-HLA holdout: AUPRC 0.587, top10 0.700
- repeated internal: AUPRC 0.454, top10 0.600

Key source-stress results:

- NEPdb: AUPRC 0.426, top10 0.700
- TESLA_mmc4: AUPRC 0.119, top10 0.300
- TESLA_mmc7: AUPRC 0.083, top10 0.100

Interpretation: hard decoys now help top-k triage, especially exact/near and
NEPdb/TESLA_mmc4 stress, but this is still not an external-valid model.

## Business Package

Directory:

`project/results/cross_neo_v1/business_package/`

Files:

- `candidate_input_template.tsv`
- `candidate_input_schema.tsv`
- `sequence_only_business_triage_demo.tsv`
- `internal_locked_demo_priority_queue.tsv`
- `commercialization_readiness_checklist.tsv`
- `BUSINESS_MODEL_CARD.md`

Intended use:

1. Use the template for patient candidate bags.
2. Run sequence-only triage only as intake screening.
3. Run full CROSS-Neo feature generation for serious candidates.
4. Report per-patient top-k with missing-context abstention and self-risk.
5. Keep wet-lab validation and public-overlap audit as business gates.

## Re-run Commands

Full v1:

```bash
bash scripts/run_cross_neo_v1_all.sh
```

Hard-decoy plus business package only:

```bash
bash scripts/run_cross_neo_v1_hard_decoy_business_only.sh
```

## Transfer Package

Zip:

`project/results/transfer_packages/cross_neo_v1_hard_decoy_business_2026_05_10.zip`

