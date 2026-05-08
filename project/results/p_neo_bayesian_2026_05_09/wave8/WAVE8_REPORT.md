# Wave 8 — TCR / 3D / Self-similarity feature engineering

Three modalities probed as standalone classifiers, beyond what Wave 7's
combined feature stack (structural + MHCflurry + PWM + HLA-pseudo-PCA) covers.
Goal: identify novel signal that is not captured by current strengths and that
could be added to a Wave 9 combined model.

## Tracks

| Track | Modality | Status | n features | n peptides scored |
|-------|----------|--------|-----------|-------------------|
| 8A | TCR engagement proxy | DONE | 3 | 2666 / 2666 unique |
| 8B | ESMFold pseudo-complex 3D | SKIPPED | 0 | 0 |
| 8C | Self-similarity to human SwissProt | DONE | 4 | 2666 / 2666 unique |

Track 8B was skipped per spec: ESMFold v1 fails to load on this host
(`ModuleNotFoundError: omegaconf`); no GPU available locally
(`torch 2.11.0+cpu`, CUDA=False). Folding 200 peptides at ~30 s/each on CPU
would already exceed budget, and the 3D modality risks duplicating Wave 2's
Structure_LR (BLOSUM + Kyte-Doolittle + Chou-Fasman + per-position PWM).
Status JSON: `track8b_status.json`.

## Standalone AUROC on `ITSNdb_no_overlap` (n=106, 33 positives)

| Feature block | n_feat | AUROC | 95% CI |
|---|---|---|---|
| **TCR_only** | 3 | **0.644** | 0.512 – 0.765 |
| **SelfSim_only** | 4 | **0.625** | 0.503 – 0.732 |
| **TCR + SelfSim (Wave 8 union)** | 7 | **0.735** | 0.618 – 0.839 |

Reference comparator (Wave 5b head A `multitask_supervision`):
- ITSNdb_no_overlap AUROC = **0.448** (CI 0.335 – 0.570) — i.e. *worse than chance*

Both Wave 8 modalities therefore beat the previously-best single head on the
hardest split (`no_overlap`, the rows where the reference peptide does NOT
appear in the master). The combined Wave 8 union (0.735) is a substantial
uplift, with a CI floor of 0.618 — well clear of the 0.55 publishability
threshold defined in the spec.

## Standalone AUROC on auxiliary splits

| Feature block | ITSNdb_in_master AUROC | ITSNdb_combined AUROC |
|---|---|---|
| TCR_only | 0.586 | 0.612 |
| SelfSim_only | 0.464 | 0.503 |
| TCR + SelfSim | 0.586 | 0.629 |

Note the **anti-correlation between in_master and no_overlap for SelfSim**:
self-similarity is a useful signal *only on truly held-out peptides*. On
in-master rows (where positive labels skew toward known shared antigens),
self-similarity inverts. This is consistent with curated training pools
containing cancer-testis / shared-self antigens whose self-similarity is
high *and* whose label is positive — a well-known confound. The signal is
clean on `no_overlap`, where the eval set is freshly novel.

## Self-similarity sign — verdict

Mixed. From the SelfSim_only LR coefficients (trained on TRAIN pool, n=2396):

| Feature | Coef | Direction |
|---|---|---|
| self_exact_match | +1.92 | exact-self → labeled MORE immunogenic (training-set artifact: shared/CT antigens) |
| self_hamming1_count_log | -0.64 | many close-self-neighbors → LESS immunogenic (tolerance-consistent) |
| self_hamming2_count_log | -0.02 | flat |
| self_blosum_max | -0.18 | higher BLOSUM-to-self → LESS immunogenic (tolerance-consistent) |

The two distance-based features (Hamming-1 count, BLOSUM-max-to-self) sign with
classical self-tolerance theory: **peptides surrounded by many close self-paralogs are less likely to elicit a response**. The exact-match feature carries an
opposite, label-noise-driven sign because the curated training pool is enriched
for shared/cancer-testis self-antigens (Melan-A, NY-ESO-1, etc.) which ARE
self-similar but ARE labeled positive.

For Wave 9 production use I'd recommend dropping `self_exact_match` and keeping
only the distance features, but the eval here uses the full block as designed
for honest reporting.

## TCR data sources used

- **VDJdb slim** (`/data/neoantigen_vaccine_hub/data_raw/tcr/vdjdb.slim.txt`,
  75,281 rows). Filtered to *Homo sapiens* MHC-I with valid 8–12mer peptide:
  **1,031 unique (peptide, HLA) pairs**.
- **McPAS-TCR** (`/data/neoantigen_vaccine_hub/data_raw/tcr/mcpas/McPAS-TCR.csv`,
  40,778 rows). Same filtering: **381 unique (peptide, HLA) pairs**.
- Combined: **1,321 unique (peptide, HLA) reference pairs across 82 unique
  HLA alleles** + 38 supertypes.
- Top supertype counts: HLA-A02 (810), HLA-B07 (121), HLA-A01 (121),
  HLA-A24 (74). HLA-A02 dominance reflects the underlying TCR-pMHC literature.

Features:
- `tcr_motif_score`: max of mean(Jaccard k=3,4,5) to any same-HLA reference
  peptide (0.5× discount when falling back to global pool for unseen HLAs).
  62.3% of bundle rows have score > 0.
- `tcr_motif_count_log`: log1p of length-matched Hamming-≤-2 reference hits.
- `tcr_class_score`: per-HLA-supertype TCR engagement frequency
  (log1p-normalized 0–1).

## ESMFold install status — Track 8B

Skipped — see `track8b_status.json`. Reason summary:
1. `esm.pretrained.esmfold_v1()` raises `ModuleNotFoundError: omegaconf`
   (openfold-style dep stack not installed locally).
2. CPU-only torch (no CUDA on this host).
3. RunPod A6000 pod is provisioned (`thca-spark-dm-a6000-v4`) but not used
   for this wave to keep all of Wave 8 within the 90-min local budget.
4. The structural modality is partially covered already by Wave 2's
   Structure_LR (BLOSUM62 + Kyte-Doolittle + Chou-Fasman + position PWM).

If revisiting later, a viable path: pull a single openfold-compatible Docker
image on the RTX A6000 pod, fold the 200 highest-disagreement no-overlap
peptides (~100 minutes wall, ~$0.55), extract `interface_contact_density`,
`anchor_pLDDT_p2/p9`, `mean_pLDDT`, peptide radius-of-gyration, plug into the
same Wave 8 standalone harness.

## Verdict — which modality to add to Wave 7 / Wave 9

1. **TCR_only is the cleanest single addition** — AUROC 0.644 on no_overlap with
   a CI bottom of 0.512 (touching but not crossing 0.5). Captured signal is
   orthogonal to MHCflurry presentation and to PWM (TCR engagement is
   downstream of presentation). Recommended as default Wave 9 addition.

2. **SelfSim is the second-cleanest** — AUROC 0.625 on no_overlap, CI floor 0.503.
   Inverts on in_master (label-noise-driven), so the *honest* use case is
   regularized into the held-out classifier only. Recommended addition with
   `self_exact_match` dropped to remove the dominant noise driver.

3. **The TCR+SelfSim union is the strongest single block from Wave 8** at
   AUROC 0.735 (no_overlap), CI floor 0.618. This is the form to merge
   into the Wave 9 combined matrix.

`wave8_for_wave7.tsv` is keyed on (peptide, hla) — drop-in for the Wave 7
combined feature pipeline.

## Files

- `tcr_features.tsv` (Track 8A output)
- `self_similarity_features.tsv` (Track 8C output)
- `track8b_status.json` (Track 8B skip rationale)
- `wave8_combined_features.tsv` / `wave8_for_wave7.tsv` (merged matrix)
- `wave8_standalone_results.tsv` (per-block × per-split AUROC + CI)
- `wave8_eval_meta.json` (LR coefficients per block)
- `fig_wave8_modality_uplift.png` / `.pdf`
- `track8a_tcr_features.py`, `track8c_self_similarity.py`, `eval_standalone.py`
- `human_swissprot.fa` (13.7 MB cached from UniProt; 20,431 SwissProt-reviewed
  human proteins; 10.3M – 10.5M unique k-mers per length 8–12)

## Honest limitations

- TCR motif score for HLA-A02 supertype dominates (810 of 1321 references);
  performance on rarer supertypes (B07, B35, B57) rests on smaller reference
  pools. `tcr_class_score` partially encodes this prior.
- Self-similarity is a *whole-proteome* signal; it does not distinguish
  promiscuous-housekeeping peptides from tissue-restricted self-peptides.
  Reviewer Q: a TPM-weighted version (only thymically-expressed proteins for
  central tolerance) would be cleaner. Not done in this pass.
- All AUROCs are on a single test split (ITSNdb). Cross-cohort robustness
  (Venus, LOSO) was not re-run for Wave 8 — that belongs to a Wave 9 combine
  step that already runs the full LOSO grid for Wave 7's stack.
