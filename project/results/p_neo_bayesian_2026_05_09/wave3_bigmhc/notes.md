# BigMHC IM on leakage-stratified ITSNdb — wave3

## Algorithm
- **BigMHC IM** (Albert et al., *Nat Mach Intell* 2023; doi:10.1038/s42256-023-00694-6).
  Repo: https://github.com/KarchinLab/bigmhc (commit at clone time: 2026-05-09).
- Architecture: bidirectional LSTM + multi-head attention + transfer-learned EL→IM head; ensemble of 7 batch-size variants (`bat512` … `bat32768`), each with an `im` subfolder applying transfer learning from the EL backbone.
- Training corpus per the paper: presentation pre-training on **MHCflurry 2.0 mass-spec / IEDB-derived BA+EL data**; immunogenicity transfer-learning on **PRIME 2.0 / Wells / Tesla** assayed peptides plus authors' curated additions.

## Install
- `git clone https://github.com/KarchinLab/bigmhc /tmp/bigmhc` — repo ships pre-trained weights (≈345 MB / model dir × 7).
- No `requirements.txt` in repo; deps already satisfied (torch 2.11.0+cpu, numpy, pandas, psutil).
- Forced CPU mode (`-d cpu`). Inference time on n=319: **~18 s wall** (8 worker jobs, 2 batches).

## Input
- Source: `bundle.tsv` filter `split ∈ {ext_itsndb_main, ext_itsndb_val}`.
- Actual n = **319** (note: task brief said 311; bundle now has 199 main + 120 val = 319).
- Columns supplied: `mhc, pep, tgt` (BigMHC schema).
- Peptide lengths: 9-mer (270), 10-mer (49). All standard 20 AAs.
- 35 unique HLA class-I alleles, all conform to `HLA-A*XX:YY` notation; BigMHC's fuzzy matcher accepted all.

## Skipped
- **0 rows skipped** — all 319 received a `BigMHC_IM` score in `[0,1]`.

## Stratified AUROC (1000-resample percentile bootstrap)

| testset            | in_master | n   | n_pos | AUROC  | 95% CI         |
|--------------------|-----------|-----|-------|--------|----------------|
| ITSNdb_combined    | all       | 319 | 136   | 0.7477 | 0.6962-0.8028  |
| ITSNdb_no_overlap  | False     | 106 |  33   | 0.6260 | 0.5139-0.7308  |
| ITSNdb_in_master   | True      | 213 | 103   | 0.8416 | 0.7869-0.8920  |

**Δ (in_master − no_overlap) = +0.2156** — large gap consistent with a leakage-aware inflation when peptides overlap our master training pool.

## Training-overlap caveat (IMPORTANT)

BigMHC IM was transfer-learned on PRIME 2.0 + Wells + IEDB-derived assays. **ITSNdb itself was a primary external benchmark in the BigMHC paper** (the authors compared BigMHC against MHCflurry/PRIME/NetMHCpan on ITSNdb). They report having explicitly **excluded ITSNdb peptides from training**, but:

1. Many ITSNdb-positive peptides also appear in IEDB / PRIME with the same HLA pairing under different study IDs, so de-duplication is not airtight.
2. Our `in_master = True` slice (n=213) was defined by overlap with **our** master pool (mostly IEDB + Tesla + curated cohorts); BigMHC's training corpus has heavy overlap with the same upstream sources, so `in_master = True` ≈ "BigMHC has likely seen this peptide-or-HLA-context."
3. The 0.842 → 0.626 collapse for the truly-external slice is therefore expected if BigMHC's training set leaks from the same upstream registries that defined our master pool.

We did NOT cross-check BigMHC's bundled training peptide list against ITSNdb at the row level (would require pulling the Mendeley archive at https://doi.org/10.17632/dvmz6pkzvb, ~5 GB). This should be flagged as a limitation in the re-benchmarking paper, and the no_overlap slice (AUROC 0.626) treated as the only honest external estimate.

## Output files
- `predictions.tsv` — 319 rows × 7 cols (peptide, hla, label, in_master, score_bigmhc_im, source, split).
- `auroc_summary.tsv` — 3 rows (combined, no_overlap, in_master).
- `itsndb_input.csv` / `itsndb_bigmhc_im.csv` — intermediate BigMHC I/O kept for audit.
- `build_outputs.py` — reproducible build script.

## Runtime
- Clone: ~30 s (5 GB repo with model weights).
- Inference: ~18 s wall on CPU (laptop-class, ~8 cores).
- Total wall: <2 min including I/O and bootstrap.
