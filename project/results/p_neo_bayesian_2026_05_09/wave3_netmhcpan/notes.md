# NetMHCpan-4.1 — leakage-stratified ITSNdb scoring

## Install path used

**Path 3 (Docker, community image).** Tried in order:

- Path 1 (DTU web service `https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/`): reachable (HTTP 200) but rate-limited form-POST workflow — skipped in favor of Path 3 since Docker hit immediately on the first community image.
- Path 2 (Bioconda): no `mamba`/`conda` available in this environment (.venv-only), skipped.
- Path 3 (Docker): pulled `kevinr9525/netmhcpan-4.1b.linux:latest` (digest `sha256:ef498936...`). Image ships the official DTU NetMHCpan-4.1b binary at `/netMHCpan-4.1/Linux_x86_64/bin/netMHCpan` plus the `data/` dir. The wrapper script `/netMHCpan-4.1/netMHCpan` had a hard-coded path bug (`data//Linux_x86_64/bin/netMHCpan`), so the script invokes the inner binary directly with `NETMHCpan=/netMHCpan-4.1` and `TMPDIR=/tmp` env vars.

Also pulled `danilotat/netmhcpan-minimal:latest` as a fallback — its entrypoint is a VCF-driven launcher (not raw netMHCpan), so it was not used.

Reported version string: `NetMHCpan version 4.1b` (matches the spec target NetMHCpan-4.1).

## HLA format conversion

Bundle format: `HLA-A*02:01`. NetMHCpan format: `HLA-A02:01` (asterisk dropped). Conversion is a single `s.replace("*", "")` in `run_score.py::hla_to_netmhc`. Allele list returned by `netMHCpan -listMHC` (n=11,858) covers all 35 unique HLAs in the ITSNdb subset → **0 rows skipped for unsupported HLA.**

## Bundle filter

Source: `project/results/p_neo_bayesian_2026_05_09/bundle.tsv`

| split             | n   |
|-------------------|-----|
| ext_itsndb_main   | 199 |
| ext_itsndb_val    | 120 |
| **total**         | **319** |

Length filter (8–11mer) and standard-AA filter dropped **0 rows** (all peptides are 9–10mers of canonical AAs).

## Stratification

| stratum             | n   | n_pos | description                          |
|---------------------|-----|-------|--------------------------------------|
| ITSNdb_combined     | 319 | 136   | full external set                    |
| ITSNdb_no_overlap   | 106 |  33   | `in_master == False` (truly external)|
| ITSNdb_in_master    | 213 | 103   | `in_master == True` (training-overlap)|

## Score definition

Primary: `score = -%Rank_EL` (sign-flipped so higher = better). All 319 rows received a non-NaN `%Rank_EL`, so the `%Rank_BA` fallback was unused (n=0).

`-BA` flag was passed at runtime so both EL and BA columns are populated in the netMHCpan output for downstream/audit use.

## Results (1000-bootstrap 95% CI)

| testset             | n   | n_pos | AUROC | 95% CI lo | 95% CI hi |
|---------------------|-----|-------|-------|-----------|-----------|
| ITSNdb_combined     | 319 | 136   | 0.567 | 0.501     | 0.632     |
| **ITSNdb_no_overlap** | 106 |  33   | **0.550** | 0.421     | 0.678     |
| ITSNdb_in_master    | 213 | 103   | 0.571 | 0.493     | 0.649     |

**Δ (in_master − no_overlap) = +0.021**, CIs wide and overlapping. NetMHCpan-4.1 shows a tiny in-master vantage that is *not* significant under the bootstrap on this stratification — but see the training-overlap caveat below.

## Skipped / failures

- `n_skipped_unsupported_hla` = **0**
- `n_dropped_pep_len` = **0**
- `n_dropped_non_std_aa` = **0**
- `n_no_score_after_predict` = **0**
- Per-allele runs that returned 0 parsed rows or non-zero exit codes: **0**
- `n_alleles_run` = **35** (one docker invocation per unique allele)

## Runtime

- Docker pull (one-time): ~30 s
- Allele-list fetch (`-listMHC`): ~3 s
- Per-allele predict (35 dockerized runs): **27.0 s** total (~0.77 s/allele incl. container cold start)
- Total wall (after pull): **33.0 s** on CPU; no GPU.

## Training-overlap caveat (CRITICAL for the leakage paper)

NetMHCpan-4.1 was trained on **>850,000 peptide instances** combining IEDB binding-affinity (BA) data and mass-spectrometry eluted-ligand (EL) data — both single-allele and multi-allele MS studies — covering 170 MHC-I molecules across HLA-A/B/C/E and non-human species. ITSNdb itself was assembled from public IEDB / mass-spec sources, so:

1. Many ITSNdb peptides are **highly likely** to be in NetMHCpan-4.1's own training corpus, regardless of our `in_master` flag (which tracks overlap with the *in-house* bayesian model's master pool, NOT NetMHCpan's training set).
2. The `ITSNdb_no_overlap = 0.550` AUROC on the truly-external-to-our-master split is therefore **still upper-bounded** for NetMHCpan, since those peptides may have been seen by NetMHCpan during its own training. A "true" no-overlap evaluation would require deduplicating against the NetMHCpan-4.1 IEDB / MS training cuts (training-data list is published under the DTU release notes).
3. The near-zero Δ (0.021) here can mean either (a) NetMHCpan-4.1 has saturated this peptide population so the in_master flag carries no extra signal, or (b) NetMHCpan-4.1 has memorized large portions of *both* strata, which would also collapse the contrast. Distinguishing requires the deduplication step in (2).
4. The relevant comparison for the leakage-aware paper is the **algorithm × stratum interaction** across the 6 algorithms — NetMHCpan's flat profile (Δ=+0.021) versus algorithms with steeper Δ (training-set leakage signature) is itself the finding, regardless of whether NetMHCpan's absolute AUROC reflects real generalization or memorization of the IEDB superset.

The mid-50s combined-set AUROC (0.567) is on the low side for a state-of-the-art predictor on benchmark MHC-I data, suggesting ITSNdb is by design weighted toward **immunogenicity** (peptide → T-cell response) rather than **presentation** (peptide → MHC binding), where NetMHCpan is calibrated. That gap between presentation and immunogenicity is the headline ITSNdb finding (Bonsack 2019, Schaap-Johansen 2021) and is consistent with what we see here.

## Scripts + outputs

- `run_score.py` — entry point (one docker call per unique HLA, parses stdout via regex, joins back to bundle by peptide).
- `predictions.tsv` — 319 rows: `peptide | hla | label | in_master | score_netmhcpan | source | split`. `score_netmhcpan = -%Rank_EL`.
- `auroc_summary.tsv` — 3 rows (combined / no_overlap / in_master).
- `_run_meta.json` — runtime + skip stats + install path.
