# Wave 9 — three additional 2024-2025 peptide-MHC immunogenicity / stability tools

**Goal**: extend the leakage-stratified ITSNdb forest from 6 algorithms to 9, adding more recent tools to strengthen the leakage-aware re-benchmarking paper.

**Test set** (frozen): `bundle.tsv` filtered to `split ∈ {ext_itsndb_main, ext_itsndb_val}`, n = 319 (213 in_master / 106 no_overlap).

**Compute**: local CPU. Runtime per tool: NetMHCstabpan 69s · T-SCAPE 118s · MHCnuggets ~3 min (per-allele loop). Total wall time well under the 90 min budget.

## Tools installed and run

| Algorithm | Source | Year | Score column | n_input | n_scored | n_skipped (HLA) | n_skipped (runtime) |
|-----------|--------|------|--------------|--------:|---------:|----------------:|--------------------:|
| **MHCnuggets-2** | KarchinLab/mhcnuggets (Shao 2020, maintenance through 2024) | 2020/2024 | -IC50 | 319 | 317 | 2 (B*56:01, B*41:02) | 0 (after A*01:01 fix) |
| **T-SCAPE / TITANiAN** | seoklab/TITANiAN (Kim et al, [bioRxiv 2025.05.11](https://www.biorxiv.org/content/10.1101/2025.05.11.653308v1)) — `pmhc_im_neo` cancer-immunotherapy variant | 2025 | model output (sigmoid, P(immunogenic)) | 319 | 319 | 0 | 0 |
| **NetMHCstabpan-1.0** | Docker community image `kevinr9525/netmhcstabpan-1.0a.linux:latest`; DTU Rasmussen 2016 | 2016 (still current) | -%Rank_Stab | 319 | 319 | 0 | 0 |

A2 fallback: original priority list was Anthem (2020) → MHCSeqNet (2019). Both Anthem candidate URLs (mleming/anthem, zjupgx/anthem, Mengmeng-NING/Anthem, ZhuLab-Fudan/Anthem) returned 404. Substituted **T-SCAPE / TITANiAN (Kim 2025 bioRxiv, Seoklab)**, which is explicitly a **2025** T-cell epitope immunogenicity predictor with a "cancer immunotherapy" model variant — closer to the paper's framing than Anthem (presentation) or MHCSeqNet (binding) would have been.

A3 NetMHCstabpan installed via Docker community image (kevinr9525/netmhcstabpan-1.0a.linux). Same-pattern install as the wave3 NetMHCpan-4.1b run (also via kevinr9525/...). Bioconda was not used because no `conda`/`mamba` is on this venv-only host.

## 9-algorithm forest table (rows ordered by no_overlap AUROC descending)

| Algorithm | combined AUROC | in_master AUROC | **no_overlap AUROC** | no_overlap 95% CI |
|-----------|---------------:|----------------:|---------------------:|:------------------|
| **MHCflurry** (prior) | 0.637 | 0.642 | **0.668** | [0.536, 0.787] |
| **BigMHC** (prior) | 0.748 | 0.842 | **0.626** | [0.514, 0.731] |
| **T-SCAPE / TITANiAN** *(new)* | 0.602 | 0.640 | **0.605** | [0.478, 0.739] |
| **DeepImmuno** (prior) | 0.702 | 0.803 | **0.579** | [0.464, 0.698] |
| **PRIME** (prior) | 0.567 | 0.564 | **0.572** | [0.452, 0.686] |
| **TransPHLA** (prior) | 0.580 | 0.583 | **0.555** | [0.440, 0.673] |
| **NetMHCpan** (prior) | 0.567 | 0.571 | **0.550** | [0.421, 0.678] |
| **MHCnuggets-2** *(new)* | 0.647 | 0.710 | **0.521** | [0.405, 0.639] |
| **NetMHCstabpan** *(new)* | 0.509 | 0.534 | **0.458** | [0.341, 0.572] |

**Δ (in_master − no_overlap)** — magnitude of the leakage signature per algorithm:
- DeepImmuno +0.224, BigMHC +0.216, MHCnuggets +0.189 (largest leakage signatures)
- T-SCAPE +0.035, NetMHCstabpan +0.076, PRIME −0.008, TransPHLA +0.028, NetMHCpan +0.021 (small / null)
- MHCflurry −0.026 (no_overlap actually beats in_master)

The new MHCnuggets-2 column joins the **leakage-suspect cluster** (DeepImmuno / BigMHC) with one of the largest in_master vs. no_overlap deltas (+0.189). T-SCAPE behaves more like the **flat-profile cluster** (PRIME, NetMHCpan), with a small Δ around the null. NetMHCstabpan is the only sub-50 AUROC entry on the no_overlap stratum, consistent with its design target (pMHC stability) being one inferential step further from immunogenicity than the binding-affinity tools.

## Best of the 3 new vs. MHCflurry 0.668

- T-SCAPE 0.605 — **−0.063** below MHCflurry no_overlap; best of the three new tools
- MHCnuggets-2 0.521 — −0.147
- NetMHCstabpan 0.458 — −0.210

None of the three new entries surpass the prior leader (MHCflurry no_overlap 0.668). The 95% CIs for T-SCAPE [0.478, 0.739] vs. MHCflurry [0.536, 0.787] **overlap heavily**, so on this n=106 stratum we cannot reject the null that T-SCAPE matches MHCflurry's clean-set performance. The interesting finding is the **rank reordering between strata**: BigMHC and DeepImmuno drop substantially when switching from in_master → no_overlap (BigMHC 0.842 → 0.626, DeepImmuno 0.803 → 0.579), while MHCflurry barely moves (0.642 → 0.668) and PRIME / NetMHCpan / TransPHLA stay flat in the mid-50s.

## Tools that failed and why

- **Anthem** (Mei 2020 Nat Commun) — none of the suggested GitHub paths resolved (mleming/anthem, zjupgx/anthem, Mengmeng-NING/Anthem, ZhuLab-Fudan/Anthem all 404). Substituted with T-SCAPE.
- **MHCSeqNet** (Phloyphisut 2019) — repo present but ships TF 1.x era code (`from keras` vs. modern `from tensorflow.keras`); was not attempted in favor of the working T-SCAPE pipeline.
- **NetTCR-2.2** (Montemurro 2024) — repo cloned, but the model is for **TCR-pMHC** binding (requires CDR3β input), not pure pMHC immunogenicity, so it is the wrong scope for the leakage-aware paper. Skipped.

The HLA-A*01:01 batch (42 ITSNdb peptides) crashed inside MHCnuggets's `rank_output=True` branch with `np.float32(5447.98)` numerical comparison error against a precomputed human-proteome IC50 ranking pickle. **Fix**: re-ran HLA-A*01:01 alone with `rank_output=False`, recovered the 42 ic50 values, sign-flipped to `-IC50` like the rest. n_skipped final = 0 runtime, 2 unsupported HLA (B*56:01, B*41:02).

## Training-overlap caveat (for the leakage-aware paper)

All three new tools are **trained on IEDB-derived data** (MHCnuggets explicitly so; T-SCAPE includes IEDB and BigMHC pre-training; NetMHCstabpan is trained on IEDB + DTU stability assays). Our `in_master` flag tracks overlap with the **in-house bayesian model's master pool**, not these tools' training corpora. So the in_master AUROC inflation seen here is conservative — true overlap with each tool's own training set may be higher, which would further inflate the in_master column relative to the no_overlap column. The Δ (in_master − no_overlap) per-algorithm pattern remains the **algorithm × stratum interaction** finding the paper relies on, regardless.

## Files

- `predictions_mhcnuggets.tsv`, `predictions_tscape.tsv`, `predictions_netmhcstabpan.tsv` — per-peptide scores
- `auroc_mhcnuggets.tsv`, `auroc_tscape.tsv`, `auroc_netmhcstabpan.tsv` — 3-row AUROC summaries (combined / no_overlap / in_master) with 1000-bootstrap 95% CI
- `wave9_combined_results.tsv` — long-form 27-row table (9 algorithms × 3 strata) merging wave3 + wave9
- `wave9_forest_pivot.tsv` — wide pivot for plotting
- `wave9_summary_table.tsv` — 9-row × 3-column AUROC summary
- `fig_wave9_9algorithm_forest.png/pdf` — 9-algorithm horizontal-error-bar forest plot, rows ordered by no_overlap descending, three colored series (in_master / combined / no_overlap)
- `_meta_mhcnuggets.json`, `_meta_tscape.json`, `_meta_netmhcstabpan.json` — runtime + skip stats
- `run_mhcnuggets.py`, `run_tscape.py`, `run_netmhcstabpan.py`, `fix_mhcnuggets_a0101.py`, `merge_and_plot.py` — reproducible scripts
- `mhcnuggets.log`, `tscape.log`, `stabpan.log` — per-tool runtime logs

## Bootstrap

1000 resamples, seed=42, sklearn `roc_auc_score` (skips bootstraps that lose class balance). Same boot configuration as wave3.
