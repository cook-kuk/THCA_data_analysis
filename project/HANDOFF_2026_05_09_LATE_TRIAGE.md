# Handoff addendum - 2026-05-09 late triage

This addendum summarizes the worktree after `project/HANDOFF_2026_05_09.md`.
It is an operating map only: no voice-protected manuscript prose.

## Immediate State

- Branch: `paper9-perturbation-extension-20260506`.
- Latest local commit: `e20eae4` (`wave 7 - synthesis classifier: combined-feature heads do not beat MHCflurry alone`).
- Worktree is dirty: about 299 changed or untracked entries.
- Storage check: `project/results` and `project/data` are correctly bind-mounted to `/data`; `/data` is healthy, but root `/` is high at 88% used. Avoid new heavy outputs outside bind-mounted paths.
- Do not use `git add -A`. Stage by topic.

## Track Decisions

| Track | Current Read | Disposition |
|---|---|---|
| Paper 1 Fig 8 / deconvolution | v2/v3/v5 and v13/v14 are usable; v15A K2 and v15C PRISM/DepMap are reserve; spatial v15B-v18 did not validate the positive mechanism chain. | Use v13/v14 as positive mechanism support. Keep v15A/v15C as reviewer reserve. Use v16-v18 only as caveat/autopsy. Stop spatial extensions. |
| Paper 2 image-DM1 | ViT-L/early CLAM RAS-like AUC=1.000 is a TSS batch-effect artifact; LOTO RAS-like AUC drops to 0.308. UNI-final is stronger (LOTO overall 0.852, RAS-like 0.718) but the TSS-only within-RAS shortcut remains. | Demote the old RAS-like-perfect claim. UNI can stay as a cautious pilot, but no center-generalized subgroup claim before K2 H&E or another external WSI cohort. |
| Paper 3 ICI | Evidence supports ICI vulnerability/readiness, not thyroid ICI response prediction. HLA-I carries the modest response signal externally; HLA-II is biology/context, not response predictor. | Keep Track A frozen. Do not start Track B without explicit author command after Paper 1 bioRxiv and Paper 2 closure. |
| HLA two-paper package | Paper 4 Graves/AITD HLA is the stronger manuscript-ready package. Paper 2 HT-overlap HLA should be transcriptomic AP/TLS/IFN biology, not allele association. | Keep separate. Do not claim DPB1*05:01 as HT-PTC finding or DRB1*04:05 as validated. Korean adult GD NGS cohort is the NComm blocker. |
| Neoantigen / cancer vaccine | Wave10 mega-table says MHCflurry and Structure_LR lead on clean ITSNdb no-overlap; many custom models inflate on in-master overlap. CROSS-Neo v1 remains HOLD. | Frame as leakage-stratified benchmark / uncertainty / abstention, not external clinical predictor or quantum advantage. |
| Synthetic lethality | PRISM MAPK and MYC/NAMPT dependency signals are useful, but not validated synthetic lethality. | Reviewer-reserve actionability only. Use "vulnerability prioritization", not treatment claim. |

## Key Numbers To Preserve

### Paper 1 / deconvolution

- v14 pooled MAPK x Panel-8 rho: `-0.327 [95% CI -0.376, -0.278]`, n=1,287.
- v13 TDS-16 vs Panel-8: BRAF V600E vs RAS Cohen's d `-1.6153` vs `-1.6163`; delta AUC TDS-16 - Panel-8 only `+0.007 / +0.012`, NS.
- v15A K2: TCGA centered-panel OOF AUC `0.960`; K2 `246/260` called DM2, median p_DM2 `0.978`.
- v15C PRISM: FDR<0.05 hits `11`; canonical MAPK hits `7/11`; hypergeometric p `1.08e-14`.
- Spatial v15B-v18: no positive support; final residual pooled rho approximately zero and KNN lag does not rescue.

### Paper 2 image audit

- Original RAS-like AUC: `1.000`.
- Pooled leave-one-TSS-out RAS-like AUC: `0.308`.
- TSS-only LR within RAS-like: `1.000`.
- Clinical-only LR AUC: `0.768`; multimodal LR AUC: `0.756`; CLAM-only lower.
- Raw male AUC: `0.35`.
- Honest TCGA-only ceiling: roughly `0.65-0.80`, seed- and split-dependent.
- UNI-final prediction-only audit: overall AUC `0.874`, clinical-only OOF AUC `0.735`, RAS-like AUC `0.795`, TSS-only within-RAS AUC `1.000`.
- UNI-final leave-one-TSS-out audit: pooled overall AUC `0.852`, RAS-like AUC `0.718`, BRAF-like AUC `0.855`, Female/Male AUC `0.862/0.875`.
- LUAD/BRCA multicohort audit scaffold refreshed, but both are NOT READY: missing per-slide UNI `.pt` feature files, slide manifests, cohort metadata TSVs, and OOF predictions.

### HLA package

- Paper 4 DPB1*05:01 pooled OR: `2.10 [1.70, 2.60]`, k=4, I2=55.16%.
- Paper 4 B*46:01 OR: `2.17 [1.39, 3.36]`, I2=85.09%; high heterogeneity.
- Paper 4 A*02:07 OR: `2.12 [1.73, 2.61]`.
- Paper 2 candidate DRB1*04:05: OR `3.93 [1.06, 14.53]`, p=0.0866, FDR=0.346; prospective target only.
- TCGA DM1-HLA partial rho after immune+stromal adjustment: HLA-I `0.266`, HLA-II `0.390`.

### Neoantigen / cancer vaccine

- Wave10 top clean ITSNdb_no_overlap: MHCflurry AUROC `0.668`, Structure_LR `0.653`, BigMHC_IM `0.626`.
- Wave10 large inflation examples: MIRO gap `+0.579`, GroupDRO `+0.543`, MoLE `+0.530`, ESM2-Bayesian `+0.528`.
- Quantum methods: near zero inflation but no quantum advantage; VQC AUROC `0.597`, QK_SVM `0.583` on no_overlap.
- Cancer-vaccine robustness: ITSNdb combined AUROC `0.734`, but no-overlap AUROC `0.431`; VenusVaccine top10_mean test AUROC `0.779`.
- CROSS-Neo v1 lockdown decision: HOLD; source-heldout collapse and public-overlap audit remain unresolved.

## Commit / Staging Plan

Recommended order, with separate commits:

1. **Paper 1 deconvolution rollup and v15-v18 closure**
   - Stage `project/results/p_deconv_2026_05_08/`, relevant Paper 1 assets under `project/papers_hub_2026_05_04/assets/paper1/`, and Paper 1 hub pages.
   - Review `project/manuscript_v8/` diffs carefully before staging because compiled manuscript files can include voice-protected text. Factual captions/cross-refs are allowed; protected prose rewrites are not.

2. **Paper 2 image-DM1 audit**
   - Stage `analysis_supp/audit_ras_auc100/`, audit scripts, audit figures, `REPORT_YU_REVIEW_PAPER2_2026_05_08.html`, and the Paper 2 audit dossier if present.
   - Commit message should emphasize confound audit and demotion, not a new positive image claim.

3. **HLA two-paper synthesis**
   - Stage `project/results/hla_two_paper_synthesis_2026_05_09/`, HLA reports, HLA hub pages/assets, and HLA deepdive scripts.
   - Keep raw downloaded external data out of git unless explicitly curated and small.

4. **Paper 3 ICI status**
   - Stage `project/reports/paper3_ici/`, `paper3.html`, `portfolio_paper3.html`, and Paper 3 asset updates.
   - Commit as status/dossier, not Track B launch.

5. **Neoantigen / cancer-vaccine benchmark**
   - Stage `project/results/p_neo_bayesian_2026_05_09/`, `cross_neo_*`, `p_cancer_vaccine_strategy_2026_05_09/`, and associated scripts/reports.
   - This is large conceptually but modest on disk. Consider separate commits for Wave10 mega-table, CROSS-Neo lockdown, and vaccine triage.

6. **Raw external data**
   - `project/data/external` is about 3.0 GB. Do not stage raw data by default.
   - Prefer manifests, checksums, and analysis reports unless the user explicitly wants data committed.

## Stop Rules

- Do not extend GSE250521 spatial MAPK x Panel work again unless a reviewer asks.
- Do not promote Paper 2 image-DM1 as a main-text positive claim until K2 H&E or another external WSI cohort is available.
- Do not merge HLA Paper 4 Graves/AITD allele claims with Paper 2 HT-overlap PTC expression biology.
- Do not claim quantum advantage, clinical vaccine selection, or validated synthetic lethality from the current neoantigen / vulnerability work.
- Do not edit Hook, Aim, Discussion 3.1, Limitations 3.4, Cover Letter paragraph 1, or Q9 beyond mechanical cross-reference scaffolding.

## Fast Next Move

The highest-value cleanup is a staged commit pass, not more analysis:

1. Lock Paper 1 deconvolution/v15-v18 with the caveat boundary.
2. Lock Paper 2 audit as a demotion/confound package.
3. Lock HLA two-paper synthesis separately.
4. Freeze neoantigen/CROSS-Neo as benchmark/status work unless the author explicitly prioritizes that manuscript.

## Post-Triage Runs Completed

Commands run after this addendum was created:

- Paper 2 audit smoke: `audit_ras_auc100.py` rerun successfully; key outputs match handoff (`RAS_like ∩ FVPTC = 16/16`, RAS-like AUC `1.000`, permutation p `0.0030`, folds 1/2/4 with zero RAS-like positives).
- LUAD/BRCA scaffold refresh: `audit_tcga_luad.py` and `audit_tcga_brca.py` rerun; both remain NOT READY because per-slide UNI features, slide manifests, cohort metadata TSVs, and OOF predictions are missing.
- Paper 2 audit dossier: `audit_make_dossier.py` regenerated `AUDIT_DOSSIER.html` and `paper2_image_dm1_audit_dossier.html`; deployed to `/var/www/papers/papers_hub_2026_05_04/`.
- Paper 1 methods reproducibility dossier: `build_paper1_methods_reproducibility_dossier.py` rerun and deployed to `/var/www/papers/papers_hub_2026_05_04/`.
- Hub deploy sanity: repo and live copies match for `index.html`, `paper1_deconvolution_rollup_v18.html`, `paper1_methods_reproducibility_dossier.html`, `paper2_image_dm1_audit_dossier.html`, `paper1_fig8_mechanism_dossier.html`, and `paper1_reviewer_defense_dashboard.html`.
