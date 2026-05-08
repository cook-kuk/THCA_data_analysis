# Paper 1 one-page web deployment — result

**Date:** 2026-05-07
**URL:** http://40.82.129.113:8012/manuscript_v8/p1_onepage_audit.html
**Server:** THCA-Secure/1.0 (PID 1690659), docroot = `/home/seungho/personal/THCA_data_analysis/project/`, port 8012
**HTTP status:** 200 OK ✅

---

## 1. Files created

### Reports / registries (`project/reports/`)
- `2026_05_06_paper1_inventory_files.txt` — 451-line result-file inventory
- `2026_05_06_paper1_dataset_suitability_registry.tsv` — 21 datasets × 19 columns
- `2026_05_06_paper1_gene_set_hierarchy.tsv` — 17 entities
- `2026_05_06_paper1_figure_logic_registry.tsv` — 25 figures
- `2026_05_06_paper1_claim_boundary_and_reviewer_risk.tsv` — 21 topics
- `2026_05_06_paper1_comprehensive_audit_master.md` — master audit
- `2026_05_06_paper1_onepage_web_deployment_result.md` — this file

### Web page + assets (`project/manuscript_v8/`)
- `p1_onepage_audit.html` (40,797 bytes; ~16 sections, sticky TOC)
- `assets/p1_onepage_audit/` (19 PNG figures + 3 Python scripts + panel_combos.tsv)

---

## 2. Figures created (19)

| # | File | Type |
|---|------|------|
| 01 | `fig01_clinical_workflow_current_vs_future.png` | schematic |
| 02 | `fig02_dataset_suitability_matrix.png` | schematic |
| 03 | `fig03_gene_set_hierarchy.png` | schematic |
| 04 | `fig04_main_evidence_chain.png` | schematic |
| 05 | `fig05_external_validation_summary.png` | data bar |
| 06 | `fig06_aggressive_cohort_context.png` | data bar |
| 07 | `fig07_main_vs_supplement_map.png` | registry view |
| 08 | `fig08_claim_boundary_board.png` | schematic |
| 09 | `fig09_reviewer_risk_heatmap.png` | schematic |
| 10 | `fig10_keep_move_drop_figure_registry.png` | registry view |
| 11 | `fig11_narrative_storyline.png` | storyline arc |
| 12 | `fig12_rai_biology_pathway.png` | thyroid biology pathway |
| 13 | `fig13_mechanism_cascade.png` | mechanism cascade |
| 14 | `fig14_clinical_timeline.png` | clinical timeline |
| 15 | `fig15_korean_cohort_applicability.png` | Korean board |
| 16 | `fig16_two_layer_architecture.png` | discovery vs deployment |
| 17 | `fig17_future_validation_board.png` | missing pieces board |
| 18 | `fig18_main_figures_meaning.png` | 5 main figures meaning |
| ★ | `fig_panel_combos.png` | NEW gene panel combination scan |

Plus reused existing figures from `project/results/p_external_expression_validation/` and `project/results/dm1_subcluster_diagnosis_2026_05_07/`.

---

## 3. Tables included (in HTML page)

| Table | Section |
|-------|---------|
| Main 11 datasets summary | §4.1 |
| Supp 6 datasets | §4.2 |
| Drop 4 datasets | §4.3 |
| Allowed vs Forbidden wording | §11 |
| Top 10 reviewer risks + answer direction | §12 |
| Missing pieces + future validation | §13 |
| Title options ranking | §15 |
| **Panel combination scan results (NEW)** | §10 |

---

## 4. Key interpretive conclusions

- **Paper 1 = compact RAI-lineage transcriptomic readout paper** (NOT a "novel 8-gene discovery paper")
- **2-layer architecture**: discovery axis (TIERA67/pan-genome ARI 0.90/0.92) + deployment readout (8-gene AUC 0.962, ARI 0.49 by design)
- **Storyline spine (6 beats)**: HOOK → AXIS → READOUT → REPLICATION → MECHANISM → REFRAME
- **Clinical reframe**: candidate triage scaffold for early flag of post-surgery RAI-failure biology (HYPOTHESIS only)
- **External validation core**: GPL570 4-cohort ρ −0.84 to −0.94 with zero-overlap NONOVERLAP panel; Korean K2 + Lee 2024 East-Asian generalizability
- **Mechanism cascade**: lineage TF backbone collapse → DNMT/STAT3/AP-1 activation → promoter hypermethylation → RAI machinery silencing
- **Survival audit**: OS event 3.4% sparse → PFI primary; age confound dominant; stage ⊥ DM axis
- **Sub-A/B versioning gap**: 17/110 overlap; re-derivation queued
- **Spatial cohort GSE250521**: marginal stage trend (ρ ≈ 0); supportive supp only; "PT" = Para-Tumor (not Papillary Tumor)
- **★ Panel combination scan finding**: 6-gene "Lineage TF + Effector minimal" (FOXE1, NKX2-1, PAX8, TG, TPO, DIO1) outperforms canonical RAI_8 (AUC 0.893 vs 0.861); recommend supp inclusion

---

## 5. What was explicitly excluded

- ✅ No voice-protected manuscript prose written (Hook / Aim / Discussion §3.1 / Limitations / Cover / Q9)
- ✅ No new GEO datasets / search
- ✅ No raw FASTQ / CEL / WSI / methylation download
- ✅ No GPU / RunPod
- ✅ No H&E / pathology DM1 retry
- ✅ No TROP2 main claim / title restoration
- ✅ No Paper 2 / Paper 3 / Paper 4 territory touch
- ✅ No clinical treatment recommendation claim
- ✅ No validated RAI response predictor claim
- ✅ No prospective clinical utility claim
- ✅ No commit performed

---

## 6. URL

**Primary:** http://40.82.129.113:8012/manuscript_v8/p1_onepage_audit.html

**Assets:** http://40.82.129.113:8012/manuscript_v8/assets/p1_onepage_audit/

**Reused figures:** http://40.82.129.113:8012/results/p_external_expression_validation/ + http://40.82.129.113:8012/results/dm1_subcluster_diagnosis_2026_05_07/

---

## 7. HTTP check

```
$ curl -sI http://40.82.129.113:8012/manuscript_v8/p1_onepage_audit.html
HTTP/1.0 200 OK
Server: THCA-Secure/1.0 Python/3.12.3
```

---

## 8. Commit proposal (NOT executed)

### P1-W1 audit registries
```
project/reports/2026_05_06_paper1_dataset_suitability_registry.tsv
project/reports/2026_05_06_paper1_gene_set_hierarchy.tsv
project/reports/2026_05_06_paper1_figure_logic_registry.tsv
project/reports/2026_05_06_paper1_claim_boundary_and_reviewer_risk.tsv
project/reports/2026_05_06_paper1_comprehensive_audit_master.md
project/reports/2026_05_06_paper1_inventory_files.txt
```
**Suggested message:** `audit: paper 1 comprehensive registries (datasets/genes/figures/claims) + master MD`

### P1-W2 web page
```
project/manuscript_v8/p1_onepage_audit.html
project/manuscript_v8/assets/p1_onepage_audit/
```
**Suggested message:** `web: paper 1 one-page comprehensive audit page (19 figures + panel combo scan + storyline arc)`

### P1-W3 deployment report
```
project/reports/2026_05_06_paper1_onepage_web_deployment_result.md
```
**Suggested message:** `docs: paper 1 one-page web deployment result (URL + HTTP check)`

### Excluded from commit
- voice-protected manuscript prose (none written)
- Paper 2 / 3 / 4 files (none touched)
- raw data (none)
- H&E / WSI / pathology branches (none)
- speculative new analysis outputs without source backing

---

**Final line:** Paper 1 one-page comprehensive audit page deployed. No voice-protected prose written.
