# Phase B — Aggressive Drug Discovery (RunPod-ready)

**Date:** 2026-05-04 (post-marathon scaffolding · author keyboard dispatch)
**Author:** Seungho Cook
**Status:** Scaffold complete — 10 sub-modules ready for RunPod dispatch.
**Cost envelope:** $80–$420 (cheapest 5 modules) → $800–$2,200 (full 10) depending on GPU choice + compound library size.
**Time envelope:** 6–18 hours per module (sequential) → 2–4 days (parallel pods).

---

## 0. Scope + binding rules

This is the **next-phase aggressive drug development pipeline** built on top of v13_drug_discovery (8 priority targets, 15 ChEMBL compounds, 8 PDB structures, 8 pocket detections, 4 docking subsets, 5 LLM drug summaries, 14 sacituzumab × thyroid clinical-trial hits).

- **Substrate:** Paper 1 v13_drug_discovery outputs (read-only).
- **Compute:** RunPod A100 / A40 / H100 (user dispatches; Claude scaffolds only).
- **Dispatch authority:** **author keyboard explicit** — `bash run_phaseB_gpu.sh --module b2 --pod $POD_ID` or per-module SSH.
- **Hard gates G1–G6** must close before each module (see §3).
- **Kill switches K1–K6** trigger immediate halt (see §3).
- **No autonomous spend** — Claude does not call RunPod GraphQL / SSH / launch.
- **Marathon discipline preserved:** voice-protected sections (Hook / Aim / Discussion / Limitations / Cover / Q9) untouched. New compute does not equal new analysis if the user decides to dispatch — it is *author-authorized* compute spend.

---

## 1. 10 Sub-modules

| Module | Name | GPU req | Time | Cost (USD) | Output |
|---|---|---|---|---:|---|
| **B1** | Structure refinement (AlphaFold3 / ESMFold rerun) | A100 80GB | 4–8 h | $30–60 | 8 refined PDB + per-residue pLDDT |
| **B2** | DiffDock pose generation (8 targets × top 50 compounds = 400 poses) | A100 40GB | 3–6 h | $20–40 | 400 docked poses + confidence |
| **B3** | GNINA rescoring + binding affinity | A40 | 2–4 h | $10–20 | Vina + CNN_score per pose |
| **B4** | KDeep / DeepPurpose binding affinity (deep learning) | A100 40GB | 4–6 h | $25–40 | pKd / pIC50 prediction per compound |
| **B5** | ChemBERTa / Mol-Former scaffold hopping (top 100 ChEMBL → 10K analogs) | A100 40GB | 6–10 h | $40–70 | 10,000 scaffold-hopped analogs + similarity |
| **B6** | PROTAC-DB matching (intracellular: CYP1B1 / PLEKHA6 / PTPRE) | A40 | 2–3 h | $10–15 | PROTAC E3-ligase × warhead matrix |
| **B7** | ADC linker / payload optimization (surface: TACSTD2 / TMPRSS4 / B3GNT3) | A100 40GB | 4–6 h | $25–40 | DAR optimization · linker stability prediction · 4 ADC head-to-head (sacituzumab govitecan / SKB264 / Dato-DXd / BNT323) |
| **B8** | Off-target profiling (SwissTargetPrediction / SEA / molecular similarity) | A40 | 2–3 h | $10–15 | Top 30 off-targets per compound + safety profile |
| **B9** | Foundation-model ADMET (Mol-Former / Chemprop / ImageMol) | A100 80GB | 6–8 h | $40–60 | Updated ADMET (Lipinski / Veber / BBB / hERG / hepatotox / clearance / volume of distribution) |
| **B10** | **FEP / MD relative binding free energy** (★ highest cost; Tier-A top 10) | H100 / A100 80GB | 12–24 h | $200–400 | ΔΔG per protein-ligand pair · pose stability · binding kinetics |

**Sequential cost (full 10 modules):** ~$410–$760 (B1–B9) + $200–400 (B10) = **$610–$1,160**
**Parallel-pod cost (3 pods × 4 modules):** **$520–$960** (faster, similar total)
**Cheapest 5 modules (B2 + B3 + B6 + B8 + B9):** **$95–$155** (dispatch in 1 day)

---

## 2. Sub-module deep dive

### B1 — Structure refinement (AlphaFold3 / ESMFold)

**Why:** v13 PDB structures heterogeneous — some from PDB direct (TACSTD2 P09758), some homology-modeled, some predicted. Refresh with AlphaFold3 (multimer + ligand-aware) for consistent quality.

**Outputs:**
- 8 refined PDB (per target)
- Per-residue pLDDT confidence
- Pocket re-detection on refined structures (compare to v13 pockets)
- Top-pocket displacement / volume change

**Critical for:** B2 (docking on better structure) + B7 (ADC binding epitope localization on TACSTD2 EC domain) + B10 (FEP requires high-quality structure).

**Tool:** AlphaFold3 (DeepMind) or AlphaFold-Multimer + ESM-IF1.

### B2 — DiffDock pose generation

**Why:** Current docking is AutoDock Vina (limited sampling, classical scoring). DiffDock = diffusion-based generative model, often 10× better pose accuracy.

**Targets × compounds:** 8 priority × top 50 (per-target) = 400 poses. Plus extension: top 30 ChemBERTa scaffold-hopped (B5) → 8 × 30 = 240 additional poses.

**Output:** 400 + 240 poses with confidence score, RMSD to reference (when available), binding mode visualization.

### B3 — GNINA rescoring

**Why:** DiffDock pose, AutoDock Vina pose, and KDeep prediction all need consistent scoring layer. GNINA = CNN-based, validated on PDBbind.

### B4 — KDeep / DeepPurpose

**Why:** Direct binding-affinity prediction from sequence + ligand SMILES. Independent of docking pose — orthogonal validation.

**Models:** KDeep (Bioinformatics 2018), DeepPurpose (CSBJ 2020), MolBERT-AffinityNet.

### B5 — Scaffold hopping

**Why:** Current ChEMBL 15 compounds include some redundant scaffolds. Use ChemBERTa embedding to find structurally diverse analogs of the top binders.

**Output:** 10K analogs with chemical similarity score · bioavailability · ADMET prediction.

### B6 — PROTAC matching (intracellular targets)

**Why:** CYP1B1 / PLEKHA6 / PTPRE are intracellular — ADC inaccessible. PROTAC = bivalent ligand (warhead + E3 ligase recruiter via linker). PROTAC-DB has 5,000+ literature PROTACs to match against pocket.

**Output:** Per-target E3 × warhead matrix · feasibility ranking.

### B7 — ADC linker / payload (★ Sacituzumab Phase B deep dive)

**Why:** Paper 1 R5 의 핵심. Sacituzumab govitecan vs SKB264 / Dato-DXd / BNT323 head-to-head computational comparison on TACSTD2 (TROP2):

- **Sacituzumab govitecan**: SN-38 (Topo-I inhibitor), DAR 7.6, hydrolysable CL2A linker
- **SKB264 / sacituzumab tirumotecan (MK-2870)**: T-030 (Topo-I), DAR 7.4, tumor-selective cleavable linker
- **Datopotamab deruxtecan (Dato-DXd)**: DXd (exatecan analog Topo-I), DAR 4, cleavable tetrapeptide
- **BNT323 / DB-1303**: TROP2 × exatecan, DAR ~8, stable peptide linker

**Computational comparison:**
- TROP2 epitope binding affinity (per ADC monoclonal Fab)
- Internalization t½ prediction (rate-limiting for SN-38 release)
- Payload release kinetics (linker pH / cathepsin sensitivity)
- TROP2 PTC expression heterogeneity (DM1 high vs low → DAR-adjusted dosing)

**Why aggressive:** 4 thyroid Phase 2 trials currently active (NCT06235216 / NCT06923826 / NCT07521670 STRAP / NCT07068542). DM1-high prospective biomarker layer = differentiation point for Paper 1.

### B8 — Off-target profiling

**Why:** Sacituzumab govitecan TNBC has neutropenia + nausea AE. Off-target prediction for thyroid context (different tissue distribution).

**Tools:** SwissTargetPrediction · SEA (Similarity Ensemble Approach) · MolNet × ChEMBL similarity.

### B9 — Foundation-model ADMET

**Why:** v13 ADMET = Lipinski + Veber + RDKit-based heuristics. Foundation models (Mol-Former, Chemprop) give more accurate predictions on:
- BBB permeability (CNS off-target)
- hERG (cardiotoxicity)
- Hepatotoxicity (CYP-mediated)
- Aqueous solubility
- Plasma protein binding
- Renal clearance
- Volume of distribution
- Half-life prediction

### B10 — FEP / MD relative binding free energy (★ Tier-A top 10)

**Why:** Final validation step. ΔΔG (relative binding free energy) is the gold-standard for ranking compounds.

**Compute:** OpenFF + GROMACS or OpenFE / Schrödinger FEP+. ~24h per protein-ligand pair on H100.

**Output:** ΔΔG per Tier-A top 10 (e.g., propofol-GABRB2, simvastatin-LDLR, etomidate-GABRB2) ± uncertainty.

**Cost-control:** Run only on **Tier-A top 10** (not all 59) — that's 10 H100-hours = ~$200–400.

---

## 3. Hard gates G1–G6 + Kill switches K1–K6

### G1–G6 (must close before each module)

| Gate | Condition | Action if not met |
|---|---|---|
| **G1** | RunPod credits ≥ $200 | Skip B10 (FEP); run B1–B9 only |
| **G2** | All 8 PDB structures available + verified (sha256) | Skip B2/B3/B7 for missing target |
| **G3** | Top 50 ChEMBL compounds per target available (per_target_top10 + ranked_candidates) | B2/B3 narrow scope |
| **G4** | TIER-A 'Go Now' top 30 verified | Skip B10 (no Tier-A → no FEP candidates) |
| **G5** | Internet access for ChemBERTa / Mol-Former / DiffDock model download | Skip B5/B9 (offline-impossible modules) |
| **G6** | User explicit `bash run_phaseB_gpu.sh --module bN --pod $ID` | Pipeline blocked. "고고" / "다 해줘" / "faster" do NOT count. |

### K1–K6 (immediate halt)

| Kill | Trigger | Reason |
|---|---|---|
| **K1** | "ADC therapy approved for thyroid" claim | No thyroid ICI approval; sacituzumab thyroid trials unblinded outcome unknown |
| **K2** | "Wet-lab IC50 derived from this pipeline" | This is in-silico only — no IC50 spectroscopy |
| **K3** | "PTC patient stratification by computational model" | requires wet-lab + clinical trial validation |
| **K4** | DiffDock confidence < 0.5 for all 50 poses on a target | Pose unreliable → re-run with refined PDB (B1) |
| **K5** | FEP ΔΔG > 5 kcal/mol uncertainty | unreliable; re-run with longer simulation or alternative force field |
| **K6** | Compute spend > $500 without user re-authorization | Halt; new authorization required |

---

## 4. RunPod dispatch (author authority)

```bash
# 0. SSH into RunPod (user-side)
ssh root@$POD_IP -p $POD_PORT

# 1. Setup environment (~ 15 min)
cd /workspace
git clone https://github.com/seunghocook/thyca-paper-2026 .
cd azure_gpu_phaseB_drug
bash scripts/setup.sh

# 2. Per-module dispatch (run in parallel via tmux / different pods)
bash scripts/b1_structure_refinement.sh    # 4-8h, $30-60
bash scripts/b2_diffdock_pose.sh            # 3-6h, $20-40
bash scripts/b3_gnina_rescore.sh            # 2-4h, $10-20
bash scripts/b4_kdeep_affinity.sh           # 4-6h, $25-40
bash scripts/b5_chemberta_scaffold.sh       # 6-10h, $40-70
bash scripts/b6_protac_db_match.sh          # 2-3h, $10-15
bash scripts/b7_adc_optimization.sh         # 4-6h, $25-40
bash scripts/b8_off_target.sh               # 2-3h, $10-15
bash scripts/b9_admet_foundation.sh         # 6-8h, $40-60
bash scripts/b10_fep_md.sh                  # 12-24h, $200-400 (★ Tier-A top 10 only)

# 3. SCP results back to local (after each module)
scp -P $POD_PORT root@$POD_IP:/workspace/azure_gpu_phaseB_drug/results/* \
    /home/seungho/personal/THCA_data_analysis/project/results/v14_drug_phaseB/

# 4. Update Paper 1 hub (Claude scaffolds; user dispatches)
# After all modules complete, regenerate Paper 1 supplement section with new TSVs
```

---

## 5. Output structure (post-Phase B)

```
project/results/v14_drug_phaseB/
├── b1_structures/
│   ├── refined/{TACSTD2,TMPRSS4,PLEKHA6,CYP1B1,LDLR,GABRB2,B3GNT3,PTPRE}.pdb
│   ├── plddt_per_residue.tsv
│   └── pocket_redetect.tsv
├── b2_diffdock/
│   ├── poses/{target}_{compound}_pose_{N}.pdb
│   ├── pose_confidence.tsv
│   └── rmsd_to_reference.tsv
├── b3_gnina/
│   └── rescore.tsv
├── b4_kdeep/
│   └── pkd_predictions.tsv
├── b5_chemberta/
│   ├── scaffold_hopped.tsv (10K analogs)
│   └── chemberta_embedding.npz
├── b6_protac/
│   └── e3_warhead_matrix.tsv
├── b7_adc_optimization/
│   ├── trop2_4adc_compare.tsv (sacituzumab gov / SKB264 / Dato-DXd / BNT323)
│   ├── linker_stability.tsv
│   └── dar_optimization.tsv
├── b8_off_target/
│   └── per_compound_top30.tsv
├── b9_admet_foundation/
│   └── admet_v2.tsv (adds clearance / Vd / half-life / PPB)
└── b10_fep_md/
    ├── ddg_tier_a_top10.tsv
    └── md_trajectory_summary.tsv
```

---

## 6. Cross-paper bridge

Phase B outputs feed:
- **Paper 1 R5 supplement** — sacituzumab govitecan vs 3 alternative ADCs head-to-head + PTC TROP2 heterogeneity
- **Paper 1 future plan Tier 2** — wet-lab IC50 (sacituzumab in PTC cell lines) is downstream after B7 + B10
- **Paper 1 future plan Tier 3** — single-cell scRNA validation of TROP2 expression heterogeneity
- **Paper 7 (agentic_research)** — Phase B is a Pattern 4 (pre-registered failed hypothesis) candidate scaffold (if FEP fails for top compounds, that's a publishable negative result)

---

## 7. Marathon attestation (Phase B only — not yet executed)

- ✅ Voice-protected sections untouched (this scaffold doesn't write Hook / Aim / Disc / Limitations / Cover / Q9)
- ✅ Paper 3/4 design bundles untouched
- ✅ No new analysis launched by Claude (RunPod dispatch = author keyboard)
- ✅ No new data download by Claude
- ✅ chmod 444 freezes preserved
- ✅ K1 (no "ICI predictor" claim) preserved across modules

**At dispatch time:**
- ✅ G6 explicit author command required
- ✅ K1–K6 kill switches monitored
- ✅ Cost cap $500 default

---

## 8. Pre-flight checklist

- [ ] RunPod credits checked (≥ $200 for cheap 5-module run; ≥ $500 for full 10)
- [ ] G1–G5 hard gates verified (PDB / ChEMBL / Tier-A / internet ready)
- [ ] G6 explicit author authorization
- [ ] tmux on pod for parallel module run
- [ ] SCP credentials verified
- [ ] Local `project/results/v14_drug_phaseB/` directory created
- [ ] Marathon discipline reaffirmed — Phase B is **author-authorized one-shot compute spend**, not a return to continuous compute mode.
