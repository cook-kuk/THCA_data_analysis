# Paper 2 BRAF Nature Sprint — Synthesis v2 / Codex takeover

Build: 2026-05-09 KST. Scope: factual synthesis and operational handoff only. Voice-protected manuscript sections remain author-only.

## Immediate Decision

Use **Option C** now.

- `thca-img-dm1-a6000` is already `EXITED`.
- Idle `thca-neo-bayesian-aux` was stopped by Codex (`r8v801csxxscw0` -> `EXITED`) to stop spend.
- `thca-spark-dm-a6000-v4` finished the active GPU jobs; Codex recovered the strict ESMFold outputs and stopped the pod to `EXITED` to cut spend.
- Do not restart H15v2/foundation-image rescue tonight. No usable H15 result was present in the local `h15_foundation_image/` directory at takeover, and the current image conclusion is already sufficient: generic CLAM does not rescue BRAF/cPTC; Paper 2 should be reframed as stratum-complementary, not as an immediate Nature-grade image model.

## One-Line v2 Conclusion

The new sprint converts the story from "image-DM1 rescue" into a stronger **Paper 1 biology/clinical axis**: in BRAF-driven cPTC, DM1 is an immune-rich, thyroid-lipid/metabolic compartment that is mostly progression-protective, while **DM2 is the BRAF-context aggressive/RAI-refractory arm**; the exception/escape is **DM1 x TERT**, which now replicates externally in Landa and pools with TCGA at HR about 19. Paper 2 survives as a stratum-complementary pathology pilot: image works in FVPTC/RAS-like morphology, but not in well-differentiated BRAF/cPTC immune biology.

## What Changed Since v1

| Claim | v1 status | v2 update | Disposition |
|---|---|---|---|
| DM1 x TERT is the aggressive niche | Replication-gated | H8/R1: Landa OS HR 12.78; TCGA+Landa REML HR 19.2 [7.16, 51.3], p=4.2e-9, I2=26% | **Hard gate passes** |
| DM2 is the BRAF bad-actor | TCGA H6 only | R1: BRAF-like PFI HR 4.67 [1.41, 15.5], p=0.012; H12: DM2 directionally more RAI-refractory and PFI HR 5.68 | **Main clinical frame** |
| Lee 2024 can replicate survival | Unknown | R4/H30: no per-sample survival, driver, or TERT; expression-side only | **Do not claim survival replication** |
| HT-13 is unique | Weak | R3: 8/13 genes unique; scalar score co-varies with TIS/TLS; Korean hashi-AUC 0.910 > TIS 0.857 | **Thyroid-Hashimoto-tuned, not novel** |
| Immune axis is purity artifact | Open | R6: HT-13 residual d remains +0.88 after leukocyte residualization; high-purity tertile d stays positive | **Reviewer defense locked** |
| Bulk HT source | Open | H10/R7: HT-13 peaks in myeloid/B/T; FA-12 peaks in malignant cells; MAPK-9 peaks in BRAF malignant | **Compartment model locked** |
| Methylation mechanism | Correlative | R9/H13: 8-gene HM450 captures dominant methylation layer; immune RNA collapses after residualizing 8-gene beta; mediation is stratum-dependent | **Keep Q9 boundary** |
| ICI translation | TCGA-internal | R5/R10: HT-13 x HLA-I pooled OR 1.67 for response, but PD-L1 IHC beats HT-13 in IMvigor | **Paper 3 reserve** |
| H&E rescue | Hopeful | H7 negative; R8 FVPTC saliency only trend; no usable local H15 output at takeover | **Future work, not load-bearing** |
| UNI image audit | Pending center-holdout | GPU LOTO complete: overall AUC 0.852, BRAF-like 0.855, RAS-like 0.718 | **Paper 2 pilot survives, not main claim** |
| K2 full-transcriptome pilot | Running | Completed 9/9 local FASTQs; all samples DM2; HT tumor-vs-normal AUC 0.40 | **Negative-control reserve** |

## Evidence Stack To Save

| Layer | Lead result | Source |
|---|---|---|
| Clinical BRAF context | DM2 vs not-DM PFI HR 4.67 [1.41, 15.5], p=0.012 | `r1_meta_survival/R1_REPORT.md` |
| Escape switch | DM1 x TERT TCGA+Landa OS REML HR 19.2 [7.16, 51.3], p=4.2e-9 | `r1_meta_survival/R1_REPORT.md` |
| External TERT replication | Landa DM1+/TERT+ OS HR 12.78, 14/14 deaths | `h8_tert_replication/H8_REPORT.md` |
| RAI/aggressive arm | BRAF-cPTC DM2 PFI HR 5.68; DM2 refractory rate 16% vs DM1 6.9% | `h12_rai_response/H12_REPORT.md` |
| Compartment biology | Myeloid/B/T carry HT-13; malignant cells carry FA-12; BRAF malignant carries MAPK-9 | `r7_singlecell/R7v2_REPORT.md` |
| Deconvolution | M2 macrophage d=+1.83, DC +1.74, Treg +1.43, TLS +0.96 in DM1 | `h10_celltype_deconv/H10_REPORT.md` |
| Purity defense | HT-13 d 1.42 -> 0.88 after leukocyte residualization; high-purity signal remains | `r6_purity_audit/R6_REPORT.md` |
| Spatial support | MAPK-active spatial subset HLA-DRA HIGH-vs-LOW delta +0.88, p=0.0079 | `r2_spatial_driver_inferred/R2_REPORT.md` |
| Korean expression support | Lee 2024 n=632 has full HT/FA/MAPK coverage; HT x FA r=-0.006; HT hashi-AUC 0.910 | `h30_lee2024_full/H30_REPORT.md`, `r3_panel_disentanglement/R3_REPORT.md` |
| Image pilot | FVPTC image AUC works; BRAF/cPTC CLAM rescue fails | `h1_symmetric_stratum/H1_REPORT.md`, `h7_perm_null/H7_REPORT.md`, `r8_fvptc_saliency/R8_REPORT.md` |

## What To Kill Or Demote

- Kill: "DM1 is aggressive" as a generic claim. Correct frame: DM1 is usually protective in BRAF-cPTC; **DM1 x TERT** is the aggressive escape.
- Kill: "multimodal image+RNA rescue" as a main claim. Current image channel adds no reliable BRAF/cPTC immune-axis readout.
- Kill: Lee 2024 survival/TERT/driver replication. Public data do not contain those fields.
- Demote: HT-13 as a unique immune signature. Use "thyroid-Hashimoto-tuned inflamed-tumor signature family member."
- Demote: K2 full-transcriptome re-quant. The local run completed, but it is n=9 with 5 tumor / 4 matched-normal samples, all prior calls DM2. It is a completed technical pilot, not Korean DM1-vs-DM2 replication.

## Figure Map v2

| Figure | Content | Status |
|---|---|---|
| Fig 1 | BRAF-cPTC DM1/DM2 clinical split + DM1 x TERT escape schematic | ready from H6/H8/R1 |
| Fig 2 | HT immune x FA thyroid-lipid orthogonal biology, with NMF/data-driven support | ready from H2/H5 |
| Fig 3 | Compartment source model: bulk deconv + scRNA + spatial HLA-DRA | ready from H10/R7/R2/H11 |
| Fig 4 | Cross-cohort expression/protein validation including Lee 2024 and Mun | ready from H4/H30 |
| Fig 5 | Clinical translation: DM2 PFI/RAI + DM1 x TERT pooled forest | ready from R1/H8/H12 |
| Fig 6 | H&E stratum-complementarity: FVPTC image yes, BRAF-cPTC image no | ready from H1/H7/R8 |
| Supp | Purity, HT-13 vs TIS/TLS, methylation boundary, ICI reserve | ready from R3/R6/R9/R5/R10 |

## Operational State

### RunPod

| Pod | State | Action |
|---|---|---|
| `thca-img-dm1-a6000` (`ytnsyachant0a4`) | `EXITED` | leave stopped |
| `thca-neo-bayesian-aux` (`r8v801csxxscw0`) | stopped to `EXITED` | leave stopped unless explicitly needed |
| `thca-spark-dm-a6000-v4` (`uvp9i2r9s6l85y`) | stopped to `EXITED` | GPU jobs finished; strict ESMFold outputs recovered |

### Local K2 quant

Working dir: `project/results/p2_braf_nature_sprint_2026_05_09/h_k2_full_quant/`.

- Runner PID at takeover: `3965359 bash run_kallisto_quant.sh`; completed at 03:21 KST.
- Script saw `N samples with FASTQ: 9`; aggregate exited 0.
- Metadata has 262 runs, but local FASTQ availability is: 253 no FASTQ, 5 R1+R2, 4 R1-only.
- Output summary: 61,228 gene symbols, HT-13/FA-12/MAPK-9 coverage 100%, 5 tumors / 4 matched normals, all prior 8-gene calls DM2.
- Tumor-vs-normal panel tests are not deployable: HT-13 AUC 0.40 (MW p=0.730); FA-12 AUC 0.75 and MAPK-9 AUC 0.75 (both MW p=0.286).
- Treat outputs as **negative-control/reviewer-reserve only** unless the missing 253 runs are fetched and paired FASTQs are complete.

Monitor:

```bash
cat project/results/p2_braf_nature_sprint_2026_05_09/h_k2_full_quant/K2_PILOT_VERDICT.md
cat project/results/p2_braf_nature_sprint_2026_05_09/h_k2_full_quant/k2_panel_aucs.tsv
```

Aggregate after quant finishes:

```bash
python3 project/results/p2_braf_nature_sprint_2026_05_09/h_k2_full_quant/summarize_k2_pilot.py
```

## Next Action

For manuscript work, do not write voice-protected sections. The highest-ROI safe next step is factual scaffolding only:

1. Update figure/caption scaffolds with the v2 map above.
2. Build a Paper 1/2 advisor dossier page from the v2 synthesis and existing figures.
3. Leave K2 and H15v2 as reviewer-reserve/future-work unless full paired FASTQs or K2 H&E become available.
