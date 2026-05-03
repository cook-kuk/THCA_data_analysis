# Paper 3 ICI — Track A Freeze Notice (2026-05-04)

**Status:** Track A FROZEN. Track B BLOCKED. Marathon mode (Paper 1 + Paper 2 writing, 5/4–6/13) preserved.
**Author:** Seungho Cook
**Freeze date:** 2026-05-04
**Single source of truth:** this notice + the design bundle at `project/reports/paper3_ici/PAPER3_ICI_TRACK_A_BUNDLE.md`.

---

## 1. Paper 3 freeze state

| Field | Value |
|---|---|
| Paper | Paper 3 |
| Working title | HLA loss, neoantigen architecture, and immune ecotypes define ICI vulnerability in molecularly dark thyroid cancer |
| One-line scope | Dark-matter thyroid cancer ICI vulnerability / immunogenomics |
| Allowed claim language | "ICI vulnerability", "ICI-readiness", "immunogenomic prioritization" |
| Forbidden claim language | "ICI response predictor" — without thyroid ICI-treated raw RNA-seq |
| Track A status | COMPLETED + FROZEN (chmod 444 on design bundle) |
| Track B status | BLOCKED |
| Track B unlock conditions | (1) Paper 1 bioRxiv / submission complete + (2) Paper 2 Task A/B/C closure + (3) explicit user command "Paper 3 Track B 시작" |
| Memory anchor | `~/.claude/projects/-home-seungho-personal-THCA-data-analysis/memory/v19_paper3_ici_track_a.md` |
| Canonical numbering | `paper_numbering_2026_05_04.md` (Paper 1=DM1 molecular / Paper 2=H&E→DM1 pathology+TCGA validation / Paper 3=ICI dark matter / Paper 4 backlog=Korean GD HLA) |

---

## 2. Tar.gz inventory

**Archive:** `project/reports/paper3_ici/paper3_ici_track_a.tar.gz`
**Compressed size:** 51,023 bytes
**Sha256:** `290349f6ebed7de64c448afbebf633baa05f1845e0685f7e7a2df9fb3d8936aa`
**Member count:** 12

| Path inside archive | Bytes | Mtime | Mode | Role |
|---|---:|---|---|---|
| `PAPER3_ICI_TRACK_A_BUNDLE.md` | 71,684 | 2026-05-04 00:16 | 444 | Single-file consolidated design (FROZEN marker prepended) |
| `paper3_ici_dataset_registry.md` | 9,947 | 2026-05-04 00:04 | 444 | Deliverable 1 — bulk / scRNA / spatial / pan-cancer ICI dataset registry |
| `paper3_ici_signature_registry.md` | 9,887 | 2026-05-04 00:05 | 444 | Deliverable 2 — signature registry + integrated score composition |
| `paper3_ici_scRNA_reference_atlas_summary.md` | 7,958 | 2026-05-04 00:06 | 444 | Deliverable 3 — scRNA atlas plan |
| `paper3_ici_HLA_neoantigen_feasibility.md` | 8,717 | 2026-05-04 00:07 | 444 | Deliverable 4 — HLA typing + LOHHLA + NetMHCpan feasibility |
| `paper3_ici_DIAL_audit_plan.md` | 8,383 | 2026-05-04 00:08 | 444 | Deliverable 5 — DIAL audit plan against pan-cancer ICI cohorts |
| `paper3_ici_figure_plan.md` | 8,570 | 2026-05-04 00:09 | 444 | Deliverable 6 — 6 main + 14 supp figure plan |
| `paper3_ici_go_no_go_verdict.md` | 5,794 | 2026-05-04 00:10 | 444 | Deliverable 7 — Conditional GO with G1–G6 gates and K1–K6 kill switches |
| `paper3_ici_12week_execution_plan.md` | 11,602 | 2026-05-04 00:11 | 444 | Deliverable 8 — 12-week Track B execution schedule |
| `todo/paper3_ici_track_b_kickoff_checklist.md` | 4,895 | 2026-05-04 00:16 | 644 | Track B kickoff verification checklist (editable) |
| `todo/paper3_ici_controlled_access_checklist.md` | 4,708 | 2026-05-04 00:17 | 644 | Controlled-access application checklist (editable) |
| `todo/paper3_ici_dataset_accession_verification_checklist.md` | 6,234 | 2026-05-04 00:18 | 644 | Dataset accession verification checklist (editable) |

**Per-file integrity (host filesystem):**
- `PAPER3_ICI_TRACK_A_BUNDLE.md` sha256 `95e2a59b5a9e0f28e2097bcf13abc94aa9b999ee31515fc4aec2f62bbf5dbdf1`

**Mode convention:** `444` = read-only frozen design artifact. `644` = editable todo (intended to receive checkmarks during Track B kickoff and beyond).

---

## 3. Design bundle index (one-line each)

1. **Dataset Registry** — bulk ≤9 cohorts, scRNA ≤6, spatial ≤4, pan-cancer ICI ≤7; all `to_verify` until Track B Wk 1; cross-paper boundary discipline annotated per cohort.
2. **Signature Registry** — IFNγ/TIS/Cytolytic/MHC-I/MHC-II/TLS/CXCL13/Myeloid/Treg/TIDE-like/IMPRES + tumor-intrinsic (BRS/ERK/TDS/lineage) + sc-derived placeholders; integrated score composition with weights deferred to Module E.
3. **scRNA Atlas Plan** — scVI/scANVI/Harmony integration, Level-1/Level-2 taxonomy, sc-derived signature lock for bulk projection; PDTC/ATC <5K caveat.
4. **HLA & Neoantigen Feasibility** — OptiType/Polysolver/arcasHLA/xHLA + LOHHLA + NetMHCpan/NetMHCIIpan; cookHLA reserved for Paper 4; degraded-scope contingencies if dbGaP not granted.
5. **DIAL Audit Plan** — sign-consistency × effect-size × tissue-transfer; PASS/FLIP/COLLAPSE/AMBIGUOUS verdict per signature; thyroid-adjusted readiness score restricted to PASS signatures.
6. **Figure Plan** — F1 cohort/dark matter / F2 ecotype NMF / F3 sc atlas / F4 HLA+neoantigen / F5 DIAL / F6 integrated score; 14 supp; cross-paper boundary lines in every relevant caption.
7. **Go/No-Go Verdict** — Conditional GO. Hard gates G1–G6. Kill switches K1–K6 (esp. K1 = no "ICI response predictor" claim).
8. **12-Week Execution Plan** — Wk 1 verify+ETL → Wk 10 integrated score → Wk 11–12 manuscript draft + cross-paper boundary sweep.

---

## 4. Checklist files (editable, kept outside the freeze)

| Path | Purpose | Editable | First-touch trigger |
|---|---|---|---|
| `paper3_ici/todo/paper3_ici_track_b_kickoff_checklist.md` | G1–G6 gate confirmation; design-bundle re-read; environment readiness | yes | At Track B kickoff only |
| `paper3_ici/todo/paper3_ici_controlled_access_checklist.md` | dbGaP / EGA / institutional applications; marathon-displacement guard (≤30 min/session, no Paper 1/2 writing-block displacement) | yes | Background-only during marathon |
| `paper3_ici/todo/paper3_ici_dataset_accession_verification_checklist.md` | Portal-lookup verification of all `to_verify` accessions; default = defer until Track B Wk 1 | yes | At Track B kickoff (or trivially-fast inline lookups) |

---

## 5. SCP retrieval (Azure VM 40.82.129.113)

```bash
# tarball (51 KB; design bundle + 8 design docs + 3 todo checklists)
scp seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/project/reports/paper3_ici/paper3_ici_track_a.tar.gz ./

# integrity check after transfer
sha256sum -c <(echo "290349f6ebed7de64c448afbebf633baa05f1845e0685f7e7a2df9fb3d8936aa  paper3_ici_track_a.tar.gz")

# extract
tar xzf paper3_ici_track_a.tar.gz
```

---

## 6. Boundary statement (binding)

- Paper 3 design bundle is **read-only** until Track B unlock. Any edit to a chmod-444 file requires explicit unlock + reason recorded here.
- Paper 1 (DM1 molecular) and Paper 2 (H&E → DM1 / pathology projection / TCGA validation) marathon work is the **active** focus.
- Paper 4 (Korean GD HLA backlog) remains gated per `v19_paper4_GD_backlog.md`; no work during Paper 3 freeze either.
- "고고" / "다 해줘" / "faster" / "추가 분석 더 해줘" do **not** override Track B block per `v17_sprint_vs_marathon_violation.md`. Unlock requires the literal phrase "Paper 3 Track B 시작".

---

Paper 3 Track A frozen. Return to Paper 1/2 marathon.
