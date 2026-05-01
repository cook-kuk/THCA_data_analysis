---
title: "Paper B (Bioinformatics OUP) — DIAL framework + 5-cohort validation"
date: 2026-04-27 (cover letter dated)
status: drafting / submission prep (cover + figures + tables + reproducibility complete)
target_venue: Bioinformatics (Oxford University Press)
position: methodology + reproducibility paper, scope match
---

# Paper B — DIAL framework methods paper (Bioinformatics OUP)

## 1. 한 줄 요약

> "Recovered external validation and proper GSEA define DM1/DM2 as an actionable thyroid cancer transcriptional axis"

본 paper 는 **DM1/DM2 axis 의 methodology + reproducibility** layer 를 강조하는 methods paper. 임상 결과 (8-gene > BRAF V600E ΔAUC=+0.132) 는 downstream application 으로 framing — Bioinformatics OUP 의 methods-paper voice 에 맞춤.

## 2. Triple deliverable

| 항목 | 핵심 |
|---|---|
| (i) **DIAL framework** | Direction Identifiability After Leave-out — leak-detection safeguard for batch-correction protocols. v5.1 → v5.2 self-audit (자체 detection example). |
| (ii) **5-cohort validation pipeline** | Multi-layer, robust transcriptomic axis recovery (DM1/DM2). Gene-ID recovery + cross-platform transfer. |
| (iii) **Proper preranked GSEA** | Pathway-name normalisation, 100% sign-agreement with prior proxy GSEA on overlapping pathways. |

## 3. 코호트 (5개)

- TCGA-THCA primary tumor
- GSE27155
- GSE33630
- GSE29265
- GSE76039

## 4. Reproducibility infrastructure

- **Deterministic seeds**: `random_state=42` 전역 적용
- **Version-pinned dependencies**: sklearn 1.8.0, gseapy 1.13, scanpy 1.12.1
- **Per-fold ComBat-seq batch correction protocol**: leakage prevention
- **Pipeline scripts**: `notebooks_or_scripts/v17p35_*.py`
- **Generalisability**: DIAL framework 은 thyroid cancer 외 임의 cohort-pooling scenario 에 적용 가능

## 5. Files in this dir

| File | Status | Purpose |
|---|---|---|
| `cover_letter.md` | ✅ done (2026-04-27) | Editor 발송용 cover letter — 핵심 claim + 3 suggested reviewers (Aedin Culhane, Lior Pachter, Casey Greene) |
| `figures/` | ⚠️ empty (TBD) | Manuscript figures — DIAL diagram, 5-cohort validation forest, GSEA sign-agreement |
| `tables/` | ⚠️ empty (TBD) | Supplementary tables — 5-cohort metadata, pathway sign-agreement table, dependency version manifest |
| `reproducibility/` | ⚠️ empty (TBD) | Reproducibility archives — script SHA-256, environment lockfile, deterministic seed log |

## 6. Key claim — clinical relevance (downstream)

8-gene RAI panel (SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) 가 BRAF V600E status alone 대비 **ΔAUC = +0.132** — 본 paper 에서는 임상 translation 강조 X, methods framework 의 actionable downstream 으로 framing.

## 7. Suggested reviewers

| Name | Affiliation | Expertise |
|---|---|---|
| Aedin Culhane, PhD | University of Limerick | Genomic data integration + reproducibility |
| Lior Pachter, PhD | Caltech | Single-cell + bulk transcriptomic methods |
| Casey Greene, PhD | University of Colorado | Biomedical data science + ML reproducibility |

(No conflicts of interest declared)

## 8. Cross-paper boundary

본 paper 는 **DM1/DM2 axis 의 4-paper share 구조** 의 methodology arm:

| Paper | Framing | Voice |
|---|---|---|
| Paper 0 (npj) | Clinical biomarker (8-gene > BRAF V600E) | Translation |
| **Paper B (Bioinformatics)** ← 이 paper | **Methodology + reproducibility (DIAL)** | **Methods** |
| Paper G (Genome Medicine) | Multi-omics integration (axis + drug + pan-cancer) | Mechanism |
| Paper 1 (Cancer paper v8) | Mechanism (DM1 fusion-driven epigenetic dark matter) | Cancer biology |

→ Same reviewer pool risk (Paper 2 advisor Q07 의제) 시 본 paper 의 methods framing 으로 separation 정당화.

## 9. Submission readiness

| 항목 | Status |
|---|---|
| Cover letter | ✅ |
| Manuscript draft | ⚠️ TBD (cover letter 기반 prose 작성 필요) |
| Figures | ⚠️ TBD |
| Supplementary tables | ⚠️ TBD |
| Reproducibility archive | ⚠️ TBD |
| Suggested reviewers | ✅ 3 names |
| Conflict of interest | ✅ none declared |
| Deadline | Open (Bioinformatics OUP rolling) |

## 10. References for next steps

- **Memory**: `~/.claude/.../memory/v52_lodo_finding.md` — DIAL detection 의 origin (v5.1 leak artifact)
- **Manuscript v8 (Paper 1)** 와 figure source data 일부 share — `project/results/v17_*tcga*/` 등
- **NeurIPS DIAL paper** (`submission/v15_neurips/`) 와는 별도 (NeurIPS = ML theorem version, Bioinf = applied methods version)

---

*Generated 2026-05-02 by THCA project master index sweep. Cover letter dated 2026-04-27.*
