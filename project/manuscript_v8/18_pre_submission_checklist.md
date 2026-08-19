# Paper 1 사전 제출 체크리스트 — 2026-05-13

> Yu 미팅 (2026-05-14) 준비용. Marathon mode 9일차 / 30일 잔여 → 2026-06-13.
> Audit GREEN (2026-05-13): boundary clean (0 Paper 2 leaks), n=865, key numbers locked (TERT 36 muts, fusion 76.8% OR 7.41, pooled HR 2.53, TPO d=2.30, β 0.385 vs 0.253, ARI 0.49/0.90/0.92).
> Target venue: **Cell Reports Medicine 1순위** → JCI Insight / Nat Commun dual reach → npj Precision Oncology fallback.

---

## A. Voice-protected (본인 키보드 필수, Claude 생성 금지)

- [ ] `03_introduction.md` §1.1 hook — **Alt B refined** ("intermediate-risk ~20%" + "~23% BRAF/RAS-neg") 또는 **hybrid** ("compass" metaphor 보존)
- [ ] `04_intro_1_1_hook.md` author hook anchor 확정 (Alt B vs hybrid 결정 + 1문장 lock)
- [ ] `06_discussion.md` §3.1 opening paragraph — **Krishnamoorthy 2025 → Landa 2016 JCI 정정 cite** (8-gene 5/8 ATC overlap reverse-causality 차단 3-layer 본인 표현)
- [ ] `06_discussion.md` §3.4 limitations — 4 main (cohort/biology/clinical/methylation) + 2 supp (image DM1 confound + Track B-lite RAI)
- [ ] `08_cover_letter.md` paragraph 1 — Cell Rep Med editor hook 본인 voice
- [ ] `09_reviewer_qa.md` Q9 — reviewer hardest question 본인 답변

## B. Figure build (technical, Claude 보조 가능)

- [ ] **Fig 1** (4 panels — Heatmap C → S1 demoted v3): TIERA67 7-category bar / 8-gene PCA / DM1-DM2 split / Driver_anchor leak-control
- [ ] **Fig 2** (clinical aggressiveness): TCGA HR 2.30 forest + MSK HR 2.67 forest
- [ ] **Fig 3** (sc external): Pu 2021 + Lu 2023 + FFPE KS p=0.44
- [ ] **Fig 4** (BRAF/TERT orthogonality)
- [ ] **Fig 5** (Korean K2+Lee n=865 validation + reflex algorithm)
- [ ] **Fig 6** (pooled HR 2.53 meta-analysis forest)
- [ ] **Fig 7** A-C (fusion 76.8% + OR 7.41) + 7D (TPO/DIO1/TSHR β bar → S6 demoted)
- [ ] **Fig 8** (3 panels): methylation heatmap + mean β bar + fusion-independent

## C. Citation lock-in

- [ ] **Landa 2016 JCI 126(3):1052-1066** — bib 추가 + Disc 3.1 정정 (Krishnamoorthy 2025 misattribution 차단 핵심)
- [ ] **Yoo 2016 PLOS Genet 12(8):e1006239** — panel origin cite Methods + Intro 1.2 + Disc 3.1
- [ ] **TCGA 2014 Cell 159(3):676-690** — Intro 1.2 canonical PTC genomic landscape
- [ ] **Pu 2021 Nat Commun** — Methods sc + Disc 3.3
- [ ] **Bradley 2010** — Disc 3.1 BRAF V600E HLA-I counter-intuitive anchor
- [ ] **Wirth 2020 NEJM** — Disc 3.2 LIBRETTO-001 (selpercatinib RET)
- [ ] **Haugen 2016 + Ringel 2025 Thyroid** — ATA 2015 + 2025 dual cite (reflex algorithm 정당화)

## D. Submission infra

- [ ] **Public code repo URL** — Zenodo DOI 확정 (GitHub release tag → Zenodo mint)
- [ ] **STAR Methods data availability** — GEO/dbGaP/cBioPortal IDs 최종 (Korean K2 PRJEB11591, GSE286332, MSK-IMPACT cBioPortal study ID)
- [ ] **Cover letter** 소속/이메일/corresponding 채우기
- [ ] **Title 최종 lock**: "An 8-gene panel reveals fusion-driven, epigenetically silenced dark matter in BRAF/RAS-negative thyroid cancer"
- [ ] **Reviewer suggestions** 3-5명 (BRAF/RAS-neg 전문 + epigenetic silencing + thyroid driver landscape)

## E. Bridge to bioRxiv

- [ ] **bioRxiv 시점 결정** (Yu 합의 — Cell Rep Med 동시? 사전?)
- [ ] **preprint server**: bioRxiv vs medRxiv (clinical 비중 따라 — 본 paper는 mechanism 비중 높으므로 bioRxiv 권장)
- [ ] **preprint license**: CC BY (default)
- [ ] **ORCID + corresponding author lock**: Seungho Cook ORCID 확인, Yu lab corresponding

## F. Marathon 잔여 (9일차 / 30일 잔여 → 2026-06-13)

| 날짜 | 마일스톤 |
|---|---|
| **5/14 (내일)** | Yu 미팅 + voice-protected 4섹션 (Hook/Disc 3.1/Limitations/Cover ¶1) 착수 |
| 5/15-5/20 | 나머지 voice 섹션 (Q9 + hook anchor lock) + bib 최종 정리 |
| 5/21-5/27 | Figure build 8종 (외부 협력자 필요 여부 5/14 미팅서 확정) |
| 5/28-6/05 | **bioRxiv 제출 (목표)** + Zenodo DOI mint |
| 6/06-6/13 | **Cell Rep Med 제출 마무리** + reviewer suggestions 확정 |

---

## Yu 미팅 5/14 어젠다 (decision-ready)

1. **Hook 결정**: Alt B refined vs hybrid (5분)
2. **bioRxiv 시점**: Cell Rep Med 동시 vs 사전 (5분)
3. **Figure build 분담**: 외부 협력자 (illustrator?) 필요 여부 (5분)
4. **Reviewer suggestions 3-5명** 후보 brainstorm (5분)
5. **Cover letter editor hook 방향** 합의 (5분)
