# K2 H&E retrieval — push email draft

**Audience:** the K2 cohort PI / SNUH-GMI / whoever currently holds K2 H&E slides
**Goal:** establish a concrete timeline (or block) for K2 H&E retrieval
**Why now:** Paper 2 image-DM1 audit (2026-05-09) confirmed TCGA-only signal is ≥50% TSS-confound — K2 with different TSS structure is the only path to a high-impact image-DM1 paper

---

## Korean draft (formal, single-paragraph, ~150 words)

> 안녕하세요 [받는 분 호칭],
>
> 지난 4월 30일 K2 코호트 PRJEB11591 RNA-seq 분석 (n=260, arcasHLA HLA imputation 포함) 까지는 잘 진행되어 결과 정리 중에 있습니다. 다만 Paper 2 (H&E pathology → DM1 분자아형 분류) 의 TCGA pilot 단계에서 site-of-origin (TSS) 배치 효과가 우세하게 작용하여 단일 코호트에서는 정직한 검증이 어렵다는 것이 audit 결과 확인되었습니다 (보고서: REPORT_YU_REVIEW_PAPER2_2026_05_08.html, audit_ras_auc100/AUDIT_REPORT_v2.md). 이를 해결하기 위해 K2 cohort의 H&E whole-slide image (WSI) retrieval가 필요한 상황입니다. 현재 K2 H&E 슬라이드의 가용성과 retrieval 일정 관련하여 (i) 보유 여부와 디지털 스캔 상태, (ii) 데이터 사용 동의 범위, (iii) 예상 timeline (Q3 2026 vs 그 이후) 세 가지를 알고 싶습니다. 빠른 답변 부탁드리며, 필요하시면 화상 미팅 잡겠습니다.
>
> 감사합니다.
> Seungho Cook

## English draft (formal, single-paragraph, ~150 words)

> Dear [recipient],
>
> Following the successful K2 cohort RNA-seq analysis (PRJEB11591, n=260, arcasHLA-imputed HLA typing) on April 30, our Paper 2 image-DM1 work has hit a confound that requires K2 H&E to resolve. A 12-test audit (2026-05-09) of the TCGA-THCA pilot (n=59, UNI foundation model + CLAM) showed that the apparent subgroup AUC = 1.000 in RAS-like cases collapses to 0.308 under leave-one-TSS-out cross-validation, and that a clinical-only logistic regression (histology + sex + TSS) attains AUC = 0.768, exceeding the foundation model. The image signal is structurally entangled with TCGA's tissue-source-site sampling and cannot be honestly validated on a single cohort. K2 H&E with a different TSS distribution would provide the external validation needed for a high-impact image-DM1 paper. Could you confirm (i) whether K2 H&E whole-slide images (or FFPE blocks) are available, (ii) the digital scanning status, (iii) the data-use consent scope, and (iv) the expected retrieval timeline (Q3 2026 vs later)? Happy to set up a video call if helpful.
>
> Best regards,
> Seungho Cook

---

## Audit summary attachment (1-page, optional)

If the recipient wants quick context, attach a 1-pager from `AUDIT_REPORT_v2.md` Sections A + I (data structure crosstab + final synthesis table).

---

## Send timing

- **Don't send** until user confirms the right recipient and reads the drafts.
- Marathon mode: outreach actions are author-keyboard. Drafts here are for user to edit/send themselves.
- Memory `v17_npj_ship_status` precedent: 4 outreach drafts written but NOT sent (user sends).
