# EGA DAR Application — EGAD00001004845 (Yoo 2019 Advanced Korean Thyroid)

## 신청 가이드 (한국어)

### Step 1: EGA 계정 등록 (15분)
1. https://ega-archive.org/register 접속
2. 계정 정보 입력:
   - First name: Seungho
   - Last name: Cook
   - Email: kukshomr@gmail.com
   - Affiliation: Independent Researcher, Seoul, Republic of Korea (Graduate School of Convergence Science and Technology, Seoul National University)
3. 이메일 verification 후 활성화 (5분-1시간)

### Step 2: PI 동의서 (Yu Hyeong Won 교수님)
- 유형원 교수님 명의 PI consent 필요
- 본인 = co-investigator 등록
- PDF template: 아래 본문 활용

### Step 3: DAR (Data Access Request) submission
- URL: https://ega-archive.org/datasets/EGAD00001004845
- "Request Access" 버튼 → 아래 form 영문 본문 paste
- 24-48시간 내 DAC (Data Access Committee) 접수 confirm
- 4-8주 후 approval 또는 추가 정보 요청

---

## DAR Form Body (영문, EGA 제출용)

**Project Title**:
Korean Thyroid Cancer Cross-Population Validation of an 8-gene RAI-Responsiveness Decision Panel and Differentiation–Immune Sub-Stratification Framework

**Principal Investigator**:
Yu Hyeong Won, M.D., Ph.D.
Department of Surgery, Seoul National University Bundang Hospital
Seongnam, Republic of Korea
Email: [Yu Hyeong Won email TBD]

**Co-Investigator**:
Seungho Cook
Independent Researcher, Seoul, Republic of Korea
Graduate School of Convergence Science and Technology, Seoul National University (part-time PhD candidate)
Email: kukshomr@gmail.com

**Dataset Requested**:
EGAD00001004845 — Yoo SK et al. Nat Commun 2019 advanced Korean thyroid cancer cohort (n=25 RNA-seq, n=113 DNA-seq spanning PDTC, ATC, refractory PTC).

**Project Description (250 words)**:

We have developed an 8-gene radioactive iodine (RAI) responsiveness biomarker panel (NIS, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) for clinical pre-RAI decision-making in papillary thyroid cancer (PTC). The panel achieves cross-validated AUC 0.962 (R1-A leak-free re-validation) on TCGA-THCA (n=513) and outperforms BRAF V600E status alone (ΔAUC = +0.113). External validation on GSE76039 (PDTC + ATC, histology ground truth) yields AUC 0.935 (95% CI 0.824–1.000), demonstrating cross-platform transferability.

To strengthen Korean cohort validation, we have already incorporated PRJEB11591 (Yoo et al. PLoS Genet 2016) primary cohort RNA-seq data. We now request EGAD00001004845 (Yoo et al. Nat Commun 2019) for the following research aims:

1. **Korean PDTC + ATC validation**: Apply 8-gene panel to advanced Korean thyroid cancers to confirm direction-correct prediction across population and histology.
2. **Cross-population trajectory analysis**: Combine TCGA + GSE76039 + Yoo 2019 to refine the PTC → PDTC → ATC dedifferentiation trajectory across Western and Korean populations.
3. **TERT promoter prevalence**: Compare advanced-disease TERT prevalence between MSK-IMPACT (n=231 thyroid, ATC 81.8%) and Korean advanced cohort.
4. **8-gene panel performance in radioiodine-refractory subset**: Evaluate panel discrimination in the refractory subset of Yoo 2019.

Our manuscript is under preparation for npj Precision Oncology submission.

**Project Duration**: 12 months (from data access date)

**Data Access Type**: Research Use only (RU)

**Data Storage**:
- Encrypted local server (LUKS encryption, /opt/thyroid-dash/secured/)
- Access restricted to PI + Co-Investigator only
- Data destroyed after 12 months
- Re-application if extension needed

**Publication Plan**:
- Primary: npj Precision Oncology (under preparation, submission Q2 2026)
- Yoo SK et al. Nat Commun 2019 will be cited as primary reference
- Both Yoo 2016 and Yoo 2019 papers will be acknowledged in Methods/References

**Investigator Commitments**:
- Will not share data with any unapproved investigator
- Will not attempt to re-identify research participants
- Will destroy all data 12 months after access date (or upon re-application)
- Will cite Yoo SK et al. Nat Commun 2019 in any publication using these data
- Will report any adverse findings or unexpected results to the DAC
- Will respect any embargo or restriction notified by the DAC

**Conflict of Interest**: The investigators have no financial conflict of interest with Yoo SK et al. or the original data providers. No commercial use is intended.

---

## 본인 액션 체크리스트

- [ ] EGA 계정 등록 (15분)
- [ ] 유형원 교수님께 PI 명의 consent + 이메일 받기
- [ ] DAR form 의 [Yu Hyeong Won email TBD] 채우기
- [ ] EGAD00001004845 페이지에서 "Request Access" → form paste
- [ ] DAC 응답 wait (4-8주)
- [ ] Approval 받으면 → data download → revision-round 분석

---

## DAR 신청 후 follow-up timeline

| Date | 액션 |
|------|------|
| Day 0 | EGA 계정 등록 + DAR submit |
| Day 1-7 | DAC 접수 confirm 이메일 |
| Day 14-28 | DAC review (필요 시 추가 질문) |
| Day 28-56 | Approval 또는 reject |
| Approval 후 | Data download (EGA 의 ega-download client 사용) |
| Approval 후 +2주 | 8-gene panel 분석 + manuscript revision-round 통합 |
