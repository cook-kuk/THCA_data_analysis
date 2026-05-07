# Paper 3 ICI — 1-page 데이터 요약 (한국어, 교수님 1분 브리핑용)

**작성일:** 2026-05-06 / **저자:** Seungho Cook

---

### 무엇을 했나

23 개 데이터셋을 1 개 마스터 TSV 에 등록 → 14 개 GEO series matrix 다운로드 (~133 MB) → 612 개 sample metadata 추출 → 7 module × 60 gene 모듈 정의 → gene coverage 검증.

### 무엇을 분석할 수 있나 (Track B 진입 시)

- **Bulk thyroid TME (~388 sample):** TCGA-THCA + GSE126698 (n=28) + GSE76039 (37) + GSE60542 (92) + GSE151179 (52) + 4 microarray (105+49+25+92).
- **scRNA (158K+ cells):** GSE184362 + GSE193581 + GSE232237 + GSE191288 + GSE148673.
- **Pan-cancer ICI DIAL:** Hugo + Riaz processed 확보; IMvigor210 + Gide + Liu + Kim 추후.

### 임상 근거 (citation only)

KEYNOTE-158, KEYNOTE-028, NCT03246958, Dierks 2021, Spartalizumab ATC, Atezo+TT ATC. 모두 raw RNA-seq 비공개 → **분석 입력 아님, 인용만**.

### Controlled access (신청 필요, 4–8 주)

- dbGaP TCGA-THCA paired BAM (LOHHLA용)
- dbGaP Liu phs000452 (DIAL Tier 2)
- EGA Yoo 2019 Korean ATC

### 변경되지 않는 한계 (정직)

> 방어 가능한 endpoint = **ICI-readiness / immunogenomic vulnerability prioritization** 까지.
> **"thyroid 에서 ICI 반응 예측"은 raw thyroid ICI-treated RNA-seq 부재로 데이터상 불가능.**
> 해제 경로 = **연세 / 서울대 / 분당 / 삼성** 비공개 cohort 협력만.

### 다음 우선순위 옵션

| | 옵션 | 마라톤 비위반 | 필요 기간 |
|---|---|---|---|
| A | Track A 추가 강화 (사전등록 / forbidden words / negative-result 시나리오) | ✅ | 1–5 일 |
| B | dbGaP / EGA 신청 sprint | ✅ | 4–8 주 lag |
| C | 비공개 thyroid ICI 협상 (교수님) | ✅ | months — K1 해제의 *유일한* 경로 |
| D | Track B 진입 (Module A bulk ecotype 시작) | Paper 1 출하 후 | 자연 트리거 |

### 가장 강한 허용 결론 (binding)

> 공개 thyroid TME 데이터와 ICI 임상 evidence 통합으로 Paper 3 를 강화할 수 있으나, 방어 가능한 endpoint 는 여전히 **ICI-readiness / immunogenomic vulnerability prioritization** 이며 **thyroid ICI response prediction 은 아니다.**

---

산출물: `project/reports/paper3_ici/` (registry / log / blockers / plan / packet) + `project/results/paper3_ici_public_data/` (133 MB 다운로드 자산) + `project/results/paper3_ici_data_registry/` (metadata / QC / coverage TSV + python script).
