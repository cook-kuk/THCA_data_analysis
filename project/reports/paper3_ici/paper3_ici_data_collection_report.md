# Paper 3 ICI — Data Collection Report (Sprint v2)

**작성일:** 2026-05-06
**작성자:** Seungho Cook
**Sprint 권한:** 사용자 명시적 권한 부여 — Track A "no public download" 일부 오버라이드 (small/moderate processed matrices만). Track B의 NMF / LOHHLA / NetMHCpan / DIAL 실행은 여전히 차단.
**Claim guard (binding):** allowed = ICI-readiness, immunogenomic vulnerability, candidate immune ecotype, hypothesis-generating prioritization. **forbidden = ICI response predictor in thyroid, clinical treatment recommendation, patient selection biomarker, validated ICI sensitivity.**

---

## 1. 결과 요약 (1줄)

23개 데이터셋을 v2 마스터 registry에 등록 — 10개 thyroid TME (public, 분석 가능), 7개 thyroid ICI 임상 evidence (citation-only), 6개 pan-cancer ICI (DIAL anchor, response modeling은 pan-cancer에서만).

---

## 2. 등록 데이터셋 카운트

| 그룹 | 수 | 분석 가능 여부 |
|---|---|---|
| Thyroid TME — 공개 RNA-seq / scRNA / microarray / spatial | 10 | ✅ 모듈 스코어링 + ecotype 후보 + dedifferentiation 축 |
| Thyroid TME — 내부 read-only (Paper 1/2 자산) | 3 (PRJEB11591, PRJKA210106, GSE250521) | △ cross-paper boundary 안에서만 read-only 참조 |
| Thyroid ICI — 임상 evidence | 7 (KEYNOTE-158, KEYNOTE-028, NCT03246958, Dierks 2021, Spartalizumab, Atezo+TT, Lenva-Pembro RAI-r) | ❌ raw RNA-seq/WES 비공개 — citation only, response prediction 클레임 불가 |
| Pan-cancer ICI — DIAL audit | 6 (IMvigor210, Hugo, Riaz, Gide, Liu, Kim GC, Cho NSCLC) | ✅ DIAL audit anchor (response modeling은 pan-cancer 한정) |

---

## 3. Track A 자산과의 차이 (v1 → v2)

| 항목 | v1 (Track A 동결) | v2 (이번 sprint) |
|---|---|---|
| 데이터셋 수 | ~16 (raw count of bulk + scRNA + spatial + ICI references) | 23 명시 등록 |
| Thyroid TME RNA-seq 추가 | TCGA-THCA + GSE76039 + Yoo 2019 | + GSE126698 (n=365) — 가장 큰 cross-subtype 추가 |
| Thyroid scRNA 추가 | GSE184362 / GSE193581 / GSE232237 / GSE191288 / GSE148673 / Han 2024 (모두 to_verify) | GSE184362 + GSE193581 우선순위 1로 명시 (cell counts + cohort 구조 확정) |
| Thyroid RAI 축 | 없음 | + GSE151179 / GSE151180 (RAI-refractory vs RAI-avid 비교) |
| Thyroid ICI 임상 evidence | 없음 | + 7개 임상 evidence — 페이퍼 Discussion에 ICI 임상 맥락 추가 |
| Microarray cross-platform | GSE33630 / GSE65144 / GSE54958 / GSE53157 / GSE29265 (Track A) | 동일 + GSE60542 추가 |

---

## 4. 분석 가능성 매트릭스 (per-dataset 요약)

### 4.1 즉시 모듈 스코어링 가능 (Track A 동결 안에서도 — 단 *지금* 실행은 마라톤 위반)

- **TCGA-THCA** — 1차 디스커버리 코호트. RNA + WES + clinical. paired BAM은 dbGaP 게이트.
- **GSE76039** — PDTC/ATC dedifferentiation anchor. n=37 작지만 TCGA aggressive subset과 페어링.
- **GSE126698** — n=365 cross-subtype. immune-module 스코어링 + ecotype replication 가능.
- **GSE184362** — scRNA n=11 patient / 158,577 cells. 1차 atlas 후보. localized + RAI-refractory metastasis 모두 포함.
- **GSE193581** — ATC transformation scRNA. 10 ATC + 7 PTC + 6 normal + 9 ATC cell lines bulk. ATC scRNA의 핵심 자원.
- **GSE151179 (+ GSE151180 supplement)** — RAI-refractory vs RAI-avid PTC. aggressive 표현형 proxy.
- **GSE33630 / GSE29265 / GSE65144 / GSE60542** — microarray cross-platform replication.
- **IMvigor210 / Hugo / Riaz / Gide / Kim GC / Cho NSCLC** — pan-cancer DIAL audit anchors. **response modeling은 여기서만 허용**.

### 4.2 분석 차단 또는 readiness만 가능

- **Liu melanoma (phs000452)** — dbGaP 신청 필요. controlled.
- **TCGA-THCA paired BAM (LOHHLA용)** — dbGaP 신청 필요.
- **Yoo 2019 Korean ATC** — EGA-controlled 가능성. v1 to_verify 유지.

### 4.3 Citation only (raw data 비공개)

- **KEYNOTE-158 thyroid cohort** — Maio 2022 J Clin Oncol. RECIST + PD-L1 TPS 보고. raw RNA-seq 없음.
- **KEYNOTE-028 thyroid cohort** — Mehnert 2019 BMC Cancer. PD-L1+ n=22, ORR ~9%. raw 비공개.
- **NCT03246958 (Sehgal/Lorch nivo+ipi)** — DTC/ATC/MTC 분리 코호트. ASCO supplement만.
- **Dierks 2021 lenva+pembro ATC/PDTC** — Thyroid 31(7):1076. WES summary supplement.
- **Spartalizumab ATC (Capdevila 2020 JCO)** — n=42 ATC, ORR ~19%, PD-L1 cutoff 분석.
- **Atezo + matched TT ATC (Cabanillas, NCT03181100)** — 코스 구분 결과.
- **Lenva+Pembro RAI-refractory DTC** — 다수 trial 진행 중.

→ 이들은 Discussion 1.1 hook + ICI 임상 맥락 cite 용도. **response prediction 클레임의 evidence 아님.**

---

## 5. Cross-paper boundary 재확인 (v2 sprint 이후)

| 데이터셋 | Paper 1 (DM1) | Paper 2 (HT-isolated) | Paper 3 (ICI vulnerability) | Paper 4 (GD HLA backlog) |
|---|---|---|---|---|
| TCGA-THCA | DM1 driver/cluster 본진 | HT-overlap subset | dark-matter (BRAF/RAS-neg) + HLA LOH + neoantigen | 인용만 |
| PRJEB11591 | 사용 | Paper 2 Pillar 1 v2 base | Korean HLA frequency anchor read-only | 인용만 |
| PRJKA210106 | 사용 | 사용 | HLA typing reference만 | 미사용 |
| GSE126698 | dark matter prevalence cite 가능 | 미사용 | ecotype replication + dedifferentiation cross-subtype | 미사용 |
| GSE184362 / GSE193581 (scRNA) | 미사용 | 미사용 | scRNA atlas + sub-state + CXCL13/TLS | 미사용 |
| GSE151179 | 미사용 | 미사용 | RAI-refractory aggressive phenotype proxy | 미사용 |
| Microarray (GSE33630 외) | 미사용 | 미사용 | cross-platform replication | 미사용 |
| GSE250521 spatial | Paper 1 supplement 본진 | 미사용 | TLS/CXCL13 overlay only | 미사용 |
| Pan-cancer ICI (6) | 미사용 | 미사용 | DIAL audit anchors | 미사용 |
| Thyroid ICI 임상 evidence (7) | 미사용 | 미사용 | citation only — Discussion + Fig 1A 임상 맥락 | 미사용 |

**Paper 1 / Paper 2 분석은 이번 sprint에서 건드리지 않음.** Paper 1/2 ETL outputs은 read-only 참조만.

---

## 6. v2 sprint가 페이퍼에 가져오는 강화 포인트

1. **GSE126698 (n=365)** — 단일 가장 큰 cross-subtype thyroid TME RNA-seq 추가. ecotype 발견의 statistical power ↑, microarray + RNA-seq 양쪽 cross-platform replication 가능.
2. **GSE193581 ATC transformation scRNA** — Track A v1에 없던 ATC scRNA 자원. dedifferentiation trajectory 의 정성적 증거 강화.
3. **RAI-refractory 축 (GSE151179)** — aggressive 표현형 proxy 추가. Module A ecotype × RAI-status 교차 분석 가능 → Module E readiness score의 RAI-refractory subset 검증 layer.
4. **Thyroid ICI 임상 evidence 7개** — Paper 3 Introduction hook + Discussion 임상 맥락 cite. "thyroid에서 ICI는 어떤 효과를 보이는가" 의 baseline benchmark.
5. **DIAL anchor 6개 명시** — pan-cancer response modeling 의 정량적 근거. Paper 3 Module D의 statistical machinery 입력 명확화.

→ **모든 강화는 readiness/vulnerability/prioritization framing 안에서만.** 어떤 강화도 "ICI response predictor in thyroid" 클레임을 정당화하지 않음.

---

## 7. 결론 (강화 가능성 + 한계 동시 명시)

> Paper 3는 v2 sprint를 통해 등록 데이터의 폭과 깊이가 모두 확장되었다. 그러나 **현재 방어 가능한 endpoint는 여전히 "ICI-readiness / immunogenomic vulnerability prioritization" 이며, "thyroid ICI response prediction" 은 아니다.** 이 한계는 raw thyroid ICI-treated RNA-seq 의 부재라는 단일 근본 원인에 의해 결정되며, sprint v2의 어떤 추가 데이터도 이 한계를 변경하지 않는다.

---

Track A 동결 (Paper 1/2 분석) 유지. 마라톤 비위반 (registry + manifest + plan만; Track B 분석 실행 없음).
