# GPU Needs Assessment — RunPod 사용 여부

**Date:** 2026-05-04 (W1 진입 직전)
**Marathon mode rule (from memory):** "새 분석은 paper-blocking 만". GPU work add 시 paper-blocking 여부 확인 필수.

---

## ★ 정직한 진단

**Phase 0 paper (Cell Rep Med / JCI Insight target) 의 모든 5-Pillar 분석은 CPU 만으로 완성됨.**

- Total Azure cost: $5.80 (CPU burst만)
- All 11 smoke tests pass (CPU)
- All 5 main + 6 suppl figures rendered (matplotlib, CPU)
- All 11 suppl TSVs built (CPU)

**현재 단계 (manuscript writing 5/4-6/13):**
→ **GPU 필요 없음.** 본인 voice-protected sections 작성 + Yu professor 미팅 + bioRxiv 6/13 submission 까지 GPU work 0건.

---

## 1. Phase 0 paper-blocking GPU 작업 — 없음

| Component | Status | Compute |
|---|---|---|
| arcasHLA RNA-seq HLA imputation | ✅ done (5/2 burst $1) | CPU (32-core, Azure D32s_v5) |
| PyDESeq2 differential expression | ✅ done | CPU |
| gseapy GSEA pre-ranked | ✅ done | CPU |
| KMeans cluster + ARI | ✅ done | CPU |
| Random-effects DerSimonian-Laird forest | ✅ done | CPU |
| Baron-Kenny mediation + bootstrap | ✅ done | CPU |
| BCR repertoire gene-level proxy | ✅ done | CPU |
| 5 main + 6 suppl figures | ✅ done | CPU (matplotlib) |
| Smoke tests 11/11 | ✅ pass | CPU |

→ **결론: Phase 0 manuscript 6/13 bioRxiv 제출 까지 GPU instance 필요 없음.**

---

## 2. GPU 가 *유의미하게 도와줄 수 있는* 작업 (paper-strengthening, NOT paper-blocking)

마라톤 rule: 모두 추가 가능하지만 **paper-blocking 여부 본인 확인 후만 진행**.

### 2.1 ★ 가장 likely 후보 — TCGA-THCA Pathology Slide Multimodal Layer

**Why genuine GPU need:**
- TCGA-THCA 에 ~500 WSI (whole-slide images) 공개 (~50-100 GB)
- Hashimoto-overlap 탐지 (lymphocytic infiltration pattern) 를 deep learning 으로 정량
- 본인 메모리 `user_image_skills.md`: "User does medical image processing including pathology"
- 본인 자산: medical imaging expertise + cookHLA Nat Commun first-author identity

**Could strengthen Pillar 5:**
- 현재 Pillar 5 mechanism: transcriptomic + HLA-II + BCR/TLS
- 추가: H&E pathology image based Hashimoto-overlap 정량 (orthogonal modality)
- "Multimodal Hashimoto-overlap detection" claim — Cell Rep Med fit 더 강해짐

**Resource estimate:**
- WSI download: ~50-100 GB from GDC (CPU/network only)
- WSI tiling + feature extraction (CTransPath / UNI / etc): RunPod A100 ~10-20 hr
- Hashimoto classifier training: ~5-10 hr A100
- Inference on 500 TCGA samples: ~2-4 hr A100
- **Total: ~$30-60 RunPod A100 ($1.69-2.89/hr)** + ~1 week wall-clock

**Decision:** **Phase 0 paper-blocking 아님.** 본인이 "Cell Rep Med vs Nat Commun reach" 기준 무조건 strengthen 하고 싶으면 OK; 그 외에는 Phase 1 trajectory 로 미루는 게 marathon 정신.

### 2.2 AlphaFold HLA-DPB1\*05:01 + TSHR autoreactive peptide structure

**Why GPU:**
- AlphaFold 2/3 추론 시 NVIDIA A100 (40-80 GB VRAM) 필요
- HLA-DPB1\*05:01 + TSHR/Tg/TPO peptide complex structure prediction
- Pillar 1 (Korean DPB1\*05:01 53%) + Pillar 5 (autoimmune mechanism) 의 structural mechanism 보강

**Resource estimate:**
- AlphaFold MSA generation: ~2-4 hr A100 per structure
- Inference: ~30 min per structure × 5-10 structures = 5-10 hr
- **Total: ~$10-20 RunPod A100**

**Decision:** **Phase 0 paper-blocking 아님.** Phase 1 (Graves' paper) 의 mechanistic supplementary figure 로 더 적합. 마라톤 mode 정신상 제외.

### 2.3 scBCR-seq full clonotype analysis

**Why GPU:**
- 만약 single-cell BCR-seq 데이터 추가 (Bundang Scenario A response 시 가능)
- immcantation-style clonotype clustering 은 large-scale 시 GPU 가속 유의미
- 현재 paper 는 gene-level BCR proxy 만 (제한 명시: Methods M7)

**Resource estimate:**
- scBCR-seq pipeline: ~5-10 hr A100 per cohort
- 만약 Bundang n=50+ scBCR-seq → ~$15-30

**Decision:** **Phase 0 paper-blocking 아님 (data 자체 없음).** Phase 1 Bundang trajectory 로.

### 2.4 K2 STAR re-quantification (D7-P3 deferred)

**Why mentioned:**
- D7-P3 calibration FAIL 의 future task
- 메모리 `v17_korean_k2_calibration.md`: "STAR re-quant 은 future task ($5-10, post-paper draft)"

**Compute:**
- **STAR alignment 은 CPU job, NOT GPU.** 32-core 1-2 day enough.
- RunPod CPU instance ($0.5-1/hr × 24-48hr = $12-50)
- **GPU 안 필요** — STAR 는 GPU 가속 안 함

**Decision:** **Phase 0 paper-blocking 아님.** Future task 그대로. RunPod CPU spot 사용 시 가능.

---

## 3. RunPod 사용 가이드 (만약 위 중 하나 진행 결정 시)

### 3.1 RunPod overview
- Spot pricing: A100 80GB ~$1.69/hr (vs on-demand ~$2.89/hr)
- Container-based: pre-built ML images (PyTorch, TensorFlow, JAX)
- Storage: $0.10/GB/month persistent volume
- Korean network: latency 100-200 ms (vs Azure Korea Central <10 ms)

### 3.2 Account setup
```
1. https://www.runpod.io 접속
2. 신용카드 등록 (Azure 계정과 별개)
3. API key 발급 → ~/.runpod/config
4. 본인 SSH key upload
```

### 3.3 Pod 선택 (작업별)

| 작업 | Pod 추천 | 시간당 cost | 예상 시간 |
|---|---|---|---|
| **WSI pathology pipeline (CTransPath)** | A40 48GB ($0.79/hr spot) | $0.79/hr | 20-40 hr |
| **AlphaFold HLA-peptide** | A100 80GB ($1.69/hr spot) | $1.69/hr | 5-10 hr |
| **scBCR-seq immcantation** | A40 48GB | $0.79/hr | 5-10 hr |
| **STAR re-quantification** | CPU 32-core (no GPU) | $0.18/hr | 24-48 hr |

### 3.4 Spot vs On-demand
- Spot: 60% 할인. eviction risk (~5-10% 시간당)
- 이전 Azure burst pattern (auto-shutdown + checkpoint) 활용 시 spot OK
- Critical run (마감 직전) 은 on-demand

### 3.5 Cost ceiling
- Phase 0 marathon: $0 권장 (마라톤 정신)
- 만약 Pathology multimodal Phase 0 추가: $30-60 ceiling
- AlphaFold Phase 0 supplement: $10-20 ceiling
- **모두 합쳐도 $100 안에**

---

## 4. 본인 결정 후보

### 후보 A — RunPod 사용 안 함 (마라톤 정신 그대로)
- Phase 0 paper 5-pillar 모두 STRONG 이미 확정
- 6/13 bioRxiv submission 까지 무중단 voice-protected 작업
- Phase 1 Bundang trajectory 진입 시 RunPod 사용 검토

### 후보 B — Pathology multimodal Phase 0 추가 (paper-strengthening)
- TCGA-THCA WSI 500개 다운로드 + CTransPath embedding + Hashimoto classifier
- Pillar 5 의 multimodal claim 추가 (transcriptomic + HLA + image)
- ~$30-60 RunPod, 1 week wall-clock
- venue ladder: Cell Rep Med → Nat Commun reach 가능성 증가
- **Marathon rule 위반 가능성**: 새 분석. "paper-blocking 인가?" — Cell Rep Med submission 만 보면 No, Nat Commun 노린다면 Yes 가능

### 후보 C — AlphaFold structural mechanism 추가 (mechanistic supplement)
- HLA-DPB1\*05:01 + TSHR/Tg peptide structure prediction
- Pillar 1 + Pillar 5 의 structural supplementary figure
- ~$10-20 RunPod
- **Marathon rule 위반 가능성**: 새 분석. paper 의 review-defense 강화 정도.

### 후보 D — STAR re-quantification 의 K2 calibration fix (D7-P3 closure)
- D7-P3 의 future task 즉시 closure
- ~$12-50 RunPod CPU (NOT GPU)
- Pillar 1 의 K2 cohort raw-score meta 진입 가능
- **Marathon rule 위반 가능성**: 분석. 그러나 D7-P3 명시적 future task 라 trajectory 내. closure 가치 있음.

---

## 5. 권고

**마라톤 mode 정신 그대로 = 후보 A (사용 안 함).**

이유:
- Phase 0 paper 5-pillar 모두 STRONG ✅
- 모든 scaffolding/infra 100% 구축 완료 ✅
- 6/13 bioRxiv submission 까지 voice-protected 본인 작업만 남음
- 추가 GPU 작업은 W3-W6 (drafting + revision) phase 에 momentum 손실
- "marathon 정신: 새 분석은 paper-blocking 만"

**예외 — 만약 본인 며칠 전 결정한 후보 F (Hybrid: Cell Rep Med direct submission + 분당 outreach + Phase 1 outline) trajectory 와 정합 시:**
- Phase 0 단계: 후보 A (GPU 사용 안 함)
- Bundang Scenario A 응답 시: 후보 D (K2 STAR re-quant) → Phase 0 paper Pillar 1 강화 + Bundang validation
- Phase 1 entry 시: 후보 B (Pathology multimodal) + 후보 C (AlphaFold) → Phase 1 paper 강화

---

## 6. 본인이 진짜 묻고 싶은 것 (추측)

만약 "GPU 작업 없는 게 정상 맞나?" 라는 의도라면:
- ★ **YES, 정상.** Phase 0 paper 의 모든 분석은 CPU 충분.
- HLA imputation, RNA-seq DE, GSEA, mediation, KMeans, forest meta — 모두 CPU.
- GPU 가 paper-blocking 영역은 Phase 1 (Bundang scBCR-seq) + Phase 2 (DIAL audit) 트래젝토리에 있음.

만약 "WSI 또는 AlphaFold 추가 해야 venue 더 강해지나?" 라면:
- 마라톤 정신: NO. 6/13 bioRxiv → Cell Rep Med Q3 2026 submission 정렬.
- 추가 시: Nat Commun stretch 가능성 증가 but 시간 비용 큼.

---

## 7. 본인이 결정해야 할 것

- [ ] 후보 A/B/C/D 중 선택 (default = A)
- [ ] B/C/D 선택 시: marathon rule 의 "paper-blocking 여부" 본인 명시 정의
- [ ] RunPod 계정 setup (만약 B/C/D 진행 시)
- [ ] Yu professor 미팅 시 GPU 작업 진행 여부 advisor 의견 확인

---

**한 줄 요약:**
**현재 Phase 0 paper 는 GPU 필요 없음. CPU 만으로 5-pillar 모두 STRONG. 6/13 bioRxiv submission 까지 RunPod 사용 안 권장 (마라톤 정신). Phase 1 Bundang Scenario A 응답 시 후보 D (K2 STAR) 또는 B (WSI multimodal) 검토.**
