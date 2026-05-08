# Causal K-Thyro Foundation — 솔직한 약점 감사 + 다음 LLM prompt

**작성일:** 2026-05-08
**대상 과제:** 삼성미래기술육성사업 30억 Technology 트랙
**과제명:** Causal K-Thyro Foundation — 공개 고재사용 바이오데이터 사전학습과 한국 능동감시–수술 코호트 기반 반사실적 갑상선암 환자여정 파운데이션 모델
**PI 조합:** 유형원 (분당SNUH 내분비외과) × 국승호 (Co-PI, agentic AI / multi-omics / QML)

---

## 컨텍스트

ChatGPT와 길게 brainstorm 한 끝에 Causal K-Thyro Foundation으로 방향 확정 직전. ChatGPT는 선정 가능성을 45–65%로 추정. 다만 ChatGPT가 끝까지 안 짚은 진짜 약점이 4+1개 있어, 이걸 정리하고 다음 LLM에 던질 prompt까지 만든다.

본 과제 방향은 옳다 — Hernán–Robins target trial emulation은 정통 인과추론 방법이고, 한국 MAeSTro/KoMPASS는 실제 그 데이터 구조를 가진다. 그러나 아래 약점을 먼저 닫지 않으면 reviewer 2가 통째로 흔들 수 있다.

---

## 약점 1. Target trial emulation × small-N의 통계적 한계

**핵심:** Hernán–Robins target trial emulation은 보통 EHR 수만~수십만 명 + 충분한 outcome event 수에서 정당화된다.

- 저위험 PTC는 본질적으로 outcome rate가 낮음 (5년 progression ~5%, mortality ~0%)
- MAeSTro-EXP가 다기관 prospective여도 progression/recurrence/mortality event 수가 수십~수백 단위
- Doubly robust learner, causal forest, IPW 모두 분산이 폭발함
- CTC/EMT 62명으로 counterfactual? 신뢰구간이 사실상 (-∞, +∞)에 가까울 수 있음

**Reviewer 2 공격 포인트:**
> "Effective sample size for AS→progression counterfactual estimand가 몇입니까? Sensitivity analysis에서 unmeasured confounding이 어느 정도면 결과가 뒤집힙니까?"

답이 약하면 과제 통째로 흔들린다.

---

## 약점 2. "Foundation model"이라는 단어가 실제 N과 안 맞음

- TCGA-THCA n≈500, CPTAC THCA n≈300, GEO thyroid scRNA 합쳐서 수만 cells
- 진짜 foundation은 cross-cancer generalization을 보여줘야 함 (TCGA-PanCancer 11k, scRNA atlas 33M cells급)
- Single-disease(PTC)에 묶으면 "domain-adapted multimodal model"이지 foundation 아님

**Reviewer 공격:**
> "왜 이게 foundation model입니까? Pan-cancer pretraining 후 thyroid fine-tuning이라면 foundation이지만, thyroid only는 domain model입니다."

**해결:** pan-cancer pretraining (TCGA 11k pan-cancer + HTAN + CELLxGENE)을 WP1에 명시하고 "PTC는 first downstream task"로 framing. ChatGPT가 짠 WP1은 이게 약함.

---

## 약점 3. AS vs Surgery에서 "treatment"가 monitoring schedule과 confounded

가장 미묘하고 위험한 약점.

- AS군: 6개월마다 ultrasound + thyroglobulin + 환자 anxiety 관리 + follow-up 자주
- Surgery군: 즉시 수술 + post-op follow-up은 다른 패턴
- "Counterfactual progression"을 계산할 때, AS군이 progression을 더 잘 발견하는 건 "AS treatment 효과"인지 "monitoring intensity 효과"인지 구분 불가
- IV, propensity score, doubly robust 모두 안 풀림 — measurement schedule이 outcome ascertainment에 직접 영향

이건 정통 causal inference 문헌에서 **detection bias / surveillance bias**라고 부르는 고전적 함정이다. AS 연구에서 반드시 나오는 질문이고, 답하기 어렵다.

**해결:** target trial protocol에 "monitoring schedule을 둘 다에 동일하게 emulate"한다고 명시 + "informative censoring under shared monitoring"으로 sensitivity analysis. 이걸 첫 페이지에 박아야 reviewer가 안심.

---

## 약점 4. Quantum-inspired sidecar의 "왜?"가 약함

- Tensor network, MPS, quantum kernel은 small-N high-D에서 사실 classical low-rank + Gaussian kernel로 대부분 환원 가능
- Medical AI에서 quantum advantage는 아직 한 번도 입증된 적 없음
- 3억을 쓴다면 "왜 30억 과제에 필요한가"의 답이 cosmetic하면 안 됨

**솔직한 framing:**
> "임상 성능 개선 주장이 아니라, NISQ 시대 의료 AI 대비 future-proofing pilot benchmark. Classical baseline 대비 우위가 없으면 negative result도 paper로 보고."

이렇게 가야 안전. 지금 문서는 "양자영감이 small-N에서 유리한 inductive bias"라고 약한 주장을 하는데, 이게 reviewer 2한테 잡힌다.

---

## 약점 5 (메타). 데이터 자산이 self-extinguishing임

ChatGPT가 한 번도 안 짚은 큰 약점.

- 한국 PTC overdiagnosis는 의료계가 *고치고 있는* 문제
- AS가 표준이 되면 future surgery cohort 비율이 점점 줄어듦
- 즉 우리가 만든 counterfactual model의 *외부 generalization*은 시간이 지날수록 약해짐
- 모델이 "한국이 과잉진단 시대의 잔재 데이터"에 의존

**Reviewer 공격:**
> "10년 후 AS가 표준이 되면 이 모델의 임상적 가치는 어떻게 됩니까? Self-defeating asset 위에 30억을 투자합니까?"

**답변 framing:**
- (a) 현재 누적된 quasi-experimental cohort에서 *재현 불가능한* causal evidence를 추출
- (b) 이 인과 framework를 다른 overtreated cancer (low-risk prostate, DCIS, indolent lymphoma 등)로 transfer

이걸 future direction에 명시해야 함.

---

## 다음 LLM에 던질 prompt

아래를 그대로 ChatGPT/Claude에 던질 것.

```text
나는 삼성미래기술육성사업 30억 Technology 트랙 제안서를 준비 중이다.
주제: Causal K-Thyro Foundation
- 공개 고재사용 바이오데이터 (TCGA, GEO, CPTAC, DepMap, scRNA/spatial)로
  PTC disease-state foundation model 사전학습
- 한국 MAeSTro/MASTER/KoMPASS active surveillance vs surgery 코호트로
  target trial emulation 기반 counterfactual patient journey 추정
- Yu 교수 5,000례 BABA robotic thyroidectomy + 62명 CTC/EMT cohort로
  biological/surgical readout validation
- Quantum-inspired tensor/kernel은 small-N high-D representation
  benchmark sidecar (3억, claim 제한)

심사위원이 공격할 5가지 약점이 이미 식별됐다:

(1) Target trial emulation × small-N (PTC progression event ~5%):
    causal forest/DR learner의 분산이 폭발할 위험. Effective sample
    size와 unmeasured confounding sensitivity를 어떻게 답할 것인가.

(2) "Foundation model"이라는 단어가 TCGA-THCA n=500 단일 질병
    설계와 안 맞음. Pan-cancer pretraining + thyroid downstream
    구조로 다시 짜야 하는가, 아니면 "domain foundation"으로 단어
    바꿔야 하는가.

(3) AS vs Surgery에서 monitoring schedule과 treatment가 confounded
    (detection/surveillance bias). Target trial protocol에서 monitoring
    을 어떻게 emulate해야 reviewer가 안심하는가.

(4) Quantum-inspired sidecar의 정당화. "임상 성능 개선" 주장은 위험.
    NISQ future-proofing pilot으로 framing하는 게 정직한가, 아니면
    아예 빼는 게 나은가.

(5) 한국 PTC overdiagnosis가 self-extinguishing data structure임 —
    AS가 표준이 되면 미래 cohort에서 재현 불가. 이 모델의 long-term
    임상가치를 어떻게 정당화하는가.

각 약점에 대해:
- Reviewer 2가 던질 질문을 1-2문장으로
- 그에 답할 가장 강한 1문단을 제안서 본문 톤으로
- 그 답을 뒷받침할 reference 또는 method (target trial emulation,
  E-value, tipping point analysis, transfer learning evidence 등)
- 만약 답이 약하면, WP 구조를 어떻게 수정해야 하는지

5개를 다 풀어달라. 각 약점당 200-300단어. 멋진 말 말고
reviewer가 인정할 수 있는 정통 통계/causal inference 언어로.

그리고 마지막에: 이 5개 약점을 다 반영한 수정된 제안서 첫
페이지(executive summary, 250단어)를 써달라.
```

---

## 솔직한 자신감 평가

| 시나리오 | 선정 가능성 |
|---|---|
| 현재 ChatGPT 버전 그대로 제출 | 45–50% |
| 위 5개 약점 다 반영 + 표현 톤다운 | 60–65% |
| + 데이터 access 4항목 (MAeSTro/MASTER/KoMPASS/BABA/CTC) 교수님 확답 | 65–70% |
| + 외부검증 기관 1곳 사전 합의 | 70% 가능 |

70% 이상은 거짓말. 삼성육성과제는 심사위원 조합 + 경쟁과제 + PI 신뢰도에 좌우됨.

---

## 다음 행동 (우선순위)

1. 위 prompt를 ChatGPT/Claude에 던져 5개 약점 답안 받기
2. 답안을 제안서 본문/Executive Summary에 박기
3. 교수님께 1-page 보여주고 데이터 access 4항목 확인
4. 답이 "된다"면 이 방향 확정
5. 답이 애매하면 → public-data pretrained + MAeSTro-only validation 축소판으로 재설계

---

## 참고: 이 방향이 옳은 이유 (강점 재확인)

- 질문이 다름: 재발 예측이 아니라 counterfactual treatment decision
- N 경쟁 회피: 중국 대규모 thyroid AI와 정면승부 안 함
- Yu 자산 진짜: MAeSTro/MASTER/KoMPASS + 5,000례 BABA + 62명 CTC/EMT
- 양자 과장 없음: clinical claim 안 하고 representation benchmark sidecar
- 삼성 트렌드 fit: 인과 파운데이션 (2025 하반기 김광호 과제) + 디지털헬스 + Advanced AI 교차점
