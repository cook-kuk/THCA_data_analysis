# "Manual Workflow에서 Agentic Research로" — 22-slide 발표 outline

_60분 대상: 개발자 / 연구자_
_본인 case study: v3 → v17p35 6개월 sprint_

PPTX 직접 build는 PptxGenJS(JS)나 python-pptx 필요 — 이 outline은
slide 단위 메시지만 정리. 본인이 Slides / Keynote에서 30분 안에 만들 수 있게.

---

| # | Slide title | Key message | Visual |
|--:|---|---|---|
| 1 | 표지: "Manual Workflow에서 Agentic Research로" | 6개월 무의식 사용한 5 pattern을 framework화 | 본인 사진 + 발표일 |
| 2 | 발표자 소개 + sprint 통계 | Seungho Cook · ~80 phase · v3→v17p35 · npj submission | 통계 카드 4개 |
| 3 | Agent vs Agentic 정의 | Anthropic 공식: workflow = 코드 결정, agent = LLM 결정 | 2-column 비교표 |
| 4 | 5-rung spectrum | Workflow → Tool-augmented → Single → Multi → Autonomous | 사다리 그림 |
| 5 | 왜 이 구분 중요 | Cost / latency / 신뢰성 tradeoff가 rung별로 다름 | tradeoff 표 |
| 6 | 6 components of an agent | Goal · Planner · Executor · Memory · Tools · Verifier | hexagon 그림 |
| 7 | 본인 case study 흐름 | Phase 1 (discovery) → 2 (validation) → 3 (extension) → B (paper) → C (submission) | timeline |
| 8 | 6 component manual mapping | 6개월간 본인이 manual로 한 6 역할 표 | 매핑표 |
| 9 | 자동화 가능 ROI | Planner + Verifier 자동화가 가장 큰 ROI | 막대 차트 |
| 10 | Pattern 1 — Fix-Amplify-Synthesis | 한 cohort 완벽 → 더 많은 cohort 확장 → 종합 → 다시 | 사이클 도식 |
| 11 | Pattern 1 — 코드 일반화 | `FixAmplifyAgent` 30줄 | 코드 스니펫 |
| 12 | Pattern 2 — Tiered Parallel DAG | v17p35 Phase B: 13 task / 4 tier / 9h → 2.5h | DAG 그림 |
| 13 | Pattern 3 — Dual Dump | raw payload + compressed markdown 동시 작성 | 2-pane 예시 |
| 14 | Pattern 4 — External Reviewer Loop | 외부 LLM critique → defense → next goal | flow 그림 |
| 15 | Pattern 5 — Honest Failure Catalog | 11 failures logged · 7 → reviewer defense 활용 | 표 |
| 16 | Reference Architecture | 5 pattern + 6 component → 합쳐서 한 diagram | 큰 그림 |
| 17 | Tech stack | Claude SDK · MCP · asyncio · aiohttp · pdfplumber | 로고 grid |
| 18 | Anti-patterns 5개 | Workflow에 agent쓰기 / verifier 없이 / context bloat / hidden side effect / no kill switch | 빨간 ✗ 5개 |
| 19 | 본인 case ROI 정량 | Phase B: 9h → 2.5h (3.6x), TERT recovery: v1 4 → v2 36 (9x), reviewer attack 17/17 defended | 숫자 카드 |
| 20 | Live Demo — thyroid_replay.py | 200줄로 5 pattern 모두 실행 (1.7초) | terminal screenshot |
| 21 | Q&A 예상 | "왜 LangGraph 안 씀?" / "Claude 의존?" / "production?" | 3 질문 미리 답변 |
| 22 | Resources | GitHub repo / docs / references / contact | QR + 링크 |

---

## Speaker notes 핵심 포인트

- **Slide 1 (표지)**: 30초 — "여러분이 6개월간 사람으로 한 일을 framework로 만들면 어떻게 보이는지"
- **Slide 3 (정의)**: 2분 — Anthropic 공식 정의 강조. workflow ≠ agent.
- **Slide 8 (case mapping)**: 4분 — 본인이 어떤 역할을 manual로 했는지 솔직하게.
- **Slide 12 (Tiered DAG)**: 4분 — 실제 wall-time 숫자 강조 (3.6x).
- **Slide 19 (ROI)**: 5분 — 숫자가 청중을 설득. TERT 4→36 강조.
- **Slide 20 (demo)**: 7분 — 실제 demo. 실패 시 fallback = 미리 캡처.
- **Q&A**: 10분 — "왜 LangGraph 안 씀?" 답: "이건 framework가 아니라 5 pattern. LangGraph 위에 돌릴 수 있음."

---

## Demo cheat sheet (slide 20)

```bash
# 1. terminal 1: 발표 PC, screen share
$ cd ~/talk_demo/v18_agentic
$ python -m agentic_research.examples.thyroid_replay --max-iter 3

# Expected output (1.7초):
==========================================================
  v18 thyroid_replay -- agentic_research demo
==========================================================

[PATTERN 1] FixAmplifyAgent
  [FIX] iteration target: ...
  -> iterations: 1, success: True

[PATTERN 2] TieredParallelOrchestrator
  -> 6/6 OK across 3 tiers

[PATTERN 3] DualDumpRecorder
  -> raw: phase_B_summary.json, compressed: phase_B_summary.md

[PATTERN 4] ReviewerLoop
  -> 3 attacks, 3 defenses

[PATTERN 5] FailureCatalog
  -> 1 entries

# 2. 두 번째 terminal에서:
$ cat /tmp/v18_thyroid_replay/dumps/phase_B_summary.md
# (compressed dump 보여주기)

$ ls /tmp/v18_thyroid_replay/memory/audit/
# (per-iteration audit 보여주기)
```

**실패 시 fallback:**
- screenshot 5장 미리 준비 (output / dump / audit / failures.jsonl / fact store)
- "demo가 실패했다는 사실 자체가 pattern 5 (failure catalog)의 예시" — 농담으로 받기

---

## Q&A 예상 + 답변

**Q1: 왜 LangGraph / LlamaIndex 안 썼어요?**
A: 이 패키지는 framework가 아니라 5 pattern. LangGraph는 graph DSL, 우리는
loop shape. 두 개 같이 쓸 수 있음. LangGraph 안에 `FixAmplifyAgent`를
node로 넣을 수 있음.

**Q2: Claude에 의존?**
A: `Tool` protocol이 abstract이라 OpenAI / Gemini / 로컬 LLM 모두 가능.
실제 demo는 Claude가 backbone이지만 패키지는 모델 무관.

**Q3: Production에 쓸 수 있나요?**
A: 0.1.0 alpha. tutorial 5에 production concerns 정리. 핵심: cost cap,
verifier always, dual dump always. external 사용자 1명 이상 받기 전엔
semver 약속 안 함.

**Q4: 다음 step?**
A: (1) v17 paper acceptance에 따라 case study extend. (2) 한국 cohort
collaboration. (3) MCP 표준 native integration. (4) 외부 reviewer LLM
2개 이상 동시 (multi-reviewer).
