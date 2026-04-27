# Speaker notes — "Manual Workflow에서 Agentic Research로"

_60분 대상 발표 · 22 slide · 청중: 개발자 / 연구자_

각 slide별로 시간 / 핵심 메시지 / 대본 / transition을 정리. 시간은
타이트하게 잡았고, Q&A 10분 별도.

---

## Slide 1 — 표지 (30초)

**메시지:** "여러분이 6개월간 사람으로 한 일을 framework로 만들면 어떻게
보이는지"

**대본:**

> "안녕하세요. Seungho Cook입니다. 오늘은 갑상선암 transcriptomic 분석을 6개월간 하면서 무의식적으로 사용한 5가지 pattern을, 어떻게 reusable framework로 일반화했는지 공유하려 합니다. 결론부터 말씀드리면 — 여러분이 long-running research project에서 manual로 하는 orchestration 작업은 대부분 자동화 가능합니다. 단, 아무 framework나 가져다 쓰면 안 됩니다."

**Transition:** "먼저 용어부터 정리하고 가겠습니다."

---

## Slide 2 — 발표자 + sprint 통계 (1분)

**메시지:** 본인 신뢰성 확보 — 실제 multi-month project 경험자.

**대본:**

> "저는 독립 computational researcher이고, 지난 6개월간 TCGA 갑상선암 cohort를 가지고 npj Precision Oncology submission을 준비했습니다. v3에서 v17p35까지, 약 80개 phase-level sprint를 돌렸습니다. 그 과정에서 발견한 5가지 pattern이 오늘 발표의 모든 것입니다."

**카드 4개 (slide 표시):**
- ~80 phase
- 6 months
- v3 → v17p35
- npj submission ready

---

## Slide 3 — Agent vs Agentic 정의 (2분)

**메시지:** Anthropic 공식 정의. workflow ≠ agent.

**대본:**

> "Anthropic이 2024년에 'Building effective agents'라는 글을 publish했는데, 거기서 명확하게 정의합니다. Workflow는 코드가 control flow를 결정하고 LLM은 fixed branch 중에서 고르는 거. Agent는 LLM이 dynamic하게 자기 tool use와 control flow를 결정하는 거. 이 구분이 중요한 이유는 — workflow에 agent를 쓰면 비용 낭비, agent가 필요한 데 workflow를 쓰면 막히기 때문입니다."

**Transition:** "그럼 이걸 5단계로 펼쳐 봅시다."

---

## Slide 4 — 5-rung spectrum (1.5분)

**메시지:** Workflow → Tool-augmented → Single agent → Multi-agent → Autonomous.

**대본:**

> "Rung 1은 workflow — LangChain SequentialChain 같은 거. Rung 2는 function-calling. Rung 3은 ReAct/AutoGPT처럼 LLM이 loop을 직접 결정. Rung 4는 multi-agent — Anthropic의 multi-agent research system. Rung 5는 autonomous — 사실상 production deployment. 오늘 framework는 rung 3-4 target입니다. Rung 5는 production guardrail이 따로 필요해서 scope 밖."

---

## Slide 5 — 왜 이 구분 중요 (2분)

**메시지:** Cost / latency / 신뢰성 tradeoff가 rung별로 다름.

**대본:**

> "Rung 1은 싸고 빠르고 predictable. Rung 3은 비싸고 non-deterministic. Rung 5는 돈을 잃을 수 있어요. 핵심은 — 여러분이 control flow를 사전에 알 수 있으면 rung 1을 쓰세요. LLM에게 '뭘 할지 결정해' 하지 마세요. 그건 단순히 비싸고 비효율적입니다. Rung 3은 진짜로 사전에 결정 못 할 때만 쓰는 거예요."

---

## Slide 6 — 6 components of an agent (3분)

**메시지:** Goal · Planner · Executor · Memory · Tools · Verifier.

**대본:**

> "어떤 agent든 이 6가지 component가 있습니다. Goal은 '뭘 만들 거냐'. Planner는 'goal을 task로 쪼개는 일'. Executor는 'task를 실제로 돌리는 일'. Memory는 'iteration 간에 state를 보존하는 일'. Tools는 'executor가 호출하는 것들 — bash, fetch, LLM, MCP'. Verifier는 '결과가 충분히 좋은지 판단하는 것'. 이 6 component가 명시적이지 않은 agentic system은, 결국 안 보이게 똑같이 6 component를 가집니다 — 단, 그게 사람의 머릿속에 있어요."

**Transition:** "제 case가 정확히 그랬습니다."

---

## Slide 7 — case study 흐름 (2분)

**메시지:** Phase 1 (discovery) → 2 (validation) → 3 (extension) → B (paper) → C (submission).

**대본:**

> "v3에서 시작한 게 첫 번째 phase. 갑상선암 TCGA 513명을 unsupervised clustering으로 DM1/DM2 axis를 발견. 그 다음 phase는 external validation — GSE76039, GSE126698 등 5 cohort. 그 다음은 extension — 8-gene RAI panel + TERT integration. 그 다음은 phase B — 13 task을 4 tier로 parallel 돌려서 paper draft 만들기. 마지막은 phase C — submission package."

---

## Slide 8 — 6 component manual mapping (4분)

**메시지:** 본인이 6 역할 manual로 한 표.

**대본:**

> "솔직히 말씀드리면 — 6개월 동안 저는 6가지 역할을 모두 사람으로 했습니다."

| Component | 본인이 한 일 |
|---|---|
| Goal | "Defensible npj submission" — 머릿속에 |
| Planner | Phase 경계 손으로 정함 |
| Executor | 여러 Claude Code session, 여러 terminal |
| Memory | repo + `.claude/memory/MEMORY.md` |
| Tools | bash, lifelines, scipy, pdfplumber |
| Verifier | 본인이 결과 읽고 + 외부 LLM critique |

> "Planner와 Verifier가 가장 큰 manual cost였어요. Executor는 Claude가 했어요. Memory는 file system이 했어요. 즉 — 자동화 가능한 게 명확히 보입니다."

---

## Slide 9 — 자동화 가능 ROI (1.5분)

**메시지:** Planner + Verifier가 가장 큰 ROI.

**대본:**

> "Planner를 자동화하면 phase decomposition 결정 시간이 줄어요. Verifier를 자동화하면 매 iteration 끝마다 사람이 읽는 시간이 줄어요. 두 개가 가장 큰 leverage입니다. Memory와 Tools는 이미 충분히 자동화돼 있어요. Goal은 본인이 정해야 하는 거고."

---

## Slide 10-11 — Pattern 1 (Fix-Amplify-Synthesis) (4분)

**메시지:** 한 cohort 완벽 → 더 많은 cohort 확장 → 종합 → 다시.

**대본:**

> "첫 번째 pattern. v3-v17 sprint 내내 무의식적으로 했던 게 이 cycle입니다. Fix는 한 cohort에서 한 가지 defect를 정밀 수정. Amplify는 그 fix를 다른 cohort로 확장. Synthesis는 전부 audit하고 paper section으로 정리. 그리고 새로운 defect가 나오면 다시 fix로. 이걸 클래스로 만든 게 `FixAmplifyAgent`. 30줄짜리. (다음 slide로) 코드 보시면 — handler 3개를 받아서 sequential로 돌려요. Verifier가 pass라고 하면 break, 아니면 다음 iteration."

**[demo]: code snippet on screen]**

---

## Slide 12 — Pattern 2 (Tiered DAG) (4분)

**메시지:** v17p35 Phase B = 13 task / 4 tier / 9h → 2.5h.

**대본:**

> "두 번째. Phase B에서 paper 만들 때 13개 task를 4개 tier로 나눴어요. Tier 1은 figure 만들기 6개 — 서로 의존성 없으니 parallel. Tier 2는 multi-panel composite. Tier 3은 manuscript draft. Tier 4는 reviewer defense. Within-tier parallelism으로 wall time 9시간에서 2.5시간으로 줄였어요. 3.6배. 이걸 framework로 만든 게 `TieredParallelOrchestrator`. asyncio 기반."

---

## Slide 13 — Pattern 3 (Dual Dump) (3분)

**메시지:** Raw payload + compressed markdown 동시 작성.

**대본:**

> "세 번째. 매 phase마다 두 파일을 같이 씁니다. Raw는 full json/tsv — replay용. Compressed는 60초 안에 읽히는 markdown — 사람용. 이게 왜 중요하냐 — 다음 iteration이 raw를 읽으면 context window가 터집니다. Compressed만 읽고 진행하면 됩니다. 이 framework에서는 `DualDumpRecorder`로 한 줄에 둘 다 작성."

---

## Slide 14 — Pattern 4 (Reviewer Loop) (4분)

**메시지:** 외부 LLM critique → defense → next goal.

**대본:**

> "네 번째 — 가장 leverage 높은 pattern입니다. 매 phase 끝마다 외부 LLM (저는 Claude한테 던지고 같은 결과를 다른 LLM한테도 던졌어요)에게 critique 받아요. 그 critique에서 top 3 attack을 뽑아서 next iteration의 goal로 promote. 이걸 17번 attack에 대해 17번 defense 만들었어요. Reviewer attack 막을 능력이 paper acceptance probability 결정합니다."

---

## Slide 15 — Pattern 5 (Failure Catalog) (3분)

**메시지:** 11 failures logged · 7 → reviewer defense.

**대본:**

> "다섯 번째. 실패한 접근을 지우지 말고 catalog에 기록. 'PRE-2 ComBat under LODO returned DIAL=0.000 — this is true biology, not a leak', 'v1 TERT recovery missed thca_tcga_pub — check ALL cBio thyroid studies'. 이런 entry 11개 쌓였고, 그 중 7개가 paper의 limitations section과 reviewer defense에 그대로 들어갔어요. 정직하게 실패 기록하면 reviewer 공격 차단됩니다."

---

## Slide 16 — Reference Architecture (2분)

**메시지:** 5 pattern + 6 component → 합쳐서 한 diagram.

**대본:**

> "5 pattern이 6 component 위에서 어떻게 합치는지 한 그림으로 보면 — Goal이 Planner를 통해 Tasks로 decompose, Executor가 Tools 호출, 결과는 Memory에 저장, Verifier가 통과 판정, 통과 안 되면 Reviewer Loop이 돌고, 실패는 Failure Catalog로. 이게 architecture입니다."

---

## Slide 17 — Tech stack (1분)

**대본:**

> "Claude SDK를 backbone으로 썼지만 Tool protocol이 abstract라 OpenAI나 로컬 LLM 다 가능. asyncio가 orchestrator 코어. aiohttp는 multi-URL fetch. pdfplumber는 PDF table 추출. MCP는 production 통합용."

---

## Slide 18 — Anti-patterns 5개 (3분)

**대본:**

> "안 하시면 좋겠는 것 5개. (1) workflow 문제에 agent 쓰기 — control flow 정해져 있으면 그냥 코드 쓰세요. (2) Verifier 없이 loop 돌리기 — 그건 그냥 pipeline. (3) Context bloat — memory 무한 증가. recall(min_confidence='HIGH') 쓰세요. (4) Hidden side effect — 파일 쓰면서 log 안 하면 replay 불가능. dual dump 항상. (5) Kill switch 없음 — wall time / cost / SIGTERM 다 갖추세요."

---

## Slide 19 — 본인 case ROI 정량 (5분 — 가장 중요)

**대본:**

> "구체적인 숫자로 ROI 보여드리겠습니다."

| 메트릭 | Manual | Framework | 배수 |
|---|---:|---:|---:|
| Phase B wall time | 9시간 | 2.5시간 | 3.6x |
| TERT recovery (Sanger calls) | 4 (v1) | 36 (v2) | 9x |
| Reviewer attack defended | 3/3 (manual) | 17/17 (with loop) | 5.7x |
| Phase decomposition (manual planning) | ~2 hours | ~10 min | 12x |

> "특히 TERT 부분 — v1에서는 cBioPortal에서 2개 study만 봐서 4개 mutation밖에 못 찾았어요. 그런데 v2에서는 9 source를 systematic하게 sweep해서 36개 발견. 같은 데이터지만 framework 쓰니까 9배 결과. Survival logrank p = 4.92e-06. 이게 paper main figure로 들어갑니다."

---

## Slide 20 — Live Demo (7분 — 시간 가장 많이)

**Setup:**
- Terminal 1: `cd ~/talk_demo/v18_agentic`
- Terminal 2: `tail -f /tmp/v18_thyroid_replay/...`

**Script:**

```
$ python -m agentic_research.examples.thyroid_replay --max-iter 3
```

> "지금 이 200줄짜리 demo가 5 pattern 모두 실제로 실행합니다. 1.7초 걸려요. (output 보면서) Pattern 1 — fix-amplify cycle 1 iteration에 succeed. Pattern 2 — 6개 task 3 tier에서 다 성공. Pattern 3 — dual dump 작성. Pattern 4 — 3 attack, 3 defense. Pattern 5 — 1 failure 기록. (다른 terminal로) 그리고 이게 dump한 파일들 — raw payload + compressed markdown."

```
$ cat /tmp/v18_thyroid_replay/dumps/phase_B_summary.md
```

> "이런 식으로 매 phase의 결과가 사람이 읽는 형태로 남아요. 다음 phase는 이걸 읽고 진행."

**실패 시 fallback:**
- Screenshot 5장 미리 준비.
- 농담: "demo 실패한 것 자체가 pattern 5 (failure catalog)의 좋은 예시"

---

## Slide 21 — Q&A 예상 (이건 발표 안 함, 준비)

(slide deck Q&A appendix로 두고 넘어가기)

---

## Slide 22 — Resources (30초)

**대본:**

> "GitHub repo, 5개 tutorial, references.md에 prior work 다 있어요. 질문 받겠습니다."

---

## 발표 후 Q&A 답변 (10분)

### Q1: 왜 LangGraph / LlamaIndex 안 썼어요?

> "이 패키지는 framework가 아니라 5 pattern입니다. LangGraph는 graph DSL로 control flow 정의용이고 우리는 loop shape 정의용. 두 개 같이 쓸 수 있어요. LangGraph 안에 `FixAmplifyAgent`를 node로 넣을 수 있어요. 그리고 — solo researcher로서 LangGraph가 dependency로 가져오는 무게가 부담스러웠어요. 이건 600줄 안 됩니다."

### Q2: Claude에 의존?

> "Tool protocol이 abstract이라 OpenAI / Gemini / Llama 다 가능. demo만 Claude. 패키지 자체는 LLM agnostic. Anthropic이 references에 있는 건 Building effective agents 글 때문이지 dependency가 아닙니다."

### Q3: Production에 쓸 수 있나요?

> "0.1.0 alpha. tutorial 5에 production concerns 정리. 핵심: cost cap, verifier always, dual dump always. external 사용자 1명 이상 받기 전엔 semver 약속 안 합니다. Production에 가려면 (1) cost tracking, (2) MCP integration, (3) multi-reviewer support이 필요합니다."

### Q4: 왜 이 5개 pattern만? 다른 건 없나?

> "6개월 sprint 끝나고 회고할 때, 자주 나타난 게 이 5개였어요. 다른 pattern들 — 'speculative branch', 'cross-cohort vote', 'lazy verification' — 도 있는데, 한두 번씩만 나타났어요. 5개로 줄인 건 의도적입니다. Pattern은 익숙해질 때까지 안 쓰면 의미 없어요. 5개 정도가 한 번에 익숙해질 수 있는 양."

### Q5: 다음 step?

> "(1) v17 paper 받아들여지면 case study extend. (2) 한국 cohort collaboration 통해 cross-cohort validation. (3) MCP native integration. (4) 외부 reviewer 2개 이상 동시 (multi-reviewer). 가장 우선은 external user 받는 것 — 본인 case 외에 진짜로 쓸 만한지 검증 필요."

---

## 발표 후 todo

1. README의 GitHub URL placeholder 채우기
2. PyPI에 0.1.0 publish
3. 본인 trip report blog post 1편 작성 (이 발표의 후기)
