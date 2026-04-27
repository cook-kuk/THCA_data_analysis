# Demo cheat sheet — slide 20

## Pre-flight (do this BEFORE the talk)

1. Confirm working directory:
   ```bash
   cd ~/talk_demo/v18_agentic
   pwd  # must end in /v18_agentic
   ```

2. Confirm package importable:
   ```bash
   python3 -c "from agentic_research import Goal, FixAmplifyAgent; print('OK')"
   ```

3. Pre-run the demo once to warm up imports:
   ```bash
   python3 -m agentic_research.examples.thyroid_replay --max-iter 2
   rm -rf /tmp/v18_thyroid_replay   # clean for the live demo
   ```

4. Capture 5 fallback screenshots in `talk/screenshots/`:
   - `01_demo_output.png` — full terminal output of the demo run
   - `02_dump_md.png` — `cat /tmp/v18_thyroid_replay/dumps/phase_B_summary.md`
   - `03_audit_dir.png` — `ls /tmp/v18_thyroid_replay/memory/audit/`
   - `04_failures.png` — `cat /tmp/v18_thyroid_replay/failures.jsonl`
   - `05_fact_store.png` — `cat /tmp/v18_thyroid_replay/memory/facts.jsonl`

---

## Live demo script (slide 20, 7 minutes)

### Beat 1 — context (30s)

> "200줄짜리 데모가 5가지 패턴을 모두 돌리는 걸 보여드리겠습니다."

### Beat 2 — run (45s)

```bash
python3 -m agentic_research.examples.thyroid_replay --max-iter 3
```

Expected output (1.7s, deterministic enough for talk):

```
==========================================================
  v18 thyroid_replay -- agentic_research demo
==========================================================

[PATTERN 1] FixAmplifyAgent
  [FIX] iteration target: Reproduce DM1/DM2 -> 8-gene RAI panel -> 4-group survival narrative
  -> iterations: 1, success: True

[PATTERN 2] TieredParallelOrchestrator
  -> 6/6 OK across 3 tiers

[PATTERN 3] DualDumpRecorder
  -> raw: phase_B_summary.json, compressed: phase_B_summary.md

[PATTERN 4] ReviewerLoop
  -> 3 attacks, 3 defenses

[PATTERN 5] FailureCatalog
  -> 1 entries

==========================================================
  Demo complete. Workspace: /tmp/v18_thyroid_replay
==========================================================
```

### Beat 3 — show the dual dump (90s)

```bash
cat /tmp/v18_thyroid_replay/dumps/phase_B_summary.md
```

> "Compressed markdown. 다음 iteration이 이걸 읽고 진행해요. Raw JSON은 옆에 같이 있고요."

```bash
ls /tmp/v18_thyroid_replay/dumps/
```

### Beat 4 — show the audit trail (90s)

```bash
ls /tmp/v18_thyroid_replay/memory/audit/
cat /tmp/v18_thyroid_replay/memory/audit/*.json | head -30
```

> "매 iteration마다 artifacts 다 기록. 어떤 phase에서 뭐가 나왔는지 1초 안에 검색."

### Beat 5 — show the failure catalog (60s)

```bash
cat /tmp/v18_thyroid_replay/failures.jsonl
```

> "실패 기록. v1 TERT recovery가 thca_tcga_pub 놓쳤다는 entry. 이걸 안 지우고 catalog로 두니까 다음 sprint가 같은 실수 안 해요."

### Beat 6 — show the fact store (45s)

```bash
cat /tmp/v18_thyroid_replay/memory/facts.jsonl
```

> "Structured facts. 'AUC 0.954, confidence HIGH'. 다음 iteration이 raw를 다시 안 읽고 fact만 recall."

### Beat 7 — wrap (30s)

> "200줄. 5 pattern. 1.7초. 6개월 manual sprint를 이렇게 generalize했어요."

---

## Failure modes + fallback

| Failure | Fallback |
|---|---|
| `python3` 명령어 없음 | `python` 시도; pre-flight에서 alias 확인 |
| 패키지 import 실패 | 미리 `pip install -e .` 확인; pre-flight 단계 |
| Demo run 실패 | screenshot 1번 보여주기, "demo 실패도 pattern 5의 좋은 예시" 농담 |
| 시간 초과 | beat 3-5만 빠르게, beat 6-7 스킵 |
| 청중 질문 끼어듦 | beat 7까지 마치고 답변 |

---

## Cleanup (발표 후)

```bash
rm -rf /tmp/v18_thyroid_replay
```
