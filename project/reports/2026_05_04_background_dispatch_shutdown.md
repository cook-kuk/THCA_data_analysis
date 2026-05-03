# Background dispatch shutdown + image-DM1 closure archival (2026-05-04)

**Author:** Seungho Cook
**Date:** 2026-05-04 (commit `e4bb222` 직후, marathon scaffolding/infra)
**Decision (user):** A안 — kill-dispatch + pods-stop-all + commit-closure-battery-md only
**Scope:** stop background dispatch on host; archive H&E → DM1 NO-GO closure decision; defer everything else.

---

## 1. Killed processes

### 1.1 `/tmp/runpod_dispatch_v2.sh` (PID 185803) — already exited naturally before kill

When this session began the dispatch loop was running (etime 03:11). At the moment user issued the kill decision, the loop had already finished its dispatch phase:

```
[01:26:05] Starting dispatch loop
[01:26:40] C READY at 216.81.151.3:10960 (uptime 25s) — dispatching
[01:26:45] C SSH OK — uploading + executing
            Started: 164
            Auto-shutdown: 24hr forced
[01:26:53] C DISPATCHED
[01:27:04] D READY at 157.157.221.29:21123 (uptime 25s) — dispatching
[01:27:14] D SSH OK — uploading + executing
            Started: 145
            Auto-shutdown: 24hr forced
[01:27:27] D DISPATCHED
[01:27:38] iter=6, dispatched=2/3 (D=1 C=1 B=0)
```

→ The loop completed dispatch to Pod C and Pod D, then exited (Pod B is `null` so 2/3 was the achievable max). Process 185803 was **no longer alive** by the time we checked.

### 1.2 `tail -F /tmp/dispatch.log` (PID 190860) — killed

```
$ kill 190860
$ ps -eo pid,etime,cmd | grep -E 'runpod_dispatch|dispatch.log' | grep -v grep
(empty — no dispatch processes remaining ✓)
```

### 1.3 No cron / systemd schedule

`crontab -l` returns no cron. The dispatcher was a one-shot foreground bash job. **No host-side respawn risk.**

---

## 2. Dispatch log tail (read-only archive of last 30 lines)

```
[01:26:05] Starting dispatch loop
[01:26:40] C READY at 216.81.151.3:10960 (uptime 25s) — dispatching
[01:26:45] C SSH OK — uploading + executing
Started: 164
Auto-shutdown: 24hr forced
[01:26:53] C DISPATCHED
[01:27:04] D READY at 157.157.221.29:21123 (uptime 25s) — dispatching
[01:27:14] D SSH OK — uploading + executing
Started: 145
Auto-shutdown: 24hr forced
[01:27:27] D DISPATCHED
[01:27:38] iter=6, dispatched=2/3 (D=1 C=1 B=0)
```

Snapshot file: `/tmp/dispatch.log` (preserved on host filesystem, not committed).

---

## 3. ⚠ Pod C / Pod D — REQUIRES USER ACTION

**Claude does NOT have the RunPod API token and was instructed not to use it.** The dispatcher already pushed work to both pods before exiting. As of this report, Pod C and Pod D are running compute on RunPod with 24-hour auto-shutdown forced.

| Pod | Role | RunPod ID (current) | SSH | Remote PID | Auto-shutdown |
|---|---|---|---|---|---|
| **C** | AlphaFold | `158ic3wrtf2j8l` | 216.81.151.3:10960 | 164 | 24 hr forced |
| **D** | K2 STAR | `8kdsltl8s2dqbo` | 157.157.221.29:21123 | 145 | 24 hr forced |
| B | (WSI) | `null` | — | — | — |

### 3.1 User action required (per Q2 = pods-stop-all)

To honor the marathon-mode "no new analysis / no compute spend" decision, **the user must stop Pod C and Pod D from the RunPod console or API**, e.g.:

- **RunPod web console:** https://www.runpod.io/console/pods → select pod by ID → Stop
- **CLI (if installed):** `runpodctl stop pod 158ic3wrtf2j8l && runpodctl stop pod 8kdsltl8s2dqbo`
- **GraphQL via curl** (using `~/.runpod/config` token):
  ```bash
  RUNPOD_KEY=$(grep RUNPOD_API_KEY ~/.runpod/config | cut -d= -f2)
  for ID in 158ic3wrtf2j8l 8kdsltl8s2dqbo; do
    curl -s -X POST https://api.runpod.io/graphql \
      -H "Authorization: Bearer $RUNPOD_KEY" \
      -H "Content-Type: application/json" \
      -d "{\"query\":\"mutation { podStop(input: {podId: \\\"$ID\\\"}) { id desiredStatus } }\"}"
  done
  ```

If the user does NOT stop them, both pods will run until either (a) their dispatched scripts finish, or (b) the 24-hour auto-shutdown fires.

### 3.2 `.runpod_pods.json` left untouched

Per Q2 instruction "do not modify `.runpod_pods.json`". File still gitignored (`.gitignore` rule from commit `e4bb222`). Current contents:

```json
{ "D": "8kdsltl8s2dqbo", "C": "158ic3wrtf2j8l", "B": null }
```

---

## 4. Committed closure files (this session)

Per Q3 = `commit-closure-battery-md only`, only the **decision artifact** is committed. The intent is to permanently record the H&E → DM1 NO-GO verdict so the question is closed for marathon and future sessions.

| Path | Bytes | Purpose |
|---|---|---|
| `CLOSURE_BATTERY_2026_05_04.md` (root) | 7,950 | Final 1-page decision summary; verdict = D (full 폐기) |
| `project/reports/pathology_dm1_closure_battery_2026_05_04.md` | 13,927 | Full forensic battery report (5 failure-mode dissection, 80-min CPU run, 7-row decision matrix) |
| `project/results/03_pathology_poc/negative_controls_raw_summary.tsv` | 6,511 | Raw-mode negative control panel (housekeeping + 50 random gene panels), supporting the NO-GO verdict |

**Already tracked (no re-add needed; unchanged since commit `25d693c`):**
- `project/results/03_pathology_poc/closure_battery_metrics.tsv` (4,095 bytes, 00:51)
- `project/results/03_pathology_poc/closure_battery_summary.png` (129,730 bytes, 00:51)

---

## 5. Intentionally NOT committed

| Path | Reason |
|---|---|
| `project_external_st/src/12_gpu/` (8 files) | New GPU pipeline scaffold dispatched to Pod C/D. Not paper-blocking under marathon scope. Defer commit until user explicitly accepts. |
| `project/notebooks_or_scripts/post_pod_B_wsi_analysis.py` | WSI post-processing for null pod B. H&E-DM1 retry adjacent — defer. |
| `project/notebooks_or_scripts/post_pod_C_alphafold_analysis.py` | AlphaFold post-processing. Pod C status pending user decision. |
| `project/notebooks_or_scripts/post_pod_D_meta3cohort.py` | K2 STAR meta-3-cohort post-processing. Pod D status pending user decision. |
| `project_external_st/results/extra/g1_external_tile_metadata.tsv.gz` | g1_embed.py output. Pipeline output not in marathon scope. |
| `project/results/03_pathology_poc/tile_metadata.tsv.gz` (modified) | Re-recorded by background closure battery; original committed in `25d693c`. Leave unstaged. |
| `project/results/03_pathology_poc/tile_metadata_resid.tsv.gz` (modified) | Same. |
| `project/results/03_pathology_poc/tile_metadata_size672.tsv.gz` | New tile size variant from closure battery. Defer. |
| `project/.runpod_pods.json` | Per Q2 = do not modify. Already gitignored via commit `e4bb222`. |
| Embeddings / tiles / raw data / WSI / FASTQ / BAM | Per absolute prohibition list. |

---

## 6. Marathon safety attestation

| Constraint | Status |
|---|---|
| Voice-protected sections (Hook / Aim / Disc §3.1 / Limitations / Cover ¶1 / Q9) | ✓ untouched this session |
| Paper 3 design bundle (chmod 444) | ✓ untouched |
| Paper 3 Track B start | ✓ NOT started |
| Paper 4 (Korean GD HLA backlog) | ✓ untouched |
| New analysis launched by Claude | ✓ none |
| New data download by Claude | ✓ none |
| H&E-DM1 closure retry by Claude | ✓ NOT attempted (closure-of-closure was user-initiated and is now archived as NO-GO) |
| TCGA WSI download by Claude | ✓ NOT attempted; Pod B (WSI) currently `null` |
| RunPod API used by Claude | ✓ NOT used (Bearer token not held; user retains pod stop authority) |
| `.runpod_pods.json` modified by Claude | ✓ NOT modified |
| `12_gpu/` pipeline committed by Claude | ✓ NOT committed (deferred) |
| Killed processes | ✓ tail -F (PID 190860). Dispatcher (185803) had already exited on its own. |

---

## 7. Open items left for next decision turn

1. **Pod C / Pod D stop confirmation.** User to verify via RunPod console or API.
2. **`12_gpu/` pipeline scripts** — keep / delete / commit decision pending. They embody a new external-ST GPU direction not yet aligned with marathon.
3. **`post_pod_B/C/D_*.py` post-processing scripts** — keep / delete / commit pending Pod outcome.
4. **Tile metadata regeneration files** (`tile_metadata*.tsv.gz`) — pre-existing committed copies remain authoritative; re-records discardable.
5. **Resume marathon work:**
   - `voice-hook` — Hook ¶1 user keyboard (first voice-protected slot)
   - `low-batch` — L1–L9 cosmetic batch (claude scaffolding)

---

## 8. Suggested user verification steps

```bash
# 1. Confirm dispatcher fully gone
ps -eo pid,etime,cmd | grep -E 'runpod_dispatch|dispatch.log' | grep -v grep
# expected: empty

# 2. Stop pods (manual — Claude did not run this)
RUNPOD_KEY=$(grep RUNPOD_API_KEY ~/.runpod/config | cut -d= -f2)
for ID in 158ic3wrtf2j8l 8kdsltl8s2dqbo; do
  curl -s -X POST https://api.runpod.io/graphql \
    -H "Authorization: Bearer $RUNPOD_KEY" \
    -H "Content-Type: application/json" \
    -d "{\"query\":\"mutation { podStop(input: {podId: \\\"$ID\\\"}) { id desiredStatus } }\"}"
done

# 3. Verify pod stop result
for ID in 158ic3wrtf2j8l 8kdsltl8s2dqbo; do
  curl -s -X POST https://api.runpod.io/graphql \
    -H "Authorization: Bearer $RUNPOD_KEY" \
    -H "Content-Type: application/json" \
    -d "{\"query\":\"query { pod(input: {podId: \\\"$ID\\\"}) { desiredStatus runtime { uptimeInSeconds } } }\"}"
done
# expected: desiredStatus = "EXITED" (or similar terminal state)
```

---

H&E-DM1 closure NO-GO archived. Background dispatcher shut down. Pods on RunPod side require user action.
