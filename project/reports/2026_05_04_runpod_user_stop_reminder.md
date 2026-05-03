# RunPod Pod C/D — user stop reminder

**Date:** 2026-05-04 · **Owner:** Seungho Cook
**Authored by:** Claude (Opus 4.7), marathon re-entry final cleanup session.

---

## 1. Boundary statement

**Claude does NOT use the RunPod API in this session or any subsequent marathon session.** No `runpod.create_pod`, `runpod.terminate_pod`, `runpod.stop_pod`, no GraphQL pod control mutations, no `curl` against `https://api.runpod.io/graphql`. All pod start / stop actions are user-side only.

Reason: marathon-mode scope (5/4 → 6/13) restricts compute spend to paper-blocking work. Image-DM1 / G1-G2-G3 / H&E pathology angle is closed under the closure battery NO-GO (`2026_05_04_image_dm1_final_nogo_decision.md`). Continuing pod compute on related work would violate that decision.

## 2. User action required — stop pods

Two RunPod pods are still running (started ~01:26 today). User should manually stop them via the RunPod console or `runpod` CLI under the user's account.

| pod | id | role | endpoint observed | self-shutdown |
|---|---|---|---|---|
| **Pod C** | `158ic3wrtf2j8l` | AlphaFold | 216.81.151.3:10960 | 24-hour forced auto-shutdown set |
| **Pod D** | `8kdsltl8s2dqbo` | K2 STAR / meta-3-cohort | 157.157.221.29:21123 | 24-hour forced auto-shutdown set |

### Recommended user action — stop both immediately

24-hour auto-shutdown is a safety net, not a guarantee. **Immediate stop** prevents:
- needless compute spend during marathon writing window
- accidental resumption of work in the deprecated image / pod-tracked angle
- runtime-state leftovers (`project/.runpod_ssh.json`, dispatch logs) drifting further

Stop options for the user (Claude does not run these):

```bash
# Option 1: RunPod web console
# https://www.runpod.io/console/pods → select Pod C and Pod D → Stop

# Option 2: runpod CLI (user environment with API key configured)
runpod stop pod 158ic3wrtf2j8l
runpod stop pod 8kdsltl8s2dqbo

# Option 3: GraphQL via curl (user runs)
curl -X POST "https://api.runpod.io/graphql?api_key=$RUNPOD_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{"query":"mutation { podStop(input: {podId: \"158ic3wrtf2j8l\"}) { id desiredStatus } }"}'
curl -X POST "https://api.runpod.io/graphql?api_key=$RUNPOD_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{"query":"mutation { podStop(input: {podId: \"8kdsltl8s2dqbo\"}) { id desiredStatus } }"}'
```

After stop, optionally pull final outputs (if Pod C/D ran to completion already):

```bash
# Pod C (AlphaFold) outputs — if any structures of interest were predicted:
scp azureuser@216.81.151.3 -P 10960:/workspace/<output_path> ./pod_C_results/
# Pod D (K2 STAR meta-3-cohort) outputs:
scp azureuser@157.157.221.29 -P 21123:/workspace/<output_path> ./pod_D_results/
```

(Pod B was the image-DM1 pathology pod and remained `null` per dispatch log — never started, no cleanup needed there.)

## 3. After stop — what to delete locally

Once user confirms pods are stopped:

| path | action | reason |
|---|---|---|
| `project/.runpod_ssh.json` | optional delete (already gitignored) | runtime state, no longer needed once pods stopped |
| `project/.runpod_pods.json` if present | optional delete (already gitignored) | runtime state |
| `/tmp/dispatch.log`, `/tmp/runpod_dispatch_v2.sh` | optional delete | dispatch loop already exited naturally |

Claude does not delete any of these unless explicitly instructed.

## 4. Marathon compute boundary going forward

Until Paper 1 bioRxiv submission (target **2026-06-13**), additional compute beyond local CPU is permitted **only** when ALL of:

1. Block is paper-blocking for Paper 1 manuscript writing
2. Cannot be done with already-cached results on disk
3. User explicitly authorizes the spend
4. Not in any of the closed/blocked categories (image-DM1, G1/G2/G3, TCGA WSI, foundation model retry)

Default during marathon: **no new RunPod / Azure GPU spend**. Claude defaults to declining any compute request that does not satisfy all four.

## 5. Cross-references

- `2026_05_04_image_dm1_final_nogo_decision.md` (closure decision authority)
- `2026_05_04_marathon_state_after_image_dm1_nogo.md` (marathon snapshot)
- `SITUATION_BACKGROUND_ACTIVITY_2026_05_04.md` (committed in `f66a6bb` — earlier dispatch state report)
- `POD_CD_RUNNING_STATE_2026_05_04.md` (root, parallel agent's pod-state report — uncommitted; user reviews)
- `f66a6bb` commit (gitignored `.runpod_ssh.json`)

---

*Pods left running by Claude observation only. User owns the stop action under marathon-mode boundary.*
