# Pod C (AlphaFold) — Raw Recovery Report (option C1)

**Author:** Seungho Cook
**Date:** 2026-05-04 (post-closure marathon mode)
**Decision:** C1 — raw-only recovery, no Paper 4 narrative
**Status:** ❌ **NO RECOVERY POSSIBLE** — pod terminated/removed from RunPod account before monitor could SCP results.

---

## 1. Verdict

Pod C (`runpod_pod_C_AlphaFold.sh` — LocalColabFold 4 HLA-peptide complex predictions) was dispatched at **01:26:53 KST 2026-05-04** to RunPod IP `216.81.151.3:10960`. RunPod-API check at **07:52 KST 2026-05-04** (~6.4 hr later) returned only Pod D in the user's account; Pod C is no longer present.

**Consequences:**

- All 8-12 hr compute spent on Pod C is lost (not SCP'd before pod removal)
- `SUMMARY.json`, per-structure rank_001 JSONs, and per-structure PDBs unrecoverable
- No raw `structure_summary.tsv` can be built (input file gone)
- `pillar5_supplement.md` — **not generated** (per C1 directive; would have been blocked anyway by Paper 4 freeze)

**Cause (most likely):** user manually stopped/terminated Pod C between dispatch (01:27 KST) and recovery monitor start (~07:52 KST). Consistent with `.runpod_pods.json` being already updated to `"C": null` before this session began (per `SITUATION_BACKGROUND_ACTIVITY_2026_05_04.md` §2.3 — user updated pods.json mid-session). Auto-shutdown 24hr cap was not yet hit.

---

## 2. Files intentionally NOT generated (per C1 directive)

| File | Reason |
|---|---|
| `pillar5_supplement.md` | Paper 4 freeze — narrative generation blocked regardless of pod outcome |
| Graves' / GD interpretation prose | Paper 4 freeze |
| HLA-DPB1\*05:01 autoreactive epitope discussion | Paper 4 freeze |
| `pillar5_*` any output named | Paper 4 freeze |
| `structure_summary.tsv` | Pod C raw input (`SUMMARY.json`) unrecoverable — would generate empty/error |

---

## 3. Files actually generated

**None.** No SCP target exists.

`project/results/c_alphafold/` was created as empty directory by the monitor scaffold prep but contains no files.

```
$ ls -la /home/seungho/personal/THCA_data_analysis/project/results/c_alphafold/
(empty or directory does not exist)
```

---

## 4. Lock + freeze status (unchanged)

- ✅ Paper 4 (Korean GD HLA Pan-Asian backlog) **freeze maintained** — no narrative file produced
- ✅ Paper 3 ICI Track A bundle **chmod 444 freeze** unchanged
- ✅ No manuscript prose modified
- ✅ Marathon mode discipline intact

---

## 5. Compute spent (sunk cost)

- Pod C runtime: ≤8-12 hr GPU (A100 estimated)
- Cost: ~$12-15 expected per pre-dispatch estimate (`runpod_pod_C_AlphaFold.sh` line 6)
- Recoverable: **0** (pod removed before SCP)

This is a sunk cost. Marathon decision criteria are unchanged: Paper 4 stays frozen until explicit unfreeze command + 4/4 gating per `v19_paper4_GD_backlog` memory.

---

## 6. Future re-entry conditions

Re-running Pod C (HLA-peptide AlphaFold structures) requires **all** of:

1. Explicit user unfreeze of Paper 4 backlog
2. 4/4 gating per `v19_paper4_GD_backlog` (Paper 1 bioRxiv submitted + Paper 2 A/B/C complete + Paper 3 design bundle frozen ✓ + explicit "Paper 4 시작" command)
3. Re-dispatch script may need refresh (RunPod GPU type availability changes)
4. ColabFold version pinning if reproducibility is needed (current script uses HEAD of `install_colabbatch_linux.sh`)

Until all 4 hold: **do not retry**.

---

## 7. Cross-references

- Decision: `POD_CD_RUNNING_STATE_2026_05_04.md` §3.3 (option C1 chosen by user 2026-05-04)
- Audit context: `project/reports/2026_05_04_paper1_dm1_full_conflict_audit.md` §8.2
- Original Pod C dispatch context: `SITUATION_BACKGROUND_ACTIVITY_2026_05_04.md` §2.3
- Pod C run script: `project/notebooks_or_scripts/runpod_pod_C_AlphaFold.sh`
- Original (UNUSED) post-pod script: `project/notebooks_or_scripts/post_pod_C_alphafold_analysis.py` — kept on disk but **not executed** per C1
- Paper 4 freeze authority: memory `v19_paper4_GD_backlog`

---

*No further Pod C action this session. Pod D recovery proceeding separately under `poll_and_recover_pod_D.sh`. Marathon manuscript work resumes whenever the user is ready.*
