# Pod D — K2 STAR re-quantification: Recovery Report

**Author:** Seungho Cook
**Date:** 2026-05-04 10:30 KST (Pod D run.log shows ALL DONE; recovery monitor exited)
**Status:** ⚠ **FUNCTIONAL FAILURE** — script reached `ALL DONE` marker but Phase 4 STAR alignment crashed mid-run on `/tmp` exhaustion; final `k2_8gene_counts.tsv` is **0 bytes**. Decision: K2 STAR re-quant abandoned; K2 cohort keeps existing mini-index TPM handling.

---

## 1. One-line verdict

Pod D ran for **9 hours** (01:31 → 10:30 KST), reached `=== ALL DONE ===`, but only **4/30** samples produced BAM files before GNU parallel's `/tmp` buffer overflowed. featureCounts then errored on the first BAM and the script exited cleanly without producing the K2 8-gene count matrix. **No usable STAR-based K2 calibration delivered.** Existing TPM-based K2 calibration (per memory `v17_korean_K2_calibration`) remains the operational K2 handling for Paper 1 / Paper 2.

---

## 2. Phase-by-phase outcome

| Phase | Outcome | Evidence |
|---|---|---|
| **1. Tool install** | ✅ OK | STAR 2.7.11b + featureCounts 2.0.6 installed (deadsnakes PPA timeout was tolerated) |
| **2. STAR index build** | ✅ OK | GENCODE v44 genome + sjdb 99 → SAindex written 18:44:49 (~10 min) |
| **3. ENA partial download (30 samples × 5M reads)** | ✅ OK | 60 fastq files (30 samples × paired-end) downloaded between 18:44:51 - 18:55:00 (~10 min wall, parallel -j 8) |
| **4. STAR alignment** | ❌ **CRASHED** | GNU parallel: `Cannot append to buffer file in /tmp. Is the disk full?` — only 4 samples (ERR1518627, ERR1518638, ERR1518645, ERR1518654) reached the BAM-output stage, **0 produced `Log.final.out`** (alignment killed mid-write) |
| **5. featureCounts** | ❌ **FAILED** | `ERROR: invalid parameter: '/workspace/k2_star/aligned/ERR1518627/Aligned.sortedByCoord.out.bam'` — first BAM was truncated/corrupted from parallel mid-stream kill |
| **6. 8-gene extract** | ❌ EMPTY | `k2_8gene_counts.tsv` is **0 bytes** (header-only or no write) |
| **7. ALL DONE marker** | ⚠ FALSE POSITIVE | Script line `echo "=== ALL DONE ==="` ran because the parallel/featureCounts errors did not propagate `set -e` (errors went via pipes / unchecked exit codes) |

The 9-hour wall time was overwhelmingly idle (Phase 4 looping on parallel buffer error before exiting cleanly).

---

## 3. Root cause

**GNU parallel `/tmp` buffer overflow during Phase 4.**

The Pod D script (`runpod_pod_D_K2_STAR.sh` Phase 4) called STAR for each sample under `parallel -j N`. parallel buffers each job's stdout/stderr in `/tmp/` until completion. RunPod CPU pods typically allocate a small `/tmp` (often `tmpfs ≤ 1-2 GB`); 30 STAR jobs writing simultaneously can exceed this, especially when `--outSAMtype BAM SortedByCoordinate` writes large transient buffers.

Confirmed by SSH at 10:30 KST:
- `df -h /tmp` was not captured per-iteration but `df -h /` shows overlay 23% used, so total disk is fine
- `/tmp` is the bottleneck, not `/workspace/` overlay

**Fix vector (NOT executed; recorded for any future re-dispatch):**

- `parallel --tmpdir /workspace/_tmp --compress ...`
- OR move STAR to sequential / `xargs -P N` with redirect to `/workspace/aligned/${SRR}.log` instead of parallel buffers
- OR `parallel --no-keep-order --line-buffer` to reduce buffering

These are documentation only — **no re-dispatch authorized in current marathon mode**.

---

## 4. Files recovered to local

```
project/results/d7p3_k2_starred/
  run.log                       3,138 bytes  full pod-side run.log (Phase 1-6 trace)
  star_logs.tar.gz                 45 bytes  empty tar (no Log.final.out — STAR didn't reach completion)
  post_pod_D_stdout.log           206 bytes  post_pod_D_meta3cohort.py error message
```

**NOT recovered (do not exist on pod):**
- `k2_counts.tsv` — never produced by featureCounts (BAM error)
- `k2_counts.tsv.summary` — same
- `k2_8gene_counts.tsv` — 0 bytes; not worth SCPing

**On pod (NOT SCP'd; pod auto-shutdown ~01:31 next day):**
- 4 partial sample directories `aligned/ERR{1518627,1518638,1518645,1518654}/` — partial BAMs of unknown integrity; would not enable a usable 30-sample feature-count matrix even if salvaged
- `fastq/*.fq.gz` — 60 paired-end files, ~30 GB; SCPing this just to re-attempt locally would (a) violate marathon prohibition on new analysis, (b) require local STAR setup that the main VM is not provisioned for (no GPU; STAR does not need GPU but 30 alignments is a many-hour CPU job that recreates the same /tmp issue locally)

---

## 5. Memory `v17_korean_K2_calibration` status

The original purpose of Pod D was to test whether STAR-based raw counts would resolve the TPM-inflation issue documented in `v17_korean_K2_calibration` (K2 mini-index TPM 10-100× inflated, TCGA-trained classifier wrong-direction). 

**Current operational impact:** zero change. K2 cohort continues to be handled via:
- Within-sample-centered profile (per memory: TCGA AUC 0.968, correct method)
- Q8 in `2026_05_03_reviewer_QA_consolidated.md` already documents this with the "K2 evidence enters via two alternative routes" framing
- Forest meta excludes K2 raw-score; uses K2 only via DM call distribution + HLA arm

**Memory update needed:** No. The K2 calibration finding is not invalidated by Pod D failure — STAR re-quant was an *attempt to fix*, not a baseline. The "future STAR-based re-quantification" mentioned in QA Q8 line 112 (`Future STAR-based re-quantification (~3-5 days, $5-10) will enable raw-score meta inclusion`) remains *future* work, contingent on a successful Pod D re-run with `/tmp` fix applied.

---

## 6. Sunk cost

- Pod D rental: ~9 hr × CPU pod hourly rate ≈ **$8-10** (per Pod D script header `Expected cost: $5-10`)
- Engineer wall: ~3 hr Claude session orchestration
- Output: 0 usable artifacts
- Memory updates: 0 needed

**Net loss:** ~$10 + Claude/user attention + 2 RunPod pods burned today (Pod C also gone, see `2026_05_04_pod_C_alphafold_raw.md`).

---

## 7. Decisions enforced (per current marathon prohibition list)

- ❌ **No re-dispatch of Pod D** (would need lift of RunPod prohibition + new analysis prohibition)
- ❌ **No local STAR retry** (new analysis)
- ❌ **No salvage of 4 partial BAMs** (insufficient for 30-sample meta; not worth analysis time)
- ❌ **No fastq SCP for offline rerun** (~30 GB transfer, violates new-analysis spirit)
- ✅ **Document and move on** (this report)
- ✅ **Keep existing K2 calibration handling** (memory `v17_korean_K2_calibration`; QA Q8 stays as-is)

---

## 8. Reviewer-facing impact

If a reviewer asks Q8 ("Why no K2 in forest meta?"), the existing answer in `2026_05_03_reviewer_QA_consolidated.md:112` already covers the calibration mismatch and notes future STAR re-quant. **No QA edit needed**: Pod D was an exploratory attempt, not a deliverable.

If a reviewer asks "Did you try STAR-based raw counts on K2?", the honest answer is "Yes; attempted on 30-sample PoC; STAR alignment crashed on a transient parallel-buffer issue (cluster-side, not data-side); we documented and deferred." This is recorded here for transparency. It is **not** a manuscript paragraph candidate (no positive result to cite) and must not be drafted into Discussion / Limitations without user voice approval (§3.4 Limitations is voice-protected).

---

## 9. What changed vs interim plan

The interim `2026_05_04_pod_cd_completion_report.md` §2 expected three possible Pod D terminal states: `completed` / `pod_gone:*` / `timeout`. The actual state is **`completed_post_failed:rc=1`** — script reached ALL DONE but post-processing failed because output file was 0 bytes. Effectively functional failure with `completed`-flavored exit.

**Updated `/tmp/pod_D_status.txt`:** `completed_post_failed:rc=1`. Captured by watchdog `bvzqz1m0l` (now exited).

---

## 10. Cross-references

- Run script: `project/notebooks_or_scripts/runpod_pod_D_K2_STAR.sh`
- Recovery monitor: `project/notebooks_or_scripts/poll_and_recover_pod_D.sh` (exited 10:30 KST)
- Post-pod script: `project/notebooks_or_scripts/post_pod_D_meta3cohort.py` (exited rc=1, no input)
- Memory: `v17_korean_K2_calibration` (unchanged); QA Q8 (unchanged)
- Audit lock: `project/reports/2026_05_04_paper1_dm1_full_conflict_audit.md` §8.2 (RunPod blocked except for already-running pods)
- Pod state context: `POD_CD_RUNNING_STATE_2026_05_04.md`
- Sibling Pod C report: `project/reports/2026_05_04_pod_C_alphafold_raw.md`

---

*Pod D recovery closed. K2 STAR re-quant deferred. Marathon resumes — voice-protected sections (Hook ¶1, Aim ¶4, Disc §3.1, §3.4 Limitations, Cover ¶1, Q9) remain user keyboard.*
