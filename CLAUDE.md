# THCA data analysis — project context for Claude

## Disk layout (IMPORTANT — 2026-05-07 reorg)

Root disk is small (123G); a 512G Premium SSD is mounted at `/data`. All large outputs MUST end up on `/data`, not on root.

The reorg made this automatic via **bind mounts** (initially used symlinks but symlinks broke git's view of tracked files in `project/results/` — bind mount is the right primitive). Just write to `project/results/...` or `project/data/...` like normal; storage lands on `/data` while the path looks identical to git:

```
project/results  ← bind mount /data/thca/repo_results   (xfs on /dev/sdb)
project/data     ← bind mount /data/thca/repo_data
```

Persistent across reboots via `/etc/fstab`:
```
/data/thca/repo_results /home/seungho/personal/THCA_data_analysis/project/results none bind,nofail 0 0
/data/thca/repo_data    /home/seungho/personal/THCA_data_analysis/project/data    none bind,nofail 0 0
```

Verify with `findmnt project/results` — should show source `/dev/sdb[/thca/repo_results]`. If empty, bind mount didn't apply (re-run `sudo mount -a`).

Already-offloaded heavy intermediates (do not re-create on root):

```
project/results/v17_korean/arcasHLA            →  /data/thca/_repo_offload/arcasHLA              (8.4G)
project/results/v17_korean/arcasHLA_GSE213647  →  /data/thca/_repo_offload/arcasHLA_GSE213647    (24G)
```

Other `/data/thca/` contents (pre-existing, don't disturb):
- `PRJEB11591_fastq/` (33G) — Korean K2 raw FASTQ
- `data_raw/`, `data_processed/` — separate analysis tracks
- `reference_kallisto/`, `v17_korean/`, `v17_gse241184/`, ... — cohort outputs

## Rules

1. **Never write large outputs (>100MB) directly to** `/home/seungho/...` paths that aren't symlinked. Use `project/results/...` or `project/data/...` (these are symlinks to `/data`).
2. **If you create a new top-level dir under `project/` that you expect to grow >1GB**, mkdir it under `/data/thca/repo_<name>/` first and symlink: `ln -s /data/thca/repo_foo project/foo`.
3. **`/tmp` is on root disk (123G).** Don't dump multi-GB intermediates there. For pipelines like GNU parallel / STAR, use `--tmpdir /data/thca/_tmp` (per `v17_pod_D_k2_star_tmp_failure` lesson).
4. Check `df -h /` periodically — if root climbs past 75%, find what skipped the convention with `du -h --max-depth=2 ~/personal/THCA_data_analysis/project | sort -hr | head`.
5. The `/data` mount has `nofail` in `/etc/fstab` — if it ever fails to mount, symlinks break silently. Sanity check: `ls -la project/results` should show `→ /data/thca/repo_results`.

## RunPod conventions

- API key + SSH key: `~/.runpod/config.toml` and `~/.runpod/ssh/RunPod-Key-Go`
- Active pod (2026-05-07): **`thca-spark-dm-a6000-v4`** (id `uvp9i2r9s6l85y`), RTX A6000 48GB, **$0.33/hr** (community cloud).
  - SSH: `ssh -i ~/.runpod/ssh/RunPod-Key-Go -p 20788 root@135.84.176.142` (port/IP via API if changed).
  - Replaced prior `thca-spark-dm-l40s-v3` after that pod's SSH wedged from disk-full and L40S host became unavailable for resume.
- Pods listed via: `curl -X POST https://api.runpod.io/graphql -H "Authorization: Bearer $(grep apikey ~/.runpod/config.toml | cut -d'"' -f2)" -H "Content-Type: application/json" -d '{"query":"query { myself { pods { id name desiredStatus runtime { ports { ip publicPort privatePort type } } } } }"}'`
- New pod creation via runpodctl (GraphQL `podFindAndDeployOnDemand` returns 403; runpodctl wraps it correctly): `runpodctl create pod --name X --gpuType "NVIDIA RTX A6000" --imageName runpod/pytorch:2.4.0-py3.11-cuda12.4.1-devel-ubuntu22.04 --containerDiskSize 100 --volumeSize 50 --volumePath /runpod-volume --ports "8888/http,22/tcp" --startSSH --communityCloud --env PUBLIC_KEY=$(cat ~/.runpod/ssh/RunPod-Key-Go.pub)`
- 50GB `/runpod-volume` on pods is **pod-specific persistent volume**, not a network volume — survives pod stop, but lost on terminate. Container disk (`/workspace`, `/tmp`, `/root`) is ephemeral and wipes on stop. **Cannot migrate volume to a different pod or different host.**
- Stopped pod's GPU type cannot be changed on resume; if same host has no free GPU of that type, pod is effectively stuck. Workaround: terminate + create new pod (volume loss).
- For STAR / arcasHLA / GNU parallel on pods: `--tmpdir /workspace/_tmp --compress` to avoid `/tmp` overflow (per `v17_pod_D_k2_star_tmp_failure` lesson).

## Marathon mode (active 2026-05-04 → 2026-06-13)

See memory `v17_marathon_mode_post_pillar1`. During this window, NEW analyses must be paper-blocking only. Default mode = scaffolding/infra (this CLAUDE.md is allowed). Voice-protected sections (Hook/Aim/Disc 3.1/Limitations/Cover Para 1/Q9) are author-keyboard only — Claude must not generate prose for those.
