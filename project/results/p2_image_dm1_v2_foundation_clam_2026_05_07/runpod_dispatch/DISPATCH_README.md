# RunPod dispatch — Foundation-model + CLAM Image-DM1 v2

**Date:** 2026-05-07
**Sprint root:** `project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/`

---

## Quick start (after RunPod pod up)

### 0. Pod setup recommendation
- **GPU:** 1× NVIDIA A100 40 GB (Phase 1) or 1× A100 80 GB (Phase 2)
- **Image:** `runpod/pytorch:2.4.0-py3.11-cuda12.4-devel-ubuntu22.04`
- **Storage:** 50 GB persistent + 100 GB temp
- **Cost:** $1.5-3/hr typical

### 1. Local → Pod sync
On local:
```bash
cd /home/seungho/personal/THCA_data_analysis
# tar the sprint workspace + 3.4 GB GSE250521 tiles + master labels
tar czf /tmp/p2_image_dm1_v2.tar.gz \
  project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/ \
  project/data/processed/GSE250521/tiles/ \
  project/results/01_spatial_score/all_spots_scored.tsv.gz \
  project/results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv

# upload via runpodctl or scp
scp /tmp/p2_image_dm1_v2.tar.gz root@<POD_IP>:/workspace/
```

### 2. Pod environment
```bash
ssh root@<POD_IP>
cd /workspace
tar xzf p2_image_dm1_v2.tar.gz
cd project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/
bash runpod_dispatch/env_setup.sh
```

### 3. Authenticate HuggingFace (UNI gated)
Go to https://huggingface.co/MahmoodLab/UNI and request access (instant for academic).
Then:
```bash
export HF_TOKEN=hf_xxx_your_token
```

### 4. Run phases
```bash
# Phase 1 — UNI on GSE250521 (1-2 hr)
bash runpod_dispatch/run_phase1.sh

# Check Phase 1 verdict
cat phase1_gse250521/PHASE1_REPORT.md
# IF PASS or MARGINAL → Phase 2

# Phase 2 — TCGA WSI download + CLAM (4-8 hr)
bash runpod_dispatch/run_phase2.sh

# Check verdict
cat phase2_tcga_clam/PHASE2_REPORT.md
```

### 5. Sync back to local
```bash
# On pod
tar czf /tmp/p2_results_back.tar.gz phase1_gse250521/ phase2_tcga_clam/

# On local
scp root@<POD_IP>:/tmp/p2_results_back.tar.gz /tmp/
cd /home/seungho/personal/THCA_data_analysis
tar xzf /tmp/p2_results_back.tar.gz -C project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/
```

### 6. Phase 3 — locally
```bash
.venv/bin/python project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/scripts/phase3_integration.py
```

---

## Cost / time estimates

| Phase | Time | GPU | Cost |
|-------|------|-----|------|
| 1 (UNI on 16 slides ~30k tiles) | 1-2 hr | 1× A100 | $4 |
| 2 (TCGA WSI download + CLAM) | 4-8 hr | 1× A100 | $15-25 |
| **Total** | 6-10 hr | | **$20-30** |

Hard cap: **$50**.

---

## Kill-switch logic

```
Phase 1 verdict
├── PASS (≥50% slides UNI PC1↔DM1 |ρ|>0.3) → Phase 2 OK
├── MARGINAL → Phase 2 OK with warning
└── FAIL → STOP; document negative for Paper 1 supp

Phase 2 verdict (held-out AUC)
├── PASS (>0.70) → Paper 2 image-DM1 launch unlocked
├── MARGINAL (0.60-0.70) → Paper 1 spatial supp only
└── FAIL (≤0.60) → closure NO-GO re-confirmed v2
```

---

## Troubleshooting

- **UNI download fails:** check HF_TOKEN; need to accept license at https://huggingface.co/MahmoodLab/UNI
- **CUDA OOM:** reduce `--batch_size` to 64 or 32 in phase1_uni_embed.py
- **GDC slow:** use `gdc-client` instead (`pip install gdc-client`); supports parallel
- **TCGA SVS open issue:** need `openslide` installed (`apt install libopenslide-dev` + `pip install openslide-python`)
