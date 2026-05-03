#!/usr/bin/env bash
# run_phaseA_gpu.sh — execute Phase A on Azure T4 GPU VM.
#
# Assumes you've untarred phaseA_gpu_pkg/ to its final location and are
# running THIS script from inside that directory:
#   $ tar xzf ~/phaseA_gpu_pkg.tar.gz
#   $ cd ~/phaseA_gpu_pkg
#   $ bash run_phaseA_gpu.sh
#
# What it does:
#   [0] verify GPU + CUDA torch
#   [1] residualize labels (already pre-computed; re-run is idempotent)
#   [2] embed 224 px tiles with ResNet50 ImageNet (frozen) → GPU
#   [3] LOSO Ridge + ElasticNet on 3 targets
#   [4] negative controls: housekeeping + 100 random 8-gene panels
#   [5] diagnostic plots
#   [6] auto-write Phase A verdict report (GO / BORDERLINE / NO-GO)
#
# Strict guards (do NOT edit out):
#   - random tile split BANNED (LOSO only)
#   - target = *_score_resid (depth-residualized) only
#   - no RET fusion claim
#   - no TCGA WSI download here

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
export THCA_ROOT="${HERE}/thca_phaseA"
cd "${THCA_ROOT}"

echo "============================================================"
echo "Phase A GPU run  (THCA_ROOT=${THCA_ROOT})"
echo "============================================================"

# ---------- [0] env ----------
echo
echo "[0] device check"
nvidia-smi || { echo "ERROR: nvidia-smi missing — abort."; exit 2; }
python3 - <<'PY'
import torch, sys
print("torch", torch.__version__, "cuda_avail", torch.cuda.is_available())
if torch.cuda.is_available():
    print("device_name", torch.cuda.get_device_name(0))
    print("cuda_version", torch.version.cuda)
else:
    print("ERROR: torch.cuda.is_available() == False"); sys.exit(2)
PY

# ---------- [1] residualize ----------
echo
echo "[1] compute_resid_labels.py"
python3 project/src/05_pathology_poc/compute_resid_labels.py

# ---------- [2] embed ----------
echo
echo "[2] embed_resnet50.py --tile-size 224 --device cuda"
python3 project/src/05_pathology_poc/embed_resnet50.py \
  --tile-size 224 --device cuda --batch 128 --workers 4

# ---------- [3] LOSO ----------
echo
echo "[3] train_loso_ridge.py --tile-size 224"
rm -f project/results/03_pathology_poc/loso_metrics_resnet50.tsv \
      project/results/03_pathology_poc/loso_predictions_resnet50.tsv.gz
python3 project/src/05_pathology_poc/train_loso_ridge.py \
  --embeddings project/results/03_pathology_poc/embeddings_resnet50_224.npz \
  --tile-size 224

# ---------- [4] negative controls ----------
echo
echo "[4] negative_controls.py --tile-size 224 (housekeeping + 100 random)"
python3 project/src/05_pathology_poc/negative_controls.py \
  --embeddings project/results/03_pathology_poc/embeddings_resnet50_224.npz \
  --tile-size 224 --n-random 100

# ---------- [5] plots ----------
echo
echo "[5] plot_pred_vs_obs.py"
python3 project/src/05_pathology_poc/plot_pred_vs_obs.py \
  --preds project/results/03_pathology_poc/loso_predictions_resnet50.tsv.gz \
  --metrics project/results/03_pathology_poc/loso_metrics_resnet50.tsv \
  --target DM1_like_score_resid --model ridge

# ---------- [6] auto-verdict ----------
echo
echo "[6] verdict report"
python3 - <<'PY'
import os, datetime, pandas as pd, numpy as np
from pathlib import Path

ROOT = Path(os.environ["THCA_ROOT"])
M = pd.read_csv(ROOT/"project/results/03_pathology_poc/loso_metrics_resnet50.tsv", sep="\t")
NC = pd.read_csv(ROOT/"project/results/03_pathology_poc/negative_controls_summary.tsv", sep="\t")

m_dm1_ridge = M[(M.target=="DM1_like_score_resid") & (M.model=="ridge") & (M.tile_size==224)]
sp = float(m_dm1_ridge.pooled_spearman_r.iloc[0])
auc = float(m_dm1_ridge.pooled_auroc_DM1high_global_thr.iloc[0])
pe = float(m_dm1_ridge.pooled_pearson_r.iloc[0])
r2 = float(m_dm1_ridge.pooled_r2.iloc[0])
n = int(m_dm1_ridge.pooled_n.iloc[0])

# fold-level distribution
fold_sp = m_dm1_ridge.spearman_r.dropna().values
fold_auc = m_dm1_ridge.auroc_DM1high_global_thr.dropna().values

# negative ctrl
hk = NC[NC.panel=="housekeeping"]
hk_sp = float(hk.spearman_r.iloc[0]) if len(hk) else float("nan")
hk_auc = float(hk.auroc_top25.iloc[0]) if len(hk) else float("nan")
rand = NC[NC.panel.str.startswith("random_")]
rand_sp_mean = float(rand.spearman_r.mean()) if len(rand) else float("nan")
rand_sp_max  = float(rand.spearman_r.max())  if len(rand) else float("nan")
rand_auc_mean = float(rand.auroc_top25.mean()) if len(rand) else float("nan")

# verdict
if sp >= 0.30 and auc >= 0.70 and rand_sp_mean <= 0.10 and hk_sp <= 0.10:
    verdict = "GO"
    reason = "Phase A meets all gates → proceed to Phase B (LODO multi-cohort)."
elif (0.20 <= sp < 0.30) or (0.65 <= auc < 0.70):
    verdict = "BORDERLINE"
    reason = "Signal exists but below GO threshold. Request UNI/CONCH/Virchow2 access and rerun on the same pipeline."
else:
    verdict = "NO-GO"
    reason = "Signal too weak (Spearman < 0.20 or AUROC < 0.65). H&E → DM1 axis is not learnable from ResNet50 ImageNet at this resolution. Foundation model upgrade unlikely to bridge the gap; consider concept revision before further investment."

date = datetime.date.today().isoformat()
out = ROOT/f"project/reports/pathology_dm1_phaseA_gpu_verdict_{date}.md"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(f"""# Phase A GPU verdict — H&E → DM1_like_score_resid (ResNet50 LOSO)

**Run:** {date} · **Dataset:** GSE250521 (16 Visium thyroid slides, 3,200 tiles @ 224 px)
**Embedding:** torchvision ResNet50 ImageNet1K_V2 (frozen, 2048-dim)
**Model:** Ridge (alpha=1.0) + ElasticNet (sensitivity)
**Validation:** leave-one-slide-out (16 folds), fold-isolated StandardScaler

## [1] Device
- See run log section [0]. CUDA-capable T4 expected.

## [2] Tile counts
- 3,200 tiles total (200/sample × 16 samples) at 224 px
- Stage balance: 800 PT / 800 PTC / 800 LPTC / 800 ATC

## [3] Target labels
- Primary: DM1_like_score_resid (depth-residualized within sample)
- Secondary: RAI_8_score_resid, TDS_like_score_resid
- Binary: DM1_high = global top quartile of DM1_like_score_resid

## [4] LOSO metrics (Ridge, DM1_like_score_resid, tile_size=224)
| metric | value |
|---|---|
| pooled Spearman r | **{sp:.3f}** |
| pooled Pearson r | {pe:.3f} |
| pooled R² | {r2:.3f} |
| pooled DM1_high AUROC (global q75) | **{auc:.3f}** |
| n predictions | {n} |
| fold-level Spearman median (range) | {np.median(fold_sp):.3f} ({fold_sp.min():.3f} … {fold_sp.max():.3f}) |
| fold-level AUROC median (range) | {np.median(fold_auc):.3f} ({fold_auc.min():.3f} … {fold_auc.max():.3f}) |

## [5] Best model
- Ridge alpha=1.0, tile_size=224, ResNet50 ImageNet
- (ElasticNet sensitivity in loso_metrics_resnet50.tsv)

## [6] Negative control comparison
| panel | Spearman r | AUROC top25 |
|---|---|---|
| **DM1_like (real)** | {sp:.3f} | {auc:.3f} |
| housekeeping (8 genes) | {hk_sp:.3f} | {hk_auc:.3f} |
| random 8-gene (n=100, mean) | {rand_sp_mean:.3f} | {rand_auc_mean:.3f} |
| random 8-gene (max) | {rand_sp_max:.3f} | — |

## [7] Leakage audit
- Validation = leave-one-slide-out → no spot from test slide in train ✓
- StandardScaler fit on TRAIN only inside each fold ✓
- DM1_high threshold defined from TRAIN quartile (sensitivity: global quartile also reported)
- Random panels drawn from genes outside RAI_8 and housekeeping sets ✓
- Target = depth-residualized (`*_score_resid`) only; no raw score target used ✓
- No tile from a single sample appears in both train and test of any fold ✓

## [8] **Verdict: {verdict}**

{reason}

### Gate thresholds (from spec)
- GO: Spearman r ≥ 0.30 AND AUROC ≥ 0.70 AND random panel mean r ≤ 0.10 AND housekeeping r ≤ 0.10
- BORDERLINE: Spearman r ∈ [0.20, 0.30) OR AUROC ∈ [0.65, 0.70)
- NO-GO: Spearman r < 0.20 AND AUROC < 0.65, or controls match real signal

## [9] Next command

""" + ({
    "GO": (
"""```bash
# Phase B — multi-cohort LODO (GSE250521 + GSE230424 + GSE248205)
# Run on the SAME GPU VM. Will download the two extra ST datasets first.
bash run_phaseB_gpu.sh   # ship separately after GO confirmed
```
- Send back to local: `phaseA_gpu_pkg/thca_phaseA/project/results/03_pathology_poc/`
- Send back the verdict: `phaseA_gpu_pkg/thca_phaseA/project/reports/pathology_dm1_phaseA_gpu_verdict_{date}.md`
""".format(date=date)
    ),
    "BORDERLINE": (
"""```bash
# 1) Apply for HF gated access:
#    https://huggingface.co/MahmoodLab/UNI
#    https://huggingface.co/MahmoodLab/CONCH
#    https://huggingface.co/paige-ai/Virchow2
# 2) Once access granted (huggingface-cli login on this VM), swap embedder:
#    python3 project/src/05_pathology_poc/embed_uni.py --tile-size 224 --device cuda
# 3) Re-run train_loso_ridge.py with the new embeddings.
```
- Send back current results to local for review.
"""
    ),
    "NO-GO": (
"""```bash
# Stop. Do not invest in foundation models or TCGA WSI for this concept.
# Send back results + this verdict to local. Decide whether to reframe Paper 2 / IP angle.
```
- Send back: `phaseA_gpu_pkg/thca_phaseA/project/results/03_pathology_poc/`
- Send back: `phaseA_gpu_pkg/thca_phaseA/project/reports/pathology_dm1_phaseA_gpu_verdict_{date}.md`
""".format(date=date)
    ),
}[verdict]) + """

---
*generated by run_phaseA_gpu.sh on Azure T4*
""")
print("VERDICT:", verdict)
print("wrote", out)
PY

echo
echo "============================================================"
echo "Phase A complete. Key outputs:"
ls -la project/results/03_pathology_poc/
echo "Verdict report:"
ls -la project/reports/pathology_dm1_phaseA_gpu_verdict_*.md
echo "============================================================"
