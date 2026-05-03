#!/usr/bin/env bash
# Phase A orchestrator. Runs full GSE250521 pathology POC end-to-end.
# Usage:  bash project/src/05_pathology_poc/run_phaseA.sh [DEVICE]
#         DEVICE = auto | cuda | cpu  (default auto)
set -euo pipefail
DEVICE="${1:-auto}"
ROOT="/home/seungho/personal/THCA_data_analysis"
cd "$ROOT"
source project/.venv/bin/activate

echo "=== [0] env ==="
python3 -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available(),
'name', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'NO_GPU')"

echo "=== [1] residualize labels ==="
python3 project/src/05_pathology_poc/compute_resid_labels.py

echo "=== [2] extract tiles 448 + 672 (224 already exists) ==="
for SZ in 448 672; do
  if [ "$(find project/data/processed/GSE250521/tiles -name "*size${SZ}.png" | wc -l)" -lt 16 ]; then
    python3 project/src/05_pathology_poc/extract_tiles.py --tile-size ${SZ} \
      --max-per-sample 200 \
      --meta-out project/results/03_pathology_poc/tile_metadata_size${SZ}.tsv.gz
  fi
done
# concat all tile-size metadata
python3 - <<'PY'
import pandas as pd, glob
fns = ["project/results/03_pathology_poc/tile_metadata.tsv.gz"] + \
      sorted(glob.glob("project/results/03_pathology_poc/tile_metadata_size*.tsv.gz"))
dfs = [pd.read_csv(f, sep="\t") for f in fns]
m = pd.concat(dfs, ignore_index=True).drop_duplicates(["spot_id","sample_id","tile_size"])
m.to_csv("project/results/03_pathology_poc/tile_metadata.tsv.gz",
         sep="\t", index=False, compression="gzip")
print("merged tile metadata:", m.shape, "sizes:", sorted(m.tile_size.unique()))
PY
# re-residualize after re-merge
python3 project/src/05_pathology_poc/compute_resid_labels.py

echo "=== [3] embed each tile size ==="
for SZ in 224 448 672; do
  python3 project/src/05_pathology_poc/embed_resnet50.py \
    --tile-size ${SZ} --device ${DEVICE} --batch 64 --workers 4
done

echo "=== [4] train LOSO ==="
rm -f project/results/03_pathology_poc/loso_metrics_resnet50.tsv \
      project/results/03_pathology_poc/loso_predictions_resnet50.tsv.gz
for SZ in 224 448 672; do
  python3 project/src/05_pathology_poc/train_loso_ridge.py \
    --embeddings project/results/03_pathology_poc/embeddings_resnet50_${SZ}.npz \
    --tile-size ${SZ} --append
done

echo "=== [5] negative controls (use BEST tile size from step 4) ==="
BEST_SZ=$(python3 - <<'PY'
import pandas as pd
m = pd.read_csv("project/results/03_pathology_poc/loso_metrics_resnet50.tsv", sep="\t")
m = m[(m.target=="DM1_like_score_resid") & (m.model=="ridge")]
best = m.groupby("tile_size")["pooled_spearman_r"].first().idxmax()
print(int(best))
PY
)
echo "best tile size: ${BEST_SZ}"
python3 project/src/05_pathology_poc/negative_controls.py \
  --embeddings project/results/03_pathology_poc/embeddings_resnet50_${BEST_SZ}.npz \
  --tile-size ${BEST_SZ} --n-random 100

echo "=== [6] plots ==="
python3 project/src/05_pathology_poc/plot_pred_vs_obs.py \
  --preds project/results/03_pathology_poc/loso_predictions_resnet50.tsv.gz \
  --metrics project/results/03_pathology_poc/loso_metrics_resnet50.tsv \
  --target DM1_like_score_resid --model ridge

echo "=== DONE ==="
ls -la project/results/03_pathology_poc/
