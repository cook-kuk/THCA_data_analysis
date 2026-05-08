#!/usr/bin/env python3
"""
ResNet50 embedding 위 multi-classifier LOSO comparison.

Phase A 1차 결과 (Ridge/ENet) 가 chance-level 인 게 ResNet50 representation 한계인지,
classifier 한계인지 분리. RandomForest / GradientBoosting / MLP / kNN 으로 같은 LOSO
target 에 대해 4 model 비교.

Inputs:
  project/results/03_pathology_poc/embeddings_resnet50_{224,448,672}.npz
  project/results/03_pathology_poc/tile_metadata_resid.tsv.gz

Outputs:
  project/results/03_pathology_poc/multi_classifier_loso_metrics.tsv
"""
from __future__ import annotations
import gzip, sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent.parent
RES = ROOT / "project/results/03_pathology_poc"
META = RES / "tile_metadata_resid.tsv.gz"
OUT = RES / "multi_classifier_loso_metrics.tsv"

from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr, pearsonr

MODELS = {
    "rf":  RandomForestRegressor(n_estimators=200, max_depth=8, n_jobs=2, random_state=0),
    "gbr": GradientBoostingRegressor(n_estimators=150, max_depth=3, random_state=0),
    "mlp": Pipeline([("sc", StandardScaler()),
                     ("m", MLPRegressor(hidden_layer_sizes=(256,64),
                                        max_iter=120, random_state=0, early_stopping=True))]),
    "knn": Pipeline([("sc", StandardScaler()),
                     ("m", KNeighborsRegressor(n_neighbors=25, weights="distance", n_jobs=2))]),
}


def main():
    if not META.exists():
        print(f"missing {META}; run compute_resid_labels.py first")
        sys.exit(0)
    meta = pd.read_csv(META, sep="\t")
    print(f"meta rows: {len(meta)}, samples: {meta.sample_id.nunique()}")

    rows = []
    for ts in [224, 448, 672]:
        emb_path = RES / f"embeddings_resnet50_{ts}.npz"
        if not emb_path.exists():
            print(f"skip size {ts} (no {emb_path.name})")
            continue
        z = np.load(emb_path, allow_pickle=True)
        emb = z["embeddings"]
        spot_ids = z["spot_ids"]
        tile_sz = int(z.get("tile_size", ts)) if "tile_size" in z.files else ts
        print(f"size {tile_sz}: emb={emb.shape}, spots={len(spot_ids)}")

        m = meta[meta.tile_size == ts].copy()
        m["row_idx"] = m.spot_id.map({s: i for i, s in enumerate(spot_ids)})
        m = m[m.row_idx.notna()].copy()
        m["row_idx"] = m.row_idx.astype(int)

        target = "DM1_like_score_resid"
        if target not in m.columns:
            print(f"target {target} missing; cols={list(m.columns)[:10]}")
            continue
        y = m[target].values
        X = emb[m.row_idx.values]
        groups = m.sample_id.values
        unique = sorted(set(groups))
        print(f"  LOSO over {len(unique)} slides ; X={X.shape} y={y.shape}")

        for mname, model in MODELS.items():
            preds = np.full_like(y, fill_value=np.nan, dtype=float)
            for g in unique:
                tr = groups != g
                te = groups == g
                if tr.sum() < 50 or te.sum() < 5:
                    continue
                try:
                    model.fit(X[tr], y[tr])
                    preds[te] = model.predict(X[te])
                except Exception as e:
                    print(f"  {mname}/{g} failed: {e}")
            keep = ~np.isnan(preds)
            if keep.sum() < 20:
                continue
            sr = float(spearmanr(y[keep], preds[keep]).statistic)
            pr = float(pearsonr(y[keep], preds[keep])[0])
            try:
                yc = (y[keep] > np.quantile(y[keep], 0.75)).astype(int)
                au = float(roc_auc_score(yc, preds[keep]))
            except Exception:
                au = float("nan")
            rows.append({"tile_size": ts, "model": mname,
                         "n": int(keep.sum()),
                         "pooled_spearman_r": sr,
                         "pooled_pearson_r": pr,
                         "pooled_auroc_top25": au})
            print(f"  {mname:>4}: r_s={sr:+.3f}  r_p={pr:+.3f}  AUROC={au:.3f}")

    if not rows:
        print("no results"); sys.exit(0)
    pd.DataFrame(rows).to_csv(OUT, sep="\t", index=False)
    print(f"\nwrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
