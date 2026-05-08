"""
Supp #1 — Bootstrap 95% CI for pooled AUC of ViT-L CLAM predictions.

Reads: phase2 clam_per_slide_predictions.tsv
Saves: analysis_supp/bootstrap_auc.{json,txt}
"""
from __future__ import annotations
import json
import time
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

PRED = Path("project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_clam/clam_per_slide_predictions.tsv")
OUT = Path("project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/analysis_supp")

N_BOOT = 10000
RNG_SEED = 42


def main() -> None:
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(PRED, sep="\t")
    y = df["label"].to_numpy(dtype=int)
    p = df["prob_DM1"].to_numpy(dtype=float)
    n = len(y)
    point_auc = float(roc_auc_score(y, p))

    rng = np.random.default_rng(RNG_SEED)
    idx_all = np.arange(n)
    aucs = np.empty(N_BOOT, dtype=float)
    valid = 0
    for b in range(N_BOOT):
        idx = rng.choice(idx_all, size=n, replace=True)
        yb = y[idx]
        pb = p[idx]
        if len(np.unique(yb)) < 2:
            aucs[b] = np.nan
            continue
        aucs[b] = roc_auc_score(yb, pb)
        valid += 1
    aucs_valid = aucs[~np.isnan(aucs)]
    mean = float(np.mean(aucs_valid))
    std = float(np.std(aucs_valid, ddof=1))
    ci_lo = float(np.percentile(aucs_valid, 2.5))
    ci_hi = float(np.percentile(aucs_valid, 97.5))
    median = float(np.median(aucs_valid))

    res = {
        "n_slides": int(n),
        "n_pos_DM1": int((y == 1).sum()),
        "n_neg_DM2": int((y == 0).sum()),
        "point_auc": point_auc,
        "n_bootstrap": int(N_BOOT),
        "n_valid": int(valid),
        "bootstrap_mean": mean,
        "bootstrap_median": median,
        "bootstrap_std": std,
        "ci_lo_2_5": ci_lo,
        "ci_hi_97_5": ci_hi,
        "rng_seed": RNG_SEED,
        "runtime_sec": float(time.time() - t0),
    }
    (OUT / "bootstrap_auc.json").write_text(json.dumps(res, indent=2))
    log = [
        "# Bootstrap 95% CI for pooled AUC (ViT-L CLAM, DM1 vs DM2)",
        f"N slides             : {n} ({(y==1).sum()} DM1 / {(y==0).sum()} DM2)",
        f"Point AUC            : {point_auc:.4f}",
        f"Bootstrap N          : {N_BOOT} (valid: {valid})",
        f"Bootstrap mean       : {mean:.4f}",
        f"Bootstrap median     : {median:.4f}",
        f"Bootstrap std        : {std:.4f}",
        f"95% CI               : [{ci_lo:.4f}, {ci_hi:.4f}]",
        f"Runtime              : {res['runtime_sec']:.1f}s",
    ]
    (OUT / "bootstrap_auc.txt").write_text("\n".join(log))
    print("\n".join(log))


if __name__ == "__main__":
    main()
