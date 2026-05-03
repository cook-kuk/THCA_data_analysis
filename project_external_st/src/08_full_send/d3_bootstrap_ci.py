#!/usr/bin/env python3
"""D3 — Bootstrap 95% CI on cross-validation Pearson r (대박 #1)."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

OUT = Path("project_external_st/results/extra/d3_bootstrap_ci.tsv")


def main():
    rng = np.random.default_rng(42)
    s = pd.read_csv("project_external_st/results/meta/sample_level_score_summary.tsv", sep="\t")
    pairs = [
        ("all 12 external slides", "mean_DM1_like_score_resid_epi25", "mean_THYROID_NONOVERLAP_score_resid_epi25"),
        ("all 12 epi50",           "mean_DM1_like_score_resid_epi50", "mean_THYROID_NONOVERLAP_score_resid_epi50"),
        ("all 12 all_spots",       "mean_DM1_like_score_resid_all",   "mean_THYROID_NONOVERLAP_score_resid_all"),
    ]
    rows = []
    for label, xcol, ycol in pairs:
        x = s[xcol].dropna().to_numpy(); y = s[ycol].dropna().to_numpy()
        if len(x) != len(y) or len(x) < 4: continue
        r_obs, p_obs = pearsonr(x, y); rho_obs, p_rho = spearmanr(x, y)
        # bootstrap 5000
        n = len(x); bs_r = []
        for _ in range(5000):
            i = rng.integers(0, n, n)
            try: bs_r.append(pearsonr(x[i], y[i])[0])
            except Exception: bs_r.append(np.nan)
        bs_r = np.array(bs_r); bs_r = bs_r[~np.isnan(bs_r)]
        ci_lo, ci_hi = np.percentile(bs_r, [2.5, 97.5])
        rows.append({"scope": label, "n": n, "pearson_r": r_obs, "pearson_p": p_obs,
                     "spearman_rho": rho_obs, "spearman_p": p_rho,
                     "bootstrap_n": len(bs_r), "ci95_lo": ci_lo, "ci95_hi": ci_hi,
                     "ci_width": ci_hi - ci_lo})
    out = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
