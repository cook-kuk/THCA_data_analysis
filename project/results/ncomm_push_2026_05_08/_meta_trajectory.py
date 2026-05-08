#!/usr/bin/env python3
"""Cohort-level Cohen's d × DerSimonian-Laird random-effects meta-analysis
for the stage trajectory.

The earlier within-cohort-anchored pooled trajectory mixed cohorts with
different stage compositions, which artefactually inflated PDTC. The right
operation is: compute Cohen's d for each stage contrast WITHIN each cohort
that has both stages, then meta-pool the d-values across cohorts.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parent

CONTRASTS = [
    ("ATC", "PDTC"),
    ("ATC", "PTC"),
    ("ATC", "normal"),
    ("PDTC", "PTC"),
    ("PDTC", "normal"),
    ("PTC", "normal"),
]


def cohens_d(a: pd.Series, b: pd.Series) -> tuple[float, float]:
    a = a.dropna()
    b = b.dropna()
    if len(a) < 2 or len(b) < 2:
        return float("nan"), float("nan")
    sd_pool = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1))
                      / (len(a) + len(b) - 2))
    if sd_pool == 0:
        return float("nan"), float("nan")
    d = (a.mean() - b.mean()) / sd_pool
    # Hedges-Welch SE for Cohen's d
    n1, n2 = len(a), len(b)
    se = np.sqrt((n1 + n2) / (n1 * n2) + d ** 2 / (2 * (n1 + n2)))
    return float(d), float(se)


def random_effects_pool(d_list: list[float], se_list: list[float]) -> dict:
    """DerSimonian-Laird random-effects meta-analysis on Cohen's d values."""
    d = np.array(d_list, dtype=float)
    se = np.array(se_list, dtype=float)
    keep = np.isfinite(d) & np.isfinite(se) & (se > 0)
    d, se = d[keep], se[keep]
    if len(d) == 0:
        return {"n_studies": 0}
    w_fe = 1.0 / (se ** 2)
    d_fe = float((w_fe * d).sum() / w_fe.sum())
    Q = float((w_fe * (d - d_fe) ** 2).sum())
    df = len(d) - 1
    if df > 0 and Q > df:
        c = float(w_fe.sum() - (w_fe ** 2).sum() / w_fe.sum())
        tau2 = max(0.0, (Q - df) / c) if c > 0 else 0.0
    else:
        tau2 = 0.0
    w_re = 1.0 / (se ** 2 + tau2)
    d_re = float((w_re * d).sum() / w_re.sum())
    se_re = float(np.sqrt(1.0 / w_re.sum()))
    z = d_re / se_re
    # two-sided z->p
    from scipy.stats import norm
    p = float(2 * (1 - norm.cdf(abs(z))))
    i2 = float(max(0.0, (Q - df) / Q) * 100) if Q > 0 else 0.0
    return {
        "n_studies": int(len(d)),
        "d_re": d_re,
        "se_re": se_re,
        "ci95_lo": d_re - 1.96 * se_re,
        "ci95_hi": d_re + 1.96 * se_re,
        "z": float(z),
        "p": p,
        "Q": Q,
        "tau2": tau2,
        "I2_pct": i2,
    }


def main() -> None:
    pooled = pd.read_csv(ROOT / "external_pooled_scores.tsv", sep="\t", index_col=0)
    pooled = pooled.loc[~pooled.index.duplicated(keep="last")]

    geo = pooled[pooled["cohort"] != "TCGA-PTC-GDC"].copy()

    # We meta-pool across cohorts using the RAW panel_log2 column. Cohen's d
    # within a cohort cancels the cohort-specific intercept automatically.

    rows = []
    for cohort in sorted(geo["cohort"].unique()):
        sub = geo[geo["cohort"] == cohort]
        for s1, s2 in CONTRASTS:
            v1 = sub.loc[sub["stage"] == s1, "panel_log2"]
            v2 = sub.loc[sub["stage"] == s2, "panel_log2"]
            if len(v1.dropna()) < 2 or len(v2.dropna()) < 2:
                continue
            d, se = cohens_d(v1, v2)
            try:
                stat, p = mannwhitneyu(v1.dropna(), v2.dropna(), alternative="two-sided")
                p = float(p)
            except ValueError:
                p = float("nan")
            rows.append({
                "cohort": cohort,
                "contrast": f"{s1}_vs_{s2}",
                "n_first": int(len(v1.dropna())),
                "n_second": int(len(v2.dropna())),
                "mean_first_log2": float(v1.mean()),
                "mean_second_log2": float(v2.mean()),
                "delta_log2": float(v1.mean() - v2.mean()),
                "cohens_d": d,
                "se_d": se,
                "mw_p": p,
            })
    per_cohort = pd.DataFrame(rows)
    per_cohort.to_csv(ROOT / "meta_per_cohort.tsv", sep="\t", index=False)

    # Random-effects pool across cohorts for each contrast
    pool_rows = []
    for contrast in [f"{s1}_vs_{s2}" for s1, s2 in CONTRASTS]:
        sub = per_cohort[per_cohort["contrast"] == contrast]
        if sub.empty:
            continue
        meta = random_effects_pool(sub["cohens_d"].tolist(), sub["se_d"].tolist())
        meta["contrast"] = contrast
        meta["cohorts_contributing"] = sub["cohort"].tolist()
        pool_rows.append(meta)
    pool_df = pd.DataFrame(pool_rows)
    pool_df.to_csv(ROOT / "meta_pooled.tsv", sep="\t", index=False)
    (ROOT / "meta_pooled.json").write_text(json.dumps(pool_rows, indent=2, default=str))

    # Add CCLE cell-line contrasts to the same JSON for completeness
    ccle_path = ROOT / "ccle_thyroid" / "thap_vs_thpa_cohens_d.json"
    if ccle_path.exists():
        ccle = json.loads(ccle_path.read_text())
        ccle_rows = []
        for k, v in ccle.items():
            ccle_rows.append({
                "source": "ccle_broad_2025",
                "contrast": k,
                **v,
            })
        (ROOT / "ccle_contrasts.json").write_text(json.dumps(ccle_rows, indent=2, default=str))

    print("=== Per-cohort Cohen's d ===")
    print(per_cohort.to_string(index=False))
    print()
    print("=== Random-effects meta-pool ===")
    for r in pool_rows:
        print(f"  {r['contrast']:18s}  k={r['n_studies']}  "
              f"d_re={r['d_re']:+.3f} [95% CI {r['ci95_lo']:+.3f}, {r['ci95_hi']:+.3f}]  "
              f"p={r['p']:.3g}  I²={r['I2_pct']:.0f}%")

    # Forest-plot–style figure
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # focus on ATC vs PDTC and ATC vs normal/PTC — the dedifferentiation contrasts
        fig, ax = plt.subplots(figsize=(7.5, 5.0))
        y = 0
        labels, ds, lo, hi = [], [], [], []
        for contrast in ["ATC_vs_normal", "ATC_vs_PTC", "ATC_vs_PDTC", "PDTC_vs_PTC", "PDTC_vs_normal", "PTC_vs_normal"]:
            sub = per_cohort[per_cohort["contrast"] == contrast]
            for _, r in sub.iterrows():
                ax.errorbar(r["cohens_d"], y,
                            xerr=1.96 * r["se_d"],
                            fmt="s", ms=6, capsize=4, color="#3b6ea8")
                ax.text(r["cohens_d"], y, f"  {r['cohort']} (n={int(r['n_first'])}/{int(r['n_second'])})",
                        va="center", fontsize=7, color="#333")
                labels.append(f"{contrast} · {r['cohort']}")
                ds.append(r["cohens_d"])
                lo.append(r["cohens_d"] - 1.96 * r["se_d"])
                hi.append(r["cohens_d"] + 1.96 * r["se_d"])
                y -= 1
            # pooled diamond
            pool = next((p for p in pool_rows if p["contrast"] == contrast), None)
            if pool and pool["n_studies"] > 0:
                ax.plot(pool["d_re"], y, marker="D", ms=12, color="#b03a2e",
                        markeredgecolor="black")
                ax.errorbar(pool["d_re"], y, xerr=1.96 * pool["se_re"],
                            fmt="none", capsize=6, color="#b03a2e", lw=2)
                ax.text(pool["d_re"], y, f"  POOLED {contrast} d={pool['d_re']:+.2f}, p={pool['p']:.2g}, I²={pool['I2_pct']:.0f}%",
                        va="center", fontsize=8, color="#b03a2e", weight="bold")
                y -= 1
            y -= 0.5
        ax.axvline(0, ls="--", c="#888", lw=0.8)
        ax.set_xlabel("Cohen's d (panel log2 mean difference / pooled SD)")
        ax.set_yticks([])
        ax.set_xlim(-3.5, 2.5)
        ax.set_title("External GEO cohort meta-analysis: 8-gene panel suppression by stage")
        plt.tight_layout()
        (ROOT / "figures").mkdir(exist_ok=True)
        plt.savefig(ROOT / "figures" / "stage_cohens_d_forest.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        print("\nfigure saved: figures/stage_cohens_d_forest.png")
    except Exception as exc:
        print(f"figure block failed: {exc}")


if __name__ == "__main__":
    main()
