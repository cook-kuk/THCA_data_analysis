#!/usr/bin/env python3
"""DM1 cross-cohort robustness deep dive — 2026-05-08.

Pulls pre-computed 8-gene scores from three independent VERIFIED-OK cohorts:
  - GSE286332 (Korean PTC vs PTC+HT, n=18 RNA)
  - TCGA-THCA (n=~500 RNA, DM1 vs DM2 + sub-A vs sub-B)
  - Mun 2025 (n=336 protein, PTC -> PDTC -> ATC dediff axis)

Computes per-cohort group effect sizes (Cohen's d, Hedges' g, 95% CI bootstrap)
and emits:
  - dm1_robustness_v2026_05_08/effect_sizes.tsv
  - dm1_robustness_v2026_05_08/per_sample_panel.tsv (consolidated)
  - dm1_robustness_v2026_05_08/forest_plot.png
  - dm1_robustness_v2026_05_08/summary.json

Paper-1 paper-blocking only — no Track B work, no voice-protected prose.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "project" / "results" / "dm1_robustness_v2026_05_08"
OUT.mkdir(parents=True, exist_ok=True)

GSE286332_PANEL = REPO / "project" / "results" / "p3_gse286332" / "8gene_panel_per_sample.tsv"
TCGA_SUB = REPO / "project" / "results" / "d6p7_dm1_subcluster" / "subcluster_scores.tsv"
TCGA_SIG = REPO / "project" / "results" / "d4p2_tcga_hashimoto_signature" / "tcga_signature_scores.tsv"
MUN_PANEL = REPO / "project" / "results" / "proteogenomic_v1" / "paper3_mun2025_dediff_layer" / "eight_gene_protein_per_sample.tsv"


def cohens_d(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Return Cohen's d and Hedges' g (small-sample bias-corrected)."""
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return float("nan"), float("nan")
    sx2, sy2 = x.var(ddof=1), y.var(ddof=1)
    sp = math.sqrt(((nx - 1) * sx2 + (ny - 1) * sy2) / (nx + ny - 2))
    if sp == 0:
        return float("nan"), float("nan")
    d = (x.mean() - y.mean()) / sp
    j = 1 - (3 / (4 * (nx + ny) - 9))
    return d, d * j


def bootstrap_ci(x: np.ndarray, y: np.ndarray, n_boot: int = 2000, seed: int = 42) -> tuple[float, float]:
    """Return 95% CI on Cohen's d via paired-bootstrap."""
    rng = np.random.default_rng(seed)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) < 2 or len(y) < 2:
        return float("nan"), float("nan")
    boots = []
    for _ in range(n_boot):
        bx = rng.choice(x, size=len(x), replace=True)
        by = rng.choice(y, size=len(y), replace=True)
        d, _ = cohens_d(bx, by)
        if math.isfinite(d):
            boots.append(d)
    if not boots:
        return float("nan"), float("nan")
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return float(lo), float(hi)


def main() -> None:
    rows = []
    per_sample = []

    # ----- GSE286332 — PTC vs PTC+HT -----
    g = pd.read_csv(GSE286332_PANEL, sep="\t", index_col=0)
    g_ptc = g.loc[g["group"] == "PTC", "RAI_score_8gene"].to_numpy()
    g_ht = g.loc[g["group"].isin(["PTC_HT", "PTC+HT"]), "RAI_score_8gene"].to_numpy()
    d, gh = cohens_d(g_ptc, g_ht)
    lo, hi = bootstrap_ci(g_ptc, g_ht)
    rows.append({
        "cohort": "GSE286332 (Korean RNA n=18)",
        "contrast": "PTC vs PTC+HT",
        "n_a": int(len(g_ptc)),
        "n_b": int(len(g_ht)),
        "mean_a": float(g_ptc.mean()),
        "mean_b": float(g_ht.mean()),
        "cohens_d": float(d),
        "hedges_g": float(gh),
        "ci_lo": lo,
        "ci_hi": hi,
        "score_axis": "RAI_score_8gene (within-sample z, mean over 8 genes)",
    })
    for s, row in g.iterrows():
        per_sample.append({
            "cohort": "GSE286332",
            "sample": s,
            "group": str(row.get("group", "")),
            "score": float(row.get("RAI_score_8gene", float("nan"))),
        })

    # ----- TCGA — sub_A vs sub_B (DM1 sub-cluster) -----
    t = pd.read_csv(TCGA_SUB, sep="\t", index_col=0)
    t_a = t.loc[t["sub"] == "sub_A", "g8_RAI"].to_numpy()
    t_b = t.loc[t["sub"] == "sub_B", "g8_RAI"].to_numpy()
    d, gh = cohens_d(t_a, t_b)
    lo, hi = bootstrap_ci(t_a, t_b)
    rows.append({
        "cohort": "TCGA-THCA (sub-cluster n=140)",
        "contrast": "sub_A vs sub_B (DM1 only)",
        "n_a": int(len(t_a)),
        "n_b": int(len(t_b)),
        "mean_a": float(t_a.mean()),
        "mean_b": float(t_b.mean()),
        "cohens_d": float(d),
        "hedges_g": float(gh),
        "ci_lo": lo,
        "ci_hi": hi,
        "score_axis": "g8_RAI (8-gene mean z)",
    })
    for s, row in t.iterrows():
        per_sample.append({
            "cohort": "TCGA-DM1-sub",
            "sample": str(s),
            "group": str(row.get("sub", "")),
            "score": float(row.get("g8_RAI", float("nan"))),
        })

    # ----- TCGA — DM1 vs DM2 -----
    sig = pd.read_csv(TCGA_SIG, sep="\t", index_col=0)
    if "DM" in sig.columns:
        merged = t.join(sig[["DM"]], how="left")
        s_dm1 = merged.loc[merged["DM"] == "DM1", "g8_RAI"].dropna().to_numpy()
        s_dm2 = merged.loc[merged["DM"] == "DM2", "g8_RAI"].dropna().to_numpy()
        if len(s_dm1) > 1 and len(s_dm2) > 1:
            d, gh = cohens_d(s_dm1, s_dm2)
            lo, hi = bootstrap_ci(s_dm1, s_dm2)
            rows.append({
                "cohort": "TCGA-THCA (DM1+DM2 join n=" + str(len(s_dm1) + len(s_dm2)) + ")",
                "contrast": "DM1 vs DM2",
                "n_a": int(len(s_dm1)),
                "n_b": int(len(s_dm2)),
                "mean_a": float(s_dm1.mean()),
                "mean_b": float(s_dm2.mean()),
                "cohens_d": float(d),
                "hedges_g": float(gh),
                "ci_lo": lo,
                "ci_hi": hi,
                "score_axis": "g8_RAI (TCGA sub-cluster) joined to DM call",
            })

    # TCGA whole-cohort g8_RAI vs hashi_otsu (HT-overlap)
    if "hashi_otsu" in sig.columns and "g8_RAI" in t.columns:
        merged = t.join(sig[["hashi_otsu"]], how="left")
        s_ht = merged.loc[merged["hashi_otsu"] == 1, "g8_RAI"].dropna().to_numpy()
        s_no = merged.loc[merged["hashi_otsu"] == 0, "g8_RAI"].dropna().to_numpy()
        if len(s_ht) > 1 and len(s_no) > 1:
            d, gh = cohens_d(s_no, s_ht)
            lo, hi = bootstrap_ci(s_no, s_ht)
            rows.append({
                "cohort": "TCGA-THCA (sub-cluster + HT-overlap call)",
                "contrast": "hashi_otsu=0 vs hashi_otsu=1",
                "n_a": int(len(s_no)),
                "n_b": int(len(s_ht)),
                "mean_a": float(s_no.mean()),
                "mean_b": float(s_ht.mean()),
                "cohens_d": float(d),
                "hedges_g": float(gh),
                "ci_lo": lo,
                "ci_hi": hi,
                "score_axis": "g8_RAI (TCGA sub-cluster) joined to hashi_otsu",
            })

    # ----- Mun 2025 protein — PTC vs PDTC vs ATC dediff axis -----
    m = pd.read_csv(MUN_PANEL, sep="\t", index_col=0)
    panel_cols = [c for c in m.columns if c not in ("group",)]
    m["g8_protein_mean"] = m[panel_cols].mean(axis=1)
    m_ptc = m.loc[m["group"] == "PTC", "g8_protein_mean"].to_numpy()
    m_pdtc = m.loc[m["group"] == "PDTC", "g8_protein_mean"].to_numpy()
    m_atc = m.loc[m["group"] == "ATC", "g8_protein_mean"].to_numpy()
    for label, ref, alt in [
        ("PTC vs ATC", m_ptc, m_atc),
        ("PTC vs PDTC", m_ptc, m_pdtc),
        ("PDTC vs ATC", m_pdtc, m_atc),
    ]:
        if len(ref) > 1 and len(alt) > 1:
            d, gh = cohens_d(ref, alt)
            lo, hi = bootstrap_ci(ref, alt)
            rows.append({
                "cohort": "Mun 2025 (protein n=" + str(len(m_ptc) + len(m_pdtc) + len(m_atc)) + ")",
                "contrast": label,
                "n_a": int(len(ref)),
                "n_b": int(len(alt)),
                "mean_a": float(ref.mean()),
                "mean_b": float(alt.mean()),
                "cohens_d": float(d),
                "hedges_g": float(gh),
                "ci_lo": lo,
                "ci_hi": hi,
                "score_axis": "8-gene protein z mean",
            })
    for s, row in m.iterrows():
        per_sample.append({
            "cohort": "Mun2025-protein",
            "sample": str(s),
            "group": str(row.get("group", "")),
            "score": float(row.get("g8_protein_mean", float("nan"))),
        })

    # ----- Save -----
    es = pd.DataFrame(rows)
    es.to_csv(OUT / "effect_sizes.tsv", sep="\t", index=False)
    pd.DataFrame(per_sample).to_csv(OUT / "per_sample_panel.tsv", sep="\t", index=False)

    summary = {
        "generated_at": "2026-05-08",
        "n_cohorts": int(es["cohort"].nunique()),
        "n_contrasts": len(es),
        "all_d_concordant_sign": bool(((es["cohens_d"] > 0).all()) or ((es["cohens_d"] < 0).all())),
        "median_abs_d": float(es["cohens_d"].abs().median()),
        "max_abs_d": float(es["cohens_d"].abs().max()),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))

    # ----- Forest plot -----
    fig, ax = plt.subplots(figsize=(9, max(3, 0.45 * len(rows) + 1.5)))
    y = np.arange(len(es))
    for i, r in es.iterrows():
        ax.errorbar(r["cohens_d"], i, xerr=[[r["cohens_d"] - r["ci_lo"]], [r["ci_hi"] - r["cohens_d"]]],
                    fmt="s", color="#222", ecolor="#666", capsize=4, markersize=6)
    ax.axvline(0, color="#888", linestyle=":", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['cohort']}\n  {r['contrast']} (n_a={r['n_a']}, n_b={r['n_b']})" for _, r in es.iterrows()], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Cohen's d (RAI / 8-gene panel)\nleft = lower in group_a, right = higher in group_a")
    ax.set_title("DM1 8-gene panel cross-cohort robustness — 2026-05-08\n(positive d = group_a > group_b along RAI / thyroid-differentiation axis)", fontsize=10)
    ax.grid(axis="x", linestyle=":", alpha=0.3)
    plt.tight_layout()
    fig.savefig(OUT / "forest_plot.png", dpi=160)
    fig.savefig(OUT / "forest_plot.pdf")
    plt.close(fig)

    print(f"wrote {OUT.relative_to(REPO)}/effect_sizes.tsv ({len(es)} rows)")
    print(f"wrote forest_plot.png + forest_plot.pdf + summary.json")
    print(f"\nEffect sizes:")
    print(es.to_string(index=False))


if __name__ == "__main__":
    main()
