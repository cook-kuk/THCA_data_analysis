"""
Vaccine + ICI synergy stratification for Korean PDAC.

The core hypothesis:
  vaccine + ICI works best when (a) antigen presentation is active (inflamed
  axis) AND (b) effector function is suppressed (suppress axis).
  Single-agent vaccine fails in (b)-heavy patients; ICI alone fails in
  (a)-cold patients. The intersection — inflamed AND suppressed — is the
  empirical target.

Per patient on TCGA-PAAD bulk RNA, we compute:
  1) inflamed score = mean(TLS + IFNG + HLA-II)
  2) suppress score = mean(myeloid + checkpoint)
  3) 4-quadrant assignment (Q1=both-high · vaccine+ICI+myeloid-modifier
                            Q2=inflamed only · vaccine OK
                            Q3=suppressed only · ICI + chemo
                            Q4=both-low · chemo only)
  4) KRAS allele (G12D / G12V / G12R / G12C / WT)
  5) Korean off-the-shelf cassette eligibility
        = mut prevalence × Korean HLA frequency for the allele's restriction
  6) Composite "Korean-aware vaccine + ICI eligibility score"
        = is_Q1 × allele_actionability × Korean_HLA_match

Outputs:
  results/synergy/per_patient_eligibility.tsv     ← per-patient score
  results/synergy/SUMMARY.json                    ← cohort-level
  figures/fig_quadrant_kras.png                   ← 2D quadrant × allele
  figures/fig_eligibility_waterfall.png           ← cohort eligibility
  figures/fig_korean_combo_eligibility.png        ← Korean-specific
"""

import json
import math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = Path("/data/pdac_poc")
RES = ROOT / "results/synergy"
RES.mkdir(parents=True, exist_ok=True)
FIG = ROOT / "figures"

EM, RO, VI, TL, AM, BL = "#1c8e6d", "#7B1F2A", "#4f2db5", "#206a78", "#b58534", "#244e73"
plt.rcParams.update({"font.family":"sans-serif","font.size":10,
                     "axes.spines.top":False,"axes.spines.right":False,
                     "figure.facecolor":"white","axes.facecolor":"white"})

# Korean HLA actionability per KRAS allele (joint mut × HLA frequency)
# Same priors as scripts 06 + 12 — AFND v3.1
KOR_HLA = {
    "G12D": {"alleles": ["A*11:01", "A*03:01", "B*07:02", "DRB1*04"],
              "korean_pct": 31.5, "european_pct": 26.4},
    "G12V": {"alleles": ["A*11:01", "B*07:02", "DRB1*0701"],
              "korean_pct": 28.1, "european_pct": 23.0},
    "G12R": {"alleles": ["B*07:02", "DRB1*04:01"],
              "korean_pct": 8.9,  "european_pct": 16.6},
    "G12C": {"alleles": ["A*02:01", "B*07:02"],
              "korean_pct": 22.0, "european_pct": 38.0},
}


def load_data():
    moff = pd.read_csv(ROOT/"results/moffitt_calls.tsv", sep="\t", index_col=0)
    immune = pd.read_csv(ROOT/"results/pdac_immune_readiness_per_sample.tsv", sep="\t", index_col=0)
    kras = pd.read_csv(ROOT/"results/kras_allele_table.tsv", sep="\t", index_col=0)
    df = moff[["moffitt_call","mut_KRAS","mut_TP53","mut_SMAD4"]].copy()
    df = df.join(immune[["myeloid_suppressive","myCAF","iCAF","HLA_I","HLA_II",
                          "IFNG_inflamed","checkpoint_exhaustion","TLS_CXCL13"]],
                 how="left")
    df["kras_allele"] = kras["kras_allele"]
    df["os_months"] = kras["os_months"]
    df["os_event"]  = kras["event"]
    return df


def assign_quadrant(df):
    df = df.copy()
    df["score_inflamed"] = df[["TLS_CXCL13","IFNG_inflamed","HLA_II"]].mean(axis=1)
    df["score_suppress"] = df[["myeloid_suppressive","checkpoint_exhaustion"]].mean(axis=1)
    inf_mid = df["score_inflamed"].median()
    sup_mid = df["score_suppress"].median()
    df["q_inflamed"] = (df["score_inflamed"] > inf_mid).astype(int)
    df["q_suppress"] = (df["score_suppress"] > sup_mid).astype(int)
    df["quadrant"] = (
        df["q_inflamed"].astype(str) + df["q_suppress"].astype(str)
    ).map({
        "11": "Q1_both_high",          # vaccine + ICI + myeloid-modifier
        "10": "Q2_inflamed_only",      # vaccine alone OK
        "01": "Q3_suppressed_only",    # ICI + chemo
        "00": "Q4_both_low",           # chemo only
    })
    return df, inf_mid, sup_mid


def annotate_combo(df):
    """Per-patient combo recommendation from paradox quadrant + KRAS allele."""
    df = df.copy()
    rec_map = {
        "Q1_both_high": "vaccine + ICI + myeloid-modifier (★ ELI-002 + atezo + anti-CSF1R)",
        "Q2_inflamed_only": "vaccine (BNT122) + atezolizumab",
        "Q3_suppressed_only": "ICI + chemo (anti-PD-1 + mFOLFIRINOX)",
        "Q4_both_low": "chemo only (mFOLFIRINOX) — vaccine+ICI not indicated",
    }
    df["combo_recommendation"] = df["quadrant"].map(rec_map)

    # Korean cassette eligibility per allele
    def actionability_korean(al):
        return KOR_HLA.get(al, {}).get("korean_pct", 0.0)
    def actionability_european(al):
        return KOR_HLA.get(al, {}).get("european_pct", 0.0)
    df["korean_HLA_match_pct"]   = df["kras_allele"].map(actionability_korean)
    df["european_HLA_match_pct"] = df["kras_allele"].map(actionability_european)

    # Composite Korean eligibility score
    df["is_Q1"] = (df["quadrant"] == "Q1_both_high").astype(int)
    df["is_Q1or2"] = df["quadrant"].isin(["Q1_both_high","Q2_inflamed_only"]).astype(int)
    # 0..1 scale
    df["korean_combo_score"] = (
        0.5 * df["is_Q1or2"] +
        0.4 * (df["korean_HLA_match_pct"] / 31.5) +
        0.1 * df["mut_KRAS"]
    )
    return df


def cohort_summary(df):
    quad_counts = df["quadrant"].value_counts().to_dict()
    summary = {
        "n_total": int(len(df)),
        "n_with_RNA": int(df["score_inflamed"].notna().sum()),
        "quadrant_counts": quad_counts,
        "fraction_per_quadrant": {k: round(v/len(df),3) for k,v in quad_counts.items()},
    }

    # By KRAS allele × quadrant
    ct = pd.crosstab(df["kras_allele"], df["quadrant"]).to_dict()
    summary["allele_x_quadrant"] = ct

    # Per-allele: fraction in Q1 (vaccine + ICI + myeloid-modifier candidate)
    summary["fraction_Q1_per_allele"] = {}
    for al, sub in df.groupby("kras_allele"):
        if len(sub) >= 5:
            summary["fraction_Q1_per_allele"][al] = round((sub["quadrant"]=="Q1_both_high").mean(), 3)

    # Per-quadrant: median KRAS-G12D fraction
    summary["G12D_fraction_per_quadrant"] = {}
    for q, sub in df.groupby("quadrant"):
        summary["G12D_fraction_per_quadrant"][q] = round((sub["kras_allele"]=="G12D").mean(), 3)

    # Korean off-the-shelf eligibility waterfall
    n_total = len(df)
    n_Q1 = int((df["quadrant"]=="Q1_both_high").sum())
    n_Q1_or_Q2 = int(df["quadrant"].isin(["Q1_both_high","Q2_inflamed_only"]).sum())
    n_kras_mut = int(df["mut_KRAS"].sum())
    n_kras_g12_actionable = int(df["kras_allele"].isin(["G12D","G12V","G12R","G12C"]).sum())
    # Korean joint-actionable patients = sum over allele of (n_allele × Korean_HLA_pct)
    kor_actionable = 0
    for al, info in KOR_HLA.items():
        n_al = int((df["kras_allele"]==al).sum())
        kor_actionable += n_al * (info["korean_pct"]/100)
    summary["waterfall"] = {
        "all_PDAC": n_total,
        "with_RNA": int(df["score_inflamed"].notna().sum()),
        "Q1_vaccine_ICI_combo_candidates": n_Q1,
        "Q1_or_Q2_vaccine_eligible": n_Q1_or_Q2,
        "KRAS_mutant": n_kras_mut,
        "KRAS_G12_hotspot": n_kras_g12_actionable,
        "Q1_AND_KRAS_G12_hotspot": int(((df["quadrant"]=="Q1_both_high") &
                                          df["kras_allele"].isin(["G12D","G12V","G12R","G12C"])).sum()),
        "Q1_AND_G12D": int(((df["quadrant"]=="Q1_both_high") & (df["kras_allele"]=="G12D")).sum()),
        "Korean_off_the_shelf_actionable_avg": round(kor_actionable, 2),
        "Korean_off_the_shelf_pct_of_cohort": round(kor_actionable/n_total*100, 1),
    }

    # Survival per quadrant (does Q1 actually have worse OS without combo?)
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import multivariate_logrank_test
    surv = df.dropna(subset=["os_months","os_event"]).query("os_months > 0").copy()
    surv["event"] = surv["os_event"].astype(int)
    medians = {}
    for q in ["Q1_both_high","Q2_inflamed_only","Q3_suppressed_only","Q4_both_low"]:
        sub = surv[surv["quadrant"]==q]
        if len(sub) >= 5:
            kmf = KaplanMeierFitter().fit(sub["os_months"], sub["event"])
            m = kmf.median_survival_time_
            medians[q] = float(m) if not (isinstance(m, float) and np.isnan(m)) else None
    summary["median_OS_per_quadrant"] = medians
    if surv["quadrant"].nunique() >= 3:
        lr = multivariate_logrank_test(surv["os_months"], surv["quadrant"], surv["event"])
        summary["logrank_4_quadrant_p"] = float(lr.p_value)
        summary["logrank_4_quadrant_test_stat"] = float(lr.test_statistic)

    return summary


def fig_quadrant(df, inf_mid, sup_mid):
    """Figure: 2D paradox quadrant with KRAS allele overlay."""
    fig, ax = plt.subplots(figsize=(7.5, 6.0), dpi=160)
    color_map = {"G12D":RO, "G12V":AM, "G12R":EM, "G12C":VI, "WT":"#aaa",
                 "G12_other":"#888","Q61":"#666","KRAS_other":"#999"}
    marker_size = 38
    for al in ["WT","G12_other","Q61","KRAS_other","G12C","G12R","G12V","G12D"]:
        sub = df[df["kras_allele"]==al]
        if len(sub) == 0: continue
        ax.scatter(sub["score_inflamed"], sub["score_suppress"],
                   c=color_map.get(al,"#888"), s=marker_size,
                   alpha=0.85 if al=="G12D" else 0.7,
                   edgecolors="white", linewidths=0.6,
                   label=f"{al} (n={len(sub)})")
    # Median crosshair
    ax.axhline(sup_mid, color="#888", lw=0.8, linestyle="--")
    ax.axvline(inf_mid, color="#888", lw=0.8, linestyle="--")

    # Quadrant labels
    xmin, xmax = df["score_inflamed"].min()-0.1, df["score_inflamed"].max()+0.1
    ymin, ymax = df["score_suppress"].min()-0.1, df["score_suppress"].max()+0.1
    ax.text(xmax*0.9, ymax*0.9, "Q1\nvaccine + ICI\n+ myeloid-modifier",
            ha="right", va="top", fontsize=10, fontweight="bold", color=RO,
            bbox=dict(boxstyle="round,pad=0.4", fc="#fce5e3", ec=RO, lw=1.5))
    ax.text(xmax*0.9, ymin*0.9, "Q2\nvaccine alone OK",
            ha="right", va="bottom", fontsize=10, fontweight="bold", color=EM,
            bbox=dict(boxstyle="round,pad=0.4", fc="#dff0dd", ec=EM, lw=1.5))
    ax.text(xmin*0.9, ymax*0.9, "Q3\nICI + chemo",
            ha="left", va="top", fontsize=10, fontweight="bold", color=BL,
            bbox=dict(boxstyle="round,pad=0.4", fc="#dde9f4", ec=BL, lw=1.5))
    ax.text(xmin*0.9, ymin*0.9, "Q4\nchemo only",
            ha="left", va="bottom", fontsize=10, fontweight="bold", color="#666",
            bbox=dict(boxstyle="round,pad=0.4", fc="#eee", ec="#888", lw=1.5))

    ax.set_xlabel("Inflamed score (TLS + IFN-γ + HLA-II)")
    ax.set_ylabel("Suppress score (myeloid + checkpoint)")
    ax.set_title("Vaccine + ICI eligibility quadrants × KRAS allele · TCGA-PAAD")
    ax.legend(loc="lower right", frameon=False, fontsize=8.5,
              bbox_to_anchor=(1.0, -0.05))
    fig.tight_layout()
    fig.savefig(FIG/"fig_quadrant_kras.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIG/"fig_quadrant_kras.svg", bbox_inches="tight")
    plt.close(fig)


def fig_waterfall(s):
    """Eligibility waterfall."""
    w = s["waterfall"]
    stages = [
        ("All TCGA-PAAD",                      w["all_PDAC"], "#0c1422"),
        ("With bulk-RNA-scored TME",           w["with_RNA"], "#244e73"),
        ("Q1 vaccine+ICI candidates",          w["Q1_vaccine_ICI_combo_candidates"], "#7B1F2A"),
        ("Q1 OR Q2 vaccine-eligible",          w["Q1_or_Q2_vaccine_eligible"], "#206a78"),
        ("KRAS mutant (any)",                  w["KRAS_mutant"], "#b58534"),
        ("KRAS G12 hotspot",                   w["KRAS_G12_hotspot"], "#1c8e6d"),
        ("Q1 AND KRAS G12 hotspot",            w["Q1_AND_KRAS_G12_hotspot"], "#4f2db5"),
        ("★ Q1 AND G12D specifically",          w["Q1_AND_G12D"], "#7B1F2A"),
        ("Korean off-the-shelf cassette\n(Q1 OR Q2 × KRAS G12 × Korean HLA)",
                                                w["Korean_off_the_shelf_actionable_avg"], "#4f2db5"),
    ]
    fig, ax = plt.subplots(figsize=(8.6, 4.6), dpi=160)
    y = np.arange(len(stages))[::-1]
    for i, (lbl, n, c) in enumerate(stages):
        ax.barh(y[i], n, color=c, alpha=.92, edgecolor="white", linewidth=2)
        pct = n / w["all_PDAC"] * 100 if w["all_PDAC"] else 0
        ax.text(n + 2, y[i], f"  n = {n}  ({pct:.1f}% of cohort)",
                va="center", fontsize=9.5, fontweight="bold", color=c)
    ax.set_yticks(y); ax.set_yticklabels([s[0] for s in stages], fontsize=10)
    ax.set_xlabel("samples")
    ax.set_xlim(0, w["all_PDAC"]*1.4)
    ax.set_title("Vaccine + ICI eligibility waterfall · TCGA-PAAD → Korean cassette projection")
    fig.tight_layout()
    fig.savefig(FIG/"fig_eligibility_waterfall.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIG/"fig_eligibility_waterfall.svg", bbox_inches="tight")
    plt.close(fig)


def fig_korean_combo(df, s):
    """Bar plot of Korean combo eligibility per KRAS allele."""
    rows = []
    for al, info in KOR_HLA.items():
        sub = df[df["kras_allele"]==al]
        if len(sub) == 0: continue
        n_al = len(sub)
        frac_Q1 = (sub["quadrant"]=="Q1_both_high").mean()
        frac_Q1or2 = sub["quadrant"].isin(["Q1_both_high","Q2_inflamed_only"]).mean()
        kor_combo_n = n_al * frac_Q1or2 * (info["korean_pct"]/100)
        eu_combo_n = n_al * frac_Q1or2 * (info["european_pct"]/100)
        rows.append({"allele": al, "n_total": n_al,
                     "frac_Q1or2": frac_Q1or2,
                     "korean_combo_n": kor_combo_n,
                     "european_combo_n": eu_combo_n,
                     "korean_pct_HLA": info["korean_pct"],
                     "european_pct_HLA": info["european_pct"]})
    rows.sort(key=lambda r: -r["korean_combo_n"])
    fig, ax = plt.subplots(figsize=(8.5, 4.0), dpi=160)
    x = np.arange(len(rows)); w = 0.38
    ax.bar(x - w/2, [r["korean_combo_n"] for r in rows], w, color=RO, label="Korean", alpha=.92)
    ax.bar(x + w/2, [r["european_combo_n"] for r in rows], w, color=BL, label="European", alpha=.92)
    for i, r in enumerate(rows):
        ax.text(i - w/2, r["korean_combo_n"]+0.2, f"{r['korean_combo_n']:.1f}",
                ha="center", fontsize=9, color=RO, fontweight="bold")
        ax.text(i + w/2, r["european_combo_n"]+0.2, f"{r['european_combo_n']:.1f}",
                ha="center", fontsize=9, color=BL)
    ax.set_xticks(x)
    ax.set_xticklabels([r["allele"] for r in rows], fontsize=11, fontweight="bold")
    ax.set_ylabel("expected vaccine+ICI eligible patients\n(per cohort × allele × HLA)")
    ax.set_title("Korean vs European vaccine+ICI eligibility per KRAS allele · TCGA-PAAD scaled")
    ax.legend(loc="upper right", frameon=False)
    fig.tight_layout()
    fig.savefig(FIG/"fig_korean_combo_eligibility.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIG/"fig_korean_combo_eligibility.svg", bbox_inches="tight")
    plt.close(fig)


def main():
    print("[22] Vaccine + ICI synergy stratification · Korean PDAC")
    print("="*70)
    df = load_data()
    print(f"  n_total = {len(df)}  · with RNA = {df['score_inflamed'].notna().sum() if 'score_inflamed' in df.columns else 0}")
    df, inf_mid, sup_mid = assign_quadrant(df)
    print(f"  quadrant midpoints · inflamed = {inf_mid:.2f}  suppress = {sup_mid:.2f}")
    df = annotate_combo(df)

    s = cohort_summary(df)
    print(f"\n  Quadrant distribution:")
    for q, n in s["quadrant_counts"].items():
        print(f"     {q:25s}  n = {n}  ({s['fraction_per_quadrant'][q]*100:.1f}%)")
    print(f"\n  Median OS per quadrant (months):")
    for q, m in s.get("median_OS_per_quadrant", {}).items():
        print(f"     {q:25s}  {m}")
    print(f"\n  4-quadrant logrank p = {s.get('logrank_4_quadrant_p')}")
    print(f"\n  Fraction in Q1 (vaccine+ICI combo) per allele:")
    for al, frac in s["fraction_Q1_per_allele"].items():
        print(f"     {al:12s}  {frac*100:.1f}%")
    print(f"\n  Waterfall:")
    for k, v in s["waterfall"].items():
        print(f"     {k:55s}  {v}")

    df.to_csv(RES/"per_patient_eligibility.tsv", sep="\t")
    json.dump(s, open(RES/"SUMMARY.json","w"), indent=2, default=str)

    fig_quadrant(df, inf_mid, sup_mid)
    fig_waterfall(s)
    fig_korean_combo(df, s)
    print(f"\n[22] artifacts → {RES}/  + 3 figures in {FIG}/")
    print(f"     fig_quadrant_kras.png · fig_eligibility_waterfall.png · fig_korean_combo_eligibility.png")


if __name__ == "__main__":
    main()
