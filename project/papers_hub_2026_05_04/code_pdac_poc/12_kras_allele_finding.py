"""
KRAS allele × Moffitt × TME × neoantigen × HLA — finding hunt.

Tracks:
  (A) Allele-specific OS (replicate Hayashi 2021 NatCancer "G12R better OS")
  (B) Allele × Moffitt subtype 2x2 / chi-square
  (C) Allele × per-module TME Cohen's d (vs all-other-KRAS-mut)
  (D) Allele × NeoQ (Balachandran R×D)
  (E) Moffitt × TME PARADOX: basal-like simultaneously inflamed AND
      myeloid-suppressed → single-agent vaccine failure mechanism;
      quantify the paradox magnitude.

Outputs:
  results/kras_allele_finding.json  ← machine-readable finding bundle
  results/kras_allele_table.tsv     ← per-sample joined table
  figures/fig_kras_allele_*.png     ← 4 publishable figures
"""

import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, chi2_contingency, mannwhitneyu
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import multivariate_logrank_test
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/data/pdac_poc")
RAW, RES, PROC, FIG = ROOT/"raw", ROOT/"results", ROOT/"processed", ROOT/"figures"

EMERALD, ROSE, VIOLET = "#35d39d", "#ef5f79", "#9b7cff"
TEAL, AMBER, BLUE = "#32b8c6", "#f2b84b", "#65a9ff"

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "axes.facecolor": "white"
})


def load():
    mut = pd.read_csv(RAW / "mutations.tsv", sep="\t")
    moff = pd.read_csv(RES / "moffitt_calls.tsv", sep="\t", index_col=0)
    clin = pd.read_csv(RAW / "clinical.tsv", sep="\t").set_index("sampleId")
    immune = pd.read_csv(RES / "pdac_immune_readiness_per_sample.tsv",
                         sep="\t", index_col=0)
    return mut, moff, clin, immune


def assign_kras_allele(mut, moff_idx):
    """Per-sample KRAS class: G12D / G12V / G12R / G12C / Q61 / other / WT."""
    out = pd.Series("WT", index=moff_idx, name="kras_allele")
    kras = mut[mut["hugo"] == "KRAS"]
    for s, g in kras.groupby("sampleId"):
        if s not in moff_idx:
            continue
        pcs = g["proteinChange"].dropna().tolist()
        if any(pc == "G12D" for pc in pcs):
            out.loc[s] = "G12D"
        elif any(pc == "G12V" for pc in pcs):
            out.loc[s] = "G12V"
        elif any(pc == "G12R" for pc in pcs):
            out.loc[s] = "G12R"
        elif any(pc == "G12C" for pc in pcs):
            out.loc[s] = "G12C"
        elif any(pc.startswith("G12") for pc in pcs):
            out.loc[s] = "G12_other"
        elif any(pc.startswith("Q61") for pc in pcs):
            out.loc[s] = "Q61"
        else:
            out.loc[s] = "KRAS_other"
    return out


def cohens_d(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    v = ((len(a)-1)*a.var(ddof=1) + (len(b)-1)*b.var(ddof=1)) / (len(a)+len(b)-2)
    sp = np.sqrt(v) if v > 0 else np.nan
    return (a.mean() - b.mean()) / sp if sp else float("nan")


def main():
    mut, moff, clin, immune = load()
    # Build joined table
    df = moff.copy()
    df["kras_allele"] = assign_kras_allele(mut, df.index)
    df["os_months"] = pd.to_numeric(clin.get("OS_MONTHS"), errors="coerce")
    df["event"] = clin.get("OS_STATUS", pd.Series(dtype=str)).fillna("").str.startswith("1").astype(int)
    df["age"] = pd.to_numeric(clin.get("AGE"), errors="coerce")
    df["stage"] = clin.get("AJCC_PATHOLOGIC_TUMOR_STAGE", pd.Series(dtype=str))
    # join immune modules
    mod_cols = ["myeloid_suppressive", "myCAF", "iCAF", "CAF_pan",
                "HLA_I", "HLA_II", "IFNG_inflamed", "checkpoint_exhaustion",
                "TLS_CXCL13"]
    df = df.join(immune[[c for c in mod_cols if c in immune.columns]],
                 how="left")
    # NeoQ
    neoq_per = pd.read_csv(RES / "neoantigen_quality_per_sample.tsv", sep="\t").set_index("sampleId")
    df = df.join(neoq_per[["top_NeoQ", "mean_NeoQ"]], how="left")

    df.to_csv(RES / "kras_allele_table.tsv", sep="\t")

    findings = {}

    # ===== (A) Allele-specific OS =====
    surv = df.dropna(subset=["os_months"]).query("os_months > 0").copy()
    main_alleles = ["G12D", "G12V", "G12R", "G12C", "Q61", "KRAS_other", "WT"]
    surv["allele_grp"] = surv["kras_allele"].where(
        surv["kras_allele"].isin(main_alleles), "G12_other")
    counts = surv["allele_grp"].value_counts().to_dict()
    findings["A_os_counts"] = counts

    # KM medians
    medians = {}
    for al in surv["allele_grp"].unique():
        s = surv[surv["allele_grp"] == al]
        if len(s) >= 5:
            kmf = KaplanMeierFitter().fit(s["os_months"], s["event"])
            m = kmf.median_survival_time_
            medians[al] = float(m) if not (isinstance(m, float) and np.isnan(m)) else None
    findings["A_median_OS_months"] = medians

    # multi-group logrank G12D vs G12V vs G12R
    target = surv[surv["allele_grp"].isin(["G12D", "G12V", "G12R"])]
    if target["allele_grp"].nunique() == 3 and len(target) >= 30:
        lr = multivariate_logrank_test(target["os_months"],
                                        target["allele_grp"], target["event"])
        findings["A_logrank_G12DVR"] = {"p_value": float(lr.p_value),
                                        "test_statistic": float(lr.test_statistic),
                                        "n": int(len(target))}
    # pairwise G12R vs G12D (Hayashi 2021 hypothesis)
    g12r = surv[surv["allele_grp"] == "G12R"]
    g12d = surv[surv["allele_grp"] == "G12D"]
    if len(g12r) >= 5 and len(g12d) >= 10:
        from lifelines.statistics import logrank_test
        lr2 = logrank_test(g12r["os_months"], g12d["os_months"],
                            g12r["event"], g12d["event"])
        findings["A_logrank_G12R_vs_G12D"] = {
            "n_G12R": int(len(g12r)), "n_G12D": int(len(g12d)),
            "median_OS_G12R": medians.get("G12R"),
            "median_OS_G12D": medians.get("G12D"),
            "p_value": float(lr2.p_value), "test_statistic": float(lr2.test_statistic)}

    # Cox: allele + Moffitt + age + stage_simplified
    surv_cox = surv.copy()
    surv_cox["basal"] = (surv_cox["moffitt_call"] == "basal-like").astype(int)
    surv_cox["G12R"] = (surv_cox["allele_grp"] == "G12R").astype(int)
    surv_cox["G12D"] = (surv_cox["allele_grp"] == "G12D").astype(int)
    surv_cox["G12V"] = (surv_cox["allele_grp"] == "G12V").astype(int)
    cox_df = surv_cox[["os_months","event","basal","G12R","G12D","G12V","age"]].dropna()
    cox_df = cox_df[cox_df["os_months"] > 0]
    if len(cox_df) >= 30:
        cph = CoxPHFitter().fit(cox_df, duration_col="os_months", event_col="event")
        s = cph.summary[["coef","exp(coef)","p","exp(coef) lower 95%","exp(coef) upper 95%"]].round(4)
        findings["A_cox_table"] = s.reset_index().to_dict(orient="records")
        findings["A_cox_concordance"] = float(cph.concordance_index_)

    # ===== (B) Allele × Moffitt enrichment =====
    contingency = pd.crosstab(df["kras_allele"], df["moffitt_call"])
    contingency.to_csv(RES / "kras_x_moffitt_contingency.tsv", sep="\t")
    findings["B_contingency"] = contingency.to_dict()
    if contingency.shape[0] >= 2 and contingency.shape[1] == 2:
        chi2, p, dof, _ = chi2_contingency(contingency.values + 0.5)
        findings["B_chi2"] = {"chi2": round(float(chi2), 3),
                              "p_value": round(float(p), 5),
                              "dof": int(dof)}
    # G12R basal fraction
    for al in ["G12D", "G12V", "G12R", "G12C"]:
        sub = df[df["kras_allele"] == al]
        if len(sub) >= 5:
            findings[f"B_basal_fraction_{al}"] = round(
                (sub["moffitt_call"] == "basal-like").mean(), 3)

    # ===== (C) Allele × TME module Cohen's d =====
    base = df[df["kras_allele"].isin(["G12_other", "KRAS_other", "WT"])]
    tme = {}
    for al in ["G12D", "G12V", "G12R", "G12C"]:
        sub = df[df["kras_allele"] == al]
        if len(sub) < 5:
            continue
        d = {}
        for m in mod_cols:
            if m in df.columns:
                d[m] = round(cohens_d(sub[m].dropna(), base[m].dropna()), 3)
        tme[al] = d
    findings["C_module_d_vs_baseline"] = tme

    # ===== (D) Allele × NeoQ + Korean HLA targetability =====
    # Korean HLA frequencies for the dominant allele restrictions per KRAS
    KRAS_HLA = {
        "G12D": {"alleles": ["A*11:01", "A*03:01", "B*07:02", "DRB1*04"],
                 "korean_combined_pct": 31.5, "european_combined_pct": 26.4},
        "G12V": {"alleles": ["A*11:01", "B*07:02", "DRB1*0701"],
                 "korean_combined_pct": 28.1, "european_combined_pct": 23.0},
        "G12R": {"alleles": ["B*07:02", "DRB1*04:01"],
                 "korean_combined_pct": 8.9, "european_combined_pct": 16.6},
        "G12C": {"alleles": ["A*02:01", "B*07:02"],
                 "korean_combined_pct": 22.0, "european_combined_pct": 38.0},
    }
    neo = {}
    for al, info in KRAS_HLA.items():
        sub = df[df["kras_allele"] == al]
        if len(sub) >= 5:
            n_top = sub["top_NeoQ"].dropna()
            neo[al] = {
                "n": int(len(sub)),
                "median_top_NeoQ": round(float(n_top.median()), 4) if len(n_top) else None,
                "korean_HLA_combined_pct": info["korean_combined_pct"],
                "european_HLA_combined_pct": info["european_combined_pct"],
                "tcga_prevalence_pct": round(len(sub)/len(df)*100, 2),
                # joint prob: carries this mutation × carries any restricted allele
                "korean_joint_actionability_pct": round(
                    (len(sub)/len(df)) * (info["korean_combined_pct"]/100) * 100, 2),
                "european_joint_actionability_pct": round(
                    (len(sub)/len(df)) * (info["european_combined_pct"]/100) * 100, 2),
            }
    findings["D_allele_neoQ_HLA"] = neo

    # ===== (E) Moffitt × TME PARADOX =====
    # basal-like IS simultaneously T-cell-inflamed (TLS, IFNG, HLA-II up) AND
    # myeloid-suppressed (myeloid_suppressive up, checkpoint exhaustion up)
    # Quantify with paired-axis coordinate per sample.
    if all(c in df.columns for c in ["TLS_CXCL13","IFNG_inflamed","myeloid_suppressive","checkpoint_exhaustion","HLA_II"]):
        df["score_inflamed"]  = df[["TLS_CXCL13","IFNG_inflamed","HLA_II"]].mean(axis=1)
        df["score_suppress"]  = df[["myeloid_suppressive","checkpoint_exhaustion"]].mean(axis=1)
        df["paradox_score"]   = df["score_inflamed"] + df["score_suppress"]
        b = df.query("moffitt_call=='basal-like'")
        c = df.query("moffitt_call=='classical'")
        findings["E_paradox"] = {
            "basal_inflamed_d_vs_classical": round(cohens_d(b["score_inflamed"].dropna(), c["score_inflamed"].dropna()), 3),
            "basal_suppress_d_vs_classical": round(cohens_d(b["score_suppress"].dropna(), c["score_suppress"].dropna()), 3),
            "basal_paradox_d_vs_classical":  round(cohens_d(b["paradox_score"].dropna(), c["paradox_score"].dropna()), 3),
            "fraction_basal_high_BOTH_inflamed_AND_suppressed":
                round(float(((b["score_inflamed"] > 0) & (b["score_suppress"] > 0)).mean()), 3),
            "fraction_classical_high_BOTH_inflamed_AND_suppressed":
                round(float(((c["score_inflamed"] > 0) & (c["score_suppress"] > 0)).mean()), 3),
        }

    # ===== Save findings =====
    with open(RES / "kras_allele_finding.json", "w") as f:
        json.dump(findings, f, indent=2, default=str)

    # ===== Figures =====
    # Fig A: KM by allele (G12D / G12V / G12R)
    fig, ax = plt.subplots(figsize=(6.5, 4.4), dpi=140)
    cmap = {"G12D": ROSE, "G12V": AMBER, "G12R": EMERALD, "G12C": VIOLET, "WT": "#888"}
    for al in ["G12D", "G12V", "G12R", "G12C"]:
        s = surv[surv["allele_grp"] == al]
        if len(s) < 5:
            continue
        kmf = KaplanMeierFitter().fit(s["os_months"], s["event"],
                                       label=f"{al} (n={len(s)})")
        kmf.plot_survival_function(ax=ax, ci_show=False,
                                    color=cmap.get(al, "#333"), linewidth=2.2)
    ax.set_xlabel("Months from diagnosis"); ax.set_ylabel("Overall survival")
    ax.set_title(f"TCGA-PAAD · OS by KRAS allele (Hayashi 2021 NatCancer hypothesis)")
    ax.set_xlim(0, 80)
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    p_str = ""
    if "A_logrank_G12R_vs_G12D" in findings:
        p_str += f"  G12R vs G12D logrank p={findings['A_logrank_G12R_vs_G12D']['p_value']:.3f}"
    if "A_logrank_G12DVR" in findings:
        p_str += f"  ·  3-way p={findings['A_logrank_G12DVR']['p_value']:.3f}"
    if p_str:
        ax.text(0.02, 0.04, p_str.strip(), transform=ax.transAxes,
                fontsize=8.5, color="#333",
                bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#bbb"))
    fig.tight_layout()
    fig.savefig(FIG / "fig_kras_allele_KM.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIG / "fig_kras_allele_KM.svg", bbox_inches="tight")
    plt.close(fig)

    # Fig B: allele × Moffitt stacked bar
    ct = pd.crosstab(df["kras_allele"], df["moffitt_call"])
    ct = ct.reindex([a for a in ["G12D","G12V","G12R","G12C","G12_other","Q61","KRAS_other","WT"] if a in ct.index])
    pct = ct.div(ct.sum(axis=1), axis=0) * 100
    fig, ax = plt.subplots(figsize=(7.5, 4.0), dpi=140)
    pct.plot(kind="barh", stacked=True, color=[ROSE, EMERALD], ax=ax,
             edgecolor="white", linewidth=0.6)
    for i, (idx, n) in enumerate(zip(ct.index, ct.sum(axis=1))):
        ax.text(102, i, f"n={int(n)}", va="center", fontsize=8.5, color="#444")
    ax.set_xlabel("% of allele class"); ax.set_xlim(0, 116)
    ax.set_ylabel("KRAS allele")
    ax.set_title("KRAS allele × Moffitt subtype enrichment (TCGA-PAAD n="+str(int(ct.values.sum()))+")")
    ax.legend(loc="lower right", frameon=False, fontsize=9, title="Moffitt")
    fig.tight_layout()
    fig.savefig(FIG / "fig_kras_allele_moffitt.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIG / "fig_kras_allele_moffitt.svg", bbox_inches="tight")
    plt.close(fig)

    # Fig C: allele × TME module heatmap (Cohen's d)
    keys = list(tme.keys())
    if keys:
        rows = sorted({m for v in tme.values() for m in v.keys()})
        arr = np.array([[tme[al].get(m, np.nan) for al in keys] for m in rows])
        fig, ax = plt.subplots(figsize=(5.5, 4.6), dpi=140)
        im = ax.imshow(arr, cmap="RdBu_r", vmin=-1.0, vmax=1.0, aspect="auto")
        ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows, fontsize=9)
        ax.set_xticks(range(len(keys))); ax.set_xticklabels(keys, fontsize=10, rotation=0)
        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                v = arr[i, j]
                if np.isfinite(v):
                    ax.text(j, i, f"{v:+.2f}", ha="center", va="center",
                            fontsize=8, color="#0a0a0a")
        ax.set_title("KRAS allele · TME module Cohen's d (vs other-KRAS / WT baseline)")
        fig.colorbar(im, ax=ax, fraction=0.04, pad=0.04, label="Cohen's d")
        fig.tight_layout()
        fig.savefig(FIG / "fig_kras_allele_TME.png", dpi=160, bbox_inches="tight")
        fig.savefig(FIG / "fig_kras_allele_TME.svg", bbox_inches="tight")
        plt.close(fig)

    # Fig D: paradox 2D map (inflamed vs suppress, colored by Moffitt)
    if "score_inflamed" in df.columns:
        fig, ax = plt.subplots(figsize=(6.0, 5.0), dpi=140)
        for cls, col in [("basal-like", ROSE), ("classical", EMERALD)]:
            sub = df[df["moffitt_call"] == cls]
            ax.scatter(sub["score_inflamed"], sub["score_suppress"],
                       s=22, alpha=0.7, c=col, edgecolors="white", linewidths=0.5,
                       label=f"{cls} (n={len(sub)})")
        ax.axhline(0, color="#bbb", lw=0.6); ax.axvline(0, color="#bbb", lw=0.6)
        ax.set_xlabel("Inflamed score (TLS + IFNG + HLA-II)")
        ax.set_ylabel("Suppress score (myeloid + checkpoint)")
        ax.set_title("PDAC paradox · basal-like = inflamed AND suppressed (TCGA-PAAD)")
        ax.legend(loc="lower right", frameon=False, fontsize=9)
        # annotate quadrants
        x_max, y_max = df["score_inflamed"].max(), df["score_suppress"].max()
        ax.text(x_max*0.9, y_max*0.95, "vaccine + ICI\n+ myeloid-modifier\nrationale",
                ha="right", va="top", fontsize=8.5, color="#7a1e2e",
                bbox=dict(boxstyle="round,pad=0.4", fc=(1,1,1,.85), ec="#e9a3b1"))
        fig.tight_layout()
        fig.savefig(FIG / "fig_paradox.png", dpi=160, bbox_inches="tight")
        fig.savefig(FIG / "fig_paradox.svg", bbox_inches="tight")
        plt.close(fig)

    # ===== Print summary =====
    print("[12] FINDINGS SUMMARY")
    print(f"  (A) median OS per allele (months):")
    for k, v in findings.get("A_median_OS_months", {}).items():
        print(f"      {k:11s}  {v}")
    if "A_logrank_G12R_vs_G12D" in findings:
        x = findings["A_logrank_G12R_vs_G12D"]
        print(f"      G12R vs G12D logrank p = {x['p_value']:.3f}  "
              f"(n_G12R={x['n_G12R']}, n_G12D={x['n_G12D']})")
    if "A_cox_concordance" in findings:
        print(f"      Cox concordance C = {findings['A_cox_concordance']:.3f}")
    print(f"  (B) chi2 KRAS_allele × Moffitt: ", findings.get("B_chi2"))
    for al in ["G12D","G12V","G12R","G12C"]:
        k = f"B_basal_fraction_{al}"
        if k in findings: print(f"      basal fraction in {al}: {findings[k]:.0%}")
    print(f"  (C) per-allele TME Cohen's d:")
    for al, m in findings["C_module_d_vs_baseline"].items():
        sig = ", ".join(f"{k}={v:+.2f}" for k,v in m.items() if abs(v) >= 0.3)
        print(f"      {al:5s}  {sig if sig else '— small effects —'}")
    print(f"  (D) Korean joint actionability per allele:")
    for al, m in findings["D_allele_neoQ_HLA"].items():
        print(f"      {al:5s}  joint Korean: {m['korean_joint_actionability_pct']}%   "
              f"joint EU: {m['european_joint_actionability_pct']}%   "
              f"prev: {m['tcga_prevalence_pct']}%")
    if "E_paradox" in findings:
        print(f"  (E) Moffitt × TME paradox:")
        for k, v in findings["E_paradox"].items():
            print(f"      {k:55s}  {v}")
    print(f"\n[12] artifacts → results/kras_allele_finding.json + 4 figures")


if __name__ == "__main__":
    main()
