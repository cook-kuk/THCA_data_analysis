#!/usr/bin/env python3
"""Phase C v2 — per-cohort + cross-cohort statistics + figures.

For each score (DM1_inflam_composite, lineage_portable_DM1, cytolytic, IFNG, HLA_I, HLA_II, TLS):
  - Per cohort: Wilcoxon (R vs NR), AUC, logistic regression slope (within-cohort z-standardized score)
  - Across cohorts: forest plot of log OR (logistic) and standardized beta
  - Survival: KM by score-tertile (only IMvigor + GSE176307 have OS)
  - Sign consistency table

Figures (under results/figures/):
  cohort_qc_barplot_sample_counts.pdf
  score_by_response_each_cohort.pdf
  forest_logOR_response.pdf
  survival_KM_if_available.pdf
  signature_correlation_heatmap.pdf

Tables (under results/tables/):
  per_cohort_response_stats.tsv
  forest_logOR_table.tsv
  survival_cox_per_cohort.tsv
  final_phase_C_ICI_summary_table.tsv
"""
from __future__ import annotations
from pathlib import Path
import json
import pandas as pd
import numpy as np
from scipy import stats as st

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/data/thca/repo_results/paper11_pancancer/phase_C_ICI")
TABLES = ROOT / "results" / "tables"
FIGS = ROOT / "results" / "figures"
TABLES.mkdir(parents=True, exist_ok=True)
FIGS.mkdir(parents=True, exist_ok=True)

SIGS_PATH = TABLES / "signature_scores_per_sample.tsv"

SCORES = ["DM1_inflam_composite", "lineage_portable_DM1", "cytolytic_GZMA_PRF1",
          "IFNG_T_cell_inflamed", "HLA_class_I", "HLA_class_II",
          "TLS_CXCL13_like", "myeloid_suppressive", "checkpoint_exhaustion"]

PALETTE = {"IMvigor210": "#1f77b4", "GSE176307": "#17becf",
           "riaz_GSE91061": "#d62728", "MGH_GSE115821": "#e377c2"}


# ---------- helpers
def cohens_d(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2: return np.nan
    s = np.sqrt((np.var(a, ddof=1) + np.var(b, ddof=1)) / 2 + 1e-12)
    return (np.mean(a) - np.mean(b)) / s


def auc_safe(y, score):
    try:
        from sklearn.metrics import roc_auc_score
        y = np.asarray(y); score = np.asarray(score)
        m = ~(pd.isna(y) | pd.isna(score))
        if m.sum() < 5 or len(set(y[m])) < 2: return np.nan
        return float(roc_auc_score(y[m], score[m]))
    except Exception:
        return np.nan


def logistic_slope(score: pd.Series, y: pd.Series):
    """Logistic regression: y ~ z(score). Returns (log_or, se, p)."""
    try:
        from sklearn.linear_model import LogisticRegression
        m = score.notna() & y.notna()
        if m.sum() < 8 or y[m].nunique() < 2: return (np.nan, np.nan, np.nan)
        z = (score[m] - score[m].mean()) / score[m].std(ddof=1)
        clf = LogisticRegression(C=1e6, solver="lbfgs", max_iter=2000)
        clf.fit(z.values.reshape(-1,1), y[m].astype(int).values)
        b = float(clf.coef_[0,0])
        # SE via Wald: 1/var of fitted log-odds
        # use statsmodels for proper SE
        import statsmodels.api as sm
        X = sm.add_constant(z.values)
        res = sm.Logit(y[m].astype(int).values, X).fit(disp=False, method="newton")
        b2 = float(res.params[1])
        se = float(res.bse[1])
        p = float(res.pvalues[1])
        return (b2, se, p)
    except Exception as e:
        return (np.nan, np.nan, np.nan)


# ---------- load
sigs = pd.read_csv(SIGS_PATH, sep="\t")
print(f"loaded {sigs.shape} from {SIGS_PATH}")

# baseline (pre-treatment only): IMvigor + GSE176307 are pre by design;
# Riaz keep only Pre; MGH split by treatment_state — Pre only for cleanest comparison.
def is_pre(row):
    tp = str(row.get("timepoint", "")).lower()
    return tp in {"pre", ""}  # IMvigor + GSE176307 default to "pre"; missing→treated as pre


sigs["use_for_stats"] = sigs.apply(is_pre, axis=1)


# ---------- per-cohort stats
rows = []
for cohort, sub in sigs.groupby("cohort"):
    sub_p = sub[sub["use_for_stats"]].copy()
    n_total = len(sub_p)
    n_resp = sub_p["response_binary_CRPR_vs_SD_PD"].notna().sum()
    n_R = int((sub_p["response_binary_CRPR_vs_SD_PD"] == 1).sum())
    n_NR = int((sub_p["response_binary_CRPR_vs_SD_PD"] == 0).sum())

    for score_col in SCORES:
        if score_col not in sub_p.columns: continue
        s = sub_p[score_col]
        y = sub_p["response_binary_CRPR_vs_SD_PD"]
        d = cohens_d(s[y == 1], s[y == 0])
        auc = auc_safe(y, s)
        try:
            tw = st.ttest_ind(s[y == 1].dropna(), s[y == 0].dropna(), equal_var=False)
            t_p = float(tw.pvalue)
        except Exception: t_p = np.nan
        try:
            wu = st.mannwhitneyu(s[y == 1].dropna(), s[y == 0].dropna(), alternative="two-sided")
            wu_p = float(wu.pvalue)
        except Exception: wu_p = np.nan
        log_or, se, p_log = logistic_slope(s, y)

        rows.append({
            "cohort": cohort, "cancer_type": sub_p["cancer_type"].mode().iloc[0] if len(sub_p) else "",
            "score": score_col, "n_total": n_total, "n_with_resp": int(n_resp),
            "n_R": n_R, "n_NR": n_NR,
            "cohens_d": d, "auc": auc, "t_p": t_p, "wilcoxon_p": wu_p,
            "logOR_per_z": log_or, "logOR_se": se, "logOR_p": p_log,
        })
per_cohort = pd.DataFrame(rows)
per_cohort.to_csv(TABLES / "per_cohort_response_stats.tsv", sep="\t", index=False)
print(f"per_cohort stats → {TABLES/'per_cohort_response_stats.tsv'} ({per_cohort.shape})")


# ---------- meta-analysis (fixed-effect inverse-variance) of log OR per score
def fixed_effect_meta(df: pd.DataFrame):
    sub = df.dropna(subset=["logOR_per_z", "logOR_se"])
    if len(sub) < 2: return np.nan, np.nan, np.nan, len(sub)
    w = 1 / (sub["logOR_se"] ** 2)
    pooled = (sub["logOR_per_z"] * w).sum() / w.sum()
    se_pool = np.sqrt(1 / w.sum())
    z = pooled / se_pool
    p = 2 * (1 - st.norm.cdf(abs(z)))
    return pooled, se_pool, p, len(sub)


meta_rows = []
for score_col in SCORES:
    sub = per_cohort[per_cohort.score == score_col]
    pooled, se_pool, p, k = fixed_effect_meta(sub)
    sign_n = (sub["logOR_per_z"] > 0).sum()
    meta_rows.append({"score": score_col, "k_cohorts": k,
                      "pooled_logOR_per_z": pooled, "pooled_se": se_pool,
                      "pooled_OR": np.exp(pooled) if not pd.isna(pooled) else np.nan,
                      "pooled_p": p, "n_positive_direction": int(sign_n),
                      "n_total_cohorts": int(len(sub))})
meta = pd.DataFrame(meta_rows)
meta.to_csv(TABLES / "forest_logOR_table.tsv", sep="\t", index=False)
print(f"meta → {TABLES/'forest_logOR_table.tsv'}")


# ---------- survival (Cox) for cohorts with OS
cox_rows = []
try:
    from lifelines import CoxPHFitter
    for cohort, sub in sigs.groupby("cohort"):
        sub_p = sub[sub["use_for_stats"]].copy()
        if "os_time" not in sub_p.columns or "os_event" not in sub_p.columns: continue
        os_t = pd.to_numeric(sub_p["os_time"], errors="coerce")
        os_e = pd.to_numeric(sub_p["os_event"], errors="coerce")
        for score_col in SCORES:
            if score_col not in sub_p.columns: continue
            s = sub_p[score_col]
            df = pd.DataFrame({"t": os_t, "e": os_e, "z": (s - s.mean())/s.std(ddof=1)}).dropna()
            if len(df) < 15 or df["e"].sum() < 5: continue
            try:
                cox = CoxPHFitter()
                cox.fit(df, "t", "e")
                hr = float(np.exp(cox.params_["z"]))
                lo = float(np.exp(cox.confidence_intervals_.loc["z","95% lower-bound"]))
                hi = float(np.exp(cox.confidence_intervals_.loc["z","95% upper-bound"]))
                p = float(cox.summary.loc["z", "p"])
                cox_rows.append({"cohort": cohort, "score": score_col, "n": len(df),
                                 "n_event": int(df["e"].sum()), "HR_per_z": hr,
                                 "HR_lo": lo, "HR_hi": hi, "p": p})
            except Exception:
                pass
except ImportError:
    print("[warn] lifelines not available — skipping Cox")
cox_df = pd.DataFrame(cox_rows)
cox_df.to_csv(TABLES / "survival_cox_per_cohort.tsv", sep="\t", index=False)
print(f"cox → {TABLES/'survival_cox_per_cohort.tsv'} ({cox_df.shape})")


# ---------- final summary table
final_rows = []
for score_col in SCORES:
    sub = per_cohort[per_cohort.score == score_col]
    m = meta[meta.score == score_col].iloc[0] if (meta.score == score_col).any() else None
    cox_sub = cox_df[cox_df.score == score_col] if len(cox_df) else pd.DataFrame()
    sign_resp = sum(sub["cohens_d"] > 0)
    final_rows.append({
        "score": score_col,
        "k_cohorts_with_response": int(sub["n_with_resp"].gt(0).sum()),
        "n_total_R+NR": int(sub["n_R"].sum() + sub["n_NR"].sum()),
        "n_cohorts_d_positive": int(sign_resp),
        "mean_cohens_d": round(float(sub["cohens_d"].mean()), 3),
        "median_AUC": round(float(sub["auc"].median()), 3),
        "pooled_OR_per_z": round(float(m["pooled_OR"]), 3) if m is not None and not pd.isna(m["pooled_OR"]) else np.nan,
        "pooled_p": float(m["pooled_p"]) if m is not None else np.nan,
        "k_cohorts_with_OS": int(len(cox_sub)),
        "n_cohorts_HR_lt_1": int((cox_sub["HR_per_z"] < 1).sum()) if len(cox_sub) else 0,
        "median_HR_OS": round(float(cox_sub["HR_per_z"].median()), 3) if len(cox_sub) else np.nan,
    })
final = pd.DataFrame(final_rows)
final.to_csv(TABLES / "final_phase_C_ICI_summary_table.tsv", sep="\t", index=False)
print(f"final → {TABLES/'final_phase_C_ICI_summary_table.tsv'}")
print(final.to_string(index=False))


# ====================================================================== FIGURES

# Fig 1 — sample count QC barplot
counts_qc = sigs.groupby("cohort").size().to_dict()
resp_qc = sigs.dropna(subset=["response_binary_CRPR_vs_SD_PD"]).groupby("cohort").size().to_dict()
pre_qc = sigs[sigs["use_for_stats"]].groupby("cohort").size().to_dict()
fig, ax = plt.subplots(figsize=(7, 4))
labels = list(counts_qc.keys())
x = np.arange(len(labels))
w = 0.27
ax.bar(x - w, [counts_qc[c] for c in labels], w, label="Total samples", color="lightgrey", edgecolor="k")
ax.bar(x, [pre_qc.get(c, 0) for c in labels], w, label="Pre-treatment used", color=[PALETTE[c] for c in labels])
ax.bar(x + w, [resp_qc.get(c, 0) for c in labels], w, label="With binary response", color="white", edgecolor="k", hatch="///")
ax.set_xticks(x); ax.set_xticklabels(labels, rotation=15, ha="right")
ax.set_ylabel("samples")
ax.legend(loc="upper left", frameon=False)
ax.set_title("Phase C v2 — cohort QC")
plt.tight_layout()
plt.savefig(FIGS / "cohort_qc_barplot_sample_counts.pdf"); plt.savefig(FIGS / "cohort_qc_barplot_sample_counts.png", dpi=140)
plt.close()


# Fig 2 — score by response, each cohort × DM1 + IFNG + lineage_portable
target_scores = ["DM1_inflam_composite", "IFNG_T_cell_inflamed", "lineage_portable_DM1"]
n_cohorts = len(sigs["cohort"].unique())
fig, axes = plt.subplots(len(target_scores), n_cohorts, figsize=(2.6 * n_cohorts, 2.6 * len(target_scores)),
                         sharey="row")
cohorts_sorted = sorted(sigs["cohort"].unique())
for i, score in enumerate(target_scores):
    for j, cohort in enumerate(cohorts_sorted):
        ax = axes[i, j] if len(target_scores) > 1 else axes[j]
        sub = sigs[(sigs.cohort == cohort) & (sigs["use_for_stats"])].copy()
        sub = sub.dropna(subset=[score, "response_binary_CRPR_vs_SD_PD"])
        if len(sub) < 5:
            ax.text(0.5, 0.5, "n<5", ha="center", va="center", transform=ax.transAxes)
            ax.set_xticks([]); ax.set_yticks([])
            continue
        for k, (lab, mask) in enumerate([("NR", sub.response_binary_CRPR_vs_SD_PD == 0),
                                         ("R",  sub.response_binary_CRPR_vs_SD_PD == 1)]):
            vals = sub.loc[mask, score].values
            ax.boxplot(vals, positions=[k], widths=0.5, patch_artist=True,
                       boxprops=dict(facecolor=("lightgrey" if lab == "NR" else PALETTE[cohort])),
                       medianprops=dict(color="k"))
            ax.scatter(np.repeat(k, len(vals)) + np.random.uniform(-0.1, 0.1, size=len(vals)),
                       vals, s=8, alpha=0.5, color="k")
        # stat
        stat = per_cohort[(per_cohort.cohort == cohort) & (per_cohort.score == score)]
        if len(stat):
            d = stat["cohens_d"].iloc[0]; auc = stat["auc"].iloc[0]; p = stat["wilcoxon_p"].iloc[0]
            ax.set_title(f"{cohort}\nd={d:.2f}, AUC={auc:.2f}, p={p:.2g}", fontsize=8)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["NR", "R"], fontsize=8)
        if j == 0:
            ax.set_ylabel(score.replace("_", " "), fontsize=8)
plt.suptitle("Phase C v2 — score by response per cohort (pre-treatment)", y=1.02)
plt.tight_layout()
plt.savefig(FIGS / "score_by_response_each_cohort.pdf", bbox_inches="tight")
plt.savefig(FIGS / "score_by_response_each_cohort.png", dpi=140, bbox_inches="tight")
plt.close()


# Fig 3 — forest of log OR
fig, ax = plt.subplots(figsize=(8, max(4, 0.4 * (len(SCORES) + 4))))
y = 0
yticks, yticklabels = [], []
for score_col in SCORES:
    sub = per_cohort[per_cohort.score == score_col].dropna(subset=["logOR_per_z","logOR_se"])
    for _, r in sub.iterrows():
        lo = r["logOR_per_z"] - 1.96 * r["logOR_se"]
        hi = r["logOR_per_z"] + 1.96 * r["logOR_se"]
        ax.plot([lo, hi], [y, y], color=PALETTE.get(r["cohort"], "k"), lw=1.5)
        ax.scatter([r["logOR_per_z"]], [y], color=PALETTE.get(r["cohort"], "k"), s=30)
        yticks.append(y); yticklabels.append(f"{score_col} | {r['cohort']} (n_R={int(r['n_R'])}/n_NR={int(r['n_NR'])})")
        y += 1
    # pooled
    m = meta[meta.score == score_col]
    if len(m) and not pd.isna(m.pooled_logOR_per_z.iloc[0]):
        mp, ms = m.pooled_logOR_per_z.iloc[0], m.pooled_se.iloc[0]
        ax.plot([mp - 1.96*ms, mp + 1.96*ms], [y, y], color="black", lw=2.5)
        ax.scatter([mp], [y], color="black", marker="D", s=60)
        yticks.append(y); yticklabels.append(f"{score_col} | POOLED")
        y += 1
    y += 0.5  # gap
ax.axvline(0, color="grey", ls=":", lw=0.7)
ax.set_yticks(yticks); ax.set_yticklabels(yticklabels, fontsize=7)
ax.set_xlabel("log(OR) per +1 z (CR/PR vs SD/PD)")
ax.set_title("Phase C v2 — forest of logistic log(OR) per signature × cohort + pooled")
plt.tight_layout()
plt.savefig(FIGS / "forest_logOR_response.pdf"); plt.savefig(FIGS / "forest_logOR_response.png", dpi=140)
plt.close()


# Fig 4 — survival KM (median split) for IMvigor + GSE176307 by DM1_inflam_composite
try:
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for k, cohort in enumerate(["IMvigor210", "GSE176307"]):
        ax = axes[k]
        sub = sigs[(sigs.cohort == cohort) & sigs["use_for_stats"]].copy()
        sub = sub.dropna(subset=["DM1_inflam_composite", "os_time", "os_event"])
        if len(sub) < 20:
            ax.text(0.5, 0.5, f"{cohort}: n<20", ha="center", va="center", transform=ax.transAxes); continue
        med = sub["DM1_inflam_composite"].median()
        sub["grp"] = np.where(sub["DM1_inflam_composite"] >= med, "high", "low")
        for grp, color in [("high", "tab:red"), ("low", "tab:blue")]:
            ss = sub[sub.grp == grp]
            kmf = KaplanMeierFitter()
            kmf.fit(ss["os_time"], ss["os_event"], label=f"{grp} (n={len(ss)})")
            kmf.plot(ax=ax, ci_show=False, color=color)
        try:
            lr = logrank_test(sub.loc[sub.grp == "high", "os_time"], sub.loc[sub.grp == "low", "os_time"],
                              event_observed_A=sub.loc[sub.grp == "high", "os_event"],
                              event_observed_B=sub.loc[sub.grp == "low", "os_event"])
            p = lr.p_value
        except Exception: p = np.nan
        ax.set_title(f"{cohort} — DM1_inflam median split  (logrank p={p:.3g})", fontsize=10)
        ax.set_xlabel("Time"); ax.set_ylabel("Survival")
    plt.tight_layout()
    plt.savefig(FIGS / "survival_KM_if_available.pdf"); plt.savefig(FIGS / "survival_KM_if_available.png", dpi=140)
    plt.close()
except ImportError:
    print("[warn] lifelines KM skipped")


# Fig 5 — signature correlation heatmap (across all samples, pooled within cohort z)
sig_z = sigs.copy()
for c, sub in sig_z.groupby("cohort"):
    for s in SCORES:
        if s in sig_z.columns:
            v = sub[s]
            sig_z.loc[sub.index, s] = (v - v.mean()) / v.std(ddof=1)
corr = sig_z[SCORES].corr(method="spearman")
fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(corr.values, cmap="RdBu_r", vmin=-1, vmax=1)
ax.set_xticks(range(len(SCORES))); ax.set_xticklabels(SCORES, rotation=45, ha="right", fontsize=8)
ax.set_yticks(range(len(SCORES))); ax.set_yticklabels(SCORES, fontsize=8)
for i in range(len(SCORES)):
    for j in range(len(SCORES)):
        ax.text(j, i, f"{corr.values[i,j]:.2f}", ha="center", va="center",
                color="white" if abs(corr.values[i,j]) > 0.5 else "black", fontsize=7)
plt.colorbar(im, fraction=0.046, label="Spearman ρ")
plt.title("Signature correlation (within-cohort z, pooled)")
plt.tight_layout()
plt.savefig(FIGS / "signature_correlation_heatmap.pdf"); plt.savefig(FIGS / "signature_correlation_heatmap.png", dpi=140)
plt.close()


# ====================================================================== verdict
final["go_no_go"] = ""
def classify_row(r):
    pos = r["pooled_OR_per_z"]
    p = r["pooled_p"]
    sign = r["n_cohorts_d_positive"]
    k = r["k_cohorts_with_response"]
    if k >= 3 and not pd.isna(pos) and pos > 1.1 and not pd.isna(p) and p < 0.05 and sign >= 3:
        return "GREEN"
    if k >= 2 and ((not pd.isna(pos) and pos > 1) or sign >= 2):
        return "YELLOW"
    return "RED"
final["go_no_go"] = final.apply(classify_row, axis=1)
final.to_csv(TABLES / "final_phase_C_ICI_summary_table.tsv", sep="\t", index=False)


# overall verdict
green = (final.go_no_go == "GREEN").sum()
yellow = (final.go_no_go == "YELLOW").sum()
print(f"\n=== verdict per signature ===")
print(final[["score","k_cohorts_with_response","n_cohorts_d_positive","pooled_OR_per_z","pooled_p","go_no_go"]].to_string(index=False))
overall = "GREEN" if green >= 2 else ("YELLOW" if yellow + green >= 3 else "RED")
print(f"\noverall verdict: {overall}  (green={green}, yellow={yellow})")

# write json summary
summary = {
    "n_cohorts_total": int(sigs["cohort"].nunique()),
    "n_samples_total": int(len(sigs)),
    "n_pre_treatment_used": int(sigs["use_for_stats"].sum()),
    "n_with_binary_response": int(sigs.dropna(subset=["response_binary_CRPR_vs_SD_PD"])
                                    [sigs["use_for_stats"]].shape[0]),
    "scores_GREEN": list(final[final.go_no_go == "GREEN"]["score"]),
    "scores_YELLOW": list(final[final.go_no_go == "YELLOW"]["score"]),
    "scores_RED": list(final[final.go_no_go == "RED"]["score"]),
    "overall_verdict": overall,
}
(TABLES / "phase_C_ICI_overall_summary.json").write_text(json.dumps(summary, indent=2, default=str))
print(f"\nfinal artifacts in {ROOT}/results/")
