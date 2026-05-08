#!/usr/bin/env python3
"""DM1 Round-3 deep dive — TERT × DM × OS, ROC AUC, histology, meta-analysis,
phospho, BCR/TLS — all Paper 1 paper-blocking, no Track B.

Pulls existing artifacts:
  - 36 TCGA-THCA TERT promoter mutations from
    /data/thca/repo_results/v17_tert_recovery/v2/parsed/S6_cbioportal_all_promoter_mutations.tsv
  - sample_master_v17_tert.tsv (histology, driver, OS)
  - tcga_with_clinical_mutations.tsv (DM call)
  - per_sample_panel.tsv (g8 scores from Round 2)
  - phospho_target_genes_S1D.tsv (Mun 2025 phospho)
  - r6_1_TLS_BCR.tsv (TLS+BCR Cohen's d)

Outputs in project/results/dm1_robustness_v2026_05_08/round3/:
  - tert_dm_os_cox.tsv + tert_dm_os.{png,pdf}
  - roc_auc_per_cohort.tsv + roc_auc.{png,pdf}
  - histology_dm_enrichment.tsv + histology_dm.{png,pdf}
  - meta_analysis_pooled.{png,pdf} + meta_results.tsv
  - phospho_layer_summary.tsv + phospho_layer.{png,pdf}
  - bcr_tls_summary.tsv
  - round3_summary.json
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
from scipy import stats
from sklearn.metrics import roc_auc_score
from lifelines import CoxPHFitter, KaplanMeierFitter

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "project" / "results" / "dm1_robustness_v2026_05_08" / "round3"
OUT.mkdir(parents=True, exist_ok=True)

TERT_MUT = Path("/data/thca/repo_results/v17_tert_recovery/v2/parsed/S6_cbioportal_all_promoter_mutations.tsv")
SAMPLE_MASTER = Path("/data/thca/repo_results/v17/tables/sample_master_v17_tert.tsv")
DM_SIG = REPO / "project" / "results" / "d4p2_tcga_hashimoto_signature" / "tcga_signature_scores.tsv"
PER_SAMPLE = REPO / "project" / "results" / "dm1_robustness_v2026_05_08" / "per_sample_panel.tsv"
PHOSPHO = Path("/data/thca/repo_results/proteogenomic_v1/paper3_mun2025_dediff_layer/phospho_target_genes_S1D.tsv")
TLS_BCR = Path("/data/thca/repo_results/audit_2026_04_30/round6/r6_1_TLS_BCR.tsv")
EFFECT_SIZES = REPO / "project" / "results" / "dm1_robustness_v2026_05_08" / "effect_sizes.tsv"
SURVIVAL_RAW = Path("/data/thca/repo_data/raw/TCGA_pancan/survival.tsv")


def short_id(s: str) -> str:
    return "-".join(str(s).split("-")[:4])[:15] if pd.notna(s) else s


# ============== 1. TERT × DM × OS ==============
def tert_dm_os():
    print("\n[1] TERT × DM × OS Cox")
    tert_df = pd.read_csv(TERT_MUT, sep="\t")
    tcga_tert_samples = (
        tert_df.loc[tert_df["study_id"].astype(str).str.contains("thca_tcga", case=False, na=False), "sample_barcode"]
        .astype(str).str.strip().unique().tolist()
    )
    tcga_tert_samples = set(s.upper()[:15] for s in tcga_tert_samples)
    print(f"  recovered TERT+ TCGA samples: {len(tcga_tert_samples)}")

    sig = pd.read_csv(DM_SIG, sep="\t", index_col=0)
    sig["sample_short"] = sig.index.to_series().apply(lambda s: str(s).upper()[:15])
    sig["tert_pos"] = sig["sample_short"].isin(tcga_tert_samples).astype(int)
    print(f"  TERT+ overlap with DM-called samples: {int(sig['tert_pos'].sum())}")

    surv = pd.read_csv(SURVIVAL_RAW, sep="\t")
    surv["sample_short"] = surv["sample"].apply(lambda s: str(s).upper()[:15])
    df = sig.merge(surv[["sample_short", "OS", "OS.time", "DSS", "DSS.time", "PFI", "PFI.time"]],
                   on="sample_short", how="left")
    df = df.dropna(subset=["DM", "OS.time"]).copy()

    rows = []
    for outcome in ("OS", "DSS", "PFI"):
        T = df[f"{outcome}.time"].astype(float)
        E = df[outcome].astype(float)
        for cov_label, cov in [("TERT_promoter", df["tert_pos"]),
                                ("DM1_indicator", (df["DM"] == "DM1").astype(int)),
                                ("TERT_or_DM1", ((df["tert_pos"] == 1) | (df["DM"] == "DM1")).astype(int))]:
            sub = pd.DataFrame({"X": cov, "T": T, "E": E}).dropna()
            if sub["X"].nunique() < 2 or len(sub) < 30 or sub["E"].sum() < 3:
                rows.append({"outcome": outcome, "covariate": cov_label, "n": len(sub),
                             "n_events": int(sub["E"].sum()), "HR": np.nan, "p": np.nan,
                             "note": "underpowered"})
                continue
            try:
                cph = CoxPHFitter().fit(sub, duration_col="T", event_col="E")
                row = cph.summary.iloc[0]
                rows.append({
                    "outcome": outcome, "covariate": cov_label, "n": int(len(sub)),
                    "n_events": int(sub["E"].sum()),
                    "HR": float(row["exp(coef)"]),
                    "HR_lo95": float(row["exp(coef) lower 95%"]),
                    "HR_hi95": float(row["exp(coef) upper 95%"]),
                    "p": float(row["p"]),
                })
            except Exception as e:
                rows.append({"outcome": outcome, "covariate": cov_label, "n": int(len(sub)),
                             "HR": np.nan, "p": np.nan, "note": f"err: {str(e)[:50]}"})

    res = pd.DataFrame(rows)
    res.to_csv(OUT / "tert_dm_os_cox.tsv", sep="\t", index=False)

    # Forest
    plot_df = res.dropna(subset=["HR"]).copy()
    if len(plot_df):
        plot_df = plot_df[(plot_df["HR_lo95"] > 0) & (plot_df["HR_hi95"] < 1e3)]
        fig, ax = plt.subplots(figsize=(11, max(3.5, 0.55 * len(plot_df) + 1)))
        y = np.arange(len(plot_df))
        for i, (_, r) in enumerate(plot_df.iterrows()):
            color = "#d62728" if (r["p"] < 0.05 and r["HR"] > 1) else ("#1f77b4" if (r["p"] < 0.05 and r["HR"] < 1) else "#7f7f7f")
            ax.errorbar(r["HR"], i,
                        xerr=[[max(r["HR"] - r["HR_lo95"], 0)], [max(r["HR_hi95"] - r["HR"], 0)]],
                        fmt="s", color=color, ecolor=color, capsize=4, markersize=7)
        ax.axvline(1.0, color="#888", linestyle=":", linewidth=0.8)
        ax.set_xscale("log")
        ax.set_yticks(y)
        labels = []
        for _, r in plot_df.iterrows():
            p = r["p"]
            p_str = f"p={p:.2g}" if p >= 1e-3 else f"p<1e-{int(-np.log10(p))}"
            labels.append(f"{r['outcome']} | {r['covariate']}  n={r['n']}, ev={r['n_events']}, {p_str}")
        ax.set_yticklabels(labels, fontsize=9)
        ax.invert_yaxis()
        ax.set_xlabel("Hazard ratio (log scale)", fontsize=9)
        ax.set_title("TCGA-THCA TERT × DM × OS/DSS/PFI Cox HR — Round 3 (2026-05-08)\n(36 recovered TERT+ samples per memory v17_tert_recovery_v2)", fontsize=10)
        ax.grid(axis="x", linestyle=":", alpha=0.3)
        plt.tight_layout()
        fig.savefig(OUT / "tert_dm_os.png", dpi=160, bbox_inches="tight")
        fig.savefig(OUT / "tert_dm_os.pdf", bbox_inches="tight")
        plt.close(fig)

    return res


# ============== 2. ROC AUC per cohort ==============
def roc_auc_per_cohort():
    print("\n[2] ROC AUC per cohort")
    ps = pd.read_csv(PER_SAMPLE, sep="\t")
    print(f"  per_sample_panel rows: {len(ps)}")

    # Define discriminations per cohort
    auc_rows = []
    rng = np.random.default_rng(42)

    def boot_auc(y, s, n=2000):
        aucs = []
        idx = np.arange(len(y))
        for _ in range(n):
            ix = rng.choice(idx, size=len(idx), replace=True)
            try:
                aucs.append(roc_auc_score(y[ix], s[ix]))
            except Exception:
                pass
        if not aucs:
            return float("nan"), float("nan"), float("nan")
        return float(np.mean(aucs)), float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))

    discriminations = [
        ("GSE286332", "PTC vs PTC_HT", lambda d: (d["group"] == "PTC_HT").astype(int).to_numpy(),
         lambda d: -d["score"].to_numpy()),  # PTC has higher score, so PTC_HT=1 corresponds to lower score
        ("TCGA-THCA", "DM1 vs DM2", lambda d: (d["group"].str.startswith("DM2")).astype(int).to_numpy(),
         lambda d: -d["score"].to_numpy()),
        ("GSE213647", "PTC vs ATC", lambda d: (d["group"] == "ATC").astype(int).to_numpy(),
         lambda d: -d["score"].to_numpy()),
        ("GSE213647", "ATC vs Normal", lambda d: (d["group"] == "ATC").astype(int).to_numpy(),
         lambda d: -d["score"].to_numpy()),
        ("Mun2025-protein", "PTC vs ATC", lambda d: (d["group"] == "ATC").astype(int).to_numpy(),
         lambda d: -d["score"].to_numpy()),
    ]

    for cohort, label, ymask, scorer in discriminations:
        sub = ps[ps["cohort"] == cohort].copy()
        if cohort == "GSE286332":
            sub = sub[sub["group"].isin(["PTC", "PTC_HT", "PTC+HT"])]
            sub.loc[sub["group"] == "PTC+HT", "group"] = "PTC_HT"
        elif cohort == "TCGA-THCA":
            sub = sub[sub["group"].isin(["DM1", "DM1+hashi", "DM2", "DM2+hashi"])]
        elif cohort == "GSE213647":
            if label == "PTC vs ATC":
                sub = sub[sub["group"].isin(["PTC", "ATC"])]
            elif label == "ATC vs Normal":
                sub = sub[sub["group"].isin(["ATC", "Normal"])]
        elif cohort == "Mun2025-protein":
            sub = sub[sub["group"].isin(["PTC", "ATC"])]

        sub = sub.dropna(subset=["score"])
        if len(sub) < 10 or sub["group"].nunique() < 2:
            auc_rows.append({"cohort": cohort, "contrast": label, "n": len(sub),
                             "auc": np.nan, "ci_lo": np.nan, "ci_hi": np.nan})
            continue

        y = ymask(sub)
        s = scorer(sub)
        try:
            auc = roc_auc_score(y, s)
            mean, lo, hi = boot_auc(y, s)
            auc_rows.append({
                "cohort": cohort, "contrast": label, "n": int(len(sub)),
                "n_pos": int(y.sum()), "n_neg": int((1 - y).sum()),
                "auc": float(auc), "auc_boot_mean": mean,
                "ci_lo": lo, "ci_hi": hi,
            })
        except Exception as e:
            auc_rows.append({"cohort": cohort, "contrast": label, "n": int(len(sub)),
                             "auc": np.nan, "note": str(e)[:50]})

    res = pd.DataFrame(auc_rows)
    res.to_csv(OUT / "roc_auc_per_cohort.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(9, max(3, 0.55 * len(res) + 1.2)))
    plot_df = res.dropna(subset=["auc"])
    y = np.arange(len(plot_df))
    for i, (_, r) in enumerate(plot_df.iterrows()):
        color = "#d62728" if r["auc"] > 0.85 else ("#fdae61" if r["auc"] > 0.7 else "#7f7f7f")
        ax.errorbar(r["auc"], i,
                    xerr=[[r["auc"] - r["ci_lo"]], [r["ci_hi"] - r["auc"]]],
                    fmt="s", color=color, ecolor=color, capsize=4, markersize=7)
    ax.axvline(0.5, color="#888", linestyle=":", linewidth=0.8)
    ax.axvline(0.7, color="#bbb", linestyle=":", linewidth=0.5, alpha=0.5)
    ax.axvline(0.9, color="#bbb", linestyle=":", linewidth=0.5, alpha=0.5)
    ax.set_xlim(0.4, 1.05)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['cohort']}\n  {r['contrast']}  n_pos={r['n_pos']}, n_neg={r['n_neg']}" for _, r in plot_df.iterrows()], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("AUC (8-gene panel as classifier; bootstrap 95% CI)", fontsize=9)
    ax.set_title("8-gene panel ROC AUC per cohort — Round 3 (2026-05-08)\nRed AUC>0.85, Orange 0.7–0.85, Grey <0.7", fontsize=10)
    ax.grid(axis="x", linestyle=":", alpha=0.3)
    plt.tight_layout()
    fig.savefig(OUT / "roc_auc.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "roc_auc.pdf", bbox_inches="tight")
    plt.close(fig)

    return res


# ============== 3. Histology × DM enrichment ==============
def histology_dm():
    print("\n[3] Histology × DM enrichment")
    sm = pd.read_csv(SAMPLE_MASTER, sep="\t", index_col="sample_id")
    sig = pd.read_csv(DM_SIG, sep="\t", index_col=0)
    sm["sample_short"] = sm.index.to_series().apply(short_id)
    sig["sample_short"] = sig.index.to_series().apply(short_id)
    merged = sig.merge(sm[["sample_short", "histology_subtype", "driver_anchor", "ajcc_stage_group"]],
                       on="sample_short", how="left")
    merged = merged.dropna(subset=["DM", "histology_subtype"])

    # 2x2 contingency: cPTC vs FVPTC × DM1 vs DM2
    ct = pd.crosstab(merged["histology_subtype"], merged["DM"])
    ct.to_csv(OUT / "histology_dm_enrichment.tsv", sep="\t")

    # fisher
    fisher_rows = []
    for hist in ct.index:
        if "DM1" in ct.columns and "DM2" in ct.columns:
            a = int(ct.loc[hist, "DM1"]) if "DM1" in ct.columns else 0
            b = int(ct.loc[hist, "DM2"]) if "DM2" in ct.columns else 0
            other_dm1 = int(ct["DM1"].sum() - a)
            other_dm2 = int(ct["DM2"].sum() - b)
            try:
                _, p = stats.fisher_exact([[a, b], [other_dm1, other_dm2]])
                or_, _ = stats.fisher_exact([[a, b], [other_dm1, other_dm2]])
                fisher_rows.append({"histology": hist, "DM1": a, "DM2": b,
                                    "other_DM1": other_dm1, "other_DM2": other_dm2,
                                    "OR": float(or_), "p": float(p)})
            except Exception:
                fisher_rows.append({"histology": hist, "DM1": a, "DM2": b, "OR": np.nan, "p": np.nan})

    fisher = pd.DataFrame(fisher_rows)
    fisher.to_csv(OUT / "histology_dm_fisher.tsv", sep="\t", index=False)

    # Plot stacked bar
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100
    ct_pct.plot(kind="bar", stacked=True, ax=ax, color=["#fdae61", "#abd9e9"], edgecolor="#333")
    for i, hist in enumerate(ct.index):
        total = int(ct.loc[hist].sum())
        ax.text(i, 102, f"n={total}", ha="center", fontsize=9)
    ax.set_ylabel("% of samples")
    ax.set_title("TCGA-THCA DM1 vs DM2 by histology subtype — Round 3 (2026-05-08)", fontsize=10)
    ax.legend(title="DM call", fontsize=9)
    ax.grid(axis="y", linestyle=":", alpha=0.3)
    plt.tight_layout()
    fig.savefig(OUT / "histology_dm.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "histology_dm.pdf", bbox_inches="tight")
    plt.close(fig)

    return ct, fisher


# ============== 4. Random-effects meta-analysis ==============
def meta_analysis():
    print("\n[4] Random-effects meta-analysis (DerSimonian-Laird)")
    es = pd.read_csv(EFFECT_SIZES, sep="\t")
    # Compute SE from CI: SE = (ci_hi - ci_lo) / (2 * 1.96)
    es["se"] = (es["ci_hi"] - es["ci_lo"]) / (2 * 1.96)
    es = es.dropna(subset=["cohens_d", "se"])
    es = es[es["se"] > 0]

    # Fixed effect (inverse-variance)
    w_fe = 1.0 / es["se"] ** 2
    d_fe = (es["cohens_d"] * w_fe).sum() / w_fe.sum()
    se_fe = math.sqrt(1.0 / w_fe.sum())

    # DerSimonian-Laird tau^2
    Q = (w_fe * (es["cohens_d"] - d_fe) ** 2).sum()
    df = len(es) - 1
    c = w_fe.sum() - (w_fe ** 2).sum() / w_fe.sum()
    tau2 = max(0.0, (Q - df) / c) if c > 0 else 0.0
    w_re = 1.0 / (es["se"] ** 2 + tau2)
    d_re = (es["cohens_d"] * w_re).sum() / w_re.sum()
    se_re = math.sqrt(1.0 / w_re.sum())

    res = {
        "fixed_effect_d": float(d_fe), "fixed_effect_se": float(se_fe),
        "fixed_effect_lo95": float(d_fe - 1.96 * se_fe), "fixed_effect_hi95": float(d_fe + 1.96 * se_fe),
        "random_effects_d": float(d_re), "random_effects_se": float(se_re),
        "random_effects_lo95": float(d_re - 1.96 * se_re), "random_effects_hi95": float(d_re + 1.96 * se_re),
        "Q_statistic": float(Q), "df": int(df), "tau_squared": float(tau2),
        "I_squared_pct": float(max(0.0, (Q - df) / Q * 100) if Q > 0 else 0.0),
        "n_studies": int(len(es)),
    }

    with open(OUT / "meta_results.json", "w") as f:
        json.dump(res, f, indent=2)
    pd.DataFrame([res]).to_csv(OUT / "meta_results.tsv", sep="\t", index=False)

    # Plot
    n = len(es)
    fig, ax = plt.subplots(figsize=(11, max(4, 0.4 * n + 2)))
    y = np.arange(n)
    es_sorted = es.sort_values("cohens_d", key=abs, ascending=False).reset_index(drop=True)
    weights = 1.0 / (es_sorted["se"] ** 2 + tau2)
    weights_norm = weights / weights.sum()
    for i, (_, r) in enumerate(es_sorted.iterrows()):
        size = 6 + 30 * weights_norm.iloc[i]
        ax.errorbar(r["cohens_d"], i,
                    xerr=[[r["cohens_d"] - r["ci_lo"]], [r["ci_hi"] - r["cohens_d"]]],
                    fmt="s", color="#222", ecolor="#666", capsize=3, markersize=size, alpha=0.7)
    # Pooled
    ax.errorbar(d_re, n + 0.5, xerr=[[d_re - res["random_effects_lo95"]], [res["random_effects_hi95"] - d_re]],
                fmt="D", color="#d62728", ecolor="#d62728", capsize=8, markersize=14)
    ax.axvline(d_re, color="#d62728", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.axvline(0, color="#888", linestyle=":", linewidth=0.8)
    ax.set_yticks(list(y) + [n + 0.5])
    yticks_labels = [f"{r['cohort']} | {r['contrast']}  n_a={r['n_a']}, n_b={r['n_b']}" for _, r in es_sorted.iterrows()]
    yticks_labels.append(f"POOLED (random-effects DL)\n  d={d_re:+.2f} [{res['random_effects_lo95']:+.2f}, {res['random_effects_hi95']:+.2f}], τ²={tau2:.2f}, I²={res['I_squared_pct']:.0f}%")
    ax.set_yticklabels(yticks_labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Cohen's d", fontsize=9)
    ax.set_title("Random-effects meta-analysis (DerSimonian-Laird) — DM1 8-gene panel, Round 3 (2026-05-08)\nMarker size ∝ inverse-variance weight; red diamond = pooled estimate", fontsize=10)
    ax.grid(axis="x", linestyle=":", alpha=0.3)
    plt.tight_layout()
    fig.savefig(OUT / "meta_analysis_pooled.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "meta_analysis_pooled.pdf", bbox_inches="tight")
    plt.close(fig)

    return res


# ============== 5. Phospho layer ==============
def phospho_layer():
    print("\n[5] Phospho layer (Mun 2025)")
    df = pd.read_csv(PHOSPHO, sep="\t")
    # Sample columns are SAMPLE|GROUP format
    sample_cols = [c for c in df.columns if "|" in c]
    groups = pd.Series([c.split("|")[1].strip() for c in sample_cols], index=sample_cols)
    print(f"  phospho rows: {len(df)} sites; {len(sample_cols)} samples; groups: {groups.unique()}")

    # Per-gene aggregate: mean across phospho sites for that gene, per sample
    agg = df.groupby("gene")[sample_cols].mean()
    # Within-sample z (per gene), then mean
    z = agg.sub(agg.mean(axis=1), axis=0).div(agg.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    panel_score = z.mean(axis=0)
    panel_score.name = "phospho_panel_score"

    # Group-stratified means + Cohen's d
    s_df = panel_score.to_frame()
    s_df["group"] = groups.reindex(s_df.index).values

    rows = []
    for la, lb in [("PTC", "PDTC"), ("PTC", "ATC"), ("PDTC", "ATC")]:
        a = s_df.loc[s_df["group"] == la, "phospho_panel_score"].dropna().to_numpy()
        b = s_df.loc[s_df["group"] == lb, "phospho_panel_score"].dropna().to_numpy()
        if len(a) < 2 or len(b) < 2:
            continue
        sp = math.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
        d = (a.mean() - b.mean()) / sp if sp > 0 else float("nan")
        try:
            p = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
        except Exception:
            p = float("nan")
        rows.append({"contrast": f"{la} vs {lb}", "n_a": int(len(a)), "n_b": int(len(b)),
                     "mean_a": float(a.mean()), "mean_b": float(b.mean()),
                     "cohens_d": float(d), "mwu_p": p})
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "phospho_layer_summary.tsv", sep="\t", index=False)

    # Plot box per group
    fig, ax = plt.subplots(figsize=(7, 4.5))
    order = ["PTC", "PTC_P", "PDTC", "ATC", "ATC_P"]
    plot_groups = [g for g in order if (s_df["group"] == g).any()]
    data = [s_df.loc[s_df["group"] == g, "phospho_panel_score"].dropna().to_numpy() for g in plot_groups]
    bp = ax.boxplot(data, labels=[f"{g}\n(n={len(d)})" for g, d in zip(plot_groups, data)], patch_artist=True, showfliers=False)
    palette = {"PTC": "#2c7fb8", "PTC_P": "#abd9e9", "PDTC": "#fdae61", "ATC": "#d7191c", "ATC_P": "#fc8d59"}
    for patch, g in zip(bp["boxes"], plot_groups):
        patch.set_facecolor(palette.get(g, "#cccccc"))
        patch.set_edgecolor("#222")
    rng = np.random.default_rng(42)
    for i, d in enumerate(data, start=1):
        ax.scatter(rng.normal(i, 0.05, len(d)), d, s=10, color="#222", alpha=0.55, zorder=3)
    ax.set_ylabel("8-gene phospho panel score (z mean)")
    ax.set_title("Mun 2025 phospho — 8-gene panel along PTC → PDTC → ATC dediff axis\n(phospho layer of 4-pillar triangulation per memory proteogenomic_v1_2026_05_08)", fontsize=10)
    ax.grid(axis="y", linestyle=":", alpha=0.3)
    plt.tight_layout()
    fig.savefig(OUT / "phospho_layer.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "phospho_layer.pdf", bbox_inches="tight")
    plt.close(fig)

    return res


# ============== 6. BCR/TLS layer ==============
def bcr_tls():
    print("\n[6] BCR/TLS layer")
    if not TLS_BCR.exists():
        return pd.DataFrame()
    df = pd.read_csv(TLS_BCR, sep="\t")
    df.to_csv(OUT / "bcr_tls_summary.tsv", sep="\t", index=False)
    print(f"  {len(df)} rows from r6_1_TLS_BCR.tsv")
    return df


def main():
    res = {}
    res["tert"] = tert_dm_os().to_dict(orient="records")
    res["roc"] = roc_auc_per_cohort().to_dict(orient="records")
    ct, fisher = histology_dm()
    res["histology_table"] = ct.to_dict()
    res["histology_fisher"] = fisher.to_dict(orient="records")
    res["meta"] = meta_analysis()
    res["phospho"] = phospho_layer().to_dict(orient="records")
    res["bcr_tls_n_rows"] = len(bcr_tls())
    res["generated_at"] = "2026-05-08 round3"

    with open(OUT / "round3_summary.json", "w") as f:
        json.dump(res, f, indent=2, default=str)

    print(f"\nwrote {OUT.relative_to(REPO)}/")
    print("Outputs:")
    for f in sorted(OUT.iterdir()):
        print(f"  {f.relative_to(OUT)}")


if __name__ == "__main__":
    main()
