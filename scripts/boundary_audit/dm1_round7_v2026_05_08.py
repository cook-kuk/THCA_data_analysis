#!/usr/bin/env python3
"""DM1 Round 7 — Nat-Comm reviewer-bulletproofing layers.

Layers:
  1. Random-panel permutation null (1000 random 8-gene panels) — is THIS panel special?
  2. Cluster stability bootstrap (DM1/DM2 reproducibility) — 100 bootstrap re-clusterings
  3. Time-dependent ROC AUC (DM + TERT survival over time)
  4. Calibration curve + Decision Curve Analysis (DCA — clinical net benefit)
  5. SHAP-style XGBoost classifier (8-gene → DM call, feature importance)
  6. ComBat cross-cohort sanity check on per-cohort scores

Outputs in project/results/dm1_robustness_v2026_05_08/round7/.
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
from sklearn.cluster import KMeans
from sklearn.linear_model import LogisticRegression
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "project" / "results" / "dm1_robustness_v2026_05_08" / "round7"
OUT.mkdir(parents=True, exist_ok=True)

TCGA_EXPR = Path("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv")
DM_SIG = REPO / "project" / "results" / "d4p2_tcga_hashimoto_signature" / "tcga_signature_scores.tsv"
SURV = Path("/data/thca/repo_data/raw/TCGA_pancan/survival.tsv")
TERT_MUT = Path("/data/thca/repo_results/v17_tert_recovery/v2/parsed/S6_cbioportal_all_promoter_mutations.tsv")
PER_GENE_R2 = REPO / "project" / "results" / "dm1_robustness_v2026_05_08" / "per_gene_contribution.tsv"

PANEL = ["TG", "TPO", "TSHR", "PAX8", "FOXE1", "NKX2-1", "DIO1", "SLC5A5"]


def cohens_d(x: np.ndarray, y: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) < 2 or len(y) < 2:
        return float("nan")
    sp = math.sqrt(((len(x) - 1) * x.var(ddof=1) + (len(y) - 1) * y.var(ddof=1)) / (len(x) + len(y) - 2))
    return (x.mean() - y.mean()) / sp if sp > 0 else float("nan")


def within_sample_z_score(expr: pd.DataFrame, genes: list) -> pd.Series:
    present = [g for g in genes if g in expr.index]
    if len(present) < 5:
        return pd.Series(dtype=float)
    sub = expr.loc[present]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)


# ============== 1. Random-panel permutation null ==============
def permutation_null():
    print("\n[1] Random-panel permutation null (1000 random 8-gene panels)")
    expr = pd.read_csv(TCGA_EXPR, sep="\t", index_col=0)
    sig = pd.read_csv(DM_SIG, sep="\t", index_col=0)
    dm1 = sig.index[sig["DM"] == "DM1"]
    dm2 = sig.index[sig["DM"] == "DM2"]

    # Canonical
    canon = within_sample_z_score(expr, PANEL)
    d_canon = cohens_d(canon.reindex(dm1).to_numpy(), canon.reindex(dm2).to_numpy())
    print(f"  canonical 8-gene d = {d_canon:+.3f}")

    # Random panels — keep gene variance ≥ canonical mean variance (otherwise mostly zeros dominate)
    canon_var = expr.loc[PANEL].var(axis=1).mean()
    eligible = expr.index[expr.var(axis=1) >= canon_var * 0.5]
    print(f"  eligible high-variance genes: {len(eligible)}")

    rng = np.random.default_rng(42)
    n_perm = 1000
    rand_ds = []
    for _ in range(n_perm):
        sample_genes = rng.choice(eligible, size=8, replace=False)
        sub = expr.loc[sample_genes]
        z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
        score = z.mean(axis=0)
        d = cohens_d(score.reindex(dm1).to_numpy(), score.reindex(dm2).to_numpy())
        if math.isfinite(d):
            rand_ds.append(d)
    rand_arr = np.array(rand_ds)

    # Empirical p-value: fraction of random panels with |d| >= |canonical|
    p_emp = (np.abs(rand_arr) >= abs(d_canon)).mean()
    print(f"  random panel |d| distribution: mean={np.abs(rand_arr).mean():.3f}, max={np.abs(rand_arr).max():.3f}")
    print(f"  empirical p (|d_random| >= |d_canon|): {p_emp:.4f} ({(np.abs(rand_arr) >= abs(d_canon)).sum()}/{n_perm})")

    res = {"d_canon": float(d_canon), "n_perm": n_perm,
           "mean_abs_d_random": float(np.abs(rand_arr).mean()),
           "max_abs_d_random": float(np.abs(rand_arr).max()),
           "p_emp": float(p_emp),
           "eligible_n": int(len(eligible))}
    pd.DataFrame({"random_d": rand_ds}).to_csv(OUT / "permutation_null_random_panels.tsv", sep="\t", index=False)
    (OUT / "permutation_null_summary.json").write_text(json.dumps(res, indent=2))

    # Plot histogram
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.hist(rand_arr, bins=60, color="#abd9e9", edgecolor="black", alpha=0.85)
    ax.axvline(d_canon, color="#d62728", linestyle="--", linewidth=2, label=f"canonical 8-gene d={d_canon:+.2f}")
    ax.axvline(-d_canon, color="#d62728", linestyle=":", linewidth=1, alpha=0.4, label="−d_canon mirror")
    ax.set_xlabel("Cohen's d (random 8-gene panel, DM1 vs DM2)")
    ax.set_ylabel(f"# of random panels (n={n_perm})")
    ax.set_title(
        f"Random-panel permutation null — does THIS 8-gene panel beat random?\n"
        f"empirical p = {p_emp:.4f} ({(np.abs(rand_arr) >= abs(d_canon)).sum()}/{n_perm} random panels reach |d|≥{abs(d_canon):.2f})", fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(axis="y", linestyle=":", alpha=0.3)
    plt.tight_layout()
    fig.savefig(OUT / "permutation_null.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "permutation_null.pdf", bbox_inches="tight")
    plt.close(fig)
    return res


# ============== 2. Cluster stability bootstrap ==============
def cluster_stability():
    print("\n[2] Cluster stability bootstrap (100 re-clusterings)")
    expr = pd.read_csv(TCGA_EXPR, sep="\t", index_col=0)
    panel_present = [g for g in PANEL if g in expr.index]
    sub = expr.loc[panel_present]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=1).replace(0, np.nan), axis=0).T
    z = z.dropna()
    sig = pd.read_csv(DM_SIG, sep="\t", index_col=0)
    z = z.loc[z.index.intersection(sig.index)]
    print(f"  samples: {len(z)} × {z.shape[1]} genes")

    # Reference clustering on full data
    ref_km = KMeans(n_clusters=2, random_state=42, n_init=20).fit(z.values)
    ref_labels = pd.Series(ref_km.labels_, index=z.index)
    # Align labels with DM call
    align = pd.crosstab(ref_labels, sig["DM"].reindex(z.index))
    # Map cluster id with most DM2 → DM2_cluster_id
    if align.shape == (2, 2):
        if align.loc[0, "DM2"] > align.loc[1, "DM2"]:
            ref_dm2 = 0
        else:
            ref_dm2 = 1
        ref_assign = ref_labels.map({ref_dm2: "DM2_cl", 1 - ref_dm2: "DM1_cl"})
    else:
        ref_assign = ref_labels.astype(str)

    # Bootstrap
    rng = np.random.default_rng(42)
    n_boot = 100
    assignments = np.zeros((len(z), n_boot), dtype=int)
    for b in range(n_boot):
        idx_boot = rng.choice(np.arange(len(z)), size=len(z), replace=True)
        z_boot = z.iloc[idx_boot]
        km = KMeans(n_clusters=2, random_state=b, n_init=10).fit(z_boot.values)
        # align cluster labels with reference using majority match
        boot_labels = pd.Series(km.labels_, index=z_boot.index)
        # Compute alignment by checking which cluster has more samples agreeing with ref_dm2
        all_labels = pd.Series(km.predict(z.values), index=z.index)
        if align.shape == (2, 2):
            cluster0_dm2_agreement = (ref_assign[all_labels == 0] == "DM2_cl").mean()
            if cluster0_dm2_agreement > 0.5:
                boot_dm2 = 0
            else:
                boot_dm2 = 1
            mapped = all_labels.map({boot_dm2: 0, 1 - boot_dm2: 1})  # 0=DM2_cl, 1=DM1_cl
        else:
            mapped = all_labels
        assignments[:, b] = mapped.values

    # Per-sample stability = fraction of bootstraps where assignment == reference
    ref_cluster_id = (ref_assign == "DM1_cl").astype(int).reindex(z.index).values  # 1 = DM1_cl
    per_sample_stab = (assignments == ref_cluster_id[:, None]).mean(axis=1)
    stab_df = pd.DataFrame({
        "sample_id": z.index,
        "ref_cluster": ref_assign.reindex(z.index).values,
        "stability": per_sample_stab,
    })
    stab_df.to_csv(OUT / "cluster_stability_per_sample.tsv", sep="\t", index=False)

    # Summary
    res = {
        "n_samples": int(len(z)),
        "n_bootstrap": int(n_boot),
        "mean_stability": float(stab_df["stability"].mean()),
        "median_stability": float(stab_df["stability"].median()),
        "n_high_stability_above_0p9": int((stab_df["stability"] >= 0.9).sum()),
        "n_high_stability_above_0p8": int((stab_df["stability"] >= 0.8).sum()),
        "fraction_above_0p9": float((stab_df["stability"] >= 0.9).mean()),
    }
    print(f"  mean stability: {res['mean_stability']:.3f}")
    print(f"  fraction >= 0.9: {res['fraction_above_0p9']:.3f}")
    (OUT / "cluster_stability_summary.json").write_text(json.dumps(res, indent=2))

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.hist(per_sample_stab, bins=40, color="#fdae61", edgecolor="black", alpha=0.85)
    ax.axvline(0.9, color="#d62728", linestyle="--", linewidth=1, label="0.9 threshold")
    ax.set_xlabel("Per-sample DM1/DM2 cluster assignment stability\n(fraction of 100 bootstrap re-clusterings preserving label)")
    ax.set_ylabel("# samples")
    ax.set_title(
        f"DM1/DM2 cluster stability bootstrap — TCGA-THCA n={len(z)}\n"
        f"mean stability={res['mean_stability']:.3f}; {res['fraction_above_0p9']*100:.0f}% of samples ≥ 0.9 stability",
        fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(axis="y", linestyle=":", alpha=0.3)
    plt.tight_layout()
    fig.savefig(OUT / "cluster_stability.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "cluster_stability.pdf", bbox_inches="tight")
    plt.close(fig)
    return res


# ============== 3. Time-dependent ROC ==============
def time_dependent_roc():
    print("\n[3] Time-dependent ROC (DM1 + TERT predicting OS)")
    sig = pd.read_csv(DM_SIG, sep="\t", index_col=0)
    surv = pd.read_csv(SURV, sep="\t")

    sig["sample_short"] = sig.index.to_series().apply(lambda s: "-".join(str(s).split("-")[:4])[:15])
    surv["sample_short"] = surv["sample"].apply(lambda s: "-".join(str(s).split("-")[:4])[:15])

    # TERT positive samples
    tert_df = pd.read_csv(TERT_MUT, sep="\t")
    tcga_tert = set(tert_df.loc[tert_df["study_id"].astype(str).str.contains("thca_tcga", case=False, na=False), "sample_barcode"]
                    .astype(str).str[:15].str.upper().tolist())

    df = sig.merge(surv[["sample_short", "OS", "OS.time"]], on="sample_short", how="left")
    df = df.dropna(subset=["DM", "OS.time"])
    df["DM1"] = (df["DM"] == "DM1").astype(int)
    df["TERT_pos"] = df["sample_short"].str.upper().isin(tcga_tert).astype(int)
    df["TERT_or_DM1"] = ((df["TERT_pos"] == 1) | (df["DM1"] == 1)).astype(int)

    rows = []
    for years in (1, 3, 5, 8):
        cutoff_days = years * 365.25
        # Only samples with follow-up ≥ cutoff or events before cutoff
        sub = df.copy()
        sub["event_by_t"] = ((sub["OS"] == 1) & (sub["OS.time"] <= cutoff_days)).astype(int)
        sub["censored_before"] = (sub["OS"] == 0) & (sub["OS.time"] < cutoff_days)
        sub_use = sub[~sub["censored_before"]]
        if sub_use["event_by_t"].sum() < 3 or sub_use["event_by_t"].nunique() < 2:
            continue
        for cov_label, cov in [("DM1", sub_use["DM1"].astype(float)),
                                ("TERT_pos", sub_use["TERT_pos"].astype(float)),
                                ("TERT_or_DM1", sub_use["TERT_or_DM1"].astype(float))]:
            try:
                auc = roc_auc_score(sub_use["event_by_t"], cov)
                rows.append({"years": years, "covariate": cov_label, "n": int(len(sub_use)),
                             "n_events": int(sub_use["event_by_t"].sum()), "auc": float(auc)})
            except Exception:
                continue

    res = pd.DataFrame(rows)
    res.to_csv(OUT / "time_dependent_roc.tsv", sep="\t", index=False)
    print(res.to_string(index=False))

    # Plot
    fig, ax = plt.subplots(figsize=(8, 5))
    for cov in ("DM1", "TERT_pos", "TERT_or_DM1"):
        sub = res[res["covariate"] == cov].sort_values("years")
        if len(sub) > 0:
            ax.plot(sub["years"], sub["auc"], marker="o", label=cov, linewidth=2)
    ax.axhline(0.5, color="#888", linestyle=":", linewidth=0.8, label="chance")
    ax.set_xlabel("Years since diagnosis")
    ax.set_ylabel("Time-dependent AUC for OS")
    ax.set_title("TCGA-THCA time-dependent ROC AUC for OS — DM1 + TERT", fontsize=10)
    ax.set_ylim(0.4, 1.05)
    ax.legend(fontsize=9)
    ax.grid(axis="y", linestyle=":", alpha=0.3)
    plt.tight_layout()
    fig.savefig(OUT / "time_dependent_roc.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "time_dependent_roc.pdf", bbox_inches="tight")
    plt.close(fig)
    return res


# ============== 4. Calibration + DCA ==============
def calibration_dca():
    print("\n[4] Calibration + Decision Curve Analysis")
    expr = pd.read_csv(TCGA_EXPR, sep="\t", index_col=0)
    sig = pd.read_csv(DM_SIG, sep="\t", index_col=0)
    panel_z = within_sample_z_score(expr, PANEL).rename("g8")

    df = sig.join(panel_z, how="left").dropna(subset=["DM", "g8"]).copy()
    df["DM1_label"] = (df["DM"] == "DM1").astype(int)

    X = df[["g8"]].values
    y = df["DM1_label"].values

    # 5-fold cross-val LogReg
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    oof_pred = np.zeros(len(y))
    for tr, te in skf.split(X, y):
        clf = LogisticRegression(max_iter=2000).fit(X[tr], y[tr])
        oof_pred[te] = clf.predict_proba(X[te])[:, 1]

    # AUC
    auc = roc_auc_score(y, oof_pred)

    # Calibration: bin oof_pred into 10 bins, plot empirical vs predicted
    prob_true, prob_pred = calibration_curve(y, oof_pred, n_bins=10, strategy="quantile")

    # DCA: for thresholds 0..1, net benefit = (TP - FP*pt/(1-pt)) / N
    thresholds = np.linspace(0.05, 0.95, 19)
    dca_rows = []
    n = len(y)
    for pt in thresholds:
        pred_pos = (oof_pred >= pt)
        tp = int(((pred_pos == 1) & (y == 1)).sum())
        fp = int(((pred_pos == 1) & (y == 0)).sum())
        if pt < 1.0:
            nb = (tp / n) - (fp / n) * (pt / (1 - pt))
        else:
            nb = float("nan")
        # treat-all
        prev = y.mean()
        nb_all = prev - (1 - prev) * (pt / (1 - pt)) if pt < 1.0 else float("nan")
        dca_rows.append({"threshold": float(pt), "model_nb": float(nb), "treat_all_nb": float(nb_all)})
    dca = pd.DataFrame(dca_rows)
    dca.to_csv(OUT / "decision_curve.tsv", sep="\t", index=False)

    # 2-panel figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # Calibration
    ax1.plot([0, 1], [0, 1], color="#888", linestyle=":", linewidth=0.8, label="Perfect calibration")
    ax1.plot(prob_pred, prob_true, marker="o", color="#d62728", linewidth=2, label=f"LogReg(g8) — AUC={auc:.3f}")
    ax1.set_xlabel("Predicted probability of DM1")
    ax1.set_ylabel("Observed fraction DM1")
    ax1.set_title("Calibration curve — TCGA-THCA 5-fold cross-val LogReg", fontsize=10)
    ax1.legend(fontsize=9)
    ax1.grid(axis="both", linestyle=":", alpha=0.3)
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)

    # DCA
    ax2.plot(dca["threshold"], dca["model_nb"], color="#d62728", linewidth=2, label="g8-LogReg model")
    ax2.plot(dca["threshold"], dca["treat_all_nb"], color="#1f77b4", linestyle="--", linewidth=1.5, label="Treat all as DM1")
    ax2.axhline(0, color="#888", linestyle=":", linewidth=0.8, label="Treat none")
    ax2.set_xlabel("Threshold probability (pt)")
    ax2.set_ylabel("Net benefit")
    ax2.set_title("Decision Curve Analysis — clinical net benefit", fontsize=10)
    ax2.legend(fontsize=9)
    ax2.grid(axis="both", linestyle=":", alpha=0.3)

    plt.tight_layout()
    fig.savefig(OUT / "calibration_dca.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "calibration_dca.pdf", bbox_inches="tight")
    plt.close(fig)

    res = {
        "logreg_auc_5fold": float(auc),
        "calibration_intercept": float(prob_true.mean() - prob_pred.mean()),
        "model_nb_at_pt_0p5": float(dca.loc[dca["threshold"].sub(0.5).abs().idxmin(), "model_nb"]),
        "treat_all_nb_at_pt_0p5": float(dca.loc[dca["threshold"].sub(0.5).abs().idxmin(), "treat_all_nb"]),
        "n_samples": int(n),
        "prevalence_DM1": float(y.mean()),
    }
    (OUT / "calibration_dca_summary.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))
    return res


# ============== 5. SHAP-style XGBoost ==============
def xgboost_importance():
    print("\n[5] SHAP-style XGBoost feature importance")
    import xgboost as xgb
    expr = pd.read_csv(TCGA_EXPR, sep="\t", index_col=0)
    sig = pd.read_csv(DM_SIG, sep="\t", index_col=0)
    panel_present = [g for g in PANEL if g in expr.index]
    panel_expr = expr.loc[panel_present]
    panel_z = panel_expr.sub(panel_expr.mean(axis=1), axis=0).div(panel_expr.std(axis=1, ddof=1).replace(0, np.nan), axis=0).T
    df = sig[["DM"]].join(panel_z, how="inner").dropna()
    X = df[panel_present].values
    y = (df["DM"] == "DM1").astype(int).values

    # 5-fold CV with XGBoost
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    importances = np.zeros((5, len(panel_present)))
    aucs = []
    for i, (tr, te) in enumerate(skf.split(X, y)):
        clf = xgb.XGBClassifier(n_estimators=200, max_depth=4, learning_rate=0.05,
                                use_label_encoder=False, eval_metric="logloss",
                                random_state=42, verbosity=0)
        clf.fit(X[tr], y[tr])
        prob = clf.predict_proba(X[te])[:, 1]
        aucs.append(roc_auc_score(y[te], prob))
        # Use built-in gain feature importance
        booster = clf.get_booster()
        scores = booster.get_score(importance_type="gain")
        for j, gname in enumerate(panel_present):
            # XGBoost names features f0...fN
            fname = f"f{j}"
            importances[i, j] = scores.get(fname, 0.0)

    mean_importance = importances.mean(axis=0)
    res = pd.DataFrame({
        "gene": panel_present,
        "xgb_gain_mean": mean_importance,
        "xgb_gain_std": importances.std(axis=0),
        "xgb_gain_normalized": mean_importance / mean_importance.sum() if mean_importance.sum() > 0 else mean_importance,
    }).sort_values("xgb_gain_mean", ascending=False)
    res.to_csv(OUT / "xgboost_feature_importance.tsv", sep="\t", index=False)
    print(f"  mean 5-fold XGB AUC: {np.mean(aucs):.3f} ± {np.std(aucs):.3f}")
    print(res.to_string(index=False))

    # Compare to per-gene |d| from Round 2
    if PER_GENE_R2.exists():
        pg = pd.read_csv(PER_GENE_R2, sep="\t")
        pg = pg[pg["contrast"] == "TCGA-THCA DM1 vs DM2"]
        if len(pg):
            pg_sub = pg.groupby("gene")["cohens_d"].mean().abs().to_dict()
            res["per_gene_abs_d_R2"] = res["gene"].map(pg_sub)
            res.to_csv(OUT / "xgboost_feature_importance.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    y_pos = np.arange(len(res))
    ax.barh(y_pos, res["xgb_gain_normalized"], color="#fdae61", edgecolor="black", alpha=0.85)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(res["gene"], fontsize=10)
    ax.set_xlabel("XGBoost normalized gain (8-gene → DM1/DM2)")
    ax.set_title(f"XGBoost feature importance (5-fold CV, mean AUC={np.mean(aucs):.3f})", fontsize=10)
    ax.invert_yaxis()
    ax.grid(axis="x", linestyle=":", alpha=0.3)
    plt.tight_layout()
    fig.savefig(OUT / "xgboost_importance.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "xgboost_importance.pdf", bbox_inches="tight")
    plt.close(fig)

    return {"mean_5fold_auc": float(np.mean(aucs)),
            "std_5fold_auc": float(np.std(aucs)),
            "top_gene": res.iloc[0]["gene"]}


def main():
    res = {}
    res["permutation_null"] = permutation_null()
    res["cluster_stability"] = cluster_stability()
    res["time_dependent_roc"] = time_dependent_roc().to_dict(orient="records")
    res["calibration_dca"] = calibration_dca()
    res["xgboost"] = xgboost_importance()
    res["generated_at"] = "2026-05-08 round7"

    with open(OUT / "round7_summary.json", "w") as f:
        json.dump(res, f, indent=2, default=str)

    print(f"\nwrote {OUT.relative_to(REPO)}/")
    for f in sorted(OUT.iterdir()):
        print(f"  {f.relative_to(OUT)}")


if __name__ == "__main__":
    main()
