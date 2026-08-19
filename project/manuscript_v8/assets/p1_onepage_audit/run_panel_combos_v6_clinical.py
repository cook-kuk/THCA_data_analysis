"""
v6 CLINICAL — paper-grade additional analyses.

Adds:
  - **Permutation null distribution** (5000 shuffles) → empirical p-value vs random
  - **Cox PH survival analysis** (TCGA OS) — does each panel score predict mortality?
  - **BRAF / RAS / unknown driver-stratified AUC** — biology subset
  - **NRI / IDI vs RAI_8 reference** — Net Reclassification Improvement
  - **ROC curves overlay** — all panels in one figure, TCGA
  - **Calibration curves** (reliability diagrams)
  - **Decision curve analysis (Vickers 2006)** — clinical net benefit

Outputs:
  - panel_combos_v6_clinical.tsv
  - fig_roc_overlay.png
  - fig_survival_km.png (Kaplan-Meier by panel-score median split)
  - fig_calibration.png
  - fig_decision_curve.png
"""
from __future__ import annotations
import gzip
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve, brier_score_loss
from sklearn.calibration import calibration_curve
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

OUT  = Path(__file__).resolve().parent
ROOT = Path("/home/seungho/personal/THCA_data_analysis")

PANELS = {
    "RAI_8 (canonical)":  ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"],
    "TF_3_within_RAI8":   ["FOXE1","NKX2-1","PAX8"],
    "NONOVERLAP_8":       ["SLC26A4","IYD","DUOX1","DUOX2","TFF3","HHEX","GLIS3","DIO2"],
    "TDS-16 (full)":      sorted(set(["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]
                                      + ["SLC26A4","IYD","DUOX1","DUOX2","TFF3","HHEX","GLIS3","DIO2"])),
    "Lineage TF + Effector minimal 6": ["FOXE1","NKX2-1","PAX8","TG","TPO","DIO1"],
    "Mechanism arm 5":    ["STAT3","FOSL1","JUNB","DNMT1","DNMT3B"],
}

def load_tcga():
    rows = {}
    needed = set([g for gl in PANELS.values() for g in gl]) | {"TITF1","NKX2_1"}
    with gzip.open(ROOT/"project/data/raw/TCGA_pancan/pancan_geneExp.gz","rt") as f:
        header = next(f).strip().split("\t"); samples = header[1:]
        for line in f:
            parts = line.rstrip("\n").split("\t"); sym = parts[0].strip()
            if sym in needed:
                rows[sym] = [float(x) if x not in ("","NA","NaN") else np.nan for x in parts[1:]]
    expr = pd.DataFrame(rows, index=samples).T
    if "NKX2-1" not in expr.index:
        if "NKX2_1" in expr.index: expr = expr.rename(index={"NKX2_1":"NKX2-1"})
        elif "TITF1" in expr.index: expr = expr.rename(index={"TITF1":"NKX2-1"})
    expr = expr[~expr.index.duplicated(keep="first")]
    master = pd.read_csv(ROOT/"project/results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv", sep="\t")
    master["short"] = master["sample_id"].str[:15]
    expr.columns = [c[:15] for c in expr.columns]
    common = sorted(set(expr.columns) & set(master["short"]))
    expr_t = expr[common]
    m = master.set_index("short").loc[common]
    return expr_t, m

def z_score(expr_sub):
    sub = expr_sub.replace([np.inf,-np.inf], np.nan).ffill(axis=1).bfill(axis=1).fillna(0)
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    return z.mean(axis=0)   # higher = differentiated

def main():
    expr, m = load_tcga()
    print(f"TCGA: {expr.shape}, DM1={int((m['dm']=='DM1').sum())} DM2={int((m['dm']=='DM2').sum())}")
    keep = m['dm'].isin(["DM1","DM2"])
    idx_keep = keep[keep].index
    ref = (m.loc[idx_keep, 'dm'] == 'DM1').astype(int).values  # 1=DM1

    # ---- compute per-panel score (higher = differentiated; we use -score for DM1) ----
    panel_scores = {}
    for pname, gene_list in PANELS.items():
        found = [g for g in gene_list if g in expr.index]
        if len(found) < 2: continue
        sub = expr.loc[found, idx_keep]
        s = z_score(sub)
        panel_scores[pname] = -s.values   # higher score → more likely DM1
    # for survival/Cox, also compute on FULL TCGA cohort (not just DM-labeled)
    panel_scores_all = {}
    for pname, gene_list in PANELS.items():
        found = [g for g in gene_list if g in expr.index]
        if len(found) < 2: continue
        sub = expr.loc[found, :]
        s = z_score(sub)
        panel_scores_all[pname] = -s

    # ============================================================
    # 1) Permutation null + empirical p
    # ============================================================
    print("\n=== permutation null (5000 shuffles) ===")
    n_perm = 5000
    rng = np.random.default_rng(42)
    perm_results = {}
    for pname, sc in panel_scores.items():
        true_auc = roc_auc_score(ref, sc)
        null_aucs = np.empty(n_perm)
        for i in range(n_perm):
            null_aucs[i] = roc_auc_score(rng.permutation(ref), sc)
        # two-sided empirical p
        p = (np.sum(np.abs(null_aucs - 0.5) >= abs(true_auc - 0.5)) + 1) / (n_perm + 1)
        perm_results[pname] = dict(
            true_auc=round(float(true_auc),3),
            null_mean=round(float(null_aucs.mean()),3),
            null_p05=round(float(np.quantile(null_aucs, 0.025)),3),
            null_p95=round(float(np.quantile(null_aucs, 0.975)),3),
            empirical_p=round(float(p),5)
        )
        print(f"  {pname}: AUC={true_auc:.3f}, null_mean={null_aucs.mean():.3f}, p={p:.4f}")

    # ============================================================
    # 2) Cox PH OS survival (TCGA, all samples with os_days+os_event)
    # ============================================================
    print("\n=== Cox PH OS survival ===")
    surv_df = m[["os_days","os_event","age","sex","driver_anchor"]].copy()
    surv_df["age_num"] = pd.to_numeric(surv_df["age"], errors="coerce")
    surv_df["sex_M"]   = (surv_df["sex"].astype(str).str.lower()=="male").astype(int)
    surv_df["BRAF"]    = (surv_df["driver_anchor"]=="BRAF").astype(int)
    surv_df["RAS"]     = (surv_df["driver_anchor"]=="RAS").astype(int)
    surv_df = surv_df.dropna(subset=["os_days","os_event"])
    surv_df["os_event"] = pd.to_numeric(surv_df["os_event"], errors="coerce")
    surv_df["os_days"]  = pd.to_numeric(surv_df["os_days"],  errors="coerce")
    surv_df = surv_df.dropna(subset=["os_days","os_event"])
    surv_df = surv_df[surv_df["os_days"]>0]
    cox_rows = []
    surv_df_idx = surv_df.index
    for pname, s in panel_scores_all.items():
        common = list(set(s.index) & set(surv_df_idx))
        if len(common) < 30: continue
        df_c = surv_df.loc[common, ["os_days","os_event","age_num","sex_M","BRAF","RAS"]].copy()
        df_c["score"] = pd.to_numeric(s.loc[common], errors="coerce")
        df_c = df_c.dropna()
        if df_c["os_event"].sum() < 5: continue
        # univariate Cox
        try:
            cph_u = CoxPHFitter().fit(df_c[["score","os_days","os_event"]], duration_col="os_days", event_col="os_event")
            hr_u = float(np.exp(cph_u.params_["score"]))
            p_u = float(cph_u.summary.loc["score","p"])
        except Exception:
            hr_u = np.nan; p_u = np.nan
        # multivariable (adjust for age + sex + BRAF + RAS)
        try:
            df_mv = df_c[["score","age_num","sex_M","BRAF","RAS","os_days","os_event"]].copy()
            df_mv = df_mv.dropna()
            cph_mv = CoxPHFitter().fit(df_mv, duration_col="os_days", event_col="os_event")
            hr_mv = float(np.exp(cph_mv.params_["score"]))
            p_mv = float(cph_mv.summary.loc["score","p"])
            c_idx = float(cph_mv.concordance_index_)
        except Exception:
            hr_mv = np.nan; p_mv = np.nan; c_idx = np.nan
        cox_rows.append({"panel": pname, "n_surv": len(df_c), "n_events": int(df_c["os_event"].sum()),
                         "HR_per_z_unadj": round(hr_u,3), "p_unadj": round(p_u,4),
                         "HR_per_z_adj_age_sex_BRAF_RAS": round(hr_mv,3),
                         "p_adj": round(p_mv,4), "C_index_mv": round(c_idx,3)})
        print(f"  {pname}: HR_unadj={hr_u:.3f} p={p_u:.4f} | HR_adj={hr_mv:.3f} p={p_mv:.4f} | C={c_idx:.3f}")

    # ============================================================
    # 3) BRAF / RAS / unknown driver-stratified AUC
    # ============================================================
    print("\n=== driver-stratified AUC (DM1 vs DM2 within driver subset) ===")
    strat_rows = []
    for driver in ["BRAF","RAS","unknown"]:
        sub_idx = m.index[(m["driver_anchor"]==driver) & m["dm"].isin(["DM1","DM2"])].intersection(idx_keep)
        sub_ref = (m.loc[sub_idx,"dm"]=="DM1").astype(int).values
        if len(sub_idx) < 10 or sub_ref.sum() < 3 or sub_ref.sum() == len(sub_ref): continue
        for pname, sc_all in panel_scores_all.items():
            sub_score = sc_all.loc[sub_idx].values
            try: auc = roc_auc_score(sub_ref, sub_score)
            except: continue
            strat_rows.append({"driver": driver, "n_subset": len(sub_idx),
                                "n_DM1": int(sub_ref.sum()), "n_DM2": int(len(sub_ref)-sub_ref.sum()),
                                "panel": pname, "AUC": round(float(auc),3)})
    print(f"  {len(strat_rows)} driver-stratified panel evals")

    # ============================================================
    # 4) NRI / IDI vs RAI_8 reference
    # ============================================================
    print("\n=== NRI / IDI vs RAI_8 reference (cross-validated probabilities) ===")
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    nri_idi_rows = []
    def cv_probs(X, y):
        out = np.zeros(len(y))
        for tr, te in skf.split(X, y):
            mdl = LogisticRegression(max_iter=2000, C=1.0).fit(X[tr], y[tr])
            out[te] = mdl.predict_proba(X[te])[:, 1]
        return out
    # RAI_8 probs as reference
    rai_genes = [g for g in PANELS["RAI_8 (canonical)"] if g in expr.index]
    sub_ref_panel = expr.loc[rai_genes, idx_keep]
    X_ref = z_score_per_gene_matrix(sub_ref_panel)
    p_ref = cv_probs(X_ref, ref)
    for pname, gene_list in PANELS.items():
        if pname == "RAI_8 (canonical)": continue
        found = [g for g in gene_list if g in expr.index]
        if len(found) < 2: continue
        X = z_score_per_gene_matrix(expr.loc[found, idx_keep])
        p_new = cv_probs(X, ref)
        # IDI = mean (p_new - p_ref | y=1) - mean (p_new - p_ref | y=0)
        idi = float(np.mean((p_new - p_ref)[ref==1]) - np.mean((p_new - p_ref)[ref==0]))
        # NRI (continuous; Pencina 2008)
        events_up = np.sum((p_new > p_ref)[ref==1]); events_down = np.sum((p_new < p_ref)[ref==1])
        nonev_down = np.sum((p_new < p_ref)[ref==0]); nonev_up = np.sum((p_new > p_ref)[ref==0])
        n_ev = max(1, ref.sum()); n_nev = max(1, len(ref)-ref.sum())
        nri = (events_up - events_down)/n_ev - (nonev_up - nonev_down)/n_nev
        # bootstrap CI on NRI
        rng2 = np.random.default_rng(7)
        nris = []
        for _ in range(500):
            idx = rng2.integers(0, len(ref), len(ref))
            r2 = ref[idx]; pn = p_new[idx]; pr = p_ref[idx]
            e_up = np.sum((pn>pr)[r2==1]); e_dn = np.sum((pn<pr)[r2==1])
            n_up = np.sum((pn>pr)[r2==0]); n_dn = np.sum((pn<pr)[r2==0])
            n_e = max(1, r2.sum()); n_ne = max(1, len(r2)-r2.sum())
            nris.append((e_up-e_dn)/n_e - (n_up-n_dn)/n_ne)
        nri_lo, nri_hi = np.quantile(nris, [0.025, 0.975])
        nri_idi_rows.append({"panel": pname, "NRI_vs_RAI8": round(nri,3),
                             "NRI_CI_low": round(nri_lo,3), "NRI_CI_high": round(nri_hi,3),
                             "IDI_vs_RAI8": round(idi,4)})
        print(f"  {pname}: NRI={nri:+.3f}[{nri_lo:.3f},{nri_hi:.3f}] IDI={idi:+.4f}")

    # ============================================================
    # 5) ROC curves overlay figure
    # ============================================================
    fig, ax = plt.subplots(figsize=(8, 7))
    colors = plt.cm.tab10(np.linspace(0,1,len(panel_scores)))
    for (pname, sc), c in zip(panel_scores.items(), colors):
        fpr, tpr, _ = roc_curve(ref, sc)
        a = roc_auc_score(ref, sc)
        ax.plot(fpr, tpr, color=c, lw=2, label=f"{pname} (AUC={a:.3f})")
    ax.plot([0,1],[0,1], color="grey", linestyle="--")
    ax.set_xlabel("False positive rate"); ax.set_ylabel("True positive rate")
    ax.set_title("TCGA-THCA ROC overlay (z-mean score, DM1 vs DM2, n=179)", fontsize=11)
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(alpha=0.3)
    plt.tight_layout(); plt.savefig(OUT/"fig_roc_overlay.png", dpi=140, bbox_inches="tight"); plt.close()
    print(f"[saved] {OUT}/fig_roc_overlay.png")

    # ============================================================
    # 6) KM survival by panel-score median split (top 3 panels)
    # ============================================================
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    surv_df_keep = surv_df.dropna(subset=["os_days","os_event"])
    for ax, pname in zip(axes, ["RAI_8 (canonical)", "NONOVERLAP_8", "TF_3_within_RAI8"]):
        sc = panel_scores_all[pname]
        common = list(set(sc.index) & set(surv_df_keep.index))
        if len(common) < 20: continue
        d = surv_df_keep.loc[common, ["os_days","os_event"]].copy()
        d["score"] = pd.to_numeric(sc.loc[common], errors="coerce")
        d = d.dropna()
        m_med = d["score"].median()
        d["grp"] = np.where(d["score"]>=m_med, "high", "low")
        # higher score = more differentiated = lower DM1 → expect LONGER survival
        for grp, color in [("low","#c0392b"),("high","#244e73")]:
            sub = d[d["grp"]==grp]
            kmf = KaplanMeierFitter().fit(sub["os_days"], sub["os_event"], label=f"{grp} (n={len(sub)})")
            kmf.plot_survival_function(ax=ax, color=color, ci_show=False)
        # logrank
        lo = d[d["grp"]=="low"]; hi = d[d["grp"]=="high"]
        lr = logrank_test(lo["os_days"], hi["os_days"], lo["os_event"], hi["os_event"])
        ax.set_title(f"{pname}\nlogrank p={lr.p_value:.4f}", fontsize=10)
        ax.set_xlabel("OS days"); ax.set_ylim(0, 1.05); ax.grid(alpha=0.3)
    plt.tight_layout(); plt.savefig(OUT/"fig_survival_km.png", dpi=140, bbox_inches="tight"); plt.close()
    print(f"[saved] {OUT}/fig_survival_km.png")

    # ============================================================
    # 7) Calibration curves (5-fold CV LogReg)
    # ============================================================
    fig, ax = plt.subplots(figsize=(7, 6))
    for (pname, gene_list), c in zip(PANELS.items(), colors):
        found = [g for g in gene_list if g in expr.index]
        if len(found) < 2: continue
        X = z_score_per_gene_matrix(expr.loc[found, idx_keep])
        p = cv_probs(X, ref)
        try:
            frac_pos, mean_pred = calibration_curve(ref, p, n_bins=10, strategy="quantile")
            brier = brier_score_loss(ref, p)
            ax.plot(mean_pred, frac_pos, "o-", color=c, label=f"{pname} (Brier={brier:.3f})", markersize=5)
        except Exception: pass
    ax.plot([0,1],[0,1], color="grey", linestyle="--", label="perfect")
    ax.set_xlabel("Mean predicted probability"); ax.set_ylabel("Fraction of positives (DM1)")
    ax.set_title("Calibration curves (5-fold CV LogReg, TCGA)", fontsize=11)
    ax.legend(fontsize=9, loc="upper left"); ax.grid(alpha=0.3)
    plt.tight_layout(); plt.savefig(OUT/"fig_calibration.png", dpi=140, bbox_inches="tight"); plt.close()
    print(f"[saved] {OUT}/fig_calibration.png")

    # ============================================================
    # 8) Decision curve analysis (Vickers 2006)
    # ============================================================
    fig, ax = plt.subplots(figsize=(8, 6))
    thresholds = np.linspace(0.05, 0.95, 91)
    prev = ref.mean()
    def net_benefit(p, ref, thr):
        pred = (p >= thr).astype(int)
        tp = np.sum((pred==1) & (ref==1)); fp = np.sum((pred==1) & (ref==0))
        n = len(ref)
        return tp/n - fp/n * (thr/(1-thr))
    # treat-all and treat-none
    treat_all = [(prev - (1-prev)*(t/(1-t))) for t in thresholds]
    ax.plot(thresholds, treat_all, color="grey", linestyle=":", label="treat all (DM1)")
    ax.axhline(0, color="black", linestyle="--", lw=0.6, label="treat none")
    for (pname, gene_list), c in zip(PANELS.items(), colors):
        found = [g for g in gene_list if g in expr.index]
        if len(found) < 2: continue
        X = z_score_per_gene_matrix(expr.loc[found, idx_keep])
        p = cv_probs(X, ref)
        nb = [net_benefit(p, ref, t) for t in thresholds]
        ax.plot(thresholds, nb, color=c, label=pname, lw=1.6)
    ax.set_xlabel("Threshold probability"); ax.set_ylabel("Net benefit")
    ax.set_title("Decision curve analysis (Vickers 2006) — TCGA DM1\n[higher = more clinical utility at that threshold]", fontsize=11)
    ax.legend(fontsize=8, loc="upper right"); ax.grid(alpha=0.3); ax.set_ylim(-0.15, 0.55)
    plt.tight_layout(); plt.savefig(OUT/"fig_decision_curve.png", dpi=140, bbox_inches="tight"); plt.close()
    print(f"[saved] {OUT}/fig_decision_curve.png")

    # ---- aggregate output TSV ----
    rows = []
    cox_map = {r["panel"]: r for r in cox_rows}
    nri_map = {r["panel"]: r for r in nri_idi_rows}
    strat_map = {(r["driver"], r["panel"]): r for r in strat_rows}
    for pname in PANELS.keys():
        p = perm_results.get(pname, {})
        c = cox_map.get(pname, {})
        n = nri_map.get(pname, {})
        rows.append({
            "panel": pname,
            "n_panel": len(PANELS[pname]),
            "TCGA_AUC": p.get("true_auc"),
            "perm_null_p": p.get("empirical_p"),
            "perm_null_mean": p.get("null_mean"),
            "perm_null_95CI_low": p.get("null_p05"),
            "perm_null_95CI_high": p.get("null_p95"),
            "CoxOS_HR_unadj": c.get("HR_per_z_unadj"),
            "CoxOS_p_unadj":  c.get("p_unadj"),
            "CoxOS_HR_adj":   c.get("HR_per_z_adj_age_sex_BRAF_RAS"),
            "CoxOS_p_adj":    c.get("p_adj"),
            "CoxOS_Cindex":   c.get("C_index_mv"),
            "AUC_in_BRAF":    strat_map.get(("BRAF",pname),{}).get("AUC"),
            "AUC_in_RAS":     strat_map.get(("RAS",pname),{}).get("AUC"),
            "AUC_in_unknown": strat_map.get(("unknown",pname),{}).get("AUC"),
            "NRI_vs_RAI8":    n.get("NRI_vs_RAI8"),
            "NRI_CI_low":     n.get("NRI_CI_low"),
            "NRI_CI_high":    n.get("NRI_CI_high"),
            "IDI_vs_RAI8":    n.get("IDI_vs_RAI8"),
        })
    df = pd.DataFrame(rows)
    df.to_csv(OUT/"panel_combos_v6_clinical.tsv", sep="\t", index=False)
    print("\n=== FINAL CLINICAL TABLE ===")
    print(df.to_string(index=False))

def z_score_per_gene_matrix(sub):
    """returns samples × genes matrix (np.array) of within-cohort gene z-scores."""
    sub = sub.replace([np.inf,-np.inf], np.nan).ffill(axis=1).bfill(axis=1).fillna(0)
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    return z.T.values

if __name__ == "__main__":
    main()
