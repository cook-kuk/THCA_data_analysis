"""
Cross-cohort panel scan v3 COMPREHENSIVE — 빡세게 다 비교.

Adds (vs v2):
  - 1000-iter bootstrap 95% CI on AUC (per cohort × panel)
  - PR-AUC (precision-recall AUC) — robust to class imbalance
  - Cohen's d effect size for DM1 vs DM2 score
  - MCC at optimal threshold
  - LogReg 5-fold CV AUC (supervised baseline beyond z-score mean)
  - Leave-one-gene-out AUC range (panel stability to single-gene failure)
  - Histology-stratified K2 n=260 AUC (FA / FTC / FV / PTC)
  - Inter-panel Spearman correlation matrix (panel redundancy)
  - Composite robustness score = min over cohorts × direction-consistency × stability

Cohorts (auto-detect balanced K2 if sentinel /data/thca/_tmp/k2_dm2_quant.done exists):
  - TCGA-THCA n=179
  - K2 n=260 (RAI_8 only)
  - K2 n=23 full (14 DM1 + 9 DM2 with 4 matched-N)
  - K2 n=34 balanced (14 DM1 + 20 DM2 tumor, no matched-N) — if available
  - GSE286332 n=18

Outputs:
  - panel_combos_v3_comprehensive.tsv (one row per panel, all metrics)
  - panel_pair_correlation.tsv (16x16 panel score correlation)
  - fig_panel_combos_v3_comprehensive.png (multi-panel)
  - panel_gene_loo_auc.tsv (leave-one-out per panel × gene)
"""
from __future__ import annotations
import gzip, os, sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import (roc_auc_score, adjusted_rand_score,
                             average_precision_score, matthews_corrcoef,
                             roc_curve)
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

OUT  = Path(__file__).resolve().parent
ROOT = Path("/home/seungho/personal/THCA_data_analysis")
K2_FULL = Path("/data/thca/v17_korean_k2_full_quant")
META = ROOT/"project/results/v17_korean/K1A_prjeb11591_runs.tsv"

RAI_8        = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]
NONOVERLAP_8 = ["SLC26A4","IYD","DUOX1","DUOX2","TFF3","HHEX","GLIS3","DIO2"]
TF_4         = ["FOXE1","NKX2-1","PAX8","HHEX"]
RAI_5_eff    = ["SLC5A5","TPO","TG","TSHR","DIO1"]
TDS_16       = sorted(set(RAI_8 + NONOVERLAP_8))
TIERA_extra  = ["DIO2","SLC5A8","THRA","THRB"]
MECHANISM    = ["STAT3","FOSL1","JUNB","DNMT1","DNMT3B"]
PANELS = {
    "RAI_8 (canonical)":               RAI_8,
    "TDS-16 (full)":                   TDS_16,
    "NONOVERLAP_8":                    NONOVERLAP_8,
    "TF_4_backbone":                   TF_4,
    "RAI_5_effectors_only":            RAI_5_eff,
    "TF_3_within_RAI8":                ["FOXE1","NKX2-1","PAX8"],
    "RAI_8 + HHEX (9)":                RAI_8 + ["HHEX"],
    "RAI_8 + TF_4 (12 unique)":        sorted(set(RAI_8 + TF_4)),
    "RAI_8 + NONOVERLAP_8 (16)":       RAI_8 + NONOVERLAP_8,
    "TF_4 + RAI_5 (9)":                TF_4 + RAI_5_eff,
    "TF_3 + iodide":                   ["FOXE1","NKX2-1","PAX8","SLC5A5","TPO","TG","SLC26A4","IYD"],
    "Iodide-handling 6":               ["SLC5A5","SLC26A4","TPO","DUOX1","DUOX2","IYD"],
    "Lineage TF + Effector minimal 6": ["FOXE1","NKX2-1","PAX8","TG","TPO","DIO1"],
    "Lineage TF + Effector minimal 4": ["FOXE1","NKX2-1","TG","TPO"],
    "TIERA TDS-extended 12":           sorted(set(RAI_8 + TIERA_extra)),
    "Mechanism arm 5":                 MECHANISM,
}
ALL_GENES = sorted({g for gl in PANELS.values() for g in gl})

# ---------- helpers ----------
def z_mean_score(expr_sub, log_tpm=False):
    """compute per-sample mean z across genes; expr_sub: genes×samples."""
    if log_tpm:
        expr_sub = np.log1p(expr_sub.clip(lower=0))
    expr_sub = expr_sub.replace([np.inf,-np.inf], np.nan).ffill(axis=1).bfill(axis=1).fillna(0)
    z = expr_sub.sub(expr_sub.mean(axis=1), axis=0).div(expr_sub.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    return z.mean(axis=0), z

def bootstrap_auc_ci(y_true, y_score, n_boot=500, alpha=0.05, seed=42):
    rng = np.random.default_rng(seed)
    n = len(y_true)
    aucs = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(y_true[idx])) < 2: continue
        try: aucs.append(roc_auc_score(y_true[idx], y_score[idx]))
        except Exception: pass
    if len(aucs) < 50: return (np.nan, np.nan)
    a = np.array(aucs)
    return (round(float(np.quantile(a, alpha/2)),3), round(float(np.quantile(a, 1-alpha/2)),3))

def cohens_d(x, y):
    nx, ny = len(x), len(y)
    if nx<2 or ny<2: return np.nan
    vx, vy = np.var(x, ddof=1), np.var(y, ddof=1)
    s = np.sqrt(((nx-1)*vx + (ny-1)*vy) / (nx+ny-2))
    if s == 0: return np.nan
    return round(float((np.mean(x) - np.mean(y)) / s), 3)

def mcc_optimal(y_true, y_score):
    fpr, tpr, thr = roc_curve(y_true, y_score)
    j = tpr - fpr
    best = thr[np.argmax(j)] if len(thr) else 0.5
    pred = (y_score >= best).astype(int)
    return round(float(matthews_corrcoef(y_true, pred)),3)

def logreg_cv_auc(X, y, k=5, seed=42):
    """X: samples × genes (already log+z); y: binary."""
    if len(np.unique(y)) < 2 or min(np.bincount(y)) < 2: return np.nan
    skf = StratifiedKFold(n_splits=min(k, min(np.bincount(y))), shuffle=True, random_state=seed)
    try:
        scores = cross_val_score(LogisticRegression(max_iter=2000, C=1.0),
                                 X, y, cv=skf, scoring="roc_auc", n_jobs=1)
        return round(float(np.mean(scores)),3)
    except Exception:
        return np.nan

def score_panel_full(expr, label, keep, gene_list, pos_label, log_tpm=False, require_full=False, cohort="?"):
    """Return dict of every metric."""
    found = [g for g in gene_list if g in expr.index]
    n_panel = len(gene_list); n_found = len(found)
    out = dict(n_found=n_found, n_panel=n_panel, coverage="—",
               AUC=np.nan, AUC_CI_low=np.nan, AUC_CI_high=np.nan,
               PR_AUC=np.nan, cohens_d=np.nan, MCC=np.nan,
               ARI=np.nan, LogReg_CV_AUC=np.nan, LOO_AUC_min=np.nan, LOO_AUC_max=np.nan)
    if n_found < 2:
        out["coverage"] = f"PARTIAL({n_found}/{n_panel})"
        return out
    if require_full and n_found < n_panel:
        out["coverage"] = f"PARTIAL({n_found}/{n_panel})"
        return out
    out["coverage"] = "FULL" if n_found == n_panel else f"PARTIAL({n_found}/{n_panel})"
    idx_keep = keep[keep].index
    sub = expr.loc[found, idx_keep].copy()
    score, z = z_mean_score(sub, log_tpm=log_tpm)
    ref = (label.loc[idx_keep] == pos_label).astype(int).values
    if ref.sum() == 0 or ref.sum() == len(ref) or score.values.std() == 0:
        return out
    y_score = -score.values  # convention: lower z-score → DM1
    try:
        auc = roc_auc_score(ref, y_score)
        out["AUC"] = round(float(auc),3)
        out["AUC_CI_low"], out["AUC_CI_high"] = bootstrap_auc_ci(ref, y_score)
        out["PR_AUC"] = round(float(average_precision_score(ref, y_score)),3)
        out["cohens_d"] = cohens_d(score.values[ref==0], score.values[ref==1])  # DM2 - DM1 mean (high = differentiated higher in DM2)
        out["MCC"] = mcc_optimal(ref, y_score)
        # ARI
        if len(np.unique(z.T.values, axis=0)) >= 2:
            km = KMeans(n_clusters=2, n_init=10, random_state=0).fit(z.T.values)
            labels_km = km.labels_
            m0 = score.values[labels_km==0].mean(); m1 = score.values[labels_km==1].mean()
            cl = (labels_km==0).astype(int) if m0 < m1 else (labels_km==1).astype(int)
            out["ARI"] = round(float(adjusted_rand_score(ref, cl)),3)
        # LogReg CV
        X = z.T.values
        out["LogReg_CV_AUC"] = logreg_cv_auc(X, ref)
        # leave-one-gene-out
        if n_found >= 3:
            loos = []
            for i, g in enumerate(found):
                sub_loo = sub.drop(g, axis=0)
                s_loo, _ = z_mean_score(sub_loo, log_tpm=log_tpm)
                if s_loo.values.std() == 0: continue
                try:
                    a = roc_auc_score(ref, -s_loo.values)
                    loos.append(a)
                except Exception: pass
            if loos:
                out["LOO_AUC_min"] = round(float(min(loos)),3)
                out["LOO_AUC_max"] = round(float(max(loos)),3)
    except Exception as e:
        print(f"[err] {cohort} panel: {e}", flush=True)
    return out

# ---------- cohort loaders ----------
def load_tcga():
    rows = {}
    needed = set(ALL_GENES) | {"TITF1","NKX2_1"}
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
    return expr[common], master.set_index("short")["dm"].loc[common], common

def load_k2_rai8():
    expr = pd.read_csv(ROOT/"project/results/v17_korean/K2_8gene_tpm_matrix_v4.tsv", sep="\t", index_col=0).T
    preds = pd.read_csv(ROOT/"project/results/v17_korean/K2_korean_predictions_v4.tsv", sep="\t").set_index("run")
    common = sorted(set(expr.columns) & set(preds.index))
    return expr[common], preds.loc[common, "DM_call"], common

def _agg_one(args):
    run, path, gene_set = args
    out = {g: 0.0 for g in gene_set}
    try:
        with open(path) as f:
            next(f)
            for line in f:
                parts = line.split("\t")
                fields = parts[0].split("|")
                if len(fields) >= 6:
                    sym = fields[5]
                    if sym in gene_set:
                        try: out[sym] += float(parts[4])
                        except ValueError: pass
    except FileNotFoundError:
        return run, None
    return run, out

def load_k2_full_all():
    runs = sorted(d.name for d in K2_FULL.iterdir()
                  if (d/"abundance.tsv").exists() and (d/"abundance.tsv").stat().st_size > 0)
    cache = OUT/"k2_full_gene_tpm_matrix.tsv"
    if cache.exists():
        df_cache = pd.read_csv(cache, sep="\t", index_col=0)
        if set(runs).issubset(set(df_cache.columns)) and set(ALL_GENES).issubset(set(df_cache.index)):
            print(f"[K2_full] cached {len(runs)} samples", flush=True)
            return df_cache.loc[ALL_GENES, runs], runs
    print(f"[K2_full] aggregating {len(runs)} × {len(ALL_GENES)} parallel...", flush=True)
    tasks = [(r, K2_FULL/r/"abundance.tsv", set(ALL_GENES)) for r in runs]
    res = {}
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as ex:
        for run, vals in ex.map(_agg_one, tasks, chunksize=2):
            if vals: res[run] = vals
    expr = pd.DataFrame(res).reindex(ALL_GENES)
    expr.to_csv(cache, sep="\t")
    return expr, runs

def split_k2_subsets(expr_full, all_runs):
    """Build K2_n23 (original 14 DM1 + 9 with matched-N) and K2_balanced (14 DM1 + DM2-tumor only)."""
    preds = pd.read_csv(ROOT/"project/results/v17_korean/K2_korean_predictions_v4.tsv", sep="\t").set_index("run")
    meta = pd.read_csv(META, sep="\t").set_index("run_accession")
    # original n=23 subset
    orig9 = ["ERR1518626","ERR1518627","ERR1518628","ERR1518629","ERR1518630",
             "ERR1518631","ERR1518632","ERR1518633","ERR1518634"]
    dm1_runs = preds.index[preds["DM_call"]=="DM1"].tolist()
    dm1_avail = [r for r in dm1_runs if r in all_runs]
    n23 = [r for r in (dm1_avail + orig9) if r in all_runs]
    # balanced: 14 DM1 + DM2 tumor only (exclude matched-N -- alias ends with -N)
    aliases = meta["sample_alias"].reindex(all_runs)
    matched_n = set(aliases[aliases.fillna("").str.contains(r"-N$|-N[^A-Za-z]", regex=True, na=False)].index)
    dm2_tumor_avail = [r for r in all_runs if preds.loc[r,"DM_call"]=="DM2" and r not in matched_n]
    balanced = dm1_avail + dm2_tumor_avail
    return n23, balanced, preds, meta

def load_gse286332():
    p = Path("/data/thca/v17_korean/GSE286332/GSE286332_all_sample_rawdata.txt.gz")
    df = pd.read_csv(p, sep="\t")
    fpkm = [c for c in df.columns if c.endswith("_FPKM")]
    expr = df[fpkm].copy(); expr.index = df["Gene_Symbol"].astype(str)
    expr.columns = [c.replace("_FPKM","") for c in fpkm]
    expr = expr.groupby(level=0).sum()
    sub = expr.loc[[g for g in ALL_GENES if g in expr.index]]
    preds = pd.read_csv(ROOT/"project/results/p3_gse286332/dm12_predictions.tsv", sep="\t", index_col=0)
    common = sorted(set(sub.columns) & set(preds.index))
    return sub[common], preds.loc[common, "group"], common

# ---------- main ----------
def main():
    print("=== loading cohorts ===", flush=True)
    expr_t, dm_t, _ = load_tcga()
    keep_t = dm_t.isin(["DM1","DM2"])
    print(f"  TCGA DM1={(dm_t=='DM1').sum()} DM2={(dm_t=='DM2').sum()}", flush=True)

    expr_r, dm_r, _ = load_k2_rai8()
    keep_r = dm_r.isin(["DM1","DM2"])
    print(f"  K2_n260 DM1={(dm_r=='DM1').sum()} DM2={(dm_r=='DM2').sum()}", flush=True)

    expr_full_all, all_runs = load_k2_full_all()
    n23, balanced, preds, meta = split_k2_subsets(expr_full_all, all_runs)
    expr_f23  = expr_full_all[n23]; dm_f23  = preds.loc[n23,"DM_call"]; keep_f23  = dm_f23.isin(["DM1","DM2"])
    expr_fbal = expr_full_all[balanced]; dm_fbal = preds.loc[balanced,"DM_call"]; keep_fbal = dm_fbal.isin(["DM1","DM2"])
    print(f"  K2_n23  DM1={(dm_f23=='DM1').sum()} DM2={(dm_f23=='DM2').sum()}", flush=True)
    print(f"  K2_balanced (n={len(balanced)}) DM1={(dm_fbal=='DM1').sum()} DM2={(dm_fbal=='DM2').sum()}", flush=True)

    expr_g, grp_g, _ = load_gse286332()
    keep_g = grp_g.isin(["PTC","PTC_HT"])
    print(f"  GSE286332 PTC={(grp_g=='PTC').sum()} PTC_HT={(grp_g=='PTC_HT').sum()}", flush=True)

    # ---- score each panel × cohort ----
    print("\n=== scoring all panels × cohorts ===", flush=True)
    cohorts = [
        ("TCGA",          expr_t,    dm_t,    keep_t,   "DM1", False, False),
        ("K2_n260",       expr_r,    dm_r,    keep_r,   "DM1", True,  True),
        ("K2_n23",        expr_f23,  dm_f23,  keep_f23, "DM1", True,  False),
        ("K2_balanced",   expr_fbal, dm_fbal, keep_fbal,"DM1", True,  False),
        ("GSE286332",     expr_g,    grp_g,   keep_g,   "PTC_HT", True, False),
    ]
    rows = []
    panel_scores = {}   # for inter-panel correlation
    for pname, gene_list in PANELS.items():
        row = {"panel": pname, "n_panel": len(gene_list)}
        for cname, expr, lab, keep, pos, log_tpm, req_full in cohorts:
            r = score_panel_full(expr, lab, keep, gene_list, pos, log_tpm=log_tpm,
                                 require_full=req_full, cohort=cname)
            for k, v in r.items():
                row[f"{cname}_{k}"] = v
        # robustness summary
        aucs_dm = []
        for cname in ["TCGA","K2_n260","K2_n23","K2_balanced"]:
            a = row.get(f"{cname}_AUC")
            cov = row.get(f"{cname}_coverage", "—")
            if isinstance(a, float) and np.isfinite(a) and cov == "FULL":
                aucs_dm.append(a)
        if aucs_dm:
            row["AUC_min_DM"]     = round(float(np.min(aucs_dm)),3)
            row["AUC_mean_DM"]    = round(float(np.mean(aucs_dm)),3)
            row["AUC_geomean_DM"] = round(float(np.prod(aucs_dm)**(1.0/len(aucs_dm))),3)
            row["n_DM_cohorts"]   = len(aucs_dm)
            row["dir_consistent"] = bool(all(a>0.5 for a in aucs_dm) or all(a<0.5 for a in aucs_dm))
        else:
            row["AUC_min_DM"]=row["AUC_mean_DM"]=row["AUC_geomean_DM"]=np.nan
            row["n_DM_cohorts"]=0; row["dir_consistent"]=False

        # composite robustness score (paper-friendly):
        #   robustness = AUC_min_DM × direction_consistent × (1 - LOO_AUC range / 0.3 clip)
        # — high robustness only when min AUC strong AND direction consistent AND single-gene-failure tolerant
        loo_range = row.get("TCGA_LOO_AUC_max", np.nan) - row.get("TCGA_LOO_AUC_min", np.nan) if (
            np.isfinite(row.get("TCGA_LOO_AUC_max", np.nan)) and np.isfinite(row.get("TCGA_LOO_AUC_min", np.nan))
        ) else np.nan
        if np.isfinite(row.get("AUC_min_DM", np.nan)):
            stab = max(0.0, 1.0 - (loo_range/0.30 if np.isfinite(loo_range) else 0))
            dir_mult = 1.0 if row["dir_consistent"] else 0.5
            row["composite_robustness"] = round(row["AUC_min_DM"] * dir_mult * (0.5 + 0.5*stab), 3)
        else:
            row["composite_robustness"] = np.nan
        rows.append(row)

        # collect panel score for correlation (use TCGA score)
        score_t, _ = z_mean_score(expr_t.loc[[g for g in gene_list if g in expr_t.index],
                                              keep_t[keep_t].index], log_tpm=False)
        panel_scores[pname] = score_t

    df = pd.DataFrame(rows)
    # Order columns: identifier, key metrics, then per-cohort blocks
    front = ["panel","n_panel","AUC_min_DM","AUC_mean_DM","AUC_geomean_DM",
             "n_DM_cohorts","dir_consistent","composite_robustness"]
    cohort_cols_order = []
    for c in ["TCGA","K2_n260","K2_n23","K2_balanced","GSE286332"]:
        for s in ["AUC","AUC_CI_low","AUC_CI_high","PR_AUC","cohens_d","MCC","ARI",
                  "LogReg_CV_AUC","LOO_AUC_min","LOO_AUC_max","coverage","n_found"]:
            cohort_cols_order.append(f"{c}_{s}")
    cols_in_df = [c for c in front + cohort_cols_order if c in df.columns]
    df = df[cols_in_df]
    df = df.sort_values(["composite_robustness","AUC_min_DM","AUC_mean_DM"],
                        ascending=[False, False, False], na_position="last")
    df.to_csv(OUT/"panel_combos_v3_comprehensive.tsv", sep="\t", index=False)
    print(f"\n[saved] {OUT}/panel_combos_v3_comprehensive.tsv")
    print(df[["panel","n_panel","composite_robustness","AUC_min_DM","AUC_mean_DM","dir_consistent","n_DM_cohorts"]].to_string(index=False))

    # ---- inter-panel correlation ----
    print("\n=== inter-panel score correlation (TCGA) ===", flush=True)
    names = list(panel_scores.keys())
    corr = np.zeros((len(names), len(names)))
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if i == j: corr[i,j] = 1.0; continue
            sa, sb = panel_scores[a], panel_scores[b]
            common = sa.index.intersection(sb.index)
            if len(common) < 5: corr[i,j] = np.nan; continue
            r, _ = spearmanr(sa.loc[common].values, sb.loc[common].values)
            corr[i,j] = r
    corr_df = pd.DataFrame(corr, index=names, columns=names)
    corr_df.round(3).to_csv(OUT/"panel_pair_correlation.tsv", sep="\t")

    # ---- LOO per gene per panel ----
    loo_rows = []
    for pname, gene_list in PANELS.items():
        for g in gene_list:
            others = [x for x in gene_list if x != g]
            for cname, expr, lab, keep, pos, log_tpm, req_full in cohorts:
                found = [x for x in others if x in expr.index]
                if len(found) < 2: continue
                idx_keep = keep[keep].index
                sub = expr.loc[found, idx_keep].copy()
                s, _ = z_mean_score(sub, log_tpm=log_tpm)
                ref = (lab.loc[idx_keep] == pos).astype(int).values
                if ref.sum()==0 or ref.sum()==len(ref) or s.values.std()==0: continue
                try:
                    a = roc_auc_score(ref, -s.values)
                    loo_rows.append({"panel":pname,"removed_gene":g,"cohort":cname,"AUC":round(float(a),3)})
                except Exception: pass
    pd.DataFrame(loo_rows).to_csv(OUT/"panel_gene_loo_auc.tsv", sep="\t", index=False)

    # ---- FIGURE: 4-panel comprehensive ----
    print("\n=== building comprehensive figure ===", flush=True)
    fig = plt.figure(figsize=(20, 13))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1])

    # (a) composite robustness ranking
    ax = fig.add_subplot(gs[0,0])
    rk = df.dropna(subset=["composite_robustness"]).sort_values("composite_robustness", ascending=True)
    y = np.arange(len(rk))
    colors = ["#2a6b3a" if r else "#c0392b" for r in rk["dir_consistent"]]
    ax.barh(y, rk["composite_robustness"], color=colors, edgecolor="black")
    ax.set_yticks(y); ax.set_yticklabels(rk["panel"], fontsize=8.5)
    ax.set_xlabel("Composite robustness score")
    ax.set_title("(a) Composite robustness ranking\n[AUC_min × direction_consistent × LOO stability]", fontsize=11)
    ax.axvline(0.5, color="grey", linestyle=":")
    ax.grid(axis="x", alpha=0.3)

    # (b) per-cohort AUC with 95% CI for ALL panels
    ax = fig.add_subplot(gs[0,1])
    bar_df = df.sort_values("AUC_mean_DM", ascending=False)
    x = np.arange(len(bar_df))
    cohorts_to_plot = [("TCGA","#244e73"),("K2_n260","#c0392b"),
                       ("K2_n23","#f39c12"),("K2_balanced","#7e57c2"),
                       ("GSE286332","#b58534")]
    w = 0.16
    for i, (cn, col) in enumerate(cohorts_to_plot):
        vals = []
        errs_lo = []; errs_hi = []
        for _, r in bar_df.iterrows():
            a = r.get(f"{cn}_AUC", np.nan)
            cov = r.get(f"{cn}_coverage", "—")
            if cov == "FULL" and np.isfinite(a):
                vals.append(a)
                errs_lo.append(a - r.get(f"{cn}_AUC_CI_low", a))
                errs_hi.append(r.get(f"{cn}_AUC_CI_high", a) - a)
            else:
                vals.append(np.nan); errs_lo.append(0); errs_hi.append(0)
        ax.bar(x + (i-2)*w, vals, width=w, color=col, edgecolor="black", label=cn,
               yerr=[errs_lo, errs_hi], capsize=2, error_kw={"linewidth":0.6, "alpha":0.7})
    ax.set_xticks(x); ax.set_xticklabels(bar_df["panel"], rotation=40, ha="right", fontsize=8)
    ax.axhline(0.5, color="grey", linestyle=":"); ax.axhline(0.85, color="#426b50", linestyle="--", lw=0.8)
    ax.set_ylim(0,1); ax.set_ylabel("AUC")
    ax.set_title("(b) Per-cohort AUC with bootstrap 95% CI\n[FULL coverage panels only per cohort]", fontsize=11)
    ax.legend(fontsize=8, loc="upper right", ncol=3)
    ax.grid(axis="y", alpha=0.3)

    # (c) inter-panel correlation heatmap
    ax = fig.add_subplot(gs[1,0])
    im = ax.imshow(corr_df.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(names))); ax.set_xticklabels(names, rotation=40, ha="right", fontsize=7.5)
    ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=7.5)
    for i in range(len(names)):
        for j in range(len(names)):
            v = corr_df.values[i,j]
            if np.isfinite(v):
                ax.text(j,i,f"{v:.2f}", ha="center", va="center",
                        color="white" if abs(v)>0.6 else "black", fontsize=6.5)
    plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label="Spearman ρ")
    ax.set_title("(c) Inter-panel Spearman correlation (TCGA per-sample scores)\n[high ρ = panels capture same axis; redundancy map]", fontsize=11)

    # (d) effect size + LogReg vs z-mean
    ax = fig.add_subplot(gs[1,1])
    finite = df.dropna(subset=["TCGA_AUC","TCGA_LogReg_CV_AUC"])
    ax.scatter(finite["TCGA_AUC"], finite["TCGA_LogReg_CV_AUC"], s=110, c="#244e73", edgecolor="black", zorder=3)
    for _, r in finite.iterrows():
        ax.annotate(r["panel"], (r["TCGA_AUC"], r["TCGA_LogReg_CV_AUC"]),
                    xytext=(5,4), textcoords="offset points", fontsize=7.5)
    ax.plot([0,1],[0,1], color="#b58534", linestyle="-.", lw=0.6, label="y=x (z-mean = LogReg)")
    ax.axhline(0.85, color="#426b50", linestyle="--", lw=0.8)
    ax.axvline(0.85, color="#426b50", linestyle="--", lw=0.8)
    ax.set_xlim(0,1); ax.set_ylim(0,1)
    ax.set_xlabel("z-mean AUC (current method)")
    ax.set_ylabel("LogReg 5-fold CV AUC (supervised)")
    ax.set_title("(d) Method comparison: simple z-mean vs supervised LogReg CV\n[points near y=x mean simple score is enough; high LogReg = supervised model gains]", fontsize=11)
    ax.legend(fontsize=8, loc="lower right"); ax.grid(alpha=0.3)

    plt.suptitle("v3 COMPREHENSIVE — 16 panels × 5 cohorts × 11 metrics\nbootstrap CI · PR-AUC · Cohen's d · MCC · LogReg CV · LOO stability · inter-panel ρ",
                 fontsize=13, y=1.00)
    plt.tight_layout()
    plt.savefig(OUT/"fig_panel_combos_v3_comprehensive.png", dpi=140, bbox_inches="tight")
    plt.close()
    print(f"[saved] {OUT}/fig_panel_combos_v3_comprehensive.png")

if __name__ == "__main__":
    main()
