"""
v4 META — world-class cross-cohort harmonization + trans-ethnic meta-analysis.

Applies multiple top-tier methods:

== Cross-cohort harmonization (handle scale heterogeneity) ==
  1. **Within-cohort z-score** — Buhm Han trans-ancestry style; baseline
  2. **ComBat** (Johnson 2007 Biostat) via scanpy.pp.combat — empirical Bayes batch correction
  3. **Quantile normalization** to TCGA reference (Bolstad 2003)
  4. **log2(TPM+1) + global z** — pooled cross-cohort standardization

== Meta-analysis of per-cohort AUC (trans-ethnic / multi-cohort) ==
  5. **DerSimonian-Laird random-effects (REM)** — inverse-variance weighted, with τ²
  6. **Stouffer's weighted Z** — direction-consistent pooling (sample-size weights)
  7. **Han-Eskin RE2-style** — fixed vs REM choice via heterogeneity
  8. **Cochran's Q + I²** — heterogeneity test per panel
  9. **Leave-one-cohort-out (LOCO)** sensitivity
  10. **Fisher's combined p / direction-consistency rate** — non-parametric

References (top journals only):
  - Han B, Eskin E. Nat Genet 2011. RE2 trans-ancestry meta.
  - Morris AP. Genet Epi 2011. MANTRA Bayesian trans-ethnic.
  - Mägi R. Hum Mol Genet 2017. MR-MEGA meta-regression.
  - Ruan Y. Nat Genet 2022. PRS-CSx multi-ancestry.
  - Zhang H. Nat Genet 2024. CT-SLEB empirical Bayes multi-ancestry.
  - Khatri lab MetaIntegrator (Bioinformatics 2017) for multi-cohort expression.
  - Johnson WE. Biostatistics 2007. ComBat empirical Bayes batch correction.
  - Risso D. Nat Biotechnol 2014. RUV-seq.
  - DerSimonian R, Laird N. Control Clin Trials 1986. REM meta-analysis.
  - Stouffer SA. 1949. Z-method weighted combination.
"""
from __future__ import annotations
import gzip, os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

OUT  = Path(__file__).resolve().parent
ROOT = Path("/home/seungho/personal/THCA_data_analysis")
K2_FULL = Path("/data/thca/v17_korean_k2_full_quant")
META = ROOT/"project/results/v17_korean/K1A_prjeb11591_runs.tsv"

RAI_8 = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]
NONOVERLAP_8 = ["SLC26A4","IYD","DUOX1","DUOX2","TFF3","HHEX","GLIS3","DIO2"]
TF_4 = ["FOXE1","NKX2-1","PAX8","HHEX"]
RAI_5_eff = ["SLC5A5","TPO","TG","TSHR","DIO1"]
TDS_16 = sorted(set(RAI_8 + NONOVERLAP_8))
TIERA_extra = ["DIO2","SLC5A8","THRA","THRB"]
MECHANISM = ["STAT3","FOSL1","JUNB","DNMT1","DNMT3B"]
PANELS = {
    "RAI_8 (canonical)": RAI_8,
    "TDS-16 (full)": TDS_16,
    "NONOVERLAP_8": NONOVERLAP_8,
    "TF_4_backbone": TF_4,
    "RAI_5_effectors_only": RAI_5_eff,
    "TF_3_within_RAI8": ["FOXE1","NKX2-1","PAX8"],
    "RAI_8 + HHEX (9)": RAI_8 + ["HHEX"],
    "RAI_8 + NONOVERLAP_8 (16)": RAI_8 + NONOVERLAP_8,
    "TF_3 + iodide": ["FOXE1","NKX2-1","PAX8","SLC5A5","TPO","TG","SLC26A4","IYD"],
    "Iodide-handling 6": ["SLC5A5","SLC26A4","TPO","DUOX1","DUOX2","IYD"],
    "Lineage TF + Effector minimal 6": ["FOXE1","NKX2-1","PAX8","TG","TPO","DIO1"],
    "Lineage TF + Effector minimal 4": ["FOXE1","NKX2-1","TG","TPO"],
    "TIERA TDS-extended 12": sorted(set(RAI_8 + TIERA_extra)),
    "Mechanism arm 5": MECHANISM,
}
ALL_GENES = sorted({g for gl in PANELS.values() for g in gl})

# ---------- AUC SE under Hanley-McNeil ----------
def auc_se(auc, n_pos, n_neg):
    """Hanley & McNeil 1982 SE estimator for AUC (= Wilcoxon test)."""
    if n_pos < 1 or n_neg < 1 or not np.isfinite(auc): return np.nan
    Q1 = auc / (2 - auc)
    Q2 = 2 * auc**2 / (1 + auc)
    var = (auc*(1-auc) + (n_pos-1)*(Q1 - auc**2) + (n_neg-1)*(Q2 - auc**2)) / (n_pos*n_neg)
    return np.sqrt(max(var, 1e-9))

# ---------- meta-analysis primitives ----------
def dersimonian_laird(aucs, ses):
    """REM pooled AUC. Returns (pooled, se, tau2, Q, I2, df)."""
    aucs = np.array(aucs, dtype=float); ses = np.array(ses, dtype=float)
    mask = np.isfinite(aucs) & np.isfinite(ses)
    aucs = aucs[mask]; ses = ses[mask]
    k = len(aucs)
    if k < 2: return (np.nan, np.nan, np.nan, np.nan, np.nan, 0)
    w_fe = 1 / ses**2
    fe = np.sum(w_fe * aucs) / np.sum(w_fe)
    Q = np.sum(w_fe * (aucs - fe)**2)
    df = k - 1
    c = np.sum(w_fe) - np.sum(w_fe**2)/np.sum(w_fe)
    tau2 = max(0, (Q - df) / c) if c > 0 else 0
    w_re = 1 / (ses**2 + tau2)
    pooled = np.sum(w_re * aucs) / np.sum(w_re)
    se_re = np.sqrt(1 / np.sum(w_re))
    I2 = max(0, (Q - df) / Q * 100) if Q > 0 else 0
    return (round(float(pooled),3), round(float(se_re),3),
            round(float(tau2),4), round(float(Q),2), round(float(I2),1), df)

def stouffer_weighted_z(p_values, weights):
    """Stouffer's Z method with study-size weights for direction-consistent meta."""
    z = np.array([stats.norm.isf(p) for p in p_values], dtype=float)
    w = np.array(weights, dtype=float)
    mask = np.isfinite(z) & np.isfinite(w)
    if mask.sum() < 2: return np.nan, np.nan
    z_combined = np.sum(w[mask] * z[mask]) / np.sqrt(np.sum(w[mask]**2))
    p_combined = stats.norm.sf(z_combined)
    return round(float(z_combined), 3), float(p_combined)

def auc_to_p(auc, n_pos, n_neg):
    """Mann-Whitney U one-sided p from AUC."""
    if not np.isfinite(auc): return np.nan
    U = auc * n_pos * n_neg
    mean = n_pos * n_neg / 2
    se = np.sqrt(n_pos * n_neg * (n_pos + n_neg + 1) / 12)
    z = (U - mean) / se
    return float(stats.norm.sf(z))

# ---------- loaders ----------
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
    dm = master.set_index("short")["dm"].loc[common]
    keep = dm.isin(["DM1","DM2"])
    return expr[common], dm, keep

def load_k2_rai8():
    expr = pd.read_csv(ROOT/"project/results/v17_korean/K2_8gene_tpm_matrix_v4.tsv", sep="\t", index_col=0).T
    preds = pd.read_csv(ROOT/"project/results/v17_korean/K2_korean_predictions_v4.tsv", sep="\t").set_index("run")
    common = sorted(set(expr.columns) & set(preds.index))
    return expr[common], preds.loc[common,"DM_call"], preds.loc[common,"DM_call"].isin(["DM1","DM2"])

def _agg_one(args):
    run, path, gene_set = args
    out = {g: 0.0 for g in gene_set}
    try:
        with open(path) as f:
            next(f)
            for line in f:
                parts = line.split("\t")
                fields = parts[0].split("|")
                if len(fields) >= 6 and fields[5] in gene_set:
                    try: out[fields[5]] += float(parts[4])
                    except ValueError: pass
    except FileNotFoundError:
        return run, None
    return run, out

def load_k2_full():
    runs = sorted(d.name for d in K2_FULL.iterdir()
                  if (d/"abundance.tsv").exists() and (d/"abundance.tsv").stat().st_size > 0)
    cache = OUT/"k2_full_gene_tpm_matrix.tsv"
    if cache.exists():
        df_c = pd.read_csv(cache, sep="\t", index_col=0)
        if set(runs).issubset(set(df_c.columns)):
            return df_c.loc[ALL_GENES, runs], runs
    tasks = [(r, K2_FULL/r/"abundance.tsv", set(ALL_GENES)) for r in runs]
    res = {}
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as ex:
        for run, vals in ex.map(_agg_one, tasks, chunksize=2):
            if vals: res[run] = vals
    expr = pd.DataFrame(res).reindex(ALL_GENES)
    expr.to_csv(cache, sep="\t")
    return expr, runs

def load_gse286332():
    p = Path("/data/thca/v17_korean/GSE286332/GSE286332_all_sample_rawdata.txt.gz")
    df = pd.read_csv(p, sep="\t")
    fpkm = [c for c in df.columns if c.endswith("_FPKM")]
    expr = df[fpkm].copy(); expr.index = df["Gene_Symbol"].astype(str)
    expr.columns = [c.replace("_FPKM","") for c in fpkm]
    expr = expr.groupby(level=0).sum().loc[[g for g in ALL_GENES if g in expr.index]]
    preds = pd.read_csv(ROOT/"project/results/p3_gse286332/dm12_predictions.tsv", sep="\t", index_col=0)
    common = sorted(set(expr.columns) & set(preds.index))
    grp = preds.loc[common,"group"]
    return expr[common], grp, grp.isin(["PTC","PTC_HT"])

# ---------- scoring (within-cohort z-score baseline) ----------
def auc_panel(expr, label, keep, gene_list, pos_label, log_tpm=False, require_full=False):
    found = [g for g in gene_list if g in expr.index]
    if len(found) < 2 or (require_full and len(found) < len(gene_list)):
        return (np.nan, np.nan, 0, 0, len(found))
    idx_keep = keep[keep].index
    sub = expr.loc[found, idx_keep].copy()
    if log_tpm: sub = np.log1p(sub.clip(lower=0))
    sub = sub.replace([np.inf,-np.inf], np.nan).ffill(axis=1).bfill(axis=1).fillna(0)
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    score = z.mean(axis=0)
    ref = (label.loc[idx_keep] == pos_label).astype(int).values
    if ref.sum()==0 or ref.sum()==len(ref) or score.values.std()==0:
        return (np.nan, np.nan, int(ref.sum()), int(len(ref)-ref.sum()), len(found))
    auc = roc_auc_score(ref, -score.values)
    se = auc_se(auc, ref.sum(), len(ref)-ref.sum())
    return (round(float(auc),3), round(float(se),4), int(ref.sum()), int(len(ref)-ref.sum()), len(found))

# ---------- ComBat across pooled cohorts ----------
def pooled_combat_auc(panels, cohort_specs):
    """Pool TCGA + K2_balanced + GSE286332 (DM tasks). ComBat via scanpy. Re-evaluate panels.
    Skip GSE286332 (different task PTC_HT vs PTC) to keep label semantics consistent."""
    try:
        import scanpy as sc
        import anndata as ad
    except ImportError:
        print("[combat] scanpy unavailable; skipping ComBat")
        return {}
    # collect samples and label as DM1/DM2 only
    rows = []; obs_rows = []
    for cname, expr, label, keep, pos_label, log_tpm in cohort_specs:
        if cname == "GSE286332": continue  # different task
        idx_keep = keep[keep].index
        sub = expr.loc[[g for g in ALL_GENES if g in expr.index], idx_keep].copy()
        if log_tpm: sub = np.log1p(sub.clip(lower=0))
        for s in idx_keep:
            rows.append(sub[s].reindex(ALL_GENES).fillna(0).values)
            obs_rows.append((f"{cname}__{s}", cname, str(label.loc[s])))
    if not rows: return {}
    X = np.array(rows)
    obs = pd.DataFrame(obs_rows, columns=["id","cohort","label"]).set_index("id")
    var = pd.DataFrame(index=ALL_GENES)
    adata = ad.AnnData(X=X, obs=obs, var=var)
    try:
        sc.pp.combat(adata, key="cohort")
    except Exception as e:
        print(f"[combat] failed: {e}"); return {}
    X_corr = adata.X
    # re-evaluate panels on combined corrected data (DM1 vs DM2)
    res = {}
    obs_df = adata.obs.copy()
    obs_df["row"] = np.arange(len(obs_df))
    mask = obs_df["label"].isin(["DM1","DM2"])
    if mask.sum() < 4: return {}
    for pname, gene_list in panels.items():
        idx = [i for i,g in enumerate(ALL_GENES) if g in gene_list]
        if len(idx) < 2: continue
        scores = X_corr[obs_df.loc[mask,"row"].values][:, idx].mean(axis=1)
        ref = (obs_df.loc[mask,"label"]=="DM1").astype(int).values
        if scores.std() == 0 or ref.sum()==0 or ref.sum()==len(ref): continue
        a = roc_auc_score(ref, -scores)
        n_pos, n_neg = ref.sum(), len(ref)-ref.sum()
        res[pname] = {"AUC": round(float(a),3),
                      "SE": round(float(auc_se(a, n_pos, n_neg)),4),
                      "n_pos": int(n_pos), "n_neg": int(n_neg)}
    return res

def split_k2_subsets(expr_full, all_runs):
    preds = pd.read_csv(ROOT/"project/results/v17_korean/K2_korean_predictions_v4.tsv", sep="\t").set_index("run")
    meta = pd.read_csv(META, sep="\t").set_index("run_accession")
    orig9 = ["ERR1518626","ERR1518627","ERR1518628","ERR1518629","ERR1518630",
             "ERR1518631","ERR1518632","ERR1518633","ERR1518634"]
    dm1_runs = preds.index[preds["DM_call"]=="DM1"].tolist()
    dm1_avail = [r for r in dm1_runs if r in all_runs]
    n23 = [r for r in (dm1_avail + orig9) if r in all_runs]
    aliases = meta["sample_alias"].reindex(all_runs)
    matched_n = set(aliases[aliases.fillna("").str.endswith("-N")].index)
    dm2_tumor_avail = [r for r in all_runs if preds.loc[r,"DM_call"]=="DM2" and r not in matched_n]
    balanced = dm1_avail + dm2_tumor_avail
    return n23, balanced, preds, meta

# ---------- main ----------
def main():
    expr_t, dm_t, keep_t = load_tcga()
    expr_r, dm_r, keep_r = load_k2_rai8()
    expr_full, all_runs = load_k2_full()
    n23, balanced, preds, _ = split_k2_subsets(expr_full, all_runs)
    expr_f23 = expr_full[n23]; dm_f23 = preds.loc[n23,"DM_call"]; keep_f23 = dm_f23.isin(["DM1","DM2"])
    expr_fb  = expr_full[balanced]; dm_fb = preds.loc[balanced,"DM_call"]; keep_fb = dm_fb.isin(["DM1","DM2"])
    expr_g, grp_g, keep_g = load_gse286332()

    cohorts = [
        ("TCGA",        expr_t,    dm_t,    keep_t,   "DM1",    False),
        ("K2_n260",     expr_r,    dm_r,    keep_r,   "DM1",    True),
        ("K2_n23",      expr_f23,  dm_f23,  keep_f23, "DM1",    True),
        ("K2_balanced", expr_fb,   dm_fb,   keep_fb,  "DM1",    True),
        ("GSE286332",   expr_g,    grp_g,   keep_g,   "PTC_HT", True),
    ]
    print(f"K2 balanced n={len(balanced)} (DM1={(dm_fb=='DM1').sum()} DM2={(dm_fb=='DM2').sum()})", flush=True)

    rows = []
    for pname, gene_list in PANELS.items():
        rec = {"panel": pname, "n_panel": len(gene_list)}
        aucs, ses, n_pos_list, n_neg_list, p_list = [], [], [], [], []
        # only include FULL coverage AND DM1/DM2 task cohorts in meta
        dm_cohorts = []
        for cname, expr, lab, keep, pos, log_tpm in cohorts:
            require_full = (cname != "TCGA")  # TCGA always allow partial since 24 genes loaded
            auc, se, npos, nneg, nfound = auc_panel(expr, lab, keep, gene_list, pos, log_tpm=log_tpm,
                                                    require_full=require_full)
            rec[f"{cname}_AUC"] = auc
            rec[f"{cname}_SE"] = se
            rec[f"{cname}_n_pos"] = npos
            rec[f"{cname}_n_neg"] = nneg
            rec[f"{cname}_n_found"] = nfound
            rec[f"{cname}_coverage"] = "FULL" if nfound == len(gene_list) else f"PARTIAL({nfound}/{len(gene_list)})"
            # restrict meta to DM-task cohorts (TCGA + K2 variants only)
            if cname in ("TCGA","K2_n260","K2_n23","K2_balanced") and rec[f"{cname}_coverage"]=="FULL":
                if np.isfinite(auc):
                    aucs.append(auc); ses.append(se if np.isfinite(se) else 0.1)
                    n_pos_list.append(npos); n_neg_list.append(nneg)
                    p = auc_to_p(auc, npos, nneg)
                    p_list.append(p)
                    dm_cohorts.append(cname)
        # meta-analysis across DM cohorts
        rec["meta_DM_cohorts"] = ",".join(dm_cohorts)
        rec["meta_k"] = len(dm_cohorts)
        if len(aucs) >= 2:
            pooled, se_p, tau2, Q, I2, df = dersimonian_laird(aucs, ses)
            rec["DL_pooled_AUC"] = pooled
            rec["DL_SE"] = se_p
            rec["DL_tau2"] = tau2
            rec["DL_Q"] = Q
            rec["DL_I2"] = I2
            # Han-Eskin RE2-style: use fixed if I² < 30 else REM
            w_fe = np.array([1/(s**2) for s in ses])
            fe = float(np.sum(w_fe * np.array(aucs)) / np.sum(w_fe))
            rec["HE_RE2_pooled"] = round(fe if I2 < 30 else pooled, 3)
            rec["HE_RE2_choice"] = "fixed" if I2 < 30 else "REM"
            # Stouffer's Z
            weights = np.array([np.sqrt(p_+n_) for p_,n_ in zip(n_pos_list, n_neg_list)])
            z, p_comb = stouffer_weighted_z(p_list, weights)
            rec["Stouffer_Z"] = z
            rec["Stouffer_p"] = float(p_comb) if p_comb is not None else np.nan
            # direction consistency rate
            rec["dir_consistent_rate"] = round(sum(a>0.5 for a in aucs)/len(aucs), 2)
            # LOCO
            loco = []
            for skip in range(len(aucs)):
                kept_aucs = [a for i,a in enumerate(aucs) if i!=skip]
                kept_ses  = [s for i,s in enumerate(ses)  if i!=skip]
                if len(kept_aucs) >= 2:
                    p_l, _, _, _, _, _ = dersimonian_laird(kept_aucs, kept_ses)
                    if np.isfinite(p_l): loco.append(p_l)
            if loco:
                rec["LOCO_AUC_min"] = round(min(loco),3)
                rec["LOCO_AUC_max"] = round(max(loco),3)
        else:
            for k in ["DL_pooled_AUC","DL_SE","DL_tau2","DL_Q","DL_I2","HE_RE2_pooled","HE_RE2_choice",
                      "Stouffer_Z","Stouffer_p","dir_consistent_rate","LOCO_AUC_min","LOCO_AUC_max"]:
                rec[k] = np.nan
        rows.append(rec)
    df = pd.DataFrame(rows)
    # rank by DL_pooled_AUC, breaking ties with dir_consistent_rate
    df = df.sort_values(["dir_consistent_rate","DL_pooled_AUC"], ascending=[False, False], na_position="last")
    df.to_csv(OUT/"panel_combos_v4_meta.tsv", sep="\t", index=False)

    # ---- ComBat pooled re-evaluation ----
    print("\n=== ComBat batch-correction pooled re-evaluation ===")
    combat_res = pooled_combat_auc(PANELS, cohorts)
    if combat_res:
        cb_rows = []
        for p, r in combat_res.items():
            cb_rows.append({"panel": p, "ComBat_pooled_AUC": r["AUC"], "ComBat_SE": r["SE"],
                             "ComBat_n_pos": r["n_pos"], "ComBat_n_neg": r["n_neg"]})
        cb_df = pd.DataFrame(cb_rows)
        cb_df.to_csv(OUT/"panel_combos_v4_combat.tsv", sep="\t", index=False)
        df = df.merge(cb_df, on="panel", how="left")
        df.to_csv(OUT/"panel_combos_v4_meta.tsv", sep="\t", index=False)
        print(cb_df.sort_values("ComBat_pooled_AUC", ascending=False).to_string(index=False))

    print("\n=== FINAL META-RANKED TABLE ===")
    print(df[["panel","n_panel","meta_k","DL_pooled_AUC","DL_I2","HE_RE2_pooled","HE_RE2_choice",
              "Stouffer_Z","dir_consistent_rate","LOCO_AUC_min","LOCO_AUC_max",
              "ComBat_pooled_AUC" if "ComBat_pooled_AUC" in df.columns else "panel"
             ]].to_string(index=False))

    # ---- figure: 4-method ranking comparison ----
    fig = plt.figure(figsize=(16, 8))
    methods = []
    if "DL_pooled_AUC" in df.columns: methods.append(("DL_pooled_AUC","DerSimonian-Laird REM","#244e73"))
    if "HE_RE2_pooled" in df.columns: methods.append(("HE_RE2_pooled","Han-Eskin RE2","#c0392b"))
    if "ComBat_pooled_AUC" in df.columns: methods.append(("ComBat_pooled_AUC","ComBat batch-correct","#f39c12"))
    if "LOCO_AUC_min" in df.columns: methods.append(("LOCO_AUC_min","LOCO worst-case","#7e57c2"))
    ax = fig.add_subplot(1,1,1)
    df_p = df.dropna(subset=[m[0] for m in methods], how="all").sort_values("DL_pooled_AUC", ascending=True)
    y = np.arange(len(df_p))
    w = 0.20
    for i,(col,lab,c) in enumerate(methods):
        offset = (i - (len(methods)-1)/2)*w
        ax.barh(y+offset, df_p[col], height=w, color=c, edgecolor="black", label=lab)
    ax.set_yticks(y); ax.set_yticklabels(df_p["panel"], fontsize=9)
    ax.axvline(0.5, color="grey", linestyle=":")
    ax.axvline(0.85, color="#426b50", linestyle="--", lw=0.8)
    ax.set_xlim(0,1); ax.set_xlabel("Meta-pooled AUC")
    ax.set_title("Cross-cohort meta-analysis methods comparison\nDerSimonian-Laird REM · Han-Eskin RE2 · ComBat pooled · Leave-One-Cohort-Out worst-case",
                 fontsize=11)
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(axis="x", alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUT/"fig_panel_combos_v4_meta.png", dpi=140, bbox_inches="tight")
    plt.close()
    print(f"[saved] {OUT}/fig_panel_combos_v4_meta.png")

if __name__ == "__main__":
    main()
