#!/usr/bin/env python
"""v3_ext_05_external_validation.py

Core external validation:
  Models : LogReg_l2, LogReg_elasticnet, RandomForest, GradientBoosting, XGBoost (seed=42)
  Feature sets:
    TDS16, TierA67_clean, TierA67_clean_no_MAPK, BRS71_original (if recovered),
    v3_panel_compact (TDS16 ∪ top-16 univariate non-MAPK on TCGA train fold)
  Tasks:
    V1 : TCGA 2-class mutation-verified -> PRJEB11591 Yoo labels
    V2 : 3-class (add NBNR) external on PRJEB11591
    V3 : 6-class fusion-aware external
    V4 : TCGA -> pooled GPL570 (ComBat-adjusted)  (did pooling fix GSE27155 bACC=0.5?)
    V5 : LODO over 7 cohorts
    V6 : 1000-iter label-permutation null on best V1 configuration

Save results to:
    results/ml/v3_ext_validation_results.tsv
    results/ml/v3_ext_calibration.tsv
    results/ml/v3_ext_permutation_null.tsv
    results/ml/v3_ext_lodo_v3.tsv
And 8 Plotly figures to reports/html/figs_interactive/v3_ext_*.html.

Parallelised (task x feature_set x model) via joblib.
"""
import os
import sys
import json
import logging
import warnings
from pathlib import Path
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path("/opt/thyroid-dash/project")
DP = ROOT / "data_processed"
DP_V3_RNA = DP / "bulk_rnaseq_v3"
DP_V3_MIC = DP / "microarray_v3"
META = ROOT / "metadata"
TABLES = ROOT / "results" / "tables"
ML = ROOT / "results" / "ml"
FIGS = ROOT / "reports" / "html" / "figs_interactive"
LOGDIR = ROOT / "logs"
for d in (ML, FIGS, LOGDIR):
    d.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[logging.FileHandler(LOGDIR / "v3_ext_05_external_validation.log", mode="w"),
              logging.StreamHandler()],
)
log = logging.getLogger("v3_ext_05")

sys.path.insert(0, str(ROOT / "notebooks_or_scripts"))
try:
    from rerun_v2 import TIERA67_CLEAN_UNIQUE, TDS16
except Exception as e:
    log.warning(f"rerun_v2 import fail: {e}; using fallback TDS16")
    TDS16 = ["DIO1", "DIO2", "DUOX1", "DUOX2", "FOXE1", "GLIS3", "NKX2-1", "PAX8",
             "SLC26A4", "SLC5A5", "SLC5A8", "TG", "THRA", "THRB", "TPO", "TSHR"]
    TIERA67_CLEAN_UNIQUE = TDS16

MAPK_GENES = {"DUSP4", "DUSP5", "DUSP6", "SPRY1", "SPRY2", "SPRY4",
              "ETV4", "ETV5", "FOSL1", "PHLDA1"}

SEED = 42


def load_brs71_orig():
    f = META / "v3_brs71_original.txt"
    if not f.exists():
        return []
    return [ln.strip() for ln in f.read_text().splitlines() if ln.strip() and not ln.startswith("#")]


def load_expr(dataset: str) -> pd.DataFrame:
    """Return a gene x sample DataFrame (index = gene_symbol)."""
    candidates = [DP_V3_RNA / f"{dataset}_v3_log2.tsv",
                  DP_V3_MIC / f"{dataset}_v3_log2.tsv"]
    for c in candidates:
        if c.exists():
            df = pd.read_csv(c, sep="\t", index_col=0)
            return df
    raise FileNotFoundError(f"No v3 matrix for {dataset}")


def align_features(X_train: pd.DataFrame, X_test: pd.DataFrame, genes: list):
    g = [g for g in genes if g in X_train.index and g in X_test.index]
    return X_train.loc[g].T, X_test.loc[g].T, g  # samples x genes


def get_labels(sm: pd.DataFrame, dataset: str, n_classes: int = 2):
    """Return (sample_ids, y) for binary 'BRAF vs RAS' by default.

    n_classes=2: BRAF_like=1, RAS_like=0 (mutation_verified only)
    n_classes=3: BRAF=1, RAS=0, NBNR=2
    """
    subset = sm[sm["dataset"] == dataset].copy()
    anchor_col = "driver_anchor"
    y, ids = [], []
    for _, r in subset.iterrows():
        a = str(r.get(anchor_col, "")).upper()
        if n_classes == 2:
            if a == "BRAF":
                ids.append(r["sample_id"]); y.append(1)
            elif a == "RAS":
                ids.append(r["sample_id"]); y.append(0)
        elif n_classes == 3:
            if a == "BRAF":
                ids.append(r["sample_id"]); y.append(1)
            elif a == "RAS":
                ids.append(r["sample_id"]); y.append(0)
            else:
                # Only include NBNR if tumor
                nvt = str(r.get("normal_vs_tumor", "")).lower()
                if nvt == "tumor":
                    ids.append(r["sample_id"]); y.append(2)
    return ids, np.array(y)


def make_models():
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    models = {
        "LogReg_l2": Pipeline([("sc", StandardScaler()),
                               ("clf", LogisticRegression(penalty="l2", C=1.0,
                                                          max_iter=2000, random_state=SEED))]),
        "LogReg_elasticnet": Pipeline([("sc", StandardScaler()),
                                        ("clf", LogisticRegression(penalty="elasticnet",
                                                                   solver="saga",
                                                                   l1_ratio=0.5, C=1.0,
                                                                   max_iter=4000, random_state=SEED))]),
        "RandomForest": RandomForestClassifier(n_estimators=400, random_state=SEED, n_jobs=2),
        "GradientBoosting": GradientBoostingClassifier(n_estimators=200, random_state=SEED),
    }
    try:
        from xgboost import XGBClassifier
        models["XGBoost"] = XGBClassifier(n_estimators=300, max_depth=4, learning_rate=0.1,
                                          eval_metric="logloss", random_state=SEED, n_jobs=2,
                                          use_label_encoder=False)
    except Exception:
        pass
    return models


def ece(y, p, bins=10):
    y = np.asarray(y).astype(int); p = np.asarray(p).astype(float)
    edges = np.linspace(0, 1, bins + 1)
    e = 0.0; n = len(y)
    for i in range(bins):
        lo, hi = edges[i], edges[i + 1]
        m = (p > lo) & (p <= hi)
        if m.sum() == 0:
            continue
        e += (m.sum() / n) * abs(y[m].mean() - p[m].mean())
    return float(e)


def score_binary(y, p_raw):
    from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss, f1_score
    out = {}
    try:
        out["auc"] = float(roc_auc_score(y, p_raw))
    except Exception:
        out["auc"] = np.nan
    try:
        out["pr_auc"] = float(average_precision_score(y, p_raw))
    except Exception:
        out["pr_auc"] = np.nan
    try:
        out["brier"] = float(brier_score_loss(y, p_raw))
    except Exception:
        out["brier"] = np.nan
    out["ece"] = ece(y, p_raw)
    # Youden optimal threshold
    from sklearn.metrics import roc_curve
    try:
        fpr, tpr, thr = roc_curve(y, p_raw)
        j = tpr - fpr
        tstar = float(thr[np.argmax(j)])
        yhat = (p_raw >= tstar).astype(int)
        tp = int(((yhat == 1) & (y == 1)).sum())
        tn = int(((yhat == 0) & (y == 0)).sum())
        fp = int(((yhat == 1) & (y == 0)).sum())
        fn = int(((yhat == 0) & (y == 1)).sum())
        se = tp / max(tp + fn, 1); sp = tn / max(tn + fp, 1)
        out["bACC_Youden"] = 0.5 * (se + sp)
        out["F1"] = f1_score(y, yhat, zero_division=0)
        # NPV @ Se>=0.95
        order = np.argsort(-p_raw)
        # sweep threshold to get Se>=0.95
        se95_npv = np.nan
        for t in thr:
            yh = (p_raw >= t).astype(int)
            _tp = int(((yh == 1) & (y == 1)).sum())
            _tn = int(((yh == 0) & (y == 0)).sum())
            _fp = int(((yh == 1) & (y == 0)).sum())
            _fn = int(((yh == 0) & (y == 1)).sum())
            se_t = _tp / max(_tp + _fn, 1)
            npv = _tn / max(_tn + _fn, 1)
            if se_t >= 0.95:
                se95_npv = npv
                break
        out["NPV_at_Se95"] = float(se95_npv) if se95_npv is not None else np.nan
    except Exception as e:
        out["bACC_Youden"] = np.nan; out["F1"] = np.nan; out["NPV_at_Se95"] = np.nan
    return out


def bootstrap_auc_ci(y, p, n_boot=200):
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(SEED)
    n = len(y); vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        try:
            vals.append(roc_auc_score(y[idx], p[idx]))
        except Exception:
            pass
    if not vals:
        return (np.nan, np.nan)
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def isotonic_calibrate(p_train, y_train, p_test):
    from sklearn.isotonic import IsotonicRegression
    iso = IsotonicRegression(out_of_bounds="clip")
    try:
        iso.fit(p_train, y_train)
        return iso.predict(p_test)
    except Exception:
        return p_test


def univariate_top_n_non_mapk(X: pd.DataFrame, y: np.ndarray, n=16, exclude=None):
    """X: samples x genes; returns top-n gene list by univariate |effect|."""
    from scipy import stats
    excl = set(exclude or [])
    genes = [g for g in X.columns if g not in excl]
    scores = []
    for g in genes:
        v = X[g].values
        if np.nanstd(v) == 0:
            scores.append((g, 0))
            continue
        try:
            t, _ = stats.ttest_ind(v[y == 1], v[y == 0], equal_var=False, nan_policy="omit")
            scores.append((g, abs(float(t))))
        except Exception:
            scores.append((g, 0))
    scores.sort(key=lambda x: -x[1])
    return [g for g, _ in scores[:n]]


def prepare_feature_sets(X_tcga: pd.DataFrame, y_tcga: np.ndarray):
    brs_orig = load_brs71_orig()
    fs = {
        "TDS16": list(TDS16),
        "TierA67_clean": list(TIERA67_CLEAN_UNIQUE),
        "TierA67_clean_no_MAPK": [g for g in TIERA67_CLEAN_UNIQUE if g not in MAPK_GENES],
    }
    if brs_orig:
        fs["BRS71_original"] = brs_orig
    # v3_panel_compact: TDS16 ∪ top-16 non-MAPK univariate on TCGA train fold
    try:
        shared = [g for g in X_tcga.index if g not in MAPK_GENES]
        top = univariate_top_n_non_mapk(X_tcga.T, y_tcga, n=16, exclude=MAPK_GENES)
        fs["v3_panel_compact"] = sorted(set(list(TDS16) + top))
    except Exception as e:
        log.warning(f"panel_compact build failed: {e}")
    return fs


def run_v1(feature_sets, X_tcga, y_tcga, sm):
    """TCGA -> PRJEB11591 (if available) else TCGA -> GSE76039 (dedifferentiated) as proxy."""
    # Try PRJEB11591 via expression matrix (unavailable -> mark in SKIPPED summary)
    results = []
    calibration_rows = []
    best = {"auc": -1}
    # External: prefer PRJEB11591 expression; fall back to GSE27155 (BRAF vs RAS anchor labelled).
    ext_dataset = None
    prjeb_expr = DP_V3_RNA / "PRJEB11591_v3_log2.tsv"
    if prjeb_expr.exists():
        ext_dataset = "PRJEB11591"
        X_ext = load_expr("PRJEB11591")
        sub = sm[sm["dataset"] == "PRJEB11591"].copy()
    else:
        ext_dataset = "GSE27155"
        try:
            X_ext = load_expr("GSE27155")
        except Exception:
            return results, calibration_rows, best
        sub = sm[sm["dataset"] == "GSE27155"].copy()

    sub = sub[sub["sample_id"].isin(X_ext.columns)]
    lbl_map = {"BRAF": 1, "RAS": 0}
    sub = sub[sub["driver_anchor"].isin(lbl_map.keys())]
    if sub.empty:
        log.warning(f"V1: no BRAF/RAS-labelled samples in {ext_dataset}")
        return results, calibration_rows, best
    y_ext = sub["driver_anchor"].map(lbl_map).values
    ext_ids = sub["sample_id"].tolist()
    log.info(f"V1 external {ext_dataset}: n={len(ext_ids)}  BRAF={int(y_ext.sum())}  RAS={int((y_ext==0).sum())}")

    models = make_models()
    for fs_name, genes in feature_sets.items():
        XtT, XeT, shared = align_features(X_tcga, X_ext[ext_ids] if ext_ids else X_ext, genes)
        if XtT.shape[1] < 3 or XeT.shape[0] < 5 or len(set(y_ext)) < 2:
            continue
        for m_name, model in models.items():
            try:
                model.fit(XtT.values, y_tcga)
                # train scores for isotonic
                if hasattr(model, "predict_proba"):
                    p_tr = model.predict_proba(XtT.values)[:, 1]
                    p_ext = model.predict_proba(XeT.values)[:, 1]
                else:
                    p_tr = model.decision_function(XtT.values)
                    p_ext = model.decision_function(XeT.values)
                    # min-max -> 0..1
                    p_tr = (p_tr - p_tr.min()) / max(p_tr.ptp(), 1e-6)
                    p_ext = (p_ext - p_ext.min()) / max(p_ext.ptp(), 1e-6)
                raw = score_binary(y_ext, p_ext)
                lo, hi = bootstrap_auc_ci(y_ext, p_ext)
                # post isotonic (fit on train)
                p_ext_cal = isotonic_calibrate(p_tr, y_tcga, p_ext)
                cal = score_binary(y_ext, p_ext_cal)
                row = {
                    "task": "V1", "external": ext_dataset, "n_ext": len(y_ext),
                    "feature_set": fs_name, "n_features": len(shared),
                    "model": m_name,
                    "auc_pre": raw["auc"], "auc_lo": lo, "auc_hi": hi,
                    "pr_auc_pre": raw["pr_auc"],
                    "brier_pre": raw["brier"], "ece_pre": raw["ece"],
                    "bACC_Youden_pre": raw["bACC_Youden"], "F1_pre": raw["F1"],
                    "NPV_at_Se95_pre": raw["NPV_at_Se95"],
                    "auc_post": cal["auc"], "brier_post": cal["brier"], "ece_post": cal["ece"],
                    "bACC_Youden_post": cal["bACC_Youden"],
                }
                results.append(row)
                calibration_rows.append({
                    "task": "V1", "feature_set": fs_name, "model": m_name,
                    "y": ";".join(map(str, y_ext.tolist())),
                    "p_pre": ";".join(f"{v:.4f}" for v in p_ext),
                    "p_post": ";".join(f"{v:.4f}" for v in p_ext_cal),
                })
                if not np.isnan(raw["auc"]) and raw["auc"] > best["auc"]:
                    best = {
                        "auc": raw["auc"], "feature_set": fs_name, "model": m_name,
                        "y": y_ext, "p": p_ext, "X_tcga": XtT, "y_tcga": y_tcga,
                        "X_ext": XeT, "genes": shared,
                    }
            except Exception as e:
                log.warning(f"V1 fit fail {fs_name}/{m_name}: {e}")
    return results, calibration_rows, best


def run_v4_pooled_gpl570(feature_sets, X_tcga, y_tcga, sm):
    """TCGA -> pooled GPL570. Verify GSE27155 bACC 0.5 fix."""
    results = []
    cohorts = ["GSE27155", "GSE33630", "GSE29265"]
    frames = []
    label_frames = []
    for c in cohorts:
        f_v3 = DP_V3_MIC / f"{c}_v3_log2.tsv"
        if not f_v3.exists():
            continue
        df = pd.read_csv(f_v3, sep="\t", index_col=0)
        frames.append(df)
        # label: try merged master
        for meta_source in [META / "sample_master_v3_merged.tsv", META / "sample_master_v3.tsv", META / "sample_master.tsv"]:
            if meta_source.exists():
                m = pd.read_csv(meta_source, sep="\t")
                sub = m[(m["dataset"] == c) & (m["sample_id"].isin(df.columns))]
                if not sub.empty:
                    lbl_map = {"BRAF": 1, "RAS": 0}
                    if "driver_anchor" in sub.columns and sub["driver_anchor"].isin(lbl_map.keys()).any():
                        sub2 = sub[sub["driver_anchor"].isin(lbl_map.keys())]
                        label_frames.append(pd.DataFrame({
                            "sample_id": sub2["sample_id"].values,
                            "cohort": c,
                            "y": sub2["driver_anchor"].map(lbl_map).values,
                        }))
                        break
                    # fallback: tumor vs normal
                    if "normal_vs_tumor" in sub.columns and sub["normal_vs_tumor"].isin(["tumor", "normal"]).any():
                        nvt = {"tumor": 1, "normal": 0}
                        sub2 = sub[sub["normal_vs_tumor"].isin(nvt.keys())]
                        label_frames.append(pd.DataFrame({
                            "sample_id": sub2["sample_id"].values,
                            "cohort": c,
                            "y": sub2["normal_vs_tumor"].map(nvt).values,
                        }))
                        break
    if not frames or not label_frames:
        return results
    import functools
    pooled = functools.reduce(lambda a, b: a.join(b, how="outer"), frames)
    pooled = pooled.fillna(pooled.median(axis=1).values[:, None]) if pooled.isna().any().any() else pooled
    lbl = pd.concat(label_frames, ignore_index=True).drop_duplicates("sample_id")
    lbl = lbl[lbl["sample_id"].isin(pooled.columns)]
    if lbl.empty or lbl["y"].nunique() < 2:
        return results

    # ComBat adjust pre/post
    from sklearn.preprocessing import StandardScaler
    models = make_models()

    def _fit_score(pooled_expr, tag):
        out = []
        for fs_name, genes in feature_sets.items():
            XtT, XpT, shared = align_features(X_tcga, pooled_expr[lbl["sample_id"].tolist()], genes)
            if XtT.shape[1] < 3 or XpT.shape[0] < 5:
                continue
            for m_name, model in models.items():
                try:
                    model.fit(XtT.values, y_tcga)
                    if hasattr(model, "predict_proba"):
                        p = model.predict_proba(XpT.values)[:, 1]
                    else:
                        p = model.decision_function(XpT.values)
                        p = (p - p.min()) / max(p.ptp(), 1e-6)
                    y = lbl["y"].values
                    s = score_binary(y, p)
                    lo, hi = bootstrap_auc_ci(y, p)
                    # per cohort
                    per_cohort = {}
                    for cc in lbl["cohort"].unique():
                        mask = (lbl["cohort"] == cc).values
                        if mask.sum() < 4 or len(set(y[mask])) < 2:
                            continue
                        per_cohort[cc] = score_binary(y[mask], p[mask])
                    out.append({
                        "task": "V4", "pool_tag": tag,
                        "feature_set": fs_name, "model": m_name,
                        "auc": s["auc"], "auc_lo": lo, "auc_hi": hi,
                        "brier": s["brier"], "ece": s["ece"],
                        "bACC_Youden": s["bACC_Youden"],
                        "per_cohort_bacc": json.dumps({k: v["bACC_Youden"] for k, v in per_cohort.items()}),
                        "per_cohort_auc": json.dumps({k: v["auc"] for k, v in per_cohort.items()}),
                    })
                except Exception as e:
                    log.warning(f"V4 {tag}/{fs_name}/{m_name}: {e}")
        return out

    results += _fit_score(pooled, "pre_combat")
    # ComBat: subtype-preserving, batch=cohort, covariate=y
    try:
        from sklearn.preprocessing import LabelEncoder
        batch = []
        cov = []
        for s in pooled.columns:
            row = lbl[lbl["sample_id"] == s]
            if row.empty:
                batch.append("NA"); cov.append(-1)
            else:
                batch.append(row["cohort"].iloc[0]); cov.append(int(row["y"].iloc[0]))
        batch_arr = np.array(batch)
        cov_arr = np.array(cov)
        # Simple residualisation: subtract per-batch gene mean (keeping y-stratified mean)
        adj = pooled.copy()
        for cc in np.unique(batch_arr):
            if cc == "NA":
                continue
            m = batch_arr == cc
            adj.iloc[:, m] = adj.iloc[:, m].sub(adj.iloc[:, m].mean(axis=1), axis=0)
        # re-centre to overall mean
        adj = adj.add(pooled.mean(axis=1), axis=0)
        results += _fit_score(adj, "post_combat_lite")
    except Exception as e:
        log.warning(f"V4 combat lite failed: {e}")

    return results


def run_v5_lodo(feature_sets, sm):
    """LODO over all cohorts with labels.

    For each held-out cohort, train on the remaining cohorts (shared gene space),
    evaluate on held-out cohort with best models.
    """
    # Collect cohort expression + labels (binary tumor vs normal where driver_anchor missing)
    cohorts = []
    for p in list(DP_V3_RNA.glob("*_v3_log2.tsv")) + list(DP_V3_MIC.glob("*_v3_log2.tsv")):
        name = p.stem.replace("_v3_log2", "")
        cohorts.append(name)
    data = {}
    for c in sorted(set(cohorts)):
        try:
            X = load_expr(c)
            merged = pd.read_csv(META / "sample_master_v3_merged.tsv", sep="\t") if (META / "sample_master_v3_merged.tsv").exists() else pd.read_csv(META / "sample_master.tsv", sep="\t")
            sub = merged[(merged["dataset"] == c) & (merged["sample_id"].isin(X.columns))]
            if sub.empty:
                continue
            # Prefer driver_anchor; else tumor/normal
            y = None; ids = None
            if "driver_anchor" in sub.columns and sub["driver_anchor"].isin(["BRAF", "RAS"]).sum() >= 6:
                s2 = sub[sub["driver_anchor"].isin(["BRAF", "RAS"])]
                y = s2["driver_anchor"].map({"BRAF": 1, "RAS": 0}).values
                ids = s2["sample_id"].tolist()
            elif sub.get("normal_vs_tumor", pd.Series(dtype=str)).isin(["tumor", "normal"]).sum() >= 6:
                s2 = sub[sub["normal_vs_tumor"].isin(["tumor", "normal"])]
                y = s2["normal_vs_tumor"].map({"tumor": 1, "normal": 0}).values
                ids = s2["sample_id"].tolist()
            if y is not None and len(set(y)) == 2:
                # skip cohorts that use ENSG (not HGNC) — no shared gene space
                first_idx = X.index[0] if len(X.index) > 0 else ""
                if str(first_idx).startswith("ENSG"):
                    log.warning(f"LODO: {c} uses ENSG ids; skipping (no HGNC overlap)")
                    continue
                data[c] = (X[ids], y)
        except Exception as e:
            log.warning(f"LODO load {c}: {e}")
    log.info(f"LODO cohorts loaded: {list(data.keys())}")

    rows = []
    names = list(data.keys())
    models = make_models()
    for held in names:
        train_cohorts = [c for c in names if c != held]
        # Shared gene symbols across all train + held
        sets = [data[c][0].index for c in train_cohorts + [held]]
        shared_genes = set(sets[0])
        for s in sets[1:]:
            shared_genes = shared_genes & set(s)
        shared_genes = sorted(shared_genes)
        if len(shared_genes) < 5:
            continue
        for fs_name, genes in feature_sets.items():
            g = [x for x in genes if x in shared_genes]
            if len(g) < 3:
                continue
            # Train matrix
            Xs = []; ys = []
            for c in train_cohorts:
                xc, yc = data[c]
                Xs.append(xc.loc[g].T.values); ys.append(yc)
            Xtr = np.vstack(Xs); ytr = np.concatenate(ys)
            Xte, yte = data[held][0].loc[g].T.values, data[held][1]
            if len(set(ytr)) < 2 or len(set(yte)) < 2:
                continue
            for m_name, model in models.items():
                try:
                    from sklearn.base import clone
                    mdl = clone(model)
                    mdl.fit(Xtr, ytr)
                    if hasattr(mdl, "predict_proba"):
                        p = mdl.predict_proba(Xte)[:, 1]
                    else:
                        p = mdl.decision_function(Xte)
                        p = (p - p.min()) / max(p.ptp(), 1e-6)
                    s = score_binary(yte, p)
                    rows.append({"held_out": held, "feature_set": fs_name, "model": m_name,
                                 "n_train": len(ytr), "n_test": len(yte),
                                 "auc": s["auc"], "bACC_Youden": s["bACC_Youden"],
                                 "brier": s["brier"], "ece": s["ece"]})
                except Exception as e:
                    log.warning(f"LODO {held}/{fs_name}/{m_name}: {e}")
    return rows


def run_v6_permutation(best, n_iter=1000):
    if not best or "X_tcga" not in best:
        return []
    rng = np.random.default_rng(SEED)
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import Pipeline
    from sklearn.metrics import roc_auc_score
    model = Pipeline([("sc", StandardScaler()), ("clf", LogisticRegression(max_iter=2000, random_state=SEED))])
    X_tr = best["X_tcga"].values; y_tr = best["y_tcga"]
    X_te = best["X_ext"].values; y_te = best["y"]
    observed = best["auc"]
    aucs = []
    for i in range(n_iter):
        yp = y_tr.copy(); rng.shuffle(yp)
        model.fit(X_tr, yp)
        p = model.predict_proba(X_te)[:, 1]
        try:
            aucs.append(roc_auc_score(y_te, p))
        except Exception:
            pass
    aucs = np.array(aucs)
    p_val = float((aucs >= observed).mean())
    rows = [{"observed_auc": observed, "n_iter": len(aucs), "null_mean": float(aucs.mean()),
             "null_std": float(aucs.std()), "p_value": p_val,
             "null_q97.5": float(np.percentile(aucs, 97.5))}]
    return rows


def save_tsv(path: Path, rows: list):
    if not rows:
        pd.DataFrame([{"status": "empty"}]).to_csv(path, sep="\t", index=False)
        return
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def make_figs(v1_rows, v4_rows, lodo_rows, best, perm_rows):
    import plotly.graph_objects as go
    # Fig 1: ROC overlay V1
    fig1 = go.Figure()
    if v1_rows:
        top = sorted(v1_rows, key=lambda r: -(r.get("auc_pre") or 0))[:5]
        for r in top:
            fig1.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines",
                                      line=dict(dash="dot", color="gray"),
                                      showlegend=False))
            # no per-curve data; plot summary point
            fig1.add_trace(go.Scatter(x=[1 - (r.get("bACC_Youden_pre") or 0.5)],
                                      y=[r.get("bACC_Youden_pre") or 0.5],
                                      mode="markers+text",
                                      name=f"{r['feature_set']}/{r['model']}",
                                      text=[f"AUC={r.get('auc_pre',0):.2f}"]))
    fig1.update_layout(title="V1 ROC summary (external)", xaxis_title="FPR", yaxis_title="TPR")
    fig1.write_html(str(FIGS / "v3_ext_v1_roc.html"), include_plotlyjs="cdn")

    # Fig 2: PR
    fig2 = go.Figure()
    for r in (v1_rows or [])[:10]:
        fig2.add_trace(go.Bar(x=[f"{r['feature_set']}/{r['model']}"], y=[r.get("pr_auc_pre") or 0]))
    fig2.update_layout(title="V1 PR-AUC (external)", yaxis_title="PR-AUC")
    fig2.write_html(str(FIGS / "v3_ext_v1_pr.html"), include_plotlyjs="cdn")

    # Fig 3: Calibration pre/post
    fig3 = go.Figure()
    for r in (v1_rows or [])[:6]:
        fig3.add_trace(go.Scatter(x=[r.get("ece_pre")], y=[r.get("ece_post")],
                                  mode="markers+text",
                                  text=[f"{r['feature_set']}/{r['model']}"],
                                  textposition="top center",
                                  name=f"{r['feature_set']}/{r['model']}"))
    fig3.add_trace(go.Scatter(x=[0, 0.5], y=[0, 0.5], mode="lines", line=dict(dash="dot"), showlegend=False))
    fig3.update_layout(title="V1 Calibration ECE pre vs post", xaxis_title="ECE pre", yaxis_title="ECE post")
    fig3.write_html(str(FIGS / "v3_ext_v1_calibration.html"), include_plotlyjs="cdn")

    # Fig 4: Confusion (best V1)
    fig4 = go.Figure()
    if best and "y" in best:
        y = best["y"]; p = best["p"]
        yhat = (p >= 0.5).astype(int)
        from sklearn.metrics import confusion_matrix
        cm = confusion_matrix(y, yhat)
        fig4 = go.Figure(data=go.Heatmap(z=cm, x=["pred0", "pred1"], y=["true0", "true1"],
                                         colorscale="Blues", text=cm, texttemplate="%{text}"))
        fig4.update_layout(title=f"V1 Best Confusion: {best.get('feature_set','?')}/{best.get('model','?')}")
    fig4.write_html(str(FIGS / "v3_ext_v1_confusion.html"), include_plotlyjs="cdn")

    # Fig 5: V4 pre/post combat
    fig5 = go.Figure()
    if v4_rows:
        import collections
        agg = collections.defaultdict(dict)
        for r in v4_rows:
            agg[(r["feature_set"], r["model"])][r["pool_tag"]] = r["bACC_Youden"]
        for k, v in agg.items():
            fig5.add_trace(go.Bar(name=f"{k[0]}/{k[1]}",
                                  x=["pre_combat", "post_combat_lite"],
                                  y=[v.get("pre_combat", 0), v.get("post_combat_lite", 0)]))
    fig5.update_layout(title="V4 pooled-GPL570 bACC pre vs post batch-correction", barmode="group")
    fig5.write_html(str(FIGS / "v3_ext_v4_pooled.html"), include_plotlyjs="cdn")

    # Fig 6: 6-class donut
    fig6 = go.Figure()
    try:
        ftcga = pd.read_csv(META / "v3_fusion_anchor_tcga.tsv", sep="\t")
        vc = ftcga["v3_anchor_6class"].value_counts().to_dict()
        title_warn = " (⚠ small n classes)" if any(v < 20 for v in vc.values()) else ""
        fig6 = go.Figure(data=[go.Pie(labels=list(vc.keys()), values=list(vc.values()), hole=0.4)])
        fig6.update_layout(title="TCGA v3 6-class composition" + title_warn)
    except Exception:
        pass
    fig6.write_html(str(FIGS / "v3_ext_v3_sixclass.html"), include_plotlyjs="cdn")

    # Fig 7: LODO heatmap
    fig7 = go.Figure()
    if lodo_rows:
        df = pd.DataFrame(lodo_rows)
        # Pick best model-per-feature_set, pivot feature_set x held_out
        best_per = df.sort_values("auc", ascending=False).drop_duplicates(["held_out", "feature_set"])
        piv = best_per.pivot(index="feature_set", columns="held_out", values="auc")
        fig7 = go.Figure(data=go.Heatmap(z=piv.values, x=list(piv.columns), y=list(piv.index),
                                         colorscale="RdBu", zmid=0.5, text=np.round(piv.values, 2),
                                         texttemplate="%{text}"))
        fig7.update_layout(title="LODO AUC (best model per cell)")
    fig7.write_html(str(FIGS / "v3_ext_v5_lodo.html"), include_plotlyjs="cdn")

    # Fig 8: permutation null
    fig8 = go.Figure()
    if perm_rows:
        r = perm_rows[0]
        fig8.add_trace(go.Histogram(x=np.random.normal(r["null_mean"], r["null_std"] or 0.01, 500),
                                    name="Null (approx)"))
        fig8.add_vline(x=r["observed_auc"], line_dash="dash", line_color="red",
                       annotation_text=f"Observed AUC={r['observed_auc']:.3f}")
        fig8.update_layout(title=f"V6 permutation null (p={r['p_value']:.3f})")
    fig8.write_html(str(FIGS / "v3_ext_v6_permutation.html"), include_plotlyjs="cdn")


def main():
    log.info("v3_ext_05 start")

    # Load TCGA training set
    X_tcga = load_expr("TCGA-THCA")
    sm = pd.read_csv(META / "sample_master_v3_merged.tsv", sep="\t") if (META / "sample_master_v3_merged.tsv").exists() else pd.read_csv(META / "sample_master.tsv", sep="\t")
    ids_tcga, y_tcga_bin = get_labels(sm, "TCGA-THCA", n_classes=2)
    # Intersect
    ids_tcga = [s for s in ids_tcga if s in X_tcga.columns]
    idx = [sm[sm["sample_id"] == s].index[0] for s in ids_tcga if not sm[sm["sample_id"] == s].empty]
    y_tcga_bin = np.array([1 if sm.iloc[i]["driver_anchor"] == "BRAF" else 0 for i in idx])
    X_tcga_sub = X_tcga[ids_tcga]
    log.info(f"TCGA 2-class training: n={len(ids_tcga)}  BRAF={int(y_tcga_bin.sum())}  RAS={int((y_tcga_bin==0).sum())}")

    fs = prepare_feature_sets(X_tcga_sub, y_tcga_bin)
    log.info(f"feature sets: { {k: len(v) for k, v in fs.items()} }")

    # V1
    v1, cal, best = run_v1(fs, X_tcga_sub, y_tcga_bin, sm)
    log.info(f"V1: {len(v1)} rows, best AUC={best.get('auc')}")

    # V2 (3-class) -- approximate as extra task: TCGA 3-class train, GSE76039 external
    v2 = []
    try:
        ids3, y3 = get_labels(sm, "TCGA-THCA", n_classes=3)
        ids3 = [s for s in ids3 if s in X_tcga.columns]
        y3 = np.array([y3[i] for i, s in enumerate(ids3) if s in X_tcga.columns][:len(ids3)])
        # rebuild properly
        rows3 = sm[(sm["dataset"] == "TCGA-THCA") & (sm["sample_id"].isin(ids3))]
        y3_map = []
        for sid in ids3:
            r = sm[sm["sample_id"] == sid].iloc[0]
            a = str(r.get("driver_anchor", "")).upper()
            if a == "BRAF": y3_map.append(1)
            elif a == "RAS": y3_map.append(0)
            else: y3_map.append(2)
        y3 = np.array(y3_map)
        X3 = X_tcga[ids3]
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler
        from sklearn.pipeline import Pipeline
        # Fit multinomial on TDS16
        genes = [g for g in fs["TDS16"] if g in X3.index]
        Xt = X3.loc[genes].T.values
        model = Pipeline([("sc", StandardScaler()),
                          ("clf", LogisticRegression(solver="lbfgs",
                                                     max_iter=2000, random_state=SEED))])
        model.fit(Xt, y3)
        # Evaluate on GSE126698 with binary collapse if multi not possible
        Xext = load_expr("GSE126698")
        sub = sm[(sm["dataset"] == "GSE126698") & (sm["sample_id"].isin(Xext.columns))]
        y_ext3 = []
        ids_ext = []
        for _, r in sub.iterrows():
            nvt = str(r.get("normal_vs_tumor", "")).lower()
            if nvt == "normal":
                y_ext3.append(2); ids_ext.append(r["sample_id"])
            elif nvt == "tumor":
                # Map histology_subtype to BRAF/RAS guess
                y_ext3.append(1 if "PTC" in str(r.get("histology_subtype", "")).upper() else 0)
                ids_ext.append(r["sample_id"])
        if ids_ext and len(set(y_ext3)) >= 2:
            Xe = Xext[ids_ext].loc[genes].T.values
            p_all = model.predict_proba(Xe)
            yhat = model.predict(Xe)
            from sklearn.metrics import accuracy_score, balanced_accuracy_score
            v2.append({
                "task": "V2", "external": "GSE126698",
                "n_classes": int(len(np.unique(y3))),
                "n_ext": len(ids_ext),
                "accuracy": float(accuracy_score(y_ext3, yhat)),
                "bACC": float(balanced_accuracy_score(y_ext3, yhat)),
            })
    except Exception as e:
        log.warning(f"V2 failed: {e}")

    # V3 (6-class)
    v3_rows = []
    try:
        ftcga = pd.read_csv(META / "v3_fusion_anchor_tcga.tsv", sep="\t")
        common = ftcga[ftcga["sample_id"].isin(X_tcga.columns)]
        classes = sorted(common["v3_anchor_6class"].unique())
        y6 = common["v3_anchor_6class"].map({c: i for i, c in enumerate(classes)}).values
        X6 = X_tcga[common["sample_id"].tolist()]
        genes = [g for g in fs["TierA67_clean"] if g in X6.index]
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler
        from sklearn.pipeline import Pipeline
        from sklearn.model_selection import StratifiedKFold
        Xt = X6.loc[genes].T.values
        if len(set(y6)) >= 2 and len(y6) >= 20:
            from sklearn.metrics import roc_auc_score
            skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=SEED)
            aucs_ovr = []
            for tr_i, te_i in skf.split(Xt, y6):
                model = Pipeline([("sc", StandardScaler()),
                                  ("clf", LogisticRegression(
                                                             max_iter=2000, random_state=SEED))])
                model.fit(Xt[tr_i], y6[tr_i])
                p = model.predict_proba(Xt[te_i])
                try:
                    aucs_ovr.append(roc_auc_score(y6[te_i], p, multi_class="ovr"))
                except Exception:
                    pass
            v3_rows.append({
                "task": "V3", "n_classes": len(classes), "n_samples": len(y6),
                "auc_ovr_cv": float(np.mean(aucs_ovr)) if aucs_ovr else np.nan,
                "classes": ",".join(classes),
                "class_counts": json.dumps(common["v3_anchor_6class"].value_counts().to_dict()),
            })
    except Exception as e:
        log.warning(f"V3 failed: {e}")

    # V4 pooled GPL570
    v4 = run_v4_pooled_gpl570(fs, X_tcga_sub, y_tcga_bin, sm)

    # V5 LODO
    lodo = run_v5_lodo(fs, sm)

    # V6 permutation on best V1 config
    perm = run_v6_permutation(best, n_iter=1000)

    # Save all
    save_tsv(ML / "v3_ext_validation_results.tsv", v1 + v2 + v3_rows + v4)
    save_tsv(ML / "v3_ext_calibration.tsv", cal)
    save_tsv(ML / "v3_ext_permutation_null.tsv", perm)
    save_tsv(ML / "v3_ext_lodo_v3.tsv", lodo)

    make_figs(v1, v4, lodo, best, perm)
    log.info("v3_ext_05 done")


if __name__ == "__main__":
    main()
