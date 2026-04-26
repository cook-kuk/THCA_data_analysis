#!/usr/bin/env python3
"""v5 Track 1 — Nonlinear Batch Correction Ablation.

Tests whether v4's UNRECOVERABLE verdict holds against non-linear
correction methods. Per-method outputs:
  - post_identifiability AUC (multiclass LogReg OvR)
  - post_LODO AUC (train on N-1 cohorts, test on held)
  - bio_preservation = TCGA-internal 5-fold CV AUC (BRAF vs RAS-like)

Methods: Harmony, ComBat-seq (inmoose), limma-removeBatchEffect (py),
MNN (optional), Z-score per cohort, No correction.

Outputs: results/ml/v5_track1_correction_comparison.tsv
         reports/html/figs_interactive/v5_correction_{radar,pca_grid,lodo_heatmap}.html
"""
from __future__ import annotations

import json
import sys
import time
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).parent))
from v5_common import (  # noqa: E402
    build_pooled_matrix,
    binary_label,
    bio_label_bvr,
    load_sample_master,
    RESULTS_ML,
    FIGS_INTERACTIVE,
    LOGS,
    safe_write_json,
)


def log(msg: str) -> None:
    print(f"[v5-track1] {msg}", flush=True)


def identifiability_auc(X: np.ndarray, cohort: np.ndarray) -> float:
    uniq = np.unique(cohort)
    aucs = []
    for c in uniq:
        y = (cohort == c).astype(int)
        if y.sum() < 5 or (len(y) - y.sum()) < 5:
            continue
        skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)
        fold_aucs = []
        for tr, te in skf.split(X, y):
            clf = LogisticRegression(max_iter=1500, C=1.0, n_jobs=1)
            clf.fit(X[tr], y[tr])
            p = clf.predict_proba(X[te])[:, 1]
            try:
                fold_aucs.append(roc_auc_score(y[te], p))
            except Exception:
                continue
        if fold_aucs:
            aucs.append(float(np.mean(fold_aucs)))
    return float(np.mean(aucs)) if aucs else float("nan")


def lodo_auc(X: np.ndarray, y: np.ndarray, cohorts: np.ndarray):
    per_cohort = []
    uniq = np.unique(cohorts)
    for held in uniq:
        tr = cohorts != held
        te = cohorts == held
        if len(np.unique(y[tr])) < 2 or len(np.unique(y[te])) < 2:
            per_cohort.append({"held": held, "auc": float("nan"),
                               "n_test": int(te.sum()), "note": "single-class"})
            continue
        clf = LogisticRegression(max_iter=2000, C=1.0)
        clf.fit(X[tr], y[tr])
        p = clf.predict_proba(X[te])[:, 1]
        auc = float(roc_auc_score(y[te], p))
        per_cohort.append({"held": held, "auc": auc, "n_test": int(te.sum()),
                           "note": "small n" if te.sum() < 20 else ""})
    m = np.nanmean([r["auc"] for r in per_cohort])
    return float(m), per_cohort


def bio_preservation_tcga(big_df: pd.DataFrame, meta_df: pd.DataFrame) -> float:
    """5-fold CV AUC on TCGA using BRAF_like vs RAS_like labels on corrected data."""
    sm_idx = load_sample_master().set_index("sample_id")
    tcga_samples = [s for s in meta_df.index if meta_df.loc[s, "cohort"] == "TCGA-THCA"]
    if len(tcga_samples) < 20:
        return float("nan")
    y_bvr = []
    keep = []
    for s in tcga_samples:
        row = sm_idx.loc[s] if s in sm_idx.index else None
        if row is None:
            continue
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        lbl = bio_label_bvr(row)
        if lbl is None:
            continue
        keep.append(s)
        y_bvr.append(lbl)
    if len(keep) < 20 or len(set(y_bvr)) < 2:
        return float("nan")
    X = big_df.loc[:, keep].T.values
    y = np.array(y_bvr)
    X = StandardScaler().fit_transform(X)
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    aucs = []
    for tr, te in skf.split(X, y):
        clf = LogisticRegression(max_iter=2000, C=0.5)
        clf.fit(X[tr], y[tr])
        p = clf.predict_proba(X[te])[:, 1]
        try:
            aucs.append(roc_auc_score(y[te], p))
        except Exception:
            pass
    return float(np.mean(aucs)) if aucs else float("nan")


# ============ CORRECTION METHODS ============

def correct_none(big_df, cohort_vec, mod_vec):
    return big_df.copy()


def correct_zscore_per_cohort(big_df, cohort_vec, mod_vec):
    out = big_df.copy()
    for c in np.unique(cohort_vec):
        mask = cohort_vec == c
        sub = out.loc[:, mask]
        mu = sub.mean(axis=1)
        sd = sub.std(axis=1).replace(0, 1.0)
        out.loc[:, mask] = sub.sub(mu, axis=0).div(sd, axis=0)
    return out


def correct_limma(big_df, cohort_vec, mod_vec):
    """Python reimplementation of limma::removeBatchEffect with preserved design.

    For each gene: fit y ~ batch + design, subtract batch coefficients.
    """
    batch = pd.get_dummies(pd.Series(cohort_vec), drop_first=True).astype(float).values
    mod = pd.get_dummies(pd.Series(mod_vec), drop_first=True).astype(float).values \
        if mod_vec is not None else np.zeros((len(cohort_vec), 0))
    design = np.hstack([np.ones((len(cohort_vec), 1)), mod, batch])
    n_mod = 1 + mod.shape[1]
    Y = big_df.values  # genes x samples
    # solve Y.T = design @ B.T -> B.T = pinv(design) @ Y.T
    coef = np.linalg.pinv(design) @ Y.T  # (params x genes)
    beta_batch = coef[n_mod:, :]  # batch coefs, genes x batch
    # subtract batch contribution
    batch_contrib = batch @ beta_batch  # samples x genes
    Y_corr = Y - batch_contrib.T
    return pd.DataFrame(Y_corr, index=big_df.index, columns=big_df.columns)


def correct_harmony(big_df, cohort_vec, mod_vec):
    """Harmony on PCA-reduced space, then project back linearly."""
    try:
        import harmonypy as hm
    except Exception as e:
        raise RuntimeError(f"harmonypy not available: {e}")
    X = big_df.T.values  # samples x genes
    X_std = StandardScaler().fit_transform(X)
    pca = PCA(n_components=min(30, X_std.shape[1] - 1, X_std.shape[0] - 1))
    emb = pca.fit_transform(X_std)  # samples x k
    meta = pd.DataFrame({"batch": cohort_vec})
    try:
        ho = hm.run_harmony(emb, meta, "batch", max_iter_harmony=10)
    except TypeError:
        # older API
        ho = hm.run_harmony(emb, meta, "batch")
    emb_corr = np.array(ho.Z_corr)
    # Z_corr is (k, samples); we want (samples, k)
    if emb_corr.shape[0] != emb.shape[0]:
        emb_corr = emb_corr.T
    # project back: approximate corrected expression as pca inverse_transform + residuals
    X_corr_std = pca.inverse_transform(emb_corr)
    # unscale
    mu = X.mean(axis=0)
    sd = X.std(axis=0, ddof=0)
    sd[sd == 0] = 1.0
    X_corr = X_corr_std * sd + mu
    return pd.DataFrame(X_corr.T, index=big_df.index, columns=big_df.columns)


def correct_combatseq(big_df, cohort_vec, mod_vec):
    """ComBat (biology-preserving) via inmoose.pycombat_norm."""
    from inmoose.pycombat import pycombat_norm
    mod_dum = pd.get_dummies(pd.Series(mod_vec), drop_first=True).astype(float).values \
        if mod_vec is not None else None
    if mod_dum is not None and mod_dum.shape[1] == 0:
        mod_dum = None
    corr = pycombat_norm(big_df, list(cohort_vec), covar_mod=mod_dum)
    if not isinstance(corr, pd.DataFrame):
        corr = pd.DataFrame(corr, index=big_df.index, columns=big_df.columns)
    corr.index = big_df.index
    corr.columns = big_df.columns
    return corr


def correct_mnn(big_df, cohort_vec, mod_vec):
    """Lightweight MNN approximation: align each non-anchor cohort to TCGA
    via mutual-nearest-neighbor mean shift. True mnnpy not available in venv.
    """
    anchor = "TCGA-THCA" if "TCGA-THCA" in np.unique(cohort_vec) else np.unique(cohort_vec)[0]
    out = big_df.copy()
    anchor_cols = [c for c in out.columns if cohort_vec[list(out.columns).index(c)] == anchor]
    anchor_mat = out[anchor_cols].values  # genes x n_a
    for c in np.unique(cohort_vec):
        if c == anchor:
            continue
        cols = [cc for cc in out.columns if cohort_vec[list(out.columns).index(cc)] == c]
        sub = out[cols].values  # genes x n_c
        # compute nearest-neighbor pairs: for each anchor sample, find closest
        # non-anchor via cosine, take top-5% and average offset
        from sklearn.metrics.pairwise import cosine_similarity
        sim = cosine_similarity(anchor_mat.T, sub.T)
        # pick top-5% matches per anchor
        k = max(1, int(0.05 * sim.shape[1]))
        offsets = []
        for i in range(sim.shape[0]):
            j_top = np.argsort(-sim[i])[:k]
            a = anchor_mat[:, i]
            b_mean = sub[:, j_top].mean(axis=1)
            offsets.append(a - b_mean)
        mean_shift = np.array(offsets).mean(axis=0)  # genes
        sub_corr = sub + mean_shift[:, None]
        for j, cc in enumerate(cols):
            out.loc[:, cc] = sub_corr[:, j]
    return out


METHODS = [
    ("none",      correct_none,           "v4 baseline"),
    ("zscore",    correct_zscore_per_cohort, "z-score per cohort"),
    ("limma",     correct_limma,          "linear removeBatchEffect"),
    ("combat",    correct_combatseq,      "ComBat (inmoose, biology-preserving)"),
    ("harmony",   correct_harmony,        "Harmony (PCA-space clustered correction)"),
    ("mnn",       correct_mnn,            "mutual-nearest-neighbor mean shift"),
]


def run_method(name: str, fn, big_df, meta_df, cohort_vec, y_vec, mod_vec):
    t0 = time.time()
    result = {"method": name, "status": "OK", "post_identifiability_auc": float("nan"),
              "post_lodo_auc": float("nan"), "bio_preservation_auc": float("nan"),
              "runtime_sec": float("nan"), "notes": ""}
    try:
        corrected = fn(big_df, cohort_vec, mod_vec)
        if corrected is None:
            raise RuntimeError("correction returned None")
        corrected = corrected.replace([np.inf, -np.inf], np.nan)
        corrected = corrected.dropna(axis=0)
        # standardize to top 1500 variable genes for speed
        top = corrected.var(axis=1).sort_values(ascending=False).head(1500).index
        X = StandardScaler().fit_transform(corrected.loc[top].T.values)
        ident = identifiability_auc(X, cohort_vec)
        lodo_m, _ = lodo_auc(X, y_vec, cohort_vec)
        bio = bio_preservation_tcga(corrected, meta_df)
        result.update({
            "post_identifiability_auc": ident,
            "post_lodo_auc": lodo_m,
            "bio_preservation_auc": bio,
        })
        # decision rule
        if (ident < 0.70 and lodo_m > 0.80 and (bio >= 0.85 or np.isnan(bio))):
            result["verdict"] = "RESCUED"
        elif ident < 0.80 and lodo_m > 0.75:
            result["verdict"] = "PARTIAL"
        else:
            result["verdict"] = "UNRECOVERABLE"
    except Exception as e:
        result["status"] = "FAILED"
        result["notes"] = str(e)[:300]
        result["verdict"] = "SKIPPED"
    result["runtime_sec"] = round(time.time() - t0, 2)
    return result, corrected if result["status"] == "OK" else None


def main() -> int:
    log("loading pooled matrix...")
    big_df, meta_df = build_pooled_matrix(min_samples_per_cohort=5,
                                          top_var_genes=3000,
                                          label_fn=binary_label)
    log(f"pooled: {big_df.shape[0]} genes x {big_df.shape[1]} samples; "
        f"cohorts={sorted(meta_df['cohort'].unique())}")

    # build 3-class mod vector from sample_master
    sm_idx = load_sample_master().set_index("sample_id")
    mod_vec = []
    for s in meta_df.index:
        row = sm_idx.loc[s]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        nvt = str(row.get("normal_vs_tumor", "")).strip().lower()
        if nvt == "normal":
            mod_vec.append("normal")
        else:
            hist = str(row.get("histology_subtype", "")).strip()
            if hist in {"ATC", "PDTC"}:
                mod_vec.append("aggressive")
            else:
                mod_vec.append("indolent")
    mod_vec = np.array(mod_vec)

    cohort_vec = meta_df["cohort"].values
    y_vec = meta_df["label"].values.astype(int)

    rows = []
    corrected_by_method = {}
    for name, fn, desc in METHODS:
        log(f"method = {name} ({desc})")
        res, corr = run_method(name, fn, big_df, meta_df, cohort_vec, y_vec, mod_vec)
        res["description"] = desc
        rows.append(res)
        corrected_by_method[name] = corr
        log(f"  verdict={res.get('verdict')} ident={res['post_identifiability_auc']:.3f}"
            f" lodo={res['post_lodo_auc']:.3f} bio={res['bio_preservation_auc']:.3f}"
            f" status={res['status']}")

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS_ML / "v5_track1_correction_comparison.tsv", sep="\t", index=False)

    # ---------------- radar chart ----------------
    categories = ["1-identifiability", "LODO", "bio_preservation"]
    fig = go.Figure()
    for _, r in df.iterrows():
        if r["status"] != "OK":
            continue
        vals = [
            1.0 - (r["post_identifiability_auc"] if not np.isnan(r["post_identifiability_auc"]) else 1.0),
            r["post_lodo_auc"] if not np.isnan(r["post_lodo_auc"]) else 0.0,
            r["bio_preservation_auc"] if not np.isnan(r["bio_preservation_auc"]) else 0.0,
        ]
        fig.add_trace(go.Scatterpolar(r=vals + [vals[0]],
                                      theta=categories + [categories[0]],
                                      fill="toself", name=r["method"]))
    fig.update_layout(title="v5 Track 1 — correction method radar (higher=better on all 3 axes)",
                      polar=dict(radialaxis=dict(range=[0, 1])))
    fig.write_html(FIGS_INTERACTIVE / "v5_correction_radar.html",
                   include_plotlyjs="cdn")

    # ---------------- 6-panel PCA grid ----------------
    fig = make_subplots(rows=2, cols=3, subplot_titles=[m[0] for m in METHODS])
    uniq_c = list(meta_df["cohort"].unique())
    palette = ["#e74c3c", "#16a085", "#3498db", "#f39c12", "#9b59b6", "#2c3e50"]
    for i, (name, _, _) in enumerate(METHODS):
        r, c = divmod(i, 3)
        corr = corrected_by_method.get(name)
        if corr is None:
            continue
        top = corr.var(axis=1).sort_values(ascending=False).head(1500).index
        Xp = StandardScaler().fit_transform(corr.loc[top].T.values)
        try:
            coords = PCA(n_components=2).fit_transform(Xp)
        except Exception:
            continue
        for j, cc in enumerate(uniq_c):
            m = meta_df["cohort"].values == cc
            fig.add_trace(go.Scatter(x=coords[m, 0], y=coords[m, 1],
                                     mode="markers", name=cc,
                                     marker=dict(size=4, color=palette[j % len(palette)]),
                                     showlegend=(i == 0),
                                     legendgroup=cc),
                          row=r + 1, col=c + 1)
    fig.update_layout(height=620, title="v5 Track 1 — PCA before and after each correction",
                      showlegend=True)
    fig.write_html(FIGS_INTERACTIVE / "v5_correction_pca_grid.html",
                   include_plotlyjs="cdn")

    # ---------------- LODO heatmap ----------------
    per_method_percohort = {}
    for name, _, _ in METHODS:
        corr = corrected_by_method.get(name)
        if corr is None:
            per_method_percohort[name] = {c: float("nan") for c in uniq_c}
            continue
        top = corr.var(axis=1).sort_values(ascending=False).head(1500).index
        Xp = StandardScaler().fit_transform(corr.loc[top].T.values)
        _, per_c = lodo_auc(Xp, y_vec, cohort_vec)
        per_method_percohort[name] = {r["held"]: r["auc"] for r in per_c}
    hm_df = pd.DataFrame(per_method_percohort).T  # method x cohort
    hm_df.to_csv(RESULTS_ML / "v5_track1_lodo_percohort.tsv", sep="\t")
    fig = go.Figure(data=go.Heatmap(
        z=hm_df.values,
        x=hm_df.columns.astype(str).tolist(),
        y=hm_df.index.tolist(),
        colorscale="Viridis", zmin=0, zmax=1,
        text=np.round(hm_df.values, 2), texttemplate="%{text}"))
    fig.update_layout(title="v5 Track 1 — LODO AUC by method x held-out cohort",
                      xaxis_title="held-out cohort", yaxis_title="correction method")
    fig.write_html(FIGS_INTERACTIVE / "v5_correction_lodo_heatmap.html",
                   include_plotlyjs="cdn")

    # pick best method on combined score
    df["score"] = (
        (1 - df["post_identifiability_auc"].fillna(1.0)) * 0.35
        + df["post_lodo_auc"].fillna(0.0) * 0.35
        + df["bio_preservation_auc"].fillna(0.0) * 0.30
    )
    best = df.sort_values("score", ascending=False).iloc[0]
    safe_write_json(RESULTS_ML / "v5_track1_summary.json", {
        "n_methods_tested": int(len(METHODS)),
        "n_methods_ok": int((df["status"] == "OK").sum()),
        "best_method": str(best["method"]),
        "best_method_verdict": str(best.get("verdict", "")),
        "best_method_metrics": {
            "post_identifiability_auc": float(best["post_identifiability_auc"])
            if not pd.isna(best["post_identifiability_auc"]) else None,
            "post_lodo_auc": float(best["post_lodo_auc"])
            if not pd.isna(best["post_lodo_auc"]) else None,
            "bio_preservation_auc": float(best["bio_preservation_auc"])
            if not pd.isna(best["bio_preservation_auc"]) else None,
        },
        "methods": df.to_dict(orient="records"),
    })

    # also dump top-correction-stable genes for track 4
    if corrected_by_method.get("combat") is not None:
        combat_df = corrected_by_method["combat"]
        sm_idx = load_sample_master().set_index("sample_id")
        # pick tumor vs normal DEG on corrected data — top 200 abs(t-stat)
        tcga = [s for s in meta_df.index if meta_df.loc[s, "cohort"] == "TCGA-THCA"]
        t_y = meta_df.loc[tcga, "label"].values
        X = combat_df.loc[:, tcga].values
        mu_t = X[:, t_y == 1].mean(axis=1)
        mu_n = X[:, t_y == 0].mean(axis=1)
        sd = X.std(axis=1) + 1e-6
        t_stat = (mu_t - mu_n) / sd
        top_stable = pd.DataFrame({
            "gene": combat_df.index,
            "combat_t_stat": t_stat,
            "abs_t": np.abs(t_stat),
        }).sort_values("abs_t", ascending=False).head(500)
        top_stable.to_csv(RESULTS_ML / "v5_track1_top_corrected_genes.tsv", sep="\t", index=False)

    log(f"done; best method = {best['method']} (score={best['score']:.3f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
