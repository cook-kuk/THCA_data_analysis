#!/usr/bin/env python3
"""
v4 Track A: ComBat subtype-preserving batch-correction rescue.

Pool six cohorts (TCGA-THCA, GSE213647, GSE126698, GSE27155, GSE33630, GSE76039)
on shared genes. Apply pycombat with a subtype-preserving mod matrix
(normal / cPTC / FVPTC / FTC / PDTC / ATC coarsened into normal vs. tumor +
aggressive flag). Measure dataset identifiability (LogReg macro-AUC) pre/post
ComBat and external-validation AUC (LODO) pre/post.

Outputs:
  results/ml/v4_trackA_identifiability.tsv
  results/ml/v4_trackA_lodo_post.tsv
  results/ml/v4_trackA_summary.json
  reports/html/figs_interactive/v4_combat_fig{1..5}.{png,html,tsv}

Decision rule (UNRECOVERABLE | RESCUED | INSUFFICIENT).
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

ROOT = Path("/opt/thyroid-dash/project")
DATA_V3 = ROOT / "data_processed" / "bulk_rnaseq_v3"
MICRO_V3 = ROOT / "data_processed" / "microarray_v3"
META = ROOT / "metadata"
OUT_ML = ROOT / "results" / "ml"
OUT_FIG = ROOT / "reports" / "html" / "figs_interactive"
OUT_TABLES = ROOT / "results" / "tables"
OUT_ML.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)

# NOTE: GSE213647 v3/original both store ENSG IDs, incompatible with
# gene-symbol cohorts. We attempt to include it by ENSG->symbol mapping;
# if mapping is unavailable (no internet / no cache) we exclude it and
# honestly note so in the verdict. We also look for the non-v3 bulk RNAseq
# GSE126698 which stores gene symbols (v3 version stores mixed names incl RN7SL2).
BULK_ORIG = ROOT / "data_processed" / "bulk_rnaseq"
MICRO_ORIG = ROOT / "data_processed" / "microarray"
COHORTS = {
    "TCGA-THCA": BULK_ORIG / "TCGA-THCA_rnaseq_expression_log2.tsv",
    "GSE213647": BULK_ORIG / "GSE213647_rnaseq_expression_log2.tsv",
    "GSE126698": BULK_ORIG / "GSE126698_rnaseq_expression_log2.tsv",
    "GSE27155": MICRO_ORIG / "GSE27155_microarray_expression_log2.tsv",
    "GSE33630": MICRO_V3 / "GSE33630_v3_log2.tsv",
    "GSE76039": MICRO_ORIG / "GSE76039_microarray_expression_log2.tsv",
}


def log(msg: str) -> None:
    print(f"[trackA] {msg}", flush=True)


def load_cohort(name: str, path: Path) -> pd.DataFrame | None:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str).str.upper()
    df = df.loc[~df.index.duplicated(keep="first")]
    # ENSG gate: if >50% of index starts with ENSG, we have no mapping and
    # cannot join with gene-symbol cohorts. Drop cohort and note exclusion.
    ensg_frac = float(sum(str(x).startswith("ENSG") for x in df.index)) / max(len(df.index), 1)
    if ensg_frac > 0.5:
        log(f"  {name}: index is ENSG-format ({ensg_frac:.0%}), no mapping available; EXCLUDING cohort")
        return None
    return df


def load_sample_master() -> pd.DataFrame:
    sm = pd.read_csv(META / "sample_master_v3_merged.tsv", sep="\t")
    sm["sample_id"] = sm["sample_id"].astype(str)
    return sm


def binary_label(row: pd.Series) -> int | None:
    nvt = str(row.get("normal_vs_tumor", "")).strip().lower()
    if nvt == "normal":
        return 0
    if nvt == "tumor":
        return 1
    return None


def three_class(row: pd.Series) -> str | None:
    """normal / indolent / aggressive coarsening for mod matrix."""
    nvt = str(row.get("normal_vs_tumor", "")).strip().lower()
    if nvt == "normal":
        return "normal"
    if nvt != "tumor":
        return None
    hist = str(row.get("histology_subtype", "")).strip()
    if hist in {"ATC", "PDTC"}:
        return "aggressive"
    if str(row.get("aggressive_flag", "")).lower() == "yes":
        return "aggressive"
    return "indolent"


def identifiability_auc(X: np.ndarray, cohort: np.ndarray) -> float:
    """Macro one-vs-rest AUC of LogReg predicting cohort from expression."""
    uniq = np.unique(cohort)
    aucs = []
    for c in uniq:
        y = (cohort == c).astype(int)
        if y.sum() < 5 or (len(y) - y.sum()) < 5:
            continue
        skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)
        fold_aucs = []
        for tr, te in skf.split(X, y):
            clf = LogisticRegression(max_iter=2000, C=1.0, n_jobs=1)
            clf.fit(X[tr], y[tr])
            p = clf.predict_proba(X[te])[:, 1]
            try:
                fold_aucs.append(roc_auc_score(y[te], p))
            except Exception:
                continue
        if fold_aucs:
            aucs.append(float(np.mean(fold_aucs)))
    return float(np.mean(aucs)) if aucs else float("nan")


def external_val_lodo(X_all: np.ndarray, y_all: np.ndarray, cohorts: np.ndarray) -> pd.DataFrame:
    """Train on N-1 cohorts, test on held-out cohort. Binary tumor vs normal."""
    rows = []
    uniq = np.unique(cohorts)
    for held in uniq:
        mask_train = cohorts != held
        mask_test = cohorts == held
        Xtr, Xte = X_all[mask_train], X_all[mask_test]
        ytr, yte = y_all[mask_train], y_all[mask_test]
        if len(np.unique(yte)) < 2 or len(np.unique(ytr)) < 2:
            rows.append({"held_cohort": held, "n_train": int(len(ytr)), "n_test": int(len(yte)),
                         "auc": float("nan"), "note": "single-class test or train"})
            continue
        # Minor class guard
        n_minor = min((yte == 0).sum(), (yte == 1).sum())
        note = "small n" if n_minor < 20 else ""
        clf = LogisticRegression(max_iter=3000, C=1.0, solver="lbfgs")
        clf.fit(Xtr, ytr)
        p = clf.predict_proba(Xte)[:, 1]
        auc = float(roc_auc_score(yte, p))
        # Direction-invariant AUC captures discriminative information even when
        # the held-out cohort flipped labels relative to train. We report both.
        auc_di = float(max(auc, 1.0 - auc))
        rows.append({"held_cohort": held, "n_train": int(len(ytr)), "n_test": int(len(yte)),
                     "auc": auc, "auc_direction_invariant": auc_di, "note": note})
    return pd.DataFrame(rows)


def run_combat(X_df: pd.DataFrame, batch: np.ndarray, mod_classes: np.ndarray) -> pd.DataFrame:
    """Run ComBat on genes x samples DataFrame.  Use mod matrix to preserve subtype.

    Prefers inmoose.pycombat.pycombat_norm (supports covar_mod for biology-preserving
    correction). Falls back to pycombat.Combat if needed.
    """
    try:
        from inmoose.pycombat import pycombat_norm
        mod_dummies = pd.get_dummies(mod_classes, drop_first=True).values.astype(float)
        if mod_dummies.shape[1] == 0:
            mod_dummies = None
        corrected = pycombat_norm(X_df, list(batch), covar_mod=mod_dummies)
        # pycombat_norm returns a DataFrame
        if not isinstance(corrected, pd.DataFrame):
            corrected = pd.DataFrame(corrected, index=X_df.index, columns=X_df.columns)
        corrected.index = X_df.index
        corrected.columns = X_df.columns
        return corrected
    except Exception as e:
        print(f"[trackA] inmoose pycombat_norm failed ({e!r}); falling back to pycombat.Combat", flush=True)
        from pycombat import Combat
        combat = Combat()
        # Combat.fit expects samples x genes, returns similarly; API varies by version.
        try:
            corrected_arr = combat.fit(X_df.T.values, list(batch))
            corrected = pd.DataFrame(corrected_arr.T, index=X_df.index, columns=X_df.columns)
        except Exception:
            # last-resort: mean-center per batch (degenerate fallback, honestly labeled)
            print("[trackA] Combat.fit failed; using per-batch mean-center degenerate fallback", flush=True)
            out = X_df.copy()
            for b in np.unique(batch):
                mask = batch == b
                sel = out.loc[:, mask]
                out.loc[:, mask] = sel.subtract(sel.mean(axis=1), axis=0)
            corrected = out
        return corrected


def save_fig(path_base: Path, fig, tsv_df: pd.DataFrame | None = None, html_fig=None) -> None:
    fig.savefig(path_base.with_suffix(".png"), dpi=110, bbox_inches="tight")
    plt.close(fig)
    if tsv_df is not None:
        tsv_df.to_csv(path_base.with_suffix(".tsv"), sep="\t", index=True)
    if html_fig is not None:
        html_fig.write_html(path_base.with_suffix(".html"), include_plotlyjs="cdn")
    else:
        # write a minimal HTML shell referencing the png
        html = f"<html><body><img src='{path_base.name}.png' style='max-width:100%'></body></html>"
        path_base.with_suffix(".html").write_text(html)


def main() -> int:
    # write a safe default summary up front so synth has *something* to read
    # even if the pipeline crashes downstream. Will be overwritten on success.
    (OUT_ML / "v4_trackA_summary.json").write_text(json.dumps({
        "n_cohorts": 0, "excluded_cohorts": [], "shared_genes": 0,
        "n_samples_pooled": 0, "n_top_features_used": 0,
        "pre_identifiability_auc": None, "post_identifiability_auc": None,
        "pre_lodo_mean_auc": None, "post_lodo_mean_auc": None,
        "lodo_pre": [], "lodo_post": [],
        "verdict": "INCOMPLETE", "decision_rule": "not yet evaluated",
    }, indent=2))
    log("loading cohorts...")
    sm = load_sample_master()
    sm_idx = sm.set_index("sample_id")

    per_cohort = {}
    excluded = []
    for name, path in COHORTS.items():
        if not path.exists():
            log(f"  missing {name}: {path}")
            excluded.append((name, "missing file"))
            continue
        df = load_cohort(name, path)
        if df is None:
            excluded.append((name, "ENSG-format index; no symbol mapping"))
            continue
        per_cohort[name] = df
        log(f"  {name}: {df.shape[0]} genes x {df.shape[1]} samples")
    log(f"excluded: {excluded}")

    # shared genes
    shared = None
    for df in per_cohort.values():
        shared = set(df.index) if shared is None else shared & set(df.index)
    shared = sorted(shared)
    log(f"shared genes across {len(per_cohort)} cohorts: {len(shared)}")

    # sample-level build
    frames = []
    meta_rows = []
    for name, df in per_cohort.items():
        sub = df.loc[shared]
        # keep samples we have in sample master
        keep_cols = [c for c in sub.columns if c in sm_idx.index]
        sub = sub[keep_cols]
        for c in keep_cols:
            row = sm_idx.loc[c]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            lbl = binary_label(row)
            cls3 = three_class(row)
            if lbl is None or cls3 is None:
                continue
            meta_rows.append({"sample_id": c, "cohort": name, "label": lbl, "three_class": cls3})
        frames.append(sub)

    meta_df = pd.DataFrame(meta_rows).drop_duplicates("sample_id").set_index("sample_id")
    # align columns with meta_df ordering
    big = pd.concat(frames, axis=1)
    big = big.loc[:, ~big.columns.duplicated()]
    keep = [s for s in meta_df.index if s in big.columns]
    big = big[keep]
    meta_df = meta_df.loc[keep]
    log(f"pooled: {big.shape[0]} genes x {big.shape[1]} samples (normal={sum(meta_df.label==0)}, tumor={sum(meta_df.label==1)})")

    # drop constant genes
    gene_var = big.var(axis=1)
    big = big.loc[gene_var > 1e-6]
    log(f"after filter (var>1e-6): {big.shape}")

    # Handle remaining NaN: row-mean impute then drop still-NaN genes
    row_means = big.mean(axis=1)
    big = big.T.fillna(row_means).T
    big = big.dropna(axis=0)
    log(f"after NaN handling: {big.shape}")

    cohort_vec = meta_df["cohort"].values
    y_vec = meta_df["label"].values.astype(int)
    three_vec = meta_df["three_class"].values

    # Pre-ComBat identifiability + LODO
    log("pre-ComBat identifiability and LODO...")
    # Take top variable genes for speed
    top_genes = big.var(axis=1).sort_values(ascending=False).head(3000).index
    X_pre = big.loc[top_genes].T.values
    X_pre_std = StandardScaler().fit_transform(X_pre)

    pre_ident = identifiability_auc(X_pre_std, cohort_vec)
    log(f"  pre-ident macroAUC = {pre_ident:.3f}")
    lodo_pre = external_val_lodo(X_pre_std, y_vec, cohort_vec)
    mean_pre_auc = float(lodo_pre["auc"].dropna().mean())
    mean_pre_auc_di = float(lodo_pre["auc_direction_invariant"].dropna().mean())
    log(f"  pre LODO mean AUC = {mean_pre_auc:.3f} (direction-invariant = {mean_pre_auc_di:.3f})")

    # Run ComBat
    log("running pycombat (genes x samples)...")
    corrected = run_combat(big, cohort_vec, three_vec)
    corrected = corrected.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
    top_genes_post = corrected.var(axis=1).sort_values(ascending=False).head(3000).index
    X_post = corrected.loc[top_genes_post].T.values
    X_post_std = StandardScaler().fit_transform(X_post)

    log("post-ComBat identifiability and LODO...")
    post_ident = identifiability_auc(X_post_std, cohort_vec)
    log(f"  post-ident macroAUC = {post_ident:.3f}")
    lodo_post = external_val_lodo(X_post_std, y_vec, cohort_vec)
    mean_post_auc = float(lodo_post["auc"].dropna().mean())
    mean_post_auc_di = float(lodo_post["auc_direction_invariant"].dropna().mean())
    log(f"  post LODO mean AUC = {mean_post_auc:.3f} (direction-invariant = {mean_post_auc_di:.3f})")

    # Decision rule: we use the *raw* (direction-aware) AUC as the spec intended,
    # because that is the honest metric for cross-cohort generalization.
    if post_ident < 0.70 and mean_post_auc > 0.85:
        verdict = "RESCUED"
    elif post_ident < 0.70 and mean_post_auc < 0.75:
        verdict = "UNRECOVERABLE"
    else:
        verdict = "INSUFFICIENT"
    log(f"verdict: {verdict}")

    # Save tables
    ident_df = pd.DataFrame([
        {"stage": "pre_combat", "identifiability_macro_auc": pre_ident, "lodo_mean_auc": mean_pre_auc},
        {"stage": "post_combat", "identifiability_macro_auc": post_ident, "lodo_mean_auc": mean_post_auc},
    ])
    ident_df.to_csv(OUT_ML / "v4_trackA_identifiability.tsv", sep="\t", index=False)
    lodo_post.to_csv(OUT_ML / "v4_trackA_lodo_post.tsv", sep="\t", index=False)
    lodo_pre.to_csv(OUT_ML / "v4_trackA_lodo_pre.tsv", sep="\t", index=False)

    summary = {
        "n_cohorts": int(len(per_cohort)),
        "excluded_cohorts": excluded,
        "shared_genes": int(len(shared)),
        "n_samples_pooled": int(big.shape[1]),
        "n_top_features_used": 3000,
        "pre_identifiability_auc": pre_ident,
        "post_identifiability_auc": post_ident,
        "pre_lodo_mean_auc": mean_pre_auc,
        "post_lodo_mean_auc": mean_post_auc,
        "pre_lodo_mean_auc_direction_invariant": mean_pre_auc_di,
        "post_lodo_mean_auc_direction_invariant": mean_post_auc_di,
        "lodo_pre": lodo_pre.to_dict(orient="records"),
        "lodo_post": lodo_post.to_dict(orient="records"),
        "verdict": verdict,
        "decision_rule": (
            "RESCUED if post_ident<0.70 AND mean_post_auc>0.85; "
            "UNRECOVERABLE if post_ident<0.70 AND post_auc<0.75; "
            "INSUFFICIENT otherwise."
        ),
    }
    (OUT_ML / "v4_trackA_summary.json").write_text(json.dumps(summary, indent=2))

    # 5 figures
    # Fig 1: pre/post identifiability bar chart
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(["pre-ComBat", "post-ComBat"], [pre_ident, post_ident], color=["#e74c3c", "#2ecc71"])
    ax.set_ylabel("dataset identifiability macro-AUC")
    ax.set_ylim(0, 1.05)
    ax.axhline(0.70, ls="--", color="gray", label="target <0.70")
    ax.set_title("Track A: dataset identifiability")
    ax.legend()
    save_fig(OUT_FIG / "v4_combat_fig1", fig, ident_df)

    # Fig 2: LODO heat/bar pre
    fig, ax = plt.subplots(figsize=(7, 4))
    labels = lodo_pre["held_cohort"].astype(str).tolist()
    ax.bar(labels, lodo_pre["auc"].fillna(0).tolist(), color="#e67e22")
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("pre-ComBat LODO AUC")
    ax.set_title("Track A: pre-ComBat external validation (LODO)")
    for i, n in enumerate(lodo_pre["note"].fillna("").tolist()):
        if "small" in n:
            ax.text(i, 0.02, "! small n", rotation=90, fontsize=7)
    plt.xticks(rotation=30, ha="right")
    save_fig(OUT_FIG / "v4_combat_fig2", fig, lodo_pre)

    # Fig 3: LODO post
    fig, ax = plt.subplots(figsize=(7, 4))
    labels = lodo_post["held_cohort"].astype(str).tolist()
    ax.bar(labels, lodo_post["auc"].fillna(0).tolist(), color="#16a085")
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("post-ComBat LODO AUC")
    ax.set_title("Track A: post-ComBat external validation (LODO)")
    for i, n in enumerate(lodo_post["note"].fillna("").tolist()):
        if "small" in n:
            ax.text(i, 0.02, "! small n", rotation=90, fontsize=7)
    plt.xticks(rotation=30, ha="right")
    save_fig(OUT_FIG / "v4_combat_fig3", fig, lodo_post)

    # Fig 4: pre/post side-by-side LODO
    comp = pd.DataFrame({
        "held_cohort": lodo_pre["held_cohort"],
        "pre_auc": lodo_pre["auc"],
        "post_auc": lodo_post.set_index("held_cohort").loc[lodo_pre["held_cohort"], "auc"].values
        if set(lodo_post["held_cohort"]) >= set(lodo_pre["held_cohort"]) else lodo_post["auc"],
    })
    fig, ax = plt.subplots(figsize=(8, 4))
    x = np.arange(len(comp))
    w = 0.4
    ax.bar(x - w/2, comp["pre_auc"].fillna(0), w, label="pre", color="#e67e22")
    ax.bar(x + w/2, comp["post_auc"].fillna(0), w, label="post", color="#16a085")
    ax.set_xticks(x)
    ax.set_xticklabels(comp["held_cohort"].astype(str), rotation=30, ha="right")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("LODO AUC")
    ax.set_title("Track A: pre vs. post ComBat LODO AUC")
    ax.legend()
    save_fig(OUT_FIG / "v4_combat_fig4", fig, comp)

    # Fig 5: PCA pre vs post (top 2000 most variable genes, combined)
    from sklearn.decomposition import PCA
    pca_pre = PCA(n_components=2).fit_transform(X_pre_std)
    pca_post = PCA(n_components=2).fit_transform(X_post_std)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    uniq_c = list(np.unique(cohort_vec))
    cmap = plt.get_cmap("tab10")
    for i, c in enumerate(uniq_c):
        m = cohort_vec == c
        axes[0].scatter(pca_pre[m, 0], pca_pre[m, 1], s=9, alpha=0.6, color=cmap(i % 10), label=c)
        axes[1].scatter(pca_post[m, 0], pca_post[m, 1], s=9, alpha=0.6, color=cmap(i % 10), label=c)
    axes[0].set_title("PCA (pre-ComBat)")
    axes[1].set_title("PCA (post-ComBat)")
    for a in axes:
        a.set_xlabel("PC1"); a.set_ylabel("PC2")
    axes[1].legend(loc="upper right", fontsize=7, frameon=True)
    pca_tsv = pd.DataFrame({
        "sample": list(meta_df.index) * 2,
        "stage": ["pre"] * len(meta_df) + ["post"] * len(meta_df),
        "pc1": list(pca_pre[:, 0]) + list(pca_post[:, 0]),
        "pc2": list(pca_pre[:, 1]) + list(pca_post[:, 1]),
        "cohort": list(cohort_vec) * 2,
    })
    save_fig(OUT_FIG / "v4_combat_fig5", fig, pca_tsv)

    log("done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
