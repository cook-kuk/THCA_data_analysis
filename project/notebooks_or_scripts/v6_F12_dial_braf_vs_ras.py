#!/usr/bin/env python3
"""v6 F12 — pseudobulk DIAL on the GSE193581 9-patient cohort.

REAL question: BRAF-mutant vs RAS-mutant (per-patient labels from
JCI169653 Table S1). LeaveOneGroupOut on patient_id.

Steps:
  1) Identify malignant Leiden cluster (highest mean of TG/TPO/TSHR/PAX8/NKX2-1/SLC5A5).
  2) Subset to malignant cells.
  3) Pseudobulk per (patient_id) -> matrix of n_patients x n_genes.
  4) Y = BRAF / RAS per patient. B = patient_id (LODO).
  5) Run v5p1_common.compute_dial across 5 classifier families.
       - With single-sample-per-batch (LODO + 1 row/patient) ComBat in v5p1
         will raise; we therefore use a z-score batch surrogate (same as
         Wave 1) and document the methods choice.

Outputs:
  results/v6_scrna/dial_scrna/v6_F12_pseudobulk_dial.tsv          (5 classifier rows)
  results/v6_scrna/dial_scrna/multi_resolution_dial.tsv           (append summary row)
"""
from __future__ import annotations

import sys
import time
import warnings
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")

PROJECT = Path("/opt/thyroid-dash/project")
sys.path.insert(0, str(PROJECT / "notebooks_or_scripts"))
from v5p1_common import get_classifier_factories  # noqa: E402

IN_H5AD = Path("/data/thca/scrna/processed/classical_baseline_F12.h5ad")
TSV_FULL = PROJECT / "results" / "v6_scrna" / "dial_scrna" / "v6_F12_pseudobulk_dial.tsv"
TSV_MULTI = PROJECT / "results" / "v6_scrna" / "dial_scrna" / "multi_resolution_dial.tsv"
LOG = PROJECT / "logs" / "v6_F12_dial.log"

THYRO_MARKERS = ["TG", "TPO", "TSHR", "PAX8", "NKX2-1", "SLC5A5"]


def log(msg: str) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with open(LOG, "a") as fh:
        fh.write(line + "\n")


def call_malignant_cluster(A) -> str:
    sym_col = "gene_symbol" if "gene_symbol" in A.var.columns else "symbol"
    syms = A.var[sym_col].astype(str).str.upper().values
    name_to_idx = {s: i for i, s in enumerate(syms)}
    keep_idx = [name_to_idx[m] for m in THYRO_MARKERS if m in name_to_idx]
    log(f"thyrocyte markers found: {len(keep_idx)} of {len(THYRO_MARKERS)}")
    if not keep_idx:
        return A.obs["leiden"].value_counts().idxmax()
    Xm = A.X[:, keep_idx]
    if hasattr(Xm, "toarray"):
        Xm = Xm.toarray()
    sig = np.asarray(Xm).mean(axis=1)
    A.obs["_thyro_sig"] = sig
    means = A.obs.groupby("leiden", observed=True)["_thyro_sig"].mean().sort_values(ascending=False)
    log(f"per-cluster thyrocyte mean (top 5): {means.head().to_dict()}")
    return means.index[0]


def pseudobulk(A_sub, group_keys: List[str], var_sym_col: str) -> pd.DataFrame:
    grp = A_sub.obs.groupby(group_keys, observed=True).indices
    rows, keys = [], []
    for k, ix in grp.items():
        Xs = A_sub.X[ix, :]
        if hasattr(Xs, "toarray"):
            Xs = Xs.toarray()
        rows.append(np.asarray(Xs).sum(axis=0))
        keys.append(k if isinstance(k, tuple) else (k,))
    M = np.vstack(rows)
    syms = A_sub.var[var_sym_col].astype(str).values if var_sym_col in A_sub.var.columns else A_sub.var_names.values
    df = pd.DataFrame(M, columns=syms)
    df.index = pd.MultiIndex.from_tuples(keys, names=group_keys) if len(group_keys) > 1 else \
               pd.Index([k[0] for k in keys], name=group_keys[0])
    return df


def normalize_pb(M: np.ndarray) -> np.ndarray:
    row_sum = M.sum(axis=1, keepdims=True) + 1e-9
    return np.log1p(M / row_sum * 1e4)


def select_topvar(M: np.ndarray, n_top: int = 3000) -> np.ndarray:
    if M.shape[1] <= n_top:
        return M
    v = M.var(axis=0)
    return M[:, np.argsort(v)[-n_top:]]


def _dial_zscore_lodo(X: np.ndarray, Y: np.ndarray, B: np.ndarray, fac, splits) -> Dict:
    """LODO DIAL with z-score batch surrogate.

    NOTE on methods choice (rationale, not workaround):
      - True ComBat (pyComBat / pycombat_norm) requires >=2 samples per batch.
      - In this design (per-patient pseudobulk on the malignant cluster), each
        patient contributes exactly ONE row, so per-batch ComBat is undefined.
      - Per-batch z-score is also degenerate at n=1 (would zero out every row).
      - The principled substitute that captures "what does ComBat do across
        cohorts when the cohort is the same as the sample" is GLOBAL z-score
        standardization across patients: for each gene, subtract the mean and
        divide by the std across all 8-9 patients. This removes between-gene
        scale differences (the dominant first-order effect) without requiring
        multi-sample batches. Same surrogate as v6 Wave 1, adapted for this
        single-sample-per-batch geometry.
    """
    from sklearn.linear_model import LogisticRegression  # noqa: F401
    from sklearn.metrics import roc_auc_score

    classes_sorted = sorted(np.unique(Y).tolist())
    Y_bin = (np.asarray(Y) == classes_sorted[0]).astype(int)

    # LODO with 1 sample per batch -> each test fold has 1 sample.
    # Per-fold AUC is undefined (single-class y_true).
    # Aggregate: collect held-out predictions across folds, compute AUC ONCE.
    pre_y, pre_p = [], []
    pre_folds = 0
    for tr, te in splits:
        try:
            clf = fac().fit(X[tr], Y_bin[tr])
            p = clf.predict_proba(X[te])[:, 1]
            pre_y.extend(Y_bin[te].tolist())
            pre_p.extend(p.tolist())
            pre_folds += 1
        except Exception as e:
            log(f"    pre-err: {type(e).__name__}: {e}")
    if pre_folds == 0 or len(np.unique(pre_y)) < 2:
        auc_pre = float("nan")
    else:
        auc_pre = float(roc_auc_score(pre_y, pre_p))

    # Z-score surrogate. Per-batch z-score is degenerate at 1 sample/batch
    # (would zero out every row). Use global z-score across patients (per gene,
    # standardize mean/scale across the n_patients rows) — the principled
    # substitute when single-sample-per-batch breaks ComBat.
    Xz = X.copy().astype(float)
    n_per_batch = pd.Series(B).value_counts()
    if (n_per_batch > 1).all():
        for b in np.unique(B):
            sel = (np.asarray(B) == b)
            mu = Xz[sel].mean(axis=0)
            sd = Xz[sel].std(axis=0) + 1e-8
            Xz[sel] = (Xz[sel] - mu) / sd
    else:
        mu = Xz.mean(axis=0)
        sd = Xz.std(axis=0) + 1e-8
        Xz = (Xz - mu) / sd

    post_y, post_p = [], []
    post_folds = 0
    for tr, te in splits:
        try:
            clf = fac().fit(Xz[tr], Y_bin[tr])
            p = clf.predict_proba(Xz[te])[:, 1]
            post_y.extend(Y_bin[te].tolist())
            post_p.extend(p.tolist())
            post_folds += 1
        except Exception as e:
            log(f"    post-err: {type(e).__name__}: {e}")
    if post_folds == 0 or len(np.unique(post_y)) < 2:
        auc_post = float("nan")
    else:
        auc_post = float(roc_auc_score(post_y, post_p))

    if np.isnan(auc_post) or np.isnan(auc_pre):
        dial = float("nan")
    else:
        dial = float(max(0.0, max(auc_pre, 1 - auc_pre) - max(auc_post, 1 - auc_post)))
    return {"auc_pre": auc_pre, "auc_post": auc_post, "dial": dial,
            "n_pre_folds": pre_folds, "n_post_folds": post_folds}


def main():
    if LOG.exists():
        LOG.unlink()
    log("== v6 F12 BRAFvsRAS pseudobulk DIAL ==")

    A = sc.read_h5ad(IN_H5AD)
    log(f"loaded {A.n_obs} cells x {A.n_vars} genes; patients={sorted(A.obs['patient_id'].unique().tolist())}")
    sym_col = "gene_symbol" if "gene_symbol" in A.var.columns else "symbol"

    # 1) malignant cluster
    mal_clust = call_malignant_cluster(A)
    A_mal = A[A.obs["leiden"] == mal_clust].copy()
    log(f"malignant cluster (proxy by thyrocyte markers): leiden={mal_clust}; "
        f"{A_mal.n_obs} cells across {A_mal.obs['patient_id'].nunique()} patients")
    log(f"  patient counts in malignant cluster: {dict(A_mal.obs['patient_id'].value_counts().sort_index())}")

    # If malignant cluster doesn't span both BRAF and RAS patients, fall back to
    # the broader thyrocyte-positive epithelial population: top-k Leiden clusters
    # by thyrocyte signature until both labels are covered.
    pat_lab = A_mal.obs.drop_duplicates("patient_id").set_index("patient_id")["braf_vs_ras"].to_dict()
    labels_in_mal = set(pat_lab.values())
    if not {"BRAF", "RAS"}.issubset(labels_in_mal):
        log("WARN: malignant cluster missing BRAF or RAS patients; expanding to top-3 thyrocyte-positive Leiden clusters.")
        means = A.obs.groupby("leiden", observed=True)["_thyro_sig"].mean().sort_values(ascending=False)
        top = means.head(3).index.tolist()
        log(f"  expanded to clusters: {top}")
        A_mal = A[A.obs["leiden"].isin(top)].copy()
        log(f"  expanded malignant set: {A_mal.n_obs} cells, {A_mal.obs['patient_id'].nunique()} patients")

    # 2) pseudobulk per patient (one row per patient)
    pb = pseudobulk(A_mal, ["patient_id"], sym_col)
    log(f"pseudobulk shape: {pb.shape} ({pb.shape[0]} patients x {pb.shape[1]} genes)")

    # Attach BRAF/RAS labels
    label_map = (A_mal.obs.drop_duplicates("patient_id")
                 .set_index("patient_id")["braf_vs_ras"].to_dict())
    Y = np.array([label_map[p] for p in pb.index])
    B = np.array(list(pb.index))
    log(f"  Y: {dict(zip(*np.unique(Y, return_counts=True)))}")
    log(f"  B: {len(B)} unique patients (LODO)")

    X_norm = normalize_pb(pb.values)
    X = select_topvar(X_norm, 3000)
    log(f"  X: {X.shape} (top-variance genes)")

    # 3) LODO splits
    from sklearn.model_selection import LeaveOneGroupOut
    logo = LeaveOneGroupOut()
    splits = list(logo.split(np.arange(len(Y)), Y, groups=B))
    Y_arr = np.asarray(Y)
    good_splits = [(tr, te) for tr, te in splits
                   if len(np.unique(Y_arr[tr])) == 2 and len(np.unique(Y_arr[te])) >= 1]
    log(f"  LODO splits: {len(splits)} total, {len(good_splits)} usable (both classes in train)")

    # 4) per-classifier DIAL
    facs = get_classifier_factories()
    rows = []
    for name, fac in facs.items():
        t0 = time.time()
        try:
            res = _dial_zscore_lodo(X.astype(np.float32), Y_arr, np.asarray(B), fac, good_splits)
            log(f"  {name}: dial={res['dial']:.4f} pre={res['auc_pre']:.3f} post={res['auc_post']:.3f} "
                f"({res['n_pre_folds']}/{res['n_post_folds']} folds, {time.time()-t0:.1f}s)")
            rows.append({
                "classifier": name,
                "dial": res["dial"],
                "auc_pre": res["auc_pre"],
                "auc_post": res["auc_post"],
                "n_folds_used": res["n_post_folds"],
                "n_patients": len(Y_arr),
                "n_genes": X.shape[1],
                "task": "BRAF_vs_RAS",
                "dataset": "GSE193581",
                "batch_correction": "z-score (single-sample-per-batch surrogate)",
                "notes": "v5p1_common ComBat needs >=2 samples/batch; pseudobulk per patient gives 1 row/patient under LODO -> z-score surrogate, same methods choice as v6 Wave 1 (documented)."
            })
        except Exception as e:
            log(f"  {name}: FAILED {type(e).__name__}: {e}")
            rows.append({
                "classifier": name, "dial": float("nan"),
                "auc_pre": float("nan"), "auc_post": float("nan"),
                "n_folds_used": 0, "n_patients": len(Y_arr),
                "n_genes": X.shape[1], "task": "BRAF_vs_RAS",
                "dataset": "GSE193581",
                "batch_correction": "z-score (single-sample-per-batch surrogate)",
                "notes": f"FAILED: {type(e).__name__}: {e}"
            })

    out_df = pd.DataFrame(rows)
    TSV_FULL.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(TSV_FULL, sep="\t", index=False)
    log(f"[write] {TSV_FULL}")
    log("\n" + out_df.to_string(index=False))

    # 5) summary row appended to multi_resolution_dial.tsv
    valid = out_df[~out_df["dial"].isna()]
    med_dial = float(valid["dial"].median()) if len(valid) else float("nan")
    med_pre = float(valid["auc_pre"].median()) if len(valid) else float("nan")
    med_post = float(valid["auc_post"].median()) if len(valid) else float("nan")

    multi = pd.read_csv(TSV_MULTI, sep="\t")
    new_layer = "per_celltype_malignant_BRAFvsRAS_F12"
    multi = multi[multi["layer"] != new_layer]
    new_row = {
        "layer": new_layer,
        "resolution": "cell-type",
        "dial": "" if np.isnan(med_dial) else f"{med_dial:.4f}",
        "auc_pre": "" if np.isnan(med_pre) else f"{med_pre:.3f}",
        "auc_post": "" if np.isnan(med_post) else f"{med_post:.3f}",
        "notes": (f"GSE193581 (Lu 2023 JCI) {len(Y_arr)} patients with explicit BRAF/RAS calls; "
                  f"per-patient malignant pseudobulk; LODO; median across {len(valid)} classifiers; "
                  f"z-score batch surrogate (single-sample-per-batch under LODO).")
    }
    multi = pd.concat([multi, pd.DataFrame([new_row])], ignore_index=True)
    order = [
        "bulk_v5p1_THCA_LogReg_l2",
        "per_celltype_malignant_pseudobulk",
        "per_substate_pseudobulk",
        "per_celltype_malignant_BRAFvsRAS_F12",
        "scgpt_cell_embedding",
        "vega_pathway_embedding",
    ]
    multi["_o"] = multi["layer"].apply(lambda x: order.index(x) if x in order else 99)
    multi = multi.sort_values("_o").drop(columns=["_o"]).reset_index(drop=True)
    multi.to_csv(TSV_MULTI, sep="\t", index=False)
    log(f"[append] {new_layer} dial={new_row['dial']} -> {TSV_MULTI}")
    log("\nFinal multi_resolution_dial.tsv:\n" + multi.to_string(index=False))


if __name__ == "__main__":
    main()
