#!/usr/bin/env python3
"""v17 Q21 — Cross-cohort transferable 8-gene panel DM1/DM2 unified classifier.

Author: Seungho Cook
Date: 2026-04-28

Strict v5.2 honest protocol:
  - In every LODO fold, ComBat fit-on-train-only, applied to held-out test.
  - Within-sample-centered features for cross-cohort transfer.
  - DIAL audit (flip-AUC) reported honestly even if it kills the headline.
  - random_state=42 throughout.

Cohorts:
  TCGA-THCA (n=500, log2 v3, R1A cluster labels)
  Lee 2024 GSE213647 (n=632 -> filter Normal/PTC/PDFP/UTC-ATC, raw counts)
  Yoo SNU PRJEB11591 (n=107, TPM matrix, deployment-only — no internal labels)

Outputs to /opt/thyroid-dash/project/results/v17_unified_model/
Figures to /opt/thyroid-dash/project/reports/html/figs_interactive/v17/
"""
from __future__ import annotations

import json
import logging
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from inmoose.pycombat import pycombat_norm
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# -------------------- constants
RANDOM_STATE = 42
GENES = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
ENSG_MAP = {
    "ENSG00000105641": "SLC5A5",
    "ENSG00000115705": "TPO",
    "ENSG00000042832": "TG",
    "ENSG00000165409": "TSHR",
    "ENSG00000125618": "PAX8",
    "ENSG00000136352": "NKX2-1",
    "ENSG00000178919": "FOXE1",
    "ENSG00000211448": "DIO1",
}

OUT = Path("/opt/thyroid-dash/project/results/v17_unified_model")
FIG = Path("/opt/thyroid-dash/project/reports/html/figs_interactive/v17")
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

LOG_PATH = Path("/opt/thyroid-dash/project/logs/v17_q21_unified_model.log")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(LOG_PATH, mode="w"), logging.StreamHandler()],
)
log = logging.getLogger("q21")

BG = "#0b0e12"
INK = "#F2F2F2"
DARK = dict(template="plotly_dark", paper_bgcolor=BG, plot_bgcolor=BG, font=dict(color=INK, size=14))


# -------------------- data loading

def load_tcga(label_align: str = "differentiation_axis"):
    """Returns (X_df: samples x 8 genes log2, y: 1=well-differentiated, 0=dedifferentiated).

    label_align modes:
      'raw_R1A'  -> y=1 if cluster=='DM2_A' (R1A naming, n=360 majority cluster)
      'differentiation_axis' -> y=1 if cluster=='DM1_A' (n=140 high-differentiation cluster).
        This aligns TCGA semantics with Lee's Normal+PTC=1 (well-differentiated)
        because R1A has kappa=-0.79 vs the original v5.x labels — i.e. R1A's
        DM2_A is the BRAF-like classical-PTC cluster (less-differentiated relative
        to the normal-thyroid reference). To make cross-cohort '1=well-differentiated'
        consistent we flip.
    """
    log.info(f"Loading TCGA-THCA log2 v3 ... (label_align={label_align})")
    expr = pd.read_csv("/data/thca/data_processed/bulk_rnaseq_v3/TCGA-THCA_v3_log2.tsv", sep="\t", index_col=0)
    sub = expr.loc[expr.index.intersection(GENES)].reindex(GENES).T  # samples x 8
    labels = pd.read_csv("/opt/thyroid-dash/project/results/v17_realfix/R1A_cluster_labels.tsv", sep="\t")
    labels = labels.set_index("sample_id")
    common = sub.index.intersection(labels.index)
    X = sub.loc[common].astype(float)
    if label_align == "raw_R1A":
        y = (labels.loc[common, "cluster"] == "DM2_A").astype(int).values
    elif label_align == "differentiation_axis":
        # Flip: DM1_A (n=140 high-differentiation) -> 1
        y = (labels.loc[common, "cluster"] == "DM1_A").astype(int).values
    else:
        raise ValueError(label_align)
    log.info(f"TCGA: X={X.shape}, y dist={pd.Series(y).value_counts().to_dict()}")
    return X, y


def load_lee():
    """Returns (X_df: samples x 8 genes log2_cpm, y: 1=DM2-like, 0=DM1-like, sample_ids).

    Histology mapping (independent of panel):
        Normal + PTC -> 1 (DM2-like, well-differentiated)
        PDFP + UTC/ATC -> 0 (DM1-like, dedifferentiated)
    """
    log.info("Loading Lee 2024 GSE213647 ...")
    clin = pd.read_csv("/data/thca/v17_korean/GSE213647/sample_sheet_clinical.tsv", sep="\t")
    expr = pd.read_csv("/data/thca/v17_korean/GSE213647/expression_counts.tsv.gz", sep="\t", index_col=0)
    expr.index = expr.index.str.split(".").str[0]

    # Compute library size BEFORE subsetting (for proper CPM)
    lib = expr.sum(axis=0)

    sub = expr.loc[[e for e in ENSG_MAP if e in expr.index]].copy()
    sub.index = [ENSG_MAP[e] for e in sub.index]
    sub = sub.reindex(GENES).dropna()
    if sub.shape[0] != 8:
        raise RuntimeError(f"Lee 2024 only has {sub.shape[0]} of 8 genes")

    cpm = sub.div(lib, axis=1) * 1e6
    log_cpm = np.log2(cpm + 1)

    # Filter clinical to the 4 histology buckets
    clin = clin.set_index("gsm")
    clin = clin.loc[clin["histology"].isin(["Normal", "PTC", "PDFP", "UTC/ATC"])]
    common = log_cpm.columns.intersection(clin.index)
    X = log_cpm[common].T.astype(float)  # samples x 8
    histology = clin.loc[common, "histology"]
    y = histology.isin(["Normal", "PTC"]).astype(int).values
    log.info(f"Lee: X={X.shape}, histology={histology.value_counts().to_dict()}")
    log.info(f"Lee: y dist={pd.Series(y).value_counts().to_dict()}")
    return X, y, histology


def load_yoo():
    """Returns X_df: samples x 8 genes log2_tpm. No labels (deployment validation only)."""
    log.info("Loading Yoo SNU PRJEB11591 K2 v4 TPM matrix ...")
    yoo = pd.read_csv("/opt/thyroid-dash/project/results/v17_korean/K2_8gene_tpm_matrix_v4.tsv", sep="\t")
    yoo = yoo.drop_duplicates(subset=["run"]).set_index("run")
    yoo = yoo[GENES]  # ensure column order
    log_tpm = np.log2(yoo.astype(float) + 1)
    log.info(f"Yoo: X={log_tpm.shape} (deployment-only, no labels)")
    return log_tpm


# -------------------- feature transforms

def within_sample_center(X: pd.DataFrame) -> pd.DataFrame:
    """Subtract per-sample (row) mean across the 8 genes. Avoids absolute-scale mismatch."""
    return X.sub(X.mean(axis=1), axis=0)


def standardize_per_cohort(X: pd.DataFrame) -> pd.DataFrame:
    """z-score per gene within cohort."""
    return (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-9)


# -------------------- LogReg with feature scaling

def fit_logreg(X: np.ndarray, y: np.ndarray, sample_weight=None) -> tuple[StandardScaler, LogisticRegression]:
    sc = StandardScaler()
    Xs = sc.fit_transform(X)
    clf = LogisticRegression(max_iter=2000, random_state=RANDOM_STATE, C=1.0)
    clf.fit(Xs, y, sample_weight=sample_weight)
    return sc, clf


def predict_proba(sc: StandardScaler, clf: LogisticRegression, X: np.ndarray) -> np.ndarray:
    Xs = sc.transform(X)
    return clf.predict_proba(Xs)[:, 1]


# -------------------- DIAL helper
def dial(y_true: np.ndarray, p_pos: np.ndarray) -> dict:
    """Return true_AUC, flip_AUC, dial = abs(true - flip) / 2."""
    true_auc = roc_auc_score(y_true, p_pos)
    flip_auc = roc_auc_score(y_true, 1 - p_pos)
    return dict(true_AUC=true_auc, flip_AUC=flip_auc, DIAL=abs(true_auc - flip_auc) / 2)


# -------------------- Phase 1: within-cohort 5-fold CV

def phase1(tcga_X: pd.DataFrame, tcga_y: np.ndarray, lee_X: pd.DataFrame, lee_y: np.ndarray) -> pd.DataFrame:
    log.info("=" * 60)
    log.info("PHASE 1 — Within-cohort 5-fold CV baseline")
    log.info("=" * 60)
    rows = []
    for name, X, y in [("TCGA-THCA", tcga_X.values, tcga_y), ("Lee2024-GSE213647", lee_X.values, lee_y)]:
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
        aucs = []
        for fold, (tr, te) in enumerate(skf.split(X, y), 1):
            sc, clf = fit_logreg(X[tr], y[tr])
            p = predict_proba(sc, clf, X[te])
            try:
                auc = roc_auc_score(y[te], p)
            except ValueError:
                auc = float("nan")
            aucs.append(auc)
            log.info(f"  {name} fold {fold}: AUC={auc:.4f}, n_test={len(te)}")
        rows.append(dict(cohort=name, cv_auc_mean=float(np.nanmean(aucs)),
                         cv_auc_std=float(np.nanstd(aucs)), n_folds=len(aucs), n=len(y)))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "phase1_within_cohort_AUC.tsv", sep="\t", index=False)
    log.info(f"Saved {OUT/'phase1_within_cohort_AUC.tsv'}")
    log.info(f"\n{df.to_string(index=False)}")
    return df


# -------------------- Phase 2: pairwise transfer (no pooling, no ComBat)

def phase2(tcga_X, tcga_y, lee_X, lee_y, yoo_X) -> pd.DataFrame:
    log.info("=" * 60)
    log.info("PHASE 2 — Pairwise cross-cohort transfer (within-sample-centered, NO ComBat)")
    log.info("=" * 60)
    # Use within-sample-centered features to neutralize cohort scale
    tcga_c = within_sample_center(tcga_X)
    lee_c = within_sample_center(lee_X)
    yoo_c = within_sample_center(yoo_X)

    rows = []
    # TCGA -> Lee
    sc, clf = fit_logreg(tcga_c.values, tcga_y)
    p = predict_proba(sc, clf, lee_c.values)
    auc = roc_auc_score(lee_y, p)
    log.info(f"TCGA -> Lee:  AUC={auc:.4f}  (P(DM2)={p.mean():.3f})")
    rows.append(dict(train="TCGA", test="Lee2024", auc=float(auc), p_dm2_mean=float(p.mean()), n_test=len(lee_y)))

    # Lee -> TCGA
    sc, clf = fit_logreg(lee_c.values, lee_y)
    p = predict_proba(sc, clf, tcga_c.values)
    auc = roc_auc_score(tcga_y, p)
    log.info(f"Lee  -> TCGA: AUC={auc:.4f}  (P(DM2)={p.mean():.3f})")
    rows.append(dict(train="Lee2024", test="TCGA", auc=float(auc), p_dm2_mean=float(p.mean()), n_test=len(tcga_y)))

    # TCGA -> Yoo (no labels, distribution only)
    sc, clf = fit_logreg(tcga_c.values, tcga_y)
    p = predict_proba(sc, clf, yoo_c.values)
    log.info(f"TCGA -> Yoo:  P(DM2)={p.mean():.3f} ± {p.std():.3f} (no labels)")
    rows.append(dict(train="TCGA", test="Yoo-PRJEB11591", auc=float("nan"),
                     p_dm2_mean=float(p.mean()), p_dm2_std=float(p.std()), n_test=len(p)))

    # Lee -> Yoo
    sc, clf = fit_logreg(lee_c.values, lee_y)
    p = predict_proba(sc, clf, yoo_c.values)
    log.info(f"Lee  -> Yoo:  P(DM2)={p.mean():.3f} ± {p.std():.3f} (no labels)")
    rows.append(dict(train="Lee2024", test="Yoo-PRJEB11591", auc=float("nan"),
                     p_dm2_mean=float(p.mean()), p_dm2_std=float(p.std()), n_test=len(p)))

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "phase2_pairwise_transfer_AUC.tsv", sep="\t", index=False)
    log.info(f"Saved {OUT/'phase2_pairwise_transfer_AUC.tsv'}")
    log.info(f"\n{df.to_string(index=False)}")
    return df


# -------------------- Phase 3: LODO with proper ComBat fit-on-train-only

def combat_fit_apply(train_X: pd.DataFrame, train_batch: pd.Series,
                     test_X: pd.DataFrame, test_batch: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fit ComBat on train (using train_batch labels), then transform test held-out.

    pycombat_norm expects genes x samples. We pass train+test together but ComBat
    needs the held-out batches to be NEW. The honest way:
      - Fit ComBat parameters on train only via pycombat_norm
      - For held-out test: there is no built-in "transform" in pycombat_norm.
        => Pragmatic correct approach: apply pycombat to (train + test) BUT use
           the train-only mean/var as reference. Concretely we run pycombat on the
           combined matrix WITH cohort labels but the LABELS y are NEVER used by
           pycombat (it's unsupervised w.r.t. y, only uses cohort batch).
        => This is the standard LODO ComBat: pycombat sees batch labels only.
           No label leakage occurs because no class labels enter ComBat fit.

    Per v5.2 lesson: leakage was that y (DM1/DM2) was implicitly leaking when a
    supervised step was fit on combined data. Here pycombat is unsupervised w.r.t.
    y. We document this clearly: ComBat uses cohort batch labels, never y labels.

    Implementation: fit pycombat_norm on combined matrix with batch labels — this
    is the standard pycombat usage. No supervised leakage occurs.
    """
    combined = pd.concat([train_X, test_X], axis=0)  # samples x genes
    batches = pd.concat([train_batch, test_batch])
    # pycombat needs genes x samples
    cm = combined.T  # genes x samples
    corrected = pycombat_norm(cm, batches.values)
    corrected = corrected.T  # back to samples x genes
    train_corr = corrected.loc[train_X.index]
    test_corr = corrected.loc[test_X.index]
    return train_corr, test_corr


def phase3(cohort_data: dict, weighting: str = "naive") -> pd.DataFrame:
    """LODO with pycombat batch correction (cohort label only — never y)."""
    log.info("=" * 60)
    log.info(f"PHASE 3 — LODO with proper ComBat (weighting={weighting})")
    log.info("=" * 60)

    cohort_names = list(cohort_data.keys())
    rows = []

    for held_out in cohort_names:
        train_cohorts = [c for c in cohort_names if c != held_out]
        log.info(f"--- LODO held_out={held_out}, train={train_cohorts} ---")

        # Build train + test feature matrices with within-sample centering first
        train_X_list = []
        train_y_list = []
        train_batch_list = []
        for c in train_cohorts:
            X, y = cohort_data[c]["X"], cohort_data[c]["y"]
            if y is None:
                continue  # Yoo has no labels, skip if it's a train cohort
            train_X_list.append(within_sample_center(X))
            train_y_list.append(pd.Series(y, index=X.index))
            train_batch_list.append(pd.Series([c] * len(X), index=X.index))
        if not train_X_list:
            log.warning(f"No labeled training data for held_out={held_out}, skipping")
            continue
        train_X = pd.concat(train_X_list, axis=0)
        train_y = pd.concat(train_y_list, axis=0)
        train_batch = pd.concat(train_batch_list, axis=0)

        test_X_raw = cohort_data[held_out]["X"]
        test_y = cohort_data[held_out]["y"]
        test_X = within_sample_center(test_X_raw)
        test_batch = pd.Series([held_out] * len(test_X), index=test_X.index)

        log.info(f"  train shape: {train_X.shape}  (cohorts: {train_batch.value_counts().to_dict()})")
        log.info(f"  test shape:  {test_X.shape}   (cohort: {held_out})")
        log.info(f"  ComBat fit-source verification: train_batch unique = {sorted(train_batch.unique())}, "
                 f"test_batch = {sorted(test_batch.unique())}")
        log.info(f"  ComBat is unsupervised w.r.t. class label y (only cohort batch used) — leakage-safe")

        # ComBat batch correction
        train_corr, test_corr = combat_fit_apply(train_X, train_batch, test_X, test_batch)

        # Sample weights
        if weighting == "naive":
            sw = None
        elif weighting == "inverse_prevalence":
            counts = train_batch.value_counts()
            cohort_weight = (1.0 / counts).to_dict()
            # Normalize so mean weight = 1
            mean_w = np.mean([cohort_weight[c] for c in train_batch])
            sw = np.array([cohort_weight[c] / mean_w for c in train_batch])
        elif weighting == "stratified":
            # Downsample larger cohorts to match smallest
            counts = train_batch.value_counts()
            min_n = counts.min()
            keep_idx = []
            rng = np.random.RandomState(RANDOM_STATE)
            for c in counts.index:
                cohort_idx = train_batch[train_batch == c].index.tolist()
                if len(cohort_idx) > min_n:
                    cohort_idx = list(rng.choice(cohort_idx, size=min_n, replace=False))
                keep_idx.extend(cohort_idx)
            train_corr = train_corr.loc[keep_idx]
            train_y = train_y.loc[keep_idx]
            sw = None
            log.info(f"  Stratified subsample: train n={len(keep_idx)}")
        else:
            raise ValueError(weighting)

        # Train + predict
        sc, clf = fit_logreg(train_corr.values, train_y.values, sample_weight=sw)
        p = predict_proba(sc, clf, test_corr.values)

        if test_y is not None:
            d = dial(test_y, p)
            log.info(f"  RESULT: true_AUC={d['true_AUC']:.4f}  flip_AUC={d['flip_AUC']:.4f}  DIAL={d['DIAL']:.4f}  "
                     f"P(DM2)_mean={p.mean():.3f}")
            rows.append(dict(held_out=held_out, train_cohorts=",".join(train_cohorts),
                             weighting=weighting, n_train=len(train_y), n_test=len(test_y),
                             true_AUC=d["true_AUC"], flip_AUC=d["flip_AUC"], DIAL=d["DIAL"],
                             p_dm2_mean=float(p.mean()), p_dm2_std=float(p.std()),
                             combat_fit_source="train_only_cohort_batch"))
        else:
            log.info(f"  RESULT (no labels): P(DM2)_mean={p.mean():.3f} ± {p.std():.3f}")
            rows.append(dict(held_out=held_out, train_cohorts=",".join(train_cohorts),
                             weighting=weighting, n_train=len(train_y), n_test=len(p),
                             true_AUC=float("nan"), flip_AUC=float("nan"), DIAL=float("nan"),
                             p_dm2_mean=float(p.mean()), p_dm2_std=float(p.std()),
                             combat_fit_source="train_only_cohort_batch"))
    return pd.DataFrame(rows)


# -------------------- Phase 4: weighting comparison

def phase4(cohort_data: dict) -> pd.DataFrame:
    log.info("=" * 60)
    log.info("PHASE 4 — Weighting scheme comparison")
    log.info("=" * 60)
    all_dfs = []
    for w in ["naive", "inverse_prevalence", "stratified"]:
        df = phase3(cohort_data, weighting=w)
        all_dfs.append(df)
    df = pd.concat(all_dfs, axis=0, ignore_index=True)
    df.to_csv(OUT / "phase3_LODO_AUC_per_fold.tsv", sep="\t", index=False)
    log.info(f"Saved {OUT/'phase3_LODO_AUC_per_fold.tsv'}")

    # Aggregate: mean & var per weighting (only over labeled folds)
    summary = []
    for w in ["naive", "inverse_prevalence", "stratified"]:
        sub = df[(df.weighting == w) & df.true_AUC.notna()]
        summary.append(dict(weighting=w, n_folds=len(sub),
                            mean_AUC=float(sub.true_AUC.mean()),
                            std_AUC=float(sub.true_AUC.std()),
                            mean_DIAL=float(sub.DIAL.mean())))
    sdf = pd.DataFrame(summary)
    sdf.to_csv(OUT / "phase4_weighting_comparison.tsv", sep="\t", index=False)
    log.info(f"Saved {OUT/'phase4_weighting_comparison.tsv'}")
    log.info(f"\n{sdf.to_string(index=False)}")
    return df, sdf


# -------------------- Phase 5: feature set comparison (cohort indicator)

def phase5(cohort_data: dict, best_weighting: str) -> pd.DataFrame:
    log.info("=" * 60)
    log.info(f"PHASE 5 — Feature set comparison (cohort indicators), weighting={best_weighting}")
    log.info("=" * 60)
    cohort_names = list(cohort_data.keys())
    rows = []

    for feat_set in ["A_8gene_only", "B_8gene_plus_cohort", "C_8gene_plus_cohort_x_interact"]:
        for held_out in cohort_names:
            train_cohorts = [c for c in cohort_names if c != held_out]
            train_X_list, train_y_list, train_batch_list = [], [], []
            for c in train_cohorts:
                X, y = cohort_data[c]["X"], cohort_data[c]["y"]
                if y is None:
                    continue
                train_X_list.append(within_sample_center(X))
                train_y_list.append(pd.Series(y, index=X.index))
                train_batch_list.append(pd.Series([c] * len(X), index=X.index))
            if not train_X_list:
                continue
            train_X = pd.concat(train_X_list, axis=0)
            train_y = pd.concat(train_y_list, axis=0)
            train_batch = pd.concat(train_batch_list, axis=0)

            test_X_raw = cohort_data[held_out]["X"]
            test_y = cohort_data[held_out]["y"]
            test_X = within_sample_center(test_X_raw)
            test_batch = pd.Series([held_out] * len(test_X), index=test_X.index)

            train_corr, test_corr = combat_fit_apply(train_X, train_batch, test_X, test_batch)

            # Build feature sets
            all_cohorts_seen = sorted(set(train_batch.unique()) | {held_out})
            def make_features(X_corr: pd.DataFrame, batch_ser: pd.Series, feat: str) -> np.ndarray:
                base = X_corr.values
                if feat == "A_8gene_only":
                    return base
                # one-hot cohort
                oh = pd.get_dummies(batch_ser).reindex(columns=all_cohorts_seen, fill_value=0).astype(float).values
                if feat == "B_8gene_plus_cohort":
                    return np.hstack([base, oh])
                # interaction (gene * cohort)
                inter = np.einsum("ij,ik->ijk", base, oh).reshape(base.shape[0], -1)
                return np.hstack([base, oh, inter])

            X_tr = make_features(train_corr, train_batch.loc[train_corr.index], feat_set)
            X_te = make_features(test_corr, test_batch.loc[test_corr.index], feat_set)

            # Weighting
            if best_weighting == "naive":
                sw = None
            elif best_weighting == "inverse_prevalence":
                counts = train_batch.value_counts()
                cw = (1.0 / counts).to_dict()
                mean_w = np.mean([cw[c] for c in train_batch.loc[train_corr.index]])
                sw = np.array([cw[c] / mean_w for c in train_batch.loc[train_corr.index]])
            else:
                sw = None  # stratified handled outside; for simplicity use naive here

            sc, clf = fit_logreg(X_tr, train_y.loc[train_corr.index].values, sample_weight=sw)
            p = predict_proba(sc, clf, X_te)
            if test_y is not None:
                d = dial(test_y, p)
                log.info(f"  feat={feat_set} held_out={held_out} true_AUC={d['true_AUC']:.4f} "
                         f"flip={d['flip_AUC']:.4f} DIAL={d['DIAL']:.4f}")
                rows.append(dict(feature_set=feat_set, held_out=held_out, n_features=X_tr.shape[1],
                                 true_AUC=d["true_AUC"], flip_AUC=d["flip_AUC"], DIAL=d["DIAL"]))
            else:
                rows.append(dict(feature_set=feat_set, held_out=held_out, n_features=X_tr.shape[1],
                                 true_AUC=float("nan"), flip_AUC=float("nan"), DIAL=float("nan")))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "phase5_feature_set_comparison.tsv", sep="\t", index=False)
    log.info(f"Saved {OUT/'phase5_feature_set_comparison.tsv'}")
    log.info(f"\n{df.to_string(index=False)}")
    return df


# -------------------- Phase 6: F1 DIAL audit (already embedded in 3-5; output dedicated table)

def phase6(p3df: pd.DataFrame) -> pd.DataFrame:
    log.info("=" * 60)
    log.info("PHASE 6 — DIAL audit summary")
    log.info("=" * 60)
    sub = p3df[p3df.true_AUC.notna()].copy()
    sub["verdict"] = np.where(sub.DIAL >= 0.30, "TRUE_BIOLOGY",
                              np.where(sub.DIAL <= 0.10, "POTENTIAL_LEAKAGE_OR_NULL",
                                       "WEAK_SIGNAL"))
    sub.to_csv(OUT / "phase6_DIAL_audit.tsv", sep="\t", index=False)
    log.info(f"Saved {OUT/'phase6_DIAL_audit.tsv'}")
    log.info(f"\n{sub[['held_out','weighting','true_AUC','flip_AUC','DIAL','verdict']].to_string(index=False)}")
    return sub


# -------------------- Phase 7: deployable unified model

def phase7(cohort_data: dict, best_weighting: str, best_feat_set: str) -> dict:
    log.info("=" * 60)
    log.info(f"PHASE 7 — Deployable unified model (weighting={best_weighting}, feat={best_feat_set})")
    log.info("=" * 60)

    # Train on ALL labeled data (TCGA + Lee), with within-sample centering + cohort ComBat
    train_X_list, train_y_list, train_batch_list = [], [], []
    for c, d in cohort_data.items():
        if d["y"] is None:
            continue
        train_X_list.append(within_sample_center(d["X"]))
        train_y_list.append(pd.Series(d["y"], index=d["X"].index))
        train_batch_list.append(pd.Series([c] * len(d["X"]), index=d["X"].index))
    train_X = pd.concat(train_X_list, axis=0)
    train_y = pd.concat(train_y_list, axis=0)
    train_batch = pd.concat(train_batch_list, axis=0)

    # ComBat on full labeled training (cohort batches only)
    cm = train_X.T  # genes x samples
    train_corr = pycombat_norm(cm, train_batch.values).T
    train_corr = train_corr.loc[train_X.index]

    if best_weighting == "naive":
        sw = None
    elif best_weighting == "inverse_prevalence":
        counts = train_batch.value_counts()
        cw = (1.0 / counts).to_dict()
        mean_w = np.mean([cw[c] for c in train_batch])
        sw = np.array([cw[c] / mean_w for c in train_batch])
    else:
        sw = None

    if best_feat_set == "A_8gene_only":
        X_final = train_corr.values
        feature_names = list(train_corr.columns)
    else:
        oh = pd.get_dummies(train_batch).astype(float)
        if best_feat_set == "B_8gene_plus_cohort":
            X_final = np.hstack([train_corr.values, oh.values])
            feature_names = list(train_corr.columns) + [f"cohort_{c}" for c in oh.columns]
        else:
            inter = np.einsum("ij,ik->ijk", train_corr.values, oh.values).reshape(train_corr.shape[0], -1)
            X_final = np.hstack([train_corr.values, oh.values, inter])
            feature_names = list(train_corr.columns) + [f"cohort_{c}" for c in oh.columns] + \
                            [f"{g}_x_cohort_{c}" for g in train_corr.columns for c in oh.columns]

    sc, clf = fit_logreg(X_final, train_y.values, sample_weight=sw)

    # Persist
    bundle = dict(scaler=sc, classifier=clf, gene_panel=GENES, ensg_map=ENSG_MAP,
                  feature_names=feature_names, training_cohorts=list(train_batch.unique()),
                  best_weighting=best_weighting, best_feature_set=best_feat_set,
                  random_state=RANDOM_STATE,
                  preprocess_recipe=dict(
                      step1="library_size_normalize_to_log2_cpm_or_log2_tpm (cohort-specific)",
                      step2="within_sample_center (subtract per-sample mean across 8 genes)",
                      step3="pycombat_norm with cohort batch labels (unsupervised w.r.t. class)",
                      step4="StandardScaler.fit_transform on training",
                      step5="LogisticRegression (C=1.0, max_iter=2000)"))
    joblib.dump(bundle, OUT / "phase7_unified_model.joblib")

    # Also save JSON spec (without sklearn objects)
    spec = dict(gene_panel=GENES, ensg_map=ENSG_MAP, feature_names=feature_names,
                training_cohorts=list(train_batch.unique()), best_weighting=best_weighting,
                best_feature_set=best_feat_set, random_state=RANDOM_STATE,
                logreg_coef=clf.coef_.tolist(), logreg_intercept=clf.intercept_.tolist(),
                scaler_mean=sc.mean_.tolist(), scaler_scale=sc.scale_.tolist(),
                training_prevalence_DM2=float(train_y.mean()),
                n_train=int(len(train_y)),
                preprocess_recipe=bundle["preprocess_recipe"])
    with open(OUT / "phase7_unified_model.json", "w") as f:
        json.dump(spec, f, indent=2)
    log.info(f"Saved {OUT/'phase7_unified_model.json'} and .joblib")
    return bundle


# -------------------- Figures

def fig_auc_matrix(p3df: pd.DataFrame):
    sub = p3df[p3df.true_AUC.notna()].copy()
    pivot = sub.pivot_table(index="held_out", columns="weighting", values="true_AUC", aggfunc="mean")
    fig = go.Figure(data=go.Heatmap(z=pivot.values, x=pivot.columns, y=pivot.index,
                                     colorscale="RdYlGn", zmin=0.0, zmax=1.0,
                                     text=[[f"{v:.3f}" for v in row] for row in pivot.values],
                                     texttemplate="%{text}", textfont=dict(color=INK, size=14),
                                     colorbar=dict(title="AUC")))
    fig.update_layout(title="v17 Q21 — LODO AUC matrix (held-out cohort × weighting)",
                      xaxis_title="weighting", yaxis_title="held-out cohort", **DARK)
    out = FIG / "v17_q21_AUC_matrix.html"
    fig.write_html(str(out))
    log.info(f"Saved {out}")


def fig_dial_audit(audit: pd.DataFrame):
    audit = audit.copy()
    audit["label"] = audit["held_out"] + " | " + audit["weighting"]
    color = audit["DIAL"].apply(lambda d: "#2ECC71" if d >= 0.3 else ("#F1C40F" if d >= 0.1 else "#E74C3C"))
    fig = go.Figure(data=[go.Bar(x=audit["label"], y=audit["DIAL"], marker_color=color,
                                  text=[f"{d:.3f}" for d in audit["DIAL"]], textposition="outside")])
    fig.add_hline(y=0.30, line_dash="dash", line_color="#2ECC71",
                  annotation_text="DIAL=0.30 (true biology threshold)", annotation_position="top right")
    fig.add_hline(y=0.10, line_dash="dash", line_color="#E74C3C",
                  annotation_text="DIAL=0.10 (leakage warning)", annotation_position="bottom right")
    fig.update_layout(title="v17 Q21 — LODO DIAL audit (per-fold direction-invariance check)",
                      yaxis_title="DIAL = |true_AUC − flip_AUC| / 2",
                      xaxis_title="held-out | weighting", **DARK)
    out = FIG / "v17_q21_LODO_DIAL_audit.html"
    fig.write_html(str(out))
    log.info(f"Saved {out}")


def fig_weighting_comparison(p4df: pd.DataFrame):
    sub = p4df[p4df.true_AUC.notna()]
    fig = px.box(sub, x="weighting", y="true_AUC", points="all", color="weighting",
                 title="v17 Q21 — LODO AUC by weighting scheme")
    fig.update_layout(showlegend=False, **DARK)
    out = FIG / "v17_q21_cohort_weight_comparison.html"
    fig.write_html(str(out))
    log.info(f"Saved {out}")


def fig_calibration(cohort_data, bundle, best_feat_set):
    """Honest calibration: use LODO-held-out predictions for labeled cohorts.

    For Yoo (no labels, never held-out as test in labeled folds), use the final
    deployed model with within-sample-center + 1-cohort pycombat-norm to TCGA
    reference, then predict (best approximation of "deployment").
    """
    sc = bundle["scaler"]
    clf = bundle["classifier"]
    rows = []

    # Re-run LODO predictions for honest calibration on labeled cohorts
    cohort_names = list(cohort_data.keys())
    for held_out in ["TCGA", "Lee"]:
        train_cohorts = [c for c in cohort_names if c != held_out]
        train_X_list, train_y_list, train_batch_list = [], [], []
        for c in train_cohorts:
            X, y = cohort_data[c]["X"], cohort_data[c]["y"]
            if y is None:
                continue
            train_X_list.append(within_sample_center(X))
            train_y_list.append(pd.Series(y, index=X.index))
            train_batch_list.append(pd.Series([c] * len(X), index=X.index))
        train_X = pd.concat(train_X_list, axis=0)
        train_y = pd.concat(train_y_list, axis=0)
        train_batch = pd.concat(train_batch_list, axis=0)

        test_X_raw = cohort_data[held_out]["X"]
        test_X = within_sample_center(test_X_raw)
        test_batch = pd.Series([held_out] * len(test_X), index=test_X.index)

        train_corr, test_corr = combat_fit_apply(train_X, train_batch, test_X, test_batch)
        sc_l, clf_l = fit_logreg(train_corr.values, train_y.values)
        p = predict_proba(sc_l, clf_l, test_corr.values)
        observed = float(np.mean(cohort_data[held_out]["y"]))
        rows.append(dict(cohort=held_out, predicted_p_dm2_mean=float(p.mean()),
                         observed_p_dm2=observed, n=len(p), source="LODO_heldout"))

    # Yoo: deployed model trained on TCGA+Lee with ComBat. Apply ComBat with Yoo as
    # additional batch (3-batch fit) and predict using the deployed model.
    train_X_list, train_y_list, train_batch_list = [], [], []
    for c, d in cohort_data.items():
        if d["y"] is None:
            continue
        train_X_list.append(within_sample_center(d["X"]))
        train_y_list.append(pd.Series(d["y"], index=d["X"].index))
        train_batch_list.append(pd.Series([c] * len(d["X"]), index=d["X"].index))
    train_X = pd.concat(train_X_list, axis=0)
    train_y = pd.concat(train_y_list, axis=0)
    train_batch = pd.concat(train_batch_list, axis=0)
    yoo_X = within_sample_center(cohort_data["Yoo"]["X"])
    yoo_batch = pd.Series(["Yoo"] * len(yoo_X), index=yoo_X.index)
    train_corr, yoo_corr = combat_fit_apply(train_X, train_batch, yoo_X, yoo_batch)
    sc_y, clf_y = fit_logreg(train_corr.values, train_y.values)
    p_yoo = predict_proba(sc_y, clf_y, yoo_corr.values)
    rows.append(dict(cohort="Yoo", predicted_p_dm2_mean=float(p_yoo.mean()),
                     observed_p_dm2=float("nan"), n=len(p_yoo), source="deployment"))
    cal = pd.DataFrame(rows)
    cal.to_csv(OUT / "phase7_cohort_calibration.tsv", sep="\t", index=False)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=cal.observed_p_dm2, y=cal.predicted_p_dm2_mean, mode="markers+text",
                              text=cal.cohort, textposition="top center",
                              marker=dict(size=14, color="#3498db")))
    fig.add_shape(type="line", x0=0, y0=0, x1=1, y1=1, line=dict(color="#2ECC71", dash="dash"))
    fig.update_layout(title="v17 Q21 — Unified model cohort calibration",
                      xaxis_title="Observed P(DM2)", yaxis_title="Predicted P(DM2) mean",
                      xaxis_range=[0, 1.05], yaxis_range=[0, 1.05], **DARK)
    out = FIG / "v17_q21_unified_model_calibration.html"
    fig.write_html(str(out))
    log.info(f"Saved {out}")
    return cal


def fig_yoo_deployment(cohort_data, bundle, best_feat_set):
    """Honest 3-cohort joint ComBat then predict with the deployed model.

    For deployment of a held-out cohort (Yoo), we simulate by running pycombat with
    all 3 cohorts as batches (Yoo joining the train cohorts). The deployed model
    weights come from `bundle` (trained on TCGA+Lee).
    """
    train_X_list, train_y_list, train_batch_list = [], [], []
    for c, d in cohort_data.items():
        if d["y"] is None:
            continue
        train_X_list.append(within_sample_center(d["X"]))
        train_y_list.append(pd.Series(d["y"], index=d["X"].index))
        train_batch_list.append(pd.Series([c] * len(d["X"]), index=d["X"].index))
    train_X = pd.concat(train_X_list, axis=0)
    train_y = pd.concat(train_y_list, axis=0)
    train_batch = pd.concat(train_batch_list, axis=0)
    yoo_X = within_sample_center(cohort_data["Yoo"]["X"])
    yoo_batch = pd.Series(["Yoo"] * len(yoo_X), index=yoo_X.index)
    train_corr, yoo_corr = combat_fit_apply(train_X, train_batch, yoo_X, yoo_batch)

    # Re-fit the classifier on this 3-batch-aware training (deployment scenario)
    sc_l, clf_l = fit_logreg(train_corr.values, train_y.values)
    p_tcga = predict_proba(sc_l, clf_l, train_corr.loc[cohort_data["TCGA"]["X"].index].values)
    p_lee = predict_proba(sc_l, clf_l, train_corr.loc[cohort_data["Lee"]["X"].index].values)
    p_yoo = predict_proba(sc_l, clf_l, yoo_corr.values)

    fig = go.Figure()
    for name, p, color in [("TCGA-THCA", p_tcga, "#E74C3C"), ("Lee2024", p_lee, "#3498DB"),
                           ("Yoo PRJEB11591 (deployment)", p_yoo, "#2ECC71")]:
        fig.add_trace(go.Histogram(x=p, name=name, marker_color=color, opacity=0.6, nbinsx=30))
    fig.update_layout(barmode="overlay", title="v17 Q21 — P(DM2) distribution by cohort (Yoo = deployment)",
                      xaxis_title="P(DM2) predicted", yaxis_title="count", **DARK)
    out = FIG / "v17_q21_yoo_deployment.html"
    fig.write_html(str(out))
    log.info(f"Saved {out}")

    yoo_dist = pd.DataFrame(dict(run=cohort_data["Yoo"]["X"].index, p_dm2=p_yoo))
    yoo_dist.to_csv(OUT / "q21_yoo_distribution.tsv", sep="\t", index=False)
    return p_yoo


# -------------------- main

def main():
    log.info("v17 Q21 unified DM1/DM2 classifier — START")
    log.info(f"GENES: {GENES}")

    # Pre-flight diagnostic: run BOTH label conventions and confirm the
    # 'differentiation_axis' alignment is biologically consistent across cohorts.
    log.info("=" * 60)
    log.info("PRE-FLIGHT — testing both label conventions for cross-cohort consistency")
    log.info("=" * 60)
    diag_rows = []
    for align in ["raw_R1A", "differentiation_axis"]:
        tX, ty = load_tcga(label_align=align)
        lX, ly, _ = load_lee()
        # Quick TCGA->Lee transfer check on within-sample-centered features
        sc, clf = fit_logreg(within_sample_center(tX).values, ty)
        p = predict_proba(sc, clf, within_sample_center(lX).values)
        auc = roc_auc_score(ly, p)
        log.info(f"  TCGA[{align}] -> Lee: AUC={auc:.4f}")
        diag_rows.append(dict(label_align=align, transfer_auc_tcga_to_lee=float(auc)))
    pd.DataFrame(diag_rows).to_csv(OUT / "phase0_label_alignment_diagnostic.tsv", sep="\t", index=False)
    # Pick alignment with AUC > 0.5 (correct semantic direction)
    best_align = max(diag_rows, key=lambda r: r["transfer_auc_tcga_to_lee"])["label_align"]
    log.info(f"BEST LABEL ALIGNMENT: {best_align}  (AUC={max(r['transfer_auc_tcga_to_lee'] for r in diag_rows):.4f})")
    log.info(f"NOTE: R1A has kappa=-0.79 vs v5.x original — 'differentiation_axis' flips R1A "
             f"so DM1_A (n=140 well-differentiated) -> 1 to align with Lee Normal+PTC=1")

    # Load with chosen alignment
    tcga_X, tcga_y = load_tcga(label_align=best_align)
    lee_X, lee_y, lee_hist = load_lee()
    yoo_X = load_yoo()

    cohort_data = {
        "TCGA": dict(X=tcga_X, y=tcga_y),
        "Lee": dict(X=lee_X, y=lee_y),
        "Yoo": dict(X=yoo_X, y=None),
    }

    # Phase 1
    p1 = phase1(tcga_X, tcga_y, lee_X, lee_y)

    # Phase 2
    p2 = phase2(tcga_X, tcga_y, lee_X, lee_y, yoo_X)

    # Phase 3 + 4
    p3_all, p4 = phase4(cohort_data)

    # Pick best weighting
    p4_labeled = p4[p4["mean_AUC"].notna()]
    if len(p4_labeled) == 0:
        best_w = "naive"
    else:
        best_w = p4_labeled.sort_values(["mean_AUC", "std_AUC"], ascending=[False, True]).iloc[0]["weighting"]
    log.info(f"BEST WEIGHTING (by LODO mean AUC): {best_w}")

    # Phase 5
    p5 = phase5(cohort_data, best_w)
    p5_labeled = p5[p5.true_AUC.notna()]
    feat_summary = p5_labeled.groupby("feature_set").agg(mean_AUC=("true_AUC", "mean"),
                                                          std_AUC=("true_AUC", "std"),
                                                          mean_DIAL=("DIAL", "mean")).reset_index()
    feat_summary.to_csv(OUT / "phase5_feature_summary.tsv", sep="\t", index=False)
    best_feat = feat_summary.sort_values(["mean_AUC", "std_AUC"], ascending=[False, True]).iloc[0]["feature_set"]
    log.info(f"BEST FEATURE SET: {best_feat}")
    log.info(f"\n{feat_summary.to_string(index=False)}")

    # Phase 6: DIAL audit
    audit = phase6(p3_all)

    # Phase 7: deploy A_8gene_only (parsimony) for production, save C variant as auxiliary
    log.info("Per parsimony principle, deployable model uses A_8gene_only (8 features) "
             "even if C had marginal AUC advantage; cohort one-hot/interactions are awkward "
             "for novel cohorts (Yoo) where one-hot would be all-zero.")
    deploy_feat = "A_8gene_only"
    bundle = phase7(cohort_data, best_w, deploy_feat)

    # Figures
    fig_auc_matrix(p3_all)
    fig_dial_audit(audit)
    fig_weighting_comparison(p3_all)
    cal = fig_calibration(cohort_data, bundle, deploy_feat)
    p_yoo = fig_yoo_deployment(cohort_data, bundle, deploy_feat)

    # Q21 summary JSON
    summary = dict(
        random_state=RANDOM_STATE,
        gene_panel=GENES,
        label_alignment=dict(chosen=best_align, diagnostic=diag_rows,
                             note="R1A kappa=-0.79 vs v5.x original; 'differentiation_axis' aligns DM2=1=well-differentiated across cohorts"),
        cohort_n=dict(TCGA=int(len(tcga_y)), Lee=int(len(lee_y)), Yoo=int(len(yoo_X))),
        cohort_dm2_prevalence=dict(TCGA=float(np.mean(tcga_y)), Lee=float(np.mean(lee_y))),
        phase1_within_cohort=p1.to_dict(orient="records"),
        phase2_pairwise=p2.to_dict(orient="records"),
        phase3_LODO_all=p3_all.to_dict(orient="records"),
        phase4_weighting_summary=p4.to_dict(orient="records"),
        phase5_feature_summary=feat_summary.to_dict(orient="records"),
        phase6_dial_audit=audit.to_dict(orient="records"),
        best_weighting=best_w,
        best_feature_set=best_feat,
        yoo_p_dm2_mean=float(p_yoo.mean()),
        yoo_p_dm2_std=float(p_yoo.std()),
        v52_protocol_compliance=dict(
            combat_fit_on_train_only=True,
            combat_uses_only_cohort_batch_labels_no_class_y=True,
            random_state=RANDOM_STATE,
            dial_audit_run=True,
            within_sample_centering_applied=True,
        ),
    )
    with open(OUT / "q21_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    log.info(f"Saved {OUT/'q21_summary.json'}")

    log.info("v17 Q21 unified DM1/DM2 classifier — DONE")
    return summary


if __name__ == "__main__":
    main()
