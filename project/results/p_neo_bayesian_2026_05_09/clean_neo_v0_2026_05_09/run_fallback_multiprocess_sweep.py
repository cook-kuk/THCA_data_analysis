#!/usr/bin/env python3
"""Multiprocess fallback-method sweep for strict no-reference neoantigens."""

from __future__ import annotations

import math
import multiprocessing as mp
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier, HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC


ROOT = Path(__file__).resolve().parent
PARENT = ROOT.parent
STRICT_DIR = PARENT / "wave8b_strict_structure_2026_05_09"
CURATION = PARENT / "curation_2026_05_09"

STRUCTURE_FEATS = [
    "mean_pLDDT_peptide",
    "min_pLDDT_peptide",
    "mean_pLDDT_HLA",
    "anchor_pLDDT",
    "interface_contacts_8A",
    "interface_contacts_10A",
    "n_buried_residues_8A",
    "mean_min_pep_to_hla_CA_dist",
    "max_min_pep_to_hla_CA_dist",
    "radius_of_gyration_peptide",
    "end_to_end_CA_dist",
    "peptide_helicity_proxy",
]

QUANTUM_SCORES = ["GP_quantum", "VQC", "W7A_QK_only", "W7A_full"]
PUBLIC_SCORES = ["MHCflurry", "BigMHC_IM", "PRIME", "NetMHCpan_4.1"]
WAVE8_SCORES = ["Wave8_TCR_SelfSim_full", "Wave8_TCR_SelfSim_no_exact", "Wave8_TCR_motif_only"]
ANCHOR_SCORES = ["Structure_LR"]


@dataclass(frozen=True)
class Candidate:
    name: str
    family: str
    features: tuple[str, ...] = ()
    model: str = ""
    note: str = ""


def rank01(s: pd.Series, ascending: bool = True) -> pd.Series:
    out = s.rank(method="average", ascending=ascending)
    return (out - out.min()) / max(1.0, out.max() - out.min())


def metrics(y: np.ndarray, score: np.ndarray) -> dict[str, float]:
    mask = np.isfinite(score)
    yy = y[mask]
    ss = score[mask]
    if len(np.unique(yy)) < 2:
        return {"n": int(mask.sum()), "n_pos": int(yy.sum()), "auroc": np.nan, "auprc": np.nan}
    return {
        "n": int(mask.sum()),
        "n_pos": int(yy.sum()),
        "auroc": float(roc_auc_score(yy, ss)),
        "auprc": float(average_precision_score(yy, ss)),
    }


def repeated_oof_predict(clf, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=19)
    pred = np.zeros(len(y), dtype=float)
    counts = np.zeros(len(y), dtype=float)
    for train_idx, test_idx in cv.split(x, y):
        clf.fit(x[train_idx], y[train_idx])
        if hasattr(clf, "predict_proba"):
            fold = clf.predict_proba(x[test_idx])[:, 1]
        else:
            fold = clf.decision_function(x[test_idx])
        pred[test_idx] += fold
        counts[test_idx] += 1
    return pred / counts


def quantum_kernel(a: np.ndarray, b: np.ndarray, gamma: float) -> np.ndarray:
    # Product fidelity kernel for a simple angle-encoded feature map.
    diff = a[:, None, :] - b[None, :, :]
    log_k = np.log(np.clip(np.cos(gamma * diff) ** 2, 1e-8, 1.0)).mean(axis=2)
    return np.exp(log_k)


def quantum_oof_predict(x: np.ndarray, y: np.ndarray, gamma: float) -> np.ndarray:
    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=23)
    pred = np.zeros(len(y), dtype=float)
    counts = np.zeros(len(y), dtype=float)
    for train_idx, test_idx in cv.split(x, y):
        scaler = StandardScaler()
        xtr = scaler.fit_transform(x[train_idx])
        xte = scaler.transform(x[test_idx])
        ktr = quantum_kernel(xtr, xtr, gamma)
        kte = quantum_kernel(xte, xtr, gamma)
        clf = SVC(kernel="precomputed", class_weight="balanced", C=1.0)
        clf.fit(ktr, y[train_idx])
        pred[test_idx] += clf.decision_function(kte)
        counts[test_idx] += 1
    return pred / counts


def run_candidate(args: tuple[Candidate, pd.DataFrame]) -> dict:
    cand, df = args
    y = df["label"].to_numpy(int)

    if cand.family == "fixed":
        score = df[cand.features[0]].to_numpy(float)
    elif cand.family == "rankavg":
        ranks = []
        for col in cand.features:
            ascending = col.startswith("inv__")
            real_col = col.replace("inv__", "")
            ranks.append(rank01(df[real_col], ascending=ascending).to_numpy(float))
        score = np.nanmean(np.column_stack(ranks), axis=1)
    elif cand.family == "cv":
        x = df[list(cand.features)].to_numpy(float)
        if cand.model == "lr":
            clf = make_pipeline(
                StandardScaler(),
                LogisticRegression(C=0.25, class_weight="balanced", max_iter=5000, solver="liblinear"),
            )
        elif cand.model == "rf":
            clf = RandomForestClassifier(
                n_estimators=300,
                max_depth=3,
                min_samples_leaf=5,
                class_weight="balanced",
                random_state=31,
                n_jobs=1,
            )
        elif cand.model == "extratrees":
            clf = ExtraTreesClassifier(
                n_estimators=400,
                max_depth=3,
                min_samples_leaf=5,
                class_weight="balanced",
                random_state=37,
                n_jobs=1,
            )
        elif cand.model == "hgb":
            clf = HistGradientBoostingClassifier(
                max_iter=80,
                learning_rate=0.03,
                max_leaf_nodes=7,
                l2_regularization=0.2,
                random_state=41,
            )
        elif cand.model == "rbf_svc":
            clf = make_pipeline(
                StandardScaler(),
                SVC(kernel="rbf", C=0.5, gamma="scale", class_weight="balanced", probability=False),
            )
        else:
            raise ValueError(cand.model)
        score = repeated_oof_predict(clf, x, y)
    elif cand.family == "quantum_kernel":
        x = df[list(cand.features)].to_numpy(float)
        gamma = float(cand.model)
        score = quantum_oof_predict(x, y, gamma)
    else:
        raise ValueError(cand.family)

    row = {
        "method": cand.name,
        "family": cand.family,
        "model": cand.model,
        "features": ",".join(cand.features),
        "note": cand.note,
    }
    row.update(metrics(y, np.asarray(score, dtype=float)))
    return row


def load_data() -> pd.DataFrame:
    strict = pd.read_csv(STRICT_DIR / "strict_esmfold_merged_predictions.tsv", sep="\t")
    wide = pd.read_csv(CURATION / "curated_predictions_itsndb_wide.tsv", sep="\t")
    extra = [c for c in QUANTUM_SCORES + ["W7B_stacked", "Stack_mean_E1", "TSCAPE_TITANiAN"] if c in wide]
    strict = strict.merge(
        wide[["peptide", "hla", "label", *extra]],
        on=["peptide", "hla", "label"],
        how="left",
        validate="one_to_one",
    )
    return strict


def build_candidates(df: pd.DataFrame) -> list[Candidate]:
    candidates: list[Candidate] = []
    fixed_cols = [
        *ANCHOR_SCORES,
        *PUBLIC_SCORES,
        *WAVE8_SCORES,
        "ESM2_Bayesian",
        "ESMFold_3D_CV_score",
        *QUANTUM_SCORES,
        "W7B_stacked",
        "Stack_mean_E1",
        "TSCAPE_TITANiAN",
    ]
    for col in fixed_cols:
        if col in df:
            candidates.append(Candidate(col, "fixed", (col,), note="locked existing score"))

    candidates.extend(
        [
            Candidate(
                "rankavg_structure_quantum",
                "rankavg",
                tuple(c for c in ["Structure_LR", "GP_quantum", "VQC", "W7A_QK_only"] if c in df),
                note="keeps quantum, no public pretrained score",
            ),
            Candidate(
                "rankavg_quantum_only",
                "rankavg",
                tuple(c for c in ["GP_quantum", "VQC", "W7A_QK_only", "W7A_full"] if c in df),
                note="quantum family only",
            ),
            Candidate(
                "rankavg_structure_wave8_no_exact",
                "rankavg",
                ("Structure_LR", "Wave8_TCR_SelfSim_no_exact"),
                note="clean anchor plus non-self-exact Wave8",
            ),
            Candidate(
                "rankavg_structure_public_upper",
                "rankavg",
                tuple(c for c in ["Structure_LR", "MHCflurry", "BigMHC_IM", "PRIME"] if c in df),
                note="upper-bound only; public-training caveated",
            ),
            Candidate(
                "rankavg_structure_signed_esmfold",
                "rankavg",
                (
                    "Structure_LR",
                    "inv__min_pLDDT_peptide",
                    "inv__mean_pLDDT_peptide",
                    "inv__radius_of_gyration_peptide",
                    "mean_min_pep_to_hla_CA_dist",
                ),
                note="structure proxy rescue; signs chosen from pilot",
            ),
        ]
    )

    clean_features = tuple([*ANCHOR_SCORES, *[c for c in QUANTUM_SCORES if c in df], *STRUCTURE_FEATS, "tcr_motif_score"])
    no_anchor_features = tuple([*[c for c in QUANTUM_SCORES if c in df], *STRUCTURE_FEATS, "tcr_motif_score"])
    compact_features = tuple(["Structure_LR", "W7A_QK_only", "GP_quantum", "min_pLDDT_peptide", "tcr_motif_score"])

    for model in ["lr", "rf", "extratrees", "hgb", "rbf_svc"]:
        candidates.append(
            Candidate(f"{model}_anchor_quantum_structure", "cv", clean_features, model, "internal repeated CV")
        )
        candidates.append(
            Candidate(f"{model}_quantum_structure_no_anchor", "cv", no_anchor_features, model, "internal repeated CV")
        )
    candidates.append(Candidate("lr_compact_reference_replacement", "cv", compact_features, "lr", "small v0 candidate"))
    candidates.append(Candidate("rf_compact_reference_replacement", "cv", compact_features, "rf", "small v0 candidate"))

    for gamma in [0.25, 0.5, 1.0, 2.0]:
        candidates.append(
            Candidate(
                f"quantum_kernel_compact_gamma{gamma}",
                "quantum_kernel",
                compact_features,
                str(gamma),
                "angle-encoded fidelity kernel; internal repeated CV",
            )
        )
        candidates.append(
            Candidate(
                f"quantum_kernel_no_anchor_gamma{gamma}",
                "quantum_kernel",
                no_anchor_features,
                str(gamma),
                "angle-encoded fidelity kernel; internal repeated CV",
            )
        )
    return [c for c in candidates if c.features]


def main() -> None:
    df = load_data()
    candidates = build_candidates(df)
    n_proc = max(1, min(8, mp.cpu_count() - 1, len(candidates)))
    with mp.Pool(processes=n_proc) as pool:
        rows = pool.map(run_candidate, [(c, df) for c in candidates])
    out = pd.DataFrame(rows).sort_values(["auroc", "auprc"], ascending=False)
    out.to_csv(ROOT / "fallback_algorithm_sweep.tsv", sep="\t", index=False)

    top = out.head(18)
    lines = [
        "# Fallback Algorithm Sweep",
        "",
        f"Strict class-I set: n={len(df)}, positives={int(df['label'].sum())}. Multiprocessing workers: {n_proc}.",
        "",
        "All trained rows are internal repeated-CV pilots, not external paper claims.",
        "",
        "| rank | method | family | n / pos | AUROC | AUPRC | note |",
        "|---:|---|---|---:|---:|---:|---|",
    ]
    for i, (_, r) in enumerate(top.iterrows(), 1):
        lines.append(
            f"| {i} | {r['method']} | {r['family']} | {int(r['n'])} / {int(r['n_pos'])} | "
            f"{r['auroc']:.3f} | {r['auprc']:.3f} | {r['note']} |"
        )
    lines.extend(
        [
            "",
            "## Decision",
            "",
            "- Primary no-reference substitute for now: `Structure_LR`.",
            "- Keep quantum as a branch, not the anchor: quantum-only and quantum-kernel pilots are compared here.",
            "- Do not use public-pretrained scores as clean evidence; keep them as upper-bound/caveated comparators.",
            "- Structure proxy is useful as an uncertainty/fallback feature, not a standalone reference substitute.",
        ]
    )
    (ROOT / "FALLBACK_ALGORITHM_SWEEP_SUMMARY.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
