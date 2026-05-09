#!/usr/bin/env python3
"""Shared utilities for CLEAN-NeoBench and BAR-Neo.

The helpers in this module are intentionally conservative: they prefer explicit
metadata, keep public pretrained tools caveated unless an audit is present, and
return NA-style metrics instead of crashing on sparse or single-class splits.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    from sklearn.isotonic import IsotonicRegression
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
    from sklearn.model_selection import StratifiedKFold, train_test_split
    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier

    SKLEARN_AVAILABLE = True
except Exception:  # pragma: no cover - exercised only in minimal environments
    IsotonicRegression = None
    LogisticRegression = None
    RandomForestClassifier = None
    GradientBoostingClassifier = None
    StratifiedKFold = None
    train_test_split = None
    average_precision_score = None
    brier_score_loss = None
    roc_auc_score = None
    SKLEARN_AVAILABLE = False


SEED = 20260509
PUBLIC_METHOD_TOKENS = (
    "mhcflurry",
    "netmhc",
    "bigmhc",
    "prime",
    "mixmhc",
    "netmhcstab",
    "immunostruct",
    "neomhci",
    "unipmt",
    "deephla",
    "deephimmuno",
    "deepimmuno",
    "transphla",
    "mhcnuggets",
    "tscape",
    "titanian",
)
KOREAN_HLA_ALLELES = {
    "HLA-A*24:02",
    "HLA-A*11:01",
    "HLA-A*02:01",
    "HLA-B*15:01",
    "HLA-B*40:01",
    "HLA-C*01:02",
    "HLA-C*03:03",
    "HLA-C*07:02",
}
AA = "ACDEFGHIKLMNPQRSTVWY"
HYDROPHOBIC = set("AILMFWVY")
AROMATIC = set("FWY")
CHARGED_POS = set("KRH")
CHARGED_NEG = set("DE")


MASTER_COLUMNS = [
    "candidate_id",
    "sample_id",
    "patient_id",
    "study_id",
    "source_name",
    "source_dataset",
    "cancer_type",
    "disease_context",
    "peptide",
    "mut_peptide",
    "wt_peptide",
    "peptide_length",
    "hla",
    "hla_gene",
    "hla_allele_4digit",
    "hla_supertype",
    "gene",
    "protein_id",
    "mutation_id",
    "source_protein_window",
    "label",
    "label_type",
    "mhc_class",
    "expression_tpm",
    "mutant_expression",
    "vaf",
    "clonality",
    "hla_loh",
    "b2m_status",
    "antigen_processing_status",
    "tumor_stage",
    "treatment_context",
    "immune_context_score",
    "tls_score",
    "ifng_score",
    "cytolytic_score",
    "exact_peptide_train_overlap",
    "exact_peptide_hla_train_overlap",
    "near_peptide_train_overlap",
    "source_protein_window_train_overlap",
    "study_train_overlap",
    "patient_train_overlap",
    "public_tool_training_overlap_any",
    "public_tool_training_overlap_detail",
    "leakage_risk_level",
    "split_exact_phla",
    "split_near_peptide_cluster",
    "split_source_heldout",
    "split_hla_heldout",
    "split_supertype_heldout",
    "split_patient_heldout",
    "split_study_heldout",
    "split_low_prevalence",
    "split_korean_hla_focus",
]

METHOD_SCORE_COLUMNS = [
    "candidate_id",
    "method_name",
    "method_family",
    "method_role",
    "score_raw",
    "score_direction",
    "score_calibrated",
    "rank_global",
    "rank_within_patient",
    "rank_within_hla",
    "runtime_status",
    "caveat",
    "uses_public_pretraining",
    "training_overlap_audited",
    "clean_comparator_allowed",
]


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def read_table(path: Path, nrows: int | None = None) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep, nrows=nrows)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    df.to_csv(path, sep="\t", index=False, na_rep="NA")


def update_manifest(output_root: Path, stage: str, payload: dict[str, Any]) -> None:
    ensure_dir(output_root)
    path = output_root / "run_manifest.json"
    if path.exists():
        try:
            manifest = json.loads(path.read_text())
        except Exception:
            manifest = {}
    else:
        manifest = {}
    manifest.setdefault("stages", {})
    manifest["stages"][stage] = payload
    warnings = payload.get("warnings", [])
    if warnings:
        manifest.setdefault("warnings", [])
        manifest["warnings"].extend(str(w) for w in warnings)
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def stable_hash(value: str, modulo: int = 5) -> int:
    digest = hashlib.sha1(str(value).encode("utf-8")).hexdigest()
    return int(digest[:10], 16) % modulo


def normalize_empty(value: Any) -> str:
    if value is None:
        return ""
    s = str(value).strip()
    if s.lower() in {"", "nan", "none", "na", "n/a", "null"}:
        return ""
    return s


def normalize_hla(value: Any) -> str:
    s = normalize_empty(value).replace("HLA_", "HLA-").replace("_", "*", 1)
    if not s:
        return ""
    if s.startswith("HLA-"):
        return s
    if re.match(r"^[ABC]\*\d{2}:\d{2}", s):
        return f"HLA-{s}"
    if re.match(r"^[ABC]\d{4}$", s):
        return f"HLA-{s[0]}*{s[1:3]}:{s[3:5]}"
    return s


def hla_gene(hla: Any) -> str:
    h = normalize_hla(hla)
    if h.startswith("HLA-") and "*" in h:
        return h.split("*", 1)[0]
    return ""


def hla_allele_4digit(hla: Any) -> str:
    h = normalize_hla(hla)
    m = re.match(r"^(HLA-[A-Z]+)\*(\d{2}):?(\d{2})", h)
    if not m:
        return h
    return f"{m.group(1)}*{m.group(2)}:{m.group(3)}"


def hla_supertype(hla: Any) -> str:
    h = hla_allele_4digit(hla)
    if not h.startswith("HLA-") or "*" not in h:
        return ""
    locus = h.split("*", 1)[0].replace("HLA-", "")
    fam = h.split("*", 1)[1].split(":", 1)[0]
    return f"{locus}{fam}"


def peptide_cluster_fallback(peptide: Any) -> str:
    pep = normalize_empty(peptide).upper()
    if not pep:
        return ""
    if len(pep) < 4:
        return f"L{len(pep)}_{pep}"
    return f"L{len(pep)}_{pep[:2]}_{pep[-2:]}"


def kmer_jaccard(a: str, b: str, k: int = 3) -> float:
    a = normalize_empty(a)
    b = normalize_empty(b)
    if not a or not b:
        return 0.0
    if len(a) < k or len(b) < k:
        return 1.0 if a == b else 0.0
    sa = {a[i : i + k] for i in range(len(a) - k + 1)}
    sb = {b[i : i + k] for i in range(len(b) - k + 1)}
    return len(sa & sb) / max(1, len(sa | sb))


def normalized_identity(a: str, b: str) -> float:
    a = normalize_empty(a)
    b = normalize_empty(b)
    if not a or not b:
        return 0.0
    matches = sum(x == y for x, y in zip(a, b))
    return matches / max(len(a), len(b))


def near_similarity(a: str, b: str) -> float:
    return max(normalized_identity(a, b), kmer_jaccard(a, b, 3), kmer_jaccard(a, b, 4))


def method_family_and_role(method_name: str) -> tuple[str, str, bool, bool, bool, str]:
    name = normalize_empty(method_name)
    low = name.lower()
    uses_public = any(tok in low for tok in PUBLIC_METHOD_TOKENS)
    audited = False
    clean_allowed = False
    caveat = ""
    if name == "Structure_LR":
        return "structure_proxy", "anchor", False, True, True, "current honest local anchor"
    if "esm2" in low and "bayesian" in low:
        return "bayesian_uncertainty", "uncertainty_only", False, True, True, "ranking not promoted; uncertainty/OOD branch"
    if "wave8" in low or "tcr" in low:
        return "tcr_aware", "internal_candidate", False, True, True, "reference sensitivity requires audit"
    if low.startswith("qk") or "qk" in low or "quantum" in low or low in {"vqc", "gp_quantum"}:
        return "qk_fallback", "bounded_fallback", False, True, True, "bounded fallback/fusion component only"
    if "anchor" in low or "counterfactual" in low or name in {"RF_biophys", "LR_8d_classical"}:
        return "internal_baseline", "anchor" if "anchor" in low else "internal_candidate", False, True, True, "internal model"
    if uses_public:
        return "deep_immunogenicity", "caveated_public_comparator", True, False, False, "public pretrained comparator; training overlap unresolved"
    if "foldsafe" in low or "fusion" in low or "stack" in low or "ensemble" in low:
        return "ensemble_or_fusion", "internal_candidate", False, True, True, "internal fusion/ensemble"
    if "plddt" in low or "esmfold" in low or "structure" in low:
        return "structure_proxy", "internal_candidate", False, True, True, "structure proxy"
    if "gp" in low or "vqc" in low:
        return "qk_fallback", "bounded_fallback", False, True, True, "bounded fallback/fusion component only"
    return "classical_biophysical", "internal_candidate", False, True, True, "internal or locally curated method"


def detect_score_direction(y: pd.Series, score: pd.Series) -> str:
    yv = pd.to_numeric(y, errors="coerce")
    sv = pd.to_numeric(score, errors="coerce")
    mask = yv.notna() & sv.notna()
    if mask.sum() < 6 or yv[mask].nunique() < 2:
        return "higher_better"
    arr_y = yv[mask].astype(int).to_numpy()
    arr_s = sv[mask].astype(float).to_numpy()
    try:
        if SKLEARN_AVAILABLE:
            auc_hi = roc_auc_score(arr_y, arr_s)
            auc_lo = roc_auc_score(arr_y, -arr_s)
            return "lower_better" if auc_lo > auc_hi + 0.02 else "higher_better"
    except Exception:
        pass
    pos = np.nanmean(arr_s[arr_y == 1])
    neg = np.nanmean(arr_s[arr_y == 0])
    return "lower_better" if pos < neg else "higher_better"


def percentile_calibrate(score: pd.Series) -> pd.Series:
    s = pd.to_numeric(score, errors="coerce")
    if s.notna().sum() == 0:
        return pd.Series(np.nan, index=score.index)
    return s.rank(method="average", pct=True)


def fold_safe_calibrate(y: pd.Series, score: pd.Series) -> pd.Series:
    """Return out-of-fold calibrated probabilities when feasible.

    If the sample is too small or sklearn is unavailable, return rank
    percentiles. This avoids using test labels for a polished probability.
    """

    yv = pd.to_numeric(y, errors="coerce")
    sv = pd.to_numeric(score, errors="coerce")
    out = pd.Series(np.nan, index=score.index, dtype=float)
    mask = yv.notna() & sv.notna()
    if mask.sum() < 20 or yv[mask].nunique() < 2 or not SKLEARN_AVAILABLE:
        out.loc[mask] = percentile_calibrate(sv[mask])
        return out
    idx = np.array(mask[mask].index)
    yy = yv.loc[idx].astype(int).to_numpy()
    ss = sv.loc[idx].astype(float).to_numpy()
    n_splits = min(5, int(np.bincount(yy).min()))
    if n_splits < 2:
        out.loc[idx] = percentile_calibrate(pd.Series(ss, index=idx))
        return out
    try:
        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=SEED)
        for train_i, test_i in skf.split(ss.reshape(-1, 1), yy):
            train_s = ss[train_i]
            train_y = yy[train_i]
            test_s = ss[test_i]
            if len(np.unique(train_y)) < 2:
                pred = pd.Series(test_s).rank(pct=True).to_numpy()
            elif len(train_y) >= 30:
                iso = IsotonicRegression(out_of_bounds="clip")
                iso.fit(train_s, train_y)
                pred = iso.predict(test_s)
            else:
                lr = LogisticRegression(max_iter=500)
                lr.fit(train_s.reshape(-1, 1), train_y)
                pred = lr.predict_proba(test_s.reshape(-1, 1))[:, 1]
            out.loc[idx[test_i]] = pred
    except Exception:
        out.loc[idx] = percentile_calibrate(pd.Series(ss, index=idx))
    return out.clip(0, 1)


def ece_score(y: np.ndarray, score: np.ndarray, n_bins: int = 10) -> float:
    y = np.asarray(y, dtype=float)
    s = np.clip(np.asarray(score, dtype=float), 0, 1)
    mask = np.isfinite(y) & np.isfinite(s)
    y = y[mask]
    s = s[mask]
    if len(y) == 0:
        return math.nan
    bins = np.linspace(0, 1, n_bins + 1)
    ece = 0.0
    for lo, hi in zip(bins[:-1], bins[1:]):
        bmask = (s >= lo) & (s < hi if hi < 1 else s <= hi)
        if bmask.sum() == 0:
            continue
        ece += float(bmask.mean() * abs(s[bmask].mean() - y[bmask].mean()))
    return ece


def binary_metrics(y: pd.Series, score: pd.Series, k_values: tuple[int, ...] = (1, 5, 10, 20)) -> dict[str, Any]:
    yv = pd.to_numeric(y, errors="coerce")
    sv = pd.to_numeric(score, errors="coerce")
    mask = yv.notna() & sv.notna()
    yy = yv[mask].astype(int).to_numpy()
    ss = sv[mask].astype(float).to_numpy()
    out: dict[str, Any] = {
        "n_total": int(len(yy)),
        "n_pos": int(yy.sum()) if len(yy) else 0,
        "prevalence": float(yy.mean()) if len(yy) else math.nan,
        "metric_note": "",
    }
    if len(yy) == 0:
        out.update({"AUROC": math.nan, "AUPRC": math.nan, "calibration_brier": math.nan, "calibration_ece": math.nan})
        out["metric_note"] = "no scorable rows"
    elif len(np.unique(yy)) < 2:
        out.update({"AUROC": math.nan, "AUPRC": math.nan, "calibration_brier": math.nan, "calibration_ece": math.nan})
        out["metric_note"] = "single-class split; AUROC/AUPRC unavailable"
    elif SKLEARN_AVAILABLE:
        clipped = np.clip(ss, 0, 1)
        out.update(
            {
                "AUROC": float(roc_auc_score(yy, ss)),
                "AUPRC": float(average_precision_score(yy, ss)),
                "calibration_brier": float(brier_score_loss(yy, clipped)),
                "calibration_ece": float(ece_score(yy, clipped)),
            }
        )
    else:
        out.update({"AUROC": math.nan, "AUPRC": fallback_auprc(yy, ss), "calibration_brier": math.nan, "calibration_ece": math.nan})
        out["metric_note"] = "sklearn unavailable; AUROC/calibration unavailable"

    order = np.argsort(-ss) if len(ss) else np.array([], dtype=int)
    pos_ranks = []
    for rank, idx in enumerate(order, start=1):
        if yy[idx] == 1:
            pos_ranks.append(rank)
    out["mean_positive_rank"] = float(np.mean(pos_ranks)) if pos_ranks else math.nan
    out["median_positive_rank"] = float(np.median(pos_ranks)) if pos_ranks else math.nan
    for k in k_values:
        kk = min(k, len(yy))
        if kk == 0:
            hit = precision = recall = math.nan
        else:
            top = yy[order[:kk]]
            hit = bool(top.sum() > 0)
            precision = float(top.mean())
            recall = float(top.sum() / max(1, yy.sum()))
        out[f"top{k}_hit"] = hit
        out[f"top{k}_precision"] = precision
        out[f"top{k}_recall"] = recall

    if len(ss):
        confidence = np.abs(np.clip(ss, 0, 1) - 0.5) * 2
        covered = confidence >= 0.8
        out["coverage_at_confidence_80"] = float(covered.mean())
        if covered.sum() and len(np.unique(yy)) >= 2:
            pred = (np.clip(ss[covered], 0, 1) >= 0.5).astype(int)
            out["risk_at_coverage_80"] = float((pred != yy[covered]).mean())
        else:
            out["risk_at_coverage_80"] = math.nan
    else:
        out["coverage_at_confidence_80"] = math.nan
        out["risk_at_coverage_80"] = math.nan
    return out


def fallback_auprc(y: np.ndarray, score: np.ndarray) -> float:
    if len(y) == 0 or y.sum() == 0:
        return math.nan
    order = np.argsort(-score)
    y_sorted = y[order]
    tp = 0
    precisions = []
    for i, val in enumerate(y_sorted, start=1):
        if val == 1:
            tp += 1
            precisions.append(tp / i)
    return float(np.mean(precisions)) if precisions else math.nan


def sequence_features(peptide: str) -> dict[str, float]:
    pep = normalize_empty(peptide).upper()
    n = max(1, len(pep))
    feats = {
        "peptide_length_numeric": float(len(pep)),
        "hydrophobic_fraction": sum(a in HYDROPHOBIC for a in pep) / n,
        "aromatic_fraction": sum(a in AROMATIC for a in pep) / n,
        "charge_proxy": (sum(a in CHARGED_POS for a in pep) - sum(a in CHARGED_NEG for a in pep)) / n,
        "cysteine_count": float(pep.count("C")),
        "glycine_proline_count": float(pep.count("G") + pep.count("P")),
    }
    if 8 <= len(pep) <= 11:
        feats["anchor_nterm_hydrophobic"] = float(pep[1] in HYDROPHOBIC) if len(pep) > 1 else 0.0
        feats["anchor_cterm_hydrophobic"] = float(pep[-1] in HYDROPHOBIC) if pep else 0.0
    else:
        feats["anchor_nterm_hydrophobic"] = 0.0
        feats["anchor_cterm_hydrophobic"] = 0.0
    return feats


def value_counts_md(series: pd.Series, limit: int = 20) -> str:
    vc = series.fillna("NA").astype(str).value_counts().head(limit)
    if vc.empty:
        return "No rows."
    return "\n".join(f"- `{idx}`: {int(val)}" for idx, val in vc.items())


def dataframe_to_markdown(df: pd.DataFrame, max_rows: int = 20) -> str:
    if df is None or df.empty:
        return "No rows available."
    d = df.head(max_rows).copy()
    try:
        return d.to_markdown(index=False)
    except Exception:
        return d.to_csv(sep="\t", index=False, na_rep="NA")
