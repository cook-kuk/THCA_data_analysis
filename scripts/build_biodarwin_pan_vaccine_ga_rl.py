#!/usr/bin/env python3
"""Develop BioDarwin Pan-Vaccine: GA/RL cancer-vaccine hunting algorithm.

This is intentionally separate from the public-source scout. The scout builds
the locked industrial testbed; this script builds the actual algorithm:

* CD8/Class-I cytotoxic axis
* CD4/Class-II helper axis
* SOM niche map for search-space coverage
* GA genome for score-controller discovery
* bandit-style RL controller for mutation-operator choice
* leakage/proxy guardrails and frozen industrial holdout evaluation

Industrial public data are evaluated only after the genome is frozen.
"""

from __future__ import annotations

import json
import math
import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
ACADEMIC = (
    ROOT
    / "project/results/cross_neo_immunogenicity_algorithm_compare_2026_05_10"
    / "cross_neo_immunogenicity_comparator_scores.tsv"
)
IND_DIR = ROOT / "project/results/public_industrial_neoantigen_benchmark_scout_2026_05_10"
GA_ADAPTER = IND_DIR / "ga_public_feature_adapter_v0/industrial_locked_v0_ga_public_feature_scores.tsv"
OUT = ROOT / "project/results/biodarwin_pan_vaccine_ga_rl_2026_05_10"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")


CD8_FEATURES = [
    "bigmhc_im_norm",
    "bigmhc_el_norm",
    "cross_core_norm",
    "cross_stress_norm",
    "cross_bma_norm",
    "cross_finetuned_norm",
    "tcr_norm",
    "foreignness_norm",
    "md_norm",
    "seq_cytotoxic_prior",
    "hla_A",
    "hla_B",
    "hla_C",
]
CD4_FEATURES = [
    "seq_helper_prior",
    "helper_promiscuity_proxy",
    "slp_processability_proxy",
    "seq_diversity",
    "charged_balance",
    "aromatic",
    "hydrophobic",
]
PUBLIC_FEATURES = [
    "bigmhc_im_norm",
    "bigmhc_el_norm",
    "seq_cytotoxic_prior",
    "seq_helper_prior",
    "helper_promiscuity_proxy",
    "slp_processability_proxy",
    "seq_diversity",
    "charged_balance",
    "hydrophobic",
    "aromatic",
    "hla_A",
    "hla_B",
    "hla_C",
]
SYNERGIES = [
    "presentation_x_tcr",
    "foreignness_x_helper",
    "stress_x_md",
    "cytotoxic_x_helper",
    "claimsafe_x_helper",
]
BASELINE_SCORES = [
    ("BigMHC_IM", "bigmhc_im_score"),
    ("BigMHC_EL", "bigmhc_el_score"),
    ("CROSS_core", "crossneo_core_score_norm"),
    ("CROSS_stress", "stress_guarded_discovery_score"),
    ("CROSS_BMA", "bma_v2_discovery_score"),
    ("CROSS_finetuned", "finetuned_experiment_priority_score"),
    ("CROSS_integrated", "immunogenicity_discovery_score"),
    ("CROSS_claimsafe", "immunogenicity_claim_safe_score"),
    ("BioDarwin_PublicFallback", "biodarwin_public_fallback_score"),
    ("BioDarwin_PanVax", "biodarwin_pan_vaccine_score"),
]
TRANSFORMS = ["identity", "sqrt", "square", "sigmoid"]
OPERATORS = [
    "cd8_reweight",
    "cd4_reweight",
    "public_reweight",
    "axis_shift",
    "synergy_shift",
    "gate_shift",
    "transform_flip",
]
CD8_CAPS = {
    "tcr_norm": 0.16,
    "cross_finetuned_norm": 0.18,
    "cross_bma_norm": 0.20,
    "cross_stress_norm": 0.24,
    "md_norm": 0.16,
}
DEFAULT_CD8_CAP = 0.30


def minmax_array(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    finite = np.isfinite(x)
    out = np.zeros_like(x, dtype=float)
    if not finite.any():
        return out
    lo = float(np.min(x[finite]))
    hi = float(np.max(x[finite]))
    if hi == lo:
        out[finite] = 0.5
    else:
        out[finite] = (x[finite] - lo) / (hi - lo)
    return out


def minmax_series(s: pd.Series) -> pd.Series:
    return pd.Series(minmax_array(pd.to_numeric(s, errors="coerce").to_numpy()), index=s.index)


def normalize_weights(w: np.ndarray) -> np.ndarray:
    w = np.clip(np.asarray(w, dtype=float), 0.0, None)
    if not np.isfinite(w).all() or w.sum() <= 0:
        return np.ones(len(w), dtype=float) / len(w)
    return w / w.sum()


def cap_cd8_weights(w: np.ndarray) -> np.ndarray:
    w = normalize_weights(w)
    caps = np.array([CD8_CAPS.get(f, DEFAULT_CD8_CAP) for f in CD8_FEATURES], dtype=float)
    for _ in range(20):
        over = w > caps
        if not over.any():
            break
        excess = float((w[over] - caps[over]).sum())
        w[over] = caps[over]
        room = caps[~over] - w[~over]
        ok = room > 1e-10
        if not ok.any() or excess <= 0:
            break
        idx = np.where(~over)[0][ok]
        w[idx] += excess * room[ok] / room[ok].sum()
    return normalize_weights(np.minimum(w, caps))


def transform_values(x: np.ndarray, name: str) -> np.ndarray:
    x = np.nan_to_num(x, nan=0.0, posinf=1.0, neginf=0.0)
    x = np.clip(x, 0.0, 1.0)
    if name == "identity":
        return x
    if name == "sqrt":
        return np.sqrt(x)
    if name == "square":
        return x * x
    if name == "sigmoid":
        return 1.0 / (1.0 + np.exp(-8.0 * (x - 0.5)))
    raise ValueError(name)


def safe_ap(y: np.ndarray, s: np.ndarray) -> float:
    y = np.asarray(y, dtype=float)
    s = np.asarray(s, dtype=float)
    m = np.isfinite(y) & np.isfinite(s)
    if m.sum() == 0 or len(np.unique(y[m])) < 2:
        return np.nan
    return float(average_precision_score(y[m], s[m]))


def safe_auc(y: np.ndarray, s: np.ndarray) -> float:
    y = np.asarray(y, dtype=float)
    s = np.asarray(s, dtype=float)
    m = np.isfinite(y) & np.isfinite(s)
    if m.sum() == 0 or len(np.unique(y[m])) < 2:
        return np.nan
    return float(roc_auc_score(y[m], s[m]))


def topk_precision(y: np.ndarray, s: np.ndarray, k: int) -> float:
    y = np.asarray(y, dtype=float)
    s = np.asarray(s, dtype=float)
    m = np.isfinite(y) & np.isfinite(s)
    if m.sum() == 0:
        return np.nan
    idx = np.where(m)[0]
    idx = idx[np.argsort(s[idx])[::-1]][: min(k, len(idx))]
    return float(np.mean(y[idx])) if len(idx) else np.nan


def peptide_features(peptide: object) -> dict[str, float]:
    pep = str(peptide or "").upper()
    pep = "".join(aa for aa in pep if aa in "ACDEFGHIKLMNPQRSTVWY")
    if not pep:
        return {
            "seq_cytotoxic_prior": 0.0,
            "seq_helper_prior": 0.0,
            "helper_promiscuity_proxy": 0.0,
            "slp_processability_proxy": 0.0,
            "seq_diversity": 0.0,
            "hydrophobic": 0.0,
            "aromatic": 0.0,
            "charged_balance": 0.0,
        }
    n = len(pep)
    diversity = min(1.0, len(set(pep)) / min(12.0, max(1.0, float(n))))
    hydrophobic = sum(aa in "AILMVPFWY" for aa in pep) / n
    aromatic = sum(aa in "FWY" for aa in pep) / n
    charged = sum(aa in "KRDE" for aa in pep) / n
    charged_balance = 1.0 - min(1.0, abs(charged - 0.22) / 0.35)
    class_i_len = 1.0 if n in {9, 10, 11} else 0.72 if n in {8, 12, 13, 14, 15} else 0.20
    class_ii_len = math.exp(-abs(n - 15.0) / 5.0) if 8 <= n <= 35 else 0.0
    helper_promiscuity = 0.38 * diversity + 0.23 * charged_balance + 0.20 * hydrophobic + 0.19 * class_ii_len
    # A cheap SLP processing prior: long enough for Class-II, not too repetitive,
    # and with moderate hydrophobic/charged balance for cleavage-friendly windows.
    slp_processability = 0.35 * class_ii_len + 0.25 * diversity + 0.20 * (1.0 - abs(hydrophobic - 0.45)) + 0.20 * charged_balance
    cytotoxic = 0.40 * class_i_len + 0.20 * diversity + 0.17 * hydrophobic + 0.13 * aromatic + 0.10 * charged_balance
    helper = 0.32 * class_ii_len + 0.24 * helper_promiscuity + 0.20 * slp_processability + 0.14 * diversity + 0.10 * charged_balance
    return {
        "seq_cytotoxic_prior": float(np.clip(cytotoxic, 0, 1)),
        "seq_helper_prior": float(np.clip(helper, 0, 1)),
        "helper_promiscuity_proxy": float(np.clip(helper_promiscuity, 0, 1)),
        "slp_processability_proxy": float(np.clip(slp_processability, 0, 1)),
        "seq_diversity": float(np.clip(diversity, 0, 1)),
        "hydrophobic": float(np.clip(hydrophobic, 0, 1)),
        "aromatic": float(np.clip(aromatic, 0, 1)),
        "charged_balance": float(np.clip(charged_balance, 0, 1)),
    }


def hla_features(hla: pd.Series) -> pd.DataFrame:
    h = hla.astype(str).str.upper().str.replace("HLA-", "", regex=False)
    return pd.DataFrame(
        {
            "hla_A": h.str.startswith("A").astype(float),
            "hla_B": h.str.startswith("B").astype(float),
            "hla_C": h.str.startswith("C").astype(float),
        },
        index=hla.index,
    )


def add_common_features(df: pd.DataFrame, peptide_col: str, hla_col: str) -> pd.DataFrame:
    out = df.copy()
    pf = pd.DataFrame([peptide_features(x) for x in out[peptide_col]], index=out.index)
    hf = hla_features(out[hla_col])
    for col in pf.columns:
        out[col] = pf[col].astype(float)
    for col in hf.columns:
        out[col] = hf[col].astype(float)
    return out


def prepare_academic() -> pd.DataFrame:
    df = pd.read_csv(ACADEMIC, sep="\t")
    df = add_common_features(df, "peptide", "hla_allele_4digit")
    rename = {
        "bigmhc_im_score": "bigmhc_im_norm",
        "bigmhc_el_score": "bigmhc_el_norm",
        "crossneo_core_score_norm": "cross_core_norm",
        "stress_guarded_discovery_score": "cross_stress_norm",
        "bma_v2_discovery_score": "cross_bma_norm",
        "finetuned_experiment_priority_score": "cross_finetuned_norm",
        "tcr_recognition_score_norm": "tcr_norm",
        "fitness_foreignness_proxy": "foreignness_norm",
        "md_control_score_norm": "md_norm",
        "immunogenicity_claim_safe_score": "claimsafe_norm",
    }
    for src, dst in rename.items():
        df[dst] = minmax_series(df[src]) if src in df else 0.0
    df["label"] = df["label"].astype(int)
    df["split_biodarwin"] = np.where(
        df["source_name"].astype(str).str.contains("validation|_Val", case=False, regex=True),
        "academic_validation_like",
        "academic_train_sources",
    )
    df["low_medium_leakage"] = df["leakage_risk_level"].astype(str).str.lower().isin(["low", "medium"])
    df["internal_feature_available"] = True
    return df


def prepare_industrial() -> pd.DataFrame:
    locked = pd.read_csv(IND_DIR / "public_industrial_locked_testset_v0.tsv", sep="\t")
    im = pd.read_csv(IND_DIR / "locked_testset_bigmhc_im_predictions.csv")
    el = pd.read_csv(IND_DIR / "locked_testset_bigmhc_el_predictions.csv")
    if "candidate_id" not in im.columns:
        im.insert(0, "candidate_id", locked["industrial_candidate_id"].to_numpy())
    if "candidate_id" not in el.columns:
        el.insert(0, "candidate_id", locked["industrial_candidate_id"].to_numpy())
    df = locked.merge(
        im[["candidate_id", "BigMHC_IM"]],
        left_on="industrial_candidate_id",
        right_on="candidate_id",
        how="left",
    ).drop(columns=["candidate_id"])
    df = df.merge(
        el[["candidate_id", "BigMHC_EL"]],
        left_on="industrial_candidate_id",
        right_on="candidate_id",
        how="left",
    ).drop(columns=["candidate_id"])
    df = add_common_features(df, "peptide", "hla_allele")
    df["bigmhc_im_norm"] = minmax_series(df["BigMHC_IM"]).fillna(0.0)
    df["bigmhc_el_norm"] = minmax_series(df["BigMHC_EL"]).fillna(0.0)
    # Industrial locked rows do not have internal CROSS/TCR/MD model calls.
    # Use neutral priors for unavailable internal axes so the frozen genome can
    # still run without retraining or source-specific fitting.
    neutral_cols = [
        "cross_core_norm",
        "cross_stress_norm",
        "cross_bma_norm",
        "cross_finetuned_norm",
        "tcr_norm",
        "foreignness_norm",
        "md_norm",
        "claimsafe_norm",
    ]
    for col in neutral_cols:
        df[col] = 0.5
    df["label"] = df["label"].astype(int)
    df["split_biodarwin"] = "industrial_locked_v0"
    df["source_name"] = df["source_id"]
    df["low_medium_leakage"] = True
    df["internal_feature_available"] = False
    if GA_ADAPTER.exists():
        ga = pd.read_csv(GA_ADAPTER, sep="\t")
        df = df.merge(
            ga[["industrial_candidate_id", "ga_public_feature_adapter_score"]],
            on="industrial_candidate_id",
            how="left",
        )
    return df


def prepare_classii_scout(best_score_fn) -> pd.DataFrame:
    path = IND_DIR / "public_industrial_candidate_pairs.tsv"
    if not path.exists():
        return pd.DataFrame()
    pairs = pd.read_csv(path, sep="\t")
    if "direct_classII_score_ready" not in pairs:
        return pd.DataFrame()
    classii = pairs[pairs["direct_classII_score_ready"]].copy()
    if classii.empty:
        return classii
    classii = classii.drop_duplicates(["source_id", "peptide_or_sequence", "hla_allele"]).reset_index(drop=True)
    classii.insert(0, "classII_candidate_id", [f"INDPUBII_SCOUT_{i+1:05d}" for i in range(len(classii))])
    classii = add_common_features(classii, "peptide_or_sequence", "hla_allele")
    for col in ["bigmhc_im_norm", "bigmhc_el_norm", "cross_core_norm", "cross_stress_norm", "cross_bma_norm", "cross_finetuned_norm", "tcr_norm", "foreignness_norm", "md_norm", "claimsafe_norm"]:
        classii[col] = 0.5
    classii["internal_feature_available"] = False
    scored = best_score_fn(classii)
    classii["biodarwin_cd8_axis_score"] = scored["cd8"]
    classii["biodarwin_cd4_helper_score"] = scored["cd4"]
    classii["biodarwin_full_mode_score"] = scored["full_pan"]
    classii["biodarwin_public_fallback_score"] = scored["public_fallback"]
    classii["biodarwin_pan_vaccine_score"] = scored["pan"]
    return classii.sort_values("biodarwin_cd4_helper_score", ascending=False)


class SimpleSOM:
    def __init__(self, width: int = 8, height: int = 8, seed: int = 17) -> None:
        self.width = width
        self.height = height
        self.seed = seed
        self.weights: np.ndarray | None = None

    def fit(self, x: np.ndarray, n_iter: int = 1800) -> "SimpleSOM":
        rng = np.random.default_rng(self.seed)
        x = np.nan_to_num(np.asarray(x, dtype=float), nan=0.0)
        n_features = x.shape[1]
        self.weights = rng.random((self.width, self.height, n_features))
        coords = np.array([(i, j) for i in range(self.width) for j in range(self.height)])
        for t in range(n_iter):
            xi = x[rng.integers(0, len(x))]
            bmu = self.bmu(xi)
            frac = 1.0 - t / max(1, n_iter)
            lr = 0.35 * frac + 0.03
            radius = max(0.8, 3.8 * frac)
            dist2 = np.sum((coords - np.array(bmu)) ** 2, axis=1).reshape(self.width, self.height)
            influence = np.exp(-dist2 / (2.0 * radius * radius))[:, :, None]
            self.weights += lr * influence * (xi - self.weights)
        return self

    def bmu(self, xi: np.ndarray) -> tuple[int, int]:
        if self.weights is None:
            raise RuntimeError("SOM is not fitted")
        d = np.sum((self.weights - xi) ** 2, axis=2)
        idx = int(np.argmin(d))
        return idx // self.height, idx % self.height

    def transform(self, x: np.ndarray) -> np.ndarray:
        return np.array([self.bmu(row) for row in np.nan_to_num(x, nan=0.0)], dtype=int)


@dataclass
class Genome:
    cd8_weights: np.ndarray
    cd4_weights: np.ndarray
    public_weights: np.ndarray
    cd8_transforms: np.ndarray
    cd4_transforms: np.ndarray
    public_transforms: np.ndarray
    synergy_weights: np.ndarray
    cd8_mix: float
    full_mode_mix: float
    low_presentation_gate: float
    low_presentation_penalty: float
    high_leakage_penalty: float
    medium_leakage_penalty: float
    complexity_penalty: float


def random_genome(rng: np.random.Generator) -> Genome:
    return Genome(
        cd8_weights=cap_cd8_weights(rng.dirichlet(np.ones(len(CD8_FEATURES)) * 0.7)),
        cd4_weights=normalize_weights(rng.dirichlet(np.ones(len(CD4_FEATURES)) * 0.9)),
        public_weights=normalize_weights(rng.dirichlet(np.ones(len(PUBLIC_FEATURES)) * 0.9)),
        cd8_transforms=rng.integers(0, len(TRANSFORMS), len(CD8_FEATURES)),
        cd4_transforms=rng.integers(0, len(TRANSFORMS), len(CD4_FEATURES)),
        public_transforms=rng.integers(0, len(TRANSFORMS), len(PUBLIC_FEATURES)),
        synergy_weights=normalize_weights(rng.dirichlet(np.ones(len(SYNERGIES)) * 0.8)),
        cd8_mix=float(rng.uniform(0.68, 0.92)),
        full_mode_mix=float(rng.uniform(0.68, 0.95)),
        low_presentation_gate=float(rng.uniform(0.00, 0.32)),
        low_presentation_penalty=float(rng.uniform(0.00, 0.42)),
        high_leakage_penalty=float(rng.uniform(0.10, 0.45)),
        medium_leakage_penalty=float(rng.uniform(0.00, 0.18)),
        complexity_penalty=float(rng.uniform(0.00, 0.05)),
    )


def axis_scores(g: Genome, df: pd.DataFrame) -> dict[str, np.ndarray]:
    cd8_x = df[CD8_FEATURES].to_numpy(dtype=float)
    cd4_x = df[CD4_FEATURES].to_numpy(dtype=float)
    public_x = df[PUBLIC_FEATURES].to_numpy(dtype=float)
    cd8_z = np.zeros_like(cd8_x)
    cd4_z = np.zeros_like(cd4_x)
    public_z = np.zeros_like(public_x)
    for i, t in enumerate(g.cd8_transforms):
        cd8_z[:, i] = transform_values(cd8_x[:, i], TRANSFORMS[int(t)])
    for i, t in enumerate(g.cd4_transforms):
        cd4_z[:, i] = transform_values(cd4_x[:, i], TRANSFORMS[int(t)])
    for i, t in enumerate(g.public_transforms):
        public_z[:, i] = transform_values(public_x[:, i], TRANSFORMS[int(t)])
    cd8 = minmax_array(cd8_z @ cap_cd8_weights(g.cd8_weights))
    cd4 = minmax_array(cd4_z @ normalize_weights(g.cd4_weights))
    public_fallback = minmax_array(public_z @ normalize_weights(g.public_weights))
    presentation = 0.58 * df["bigmhc_el_norm"].to_numpy(dtype=float) + 0.42 * df["bigmhc_im_norm"].to_numpy(dtype=float)
    sy = np.vstack(
        [
            presentation * df["tcr_norm"].to_numpy(dtype=float),
            df["foreignness_norm"].to_numpy(dtype=float) * cd4,
            df["cross_stress_norm"].to_numpy(dtype=float) * df["md_norm"].to_numpy(dtype=float),
            cd8 * cd4,
            df["claimsafe_norm"].to_numpy(dtype=float) * cd4,
        ]
    ).T
    synergy = minmax_array(sy @ normalize_weights(g.synergy_weights))
    full_pan = g.cd8_mix * cd8 + (1.0 - g.cd8_mix) * cd4 + 0.18 * synergy
    low_presentation = presentation < g.low_presentation_gate
    full_pan[low_presentation] *= 1.0 - g.low_presentation_penalty
    public_fallback[low_presentation] *= 1.0 - 0.5 * g.low_presentation_penalty
    if "leakage_risk_level" in df:
        risk = df["leakage_risk_level"].astype(str).str.lower()
        full_pan[risk.eq("high").to_numpy()] *= 1.0 - g.high_leakage_penalty
        full_pan[risk.eq("medium").to_numpy()] *= 1.0 - g.medium_leakage_penalty
    full_pan = minmax_array(full_pan)
    public_fallback = minmax_array(public_fallback)
    available = df.get("internal_feature_available", pd.Series(True, index=df.index)).astype(bool).to_numpy()
    guarded = np.where(available, g.full_mode_mix * full_pan + (1.0 - g.full_mode_mix) * public_fallback, public_fallback)
    return {
        "cd8": minmax_array(cd8),
        "cd4": minmax_array(cd4),
        "synergy": synergy,
        "full_pan": full_pan,
        "public_fallback": public_fallback,
        "pan": minmax_array(guarded),
    }


def som_coverage_score(y: np.ndarray, score: np.ndarray, som_xy: np.ndarray, k: int = 96) -> float:
    idx = np.argsort(score)[::-1][: min(k, len(score))]
    if len(idx) == 0:
        return 0.0
    coords = som_xy[idx]
    occupied = len({tuple(x) for x in coords})
    coverage = occupied / min(k, 64)
    hits = float(np.mean(y[idx])) if len(idx) else 0.0
    return float(0.65 * coverage + 0.35 * hits)


def source_balanced_ap(y: np.ndarray, s: np.ndarray, source_masks: list[np.ndarray]) -> float:
    vals = [safe_ap(y[m], s[m]) for m in source_masks]
    vals = [v for v in vals if np.isfinite(v)]
    return float(np.mean(vals)) if vals else 0.0


def evaluate_genome(
    g: Genome,
    df: pd.DataFrame,
    train_mask: np.ndarray,
    low_mask: np.ndarray,
    source_masks: list[np.ndarray],
    som_xy: np.ndarray,
) -> tuple[float, dict[str, np.ndarray], dict[str, float]]:
    scores = axis_scores(g, df)
    s = scores["pan"]
    y = df["label"].to_numpy(dtype=int)
    train_ap = safe_ap(y[train_mask], s[train_mask])
    train_auc = safe_auc(y[train_mask], s[train_mask])
    low_ap = safe_ap(y[train_mask & low_mask], s[train_mask & low_mask])
    source_ap = source_balanced_ap(y, s, source_masks)
    deployment_s = scores["public_fallback"]
    deployment_ap = safe_ap(y[train_mask], deployment_s[train_mask])
    deployment_low_ap = safe_ap(y[train_mask & low_mask], deployment_s[train_mask & low_mask])
    top96 = topk_precision(y[train_mask], s[train_mask], 96)
    som_score = som_coverage_score(y[train_mask], s[train_mask], som_xy[train_mask], 96)
    base = 0.5 * df["bigmhc_im_norm"].to_numpy(dtype=float) + 0.5 * df["bigmhc_el_norm"].to_numpy(dtype=float)
    disagreement = float(np.mean(np.abs(s[train_mask] - base[train_mask])))
    novelty = 1.0 - min(1.0, abs(disagreement - 0.22) / 0.22)
    cd4_ratio = float(np.mean(scores["cd4"][train_mask]) / max(1e-6, np.mean(scores["cd8"][train_mask])))
    axis_balance = 1.0 - min(1.0, abs(cd4_ratio - 0.75) / 0.75)
    active = float(
        np.sum(cap_cd8_weights(g.cd8_weights) > 0.03)
        + np.sum(normalize_weights(g.cd4_weights) > 0.04)
        + np.sum(normalize_weights(g.public_weights) > 0.04)
    )
    complexity = active / (len(CD8_FEATURES) + len(CD4_FEATURES) + len(PUBLIC_FEATURES)) + g.complexity_penalty
    comps = {
        "train_ap": 0.0 if not np.isfinite(train_ap) else float(train_ap),
        "train_auc": 0.0 if not np.isfinite(train_auc) else float(train_auc),
        "low_medium_ap": 0.0 if not np.isfinite(low_ap) else float(low_ap),
        "source_balanced_ap": source_ap,
        "deployment_public_ap": 0.0 if not np.isfinite(deployment_ap) else float(deployment_ap),
        "deployment_low_medium_ap": 0.0 if not np.isfinite(deployment_low_ap) else float(deployment_low_ap),
        "top96_precision": 0.0 if not np.isfinite(top96) else float(top96),
        "som_niche_score": som_score,
        "novelty_disagreement": novelty,
        "axis_balance": axis_balance,
        "complexity": complexity,
    }
    fitness = (
        0.20 * comps["train_ap"]
        + 0.14 * comps["source_balanced_ap"]
        + 0.13 * comps["low_medium_ap"]
        + 0.13 * comps["deployment_public_ap"]
        + 0.10 * comps["deployment_low_medium_ap"]
        + 0.11 * comps["top96_precision"]
        + 0.07 * comps["som_niche_score"]
        + 0.05 * comps["novelty_disagreement"]
        + 0.04 * comps["axis_balance"]
        + 0.05 * comps["train_auc"]
        - 0.05 * comps["complexity"]
    )
    comps["fitness"] = float(fitness)
    return float(fitness), scores, comps


def mutate(g: Genome, rng: np.random.Generator, op: str, scale: float) -> Genome:
    cd8 = g.cd8_weights.copy()
    cd4 = g.cd4_weights.copy()
    public = g.public_weights.copy()
    cd8_t = g.cd8_transforms.copy()
    cd4_t = g.cd4_transforms.copy()
    public_t = g.public_transforms.copy()
    sy = g.synergy_weights.copy()
    cd8_mix = g.cd8_mix
    full_mode_mix = g.full_mode_mix
    gate = g.low_presentation_gate
    gate_pen = g.low_presentation_penalty
    high_pen = g.high_leakage_penalty
    med_pen = g.medium_leakage_penalty
    cx = g.complexity_penalty
    if op == "cd8_reweight":
        for i in rng.choice(len(cd8), size=max(1, len(cd8) // 4), replace=False):
            cd8[i] *= float(np.exp(rng.normal(0, scale)))
    elif op == "cd4_reweight":
        for i in rng.choice(len(cd4), size=max(1, len(cd4) // 3), replace=False):
            cd4[i] *= float(np.exp(rng.normal(0, scale)))
    elif op == "public_reweight":
        for i in rng.choice(len(public), size=max(1, len(public) // 4), replace=False):
            public[i] *= float(np.exp(rng.normal(0, scale)))
    elif op == "axis_shift":
        cd8_mix = float(np.clip(cd8_mix + rng.normal(0, scale * 0.08), 0.60, 0.95))
        full_mode_mix = float(np.clip(full_mode_mix + rng.normal(0, scale * 0.08), 0.45, 0.98))
    elif op == "synergy_shift":
        for i in rng.choice(len(sy), size=2, replace=False):
            sy[i] *= float(np.exp(rng.normal(0, scale)))
    elif op == "gate_shift":
        gate = float(np.clip(gate + rng.normal(0, scale * 0.06), 0.0, 0.55))
        gate_pen = float(np.clip(gate_pen + rng.normal(0, scale * 0.08), 0.0, 0.65))
        high_pen = float(np.clip(high_pen + rng.normal(0, scale * 0.08), 0.0, 0.70))
        med_pen = float(np.clip(med_pen + rng.normal(0, scale * 0.04), 0.0, 0.35))
    elif op == "transform_flip":
        if rng.random() < 0.65:
            cd8_t[rng.integers(0, len(cd8_t))] = rng.integers(0, len(TRANSFORMS))
        elif rng.random() < 0.82:
            cd4_t[rng.integers(0, len(cd4_t))] = rng.integers(0, len(TRANSFORMS))
        else:
            public_t[rng.integers(0, len(public_t))] = rng.integers(0, len(TRANSFORMS))
    cx = float(np.clip(cx + rng.normal(0, 0.006), 0.0, 0.08))
    return Genome(
        cd8_weights=cap_cd8_weights(cd8),
        cd4_weights=normalize_weights(cd4),
        public_weights=normalize_weights(public),
        cd8_transforms=cd8_t,
        cd4_transforms=cd4_t,
        public_transforms=public_t,
        synergy_weights=normalize_weights(sy),
        cd8_mix=cd8_mix,
        full_mode_mix=full_mode_mix,
        low_presentation_gate=gate,
        low_presentation_penalty=gate_pen,
        high_leakage_penalty=high_pen,
        medium_leakage_penalty=med_pen,
        complexity_penalty=cx,
    )


def crossover(a: Genome, b: Genome, rng: np.random.Generator) -> Genome:
    cd8_mask = rng.random(len(CD8_FEATURES)) < 0.5
    cd4_mask = rng.random(len(CD4_FEATURES)) < 0.5
    public_mask = rng.random(len(PUBLIC_FEATURES)) < 0.5
    sy_mask = rng.random(len(SYNERGIES)) < 0.5
    return Genome(
        cd8_weights=cap_cd8_weights(np.where(cd8_mask, a.cd8_weights, b.cd8_weights)),
        cd4_weights=normalize_weights(np.where(cd4_mask, a.cd4_weights, b.cd4_weights)),
        public_weights=normalize_weights(np.where(public_mask, a.public_weights, b.public_weights)),
        cd8_transforms=np.where(cd8_mask, a.cd8_transforms, b.cd8_transforms),
        cd4_transforms=np.where(cd4_mask, a.cd4_transforms, b.cd4_transforms),
        public_transforms=np.where(public_mask, a.public_transforms, b.public_transforms),
        synergy_weights=normalize_weights(np.where(sy_mask, a.synergy_weights, b.synergy_weights)),
        cd8_mix=float((a.cd8_mix + b.cd8_mix) / 2),
        full_mode_mix=float((a.full_mode_mix + b.full_mode_mix) / 2),
        low_presentation_gate=float((a.low_presentation_gate + b.low_presentation_gate) / 2),
        low_presentation_penalty=float((a.low_presentation_penalty + b.low_presentation_penalty) / 2),
        high_leakage_penalty=float((a.high_leakage_penalty + b.high_leakage_penalty) / 2),
        medium_leakage_penalty=float((a.medium_leakage_penalty + b.medium_leakage_penalty) / 2),
        complexity_penalty=float((a.complexity_penalty + b.complexity_penalty) / 2),
    )


def choose_operator(rng: np.random.Generator, pulls: dict[str, int], rewards: dict[str, float], gen: int) -> str:
    if gen < 5 or rng.random() < 0.12:
        return str(rng.choice(OPERATORS))
    total = sum(pulls.values()) + 1
    scores = {}
    for op in OPERATORS:
        n = pulls.get(op, 0) + 1
        mean = rewards.get(op, 0.0) / n
        scores[op] = mean + 0.08 * math.sqrt(math.log(total + 1) / n)
    return max(scores, key=scores.get)


def run_ga(df: pd.DataFrame, som_xy: np.ndarray, generations: int = 75, population_size: int = 144, seed: int = 29):
    rng = np.random.default_rng(seed)
    train_mask = df["split_biodarwin"].eq("academic_train_sources").to_numpy()
    low_mask = df["low_medium_leakage"].to_numpy()
    sources = df["source_name"].astype(str)
    source_masks = [train_mask & sources.eq(src).to_numpy() for src in sorted(sources[train_mask].unique())]
    pop = [random_genome(rng) for _ in range(population_size)]
    pulls = {op: 0 for op in OPERATORS}
    rewards = {op: 0.0 for op in OPERATORS}
    best_g = pop[0]
    best_scores = axis_scores(best_g, df)
    best_fit = -1e9
    best_comps: dict[str, float] = {}
    trace = []
    op_rows = []
    for gen in range(generations):
        scored = []
        for g in pop:
            fit, scores, comps = evaluate_genome(g, df, train_mask, low_mask, source_masks, som_xy)
            scored.append((fit, g, scores, comps))
        scored.sort(key=lambda x: x[0], reverse=True)
        if scored[0][0] > best_fit:
            best_fit, best_g, best_scores, best_comps = scored[0]
        trace.append(
            {
                "generation": gen,
                "best_fitness": scored[0][0],
                "mean_fitness": float(np.mean([x[0] for x in scored])),
                **scored[0][3],
            }
        )
        elites = [x[1] for x in scored[: max(6, population_size // 12)]]
        next_pop = elites.copy()
        candidate_pool = scored[: max(18, population_size // 2)]
        rate = max(0.18, 0.55 * (1.0 - gen / max(1, generations)))
        while len(next_pop) < population_size:
            pidx = rng.choice(len(candidate_pool), size=2, replace=False)
            parent_fit = max(candidate_pool[pidx[0]][0], candidate_pool[pidx[1]][0])
            child = crossover(candidate_pool[pidx[0]][1], candidate_pool[pidx[1]][1], rng)
            op = choose_operator(rng, pulls, rewards, gen)
            child = mutate(child, rng, op, rate)
            child_fit, _, _ = evaluate_genome(child, df, train_mask, low_mask, source_masks, som_xy)
            pulls[op] += 1
            rewards[op] += max(0.0, child_fit - parent_fit)
            next_pop.append(child)
        op_rows.append({"generation": gen, **{f"pulls_{op}": pulls[op] for op in OPERATORS}, **{f"reward_{op}": rewards[op] for op in OPERATORS}})
        pop = next_pop
    return best_g, best_scores, best_comps, pd.DataFrame(trace), pd.DataFrame(op_rows), pulls, rewards


def metrics_table(df: pd.DataFrame, label_col: str, score_specs: list[tuple[str, str]], split_col: str | None = None) -> pd.DataFrame:
    rows = []
    if split_col:
        splits = {str(k): df[split_col].astype(str).eq(str(k)) for k in sorted(df[split_col].astype(str).unique())}
    else:
        splits = {"all": pd.Series(True, index=df.index)}
    for split, mask in splits.items():
        d = df[mask].copy()
        y = d[label_col].astype(int).to_numpy()
        for name, col in score_specs:
            if col not in d:
                continue
            s = pd.to_numeric(d[col], errors="coerce").to_numpy()
            row = {"split": split, "algorithm": name, "n": int(len(d)), "positives": int(y.sum()), "positive_rate": float(np.mean(y))}
            row["AUPRC"] = safe_ap(y, s)
            row["AUROC"] = safe_auc(y, s)
            for k in [5, 10, 24, 96]:
                val = topk_precision(y, s, k)
                row[f"top{k}_precision"] = val
                row[f"top{k}_hits"] = int(round(val * min(k, len(d)))) if np.isfinite(val) else 0
            rows.append(row)
    return pd.DataFrame(rows)


def genome_tables(g: Genome) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cd8 = pd.DataFrame(
        {
            "axis": "CD8_ClassI",
            "feature": CD8_FEATURES,
            "weight": cap_cd8_weights(g.cd8_weights),
            "transform": [TRANSFORMS[int(i)] for i in g.cd8_transforms],
        }
    ).sort_values("weight", ascending=False)
    cd4 = pd.DataFrame(
        {
            "axis": "CD4_ClassII_helper",
            "feature": CD4_FEATURES,
            "weight": normalize_weights(g.cd4_weights),
            "transform": [TRANSFORMS[int(i)] for i in g.cd4_transforms],
        }
    ).sort_values("weight", ascending=False)
    public = pd.DataFrame(
        {
            "axis": "Deployment_public_fallback",
            "feature": PUBLIC_FEATURES,
            "weight": normalize_weights(g.public_weights),
            "transform": [TRANSFORMS[int(i)] for i in g.public_transforms],
        }
    ).sort_values("weight", ascending=False)
    ctrl = pd.DataFrame(
        [
            {"parameter": "cd8_mix", "value": g.cd8_mix},
            {"parameter": "cd4_mix", "value": 1.0 - g.cd8_mix},
            {"parameter": "full_mode_mix", "value": g.full_mode_mix},
            {"parameter": "public_fallback_mix_in_full_mode", "value": 1.0 - g.full_mode_mix},
            {"parameter": "low_presentation_gate", "value": g.low_presentation_gate},
            {"parameter": "low_presentation_penalty", "value": g.low_presentation_penalty},
            {"parameter": "high_leakage_penalty", "value": g.high_leakage_penalty},
            {"parameter": "medium_leakage_penalty", "value": g.medium_leakage_penalty},
            {"parameter": "complexity_penalty", "value": g.complexity_penalty},
            *[
                {"parameter": f"synergy_{name}", "value": val}
                for name, val in zip(SYNERGIES, normalize_weights(g.synergy_weights), strict=True)
            ],
        ]
    )
    return cd8, pd.concat([cd4, public], ignore_index=True), ctrl


def write_figures(trace: pd.DataFrame, metrics: pd.DataFrame, weights: pd.DataFrame, academic: pd.DataFrame, som_xy: np.ndarray) -> None:
    fig_dir = OUT / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    ax.plot(trace["generation"], trace["best_fitness"], label="best", color="#2f6f73")
    ax.plot(trace["generation"], trace["mean_fitness"], label="mean", color="#c59b3b")
    ax.set_xlabel("Generation")
    ax.set_ylabel("Fitness")
    ax.set_title("BioDarwin GA/RL convergence")
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig1_biodarwin_convergence.png", dpi=190)
    plt.close(fig)

    show = metrics[metrics["split"].isin(["academic_validation_like", "industrial_locked_v0"])].copy()
    show = show[
        show["algorithm"].isin(
            [
                "BigMHC_IM",
                "BigMHC_EL",
                "CROSS_stress",
                "CROSS_claimsafe",
                "GA_public_feature_adapter",
                "BioDarwin_PublicFallback",
                "BioDarwin_PanVax",
            ]
        )
    ]
    if not show.empty:
        fig, ax = plt.subplots(figsize=(10, 5.4))
        labels = show["split"] + " · " + show["algorithm"]
        ax.barh(labels, show["AUPRC"], color="#2f6f73")
        ax.scatter(show["AUROC"], labels, color="#c59b3b", s=52, label="AUROC", zorder=3)
        ax.set_xlim(0, 1)
        ax.set_xlabel("Score")
        ax.set_title("Held-out validation and public industrial OOD")
        ax.legend(loc="lower right")
        fig.tight_layout()
        fig.savefig(fig_dir / "fig2_heldout_industrial_metrics.png", dpi=190)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    wshow = weights.sort_values("weight").tail(16)
    colors = np.where(wshow["axis"].eq("CD8_ClassI"), "#2f6f73", "#8f6bb3")
    ax.barh(wshow["feature"], wshow["weight"], color=colors)
    ax.set_xlabel("Frozen genome weight")
    ax.set_title("BioDarwin feature architecture")
    fig.tight_layout()
    fig.savefig(fig_dir / "fig3_frozen_architecture_weights.png", dpi=190)
    plt.close(fig)

    y = academic["label"].to_numpy(dtype=int)
    s = academic["biodarwin_pan_vaccine_score"].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(7.6, 6.4))
    sc = ax.scatter(som_xy[:, 0], som_xy[:, 1], c=s, s=24 + 30 * y, cmap="viridis", alpha=0.78, edgecolors="none")
    ax.set_xlabel("SOM x")
    ax.set_ylabel("SOM y")
    ax.set_title("SOM niche map: score intensity and known hits")
    fig.colorbar(sc, ax=ax, label="BioDarwin score")
    fig.tight_layout()
    fig.savefig(fig_dir / "fig4_som_niche_map.png", dpi=190)
    plt.close(fig)


def table_html(df: pd.DataFrame, max_rows: int = 30) -> str:
    if df.empty:
        return "<p>No rows.</p>"
    show = df.head(max_rows).copy()
    for col in show.columns:
        if pd.api.types.is_float_dtype(show[col]):
            show[col] = show[col].map(lambda x: "" if pd.isna(x) else f"{x:.4f}")
    return show.to_html(index=False, escape=True, classes="data")


def write_report_and_html(summary: dict[str, object], metrics: pd.DataFrame, weights: pd.DataFrame, ctrl: pd.DataFrame, classii: pd.DataFrame) -> None:
    report = [
        "# BioDarwin Pan-Vaccine GA/RL algorithm",
        "",
        f"Generated: {summary['generated_at']}",
        "",
        "## Claim boundary",
        "",
        "- This is the algorithm artifact, not the public-source scout.",
        "- GA/RL fitness used academic training sources only; validation-like and industrial public locked rows are reported after freezing.",
        "- Class I is scored directly; Class II is represented as a helper axis and scout lane until curated Class-II labels are added.",
        "- Industrial/public patent rows remain do-not-train and manual-QA required.",
        "",
        "## Outputs",
        "",
        f"- Metrics: `{OUT / 'biodarwin_pan_vaccine_metrics.tsv'}`",
        f"- Scores: `{OUT / 'biodarwin_academic_candidate_scores.tsv'}`",
        f"- Industrial scores: `{OUT / 'biodarwin_industrial_locked_scores.tsv'}`",
        f"- Class-II scout: `{OUT / 'biodarwin_classII_scout_scores.tsv'}`",
        "",
        "## Frozen metrics",
        "",
        metrics.to_markdown(index=False),
        "",
    ]
    (OUT / "BIODARWIN_PAN_VACCINE_GA_RL_REPORT.md").write_text("\n".join(report), encoding="utf-8")

    asset_dir = HUB / "assets/biodarwin_pan_vaccine_ga_rl"
    asset_dir.mkdir(parents=True, exist_ok=True)
    for fig in (OUT / "figures").glob("*.png"):
        shutil.copy2(fig, asset_dir / fig.name)

    best_ind = metrics[metrics["split"].eq("industrial_locked_v0")].sort_values("AUPRC", ascending=False).head(1)
    best_val = metrics[metrics["split"].eq("academic_validation_like")].sort_values("AUPRC", ascending=False).head(1)
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>BioDarwin Pan-Vaccine GA/RL</title>
<style>
:root {{ --bg:#101418; --panel:#151c22; --ink:#e9edf0; --muted:#9aa7af; --gold:#c59b3b; --line:#29343c; --teal:#69a7a2; --violet:#8f6bb3; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font:15px/1.55 system-ui,-apple-system,Segoe UI,sans-serif; }}
header {{ padding:42px 5vw 30px; background:#0c1115; border-bottom:1px solid var(--line); }}
.kicker {{ color:var(--gold); text-transform:uppercase; letter-spacing:.12em; font-weight:700; font-size:12px; }}
h1 {{ margin:.3rem 0 .65rem; font-size:clamp(30px,4vw,54px); line-height:1.05; }}
.lead {{ max-width:1040px; color:#d0d8dd; font-size:18px; }}
.stats {{ display:grid; grid-template-columns:repeat(5,minmax(130px,1fr)); gap:10px; margin-top:22px; }}
.stat {{ border:1px solid var(--line); background:var(--panel); padding:13px 14px; border-radius:8px; }}
.stat b {{ display:block; font-size:23px; }}
.stat span {{ color:var(--muted); font-size:12px; }}
main {{ padding:30px 5vw 60px; display:grid; grid-template-columns:270px minmax(0,1fr); gap:28px; }}
nav {{ position:sticky; top:0; align-self:start; max-height:100vh; overflow:auto; padding:14px; border:1px solid var(--line); border-radius:8px; background:#111820; }}
nav a {{ display:block; color:#d6dee2; text-decoration:none; padding:7px 0; border-bottom:1px solid #202b32; }}
h2 {{ border-bottom:1px solid var(--line); padding-bottom:8px; }}
.num {{ color:var(--gold); margin-right:8px; }}
.note {{ border-left:3px solid var(--gold); background:#171f26; padding:10px 14px; color:#dce3e7; }}
table.data {{ width:100%; border-collapse:collapse; font-size:13px; margin:12px 0 28px; }}
table.data th, table.data td {{ border-bottom:1px solid var(--line); padding:7px 8px; text-align:left; vertical-align:top; }}
table.data th {{ color:#f4d891; background:#141b21; }}
.grid {{ display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:18px; }}
img {{ max-width:100%; border:1px solid var(--line); border-radius:8px; background:#fff; }}
code {{ color:#f4d891; }}
@media(max-width:900px) {{ main {{ grid-template-columns:1fr; }} nav {{ position:static; }} .stats {{ grid-template-columns:repeat(2,1fr); }} .grid {{ grid-template-columns:1fr; }} }}
</style>
</head>
<body>
<header>
  <div class="kicker">BioDarwin · pan-vaccine GA/RL algorithm</div>
  <h1>Class I + Class II cancer-vaccine hunting engine</h1>
  <p class="lead">Independent algorithm artifact with CD8/Class-I scoring, CD4/Class-II helper scoring, SOM niche coverage and bandit-controlled GA evolution. Industrial public rows are held out and do-not-train.</p>
  <div class="stats">
    <div class="stat"><b>{summary['academic_n']}</b><span>academic rows</span></div>
    <div class="stat"><b>{summary['industrial_n']}</b><span>industrial locked rows</span></div>
    <div class="stat"><b>{summary['classII_scout_n']}</b><span>Class-II scout rows</span></div>
    <div class="stat"><b>{best_val['AUPRC'].iloc[0] if not best_val.empty else float('nan'):.3f}</b><span>best validation AUPRC</span></div>
    <div class="stat"><b>{best_ind['AUPRC'].iloc[0] if not best_ind.empty else float('nan'):.3f}</b><span>best industrial AUPRC</span></div>
  </div>
</header>
<main>
<nav>
  <a href="#summary">01 Summary</a>
  <a href="#figures">02 Figures</a>
  <a href="#metrics">03 Metrics</a>
  <a href="#architecture">04 Architecture</a>
  <a href="#classii">05 Class-II Scout</a>
  <a href="#guardrails">06 Guardrails</a>
</nav>
<article>
<section id="summary">
  <h2><span class="num">01</span>Summary</h2>
  <p class="note">This is no longer a rough adapter. It is a frozen pan-vaccine controller: CD8 axis + CD4 helper axis + SOM map + GA/RL operator controller + leakage guards.</p>
</section>
<section id="figures">
  <h2><span class="num">02</span>Figures</h2>
  <div class="grid">
    <img src="assets/biodarwin_pan_vaccine_ga_rl/fig1_biodarwin_convergence.png" alt="convergence">
    <img src="assets/biodarwin_pan_vaccine_ga_rl/fig2_heldout_industrial_metrics.png" alt="metrics">
  </div>
  <div class="grid">
    <img src="assets/biodarwin_pan_vaccine_ga_rl/fig3_frozen_architecture_weights.png" alt="weights">
    <img src="assets/biodarwin_pan_vaccine_ga_rl/fig4_som_niche_map.png" alt="som">
  </div>
</section>
<section id="metrics">
  <h2><span class="num">03</span>Metrics</h2>
  {table_html(metrics, 80)}
</section>
<section id="architecture">
  <h2><span class="num">04</span>Frozen Architecture</h2>
  {table_html(weights, 30)}
  {table_html(ctrl, 30)}
</section>
<section id="classii">
  <h2><span class="num">05</span>Class-II Scout Lane</h2>
  {table_html(classii, 40)}
</section>
<section id="guardrails">
  <h2><span class="num">06</span>Guardrails</h2>
  <p>Industrial public data are not used for training. Class-II rows without curated response labels are scored for scout/plate design only. Suspiciously strong results trigger overlap and source-leakage audits before claim escalation.</p>
</section>
</article>
</main>
</body>
</html>
"""
    html_path = HUB / "biodarwin_pan_vaccine_ga_rl.html"
    html_path.write_text(html, encoding="utf-8")
    try:
        live_asset = LIVE_HUB / "assets/biodarwin_pan_vaccine_ga_rl"
        live_asset.mkdir(parents=True, exist_ok=True)
        for fig in asset_dir.glob("*.png"):
            shutil.copy2(fig, live_asset / fig.name)
        shutil.copy2(html_path, LIVE_HUB / html_path.name)
    except Exception as exc:
        print(f"deploy warning: {exc}")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    academic = prepare_academic()
    som_cols = [
        "bigmhc_im_norm",
        "bigmhc_el_norm",
        "tcr_norm",
        "foreignness_norm",
        "md_norm",
        "seq_helper_prior",
        "seq_diversity",
        "hla_A",
        "hla_B",
        "hla_C",
    ]
    som = SimpleSOM(width=8, height=8, seed=17).fit(academic[som_cols].to_numpy(dtype=float), n_iter=1800)
    som_xy = som.transform(academic[som_cols].to_numpy(dtype=float))
    best, best_scores, comps, trace, op_trace, pulls, rewards = run_ga(academic, som_xy)
    academic["biodarwin_cd8_axis_score"] = best_scores["cd8"]
    academic["biodarwin_cd4_helper_score"] = best_scores["cd4"]
    academic["biodarwin_synergy_score"] = best_scores["synergy"]
    academic["biodarwin_full_mode_score"] = best_scores["full_pan"]
    academic["biodarwin_public_fallback_score"] = best_scores["public_fallback"]
    academic["biodarwin_pan_vaccine_score"] = best_scores["pan"]
    academic["som_x"] = som_xy[:, 0]
    academic["som_y"] = som_xy[:, 1]

    industrial = prepare_industrial()
    ind_scores = axis_scores(best, industrial)
    industrial["biodarwin_cd8_axis_score"] = ind_scores["cd8"]
    industrial["biodarwin_cd4_helper_score"] = ind_scores["cd4"]
    industrial["biodarwin_synergy_score"] = ind_scores["synergy"]
    industrial["biodarwin_full_mode_score"] = ind_scores["full_pan"]
    industrial["biodarwin_public_fallback_score"] = ind_scores["public_fallback"]
    industrial["biodarwin_pan_vaccine_score"] = ind_scores["pan"]

    def frozen_score_fn(df: pd.DataFrame) -> dict[str, np.ndarray]:
        return axis_scores(best, df)

    classii = prepare_classii_scout(frozen_score_fn)

    academic_metrics = metrics_table(academic, "label", BASELINE_SCORES, split_col="split_biodarwin")
    industrial_specs = [
        ("BigMHC_IM", "BigMHC_IM"),
        ("BigMHC_EL", "BigMHC_EL"),
        ("BioDarwin_PublicFallback", "biodarwin_public_fallback_score"),
        ("BioDarwin_PanVax", "biodarwin_pan_vaccine_score"),
    ]
    if "ga_public_feature_adapter_score" in industrial:
        industrial_specs.insert(2, ("GA_public_feature_adapter", "ga_public_feature_adapter_score"))
    industrial_metrics = metrics_table(industrial, "label", industrial_specs, split_col="split_biodarwin")
    all_metrics = pd.concat([academic_metrics, industrial_metrics], ignore_index=True)

    cd8_w, cd4_w, ctrl = genome_tables(best)
    weights = pd.concat([cd8_w, cd4_w], ignore_index=True)
    operator_summary = pd.DataFrame(
        [
            {"operator": op, "pulls": pulls[op], "positive_reward": rewards[op], "mean_positive_reward": rewards[op] / max(1, pulls[op])}
            for op in OPERATORS
        ]
    ).sort_values("mean_positive_reward", ascending=False)

    academic.to_csv(OUT / "biodarwin_academic_candidate_scores.tsv", sep="\t", index=False)
    industrial.to_csv(OUT / "biodarwin_industrial_locked_scores.tsv", sep="\t", index=False)
    classii.to_csv(OUT / "biodarwin_classII_scout_scores.tsv", sep="\t", index=False)
    all_metrics.to_csv(OUT / "biodarwin_pan_vaccine_metrics.tsv", sep="\t", index=False)
    trace.to_csv(OUT / "biodarwin_generation_trace.tsv", sep="\t", index=False)
    op_trace.to_csv(OUT / "biodarwin_operator_bandit_trace.tsv", sep="\t", index=False)
    operator_summary.to_csv(OUT / "biodarwin_operator_bandit_summary.tsv", sep="\t", index=False)
    weights.to_csv(OUT / "biodarwin_frozen_feature_weights.tsv", sep="\t", index=False)
    ctrl.to_csv(OUT / "biodarwin_frozen_controller.tsv", sep="\t", index=False)

    write_figures(trace, all_metrics, weights, academic, som_xy)
    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "algorithm": "BioDarwin_PanVax_GA_RL_v2_deployment_guarded",
        "claim_boundary": "academic training only for GA/RL; validation-like and industrial public locked sets are post-freeze evaluations; Class-II scout rows are not supervised claims; public fallback mode is evolved without industrial labels",
        "academic_n": int(len(academic)),
        "industrial_n": int(len(industrial)),
        "classII_scout_n": int(len(classii)),
        "fitness_components": comps,
        "best_validation_like": all_metrics[all_metrics["split"].eq("academic_validation_like")].sort_values("AUPRC", ascending=False).head(1).to_dict("records"),
        "best_industrial": all_metrics[all_metrics["split"].eq("industrial_locked_v0")].sort_values("AUPRC", ascending=False).head(1).to_dict("records"),
        "outputs": {
            "out_dir": str(OUT),
            "metrics": str(OUT / "biodarwin_pan_vaccine_metrics.tsv"),
            "academic_scores": str(OUT / "biodarwin_academic_candidate_scores.tsv"),
            "industrial_scores": str(OUT / "biodarwin_industrial_locked_scores.tsv"),
            "classII_scout_scores": str(OUT / "biodarwin_classII_scout_scores.tsv"),
            "html": str(HUB / "biodarwin_pan_vaccine_ga_rl.html"),
            "live_html": str(LIVE_HUB / "biodarwin_pan_vaccine_ga_rl.html"),
        },
    }
    (OUT / "biodarwin_pan_vaccine_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    write_report_and_html(summary, all_metrics, weights, ctrl, classii)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
