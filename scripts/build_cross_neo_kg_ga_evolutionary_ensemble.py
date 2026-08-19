#!/usr/bin/env python3
"""Build a knowledge-graph-guided genetic ensemble for CROSS-Neo.

This is a retrospective architecture-search scaffold, not a prospective claim.
The knowledge graph stores model concepts, assay roles and guardrails; the
genetic algorithm evolves valid score-controller combinations from those
concepts. Actual 96-well results should be used as the prospective validation
layer after assay readout.
"""

from __future__ import annotations

import argparse
import html
import json
import math
import random
import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Callable

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCORES = (
    ROOT
    / "project/results/cross_neo_immunogenicity_algorithm_compare_2026_05_10"
    / "cross_neo_immunogenicity_comparator_scores.tsv"
)
DEFAULT_WETLAB = (
    ROOT
    / "project/results/cross_neo_wetlab_algorithm_comparison_v6_2026_05_10"
    / "v6_wetlab_algorithm_joined_calls.tsv"
)
DEFAULT_OUT = ROOT / "project/results/cross_neo_kg_ga_evolutionary_ensemble_2026_05_10"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")


@dataclass(frozen=True)
class FeatureSpec:
    name: str
    column: str
    concept: str
    role: str
    layer: str


FEATURES = [
    FeatureSpec("BigMHC_IM", "bigmhc_im_score", "public_immunogenicity_prior", "presentation_immunogenicity", "atomic"),
    FeatureSpec("BigMHC_EL", "bigmhc_el_score", "public_presentation_prior", "hla_presentation", "atomic"),
    FeatureSpec("CROSS_core", "crossneo_core_score_norm", "internal_ranker", "candidate_ranking", "atomic"),
    FeatureSpec("CROSS_stress", "stress_guarded_discovery_score", "stress_guarded_ranker", "robust_candidate_ranking", "atomic"),
    FeatureSpec("CROSS_BMA", "bma_v2_discovery_score", "bayesian_model_average", "ensemble_stability", "atomic"),
    FeatureSpec("CROSS_finetuned", "finetuned_experiment_priority_score", "experiment_priority", "wetlab_priority", "atomic"),
    FeatureSpec("Impact_portfolio", "impact_portfolio_score", "translation_impact", "product_value", "atomic"),
    FeatureSpec("TCR_expert", "tcr_recognition_score_norm", "tcr_recognition", "pmhc_tcr_likelihood", "atomic"),
    FeatureSpec("Foreignness", "fitness_foreignness_proxy", "mutant_foreignness", "self_nonself_gap", "atomic"),
    FeatureSpec("MD_control", "md_control_score_norm", "structure_control_readiness", "claim_safety", "atomic"),
    FeatureSpec("Fixed_integrated", "immunogenicity_discovery_score", "fixed_integrated_score", "legacy_controller", "meta"),
    FeatureSpec("Fixed_claimsafe", "immunogenicity_claim_safe_score", "fixed_claim_safe_score", "legacy_claim_guard", "meta"),
]

BASELINE_SCORES = [
    ("BigMHC_IM", "bigmhc_im_score"),
    ("BigMHC_EL", "bigmhc_el_score"),
    ("CROSS_stress", "stress_guarded_discovery_score"),
    ("CROSS_BMA", "bma_v2_discovery_score"),
    ("CROSS_finetuned", "finetuned_experiment_priority_score"),
    ("CROSS_integrated", "immunogenicity_discovery_score"),
    ("CROSS_claimsafe", "immunogenicity_claim_safe_score"),
    ("KG_GA_evolved", "kg_ga_evolved_score"),
]

TRANSFORMS = ["identity", "sqrt", "square", "sigmoid"]
SYNERGY_PAIRS = [
    ("BigMHC_EL", "TCR_expert", "presentation_x_tcr"),
    ("BigMHC_IM", "Foreignness", "immunogenicity_x_foreignness"),
    ("CROSS_stress", "BigMHC_EL", "stress_x_presentation"),
    ("CROSS_core", "TCR_expert", "ranker_x_tcr"),
    ("MD_control", "Foreignness", "claim_safety_x_foreignness"),
    ("Fixed_claimsafe", "BigMHC_IM", "claimsafe_x_public_im"),
]

FEATURE_WEIGHT_CAPS = {
    "TCR_expert": 0.18,
    "Impact_portfolio": 0.22,
    "Fixed_integrated": 0.20,
    "Fixed_claimsafe": 0.20,
}
DEFAULT_FEATURE_WEIGHT_CAP = 0.30


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t")


def minmax(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    lo = x.min(skipna=True)
    hi = x.max(skipna=True)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi == lo:
        return pd.Series(np.where(x.notna(), 0.5, np.nan), index=s.index)
    return (x - lo) / (hi - lo)


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
        return out
    out[finite] = (x[finite] - lo) / (hi - lo)
    return out


def safe_metric(y: np.ndarray, score: np.ndarray, metric: str) -> float:
    mask = np.isfinite(score) & np.isfinite(y)
    if mask.sum() == 0 or len(np.unique(y[mask])) < 2:
        return np.nan
    if metric == "ap":
        return float(average_precision_score(y[mask], score[mask]))
    if metric == "auc":
        return float(roc_auc_score(y[mask], score[mask]))
    raise ValueError(metric)


def topk_precision(y: np.ndarray, score: np.ndarray, k: int) -> float:
    mask = np.isfinite(score) & np.isfinite(y)
    if mask.sum() == 0:
        return np.nan
    k = min(k, int(mask.sum()))
    idx = np.where(mask)[0]
    order = idx[np.argsort(score[idx])[::-1]][:k]
    return float(np.mean(y[order])) if len(order) else np.nan


def source_balanced_ap(y: np.ndarray, score: np.ndarray, source_masks: list[np.ndarray]) -> float:
    vals = []
    for m in source_masks:
        val = safe_metric(y[m], score[m], "ap")
        if np.isfinite(val):
            vals.append(val)
    return float(np.mean(vals)) if vals else np.nan


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


def normalize_weights(w: np.ndarray) -> np.ndarray:
    w = np.clip(np.asarray(w, dtype=float), 0.0, None)
    if not np.isfinite(w).all() or w.sum() <= 0:
        return np.ones_like(w) / len(w)
    return w / w.sum()


def apply_feature_caps(w: np.ndarray) -> np.ndarray:
    """Normalize weights while preventing a single proxy-like concept takeover."""
    w = normalize_weights(w)
    caps = np.array(
        [FEATURE_WEIGHT_CAPS.get(spec.name, DEFAULT_FEATURE_WEIGHT_CAP) for spec in FEATURES],
        dtype=float,
    )
    for _ in range(20):
        over = w > caps
        if not over.any():
            break
        excess = float((w[over] - caps[over]).sum())
        w[over] = caps[over]
        under = ~over
        room = caps[under] - w[under]
        positive_room = room > 1e-12
        if not positive_room.any() or excess <= 0:
            break
        idx = np.where(under)[0][positive_room]
        room = room[positive_room]
        w[idx] += excess * room / room.sum()
    return normalize_weights(np.minimum(w, caps))


@dataclass
class Genome:
    weights: np.ndarray
    transforms: np.ndarray
    synergy_weights: np.ndarray
    high_leakage_penalty: float
    medium_leakage_penalty: float
    low_presentation_threshold: float
    low_presentation_penalty: float
    claim_blocker_penalty: float

    def copy(self) -> "Genome":
        return Genome(
            self.weights.copy(),
            self.transforms.copy(),
            self.synergy_weights.copy(),
            float(self.high_leakage_penalty),
            float(self.medium_leakage_penalty),
            float(self.low_presentation_threshold),
            float(self.low_presentation_penalty),
            float(self.claim_blocker_penalty),
        )


def random_genome(rng: np.random.Generator) -> Genome:
    weights = apply_feature_caps(rng.dirichlet(np.ones(len(FEATURES)) * 0.8))
    transforms = rng.integers(0, len(TRANSFORMS), size=len(FEATURES))
    synergy_weights = rng.uniform(0.0, 0.20, size=len(SYNERGY_PAIRS))
    return Genome(
        weights=weights,
        transforms=transforms,
        synergy_weights=synergy_weights,
        high_leakage_penalty=float(rng.uniform(0.0, 0.45)),
        medium_leakage_penalty=float(rng.uniform(0.0, 0.20)),
        low_presentation_threshold=float(rng.uniform(0.05, 0.35)),
        low_presentation_penalty=float(rng.uniform(0.0, 0.45)),
        claim_blocker_penalty=float(rng.uniform(0.0, 0.35)),
    )


def prepare_features(df: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray, dict[str, int]]:
    out = df.copy()
    feature_arrays = []
    feature_index = {}
    for i, spec in enumerate(FEATURES):
        if spec.column not in out.columns:
            out[spec.column] = np.nan
        norm_col = f"kg_norm_{spec.name}"
        out[norm_col] = minmax(out[spec.column]).fillna(0.0)
        feature_arrays.append(out[norm_col].to_numpy(dtype=float))
        feature_index[spec.name] = i
    return out, np.vstack(feature_arrays).T, feature_index


def score_genome(
    genome: Genome,
    x: np.ndarray,
    feature_index: dict[str, int],
    leakage: np.ndarray,
    blockers: np.ndarray,
) -> np.ndarray:
    transformed = np.zeros_like(x, dtype=float)
    for i, t_idx in enumerate(genome.transforms):
        transformed[:, i] = transform_values(x[:, i], TRANSFORMS[int(t_idx)])
    weights = apply_feature_caps(genome.weights)
    score = transformed @ weights

    synergy = np.zeros(x.shape[0], dtype=float)
    for j, (left, right, _) in enumerate(SYNERGY_PAIRS):
        li = feature_index[left]
        ri = feature_index[right]
        synergy += genome.synergy_weights[j] * transformed[:, li] * transformed[:, ri]
    if np.sum(genome.synergy_weights) > 0:
        score = 0.85 * score + 0.15 * minmax_array(synergy)

    penalty = np.ones_like(score)
    penalty[leakage == "high"] *= 1.0 - genome.high_leakage_penalty
    penalty[leakage == "medium"] *= 1.0 - genome.medium_leakage_penalty
    el = x[:, feature_index["BigMHC_EL"]]
    penalty[el < genome.low_presentation_threshold] *= 1.0 - genome.low_presentation_penalty
    penalty[blockers] *= 1.0 - genome.claim_blocker_penalty
    score = score * np.clip(penalty, 0.05, 1.0)
    return minmax_array(score)


def fitness_components(
    df: pd.DataFrame,
    y: np.ndarray,
    score: np.ndarray,
    train_mask: np.ndarray,
    lowrisk_mask: np.ndarray,
    x: np.ndarray,
    feature_index: dict[str, int],
    train_source_masks: list[np.ndarray],
) -> dict[str, float]:
    train_ap = safe_metric(y[train_mask], score[train_mask], "ap")
    train_auc = safe_metric(y[train_mask], score[train_mask], "auc")
    source_ap = source_balanced_ap(y, score, train_source_masks)
    low_mask = train_mask & lowrisk_mask
    low_ap = safe_metric(y[low_mask], score[low_mask], "ap")
    top96 = topk_precision(y[train_mask], score[train_mask], 96)
    corr_penalty = 0.0
    refs = ["BigMHC_EL", "CROSS_stress", "Fixed_claimsafe"]
    corrs = []
    for ref in refs:
        r = x[:, feature_index[ref]]
        if np.std(r[train_mask]) > 0 and np.std(score[train_mask]) > 0:
            corrs.append(abs(float(np.corrcoef(score[train_mask], r[train_mask])[0, 1])))
    if corrs:
        corr_penalty = max(corrs)
    diversity = 1.0 - corr_penalty
    pieces = {
        "train_ap": train_ap,
        "train_auc": train_auc,
        "source_balanced_ap": source_ap,
        "lowrisk_ap": low_ap,
        "top96_precision": top96,
        "diversity": diversity,
    }
    fill = {k: (0.0 if not np.isfinite(v) else float(v)) for k, v in pieces.items()}
    fill["fitness"] = (
        0.30 * fill["train_ap"]
        + 0.25 * fill["source_balanced_ap"]
        + 0.15 * fill["lowrisk_ap"]
        + 0.15 * fill["top96_precision"]
        + 0.10 * fill["train_auc"]
        + 0.05 * fill["diversity"]
    )
    return fill


def mutate(genome: Genome, rng: np.random.Generator, rate: float) -> Genome:
    g = genome.copy()
    for i in range(len(g.weights)):
        if rng.random() < rate:
            g.weights[i] *= float(np.exp(rng.normal(0.0, 0.45)))
        if rng.random() < rate * 0.35:
            g.transforms[i] = rng.integers(0, len(TRANSFORMS))
    g.weights = apply_feature_caps(g.weights)
    for j in range(len(g.synergy_weights)):
        if rng.random() < rate:
            g.synergy_weights[j] = float(np.clip(g.synergy_weights[j] + rng.normal(0.0, 0.04), 0.0, 0.35))
    for attr, lo, hi, sigma in [
        ("high_leakage_penalty", 0.0, 0.65, 0.06),
        ("medium_leakage_penalty", 0.0, 0.35, 0.04),
        ("low_presentation_threshold", 0.0, 0.50, 0.04),
        ("low_presentation_penalty", 0.0, 0.65, 0.06),
        ("claim_blocker_penalty", 0.0, 0.55, 0.05),
    ]:
        if rng.random() < rate:
            setattr(g, attr, float(np.clip(getattr(g, attr) + rng.normal(0.0, sigma), lo, hi)))
    return g


def crossover(a: Genome, b: Genome, rng: np.random.Generator) -> Genome:
    mask = rng.random(len(a.weights)) < 0.5
    weights = np.where(mask, a.weights, b.weights)
    weights = apply_feature_caps(weights + rng.normal(0.0, 0.01, size=len(weights)))
    transforms = np.where(mask, a.transforms, b.transforms)
    smask = rng.random(len(a.synergy_weights)) < 0.5
    synergy = np.where(smask, a.synergy_weights, b.synergy_weights)
    return Genome(
        weights=weights,
        transforms=transforms,
        synergy_weights=np.clip(synergy, 0.0, 0.35),
        high_leakage_penalty=float((a.high_leakage_penalty + b.high_leakage_penalty) / 2),
        medium_leakage_penalty=float((a.medium_leakage_penalty + b.medium_leakage_penalty) / 2),
        low_presentation_threshold=float((a.low_presentation_threshold + b.low_presentation_threshold) / 2),
        low_presentation_penalty=float((a.low_presentation_penalty + b.low_presentation_penalty) / 2),
        claim_blocker_penalty=float((a.claim_blocker_penalty + b.claim_blocker_penalty) / 2),
    )


def run_ga(
    df: pd.DataFrame,
    x: np.ndarray,
    feature_index: dict[str, int],
    generations: int,
    population_size: int,
    seed: int,
) -> tuple[Genome, pd.DataFrame, np.ndarray, dict[str, float]]:
    rng = np.random.default_rng(seed)
    random.seed(seed)
    y = pd.to_numeric(df["label"], errors="coerce").fillna(0).to_numpy(dtype=float)
    source = df["source_name"].astype(str)
    frozen_mask = source.str.contains("validation|_Val", case=False, regex=True).to_numpy()
    train_mask = ~frozen_mask
    train_source_masks = [
        train_mask & source.eq(src).to_numpy()
        for src in sorted(source[train_mask].unique())
    ]
    lowrisk_mask = df["leakage_risk_level"].astype(str).str.lower().isin(["low", "medium"]).to_numpy()
    leakage = df["leakage_risk_level"].astype(str).str.lower().to_numpy()
    blockers = df.get("immunogenicity_claim_blockers", "none").astype(str).ne("none").to_numpy()

    population = [random_genome(rng) for _ in range(population_size)]
    trace = []
    best_genome = population[0]
    best_score = None
    best_components: dict[str, float] = {}

    def evaluate(g: Genome) -> tuple[float, np.ndarray, dict[str, float]]:
        scores = score_genome(g, x, feature_index, leakage, blockers)
        comps = fitness_components(df, y, scores, train_mask, lowrisk_mask, x, feature_index, train_source_masks)
        return comps["fitness"], scores, comps

    scored: list[tuple[float, Genome, np.ndarray, dict[str, float]]] = []
    for gen in range(generations):
        scored = []
        for g in population:
            fit, scores, comps = evaluate(g)
            scored.append((fit, g, scores, comps))
        scored.sort(key=lambda z: z[0], reverse=True)
        if best_score is None or scored[0][0] > best_score:
            best_score = scored[0][0]
            best_genome = scored[0][1].copy()
            best_components = scored[0][3]
        trace.append(
            {
                "generation": gen,
                "best_fitness": float(scored[0][0]),
                "mean_fitness": float(np.mean([z[0] for z in scored])),
                **{f"best_{k}": v for k, v in scored[0][3].items()},
            }
        )

        elite_n = max(4, population_size // 12)
        next_population = [z[1].copy() for z in scored[:elite_n]]
        while len(next_population) < population_size:
            pool = random.sample(scored[: max(12, population_size // 2)], k=4)
            parent_a = max(pool[:2], key=lambda z: z[0])[1]
            parent_b = max(pool[2:], key=lambda z: z[0])[1]
            child = crossover(parent_a, parent_b, rng)
            rate = max(0.04, 0.20 * (1.0 - gen / max(generations, 1)))
            child = mutate(child, rng, rate)
            next_population.append(child)
        population = next_population

    final_score = score_genome(best_genome, x, feature_index, leakage, blockers)
    return best_genome, pd.DataFrame(trace), final_score, best_components


def metric_row(df: pd.DataFrame, label_col: str, score_col: str, name: str, split_name: str, mask: pd.Series) -> dict[str, object]:
    d = df.loc[mask, [label_col, score_col]].copy()
    d[label_col] = pd.to_numeric(d[label_col], errors="coerce")
    d[score_col] = pd.to_numeric(d[score_col], errors="coerce")
    d = d.dropna()
    row = {
        "split": split_name,
        "algorithm": name,
        "score_column": score_col,
        "n": int(len(d)),
        "positives": int(d[label_col].sum()) if len(d) else 0,
        "positive_rate": float(d[label_col].mean()) if len(d) else np.nan,
        "AUPRC": np.nan,
        "AUROC": np.nan,
    }
    if len(d) and d[label_col].nunique() == 2:
        row["AUPRC"] = float(average_precision_score(d[label_col], d[score_col]))
        row["AUROC"] = float(roc_auc_score(d[label_col], d[score_col]))
    ordered = d.sort_values(score_col, ascending=False)
    for k in [10, 24, 48, 96]:
        top = ordered.head(k)
        row[f"top{k}_hits"] = int(top[label_col].sum()) if len(top) else 0
        row[f"top{k}_precision"] = float(top[label_col].mean()) if len(top) else np.nan
    return row


def build_benchmark(df: pd.DataFrame) -> pd.DataFrame:
    validation_mask = df["source_name"].astype(str).str.contains("validation|_Val", case=False, regex=True)
    split_masks = {
        "all": pd.Series(True, index=df.index),
        "ga_train_nonvalidation_sources": ~validation_mask,
        "frozen_validation_like_sources": validation_mask,
        "low_medium_leakage": df["leakage_risk_level"].astype(str).str.lower().isin(["low", "medium"]),
    }
    rows = []
    for split_name, mask in split_masks.items():
        for name, col in BASELINE_SCORES:
            if col in df.columns:
                rows.append(metric_row(df, "label", col, name, split_name, mask))
    return pd.DataFrame(rows)


def build_source_benchmark(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for src in sorted(df["source_name"].astype(str).unique()):
        mask = df["source_name"].astype(str).eq(src)
        for name, col in BASELINE_SCORES:
            if col in df.columns:
                rows.append(metric_row(df, "label", col, name, src, mask))
    return pd.DataFrame(rows)


def build_wetlab_smoke(wetlab_path: Path, scores: pd.DataFrame) -> pd.DataFrame:
    if not wetlab_path.exists():
        return pd.DataFrame()
    wet = pd.read_csv(wetlab_path, sep="\t")
    cols = ["candidate_id", "wetlab_hit_binary", "wetlab_countable"]
    cols = [c for c in cols if c in wet.columns]
    if len(cols) < 3:
        return pd.DataFrame()
    merged = wet[cols].merge(
        scores[["candidate_id"] + [c for _, c in BASELINE_SCORES if c in scores.columns]],
        on="candidate_id",
        how="left",
    )
    merged["wetlab_countable"] = merged["wetlab_countable"].astype(str).isin(["True", "true", "1", "PASS", "FAIL"]) | (
        merged["wetlab_countable"] is True
    )
    rows = []
    mask = merged["wetlab_countable"].astype(bool)
    for name, col in BASELINE_SCORES:
        if col in merged.columns:
            rows.append(metric_row(merged, "wetlab_hit_binary", col, name, "v6_smoke_confirmatory", mask))
    return pd.DataFrame(rows)


def genome_to_tables(genome: Genome) -> tuple[pd.DataFrame, pd.DataFrame]:
    weights = apply_feature_caps(genome.weights)
    rows = []
    for i, spec in enumerate(FEATURES):
        rows.append(
            {
                "feature": spec.name,
                "source_column": spec.column,
                "concept": spec.concept,
                "role": spec.role,
                "layer": spec.layer,
                "weight": float(weights[i]),
                "transform": TRANSFORMS[int(genome.transforms[i])],
            }
        )
    feature_table = pd.DataFrame(rows).sort_values("weight", ascending=False)
    synergy_rows = []
    for j, (left, right, name) in enumerate(SYNERGY_PAIRS):
        synergy_rows.append({"synergy": name, "left": left, "right": right, "weight": float(genome.synergy_weights[j])})
    gates = [
        {"gate": "high_leakage_penalty", "value": genome.high_leakage_penalty},
        {"gate": "medium_leakage_penalty", "value": genome.medium_leakage_penalty},
        {"gate": "low_presentation_threshold", "value": genome.low_presentation_threshold},
        {"gate": "low_presentation_penalty", "value": genome.low_presentation_penalty},
        {"gate": "claim_blocker_penalty", "value": genome.claim_blocker_penalty},
        {"gate": "single_concept_default_weight_cap", "value": DEFAULT_FEATURE_WEIGHT_CAP},
    ]
    for feature, cap in FEATURE_WEIGHT_CAPS.items():
        gates.append({"gate": f"single_concept_weight_cap:{feature}", "value": cap})
    gates_table = pd.DataFrame(gates + synergy_rows)
    return feature_table, gates_table


def build_knowledge_graph(genome: Genome, out_dir: Path) -> tuple[nx.DiGraph, pd.DataFrame, pd.DataFrame]:
    g = nx.DiGraph()
    g.add_node("CROSS-Neo_KG_GA", type="architecture", layer="controller")
    for obj in ["known_label_hit_yield", "source_balanced_robustness", "low_leakage_generalization", "wetlab_hit_fail_interpretability"]:
        g.add_node(obj, type="objective", layer="fitness")
        g.add_edge("CROSS-Neo_KG_GA", obj, relation="optimizes")

    for spec in FEATURES:
        g.add_node(spec.concept, type="concept", layer=spec.layer, role=spec.role)
        g.add_node(spec.name, type="feature", column=spec.column, layer=spec.layer)
        g.add_edge(spec.concept, spec.name, relation="materialized_as")
        g.add_edge(spec.name, "CROSS-Neo_KG_GA", relation="candidate_input")

    for transform in TRANSFORMS:
        g.add_node(f"transform:{transform}", type="transform", layer="architecture_gene")
        g.add_edge(f"transform:{transform}", "CROSS-Neo_KG_GA", relation="allowed_gene")

    for left, right, name in SYNERGY_PAIRS:
        g.add_node(name, type="synergy", layer="architecture_gene")
        g.add_edge(left, name, relation="synergy_left")
        g.add_edge(right, name, relation="synergy_right")
        g.add_edge(name, "CROSS-Neo_KG_GA", relation="candidate_input")

    for gate in [
        "high_leakage_penalty",
        "medium_leakage_penalty",
        "low_presentation_gate",
        "claim_blocker_penalty",
        "single_concept_dominance_cap",
        "sparse_proxy_guard",
    ]:
        g.add_node(gate, type="guardrail", layer="claim_safety")
        g.add_edge(gate, "CROSS-Neo_KG_GA", relation="caps_or_penalizes")

    assay_arms = [
        "A_clean_discovery",
        "B_mechanism_TCR_MD",
        "C_label_rescue",
        "D_specificity_moat",
        "E_positive_QC_control",
        "F_model_boundary",
    ]
    for arm in assay_arms:
        g.add_node(arm, type="assay_arm", layer="wetlab")
        g.add_edge("CROSS-Neo_KG_GA", arm, relation="prioritizes_or_interprets")

    feature_table, gates_table = genome_to_tables(genome)
    for _, row in feature_table.iterrows():
        g.nodes[row["feature"]]["evolved_weight"] = float(row["weight"])
        g.nodes[row["feature"]]["evolved_transform"] = row["transform"]

    nodes = pd.DataFrame([{"node": n, **attrs} for n, attrs in g.nodes(data=True)])
    edges = pd.DataFrame([{"source": u, "target": v, **attrs} for u, v, attrs in g.edges(data=True)])
    nodes.to_csv(out_dir / "knowledge_graph_nodes.tsv", sep="\t", index=False)
    edges.to_csv(out_dir / "knowledge_graph_edges.tsv", sep="\t", index=False)
    nx.write_graphml(g, out_dir / "cross_neo_kg_ga.graphml")
    return g, nodes, edges


def plot_outputs(out_dir: Path, graph: nx.DiGraph, trace: pd.DataFrame, benchmark: pd.DataFrame, scores: pd.DataFrame) -> None:
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    ax.plot(trace["generation"], trace["best_fitness"], color="#2f6f73", label="best fitness")
    ax.plot(trace["generation"], trace["mean_fitness"], color="#c59b3b", label="mean fitness")
    ax.set_xlabel("Generation")
    ax.set_ylabel("Fitness")
    ax.set_title("KG-guided GA convergence")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(fig_dir / "fig1_ga_convergence.png", dpi=180)
    plt.close(fig)

    all_metrics = benchmark[benchmark["split"].eq("all")].copy()
    all_metrics = all_metrics.sort_values("AUPRC", ascending=True)
    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    ax.barh(all_metrics["algorithm"], all_metrics["AUPRC"], color="#2f6f73", label="AUPRC")
    ax.scatter(all_metrics["AUROC"], all_metrics["algorithm"], color="#c59b3b", s=60, label="AUROC", zorder=3)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Score")
    ax.set_title("Retrospective benchmark: evolved ensemble vs comparators")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(fig_dir / "fig2_evolved_vs_baselines.png", dpi=180)
    plt.close(fig)

    val = benchmark[benchmark["split"].eq("frozen_validation_like_sources")].sort_values("AUPRC", ascending=True)
    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    ax.barh(val["algorithm"], val["AUPRC"], color="#6f7f95")
    ax.set_xlim(0, 1)
    ax.set_xlabel("AUPRC")
    ax.set_title("Validation-like source check")
    fig.tight_layout()
    fig.savefig(fig_dir / "fig3_validation_like_auprc.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 7))
    pos = nx.spring_layout(graph, seed=11, k=0.65)
    node_colors = []
    for _, attrs in graph.nodes(data=True):
        typ = attrs.get("type", "")
        if typ == "feature":
            node_colors.append("#69a7a2")
        elif typ == "objective":
            node_colors.append("#c59b3b")
        elif typ == "guardrail":
            node_colors.append("#d16b5f")
        elif typ == "assay_arm":
            node_colors.append("#8f9cc2")
        else:
            node_colors.append("#d7dde2")
    nx.draw_networkx_nodes(graph, pos, node_size=320, node_color=node_colors, linewidths=0.6, edgecolors="#27323a", ax=ax)
    nx.draw_networkx_edges(graph, pos, arrows=False, width=0.7, alpha=0.35, edge_color="#77838c", ax=ax)
    labels = {n: n.replace("transform:", "t:") for n in graph.nodes}
    nx.draw_networkx_labels(graph, pos, labels=labels, font_size=6, ax=ax)
    ax.set_title("CROSS-Neo knowledge graph for evolutionary ensemble search")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(fig_dir / "fig4_knowledge_graph.png", dpi=220)
    plt.close(fig)

    top = scores.sort_values("kg_ga_evolved_score", ascending=False).head(96).copy()
    fig, ax = plt.subplots(figsize=(10, 4.8))
    colors = top["label"].map({1: "#2f6f73", 0: "#d16b5f"}).fillna("#777")
    ax.bar(np.arange(len(top)), top["kg_ga_evolved_score"], color=colors)
    ax.set_xlabel("Top 96 candidates by evolved score")
    ax.set_ylabel("KG-GA evolved score")
    ax.set_title("Top-96 known-label recovery")
    ax.set_xticks([])
    fig.tight_layout()
    fig.savefig(fig_dir / "fig5_top96_evolved_score.png", dpi=180)
    plt.close(fig)


def table_html(df: pd.DataFrame, max_rows: int = 30) -> str:
    if df.empty:
        return "<p>No rows.</p>"
    show = df.head(max_rows).copy()
    for col in show.columns:
        if pd.api.types.is_float_dtype(show[col]):
            show[col] = show[col].map(lambda x: "" if pd.isna(x) else f"{x:.3f}")
    return show.to_html(index=False, escape=True, classes="data")


def write_report(
    out_dir: Path,
    benchmark: pd.DataFrame,
    source_benchmark: pd.DataFrame,
    wetlab: pd.DataFrame,
    feature_table: pd.DataFrame,
    gates_table: pd.DataFrame,
) -> None:
    all_rank = benchmark[benchmark["split"].eq("all")].sort_values("AUPRC", ascending=False)
    val_rank = benchmark[benchmark["split"].eq("frozen_validation_like_sources")].sort_values("AUPRC", ascending=False)
    best_all = all_rank.iloc[0].to_dict() if len(all_rank) else {}
    best_val = val_rank.iloc[0].to_dict() if len(val_rank) else {}
    lines = [
        "# CROSS-Neo KG-GA evolutionary ensemble",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Bottom line",
        "",
        "- This run builds a concept knowledge graph and evolves a score/controller architecture from the graph.",
        "- It is retrospective architecture search. The v6 96-well assay remains the prospective validation layer.",
    ]
    if best_all:
        lines.append(
            f"- Best all-candidate retrospective AUPRC: `{best_all['algorithm']}` "
            f"(AUPRC {best_all['AUPRC']:.3f}, AUROC {best_all['AUROC']:.3f})."
        )
    if best_val:
        lines.append(
            f"- Best validation-like source AUPRC: `{best_val['algorithm']}` "
            f"(AUPRC {best_val['AUPRC']:.3f}, AUROC {best_val['AUROC']:.3f})."
        )
    lines.extend(
        [
            "",
            "## Evolved feature weights",
            "",
            feature_table.to_markdown(index=False),
            "",
            "## Evolved gates and synergies",
            "",
            gates_table.to_markdown(index=False),
            "",
            "## Benchmark",
            "",
            benchmark.to_markdown(index=False),
            "",
            "## Source-level benchmark",
            "",
            source_benchmark.to_markdown(index=False),
        ]
    )
    if not wetlab.empty:
        lines.extend(["", "## v6 smoke wetlab comparison", "", wetlab.to_markdown(index=False)])
    lines.extend(
        [
            "",
            "## Claim boundary",
            "",
            "Do not claim prospective superiority from this search alone. The claim-safe wording is: a KG-guided evolutionary controller generated a candidate architecture that is benchmarked retrospectively and is ready for prospective v6 assay validation.",
            "",
        ]
    )
    (out_dir / "KG_GA_EVOLUTIONARY_ENSEMBLE_REPORT.md").write_text("\n".join(lines), encoding="utf-8")


def write_html(
    out_dir: Path,
    benchmark: pd.DataFrame,
    source_benchmark: pd.DataFrame,
    wetlab: pd.DataFrame,
    feature_table: pd.DataFrame,
    gates_table: pd.DataFrame,
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
) -> dict[str, object]:
    asset_dir = HUB / "assets/cross_neo_kg_ga_evolutionary_ensemble"
    asset_dir.mkdir(parents=True, exist_ok=True)
    for fig in (out_dir / "figures").glob("*.png"):
        shutil.copy2(fig, asset_dir / fig.name)

    all_rank = benchmark[benchmark["split"].eq("all")].sort_values("AUPRC", ascending=False)
    val_rank = benchmark[benchmark["split"].eq("frozen_validation_like_sources")].sort_values("AUPRC", ascending=False)
    kg_all = all_rank[all_rank["algorithm"].eq("KG_GA_evolved")]
    kg_val = val_rank[val_rank["algorithm"].eq("KG_GA_evolved")]
    kg_all_ap = float(kg_all.iloc[0]["AUPRC"]) if len(kg_all) else float("nan")
    kg_val_ap = float(kg_val.iloc[0]["AUPRC"]) if len(kg_val) else float("nan")
    top_features = feature_table.head(8)
    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>CROSS-Neo KG-GA Evolutionary Ensemble</title>
  <style>
    :root {{ --bg:#101418; --panel:#151c22; --ink:#e9edf0; --muted:#9aa7af; --gold:#c59b3b; --line:#29343c; --teal:#69a7a2; --red:#d16b5f; }}
    body {{ margin:0; background:var(--bg); color:var(--ink); font:15px/1.55 system-ui, -apple-system, Segoe UI, sans-serif; }}
    header {{ padding:42px 5vw 30px; background:#0c1115; border-bottom:1px solid var(--line); }}
    .kicker {{ color:var(--gold); text-transform:uppercase; letter-spacing:.12em; font-weight:700; font-size:12px; }}
    h1 {{ margin:.3rem 0 .65rem; font-size:clamp(30px,4vw,54px); line-height:1.05; }}
    .lead {{ color:#d0d8dd; max-width:1040px; font-size:18px; }}
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
  <div class="kicker">CROSS-Neo · knowledge graph · genetic architecture search</div>
  <h1>KG-guided evolutionary immunogenicity controller</h1>
  <p class="lead">Model concepts, public predictors, internal rankers, TCR/foreignness/MD signals and claim-safety guardrails are represented as a knowledge graph. A genetic algorithm evolves a valid scoring architecture from that graph, then reports retrospective and validation-like checks.</p>
  <div class="stats">
    <div class="stat"><b>{len(nodes)}</b><span>KG nodes</span></div>
    <div class="stat"><b>{len(edges)}</b><span>KG edges</span></div>
    <div class="stat"><b>{kg_all_ap:.3f}</b><span>KG-GA all AUPRC</span></div>
    <div class="stat"><b>{kg_val_ap:.3f}</b><span>KG-GA validation-like AUPRC</span></div>
    <div class="stat"><b>{html.escape(str(top_features.iloc[0]['feature']))}</b><span>top evolved feature</span></div>
  </div>
</header>
<main>
<nav>
  <a href="#summary">01 Summary</a>
  <a href="#kg">02 Knowledge Graph</a>
  <a href="#architecture">03 Evolved Architecture</a>
  <a href="#benchmark">04 Benchmark</a>
  <a href="#source">05 Source Checks</a>
  <a href="#wetlab">06 v6 Assay Handoff</a>
</nav>
<article>
<section id="summary">
  <h2><span class="num">01</span>Summary</h2>
  <p class="note">This page is an architecture-search dossier, not a prospective superiority claim. The strong claim unlocks only after the v6 96-well assay is interpreted with the same locked scoring architecture.</p>
  <div class="grid">
    <img src="assets/cross_neo_kg_ga_evolutionary_ensemble/fig1_ga_convergence.png" alt="GA convergence">
    <img src="assets/cross_neo_kg_ga_evolutionary_ensemble/fig2_evolved_vs_baselines.png" alt="evolved vs baselines">
  </div>
  <div class="grid">
    <img src="assets/cross_neo_kg_ga_evolutionary_ensemble/fig3_validation_like_auprc.png" alt="validation-like AUPRC">
    <img src="assets/cross_neo_kg_ga_evolutionary_ensemble/fig5_top96_evolved_score.png" alt="top 96 evolved score">
  </div>
</section>
<section id="kg">
  <h2><span class="num">02</span>Knowledge Graph</h2>
  <img src="assets/cross_neo_kg_ga_evolutionary_ensemble/fig4_knowledge_graph.png" alt="knowledge graph">
  {table_html(nodes.sort_values(["type", "node"]), 35)}
</section>
<section id="architecture">
  <h2><span class="num">03</span>Evolved Architecture</h2>
  {table_html(feature_table, 20)}
  {table_html(gates_table, 20)}
</section>
<section id="benchmark">
  <h2><span class="num">04</span>Benchmark</h2>
  {table_html(benchmark.sort_values(["split", "AUPRC"], ascending=[True, False]), 60)}
</section>
<section id="source">
  <h2><span class="num">05</span>Source Checks</h2>
  {table_html(source_benchmark.sort_values(["split", "AUPRC"], ascending=[True, False]), 70)}
</section>
<section id="wetlab">
  <h2><span class="num">06</span>v6 Assay Handoff</h2>
  <p>The same locked <code>kg_ga_evolved_score</code> is merged into the 96-well interpreter output and benchmarked against BigMHC and CROSS-Neo comparators.</p>
  {table_html(wetlab.sort_values("AUPRC", ascending=False) if not wetlab.empty else wetlab, 20)}
</section>
</article>
</main>
</body>
</html>
"""
    html_path = HUB / "cross_neo_kg_ga_evolutionary_ensemble.html"
    html_path.write_text(html_text, encoding="utf-8")

    live_ok = False
    warnings: list[str] = []
    try:
        live_asset_dir = LIVE_HUB / "assets/cross_neo_kg_ga_evolutionary_ensemble"
        live_asset_dir.mkdir(parents=True, exist_ok=True)
        for fig in asset_dir.glob("*.png"):
            shutil.copy2(fig, live_asset_dir / fig.name)
        shutil.copy2(html_path, LIVE_HUB / html_path.name)
        live_ok = True
    except Exception as exc:
        warnings.append(str(exc))
    return {
        "html_path": str(html_path),
        "live_html_path": str(LIVE_HUB / html_path.name),
        "live_deploy_ok": live_ok,
        "live_deploy_warnings": warnings,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, default=DEFAULT_SCORES)
    parser.add_argument("--wetlab", type=Path, default=DEFAULT_WETLAB)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--generations", type=int, default=45)
    parser.add_argument("--population-size", type=int, default=96)
    parser.add_argument("--seed", type=int, default=11)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    scores = read_tsv(args.scores)
    scores, x, feature_index = prepare_features(scores)
    genome, trace, evolved, fitness = run_ga(
        scores,
        x,
        feature_index,
        generations=args.generations,
        population_size=args.population_size,
        seed=args.seed,
    )
    scores["kg_ga_evolved_score"] = evolved
    scores["kg_ga_evolved_rank"] = scores["kg_ga_evolved_score"].rank(method="first", ascending=False).astype(int)

    feature_table, gates_table = genome_to_tables(genome)
    graph, nodes, edges = build_knowledge_graph(genome, args.out_dir)
    benchmark = build_benchmark(scores)
    source_benchmark = build_source_benchmark(scores)
    wetlab = build_wetlab_smoke(args.wetlab, scores)

    scores.to_csv(args.out_dir / "kg_ga_evolved_candidate_scores.tsv", sep="\t", index=False)
    feature_table.to_csv(args.out_dir / "kg_ga_best_feature_weights.tsv", sep="\t", index=False)
    gates_table.to_csv(args.out_dir / "kg_ga_best_gates_and_synergies.tsv", sep="\t", index=False)
    trace.to_csv(args.out_dir / "kg_ga_generation_trace.tsv", sep="\t", index=False)
    benchmark.to_csv(args.out_dir / "kg_ga_benchmark.tsv", sep="\t", index=False)
    source_benchmark.to_csv(args.out_dir / "kg_ga_source_benchmark.tsv", sep="\t", index=False)
    wetlab.to_csv(args.out_dir / "kg_ga_v6_smoke_wetlab_comparison.tsv", sep="\t", index=False)

    architecture = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "seed": args.seed,
        "generations": args.generations,
        "population_size": args.population_size,
        "fitness_components": fitness,
        "weights": feature_table.to_dict(orient="records"),
        "gates_and_synergies": gates_table.to_dict(orient="records"),
        "claim_boundary": "retrospective KG-guided architecture search; prospective superiority requires locked v6 wetlab validation",
    }
    (args.out_dir / "kg_ga_best_architecture.json").write_text(json.dumps(architecture, indent=2), encoding="utf-8")

    plot_outputs(args.out_dir, graph, trace, benchmark, scores)
    write_report(args.out_dir, benchmark, source_benchmark, wetlab, feature_table, gates_table)
    html_info = write_html(args.out_dir, benchmark, source_benchmark, wetlab, feature_table, gates_table, nodes, edges)

    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "input_scores": str(args.scores),
        "out_dir": str(args.out_dir),
        "n_candidates": int(len(scores)),
        "n_kg_nodes": int(len(nodes)),
        "n_kg_edges": int(len(edges)),
        "ga_generations": int(args.generations),
        "ga_population_size": int(args.population_size),
        "best_fitness_components": fitness,
        "best_all_metrics": benchmark[
            benchmark["split"].eq("all") & benchmark["algorithm"].eq("KG_GA_evolved")
        ].to_dict(orient="records"),
        "best_validation_like_metrics": benchmark[
            benchmark["split"].eq("frozen_validation_like_sources") & benchmark["algorithm"].eq("KG_GA_evolved")
        ].to_dict(orient="records"),
        **html_info,
        "output_files": [
            str(args.out_dir / "kg_ga_evolved_candidate_scores.tsv"),
            str(args.out_dir / "kg_ga_best_feature_weights.tsv"),
            str(args.out_dir / "kg_ga_best_gates_and_synergies.tsv"),
            str(args.out_dir / "kg_ga_benchmark.tsv"),
            str(args.out_dir / "KG_GA_EVOLUTIONARY_ENSEMBLE_REPORT.md"),
        ],
    }
    (args.out_dir / "kg_ga_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
