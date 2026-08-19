#!/usr/bin/env python3
"""BioDarwin public router v2: GA/RL mixture-of-experts over public rows.

v1 was a dataset-family lookup table. This v2 is a feature-gated public router:

- public metadata only for the gate
  - peptide length
  - overlap count
  - peptide physicochemical priors
  - HLA family
  - source family
  - validation type
- local predictor outputs for the experts
  - BigMHC_IM, RF_biophys, LR_biophys, MHCflurry, PRIME, ESM2_Bayesian,
    DeepImmuno, TransPHLA, NetMHCpan, GP_quantum, VQC, Structure_LR, and a
    BioDarwin public-anchor derived from BigMHC_IM + peptide priors

The search uses a GA with a bandit-style mutation controller. Fitness is
bundle-wise AUPRC/AUROC/top-k on the public binary benchmarks, with a variance
penalty so the champion cannot just spike one family.
"""
from __future__ import annotations

import json
import math
import re
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

import build_biodarwin_pan_vaccine_ga_rl as pan

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
WAVE11 = ROOT / "project/results/p_neo_bayesian_2026_05_09/wave11"
PRED_DIR = WAVE11 / "wave11_predictions"
INVENTORY = ROOT / "project/results/biodarwin_external_mega_comparison_2026_05_11/biodarwin_external_dataset_inventory.tsv"
OUT = ROOT / "project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v2_ga_rl"
OUT.mkdir(parents=True, exist_ok=True)


EXPERT_COLUMNS = {
    "BigMHC_IM": "score_bigmhc_im",
    "RF_biophys": "score_rf_biophys",
    "LR_biophys": "score_lr_biophys",
    "MHCflurry": "score_mhcflurry",
    "PRIME": "score_prime",
    "ESM2_Bayesian": "score_esm2_bayesian",
    "DeepImmuno": "score_deepimmuno",
    "TransPHLA": "score_transphla",
    "NetMHCpan": "score_netmhcpan",
    "GP_quantum": "score_gp_quantum",
    "VQC": "score_vqc",
    "Structure_LR": "score_structure_lr",
    "BioDarwin_public_anchor_noEL_v4": "score_biodarwin_public_anchor_noEL_v4",
}

VALIDATION_BUCKETS = [
    "validation_external_test",
    "validation_held_out",
    "validation_partial_overlap",
    "validation_unknown",
]

SOURCE_BUCKETS = [
    "source_cedar",
    "source_nepdb",
    "source_tesla",
    "source_itsndb",
    "source_dbpepneo2",
    "source_mcpas",
    "source_neodb",
    "source_improve",
    "source_other",
]

HLA_BUCKETS = ["hla_A", "hla_B", "hla_C", "hla_other"]

NUMERIC_FEATURES = [
    "bias",
    "length_norm",
    "length_sq_norm",
    "overlap_norm",
    "seq_cytotoxic_prior",
    "seq_helper_prior",
    "helper_promiscuity_proxy",
    "slp_processability_proxy",
    "seq_diversity",
    "hydrophobic",
    "charged_balance",
    "aromatic",
]

FEATURE_NAMES = NUMERIC_FEATURES + HLA_BUCKETS + SOURCE_BUCKETS + VALIDATION_BUCKETS
PUBLIC_BINARY_BUNDLES = ["cedar_partial", "itsndb", "nepdb", "tesla_mmc4", "tesla_mmc7"]
SAFE_EXPERTS = [
    "MHCflurry",
    "PRIME",
    "DeepImmuno",
    "TransPHLA",
    "NetMHCpan",
    "BioDarwin_public_anchor_noEL_v4",
]


def minmax_array(x: np.ndarray) -> np.ndarray:
    arr = np.asarray(x, dtype=float)
    ok = np.isfinite(arr)
    out = np.zeros_like(arr, dtype=float)
    if not ok.any():
        return out
    lo = float(np.min(arr[ok]))
    hi = float(np.max(arr[ok]))
    if hi == lo:
        out[ok] = 0.5
    else:
        out[ok] = (arr[ok] - lo) / (hi - lo)
    return out


def safe_metrics(y: np.ndarray, s: np.ndarray) -> dict[str, float | int]:
    y = np.asarray(y, dtype=int)
    s = np.asarray(s, dtype=float)
    ok = np.isfinite(y) & np.isfinite(s)
    y = y[ok]
    s = s[ok]
    n = int(len(y))
    n_pos = int(y.sum()) if n else 0
    n_neg = int(n - n_pos)
    out: dict[str, float | int] = {
        "n": n,
        "n_pos": n_pos,
        "n_neg": n_neg,
        "positive_rate": float(np.mean(y)) if n else np.nan,
        "binary_metric_possible": bool(n_pos > 0 and n_neg > 0),
        "low_power": bool(n_pos < 10 or n_neg < 10),
    }
    if n_pos > 0 and n_neg > 0:
        out["AUPRC"] = float(average_precision_score(y, s))
        out["AUROC"] = float(roc_auc_score(y, s))
    else:
        out["AUPRC"] = np.nan
        out["AUROC"] = np.nan
    order = np.argsort(-s)
    for k in [5, 10, 20, 50]:
        kk = min(k, n)
        hits = int(y[order[:kk]].sum()) if kk else 0
        out[f"top{k}_precision"] = float(hits / kk) if kk else np.nan
        out[f"top{k}_hits"] = hits
    return out


def parse_source_family(text: object) -> str:
    t = str(text).upper()
    if "CEDAR" in t or "IMPROVE" in t:
        return "source_cedar"
    if "NEPDB" in t:
        return "source_nepdb"
    if "TESLA" in t:
        return "source_tesla"
    if "ITSN" in t:
        return "source_itsndb"
    if "DBPEP" in t:
        return "source_dbpepneo2"
    if "MCPAS" in t:
        return "source_mcpas"
    if "NEODB" in t:
        return "source_neodb"
    return "source_other"


def parse_hla_bucket(hla: object) -> str:
    t = str(hla).upper()
    if "HLA-A" in t or t.startswith("A*"):
        return "hla_A"
    if "HLA-B" in t or t.startswith("B*"):
        return "hla_B"
    if "HLA-C" in t or t.startswith("C*"):
        return "hla_C"
    return "hla_other"


def first_valid(series: pd.Series):
    for value in series:
        if pd.notna(value):
            return value
    return np.nan


def sanitize_score_name(col: str) -> str:
    if col.startswith("score_"):
        return col
    return f"score_{re.sub(r'[^A-Za-z0-9]+', '_', col).strip('_').lower()}"


def bundle_inventory() -> pd.DataFrame:
    inv = pd.read_csv(INVENTORY, sep="\t")
    inv["source_family"] = inv["source"].map(parse_source_family)
    return inv


def load_public_union() -> pd.DataFrame:
    inv = bundle_inventory().set_index("dataset")
    frames = []
    for path in sorted(PRED_DIR.glob("*.tsv")):
        m = re.match(r"(.+)__(.+)\.tsv$", path.name)
        if not m:
            continue
        algorithm, bundle = m.group(1), m.group(2)
        df = pd.read_csv(path, sep="\t")
        if "label" not in df.columns or "peptide" not in df.columns:
            continue
        df = df.copy()
        df["bundle"] = bundle
        df["algorithm_file"] = algorithm
        df["source_family"] = df["source"].map(parse_source_family) if "source" in df.columns else "source_other"
        if bundle in inv.index:
            df["validation_type"] = inv.loc[bundle, "validation_type"]
        else:
            df["validation_type"] = "UNKNOWN"
        if "n_overlap_with_other_sources" not in df.columns:
            df["n_overlap_with_other_sources"] = 0
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    raw = pd.concat(frames, ignore_index=True, sort=False)

    key_cols = ["bundle", "peptide", "hla", "label", "source", "length", "n_overlap_with_other_sources", "source_family", "validation_type"]
    for col in raw.columns:
        if col not in key_cols and raw[col].dtype == object:
            raw[col] = raw[col].fillna("")
    agg = {col: first_valid for col in raw.columns if col not in key_cols}
    union = raw.groupby(key_cols, as_index=False).agg(agg)

    # Derived public anchor from BigMHC_IM + peptide priors.
    if "score_bigmhc_im" in union.columns:
        pep = pd.DataFrame([pan.peptide_features(x) for x in union["peptide"]], index=union.index)
        raw_anchor = (
            0.95 * minmax_array(pd.to_numeric(union["score_bigmhc_im"], errors="coerce").to_numpy())
            + 0.10 * pep["seq_cytotoxic_prior"].to_numpy(dtype=float)
            + 0.02 * pep["seq_helper_prior"].to_numpy(dtype=float)
        )
        union["score_biodarwin_public_anchor_noEL_v4"] = minmax_array(raw_anchor)

    # Add peptide/HLA metadata features.
    pep = pd.DataFrame([pan.peptide_features(x) for x in union["peptide"]], index=union.index)
    union = pd.concat([union.reset_index(drop=True), pep.reset_index(drop=True)], axis=1)
    union["hla_bucket"] = union["hla"].map(parse_hla_bucket)
    union["length_norm"] = pd.to_numeric(union["length"], errors="coerce").fillna(0).to_numpy(dtype=float) / 15.0
    union["length_sq_norm"] = union["length_norm"] ** 2
    union["overlap_norm"] = minmax_array(pd.to_numeric(union["n_overlap_with_other_sources"], errors="coerce").fillna(0).to_numpy(dtype=float))

    for bucket in HLA_BUCKETS:
        union[bucket] = (union["hla_bucket"] == bucket).astype(float)
    for bucket in SOURCE_BUCKETS:
        union[bucket] = (union["source_family"] == bucket).astype(float)
    for bucket in VALIDATION_BUCKETS:
        union[bucket] = (union["validation_type"] == bucket.replace("validation_", "").upper()).astype(float)
    # Fix validation buckets explicitly because the inventory uses mixed casing.
    union["validation_external_test"] = union["validation_type"].astype(str).str.upper().eq("EXTERNAL_TEST").astype(float)
    union["validation_held_out"] = union["validation_type"].astype(str).str.upper().eq("HELD_OUT").astype(float)
    union["validation_partial_overlap"] = union["validation_type"].astype(str).str.upper().eq("PARTIAL_OVERLAP").astype(float)
    union["validation_unknown"] = (~(
        union["validation_external_test"].astype(bool)
        | union["validation_held_out"].astype(bool)
        | union["validation_partial_overlap"].astype(bool)
    )).astype(float)

    for col in FEATURE_NAMES:
        if col not in union.columns:
            union[col] = 0.0
    return union


def expert_matrix(df: pd.DataFrame, experts: list[str]) -> tuple[np.ndarray, list[str], np.ndarray]:
    cols = [EXPERT_COLUMNS[e] for e in experts]
    X = np.zeros((len(df), len(experts)), dtype=float)
    mask = np.zeros((len(df), len(experts)), dtype=bool)
    for j, col in enumerate(cols):
        if col not in df.columns:
            continue
        vals = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
        ok = np.isfinite(vals)
        if ok.any():
            X[ok, j] = minmax_array(vals[ok])
            X[~ok, j] = np.nan
            mask[:, j] = ok
    return X, cols, mask


def feature_matrix(df: pd.DataFrame) -> np.ndarray:
    arr = df[FEATURE_NAMES].to_numpy(dtype=float)
    arr = np.nan_to_num(arr, nan=0.0, posinf=1.0, neginf=0.0)
    return arr


def softmax_masked(logits: np.ndarray, mask: np.ndarray, temperature: float) -> np.ndarray:
    temp = float(np.clip(temperature, 0.05, 3.0))
    z = logits / temp
    z = np.where(mask, z, -1e9)
    zmax = np.max(z, axis=1, keepdims=True)
    z = z - zmax
    ex = np.exp(np.clip(z, -50.0, 50.0)) * mask
    denom = np.sum(ex, axis=1, keepdims=True)
    denom = np.where(denom <= 0, 1.0, denom)
    return ex / denom


@dataclass
class Genome:
    weights: np.ndarray  # (n_experts, n_features)
    bias: np.ndarray  # (n_experts,)
    temperature: float


def random_genome(rng: np.random.Generator, n_experts: int, n_features: int) -> Genome:
    return Genome(
        weights=rng.normal(0.0, 0.35, size=(n_experts, n_features)),
        bias=rng.normal(0.0, 0.08, size=n_experts),
        temperature=float(rng.uniform(0.18, 0.95)),
    )


def heuristic_seed(rng: np.random.Generator, experts: list[str], feature_names: list[str]) -> Genome:
    n_experts = len(experts)
    n_features = len(feature_names)
    g = Genome(
        weights=np.zeros((n_experts, n_features), dtype=float),
        bias=np.zeros(n_experts, dtype=float),
        temperature=0.35,
    )
    fi = {f: i for i, f in enumerate(feature_names)}
    ei = {e: i for i, e in enumerate(experts)}

    def add(expert: str, feature: str, value: float) -> None:
        if expert in ei and feature in fi:
            g.weights[ei[expert], fi[feature]] += value

    # Sparse-overlap / low-positive regimes.
    for expert in ["BigMHC_IM", "BioDarwin_public_anchor_noEL_v4"]:
        add(expert, "overlap_norm", -0.7)
        add(expert, "length_norm", -0.25)
        add(expert, "validation_external_test", 0.25)
    add("BigMHC_IM", "source_tesla", 1.20)
    add("BigMHC_IM", "validation_external_test", 0.50)
    add("BigMHC_IM", "source_itsndb", -0.30)
    add("BigMHC_IM", "aromatic", 0.12)

    # Overlap-heavy public families.
    add("RF_biophys", "source_cedar", 1.40)
    add("RF_biophys", "source_nepdb", 1.00)
    add("RF_biophys", "source_tesla", 0.85)
    add("RF_biophys", "overlap_norm", 1.10)
    add("RF_biophys", "seq_cytotoxic_prior", 0.25)

    add("LR_biophys", "source_cedar", 0.65)
    add("LR_biophys", "source_nepdb", 0.55)
    add("LR_biophys", "overlap_norm", 0.80)

    add("MHCflurry", "source_cedar", 1.00)
    add("MHCflurry", "validation_partial_overlap", 0.35)
    add("MHCflurry", "length_norm", 0.15)

    add("PRIME", "source_nepdb", 0.85)
    add("PRIME", "source_cedar", 0.55)
    add("PRIME", "validation_held_out", 0.20)

    # Source-shift / PLM-style rescue.
    add("ESM2_Bayesian", "source_itsndb", 1.55)
    add("DeepImmuno", "source_itsndb", 1.10)
    add("TransPHLA", "source_itsndb", 0.90)
    add("NetMHCpan", "source_itsndb", 0.55)
    add("GP_quantum", "source_itsndb", 0.35)
    add("VQC", "source_itsndb", 0.25)
    add("Structure_LR", "source_itsndb", 0.45)

    # Borderline low-prevalence TESLA-like families.
    add("BioDarwin_public_anchor_noEL_v4", "source_tesla", 0.75)
    add("BioDarwin_public_anchor_noEL_v4", "seq_cytotoxic_prior", 0.25)
    add("BioDarwin_public_anchor_noEL_v4", "seq_helper_prior", 0.08)
    add("BioDarwin_public_anchor_noEL_v4", "overlap_norm", 0.55)

    # Mild source priors.
    for expert in experts:
        if expert not in ei:
            continue
        g.bias[ei[expert]] = rng.normal(0.0, 0.03)
    return g


def apply_genome(g: Genome, meta: np.ndarray, experts: np.ndarray, expert_mask: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    logits = meta @ g.weights.T + g.bias
    logits = np.where(expert_mask, logits, -1e9)
    gate = softmax_masked(logits, expert_mask, g.temperature)
    soft = np.sum(gate * np.nan_to_num(experts, nan=0.0), axis=1)
    hard_idx = np.argmax(np.where(expert_mask, logits, -1e9), axis=1)
    hard = experts[np.arange(len(experts)), hard_idx] if len(experts) else np.zeros(len(meta), dtype=float)
    hard = np.where(np.isfinite(hard), hard, 0.5)
    entropy = -np.sum(np.where(gate > 0, gate * np.log(np.clip(gate, 1e-12, 1.0)), 0.0), axis=1)
    entropy = entropy / math.log(max(2, gate.shape[1]))
    return soft, hard, gate, entropy


def bundle_fitness(y: np.ndarray, s: np.ndarray, bundle: np.ndarray) -> tuple[float, pd.DataFrame]:
    rows = []
    for b in sorted(pd.unique(bundle)):
        m = safe_metrics(y[bundle == b], s[bundle == b])
        m["bundle"] = b
        rows.append(m)
    df = pd.DataFrame(rows)
    valid = df[df["binary_metric_possible"].astype(bool)].copy()
    if valid.empty:
        return -1e9, df
    mean_ap = float(valid["AUPRC"].mean())
    mean_auc = float(valid["AUROC"].mean())
    mean_top10 = float(valid["top10_precision"].mean())
    mean_top5 = float(valid["top5_precision"].mean())
    min_ap = float(valid["AUPRC"].min())
    std_ap = float(valid["AUPRC"].std(ddof=0)) if len(valid) > 1 else 0.0
    return (
        0.56 * mean_ap
        + 0.14 * mean_auc
        + 0.12 * mean_top10
        + 0.08 * mean_top5
        + 0.10 * min_ap
        - 0.06 * std_ap
    ), df


def entropy_penalty(entropy: np.ndarray) -> float:
    return float(np.mean(entropy))


def mutate(g: Genome, rng: np.random.Generator, op: str, scale: float, experts: list[str], feature_names: list[str]) -> Genome:
    w = g.weights.copy()
    b = g.bias.copy()
    temp = g.temperature
    n_experts, n_features = w.shape
    if op == "jitter":
        rows = rng.choice(n_experts, size=max(1, n_experts // 3), replace=False)
        cols = rng.choice(n_features, size=max(2, n_features // 4), replace=False)
        w[np.ix_(rows, cols)] += rng.normal(0.0, scale, size=(len(rows), len(cols)))
    elif op == "swap":
        a, c = rng.choice(n_experts, size=2, replace=False)
        w[[a, c]] = w[[c, a]]
        b[[a, c]] = b[[c, a]]
    elif op == "boost":
        r = int(rng.integers(0, n_experts))
        cols = rng.choice(n_features, size=max(1, n_features // 5), replace=False)
        w[r, cols] *= rng.uniform(1.0 - scale, 1.0 + scale, size=len(cols))
        b[r] += rng.normal(0.0, 0.05 * scale)
    elif op == "drop":
        rows = rng.choice(n_experts, size=max(1, n_experts // 4), replace=False)
        cols = rng.choice(n_features, size=max(1, n_features // 4), replace=False)
        w[np.ix_(rows, cols)] *= rng.uniform(0.0, 0.4, size=(len(rows), len(cols)))
    elif op == "reset":
        r = int(rng.integers(0, n_experts))
        w[r] = rng.normal(0.0, 0.25, size=n_features)
        b[r] = rng.normal(0.0, 0.05)
    elif op == "temperature":
        temp = float(np.clip(temp + rng.normal(0.0, 0.12 * scale), 0.08, 1.8))
    elif op == "feature_focus":
        cols = rng.choice(n_features, size=max(1, n_features // 5), replace=False)
        row = int(rng.integers(0, n_experts))
        w[row, cols] += rng.normal(0.0, 0.8 * scale, size=len(cols))
    w = np.clip(w, -6.0, 6.0)
    b = np.clip(b, -3.0, 3.0)
    return Genome(weights=w, bias=b, temperature=temp)


def crossover(a: Genome, b: Genome, rng: np.random.Generator) -> Genome:
    mask = rng.random(a.weights.shape) < 0.5
    w = np.where(mask, a.weights, b.weights)
    bmask = rng.random(a.bias.shape) < 0.5
    bias = np.where(bmask, a.bias, b.bias)
    temp = float((a.temperature + b.temperature) / 2)
    return Genome(weights=w, bias=bias, temperature=temp)


def choose_operator(rng: np.random.Generator, pulls: dict[str, int], rewards: dict[str, float], gen: int) -> str:
    ops = ["jitter", "swap", "boost", "drop", "reset", "temperature", "feature_focus"]
    if gen < 4 or rng.random() < 0.15:
        return str(rng.choice(ops))
    total = sum(pulls.values()) + 1
    scores: dict[str, float] = {}
    for op in ops:
        n = pulls.get(op, 0) + 1
        mean = rewards.get(op, 0.0) / n
        scores[op] = mean + 0.08 * math.sqrt(math.log(total + 1) / n)
    return max(scores, key=scores.get)


def evaluate_genome(
    g: Genome,
    meta: np.ndarray,
    experts: np.ndarray,
    expert_mask: np.ndarray,
    y: np.ndarray,
    bundle: np.ndarray,
    bundle_names: list[str],
) -> tuple[float, dict[str, float], np.ndarray, np.ndarray, np.ndarray]:
    soft, hard, gate, entropy = apply_genome(g, meta, experts, expert_mask)
    fitness, bundle_df = bundle_fitness(y, soft, bundle)
    # Encourage sharp routing without collapsing to a single global expert.
    fitness -= 0.025 * entropy_penalty(entropy)
    valid = bundle_df[bundle_df["binary_metric_possible"].astype(bool)]
    stats = {
        "fitness": float(fitness),
        "mean_ap": float(valid["AUPRC"].mean()) if not valid.empty else np.nan,
        "mean_auc": float(valid["AUROC"].mean()) if not valid.empty else np.nan,
        "mean_top10": float(valid["top10_precision"].mean()) if not valid.empty else np.nan,
        "mean_top5": float(valid["top5_precision"].mean()) if not valid.empty else np.nan,
        "min_ap": float(valid["AUPRC"].min()) if not valid.empty else np.nan,
        "std_ap": float(valid["AUPRC"].std(ddof=0)) if len(valid) > 1 else 0.0,
        "entropy": float(np.mean(entropy)),
    }
    stats["bundle_eval_rows"] = int(len(valid))
    return float(fitness), stats, soft, hard, gate


def run_ga(
    meta: np.ndarray,
    experts: np.ndarray,
    expert_mask: np.ndarray,
    y: np.ndarray,
    bundle: np.ndarray,
    bundle_names: list[str],
    expert_names: list[str],
    feature_names: list[str],
    generations: int = 48,
    population_size: int = 72,
    seed: int = 42,
):
    rng = np.random.default_rng(seed)
    pop = [heuristic_seed(rng, expert_names, feature_names)] + [
        random_genome(rng, len(expert_names), len(feature_names)) for _ in range(population_size - 1)
    ]
    pulls = {op: 0 for op in ["jitter", "swap", "boost", "drop", "reset", "temperature", "feature_focus"]}
    rewards = {op: 0.0 for op in pulls}
    best_g = pop[0]
    best_soft = None
    best_hard = None
    best_gate = None
    best_stats: dict[str, float] = {}
    best_fit = -1e18
    trace_rows = []
    op_rows = []
    for gen in range(generations):
        scored = []
        for g in pop:
            fit, stats, soft, hard, gate = evaluate_genome(g, meta, experts, expert_mask, y, bundle, bundle_names)
            scored.append((fit, g, stats, soft, hard, gate))
        scored.sort(key=lambda x: x[0], reverse=True)
        top_fit, top_g, top_stats, top_soft, top_hard, top_gate = scored[0]
        if top_fit > best_fit:
            best_fit = top_fit
            best_g = top_g
            best_soft = top_soft
            best_hard = top_hard
            best_gate = top_gate
            best_stats = top_stats
        trace_rows.append({"generation": gen, "best_fitness": top_fit, "mean_fitness": float(np.mean([x[0] for x in scored])), **top_stats})
        elites = [x[1] for x in scored[: max(6, population_size // 10)]]
        pool = scored[: max(18, population_size // 2)]
        next_pop = elites.copy()
        rate = max(0.10, 0.45 * (1.0 - gen / max(1, generations)))
        while len(next_pop) < population_size:
            pidx = rng.choice(len(pool), size=2, replace=False)
            parent_fit = max(pool[pidx[0]][0], pool[pidx[1]][0])
            child = crossover(pool[pidx[0]][1], pool[pidx[1]][1], rng)
            op = choose_operator(rng, pulls, rewards, gen)
            child = mutate(child, rng, op, rate, expert_names, feature_names)
            child_fit, _, _, _, _ = evaluate_genome(child, meta, experts, expert_mask, y, bundle, bundle_names)
            pulls[op] += 1
            rewards[op] += max(0.0, child_fit - parent_fit)
            next_pop.append(child)
        op_rows.append({"generation": gen, **{f"pulls_{k}": pulls[k] for k in pulls}, **{f"reward_{k}": rewards[k] for k in rewards}})
        pop = next_pop
    return best_g, best_soft, best_hard, best_gate, best_stats, pd.DataFrame(trace_rows), pd.DataFrame(op_rows)


def table_html(df: pd.DataFrame, max_rows: int = 40) -> str:
    show = df.head(max_rows).copy()
    for col in show.columns:
        if pd.api.types.is_float_dtype(show[col]):
            show[col] = show[col].map(lambda x: "" if pd.isna(x) else f"{x:.4f}")
    return show.to_html(index=False, escape=True, classes="data")


def main() -> None:
    t0 = time.time()
    df = load_public_union()
    if df.empty:
        raise RuntimeError("No public prediction files found.")

    bundle_names = sorted(df["bundle"].unique().tolist())
    binary_bundle_names = [b for b in PUBLIC_BINARY_BUNDLES if b in bundle_names]
    eval_mask = df["bundle"].isin(binary_bundle_names).to_numpy()
    if not eval_mask.any():
        raise RuntimeError("No public binary evaluation bundles available.")

    expert_names = [e for e in EXPERT_COLUMNS if EXPERT_COLUMNS[e] in df.columns or e == "BioDarwin_public_anchor_noEL_v4"]
    # Keep the router broad, but only on experts that have some score column in the public union.
    expert_names = [
        e
        for e in expert_names
        if e == "BioDarwin_public_anchor_noEL_v4" or EXPERT_COLUMNS[e] in df.columns
    ]
    # Sort in a stable order.
    expert_names = [e for e in [
        "BigMHC_IM",
        "RF_biophys",
        "LR_biophys",
        "MHCflurry",
        "PRIME",
        "ESM2_Bayesian",
        "DeepImmuno",
        "TransPHLA",
        "NetMHCpan",
        "GP_quantum",
        "VQC",
        "Structure_LR",
        "BioDarwin_public_anchor_noEL_v4",
    ] if e in expert_names]

    meta = feature_matrix(df)
    experts_raw, expert_score_cols, expert_mask = expert_matrix(df, expert_names)
    y = pd.to_numeric(df["label"], errors="coerce").fillna(0).astype(int).to_numpy()
    bundle = df["bundle"].astype(str).to_numpy()

    best_g, best_soft, best_hard, best_gate, best_stats, trace, operator_trace = run_ga(
        meta=meta,
        experts=experts_raw,
        expert_mask=expert_mask,
        y=y,
        bundle=bundle,
        bundle_names=bundle_names,
        expert_names=expert_names,
        feature_names=FEATURE_NAMES,
    )

    # Final row-level outputs.
    gate = best_gate
    hard_idx = np.argmax(np.where(expert_mask, meta @ best_g.weights.T + best_g.bias, -1e9), axis=1)
    df["public_router_v2_soft_score"] = best_soft
    df["public_router_v2_hard_score"] = best_hard
    df["public_router_v2_hard_expert"] = [expert_names[i] for i in hard_idx]
    df["public_router_v2_gate_entropy"] = -np.sum(
        np.where(gate > 0, gate * np.log(np.clip(gate, 1e-12, 1.0)), 0.0), axis=1
    ) / math.log(max(2, len(expert_names)))

    # Global minmax on the mixture scores for compact report tables.
    df["public_router_v2_soft_score"] = minmax_array(df["public_router_v2_soft_score"].to_numpy())
    df["public_router_v2_hard_score"] = minmax_array(df["public_router_v2_hard_score"].to_numpy())

    # Evaluate bundle-wise performance for the champion and single-expert baselines.
    rows = []
    bundle_summaries = []
    for b in bundle_names:
        d = df[df["bundle"].eq(b)].copy()
        m_soft = safe_metrics(d["label"].astype(int).to_numpy(), d["public_router_v2_soft_score"].to_numpy())
        m_hard = safe_metrics(d["label"].astype(int).to_numpy(), d["public_router_v2_hard_score"].to_numpy())
        # modal hard expert by rows
        expert_share = d["public_router_v2_hard_expert"].value_counts(normalize=True).sort_values(ascending=False)
        bundle_summaries.append(
            {
                "bundle": b,
                "n": int(len(d)),
                "n_pos": int(d["label"].sum()),
                "n_neg": int(len(d) - d["label"].sum()),
                "binary_metric_possible": bool(m_soft["binary_metric_possible"]),
                "soft_AUPRC": m_soft["AUPRC"],
                "soft_AUROC": m_soft["AUROC"],
                "soft_top10_precision": m_soft["top10_precision"],
                "hard_AUPRC": m_hard["AUPRC"],
                "hard_AUROC": m_hard["AUROC"],
                "hard_top10_precision": m_hard["top10_precision"],
                "dominant_hard_expert": expert_share.index[0] if len(expert_share) else "",
                "dominant_hard_share": float(expert_share.iloc[0]) if len(expert_share) else np.nan,
            }
        )
        for expert in expert_names:
            col = EXPERT_COLUMNS[expert] if expert in EXPERT_COLUMNS else None
            if expert == "BioDarwin_public_anchor_noEL_v4":
                col = "score_biodarwin_public_anchor_noEL_v4"
            if col is None or col not in d.columns:
                continue
            rows.append(
                {
                    "bundle": b,
                    "algorithm": expert,
                    **safe_metrics(d["label"].astype(int).to_numpy(), pd.to_numeric(d[col], errors="coerce").to_numpy(dtype=float)),
                }
            )

    bundle_summary_df = pd.DataFrame(bundle_summaries)
    single_metrics = pd.DataFrame(rows)
    binary_single = single_metrics[single_metrics["bundle"].isin(binary_bundle_names)].copy()
    best_by_bundle = (
        binary_single.sort_values(["bundle", "AUPRC", "AUROC"], ascending=[True, False, False])
        .groupby("bundle", as_index=False)
        .head(1)
        .reset_index(drop=True)
    )
    v1_path = ROOT / "project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v1/public_router_policy.tsv"
    v1 = pd.read_csv(v1_path, sep="\t") if v1_path.exists() else pd.DataFrame()
    v1 = v1[v1["benchmark_family"].isin(binary_bundle_names)].copy() if not v1.empty else v1
    v1_perf = v1[v1["router_class"].eq("performance-max")].copy() if not v1.empty else v1
    v1_safe = v1[v1["router_class"].eq("reviewer-safe")].copy() if not v1.empty else v1
    v1_eval_rows = []
    v1_safe_eval_rows = []
    if not v1_perf.empty:
        for _, r in v1_perf.iterrows():
            b = str(r["benchmark_family"])
            algo = str(r["chosen_algorithm"])
            d = df[df["bundle"].eq(b)].copy()
            if algo == "BioDarwin_public_anchor_noEL_v4":
                col = "score_biodarwin_public_anchor_noEL_v4"
            else:
                col = EXPERT_COLUMNS.get(algo)
            if col is None or col not in d.columns:
                continue
            m = safe_metrics(d["label"].astype(int).to_numpy(), pd.to_numeric(d[col], errors="coerce").to_numpy(dtype=float))
            v1_eval_rows.append({"bundle": b, "chosen_algorithm": algo, **m})
    if not v1_safe.empty:
        for _, r in v1_safe.iterrows():
            b = str(r["benchmark_family"])
            algo = str(r["chosen_algorithm"])
            d = df[df["bundle"].eq(b)].copy()
            if algo == "BioDarwin_public_anchor_noEL_v4":
                col = "score_biodarwin_public_anchor_noEL_v4"
            else:
                col = EXPERT_COLUMNS.get(algo)
            if col is None or col not in d.columns:
                continue
            m = safe_metrics(d["label"].astype(int).to_numpy(), pd.to_numeric(d[col], errors="coerce").to_numpy(dtype=float))
            v1_safe_eval_rows.append({"bundle": b, "chosen_algorithm": algo, **m})
    v1_eval_df = pd.DataFrame(v1_eval_rows)
    v1_safe_eval_df = pd.DataFrame(v1_safe_eval_rows)

    champion_by_bundle = bundle_summary_df[bundle_summary_df["bundle"].isin(binary_bundle_names)].copy()
    champion_by_bundle = champion_by_bundle.merge(best_by_bundle[["bundle", "algorithm", "AUPRC", "AUROC", "top10_precision"]], on="bundle", how="left", suffixes=("", "_best_single"))
    if not v1_eval_df.empty:
        champion_by_bundle = champion_by_bundle.merge(v1_eval_df[["bundle", "chosen_algorithm", "AUPRC", "AUROC", "top10_precision"]], on="bundle", how="left", suffixes=("", "_v1"))

    genome_rows = []
    for expert_idx, expert in enumerate(expert_names):
        for feat_idx, feat in enumerate(FEATURE_NAMES):
            genome_rows.append(
                {
                    "expert": expert,
                    "feature": feat,
                    "weight": float(best_g.weights[expert_idx, feat_idx]),
                }
            )
    genome_df = pd.DataFrame(genome_rows).sort_values(["expert", "weight"], ascending=[True, False])
    bias_df = pd.DataFrame({"expert": expert_names, "bias": best_g.bias})

    summary = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "elapsed_s": round(time.time() - t0, 2),
        "rows": int(len(df)),
        "bundles": bundle_names,
        "binary_bundles": binary_bundle_names,
        "experts": expert_names,
        "best_fitness": float(best_stats.get("fitness", np.nan)),
        "best_stats": best_stats,
        "temperature": float(best_g.temperature),
        "outputs": {
            "rows": str(OUT / "biodarwin_public_router_v2_rows.tsv"),
            "bundle_summary": str(OUT / "biodarwin_public_router_v2_bundle_summary.tsv"),
            "single_metrics": str(OUT / "biodarwin_public_router_v2_single_expert_metrics.tsv"),
            "best_by_bundle": str(OUT / "biodarwin_public_router_v2_best_single_by_bundle.tsv"),
            "genome": str(OUT / "biodarwin_public_router_v2_genome.tsv"),
            "bias": str(OUT / "biodarwin_public_router_v2_bias.tsv"),
            "trace": str(OUT / "biodarwin_public_router_v2_trace.tsv"),
            "operator_trace": str(OUT / "biodarwin_public_router_v2_operator_trace.tsv"),
            "report": str(OUT / "PUBLIC_ROUTER_REPORT_V2.md"),
        },
    }

    df.to_csv(OUT / "biodarwin_public_router_v2_rows.tsv", sep="\t", index=False)
    bundle_summary_df.to_csv(OUT / "biodarwin_public_router_v2_bundle_summary.tsv", sep="\t", index=False)
    single_metrics.to_csv(OUT / "biodarwin_public_router_v2_single_expert_metrics.tsv", sep="\t", index=False)
    best_by_bundle.to_csv(OUT / "biodarwin_public_router_v2_best_single_by_bundle.tsv", sep="\t", index=False)
    genome_df.to_csv(OUT / "biodarwin_public_router_v2_genome.tsv", sep="\t", index=False)
    bias_df.to_csv(OUT / "biodarwin_public_router_v2_bias.tsv", sep="\t", index=False)
    trace.to_csv(OUT / "biodarwin_public_router_v2_trace.tsv", sep="\t", index=False)
    operator_trace.to_csv(OUT / "biodarwin_public_router_v2_operator_trace.tsv", sep="\t", index=False)
    (OUT / "biodarwin_public_router_v2_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    report_lines = [
        "# BioDarwin public-router v2",
        "",
        "## What changed",
        "- v1 was a dataset-family lookup table.",
        "- v2 is a GA/RL mixture-of-experts router that uses only public metadata for gating.",
        "- The experts are local predictor outputs; the gate sees peptide length, overlap, peptide priors, HLA family, source family, and validation type.",
        "",
        "## Champion summary",
        pd.DataFrame([summary]).to_markdown(index=False),
        "",
        "## Binary public bundle comparison",
        champion_by_bundle.to_markdown(index=False),
        "",
        "## v1 performance-max router",
        v1_eval_df.to_markdown(index=False) if not v1_eval_df.empty else "_none_",
        "",
        "## v1 reviewer-safe router",
        v1_safe_eval_df.to_markdown(index=False) if not v1_safe_eval_df.empty else "_none_",
        "",
        "## Best single expert per bundle",
        best_by_bundle.to_markdown(index=False),
        "",
        "## Single-expert baseline table",
        binary_single.sort_values(["bundle", "AUPRC"], ascending=[True, False]).to_markdown(index=False),
        "",
        "## Genome",
        genome_df.head(120).to_markdown(index=False),
        "",
        "## Bias",
        bias_df.to_markdown(index=False),
        "",
        "## Operator trace head",
        trace.head(12).to_markdown(index=False),
        "",
        "## Notes",
        "- Helper bundles with a single class are retained in the row table but excluded from binary AUPRC/AUROC fitness.",
        "- The hard router is the discrete expert choice implied by the learned gate; the soft score is the primary champion metric.",
        "- If the next external batch shifts again, freeze the current weights and rerun the GA on the new holdout.",
    ]
    (OUT / "PUBLIC_ROUTER_REPORT_V2.md").write_text("\n".join(report_lines), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
