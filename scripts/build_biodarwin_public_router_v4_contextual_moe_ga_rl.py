#!/usr/bin/env python3
"""BioDarwin public router v4: contextual MoE over the v3 routers.

This is the next public-router step after v3 bagging. The experts are the three
strong v3 routers:

- full-data router
- LOSO-bagged router
- best contextual blend router

The gate only sees public metadata:

- peptide length
- overlap count
- peptide priors
- HLA family
- source family
- validation type

The GA searches a contextual gate over the three experts to maximize
bundle-wise AUPRC/AUROC/top-k on the binary public bundles.
"""
from __future__ import annotations

import json
import math
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

import build_biodarwin_public_router_v2_ga_rl as v2

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
V3_DIR = ROOT / "project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v3_loso_ga_rl"
OUT = ROOT / "project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v4_contextual_moe_ga_rl"
OUT.mkdir(parents=True, exist_ok=True)

V3_ROWS = V3_DIR / "biodarwin_public_router_v3_rows.tsv"
V3_SUMMARY = V3_DIR / "biodarwin_public_router_v3_bundle_summary.tsv"


FEATURES = v2.FEATURE_NAMES
EXPERTS = ["v3_full_soft_score", "v3_bagged_soft_score", "v3_blend_soft_score"]


def safe_metrics(y: np.ndarray, s: np.ndarray) -> dict[str, float | int]:
    return v2.safe_metrics(y, s)


def load_v3_rows() -> pd.DataFrame:
    if not V3_ROWS.exists():
        raise FileNotFoundError(V3_ROWS)
    df = pd.read_csv(V3_ROWS, sep="\t")
    needed = {"label", "bundle", "v3_full_soft_score", "v3_bagged_soft_score", "v3_blend_soft_score"}
    missing = needed - set(df.columns)
    if missing:
        raise RuntimeError(f"Missing columns in v3 rows: {sorted(missing)}")
    return df


def feature_matrix(df: pd.DataFrame) -> np.ndarray:
    arr = df[FEATURES].to_numpy(dtype=float)
    arr = np.nan_to_num(arr, nan=0.0, posinf=1.0, neginf=0.0)
    return arr


def expert_matrix(df: pd.DataFrame) -> np.ndarray:
    arr = df[EXPERTS].to_numpy(dtype=float)
    arr = np.nan_to_num(arr, nan=0.0, posinf=1.0, neginf=0.0)
    return arr


def softmax_masked(logits: np.ndarray, temperature: float) -> np.ndarray:
    temp = float(np.clip(temperature, 0.05, 3.0))
    z = logits / temp
    zmax = np.max(z, axis=1, keepdims=True)
    z = z - zmax
    ex = np.exp(np.clip(z, -50.0, 50.0))
    denom = np.sum(ex, axis=1, keepdims=True)
    denom = np.where(denom <= 0, 1.0, denom)
    return ex / denom


@dataclass
class Genome:
    weights: np.ndarray  # (n_experts, n_features)
    bias: np.ndarray  # (n_experts,)
    temperature: float


def normalize_rows(w: np.ndarray) -> np.ndarray:
    w = np.asarray(w, dtype=float)
    return np.clip(w, -6.0, 6.0)


def random_genome(rng: np.random.Generator, n_experts: int, n_features: int) -> Genome:
    return Genome(
        weights=rng.normal(0.0, 0.28, size=(n_experts, n_features)),
        bias=rng.normal(0.0, 0.08, size=n_experts),
        temperature=float(rng.uniform(0.18, 0.95)),
    )


def heuristic_seed(rng: np.random.Generator, n_features: int) -> Genome:
    g = Genome(
        weights=np.zeros((len(EXPERTS), n_features), dtype=float),
        bias=np.zeros(len(EXPERTS), dtype=float),
        temperature=0.30,
    )
    fi = {f: i for i, f in enumerate(FEATURES)}

    def add(expert: str, feature: str, value: float) -> None:
        if feature in fi:
            g.weights[EXPERTS.index(expert), fi[feature]] += value

    # Full-data router is best when metadata is not heavily shifted.
    add("v3_full_soft_score", "validation_external_test", 0.50)
    add("v3_full_soft_score", "source_cedar", 0.45)
    add("v3_full_soft_score", "source_nepdb", 0.40)
    add("v3_full_soft_score", "overlap_norm", -0.20)

    # Bagged router helps source-shift and overlap-heavy regimes.
    add("v3_bagged_soft_score", "source_tesla", 0.85)
    add("v3_bagged_soft_score", "source_itsndb", 0.65)
    add("v3_bagged_soft_score", "source_nepdb", 0.20)
    add("v3_bagged_soft_score", "overlap_norm", 0.40)
    add("v3_bagged_soft_score", "validation_held_out", 0.20)

    # v3 blend is the default safe expert.
    add("v3_blend_soft_score", "validation_partial_overlap", 0.55)
    add("v3_blend_soft_score", "source_cedar", 0.35)
    add("v3_blend_soft_score", "source_tesla", 0.35)
    add("v3_blend_soft_score", "seq_cytotoxic_prior", 0.15)
    add("v3_blend_soft_score", "seq_helper_prior", 0.06)
    return g


def apply_genome(g: Genome, meta: np.ndarray, experts: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    logits = meta @ g.weights.T + g.bias
    gate = softmax_masked(logits, g.temperature)
    score = np.sum(gate * experts, axis=1)
    return score, gate, logits


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
        0.55 * mean_ap
        + 0.15 * mean_auc
        + 0.12 * mean_top10
        + 0.08 * mean_top5
        + 0.10 * min_ap
        - 0.05 * std_ap
    ), df


def entropy_penalty(gate: np.ndarray) -> float:
    gate = np.clip(gate, 1e-12, 1.0)
    ent = -np.sum(gate * np.log(gate), axis=1) / math.log(gate.shape[1])
    return float(np.mean(ent))


def mutate(g: Genome, rng: np.random.Generator, op: str, scale: float) -> Genome:
    w = g.weights.copy()
    b = g.bias.copy()
    temp = g.temperature
    if op == "jitter":
        rows = rng.choice(w.shape[0], size=1, replace=False)
        cols = rng.choice(w.shape[1], size=max(2, w.shape[1] // 5), replace=False)
        w[np.ix_(rows, cols)] += rng.normal(0.0, scale, size=(len(rows), len(cols)))
    elif op == "swap":
        a, c = rng.choice(w.shape[0], size=2, replace=False)
        w[[a, c]] = w[[c, a]]
        b[[a, c]] = b[[c, a]]
    elif op == "boost":
        r = int(rng.integers(0, w.shape[0]))
        cols = rng.choice(w.shape[1], size=max(1, w.shape[1] // 4), replace=False)
        w[r, cols] += rng.normal(0.0, 0.7 * scale, size=len(cols))
        b[r] += rng.normal(0.0, 0.04)
    elif op == "drop":
        rows = rng.choice(w.shape[0], size=1, replace=False)
        cols = rng.choice(w.shape[1], size=max(1, w.shape[1] // 4), replace=False)
        w[np.ix_(rows, cols)] *= rng.uniform(0.0, 0.4, size=(len(rows), len(cols)))
    elif op == "temperature":
        temp = float(np.clip(temp + rng.normal(0.0, 0.10 * scale), 0.08, 1.8))
    elif op == "feature_focus":
        r = int(rng.integers(0, w.shape[0]))
        c = int(rng.integers(0, w.shape[1]))
        w[r, c] += rng.normal(0.0, 1.1 * scale)
    return Genome(weights=normalize_rows(w), bias=np.clip(b, -3.0, 3.0), temperature=temp)


def crossover(a: Genome, b: Genome, rng: np.random.Generator) -> Genome:
    mask = rng.random(a.weights.shape) < 0.5
    w = np.where(mask, a.weights, b.weights)
    bmask = rng.random(a.bias.shape) < 0.5
    bias = np.where(bmask, a.bias, b.bias)
    return Genome(weights=w, bias=bias, temperature=float((a.temperature + b.temperature) / 2))


def choose_operator(rng: np.random.Generator, pulls: dict[str, int], rewards: dict[str, float], gen: int) -> str:
    ops = ["jitter", "swap", "boost", "drop", "temperature", "feature_focus"]
    if gen < 4 or rng.random() < 0.15:
        return str(rng.choice(ops))
    total = sum(pulls.values()) + 1
    scores = {}
    for op in ops:
        n = pulls.get(op, 0) + 1
        mean = rewards.get(op, 0.0) / n
        scores[op] = mean + 0.08 * math.sqrt(math.log(total + 1) / n)
    return max(scores, key=scores.get)


def evaluate_genome(
    g: Genome,
    meta: np.ndarray,
    experts: np.ndarray,
    y: np.ndarray,
    bundle: np.ndarray,
) -> tuple[float, dict[str, float], np.ndarray, np.ndarray]:
    score, gate, _ = apply_genome(g, meta, experts)
    fitness, bundle_df = bundle_fitness(y, score, bundle)
    fitness -= 0.03 * entropy_penalty(gate)
    valid = bundle_df[bundle_df["binary_metric_possible"].astype(bool)]
    stats = {
        "fitness": float(fitness),
        "mean_ap": float(valid["AUPRC"].mean()) if not valid.empty else np.nan,
        "mean_auc": float(valid["AUROC"].mean()) if not valid.empty else np.nan,
        "mean_top10": float(valid["top10_precision"].mean()) if not valid.empty else np.nan,
        "mean_top5": float(valid["top5_precision"].mean()) if not valid.empty else np.nan,
        "min_ap": float(valid["AUPRC"].min()) if not valid.empty else np.nan,
        "std_ap": float(valid["AUPRC"].std(ddof=0)) if len(valid) > 1 else 0.0,
        "entropy": float(np.mean(-np.sum(np.clip(gate, 1e-12, 1.0) * np.log(np.clip(gate, 1e-12, 1.0)), axis=1) / math.log(gate.shape[1]))),
    }
    return float(fitness), stats, score, gate


def run_ga(
    meta: np.ndarray,
    experts: np.ndarray,
    y: np.ndarray,
    bundle: np.ndarray,
    generations: int = 48,
    population_size: int = 80,
    seed: int = 77,
):
    rng = np.random.default_rng(seed)
    pop = [heuristic_seed(rng, meta.shape[1])] + [random_genome(rng, len(EXPERTS), meta.shape[1]) for _ in range(population_size - 1)]
    pulls = {op: 0 for op in ["jitter", "swap", "boost", "drop", "temperature", "feature_focus"]}
    rewards = {op: 0.0 for op in pulls}
    best_g = pop[0]
    best_score = None
    best_gate = None
    best_stats: dict[str, float] = {}
    best_fit = -1e18
    trace_rows = []
    op_rows = []
    for gen in range(generations):
        scored = []
        for g in pop:
            fit, stats, score, gate = evaluate_genome(g, meta, experts, y, bundle)
            scored.append((fit, g, stats, score, gate))
        scored.sort(key=lambda x: x[0], reverse=True)
        top_fit, top_g, top_stats, top_score, top_gate = scored[0]
        if top_fit > best_fit:
            best_fit = top_fit
            best_g = top_g
            best_score = top_score
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
            child = mutate(child, rng, op, rate)
            child_fit, _, _, _ = evaluate_genome(child, meta, experts, y, bundle)
            pulls[op] += 1
            rewards[op] += max(0.0, child_fit - parent_fit)
            next_pop.append(child)
        op_rows.append({"generation": gen, **{f"pulls_{k}": pulls[k] for k in pulls}, **{f"reward_{k}": rewards[k] for k in rewards}})
        pop = next_pop
    return best_g, best_score, best_gate, best_stats, pd.DataFrame(trace_rows), pd.DataFrame(op_rows)


def main() -> None:
    t0 = time.time()
    df = load_v3_rows()
    bundle_names = sorted(df["bundle"].astype(str).unique().tolist())
    binary_bundles = [b for b in v2.PUBLIC_BINARY_BUNDLES if b in bundle_names]
    if not binary_bundles:
        raise RuntimeError("No binary public bundles available.")

    meta = feature_matrix(df)
    experts = expert_matrix(df)
    y = pd.to_numeric(df["label"], errors="coerce").fillna(0).astype(int).to_numpy()
    bundle = df["bundle"].astype(str).to_numpy()

    best_g, best_score, best_gate, best_stats, trace, op_trace = run_ga(meta, experts, y, bundle)
    df["v4_context_score"] = v2.minmax_array(best_score)
    df["v4_gate_entropy"] = -np.sum(np.clip(best_gate, 1e-12, 1.0) * np.log(np.clip(best_gate, 1e-12, 1.0)), axis=1) / math.log(best_gate.shape[1])
    gate_argmax = np.argmax(best_gate, axis=1)
    df["v4_context_expert"] = [EXPERTS[i] for i in gate_argmax]

    rows = []
    for b in bundle_names:
        d = df[df["bundle"].eq(b)].copy()
        m = safe_metrics(d["label"].astype(int).to_numpy(), d["v4_context_score"].to_numpy())
        expert_share = d["v4_context_expert"].value_counts(normalize=True).sort_values(ascending=False)
        rows.append(
            {
                "bundle": b,
                "n": int(len(d)),
                "n_pos": int(d["label"].sum()),
                "n_neg": int(len(d) - d["label"].sum()),
                "AUPRC": m["AUPRC"],
                "AUROC": m["AUROC"],
                "top10_precision": m["top10_precision"],
                "dominant_expert": expert_share.index[0] if len(expert_share) else "",
                "dominant_share": float(expert_share.iloc[0]) if len(expert_share) else np.nan,
            }
        )
    bundle_df = pd.DataFrame(rows)

    v3 = pd.read_csv(V3_SUMMARY, sep="\t") if V3_SUMMARY.exists() else pd.DataFrame()
    v3_cmp = v3[v3["bundle"].isin(binary_bundles)][["bundle", "blend_AUPRC", "blend_AUROC"]].copy() if not v3.empty else pd.DataFrame()
    v1_path = ROOT / "project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v1/public_router_policy.tsv"
    v1 = pd.read_csv(v1_path, sep="\t")
    v1_perf = v1[v1["router_class"].eq("performance-max") & v1["benchmark_family"].isin(binary_bundles)].copy()
    v1_eval_rows = []
    for _, r in v1_perf.iterrows():
        b = str(r["benchmark_family"])
        algo = str(r["chosen_algorithm"])
        col = "score_biodarwin_public_anchor_noEL_v4" if algo == "BioDarwin_public_anchor_noEL_v4" else v2.EXPERT_COLUMNS.get(algo)
        d = df[df["bundle"].eq(b)]
        if col and col in d.columns:
            m = safe_metrics(d["label"].astype(int).to_numpy(), pd.to_numeric(d[col], errors="coerce").to_numpy(dtype=float))
            v1_eval_rows.append({"bundle": b, "chosen_algorithm": algo, **m})
    v1_eval_df = pd.DataFrame(v1_eval_rows)

    summary = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "elapsed_s": round(time.time() - t0, 2),
        "rows": int(len(df)),
        "binary_bundles": binary_bundles,
        "best_fitness": float(best_stats.get("fitness", np.nan)),
        "best_stats": best_stats,
        "mean_ap": float(bundle_df[bundle_df["bundle"].isin(binary_bundles)]["AUPRC"].mean()),
        "mean_auc": float(bundle_df[bundle_df["bundle"].isin(binary_bundles)]["AUROC"].mean()),
        "mean_top10": float(bundle_df[bundle_df["bundle"].isin(binary_bundles)]["top10_precision"].mean()),
        "outputs": {
            "rows": str(OUT / "biodarwin_public_router_v4_rows.tsv"),
            "bundle_summary": str(OUT / "biodarwin_public_router_v4_bundle_summary.tsv"),
            "trace": str(OUT / "biodarwin_public_router_v4_trace.tsv"),
            "operator_trace": str(OUT / "biodarwin_public_router_v4_operator_trace.tsv"),
            "v1_eval": str(OUT / "biodarwin_public_router_v4_v1_eval.tsv"),
            "report": str(OUT / "PUBLIC_ROUTER_REPORT_V4.md"),
        },
    }

    df.to_csv(OUT / "biodarwin_public_router_v4_rows.tsv", sep="\t", index=False)
    bundle_df.to_csv(OUT / "biodarwin_public_router_v4_bundle_summary.tsv", sep="\t", index=False)
    trace.to_csv(OUT / "biodarwin_public_router_v4_trace.tsv", sep="\t", index=False)
    op_trace.to_csv(OUT / "biodarwin_public_router_v4_operator_trace.tsv", sep="\t", index=False)
    v1_eval_df.to_csv(OUT / "biodarwin_public_router_v4_v1_eval.tsv", sep="\t", index=False)
    (OUT / "PUBLIC_ROUTER_REPORT_V4.md").write_text(
        "\n".join(
            [
                "# BioDarwin public-router v4",
                "",
                "## What changed",
                "- v3 combined the best full-data, bagged, and blend routers.",
                "- v4 uses a contextual GA/RL MoE gate over those three routers.",
                "- The gate only sees public metadata; no dataset-name lookup is used.",
                "",
                "## Summary",
                pd.DataFrame([summary]).to_markdown(index=False),
                "",
                "## Bundle summary",
                bundle_df.to_markdown(index=False),
                "",
                "## v1 performance-max comparison",
                v1_eval_df.to_markdown(index=False),
                "",
                "## v3 comparison",
                v3_cmp.to_markdown(index=False) if not v3_cmp.empty else "_none_",
                "",
                "## Trace head",
                trace.head(12).to_markdown(index=False),
                "",
                "## Notes",
                "- Deployable score is `v4_context_score`.",
                "- Expert vote is inspection-only.",
            ]
        ),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
