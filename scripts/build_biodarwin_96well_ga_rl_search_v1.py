#!/usr/bin/env python3
"""GA/RL search over the BioDarwin 96-well interpreter outputs.

This is a frozen wetlab-facing search scaffold. It does not retrain BigMHC or
RF. It evolves a meta-score over the local model outputs already computed on
the 96-well plate.
"""
from __future__ import annotations

import json
import math
import shutil
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
PLATE = ROOT / "project/results/biodarwin_96well_interpreter_2026_05_11/biodarwin_96well_interpreter_v1.tsv"
OUT = ROOT / "project/results/biodarwin_96well_ga_rl_search_2026_05_11"
OUT.mkdir(parents=True, exist_ok=True)
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")


FEATURES = [
    "BioDarwin_public_anchor_v3",
    "BioDarwin_mode_bank_v4",
    "BioDarwin_regime_router_v1",
    "BigMHC_IM",
    "GA_public_feature_adapter_v0",
    "BioDarwin_RF_anchor_heuristic_v1",
    "RF_biophys_score",
    "BigMHC_EL",
    "score_consensus",
    "score_disagreement",
    "topline_vote_count",
]
TRANSFORMS = ["identity", "sqrt", "square", "sigmoid"]


def minmax(x: pd.Series | np.ndarray) -> np.ndarray:
    arr = np.asarray(pd.to_numeric(x, errors="coerce"), dtype=float)
    ok = np.isfinite(arr)
    out = np.zeros_like(arr, dtype=float)
    if not ok.any():
        return out
    lo = arr[ok].min()
    hi = arr[ok].max()
    if hi == lo:
        out[ok] = 0.5
    else:
        out[ok] = (arr[ok] - lo) / (hi - lo)
    return out


def transform(x: np.ndarray, name: str) -> np.ndarray:
    x = np.nan_to_num(np.asarray(x, dtype=float), nan=0.0, posinf=1.0, neginf=0.0)
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


@dataclass
class Genome:
    weights: np.ndarray
    transforms: np.ndarray
    consensus_floor: float
    disagreement_gate: float
    disagreement_penalty: float
    vote_boost: float
    sparsity_penalty: float


def normalize_weights(w: np.ndarray) -> np.ndarray:
    w = np.clip(np.asarray(w, dtype=float), 0.0, None)
    if not np.isfinite(w).all() or w.sum() <= 0:
        return np.ones_like(w) / len(w)
    return w / w.sum()


def random_genome(rng: np.random.Generator) -> Genome:
    return Genome(
        weights=normalize_weights(rng.dirichlet(np.ones(len(FEATURES)) * 0.8)),
        transforms=rng.integers(0, len(TRANSFORMS), size=len(FEATURES)),
        consensus_floor=float(rng.uniform(0.20, 0.65)),
        disagreement_gate=float(rng.uniform(0.08, 0.32)),
        disagreement_penalty=float(rng.uniform(0.00, 0.35)),
        vote_boost=float(rng.uniform(0.00, 0.15)),
        sparsity_penalty=float(rng.uniform(0.00, 0.10)),
    )


def score_genome(g: Genome, x: np.ndarray, feature_index: dict[str, int]) -> np.ndarray:
    z = np.zeros_like(x, dtype=float)
    for i, t in enumerate(g.transforms):
        z[:, i] = transform(x[:, i], TRANSFORMS[int(t)])
    raw = z @ normalize_weights(g.weights)
    raw += g.vote_boost * x[:, feature_index["topline_vote_count"]]
    raw += 0.10 * x[:, feature_index["score_consensus"]]
    raw -= g.disagreement_penalty * np.clip(x[:, feature_index["score_disagreement"]] - g.disagreement_gate, 0.0, 1.0)
    raw *= np.where(x[:, feature_index["score_consensus"]] < g.consensus_floor, 1.0 - 0.45 * g.disagreement_penalty, 1.0)
    return minmax(raw)


def metrics(y: np.ndarray, s: np.ndarray) -> dict[str, float]:
    out = {
        "n": int(len(y)),
        "positives": int(y.sum()),
        "positive_rate": float(np.mean(y)) if len(y) else np.nan,
    }
    if len(np.unique(y)) == 2:
        out["AUPRC"] = float(average_precision_score(y, s))
        out["AUROC"] = float(roc_auc_score(y, s))
    else:
        out["AUPRC"] = np.nan
        out["AUROC"] = np.nan
    order = np.argsort(-s)
    for k in [5, 10, 24]:
        kk = min(k, len(y))
        hits = int(y[order[:kk]].sum()) if kk else 0
        out[f"top{k}_hits"] = hits
        out[f"top{k}_precision"] = float(hits / kk) if kk else np.nan
    out["rank_sum"] = float(np.sum(np.where(y == 1)[0] + 1)) if len(y) else np.nan
    return out


def entropy_like(w: np.ndarray) -> float:
    w = normalize_weights(w)
    w = np.clip(w, 1e-9, 1.0)
    ent = -float(np.sum(w * np.log(w))) / math.log(len(w))
    active = float(np.sum(w > 0.05)) / len(w)
    return 0.5 * ent + 0.5 * active


def evaluate_genome(g: Genome, x: np.ndarray, y: np.ndarray, metric_mask: np.ndarray, feature_index: dict[str, int]) -> tuple[float, np.ndarray, dict[str, float]]:
    s = score_genome(g, x, feature_index)
    y_m = y[metric_mask]
    s_m = s[metric_mask]
    m = metrics(y_m, s_m)
    sparsity = entropy_like(g.weights)
    fitness = (
        0.42 * (0.0 if not np.isfinite(m["AUPRC"]) else m["AUPRC"])
        + 0.14 * (0.0 if not np.isfinite(m["AUROC"]) else m["AUROC"])
        + 0.18 * m["top5_precision"]
        + 0.16 * m["top10_precision"]
        + 0.05 * m["top24_precision"]
        + 0.05 * sparsity
        - 0.04 * g.sparsity_penalty
    )
    m["sparsity"] = sparsity
    m["fitness"] = float(fitness)
    return float(fitness), s, m


def mutate(g: Genome, rng: np.random.Generator, op: str, scale: float) -> Genome:
    w = g.weights.copy()
    t = g.transforms.copy()
    consensus_floor = g.consensus_floor
    disagreement_gate = g.disagreement_gate
    disagreement_penalty = g.disagreement_penalty
    vote_boost = g.vote_boost
    sparsity_penalty = g.sparsity_penalty
    if op == "reweight":
        idx = rng.choice(len(w), size=max(1, len(w) // 4), replace=False)
        w[idx] *= np.exp(rng.normal(0, scale, size=len(idx)))
    elif op == "transform":
        idx = rng.integers(0, len(t))
        t[idx] = rng.integers(0, len(TRANSFORMS))
    elif op == "gate":
        consensus_floor = float(np.clip(consensus_floor + rng.normal(0, scale * 0.08), 0.0, 0.9))
        disagreement_gate = float(np.clip(disagreement_gate + rng.normal(0, scale * 0.05), 0.0, 0.6))
        disagreement_penalty = float(np.clip(disagreement_penalty + rng.normal(0, scale * 0.05), 0.0, 0.6))
    elif op == "vote":
        vote_boost = float(np.clip(vote_boost + rng.normal(0, scale * 0.04), 0.0, 0.25))
    elif op == "sparsity":
        sparsity_penalty = float(np.clip(sparsity_penalty + rng.normal(0, 0.02), 0.0, 0.2))
    return Genome(
        weights=normalize_weights(w),
        transforms=t,
        consensus_floor=consensus_floor,
        disagreement_gate=disagreement_gate,
        disagreement_penalty=disagreement_penalty,
        vote_boost=vote_boost,
        sparsity_penalty=sparsity_penalty,
    )


def crossover(a: Genome, b: Genome, rng: np.random.Generator) -> Genome:
    mask = rng.random(len(a.weights)) < 0.5
    tmask = rng.random(len(a.transforms)) < 0.5
    return Genome(
        weights=normalize_weights(np.where(mask, a.weights, b.weights)),
        transforms=np.where(tmask, a.transforms, b.transforms),
        consensus_floor=float((a.consensus_floor + b.consensus_floor) / 2),
        disagreement_gate=float((a.disagreement_gate + b.disagreement_gate) / 2),
        disagreement_penalty=float((a.disagreement_penalty + b.disagreement_penalty) / 2),
        vote_boost=float((a.vote_boost + b.vote_boost) / 2),
        sparsity_penalty=float((a.sparsity_penalty + b.sparsity_penalty) / 2),
    )


def choose_operator(rng: np.random.Generator, pulls: dict[str, int], rewards: dict[str, float], gen: int) -> str:
    ops = ["reweight", "transform", "gate", "vote", "sparsity"]
    if gen < 5 or rng.random() < 0.15:
        return str(rng.choice(ops))
    total = sum(pulls.values()) + 1
    scores = {}
    for op in ops:
        n = pulls.get(op, 0) + 1
        mean = rewards.get(op, 0.0) / n
        scores[op] = mean + 0.08 * math.sqrt(math.log(total + 1) / n)
    return max(scores, key=scores.get)


def run_ga(x: np.ndarray, y: np.ndarray, metric_mask: np.ndarray, feature_index: dict[str, int], generations: int = 48, population_size: int = 96, seed: int = 29):
    rng = np.random.default_rng(seed)
    pop = [random_genome(rng) for _ in range(population_size)]
    pulls = {op: 0 for op in ["reweight", "transform", "gate", "vote", "sparsity"]}
    rewards = {op: 0.0 for op in pulls}
    best_g = pop[0]
    best_s = None
    best_fit = -1e9
    best_m: dict[str, float] = {}
    trace = []
    op_rows = []
    for gen in range(generations):
        scored = []
        for g in pop:
            fit, s, m = evaluate_genome(g, x, y, metric_mask, feature_index)
            scored.append((fit, g, s, m))
        scored.sort(key=lambda z: z[0], reverse=True)
        if scored[0][0] > best_fit:
            best_fit, best_g, best_s, best_m = scored[0]
        trace.append({"generation": gen, "best_fitness": scored[0][0], "mean_fitness": float(np.mean([z[0] for z in scored])), **scored[0][3]})
        elites = [z[1] for z in scored[: max(6, population_size // 12)]]
        candidate_pool = scored[: max(18, population_size // 2)]
        next_pop = elites.copy()
        rate = max(0.16, 0.48 * (1.0 - gen / max(1, generations)))
        while len(next_pop) < population_size:
            pidx = rng.choice(len(candidate_pool), size=2, replace=False)
            parent_fit = max(candidate_pool[pidx[0]][0], candidate_pool[pidx[1]][0])
            child = crossover(candidate_pool[pidx[0]][1], candidate_pool[pidx[1]][1], rng)
            op = choose_operator(rng, pulls, rewards, gen)
            child = mutate(child, rng, op, rate)
            child_fit, _, _ = evaluate_genome(child, x, y, metric_mask, feature_index)
            pulls[op] += 1
            rewards[op] += max(0.0, child_fit - parent_fit)
            next_pop.append(child)
        op_rows.append({"generation": gen, **{f"pulls_{k}": pulls[k] for k in pulls}, **{f"reward_{k}": rewards[k] for k in rewards}})
        pop = next_pop
    return best_g, best_s, best_m, pd.DataFrame(trace), pd.DataFrame(op_rows), pulls, rewards


def table_html(df: pd.DataFrame, max_rows: int = 40) -> str:
    show = df.head(max_rows).copy()
    for col in show.columns:
        if pd.api.types.is_float_dtype(show[col]):
            show[col] = show[col].map(lambda x: "" if pd.isna(x) else f"{x:.4f}")
    return show.to_html(index=False, escape=True, classes="data")


def write_html(summary: dict[str, object], metrics: pd.DataFrame, champion_rows: pd.DataFrame, plate: pd.DataFrame) -> None:
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>BioDarwin 96-well GA/RL search v1</title>
<style>
body {{ margin:0; background:#0f1418; color:#e8edf0; font:14px/1.55 system-ui,-apple-system,Segoe UI,sans-serif; }}
header {{ padding:34px 5vw 20px; background:#0c1014; border-bottom:1px solid #27313a; }}
h1 {{ margin:.2rem 0 .5rem; font-size:34px; line-height:1.05; }}
.lead {{ color:#cfd8dd; max-width:1100px; }}
main {{ padding:24px 5vw 56px; }}
.note {{ background:#171e24; border-left:3px solid #c59b3b; padding:10px 14px; margin:0 0 18px; }}
table.data {{ width:100%; border-collapse:collapse; font-size:13px; margin:10px 0 28px; }}
table.data th, table.data td {{ border-bottom:1px solid #2a3640; padding:7px 8px; text-align:left; vertical-align:top; }}
table.data th {{ color:#f1d58b; background:#141b21; position:sticky; top:0; }}
</style>
</head>
<body>
<header>
  <div style="color:#c59b3b;text-transform:uppercase;letter-spacing:.12em;font-weight:700;font-size:12px;">BioDarwin 96-well GA/RL search</div>
  <h1>Meta-score search over the 96-well plate</h1>
  <p class="lead">The genome evolves a weighted meta-score over local predictors. Fitness is measured on the 25 locked TG4050 control rows and uses AUPRC, AUROC, and top-k hit precision.</p>
</header>
<main>
  <div class="note">This is a search scaffold, not a frozen claim model. It should be re-run after the next wetlab batch or new locked plate.</div>
  <h2>Summary</h2>
  <pre>{json.dumps(summary, indent=2)}</pre>
  <h2>Locked-control comparison</h2>
  {table_html(metrics, 20)}
  <h2>Champion ranking</h2>
  {table_html(champion_rows, 30)}
  <h2>Top plate rows by champion</h2>
  {table_html(plate.sort_values("BioDarwin_wetlab_search_v1", ascending=False)[[
      "plate_well","source_id","peptide","hla_allele","analysis_role","wetlab_priority",
      "BioDarwin_wetlab_search_v1","score_consensus","score_disagreement","topline_vote_count"
  ]], 30)}
</main>
</body>
</html>
"""
    path = HUB / "biodarwin_96well_ga_rl_search_v1.html"
    path.write_text(html, encoding="utf-8")
    try:
        shutil.copy2(path, LIVE_HUB / path.name)
    except Exception as exc:
        print(f"deploy warning: {exc}")


def main() -> None:
    t0 = time.time()
    plate = pd.read_csv(PLATE, sep="\t")
    feature_cols = FEATURES
    for col in feature_cols:
        plate[f"{col}_norm"] = minmax(plate[col])
    x_cols = [f"{c}_norm" for c in feature_cols]
    x = plate[x_cols].to_numpy(dtype=float)
    feature_index = {c: i for i, c in enumerate(feature_cols)}
    y = pd.to_numeric(plate["label"], errors="coerce").fillna(0).astype(int).to_numpy()
    metric_mask = plate["analysis_role"].astype(str).eq("locked_metric_control").to_numpy()

    best_g, best_s, best_m, trace, op_rows, pulls, rewards = run_ga(x, y, metric_mask, feature_index)
    plate["BioDarwin_wetlab_search_v1"] = best_s
    plate["BioDarwin_wetlab_search_rank_v1"] = plate["BioDarwin_wetlab_search_v1"].rank(ascending=False, method="first").astype(int)

    baselines = []
    for algo, col in [
        ("BioDarwin_public_anchor_v3", "BioDarwin_public_anchor_v3"),
        ("BioDarwin_mode_bank_v4", "BioDarwin_mode_bank_v4"),
        ("BioDarwin_regime_router_v1", "BioDarwin_regime_router_v1"),
        ("BigMHC_IM", "BigMHC_IM"),
        ("GA_public_feature_adapter_v0", "GA_public_feature_adapter_v0"),
        ("BioDarwin_wetlab_search_v1", "BioDarwin_wetlab_search_v1"),
    ]:
        d = plate.loc[metric_mask, ["label", col]].copy()
        y_m = d["label"].astype(int).to_numpy()
        s_m = pd.to_numeric(d[col], errors="coerce").to_numpy(dtype=float)
        m = metrics(y_m, s_m)
        baselines.append({"algorithm": algo, **m})
    metrics_df = pd.DataFrame(baselines).sort_values(["AUPRC", "AUROC"], ascending=False)

    genome_rows = []
    for feat, w, t in zip(FEATURES, best_g.weights, best_g.transforms):
        genome_rows.append({"feature": feat, "weight": float(w), "transform": TRANSFORMS[int(t)]})
    genome_df = pd.DataFrame(genome_rows).sort_values("weight", ascending=False)

    summary = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "elapsed_s": round(time.time() - t0, 2),
        "plate_rows": int(len(plate)),
        "locked_metric_rows": int(metric_mask.sum()),
        "best_fitness": float(best_m.get("fitness", np.nan)),
        "best_metrics": best_m,
        "outputs": {
            "plate": str(OUT / "biodarwin_96well_ga_rl_search_v1.tsv"),
            "metrics": str(OUT / "biodarwin_96well_ga_rl_search_v1_metrics.tsv"),
            "genome": str(OUT / "biodarwin_96well_ga_rl_search_v1_genome.tsv"),
            "trace": str(OUT / "biodarwin_96well_ga_rl_search_v1_trace.tsv"),
            "html": str(HUB / "biodarwin_96well_ga_rl_search_v1.html"),
        },
    }

    plate.to_csv(OUT / "biodarwin_96well_ga_rl_search_v1.tsv", sep="\t", index=False)
    metrics_df.to_csv(OUT / "biodarwin_96well_ga_rl_search_v1_metrics.tsv", sep="\t", index=False)
    genome_df.to_csv(OUT / "biodarwin_96well_ga_rl_search_v1_genome.tsv", sep="\t", index=False)
    trace.to_csv(OUT / "biodarwin_96well_ga_rl_search_v1_trace.tsv", sep="\t", index=False)
    op_rows.to_csv(OUT / "biodarwin_96well_ga_rl_search_v1_operator_trace.tsv", sep="\t", index=False)
    (OUT / "biodarwin_96well_ga_rl_search_v1_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (OUT / "BIODARWIN_96WELL_GA_RL_SEARCH_V1.md").write_text(
        "\n".join(
            [
                "# BioDarwin 96-well GA/RL search v1",
                "",
                f"Generated: {summary['generated_at']}",
                "",
                "## Locked-control comparison",
                metrics_df.to_markdown(index=False),
                "",
                "## Best genome",
                genome_df.to_markdown(index=False),
                "",
                "## Summary",
                pd.DataFrame([summary]).to_markdown(index=False),
            ]
        ),
        encoding="utf-8",
    )
    write_html(summary, metrics_df, plate.sort_values("BioDarwin_wetlab_search_v1", ascending=False).head(30), plate)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
