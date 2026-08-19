#!/usr/bin/env python3
"""BioDarwin public router v3: LOSO bagged GA/RL mixture-of-experts.

This goes beyond the v2 full-data GA/RL router by doing leave-one-bundle-out
search on the binary public bundles, bagging the fold champions, and then
re-blending the bagged router with the full-data champion.

The gate still only sees public metadata:
- peptide length
- overlap count
- peptide priors
- HLA family
- source family
- validation type

The experts are the same public local predictors used in v2.
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd

import build_biodarwin_public_router_v2_ga_rl as v2

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
V2_DIR = ROOT / "project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v2_ga_rl"
OUT = ROOT / "project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v3_loso_ga_rl"
OUT.mkdir(parents=True, exist_ok=True)

FULL_GENOME_PATH = V2_DIR / "biodarwin_public_router_v2_genome.tsv"
FULL_BIAS_PATH = V2_DIR / "biodarwin_public_router_v2_bias.tsv"
FULL_SUMMARY_PATH = V2_DIR / "biodarwin_public_router_v2_summary.json"


def load_full_seed(expert_names: list[str], feature_names: list[str]) -> v2.Genome:
    weights = np.zeros((len(expert_names), len(feature_names)), dtype=float)
    if FULL_GENOME_PATH.exists():
        g = pd.read_csv(FULL_GENOME_PATH, sep="\t")
        ei = {e: i for i, e in enumerate(expert_names)}
        fi = {f: i for i, f in enumerate(feature_names)}
        for _, r in g.iterrows():
            e = str(r["expert"])
            f = str(r["feature"])
            if e in ei and f in fi:
                weights[ei[e], fi[f]] = float(r["weight"])
    bias = np.zeros(len(expert_names), dtype=float)
    if FULL_BIAS_PATH.exists():
        b = pd.read_csv(FULL_BIAS_PATH, sep="\t")
        bi = {e: i for i, e in enumerate(expert_names)}
        for _, r in b.iterrows():
            e = str(r["expert"])
            if e in bi:
                bias[bi[e]] = float(r["bias"])
    temperature = 0.6753977935013222
    if FULL_SUMMARY_PATH.exists():
        try:
            temperature = float(json.loads(FULL_SUMMARY_PATH.read_text(encoding="utf-8")).get("temperature", temperature))
        except Exception:
            pass
    return v2.Genome(weights=weights, bias=bias, temperature=temperature)


def seeded_ga(
    meta: np.ndarray,
    experts: np.ndarray,
    expert_mask: np.ndarray,
    y: np.ndarray,
    bundle: np.ndarray,
    bundle_names: list[str],
    expert_names: list[str],
    feature_names: list[str],
    seed_genomes: list[v2.Genome],
    generations: int,
    population_size: int,
    seed: int,
):
    rng = np.random.default_rng(seed)
    pop = list(seed_genomes)
    while len(pop) < population_size:
        pop.append(v2.random_genome(rng, len(expert_names), len(feature_names)))
    pulls = {op: 0 for op in ["jitter", "swap", "boost", "drop", "reset", "temperature", "feature_focus"]}
    rewards = {op: 0.0 for op in pulls}
    best_g = pop[0]
    best_soft = None
    best_hard = None
    best_gate = None
    best_stats: dict[str, float] = {}
    best_fit = -1e18
    trace_rows = []
    for gen in range(generations):
        scored = []
        for g in pop:
            fit, stats, soft, hard, gate = v2.evaluate_genome(g, meta, experts, expert_mask, y, bundle, bundle_names)
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
            child = v2.crossover(pool[pidx[0]][1], pool[pidx[1]][1], rng)
            op = v2.choose_operator(rng, pulls, rewards, gen)
            child = v2.mutate(child, rng, op, rate, expert_names, feature_names)
            child_fit, _, _, _, _ = v2.evaluate_genome(child, meta, experts, expert_mask, y, bundle, bundle_names)
            pulls[op] += 1
            rewards[op] += max(0.0, child_fit - parent_fit)
            next_pop.append(child)
        pop = next_pop
    return best_g, best_soft, best_hard, best_gate, best_stats, pd.DataFrame(trace_rows)


def metric_bundle(y: np.ndarray, s: np.ndarray, bundle: np.ndarray) -> pd.DataFrame:
    rows = []
    for b in sorted(pd.unique(bundle)):
        m = v2.safe_metrics(y[bundle == b], s[bundle == b])
        m["bundle"] = b
        rows.append(m)
    return pd.DataFrame(rows)


def main() -> None:
    t0 = time.time()
    df = v2.load_public_union()
    if df.empty:
        raise RuntimeError("No public prediction files found.")

    expert_names = [
        e for e in [
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
        ]
        if e == "BioDarwin_public_anchor_noEL_v4" or v2.EXPERT_COLUMNS[e] in df.columns
    ]
    meta = v2.feature_matrix(df)
    experts_raw, expert_score_cols, expert_mask = v2.expert_matrix(df, expert_names)
    y = pd.to_numeric(df["label"], errors="coerce").fillna(0).astype(int).to_numpy()
    bundle = df["bundle"].astype(str).to_numpy()
    bundle_names = sorted(df["bundle"].astype(str).unique().tolist())
    binary_bundles = [b for b in v2.PUBLIC_BINARY_BUNDLES if b in bundle_names]
    if not binary_bundles:
        raise RuntimeError("No binary public bundles available.")

    full_seed = load_full_seed(expert_names, v2.FEATURE_NAMES)
    heuristic_seed = v2.heuristic_seed(np.random.default_rng(7), expert_names, v2.FEATURE_NAMES)
    fold_rows = []
    fold_trace_rows = []
    fold_models: dict[str, v2.Genome] = {}
    fold_soft_scores: dict[str, np.ndarray] = {}
    fold_hard_experts: dict[str, np.ndarray] = {}

    for i, holdout in enumerate(binary_bundles):
        train_mask = bundle != holdout
        train_bundles = [b for b in bundle_names if b != holdout]
        seed_genomes = [
            full_seed,
            heuristic_seed,
            v2.mutate(full_seed, np.random.default_rng(100 + i), "jitter", 0.18, expert_names, v2.FEATURE_NAMES),
        ]
        best_g, best_soft_train, best_hard_train, best_gate_train, best_stats_train, trace = seeded_ga(
            meta=meta[train_mask],
            experts=experts_raw[train_mask],
            expert_mask=expert_mask[train_mask],
            y=y[train_mask],
            bundle=bundle[train_mask],
            bundle_names=train_bundles,
            expert_names=expert_names,
            feature_names=v2.FEATURE_NAMES,
            seed_genomes=seed_genomes,
            generations=24,
            population_size=48,
            seed=1000 + i,
        )
        fold_models[holdout] = best_g
        fold_trace_rows.append(trace.assign(holdout_bundle=holdout))
        soft_all, hard_all, gate_all, entropy_all = v2.apply_genome(best_g, meta, experts_raw, expert_mask)
        fold_soft_scores[holdout] = soft_all
        fold_hard_experts[holdout] = np.array([expert_names[j] for j in np.argmax(np.where(expert_mask, meta @ best_g.weights.T + best_g.bias, -1e9), axis=1)])

        holdout_mask = bundle == holdout
        holdout_metric = v2.safe_metrics(y[holdout_mask], soft_all[holdout_mask])
        train_metric = v2.safe_metrics(y[train_mask], soft_all[train_mask])
        fold_rows.append(
            {
                "holdout_bundle": holdout,
                "train_bundles": ",".join(train_bundles),
                "train_fitness": best_stats_train.get("fitness", np.nan),
                "train_mean_ap": best_stats_train.get("mean_ap", np.nan),
                "train_mean_auc": best_stats_train.get("mean_auc", np.nan),
                "train_min_ap": best_stats_train.get("min_ap", np.nan),
                "train_std_ap": best_stats_train.get("std_ap", np.nan),
                "holdout_n": int(holdout_metric["n"]),
                "holdout_n_pos": int(holdout_metric["n_pos"]),
                "holdout_n_neg": int(holdout_metric["n_neg"]),
                "holdout_AUPRC": holdout_metric["AUPRC"],
                "holdout_AUROC": holdout_metric["AUROC"],
                "holdout_top10_precision": holdout_metric["top10_precision"],
                "train_AUPRC": train_metric["AUPRC"],
                "train_AUROC": train_metric["AUROC"],
            }
        )

    fold_df = pd.DataFrame(fold_rows)
    fold_trace_df = pd.concat(fold_trace_rows, ignore_index=True) if fold_trace_rows else pd.DataFrame()

    fold_stack = np.vstack([fold_soft_scores[b] for b in binary_bundles])
    bagged_soft = np.mean(fold_stack, axis=0)
    full_soft, full_hard, full_gate, full_entropy = v2.apply_genome(full_seed, meta, experts_raw, expert_mask)

    # Optimize a single blend weight between the full-data champion and the LOSO bag.
    alpha_rows = []
    best_alpha = 0.5
    best_alpha_fit = -1e18
    best_blend = None
    best_blend_bundle = None
    for alpha in np.linspace(0.0, 1.0, 11):
        blend = alpha * bagged_soft + (1.0 - alpha) * full_soft
        m = v2.bundle_fitness(y[np.isin(bundle, binary_bundles)], blend[np.isin(bundle, binary_bundles)], bundle[np.isin(bundle, binary_bundles)])
        alpha_rows.append({"alpha_bagged": float(alpha), "fitness": float(m[0]), **m[1].assign(note="")[["bundle", "AUPRC", "AUROC", "top10_precision"]].head(0).iloc[:0].to_dict()})
        if m[0] > best_alpha_fit:
            best_alpha_fit = m[0]
            best_alpha = float(alpha)
            best_blend = blend
            best_blend_bundle = m[1]

    # Final scores.
    df["v3_full_soft_score"] = v2.minmax_array(full_soft)
    df["v3_bagged_soft_score"] = v2.minmax_array(bagged_soft)
    df["v3_blend_soft_score"] = v2.minmax_array(best_blend if best_blend is not None else bagged_soft)
    df["v3_full_hard_expert"] = [expert_names[i] for i in np.argmax(np.where(expert_mask, meta @ full_seed.weights.T + full_seed.bias, -1e9), axis=1)]
    bagged_hard_idx = np.argmax(np.where(expert_mask, meta @ np.mean(np.stack([g.weights for g in fold_models.values()], axis=0), axis=0).T + np.mean(np.stack([g.bias for g in fold_models.values()], axis=0), axis=0), -1e9), axis=1)
    df["v3_bagged_hard_expert"] = [expert_names[i] for i in bagged_hard_idx]

    rows = []
    for b in bundle_names:
        d = df[df["bundle"].eq(b)].copy()
        m_full = v2.safe_metrics(d["label"].astype(int).to_numpy(), d["v3_full_soft_score"].to_numpy())
        m_bag = v2.safe_metrics(d["label"].astype(int).to_numpy(), d["v3_bagged_soft_score"].to_numpy())
        m_blend = v2.safe_metrics(d["label"].astype(int).to_numpy(), d["v3_blend_soft_score"].to_numpy())
        hard_share = d["v3_bagged_hard_expert"].value_counts(normalize=True).sort_values(ascending=False)
        rows.append(
            {
                "bundle": b,
                "n": int(len(d)),
                "n_pos": int(d["label"].sum()),
                "n_neg": int(len(d) - d["label"].sum()),
                "full_AUPRC": m_full["AUPRC"],
                "full_AUROC": m_full["AUROC"],
                "bagged_AUPRC": m_bag["AUPRC"],
                "bagged_AUROC": m_bag["AUROC"],
                "blend_AUPRC": m_blend["AUPRC"],
                "blend_AUROC": m_blend["AUROC"],
                "blend_top10_precision": m_blend["top10_precision"],
                "dominant_bagged_hard_expert": hard_share.index[0] if len(hard_share) else "",
                "dominant_bagged_hard_share": float(hard_share.iloc[0]) if len(hard_share) else np.nan,
            }
        )
    bundle_df = pd.DataFrame(rows)

    # Compare against v2 and the prior v1 lookup.
    v2_rows = pd.read_csv(V2_DIR / "biodarwin_public_router_v2_bundle_summary.tsv", sep="\t")
    v2_cmp = v2_rows[v2_rows["bundle"].isin(binary_bundles)][["bundle", "soft_AUPRC", "soft_AUROC", "hard_AUPRC", "hard_AUROC"]].copy()
    v1_path = ROOT / "project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v1/public_router_policy.tsv"
    v1 = pd.read_csv(v1_path, sep="\t")
    v1_perf = v1[v1["router_class"].eq("performance-max") & v1["benchmark_family"].isin(binary_bundles)].copy()
    v1_safe = v1[v1["router_class"].eq("reviewer-safe") & v1["benchmark_family"].isin(binary_bundles)].copy()
    v1_eval_rows = []
    for _, r in v1_perf.iterrows():
        b = str(r["benchmark_family"])
        algo = str(r["chosen_algorithm"])
        col = "score_biodarwin_public_anchor_noEL_v4" if algo == "BioDarwin_public_anchor_noEL_v4" else v2.EXPERT_COLUMNS.get(algo)
        d = df[df["bundle"].eq(b)]
        if col and col in d.columns:
            m = v2.safe_metrics(d["label"].astype(int).to_numpy(), pd.to_numeric(d[col], errors="coerce").to_numpy(dtype=float))
            v1_eval_rows.append({"bundle": b, "chosen_algorithm": algo, "router_class": "performance-max", **m})
    for _, r in v1_safe.iterrows():
        b = str(r["benchmark_family"])
        algo = str(r["chosen_algorithm"])
        col = "score_biodarwin_public_anchor_noEL_v4" if algo == "BioDarwin_public_anchor_noEL_v4" else v2.EXPERT_COLUMNS.get(algo)
        d = df[df["bundle"].eq(b)]
        if col and col in d.columns:
            m = v2.safe_metrics(d["label"].astype(int).to_numpy(), pd.to_numeric(d[col], errors="coerce").to_numpy(dtype=float))
            v1_eval_rows.append({"bundle": b, "chosen_algorithm": algo, "router_class": "reviewer-safe", **m})
    v1_eval_df = pd.DataFrame(v1_eval_rows)

    summary = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "elapsed_s": round(time.time() - t0, 2),
        "rows": int(len(df)),
        "binary_bundles": binary_bundles,
        "holdout_folds": int(len(fold_df)),
        "best_alpha_bagged": best_alpha,
        "alpha_fitness": float(best_alpha_fit),
        "v2_mean_soft_ap": float(v2_cmp["soft_AUPRC"].mean()),
        "v3_mean_full_ap": float(bundle_df[bundle_df["bundle"].isin(binary_bundles)]["full_AUPRC"].mean()),
        "v3_mean_bagged_ap": float(bundle_df[bundle_df["bundle"].isin(binary_bundles)]["bagged_AUPRC"].mean()),
        "v3_mean_blend_ap": float(bundle_df[bundle_df["bundle"].isin(binary_bundles)]["blend_AUPRC"].mean()),
        "v3_mean_blend_auc": float(bundle_df[bundle_df["bundle"].isin(binary_bundles)]["blend_AUROC"].mean()),
        "outputs": {
            "rows": str(OUT / "biodarwin_public_router_v3_rows.tsv"),
            "bundle_summary": str(OUT / "biodarwin_public_router_v3_bundle_summary.tsv"),
            "fold_summary": str(OUT / "biodarwin_public_router_v3_fold_summary.tsv"),
            "alpha_grid": str(OUT / "biodarwin_public_router_v3_alpha_grid.tsv"),
            "v1_eval": str(OUT / "biodarwin_public_router_v3_v1_eval.tsv"),
            "trace": str(OUT / "biodarwin_public_router_v3_fold_trace.tsv"),
            "report": str(OUT / "PUBLIC_ROUTER_REPORT_V3.md"),
        },
    }

    df.to_csv(OUT / "biodarwin_public_router_v3_rows.tsv", sep="\t", index=False)
    bundle_df.to_csv(OUT / "biodarwin_public_router_v3_bundle_summary.tsv", sep="\t", index=False)
    fold_df.to_csv(OUT / "biodarwin_public_router_v3_fold_summary.tsv", sep="\t", index=False)
    fold_trace_df.to_csv(OUT / "biodarwin_public_router_v3_fold_trace.tsv", sep="\t", index=False)
    pd.DataFrame([{"alpha_bagged": float(x), "fitness": float(v2.bundle_fitness(y[np.isin(bundle, binary_bundles)], (x * bagged_soft + (1.0 - x) * full_soft)[np.isin(bundle, binary_bundles)], bundle[np.isin(bundle, binary_bundles)])[0])} for x in np.linspace(0.0, 1.0, 11)]).to_csv(OUT / "biodarwin_public_router_v3_alpha_grid.tsv", sep="\t", index=False)
    v1_eval_df.to_csv(OUT / "biodarwin_public_router_v3_v1_eval.tsv", sep="\t", index=False)
    (OUT / "PUBLIC_ROUTER_REPORT_V3.md").write_text(
        "\n".join(
            [
                "# BioDarwin public-router v3",
                "",
                "## What changed",
                "- v2 used a single GA/RL router on the full public union.",
                "- v3 does leave-one-bundle-out GA/RL, then bags the fold champions, then blends bagged vs full-data router.",
                "- This is the more generalizable public router scaffold.",
                "",
                "## Summary",
                pd.DataFrame([summary]).to_markdown(index=False),
                "",
                "## Fold champions",
                fold_df.to_markdown(index=False),
                "",
                "## Binary bundle summary",
                bundle_df.to_markdown(index=False),
                "",
                "## v2 comparison",
                v2_cmp.to_markdown(index=False),
                "",
                "## v1 comparison",
                v1_eval_df.to_markdown(index=False),
                "",
                "## Alpha grid",
                pd.read_csv(OUT / "biodarwin_public_router_v3_alpha_grid.tsv", sep="\t").to_markdown(index=False),
                "",
                "## Notes",
                "- The reported champion is the bagged/full blend with the best public-bundle fitness.",
                "- The hard router is only for inspection; the soft blend is the deployable score.",
            ]
        ),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
