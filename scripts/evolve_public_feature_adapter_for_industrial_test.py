#!/usr/bin/env python3
"""Evolve a public-feature adapter on academic data, then test industrial OOD.

The industrial public locked set is never used for GA fitness. The GA evolves a
small, deployable scoring adapter using only features available in both the
academic benchmark and the industrial testset: peptide sequence, HLA allele and
BigMHC IM/EL scores. The frozen adapter is then evaluated on the public
industrial locked v0 set.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
ACADEMIC = ROOT / "project/results/cross_neo_immunogenicity_algorithm_compare_2026_05_10/cross_neo_immunogenicity_comparator_scores.tsv"
IND_DIR = ROOT / "project/results/public_industrial_neoantigen_benchmark_scout_2026_05_10"
OUT = IND_DIR / "ga_public_feature_adapter_v0"


FEATURES = [
    "bigmhc_im_norm",
    "bigmhc_el_norm",
    "peptide_quality_proxy",
    "length_pref",
    "diversity",
    "hydrophobic",
    "aromatic",
    "charged_balance",
    "hla_A",
    "hla_B",
    "hla_C",
    "im_x_quality",
    "el_x_quality",
    "im_x_el",
]
TRANSFORMS = ["identity", "sqrt", "square", "sigmoid"]


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


def transform(x: np.ndarray, name: str) -> np.ndarray:
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


def peptide_features(peptide: object) -> dict[str, float]:
    pep = str(peptide or "").upper()
    if not pep:
        return {k: 0.0 for k in ["length_pref", "diversity", "hydrophobic", "aromatic", "charged_balance", "peptide_quality_proxy"]}
    n = len(pep)
    length_pref = 1.0 if n in {9, 10, 11} else 0.75 if n in {8, 12, 13, 14, 15} else 0.0
    diversity = min(1.0, len(set(pep)) / 8.0)
    hydrophobic = sum(aa in "AILMVPFWY" for aa in pep) / n
    aromatic = sum(aa in "FWY" for aa in pep) / n
    charged = sum(aa in "KRDE" for aa in pep) / n
    charged_balance = 1.0 - min(1.0, abs(charged - 0.22) / 0.35)
    quality = 0.42 * length_pref + 0.22 * diversity + 0.16 * hydrophobic + 0.10 * aromatic + 0.10 * charged_balance
    return {
        "length_pref": float(length_pref),
        "diversity": float(diversity),
        "hydrophobic": float(hydrophobic),
        "aromatic": float(aromatic),
        "charged_balance": float(charged_balance),
        "peptide_quality_proxy": float(quality),
    }


def normalize_hla_col(s: pd.Series) -> pd.Series:
    return s.astype(str).str.upper().str.replace("HLA-", "", regex=False)


def build_feature_frame(df: pd.DataFrame, peptide_col: str, hla_col: str, im_col: str, el_col: str) -> pd.DataFrame:
    out = df.copy()
    pf = pd.DataFrame([peptide_features(x) for x in out[peptide_col]], index=out.index)
    for col in pf.columns:
        out[col] = pf[col]
    out["bigmhc_im_norm"] = minmax_series(out[im_col]).fillna(0.0)
    out["bigmhc_el_norm"] = minmax_series(out[el_col]).fillna(0.0)
    hla = normalize_hla_col(out[hla_col])
    out["hla_A"] = hla.str.startswith("A").astype(float)
    out["hla_B"] = hla.str.startswith("B").astype(float)
    out["hla_C"] = hla.str.startswith("C").astype(float)
    out["im_x_quality"] = out["bigmhc_im_norm"] * out["peptide_quality_proxy"]
    out["el_x_quality"] = out["bigmhc_el_norm"] * out["peptide_quality_proxy"]
    out["im_x_el"] = out["bigmhc_im_norm"] * out["bigmhc_el_norm"]
    return out


@dataclass
class Genome:
    weights: np.ndarray
    transforms: np.ndarray
    im_floor: float
    el_floor: float
    floor_penalty: float
    complexity_penalty: float


def normalize_weights(w: np.ndarray) -> np.ndarray:
    w = np.clip(np.asarray(w, dtype=float), 0.0, None)
    if not np.isfinite(w).all() or w.sum() <= 0:
        return np.ones_like(w) / len(w)
    return w / w.sum()


def random_genome(rng: np.random.Generator) -> Genome:
    return Genome(
        weights=normalize_weights(rng.dirichlet(np.ones(len(FEATURES)) * 0.8)),
        transforms=rng.integers(0, len(TRANSFORMS), size=len(FEATURES)),
        im_floor=float(rng.uniform(0.0, 0.35)),
        el_floor=float(rng.uniform(0.0, 0.35)),
        floor_penalty=float(rng.uniform(0.0, 0.45)),
        complexity_penalty=float(rng.uniform(0.0, 0.05)),
    )


def score_genome(g: Genome, x: np.ndarray, feature_index: dict[str, int]) -> np.ndarray:
    z = np.zeros_like(x, dtype=float)
    for i, t in enumerate(g.transforms):
        z[:, i] = transform(x[:, i], TRANSFORMS[int(t)])
    score = z @ normalize_weights(g.weights)
    penalty = np.ones(len(score))
    penalty[x[:, feature_index["bigmhc_im_norm"]] < g.im_floor] *= 1.0 - g.floor_penalty
    penalty[x[:, feature_index["bigmhc_el_norm"]] < g.el_floor] *= 1.0 - 0.5 * g.floor_penalty
    return minmax_array(score * np.clip(penalty, 0.05, 1.0))


def safe_ap(y: np.ndarray, s: np.ndarray) -> float:
    return float(average_precision_score(y, s)) if len(np.unique(y)) == 2 else np.nan


def safe_auc(y: np.ndarray, s: np.ndarray) -> float:
    return float(roc_auc_score(y, s)) if len(np.unique(y)) == 2 else np.nan


def topk(y: np.ndarray, s: np.ndarray, k: int) -> float:
    idx = np.argsort(s)[::-1][: min(k, len(s))]
    return float(np.mean(y[idx])) if len(idx) else np.nan


def evaluate(
    g: Genome,
    x: np.ndarray,
    y: np.ndarray,
    train_mask: np.ndarray,
    low_mask: np.ndarray,
    source_masks: list[np.ndarray],
    feature_index: dict[str, int],
) -> tuple[float, np.ndarray, dict[str, float]]:
    s = score_genome(g, x, feature_index)
    train_ap = safe_ap(y[train_mask], s[train_mask])
    train_auc = safe_auc(y[train_mask], s[train_mask])
    source_aps = [safe_ap(y[m], s[m]) for m in source_masks]
    source_aps = [v for v in source_aps if np.isfinite(v)]
    source_balanced = float(np.mean(source_aps)) if source_aps else 0.0
    low_ap = safe_ap(y[train_mask & low_mask], s[train_mask & low_mask])
    top96 = topk(y[train_mask], s[train_mask], 96)
    active_weights = float(np.sum(normalize_weights(g.weights) > 0.03))
    complexity = active_weights / len(FEATURES) + g.complexity_penalty
    comps = {
        "train_ap": 0.0 if not np.isfinite(train_ap) else train_ap,
        "train_auc": 0.0 if not np.isfinite(train_auc) else train_auc,
        "source_balanced_ap": source_balanced,
        "low_leakage_ap": 0.0 if not np.isfinite(low_ap) else low_ap,
        "top96_precision": 0.0 if not np.isfinite(top96) else top96,
        "complexity": complexity,
    }
    fitness = (
        0.30 * comps["train_ap"]
        + 0.22 * comps["source_balanced_ap"]
        + 0.18 * comps["low_leakage_ap"]
        + 0.15 * comps["top96_precision"]
        + 0.10 * comps["train_auc"]
        - 0.05 * comps["complexity"]
    )
    comps["fitness"] = float(fitness)
    return float(fitness), s, comps


def mutate(g: Genome, rng: np.random.Generator, rate: float) -> Genome:
    weights = g.weights.copy()
    transforms = g.transforms.copy()
    for i in range(len(weights)):
        if rng.random() < rate:
            weights[i] *= float(np.exp(rng.normal(0.0, 0.4)))
        if rng.random() < rate * 0.25:
            transforms[i] = rng.integers(0, len(TRANSFORMS))
    return Genome(
        weights=normalize_weights(weights),
        transforms=transforms,
        im_floor=float(np.clip(g.im_floor + (rng.normal(0, 0.04) if rng.random() < rate else 0), 0, 0.6)),
        el_floor=float(np.clip(g.el_floor + (rng.normal(0, 0.04) if rng.random() < rate else 0), 0, 0.6)),
        floor_penalty=float(np.clip(g.floor_penalty + (rng.normal(0, 0.05) if rng.random() < rate else 0), 0, 0.7)),
        complexity_penalty=float(np.clip(g.complexity_penalty + (rng.normal(0, 0.01) if rng.random() < rate else 0), 0, 0.1)),
    )


def crossover(a: Genome, b: Genome, rng: np.random.Generator) -> Genome:
    mask = rng.random(len(a.weights)) < 0.5
    return Genome(
        weights=normalize_weights(np.where(mask, a.weights, b.weights)),
        transforms=np.where(mask, a.transforms, b.transforms),
        im_floor=float((a.im_floor + b.im_floor) / 2),
        el_floor=float((a.el_floor + b.el_floor) / 2),
        floor_penalty=float((a.floor_penalty + b.floor_penalty) / 2),
        complexity_penalty=float((a.complexity_penalty + b.complexity_penalty) / 2),
    )


def run_ga(df: pd.DataFrame, generations: int, population_size: int, seed: int) -> tuple[Genome, pd.DataFrame, np.ndarray, dict[str, float]]:
    rng = np.random.default_rng(seed)
    x = df[FEATURES].to_numpy(dtype=float)
    y = df["label"].astype(int).to_numpy()
    feature_index = {f: i for i, f in enumerate(FEATURES)}
    validation = df["source_name"].astype(str).str.contains("validation|_Val", case=False, regex=True).to_numpy()
    train_mask = ~validation
    low_mask = df["leakage_risk_level"].astype(str).str.lower().isin(["low", "medium"]).to_numpy()
    sources = df["source_name"].astype(str)
    source_masks = [train_mask & sources.eq(src).to_numpy() for src in sorted(sources[train_mask].unique())]
    pop = [random_genome(rng) for _ in range(population_size)]
    best_g = pop[0]
    best_fit = -1e9
    best_s = np.zeros(len(df))
    best_comps: dict[str, float] = {}
    trace = []
    for gen in range(generations):
        scored = []
        for g in pop:
            fit, s, comps = evaluate(g, x, y, train_mask, low_mask, source_masks, feature_index)
            scored.append((fit, g, s, comps))
        scored.sort(key=lambda z: z[0], reverse=True)
        if scored[0][0] > best_fit:
            best_fit, best_g, best_s, best_comps = scored[0][0], scored[0][1], scored[0][2], scored[0][3]
        trace.append({"generation": gen, "best_fitness": scored[0][0], "mean_fitness": float(np.mean([z[0] for z in scored])), **scored[0][3]})
        elites = [z[1] for z in scored[: max(4, population_size // 10)]]
        nxt = elites.copy()
        while len(nxt) < population_size:
            parents = rng.choice(len(scored[: max(12, population_size // 2)]), size=2, replace=False)
            child = crossover(scored[parents[0]][1], scored[parents[1]][1], rng)
            rate = max(0.05, 0.22 * (1 - gen / max(generations, 1)))
            nxt.append(mutate(child, rng, rate))
        pop = nxt
    return best_g, pd.DataFrame(trace), best_s, best_comps


def metrics_table(df: pd.DataFrame, label_col: str, score_specs: list[tuple[str, str]], split_col: str | None = None) -> pd.DataFrame:
    rows = []
    if split_col is None:
        splits = {"all": pd.Series(True, index=df.index)}
    else:
        splits = {str(k): df[split_col].astype(str).eq(str(k)) for k in sorted(df[split_col].astype(str).unique())}
    for split, mask in splits.items():
        d = df[mask].copy()
        y = d[label_col].astype(int).to_numpy()
        for name, col in score_specs:
            s = pd.to_numeric(d[col], errors="coerce").to_numpy()
            row = {"split": split, "algorithm": name, "n": int(len(d)), "positives": int(y.sum()), "positive_rate": float(np.mean(y))}
            row["AUPRC"] = safe_ap(y, s)
            row["AUROC"] = safe_auc(y, s)
            for k in [3, 5, 10, 15, 25, 96]:
                row[f"top{k}_precision"] = topk(y, s, k)
                row[f"top{k}_hits"] = int(round(row[f"top{k}_precision"] * min(k, len(d)))) if np.isfinite(row[f"top{k}_precision"]) else 0
            rows.append(row)
    return pd.DataFrame(rows)


def genome_tables(g: Genome) -> tuple[pd.DataFrame, pd.DataFrame]:
    weights = normalize_weights(g.weights)
    wt = pd.DataFrame(
        {
            "feature": FEATURES,
            "weight": weights,
            "transform": [TRANSFORMS[int(i)] for i in g.transforms],
        }
    ).sort_values("weight", ascending=False)
    gates = pd.DataFrame(
        [
            {"gate": "im_floor", "value": g.im_floor},
            {"gate": "el_floor", "value": g.el_floor},
            {"gate": "floor_penalty", "value": g.floor_penalty},
            {"gate": "complexity_penalty", "value": g.complexity_penalty},
        ]
    )
    return wt, gates


def plot_outputs(trace: pd.DataFrame, industrial_metrics: pd.DataFrame) -> None:
    fig_dir = OUT / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    ax.plot(trace["generation"], trace["best_fitness"], label="best", color="#2f6f73")
    ax.plot(trace["generation"], trace["mean_fitness"], label="mean", color="#c59b3b")
    ax.set_xlabel("Generation")
    ax.set_ylabel("Fitness")
    ax.set_title("GA evolution on academic public features")
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig1_ga_public_feature_convergence.png", dpi=180)
    plt.close(fig)

    show = industrial_metrics[industrial_metrics["split"].eq("industrial_locked_v0")].sort_values("AUPRC")
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    ax.barh(show["algorithm"], show["AUPRC"], color="#2f6f73", label="AUPRC")
    ax.scatter(show["AUROC"], show["algorithm"], color="#c59b3b", s=65, label="AUROC", zorder=3)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Score")
    ax.set_title("Frozen adapter on industrial locked v0")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(fig_dir / "fig2_industrial_locked_metrics.png", dpi=180)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    academic = pd.read_csv(ACADEMIC, sep="\t")
    academic_feat = build_feature_frame(academic, "peptide", "hla_allele_4digit", "bigmhc_im_score", "bigmhc_el_score")
    best, trace, score, comps = run_ga(academic_feat, generations=60, population_size=120, seed=23)
    academic_feat["ga_public_feature_adapter_score"] = score

    locked = pd.read_csv(IND_DIR / "public_industrial_locked_testset_v0.tsv", sep="\t")
    im = pd.read_csv(IND_DIR / "locked_testset_bigmhc_im_predictions.csv")
    el = pd.read_csv(IND_DIR / "locked_testset_bigmhc_el_predictions.csv")
    if "candidate_id" not in im.columns:
        im.insert(0, "candidate_id", locked["industrial_candidate_id"].to_numpy())
    if "candidate_id" not in el.columns:
        el.insert(0, "candidate_id", locked["industrial_candidate_id"].to_numpy())
    industrial = locked.merge(im[["candidate_id", "BigMHC_IM"]], left_on="industrial_candidate_id", right_on="candidate_id", how="left").drop(columns=["candidate_id"])
    industrial = industrial.merge(el[["candidate_id", "BigMHC_EL"]], left_on="industrial_candidate_id", right_on="candidate_id", how="left").drop(columns=["candidate_id"])
    industrial = build_feature_frame(industrial, "peptide", "hla_allele", "BigMHC_IM", "BigMHC_EL")
    x_ind = industrial[FEATURES].to_numpy(dtype=float)
    feature_index = {f: i for i, f in enumerate(FEATURES)}
    industrial["ga_public_feature_adapter_score"] = score_genome(best, x_ind, feature_index)
    industrial["split"] = "industrial_locked_v0"

    academic_metrics = metrics_table(
        academic_feat,
        "label",
        [
            ("BigMHC_IM", "bigmhc_im_score"),
            ("BigMHC_EL", "bigmhc_el_score"),
            ("GA_public_feature_adapter", "ga_public_feature_adapter_score"),
        ],
        split_col=None,
    )
    industrial_metrics = metrics_table(
        industrial,
        "label",
        [
            ("BigMHC_IM", "BigMHC_IM"),
            ("BigMHC_EL", "BigMHC_EL"),
            ("GA_public_feature_adapter", "ga_public_feature_adapter_score"),
        ],
        split_col="split",
    )
    validation_mask = academic_feat["source_name"].astype(str).str.contains("validation|_Val", case=False, regex=True)
    validation_like = academic_feat.copy()
    validation_like["split"] = np.where(validation_mask, "academic_validation_like", "academic_train_sources")
    academic_split_metrics = metrics_table(
        validation_like,
        "label",
        [
            ("BigMHC_IM", "bigmhc_im_score"),
            ("BigMHC_EL", "bigmhc_el_score"),
            ("GA_public_feature_adapter", "ga_public_feature_adapter_score"),
        ],
        split_col="split",
    )

    weights, gates = genome_tables(best)
    academic_feat.to_csv(OUT / "academic_public_feature_adapter_scores.tsv", sep="\t", index=False)
    industrial.to_csv(OUT / "industrial_locked_v0_ga_public_feature_scores.tsv", sep="\t", index=False)
    trace.to_csv(OUT / "ga_public_feature_generation_trace.tsv", sep="\t", index=False)
    weights.to_csv(OUT / "ga_public_feature_best_weights.tsv", sep="\t", index=False)
    gates.to_csv(OUT / "ga_public_feature_best_gates.tsv", sep="\t", index=False)
    pd.concat([academic_metrics, academic_split_metrics, industrial_metrics], ignore_index=True).to_csv(
        OUT / "ga_public_feature_adapter_metrics.tsv", sep="\t", index=False
    )
    plot_outputs(trace, industrial_metrics)
    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "claim_boundary": "GA evolved on academic data only; industrial locked v0 used only once as OOD test; manual QA still required",
        "fitness_components": comps,
        "industrial_best": industrial_metrics.sort_values("AUPRC", ascending=False).iloc[0].to_dict(),
        "outputs": {
            "out_dir": str(OUT),
            "weights": str(OUT / "ga_public_feature_best_weights.tsv"),
            "metrics": str(OUT / "ga_public_feature_adapter_metrics.tsv"),
            "industrial_scores": str(OUT / "industrial_locked_v0_ga_public_feature_scores.tsv"),
        },
    }
    (OUT / "ga_public_feature_adapter_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
