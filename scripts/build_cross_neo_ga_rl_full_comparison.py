#!/usr/bin/env python3
"""Build a consolidated GA/RL-vs-all comparison dossier for CROSS-Neo.

This script compares the current KG-GA evolved controller, fixed CROSS-Neo
integrated scores, public predictors, clean-track local scores, production
stack scores, and paperability-aware Darwin-RL ranking. It is a report/figure
builder, not a protected manuscript prose generator.
"""

from __future__ import annotations

import html
import json
import math
import shutil
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
KG = ROOT / "project/results/cross_neo_kg_ga_evolutionary_ensemble_2026_05_10"
PAPER = ROOT / "project/results/cross_neo_paperability_darwin_rl_2026_05_10"
NEO = ROOT / "project/results/neoimmune_stack_2026_05_10"
PUBLIC_CTX = ROOT / "project/results/cross_neo_product_demo_2026_05_10/strong_competitors/strong_competitor_context_public_benchmarks.tsv"
OUT = ROOT / "project/results/cross_neo_ga_rl_full_comparison_2026_05_10"
FIG = OUT / "figures"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET = HUB / "assets/cross_neo_ga_rl_full_comparison"
LIVE_ASSET = LIVE / "assets/cross_neo_ga_rl_full_comparison"
PAGE = HUB / "cross_neo_ga_rl_full_comparison.html"
LIVE_PAGE = LIVE / "cross_neo_ga_rl_full_comparison.html"


ALGORITHM_FAMILY = {
    "KG_GA_evolved": "new_ga_controller",
    "KG_GA_evolved_controller": "new_ga_controller",
    "CROSS_integrated": "fixed_integrated_crossneo",
    "KG_GA_fixed_integrated": "fixed_integrated_crossneo",
    "CROSS_claimsafe": "claim_safe_crossneo",
    "KG_GA_fixed_claimsafe": "claim_safe_crossneo",
    "CROSS_stress": "stress_guarded_internal",
    "CROSS_BMA": "bayesian_internal",
    "CROSS_finetuned": "finetuned_internal",
    "BigMHC_IM": "public_immunogenicity",
    "BigMHC_IM_inside_GA": "public_immunogenicity",
    "BigMHC_EL": "public_presentation",
    "BigMHC_EL_inside_GA": "public_presentation",
    "BAR-Neo_confidence": "production_reliability_ranker",
    "BAR-Neo": "production_reliability_ranker",
    "BAR-Neo_patient_gated": "production_reliability_ranker",
    "BAR-Neo-X": "production_reliability_ranker",
    "BAR-Neo-X_claim_safe": "production_reliability_ranker",
    "Wave8_TCR_SelfSim_no_exact": "clean_local_tcr",
    "Wave8_TCR_motif_only": "clean_local_tcr",
    "Wave8_TCR_SelfSim_full": "clean_local_tcr",
}


CLAIM_TRACK = {
    "KG_GA_evolved": "production/experiment-priority",
    "KG_GA_evolved_controller": "production/experiment-priority",
    "CROSS_integrated": "production/experiment-priority",
    "KG_GA_fixed_integrated": "production/experiment-priority",
    "CROSS_claimsafe": "claim-safe prioritization",
    "KG_GA_fixed_claimsafe": "claim-safe prioritization",
    "CROSS_stress": "claim-safe prioritization",
    "CROSS_BMA": "claim-safe prioritization",
    "CROSS_finetuned": "experiment-priority",
    "BigMHC_IM": "external comparator",
    "BigMHC_IM_inside_GA": "external comparator component",
    "BigMHC_EL": "external comparator",
    "BigMHC_EL_inside_GA": "external comparator component",
    "BAR-Neo_confidence": "strict no-overlap production reliability",
    "BAR-Neo": "strict no-overlap production reliability",
    "BAR-Neo_patient_gated": "strict no-overlap production reliability",
    "BAR-Neo-X": "strict no-overlap production reliability",
    "BAR-Neo-X_claim_safe": "strict no-overlap production reliability",
}


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t")


def esc(value: object) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    return html.escape(str(value))


def fmt(x: object, digits: int = 3) -> str:
    try:
        if pd.isna(x):
            return "NA"
        return f"{float(x):.{digits}f}"
    except Exception:
        return esc(x)


def safe_copy(src: Path, dst: Path) -> None:
    try:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    except PermissionError:
        print(f"[ga-rl-compare] skip permission-denied copy {dst}")


def table_html(df: pd.DataFrame, cols: list[str], max_rows: int = 30) -> str:
    if df.empty:
        return "<p class='muted'>No rows.</p>"
    keep = [c for c in cols if c in df.columns]
    if not keep:
        return "<p class='muted'>Requested columns not present.</p>"
    rows = ["<table><thead><tr>"]
    rows += [f"<th>{esc(c.replace('_', ' '))}</th>" for c in keep]
    rows.append("</tr></thead><tbody>")
    for _, r in df.head(max_rows).iterrows():
        rows.append("<tr>")
        for c in keep:
            v = r[c]
            rows.append(f"<td>{fmt(v) if isinstance(v, (int, float, np.integer, np.floating)) else esc(v)}</td>")
        rows.append("</tr>")
    rows.append("</tbody></table>")
    return "\n".join(rows)


def markdown_table(df: pd.DataFrame) -> str:
    view = df.copy()
    return view.where(pd.notna(view), "NA").to_markdown(index=False)


def metric_lookup(df: pd.DataFrame, split: str, algorithm: str, col: str) -> float:
    if df.empty:
        return math.nan
    m = df[df["split"].astype(str).eq(split) & df["algorithm"].astype(str).eq(algorithm)]
    if m.empty or col not in m.columns:
        return math.nan
    return float(pd.to_numeric(m.iloc[0][col], errors="coerce"))


def truthy(value: object) -> bool:
    if value is None:
        return False
    try:
        if pd.isna(value):
            return False
    except Exception:
        pass
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "y"}
    return bool(value)


def first_text(values: pd.Series) -> str:
    for v in values:
        try:
            if pd.isna(v):
                continue
        except Exception:
            pass
        s = str(v).strip()
        if s and s.lower() != "nan":
            return s
    return ""


def unique_join(values: pd.Series) -> str:
    seen: list[str] = []
    for v in values:
        try:
            if pd.isna(v):
                continue
        except Exception:
            pass
        s = str(v).strip()
        if s and s.lower() != "nan" and s not in seen:
            seen.append(s)
    return " + ".join(seen)


def collapse_algorithm_rows(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    text_join_cols = {"claim_track", "claim_boundary_note", "reason_clean_disallowed"}
    text_first_cols = {"algorithm_family", "track", "score_column"}
    bool_cols = {"clean_track_allowed"}
    for alg, g in df.groupby("algorithm", sort=False):
        row: dict[str, object] = {"algorithm": alg}
        for c in df.columns:
            if c == "algorithm":
                continue
            if c in text_join_cols:
                row[c] = unique_join(g[c])
            elif c in text_first_cols:
                row[c] = first_text(g[c])
            elif c in bool_cols:
                row[c] = any(truthy(v) for v in g[c])
            else:
                numeric = pd.to_numeric(g[c], errors="coerce")
                row[c] = float(numeric.max()) if numeric.notna().any() else first_text(g[c])
        rows.append(row)
    return pd.DataFrame(rows)


def build_master_table() -> pd.DataFrame:
    kg = read_tsv(KG / "kg_ga_benchmark.tsv")
    v6 = read_tsv(KG / "kg_ga_v6_smoke_wetlab_comparison.tsv")
    paper = read_tsv(PAPER / "paperability_algorithm_ranking.tsv")
    clean = read_tsv(NEO / "metrics/leaderboard_clean_track_strict_no_overlap.tsv")
    prod = read_tsv(NEO / "metrics/leaderboard_production_track_strict_no_overlap.tsv")
    ga_rl = read_tsv(NEO / "metrics/ga_rl_algorithm_metrics.tsv")

    algorithms = sorted(set(kg["algorithm"].astype(str)) | set(v6["algorithm"].astype(str)))
    rows: list[dict[str, object]] = []
    for alg in algorithms:
        prow = paper[paper["algorithm"].astype(str).eq(alg)]
        rows.append(
            {
                "algorithm": alg,
                "algorithm_family": ALGORITHM_FAMILY.get(alg, "other_crossneo"),
                "claim_track": CLAIM_TRACK.get(alg, "comparison only"),
                "all_AUPRC": metric_lookup(kg, "all", alg, "AUPRC"),
                "all_AUROC": metric_lookup(kg, "all", alg, "AUROC"),
                "all_top96_hits": metric_lookup(kg, "all", alg, "top96_hits"),
                "all_top96_precision": metric_lookup(kg, "all", alg, "top96_precision"),
                "frozen_validation_AUPRC": metric_lookup(kg, "frozen_validation_like_sources", alg, "AUPRC"),
                "frozen_validation_AUROC": metric_lookup(kg, "frozen_validation_like_sources", alg, "AUROC"),
                "frozen_top24_hits": metric_lookup(kg, "frozen_validation_like_sources", alg, "top24_hits"),
                "frozen_top24_precision": metric_lookup(kg, "frozen_validation_like_sources", alg, "top24_precision"),
                "low_medium_leakage_AUPRC": metric_lookup(kg, "low_medium_leakage", alg, "AUPRC"),
                "low_medium_leakage_AUROC": metric_lookup(kg, "low_medium_leakage", alg, "AUROC"),
                "low_medium_top96_hits": metric_lookup(kg, "low_medium_leakage", alg, "top96_hits"),
                "v6_smoke_AUPRC": metric_lookup(v6, "v6_smoke_confirmatory", alg, "AUPRC"),
                "v6_smoke_AUROC": metric_lookup(v6, "v6_smoke_confirmatory", alg, "AUROC"),
                "paperability_reward": float(prow.iloc[0]["paperability_reward"]) if not prow.empty else math.nan,
            }
        )

    for _, r in ga_rl.iterrows():
        alg = str(r["algorithm"])
        rows.append(
            {
                "algorithm": alg,
                "algorithm_family": ALGORITHM_FAMILY.get(alg, "neoimmune_ga_rl_component"),
                "claim_track": CLAIM_TRACK.get(alg, "production/experiment-priority"),
                "neo_stack_AUPRC": float(r.get("AUPRC", math.nan)),
                "neo_stack_AUROC": float(r.get("AUROC", math.nan)),
                "neo_precision20": float(r.get("Precision@20", math.nan)),
                "neo_precision34": float(r.get("Precision@34", math.nan)),
                "neo_patient_hit_rate34": float(r.get("patient_hit_rate@34", math.nan)),
                "clean_track_allowed": bool(r.get("clean_track_allowed", False)),
                "claim_boundary_note": r.get("reason_clean_disallowed", ""),
            }
        )

    for src, track, limit in [
        (clean, "strict no-overlap clean", 6),
        (prod, "strict no-overlap production", 12),
    ]:
        must_models = {"KG_GA_evolved_controller", "BAR-Neo_confidence", "BigMHC_IM", "BigMHC_EL"}
        src = pd.concat(
            [
                src.head(limit).copy(),
                src[src["model"].astype(str).isin(must_models)].copy(),
            ],
            ignore_index=True,
        ).drop_duplicates(subset=["track", "model", "split"], keep="first")
        for _, r in src.iterrows():
            alg = str(r["model"])
            rows.append(
                {
                    "algorithm": alg,
                    "algorithm_family": ALGORITHM_FAMILY.get(alg, "neoimmune_stack_baseline"),
                    "claim_track": track,
                    "neo_strict_AUPRC": float(r.get("AUPRC", math.nan)),
                    "neo_strict_AUROC": float(r.get("AUROC", math.nan)),
                    "neo_precision20": float(r.get("Precision@20", math.nan)),
                    "neo_patient_hit_rate20": float(r.get("patient_hit_rate@20", math.nan)),
                    "patients_evaluated": float(r.get("patients_evaluated", math.nan)),
                    "clean_track_allowed": track.endswith("clean"),
                }
            )

    out = collapse_algorithm_rows(pd.DataFrame(rows))
    metric_cols = [
        "all_AUPRC",
        "all_AUROC",
        "all_top96_hits",
        "all_top96_precision",
        "frozen_validation_AUPRC",
        "frozen_validation_AUROC",
        "frozen_top24_hits",
        "frozen_top24_precision",
        "low_medium_leakage_AUPRC",
        "low_medium_leakage_AUROC",
        "low_medium_top96_hits",
        "v6_smoke_AUPRC",
        "v6_smoke_AUROC",
        "paperability_reward",
        "neo_stack_AUPRC",
        "neo_stack_AUROC",
        "neo_strict_AUPRC",
        "neo_strict_AUROC",
        "neo_precision20",
        "neo_precision34",
        "neo_patient_hit_rate20",
        "neo_patient_hit_rate34",
        "patients_evaluated",
    ]
    for c in metric_cols:
        if c not in out.columns:
            out[c] = np.nan
        out[c] = pd.to_numeric(out[c], errors="coerce")
    out["headline_score"] = (
        out["all_AUPRC"].fillna(0) * 0.28
        + out["frozen_validation_AUPRC"].fillna(0) * 0.24
        + out["low_medium_leakage_AUPRC"].fillna(0) * 0.18
        + out["v6_smoke_AUPRC"].fillna(0) * 0.12
        + out["paperability_reward"].fillna(0) * 0.12
        + out["neo_strict_AUPRC"].fillna(0) * 0.06
    )
    out["primary_disposition"] = out.apply(disposition, axis=1)
    out = out.sort_values(["headline_score", "all_AUPRC", "neo_strict_AUPRC"], ascending=False)
    return out


def disposition(r: pd.Series) -> str:
    alg = str(r.get("algorithm", ""))
    if alg == "KG_GA_evolved":
        return "CHAMPION for retrospective + validation-like + low-leakage benchmark; freeze before prospective wetlab"
    if alg == "KG_GA_evolved_controller":
        return "Production GA/RL branch in NeoImmune stack; strong but not clean-track allowed"
    if alg == "CROSS_claimsafe":
        return "Claim-safe fallback when GA/RL proxy/product features are too risky"
    if alg == "CROSS_integrated":
        return "Fixed integrated comparator; strong but less novel than KG-GA"
    if alg.startswith("BigMHC"):
        return "External public comparator/component, not CROSS-Neo claim engine"
    if alg.startswith("BAR-Neo"):
        return "Best strict no-overlap production reliability baseline; compare separately from 2,715-board GA"
    if "Wave8" in alg:
        return "Clean local TCR baseline; useful ablation/control"
    return "Comparator or component"


def build_winloss(master: pd.DataFrame) -> pd.DataFrame:
    champion = master[master["algorithm"].eq("KG_GA_evolved")].iloc[0]
    comparators = master[
        master["algorithm"].isin(
            [
                "BigMHC_IM",
                "BigMHC_EL",
                "CROSS_stress",
                "CROSS_BMA",
                "CROSS_finetuned",
                "CROSS_integrated",
                "CROSS_claimsafe",
            ]
        )
    ].copy()
    rows = []
    metrics = [
        "all_AUPRC",
        "all_AUROC",
        "frozen_validation_AUPRC",
        "frozen_validation_AUROC",
        "low_medium_leakage_AUPRC",
        "low_medium_leakage_AUROC",
        "v6_smoke_AUPRC",
        "paperability_reward",
    ]
    for _, r in comparators.iterrows():
        row = {"champion": "KG_GA_evolved", "comparator": r["algorithm"]}
        wins = 0
        tested = 0
        for m in metrics:
            if m not in r or m not in champion:
                continue
            a = champion[m]
            b = r[m]
            if pd.notna(a) and pd.notna(b):
                tested += 1
                delta = float(a) - float(b)
                row[f"delta_{m}"] = delta
                if delta > 0:
                    wins += 1
        row["wins"] = wins
        row["tested_metrics"] = tested
        row["win_fraction"] = wins / tested if tested else np.nan
        if row["win_fraction"] == 1:
            row["disposition"] = "champion wins all tested axes"
        elif row["win_fraction"] >= 0.75:
            row["disposition"] = "champion wins most axes"
        else:
            row["disposition"] = "mixed; inspect caveats"
        rows.append(row)
    return pd.DataFrame(rows).sort_values("win_fraction", ascending=False)


def build_claim_matrix(master: pd.DataFrame) -> pd.DataFrame:
    prod = read_tsv(NEO / "metrics/leaderboard_production_track_strict_no_overlap.tsv")
    bar = prod[prod["model"].astype(str).eq("BAR-Neo_confidence")].head(1)
    kg_prod = prod[prod["model"].astype(str).eq("KG_GA_evolved_controller")].head(1)
    if not bar.empty and not kg_prod.empty:
        br = bar.iloc[0]
        kr = kg_prod.iloc[0]
        strict_evidence = (
            f"Strict no-overlap production AUPRC: BAR-Neo_confidence {float(br['AUPRC']):.3f} "
            f"(n={int(float(br['n']))}) vs KG_GA_evolved_controller {float(kr['AUPRC']):.3f} "
            f"(n={int(float(kr['n']))})."
        )
    else:
        strict_evidence = "Strict no-overlap production leaderboard does not support KG-GA as the reliability winner."
    rows = [
        {
            "claim_item": "Is GA/RL the best current algorithm on the 2,715-row retrospective CROSS-Neo board?",
            "answer": "YES",
            "evidence": "KG_GA_evolved all AUPRC 0.933 / AUROC 0.943 / top96 96/96.",
            "boundary": "Retrospective labels; not prospective wetlab proof.",
        },
        {
            "claim_item": "Does GA/RL beat public BigMHC in this local board?",
            "answer": "YES",
            "evidence": "KG_GA all AUPRC 0.933 vs BigMHC_IM 0.672 and BigMHC_EL 0.648.",
            "boundary": "BigMHC remains public comparator and useful feature, not an obsolete method.",
        },
        {
            "claim_item": "Does GA/RL remain strong on frozen validation-like sources?",
            "answer": "YES",
            "evidence": "KG_GA frozen validation-like AUPRC 0.978 / AUROC 0.999, top24 contains all 11 positives.",
            "boundary": "Validation-like split is still retrospective, not prospective clinical validation.",
        },
        {
            "claim_item": "Does GA/RL dominate low/medium leakage subset?",
            "answer": "YES",
            "evidence": "KG_GA low/medium leakage AUPRC 0.757; next strongest among CROSS comparators is finetuned 0.669.",
            "boundary": "Leakage risk reduced, not eliminated.",
        },
        {
            "claim_item": "Is GA/RL best on v6 smoke confirmatory board?",
            "answer": "NO",
            "evidence": "BigMHC_EL has v6 AUPRC 0.866; KG_GA has 0.812. KG_GA is still strong but not the v6 AUPRC winner.",
            "boundary": "Small n=24 smoke board; all top24 recover 17/24 because the board itself has 17 positives.",
        },
        {
            "claim_item": "Is GA/RL allowed as clean science track?",
            "answer": "NO",
            "evidence": "NeoImmune GA/RL includes public predictor/product-value components and is flagged clean_track_allowed=False.",
            "boundary": "Use as production/experiment-priority branch; keep clean local track separate.",
        },
        {
            "claim_item": "Does GA/RL beat strict no-overlap production reliability rankers?",
            "answer": "NO / context-specific",
            "evidence": strict_evidence,
            "boundary": "Strict no-overlap rows are model-availability-specific and separate from the 2,715-row GA board; use BAR-Neo for reliability and KG-GA for integrated discovery.",
        },
        {
            "claim_item": "Can we call it immunogenicity proof?",
            "answer": "NO",
            "evidence": "All current wins are label/benchmark/prioritization wins.",
            "boundary": "Prospective mutant > WT/decoy activation and killing assays are required.",
        },
    ]
    return pd.DataFrame(rows)


def build_ablation_gap() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "module": "Knowledge graph",
                "current_status": "implemented",
                "evidence": "51 KG nodes / 62 edges in KG-GA package; method-mining graph also exists.",
                "gap": "Need graph-node removal ablations for final method paper.",
                "priority": "high",
            },
            {
                "module": "Genetic algorithm controller",
                "current_status": "implemented",
                "evidence": "45 generations x 96 population; best all AUPRC 0.933.",
                "gap": "Need same-budget random search and GA-only vs guarded-GA ablations.",
                "priority": "high",
            },
            {
                "module": "RL controller",
                "current_status": "blueprint/paperability scaffold",
                "evidence": "Darwin-RL method and paperability reward graphs exist, but no trained RL policy score column yet.",
                "gap": "Train or simulate policy-learning over mutation/action choices; otherwise call it Darwin/RL blueprint, not trained RL.",
                "priority": "highest",
            },
            {
                "module": "Claim-safe fallback",
                "current_status": "implemented",
                "evidence": "CROSS_claimsafe top96 96/96, paperability 0.750, strong fallback.",
                "gap": "Use when production GA features are considered too risky.",
                "priority": "medium",
            },
            {
                "module": "Public comparator graft",
                "current_status": "implemented",
                "evidence": "BigMHC IM/EL scored 2,715 rows and used as external comparator/features.",
                "gap": "Public training overlap audit remains needed before manuscript feature claims.",
                "priority": "high",
            },
            {
                "module": "Prospective wetlab lock",
                "current_status": "not yet passed",
                "evidence": "v6 packet and assay interpreter exist; no real wetlab outcomes yet.",
                "gap": "Freeze KG-GA score/thresholds and run mutant/WT/decoy assay readout.",
                "priority": "highest",
            },
        ]
    )


def build_public_competitor_context() -> tuple[pd.DataFrame, pd.DataFrame]:
    if not PUBLIC_CTX.exists():
        return pd.DataFrame(), pd.DataFrame()
    raw = read_tsv(PUBLIC_CTX)
    ctx = raw[raw["method_family"].astype(str).eq("public_or_prior_internal_competitor")].copy()
    for c in ["n_scored", "n_pos", "AUPRC", "AUROC", "top5_precision", "top10_precision", "top20_precision", "coverage"]:
        if c in ctx.columns:
            ctx[c] = pd.to_numeric(ctx[c], errors="coerce")
    summary = (
        ctx.groupby("method", dropna=False)
        .agg(
            n_context_rows=("method", "size"),
            total_scored=("n_scored", "sum"),
            datasets_with_auroc=("AUROC", "count"),
            mean_context_AUROC=("AUROC", "mean"),
            best_context_AUROC=("AUROC", "max"),
            min_context_AUROC=("AUROC", "min"),
        )
        .reset_index()
        .sort_values(["best_context_AUROC", "mean_context_AUROC"], ascending=False)
    )
    summary["claim_use"] = "context-only, AUROC-only in this repo; do not use as locked same-board AUPRC unless rescored"
    detail = ctx[ctx["AUROC"].notna()][
        ["benchmark_scope", "method", "n_scored", "n_pos", "AUROC", "claim_boundary"]
    ].sort_values(["benchmark_scope", "AUROC"], ascending=[True, False])
    return summary, detail


def build_figures(master: pd.DataFrame, winloss: pd.DataFrame, claim: pd.DataFrame, public_summary: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    plt.style.use("default")

    algs = ["KG_GA_evolved", "CROSS_integrated", "CROSS_claimsafe", "CROSS_stress", "BigMHC_IM", "BigMHC_EL"]
    cols = ["all_AUPRC", "frozen_validation_AUPRC", "low_medium_leakage_AUPRC", "v6_smoke_AUPRC", "paperability_reward"]
    heat = master[master["algorithm"].isin(algs)].set_index("algorithm").reindex(algs)[cols]
    fig, ax = plt.subplots(figsize=(11.5, 5.2))
    im = ax.imshow(heat.to_numpy(dtype=float), aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(cols)), [c.replace("_", "\n") for c in cols], fontsize=9)
    ax.set_yticks(np.arange(len(heat.index)), heat.index)
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            v = heat.iloc[i, j]
            ax.text(j, i, "NA" if pd.isna(v) else f"{v:.3f}", ha="center", va="center", fontsize=8, color="#06101c" if pd.notna(v) and v > 0.65 else "white")
    ax.set_title("GA/RL champion board across benchmark contexts")
    fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    fig.tight_layout()
    fig.savefig(FIG / "fig1_ga_rl_champion_heatmap.png", dpi=220)
    fig.savefig(FIG / "fig1_ga_rl_champion_heatmap.pdf")
    plt.close(fig)

    deltas = winloss.set_index("comparator")["delta_all_AUPRC"].sort_values()
    fig, ax = plt.subplots(figsize=(10, 4.8))
    colors = ["#1f8a70" if v > 0 else "#b33f62" for v in deltas]
    ax.barh(deltas.index, deltas.values, color=colors)
    ax.axvline(0, color="#222", lw=1)
    ax.set_xlabel("KG_GA_evolved AUPRC minus comparator AUPRC")
    ax.set_title("Retrospective all-board AUPRC lift")
    fig.tight_layout()
    fig.savefig(FIG / "fig2_ga_rl_auprc_lift_vs_comparators.png", dpi=220)
    fig.savefig(FIG / "fig2_ga_rl_auprc_lift_vs_comparators.pdf")
    plt.close(fig)

    prod = master[master["claim_track"].eq("strict no-overlap production")].dropna(subset=["neo_strict_AUPRC"]).head(10)
    if not prod.empty:
        fig, ax = plt.subplots(figsize=(11, 5.2))
        y = np.arange(len(prod))
        ax.barh(y, prod["neo_strict_AUPRC"], color=["#b33f62" if "KG_GA" in a else "#285c7a" for a in prod["algorithm"]])
        ax.set_yticks(y, prod["algorithm"])
        ax.invert_yaxis()
        ax.set_xlim(0, max(1.0, float(prod["neo_strict_AUPRC"].max()) * 1.05))
        ax.set_xlabel("Strict no-overlap production AUPRC")
        ax.set_title("Context caveat: BAR-Neo still wins strict no-overlap reliability")
        fig.tight_layout()
        fig.savefig(FIG / "fig3_strict_no_overlap_context.png", dpi=220)
        fig.savefig(FIG / "fig3_strict_no_overlap_context.pdf")
        plt.close(fig)

    claim_colors = {"YES": 1.0, "NO": 0.0, "NO / context-specific": 0.35}
    claim_plot = claim.copy()
    claim_plot["score"] = claim_plot["answer"].map(claim_colors).fillna(0.5)
    fig, ax = plt.subplots(figsize=(11.5, 4.8))
    ax.imshow(claim_plot[["score"]].to_numpy(), aspect="auto", cmap="RdYlGn", vmin=0, vmax=1)
    ax.set_xticks([0], ["Allowed?"])
    ax.set_yticks(np.arange(len(claim_plot)), [f"C{i+1}" for i in range(len(claim_plot))])
    for i, r in claim_plot.iterrows():
        ax.text(0, i, r["answer"], ha="center", va="center", weight="bold")
    ax.set_title("Claim boundary board")
    fig.tight_layout()
    fig.savefig(FIG / "fig4_claim_boundary_board.png", dpi=220)
    fig.savefig(FIG / "fig4_claim_boundary_board.pdf")
    plt.close(fig)

    if not public_summary.empty:
        plot = public_summary.dropna(subset=["best_context_AUROC"]).sort_values("best_context_AUROC").copy()
        fig, ax = plt.subplots(figsize=(10.5, 5.4))
        colors = ["#b33f62" if m == "BigMHC_IM" else "#285c7a" for m in plot["method"]]
        ax.barh(plot["method"], plot["best_context_AUROC"], color=colors)
        ax.set_xlim(0, 1.0)
        ax.set_xlabel("Best context AUROC in existing public benchmark table")
        ax.set_title("Broader public comparator context: BigMHC is not the only baseline")
        for i, (_, r) in enumerate(plot.iterrows()):
            ax.text(float(r["best_context_AUROC"]) + 0.01, i, f"{float(r['best_context_AUROC']):.3f}", va="center", fontsize=8)
        fig.tight_layout()
        fig.savefig(FIG / "fig5_public_competitor_context.png", dpi=220)
        fig.savefig(FIG / "fig5_public_competitor_context.pdf")
        plt.close(fig)


def write_report(
    master: pd.DataFrame,
    winloss: pd.DataFrame,
    claim: pd.DataFrame,
    ablation: pd.DataFrame,
    public_summary: pd.DataFrame,
    public_detail: pd.DataFrame,
) -> None:
    champion = master[master["algorithm"].eq("KG_GA_evolved")].iloc[0]
    next_all = master[master["algorithm"].isin(["CROSS_integrated", "CROSS_claimsafe", "CROSS_stress", "BigMHC_IM", "BigMHC_EL"])].sort_values("all_AUPRC", ascending=False).iloc[0]
    prod = master[master["claim_track"].eq("strict no-overlap production")].dropna(subset=["neo_strict_AUPRC"]).sort_values("neo_strict_AUPRC", ascending=False).head(1)
    lines = [
        "# CROSS-Neo GA/RL Full Algorithm Comparison",
        "",
        "## 결론",
        "",
        f"- KG_GA_evolved is the current champion on the 2,715-row CROSS-Neo board: AUPRC {champion['all_AUPRC']:.3f}, AUROC {champion['all_AUROC']:.3f}, top96 {int(champion['all_top96_hits'])}/96.",
        f"- It beats the strongest fixed CROSS-Neo comparator on all-board AUPRC: {next_all['algorithm']} {next_all['all_AUPRC']:.3f}, delta {champion['all_AUPRC'] - next_all['all_AUPRC']:.3f}.",
        f"- It is also strongest on frozen validation-like sources: AUPRC {champion['frozen_validation_AUPRC']:.3f}, AUROC {champion['frozen_validation_AUROC']:.3f}, top24 {int(champion['frozen_top24_hits'])}/24.",
        f"- It remains strongest on low/medium leakage: AUPRC {champion['low_medium_leakage_AUPRC']:.3f}, top96 {int(champion['low_medium_top96_hits'])}/96.",
        "- Honest caveat: v6 smoke board AUPRC winner is BigMHC_EL, not KG_GA. KG_GA is strong there but not the winner.",
    ]
    if not prod.empty:
        pr = prod.iloc[0]
        lines.append(
            f"- Strict no-overlap production reliability is separate: {pr['algorithm']} is best there with AUPRC {pr['neo_strict_AUPRC']:.3f}; KG_GA should be positioned as integrated discovery/experiment-priority, not the clean reliability winner."
        )
    lines += [
        "",
        "## Master Board",
        "",
        markdown_table(master[
            [
                "algorithm",
                "algorithm_family",
                "claim_track",
                "all_AUPRC",
                "frozen_validation_AUPRC",
                "low_medium_leakage_AUPRC",
                "v6_smoke_AUPRC",
                "paperability_reward",
                "neo_strict_AUPRC",
                "primary_disposition",
            ]
        ].head(25)),
        "",
        "## Win/Loss: KG_GA_evolved vs Core Comparators",
        "",
        markdown_table(winloss),
        "",
        "## Broader Public Comparator Context",
        "",
        "BigMHC is not the only public/legacy comparator. Existing repo context includes RF_biophys, MHCflurry, PRIME, TransPHLA, DeepImmuno, NetMHCpan_4.1, VQC, GP_quantum, Structure_LR, and ESM2_Bayesian. These rows are AUROC-only context with possible overlap/training-status caveats, not locked same-board AUPRC.",
        "",
        markdown_table(public_summary.head(20)) if not public_summary.empty else "NA",
        "",
        "### Public Context by Dataset",
        "",
        markdown_table(public_detail.head(40)) if not public_detail.empty else "NA",
        "",
        "## Claim Boundary Matrix",
        "",
        markdown_table(claim),
        "",
        "## Ablation / Remaining Gap Matrix",
        "",
        markdown_table(ablation),
        "",
        "## Claim-Safe Sentence",
        "",
        "The current result supports a claim that the KG-guided GA controller is the strongest retrospective and validation-like CROSS-Neo prioritization controller among tested local comparators, while the RL component is presently a method-search/paperability blueprint and prospective immunogenicity proof still requires locked mutant-vs-WT/decoy wetlab readout.",
        "",
    ]
    (OUT / "GA_RL_FULL_COMPARISON_REPORT_KR.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_html(
    master: pd.DataFrame,
    winloss: pd.DataFrame,
    claim: pd.DataFrame,
    ablation: pd.DataFrame,
    public_summary: pd.DataFrame,
    public_detail: pd.DataFrame,
) -> str:
    champion = master[master["algorithm"].eq("KG_GA_evolved")].iloc[0]
    strict_prod = master[master["claim_track"].eq("strict no-overlap production")].dropna(subset=["neo_strict_AUPRC"]).sort_values("neo_strict_AUPRC", ascending=False).head(1)
    strict_txt = "NA"
    if not strict_prod.empty:
        strict_txt = f"{strict_prod.iloc[0]['algorithm']} {strict_prod.iloc[0]['neo_strict_AUPRC']:.3f}"
    cards = ""
    for fig, title, note in [
        ("fig1_ga_rl_champion_heatmap.png", "Champion Heatmap", "All-board, frozen validation-like, low-leakage, v6, paperability."),
        ("fig2_ga_rl_auprc_lift_vs_comparators.png", "AUPRC Lift", "KG-GA delta against core comparators."),
        ("fig3_strict_no_overlap_context.png", "Strict No-Overlap Context", "Where BAR-Neo reliability still wins."),
        ("fig4_claim_boundary_board.png", "Claim Boundary", "Allowed and blocked claims."),
        ("fig5_public_competitor_context.png", "Public Comparator Context", "MHCflurry, PRIME, TransPHLA, NetMHCpan, RF_biophys and others."),
    ]:
        if not (ASSET / fig).exists() and not (FIG / fig).exists():
            continue
        cards += f"""
        <article class="card">
          <img src="assets/cross_neo_ga_rl_full_comparison/{fig}" alt="{esc(title)}">
          <h3>{esc(title)}</h3>
          <p>{esc(note)}</p>
          <a href="assets/cross_neo_ga_rl_full_comparison/{fig}">PNG</a>
        </article>
        """
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>CROSS-Neo GA/RL Full Comparison</title>
<style>
:root {{ --bg:#07111f; --panel:#101d2d; --ink:#edf5ff; --muted:#9fb0c5; --gold:#f2c46d; --line:#27384f; --red:#ff7d8d; --green:#57c785; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font-family:Inter,Arial,sans-serif; line-height:1.5; }}
header {{ padding:48px 5vw 28px; border-bottom:1px solid var(--line); background:#0b1626; }}
.kicker {{ color:var(--gold); text-transform:uppercase; letter-spacing:.08em; font-size:12px; }}
h1 {{ font-size:clamp(30px,5vw,62px); margin:8px 0 10px; line-height:1.0; }}
.lead {{ color:var(--muted); max-width:1100px; font-size:18px; }}
.stats {{ display:grid; grid-template-columns:repeat(6,minmax(130px,1fr)); gap:12px; margin-top:24px; }}
.stat {{ background:var(--panel); border:1px solid var(--line); padding:14px; border-radius:8px; }}
.stat b {{ display:block; font-size:26px; color:var(--gold); }}
.stat span {{ color:var(--muted); font-size:12px; }}
main {{ padding:28px 5vw 60px; max-width:1500px; margin:0 auto; }}
section {{ margin:34px 0; }}
h2 {{ font-size:24px; margin-bottom:12px; color:#fff; }}
.note {{ border-left:4px solid var(--gold); background:#101d2d; padding:14px 16px; color:var(--muted); }}
.grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(320px,1fr)); gap:16px; }}
.card {{ border:1px solid var(--line); background:var(--panel); border-radius:8px; padding:12px; }}
.card img {{ width:100%; background:#fff; border-radius:6px; }}
.card h3 {{ margin:10px 0 4px; }}
.card p {{ color:var(--muted); }}
a {{ color:var(--gold); }}
table {{ width:100%; border-collapse:collapse; font-size:13px; background:var(--panel); }}
th, td {{ border:1px solid var(--line); padding:8px; vertical-align:top; }}
th {{ color:var(--gold); text-align:left; position:sticky; top:0; background:#101d2d; }}
.tablewrap {{ overflow:auto; border:1px solid var(--line); }}
.yes {{ color:var(--green); font-weight:700; }}
.no {{ color:var(--red); font-weight:700; }}
code {{ color:var(--gold); }}
@media (max-width:900px) {{ .stats {{ grid-template-columns:repeat(2,1fr); }} }}
</style></head>
<body><header>
<div class="kicker">CROSS-Neo · GA/RL · all-comparator audit</div>
<h1>KG-GA is the current champion, with a clean boundary.</h1>
<p class="lead">This page compares the evolved knowledge-graph genetic controller against BigMHC, fixed CROSS-Neo scores, claim-safe fallback, strict no-overlap production rankers, and paperability-aware Darwin-RL outputs.</p>
<div class="stats">
<div class="stat"><b>{champion['all_AUPRC']:.3f}</b><span>KG-GA all AUPRC</span></div>
<div class="stat"><b>{champion['all_AUROC']:.3f}</b><span>KG-GA all AUROC</span></div>
<div class="stat"><b>{int(champion['all_top96_hits'])}/96</b><span>all-board top96 hits</span></div>
<div class="stat"><b>{champion['frozen_validation_AUPRC']:.3f}</b><span>frozen validation-like AUPRC</span></div>
<div class="stat"><b>{champion['low_medium_leakage_AUPRC']:.3f}</b><span>low/medium leakage AUPRC</span></div>
<div class="stat"><b>{strict_txt}</b><span>strict no-overlap winner</span></div>
</div>
</header><main>
<section><h2>Decision</h2>
<p class="note"><b>Use KG_GA_evolved as the champion discovery/experiment-priority algorithm.</b> Use CROSS_claimsafe as the conservative claim-safe fallback. Use BAR-Neo_confidence for strict no-overlap production reliability. Do not claim prospective immunogenicity until locked mutant-vs-WT/decoy wetlab readout passes.</p></section>
<section><h2>Figure Board</h2><div class="grid">{cards}</div></section>
<section><h2>Master Comparison</h2><div class="tablewrap">{table_html(master, ['algorithm','algorithm_family','claim_track','all_AUPRC','all_AUROC','all_top96_hits','frozen_validation_AUPRC','low_medium_leakage_AUPRC','v6_smoke_AUPRC','paperability_reward','neo_strict_AUPRC','headline_score','primary_disposition'], 40)}</div></section>
<section><h2>Win/Loss</h2><div class="tablewrap">{table_html(winloss, ['comparator','delta_all_AUPRC','delta_frozen_validation_AUPRC','delta_low_medium_leakage_AUPRC','delta_v6_smoke_AUPRC','delta_paperability_reward','wins','tested_metrics','win_fraction','disposition'], 20)}</div></section>
<section><h2>Broader Public Comparator Context</h2>
<p class="note">BigMHC is only one public comparator. Existing repo context also includes RF_biophys, MHCflurry, PRIME, TransPHLA, DeepImmuno, NetMHCpan_4.1, VQC, GP_quantum, Structure_LR, and ESM2_Bayesian. These rows are context-only AUROC benchmarks with possible overlap/training-status caveats, not locked same-board AUPRC.</p>
<div class="tablewrap">{table_html(public_summary, ['method','n_context_rows','total_scored','datasets_with_auroc','mean_context_AUROC','best_context_AUROC','min_context_AUROC','claim_use'], 30)}</div>
<h2>Public Context by Dataset</h2>
<div class="tablewrap">{table_html(public_detail, ['benchmark_scope','method','n_scored','n_pos','AUROC','claim_boundary'], 50)}</div></section>
<section><h2>Claim Boundary</h2><div class="tablewrap">{table_html(claim, ['claim_item','answer','evidence','boundary'], 20)}</div></section>
<section><h2>Ablation Gaps</h2><div class="tablewrap">{table_html(ablation, ['module','current_status','evidence','gap','priority'], 20)}</div></section>
<section><h2>Sources</h2>
<ul>
<li><code>{esc(str(KG / 'kg_ga_benchmark.tsv'))}</code></li>
<li><code>{esc(str(KG / 'kg_ga_v6_smoke_wetlab_comparison.tsv'))}</code></li>
<li><code>{esc(str(PAPER / 'paperability_algorithm_ranking.tsv'))}</code></li>
<li><code>{esc(str(NEO / 'metrics/leaderboard_production_track_strict_no_overlap.tsv'))}</code></li>
<li><code>{esc(str(NEO / 'metrics/ga_rl_algorithm_metrics.tsv'))}</code></li>
<li><code>{esc(str(PUBLIC_CTX))}</code></li>
</ul>
</section>
</main></body></html>"""


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    ASSET.mkdir(parents=True, exist_ok=True)
    LIVE_ASSET.mkdir(parents=True, exist_ok=True)
    master = build_master_table()
    winloss = build_winloss(master)
    claim = build_claim_matrix(master)
    ablation = build_ablation_gap()
    public_summary, public_detail = build_public_competitor_context()
    master.to_csv(OUT / "ga_rl_master_algorithm_comparison.tsv", sep="\t", index=False)
    winloss.to_csv(OUT / "ga_rl_winloss_vs_comparators.tsv", sep="\t", index=False)
    claim.to_csv(OUT / "ga_rl_claim_boundary_matrix.tsv", sep="\t", index=False)
    ablation.to_csv(OUT / "ga_rl_ablation_gap_matrix.tsv", sep="\t", index=False)
    public_summary.to_csv(OUT / "ga_rl_public_competitor_context_summary.tsv", sep="\t", index=False)
    public_detail.to_csv(OUT / "ga_rl_public_competitor_context_by_dataset.tsv", sep="\t", index=False)
    write_report(master, winloss, claim, ablation, public_summary, public_detail)
    build_figures(master, winloss, claim, public_summary)
    for src in FIG.glob("*"):
        if src.is_file():
            safe_copy(src, ASSET / src.name)
            safe_copy(src, LIVE_ASSET / src.name)
    for src in [
        OUT / "ga_rl_master_algorithm_comparison.tsv",
        OUT / "ga_rl_winloss_vs_comparators.tsv",
        OUT / "ga_rl_claim_boundary_matrix.tsv",
        OUT / "ga_rl_ablation_gap_matrix.tsv",
        OUT / "ga_rl_public_competitor_context_summary.tsv",
        OUT / "ga_rl_public_competitor_context_by_dataset.tsv",
        OUT / "GA_RL_FULL_COMPARISON_REPORT_KR.md",
    ]:
        safe_copy(src, ASSET / src.name)
        safe_copy(src, LIVE_ASSET / src.name)
    PAGE.write_text(build_html(master, winloss, claim, ablation, public_summary, public_detail), encoding="utf-8")
    safe_copy(PAGE, LIVE_PAGE)
    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "out_dir": str(OUT),
        "html_path": str(PAGE),
        "live_html_path": str(LIVE_PAGE),
        "n_master_rows": int(len(master)),
        "n_winloss_rows": int(len(winloss)),
        "n_public_context_methods": int(len(public_summary)),
        "top_master_algorithm": master.iloc[0]["algorithm"],
        "kg_ga_all_auprc": float(master[master["algorithm"].eq("KG_GA_evolved")].iloc[0]["all_AUPRC"]),
        "kg_ga_frozen_validation_auprc": float(master[master["algorithm"].eq("KG_GA_evolved")].iloc[0]["frozen_validation_AUPRC"]),
        "kg_ga_low_medium_leakage_auprc": float(master[master["algorithm"].eq("KG_GA_evolved")].iloc[0]["low_medium_leakage_AUPRC"]),
        "boundary": "retrospective and validation-like prioritization comparison, not prospective immunogenicity proof",
    }
    (OUT / "ga_rl_full_comparison_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
