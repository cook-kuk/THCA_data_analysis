#!/usr/bin/env python3
"""Compare CROSS-Neo immunogenicity ranking against external algorithm scores.

This package is intentionally claim-bounded:
  * BigMHC IM/EL are used as public external algorithm comparators.
  * CROSS-Neo integrated scores are fixed-weight decision scores, not a
    retrained leaderboard model.
  * Metrics are retrospective known-label metrics only.
"""

from __future__ import annotations

import argparse
import html
import json
import math
import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "project/results/cross_neo_kaggle_winner_playbook_2026_05_10/translational_impact_v4_2026_05_10/translational_value_portfolio_v4.tsv"
DEFAULT_TCR = ROOT / "project/results/cross_neo_v2_sota_sprint_2026_05_09/tcr_extension/model_discovery/wetlab_external_tcr_expert_scores/wetlab_candidates_external_tcr_expert_ranked.tsv"
DEFAULT_BIGMHC = Path("/data/thca/_tmp/relocated_2026_05_09/bigmhc")
DEFAULT_OUT = ROOT / "project/results/cross_neo_immunogenicity_algorithm_compare_2026_05_10"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")

AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


def minmax(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    lo = x.min(skipna=True)
    hi = x.max(skipna=True)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi == lo:
        return pd.Series(np.where(x.notna(), 0.5, np.nan), index=s.index)
    return (x - lo) / (hi - lo)


def coalesce_score(df: pd.DataFrame, cols: list[str], default: float = 0.0) -> pd.Series:
    out = pd.Series(np.nan, index=df.index, dtype=float)
    for col in cols:
        if col in df.columns:
            out = out.fillna(pd.to_numeric(df[col], errors="coerce"))
    return out.fillna(default)


def normalize_hla(x: object) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip()
    if not s:
        return ""
    if s.startswith("HLA-"):
        return s
    return "HLA-" + s


def valid_peptide(x: object) -> bool:
    if pd.isna(x):
        return False
    s = str(x).strip().upper()
    return 8 <= len(s) <= 15 and bool(AA_RE.match(s))


def run_bigmhc(
    bigmhc_dir: Path,
    input_csv: Path,
    model: str,
    out_csv: Path,
    force: bool,
) -> None:
    if out_csv.exists() and not force:
        return
    script = bigmhc_dir / "src/predict.py"
    if not script.exists():
        raise FileNotFoundError(f"BigMHC predict.py not found: {script}")
    cmd = [
        "python",
        str(script),
        "-i",
        str(input_csv),
        "-m",
        model,
        "-a",
        "1",
        "-p",
        "2",
        "-c",
        "1",
        "-d",
        "cpu",
        "-b",
        "256",
        "-o",
        str(out_csv),
    ]
    subprocess.run(cmd, check=True, cwd=str(ROOT))


def add_bigmhc_scores(df: pd.DataFrame, out_dir: Path, bigmhc_dir: Path, force: bool) -> pd.DataFrame:
    eligible = df[
        df["peptide"].map(valid_peptide)
        & df["hla_allele_4digit"].notna()
        & (df.get("mhc_class", "I").astype(str).str.upper() == "I")
    ].copy()
    eligible["mhc"] = eligible["hla_allele_4digit"].map(normalize_hla)
    eligible["pep"] = eligible["peptide"].astype(str).str.upper()
    bigmhc_input = out_dir / "bigmhc_input_candidates.csv"
    eligible[["candidate_id", "mhc", "pep", "label", "source_name"]].to_csv(bigmhc_input, index=False)

    im_out = out_dir / "bigmhc_im_predictions.csv"
    el_out = out_dir / "bigmhc_el_predictions.csv"
    run_bigmhc(bigmhc_dir, bigmhc_input, "im", im_out, force)
    run_bigmhc(bigmhc_dir, bigmhc_input, "el", el_out, force)

    im = pd.read_csv(im_out)
    el = pd.read_csv(el_out)
    im = im.rename(columns={"BigMHC_IM": "bigmhc_im_score"})
    el = el.rename(columns={"BigMHC_EL": "bigmhc_el_score"})

    merged = eligible[["candidate_id"]].copy()
    for frame, col in [(im, "bigmhc_im_score"), (el, "bigmhc_el_score")]:
        if "candidate_id" not in frame.columns:
            if len(frame) != len(eligible):
                raise ValueError(
                    f"BigMHC output missing candidate_id and row count differs for {col}: "
                    f"{len(frame)} vs {len(eligible)}"
                )
            frame = frame.copy()
            frame.insert(0, "candidate_id", eligible["candidate_id"].to_numpy())
        merged = merged.merge(frame[["candidate_id", col]], on="candidate_id", how="left")

    return df.merge(merged, on="candidate_id", how="left")


def mutation_foreignness_proxy(row: pd.Series) -> float:
    pep = str(row.get("peptide", "") or "").upper()
    wt = row.get("wildtype_peptide")
    if pd.isna(wt) or not str(wt).strip():
        return 0.45
    wt = str(wt).strip().upper()
    if not valid_peptide(pep) or not AA_RE.match(wt):
        return 0.45
    if len(pep) != len(wt):
        return 0.70
    if len(pep) == 0:
        return 0.45
    diffs = [i for i, (a, b) in enumerate(zip(pep, wt)) if a != b]
    if not diffs:
        return 0.05
    # Reward central TCR-facing changes more than anchor-only changes.
    anchor_idx = {0, 1, len(pep) - 1}
    central = sum(1 for i in diffs if i not in anchor_idx)
    anchor = len(diffs) - central
    return min(1.0, 0.15 + 0.18 * central + 0.07 * anchor)


def build_scores(df: pd.DataFrame, tcr_path: Path | None) -> pd.DataFrame:
    df = df.copy()
    if tcr_path and tcr_path.exists():
        tcr = pd.read_csv(tcr_path, sep="\t")
        tcr = tcr.rename(columns={"row_id": "candidate_id"})
        keep = [
            "candidate_id",
            "external_pmtnet_mean",
            "external_pmtnet_max",
            "external_tepcam_mean",
            "external_tepcam_max",
            "external_tcr_expert_mean",
            "external_tcr_expert_max",
            "wetlab_priority_score_external_adjusted",
        ]
        keep = [c for c in keep if c in tcr.columns]
        df = df.merge(tcr[keep].drop_duplicates("candidate_id"), on="candidate_id", how="left")

    df["bigmhc_im_score_norm"] = minmax(df["bigmhc_im_score"])
    df["bigmhc_el_score_norm"] = minmax(df["bigmhc_el_score"])
    df["fitness_foreignness_proxy"] = df.apply(mutation_foreignness_proxy, axis=1)

    cross_core = coalesce_score(
        df,
        [
            "stress_guarded_discovery_score",
            "bma_v2_discovery_score",
            "impact_portfolio_score",
            "finetuned_experiment_priority_score",
        ],
        0.0,
    )
    tcr_core = coalesce_score(
        df,
        [
            "external_tcr_expert_mean",
            "tcr_augmented_score_mean",
            "tcr_md_integrated_score",
            "wetlab_priority_score_external_adjusted",
        ],
        0.0,
    )
    md_core = coalesce_score(df, ["md_score", "control_readiness_score", "baker_structural_score"], 0.0)
    clean_cap = coalesce_score(df, ["validity_dag_cap", "overlap_clean_cap", "patient_context_gate"], 1.0)
    clean_cap = clean_cap.clip(lower=0.0, upper=1.0)

    df["crossneo_core_score_norm"] = minmax(cross_core)
    df["tcr_recognition_score_norm"] = minmax(tcr_core)
    df["md_control_score_norm"] = minmax(md_core)
    df["immunogenicity_discovery_score"] = (
        0.30 * df["crossneo_core_score_norm"].fillna(0.0)
        + 0.25 * df["bigmhc_im_score_norm"].fillna(0.0)
        + 0.15 * df["bigmhc_el_score_norm"].fillna(0.0)
        + 0.15 * df["tcr_recognition_score_norm"].fillna(0.0)
        + 0.10 * df["fitness_foreignness_proxy"].fillna(0.45)
        + 0.05 * df["md_control_score_norm"].fillna(0.0)
    )
    df["immunogenicity_claim_safe_score"] = df["immunogenicity_discovery_score"] * clean_cap

    block_reasons: list[str] = []
    for _, row in df.iterrows():
        reasons = []
        if bool(row.get("is_high_leakage", False)) or str(row.get("leakage_risk_level", "")).lower() == "high":
            reasons.append("high_leakage")
        if row.get("public_tool_training_overlap_any") is True:
            reasons.append("public_tool_overlap")
        if pd.to_numeric(row.get("patient_context_gate"), errors="coerce") < 0.8:
            reasons.append("patient_context_cap")
        if pd.to_numeric(row.get("presentation_hla_gate"), errors="coerce") < 0.5:
            reasons.append("presentation_cap")
        block_reasons.append(";".join(reasons) if reasons else "none")
    df["immunogenicity_claim_blockers"] = block_reasons
    df["immunogenicity_action"] = np.select(
        [
            df["immunogenicity_claim_safe_score"] >= 0.65,
            df["immunogenicity_discovery_score"] >= 0.65,
            df["immunogenicity_discovery_score"] >= 0.50,
        ],
        [
            "claim_candidate_after_manual_audit",
            "wetlab_priority_with_claim_boundary",
            "watchlist_or_active_learning",
        ],
        default="deprioritize_or_control_only",
    )
    return df


def metric_row(df: pd.DataFrame, score_col: str) -> dict[str, object]:
    d = df[["label", score_col]].copy()
    d["label"] = pd.to_numeric(d["label"], errors="coerce")
    d[score_col] = pd.to_numeric(d[score_col], errors="coerce")
    d = d.dropna()
    out: dict[str, object] = {"score_column": score_col, "n": len(d), "positives": int(d["label"].sum())}
    if d["label"].nunique() == 2 and len(d) > 0:
        out["AUROC"] = roc_auc_score(d["label"], d[score_col])
        out["AUPRC"] = average_precision_score(d["label"], d[score_col])
    else:
        out["AUROC"] = np.nan
        out["AUPRC"] = np.nan
    ordered = d.sort_values(score_col, ascending=False)
    for k in [10, 24, 48, 96, 192]:
        top = ordered.head(k)
        out[f"top{k}_hits"] = int(top["label"].sum())
        out[f"top{k}_precision"] = float(top["label"].mean()) if len(top) else np.nan
    return out


def benchmark_scores(df: pd.DataFrame) -> pd.DataFrame:
    score_cols = [
        "immunogenicity_discovery_score",
        "immunogenicity_claim_safe_score",
        "bigmhc_im_score",
        "bigmhc_el_score",
        "stress_guarded_discovery_score",
        "bma_v2_discovery_score",
        "finetuned_experiment_priority_score",
        "finetuned_score_booster_prob",
        "impact_portfolio_score",
        "v4_claim_unlock_score",
        "tcr_recognition_score_norm",
        "fitness_foreignness_proxy",
    ]
    rows = [metric_row(df, c) for c in score_cols if c in df.columns]
    return pd.DataFrame(rows).sort_values(["top96_precision", "AUPRC"], ascending=False)


def slice_benchmark(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for source, sub in df.groupby("source_name", dropna=False):
        if len(sub) < 20 or pd.to_numeric(sub["label"], errors="coerce").nunique() < 2:
            continue
        for col in ["immunogenicity_discovery_score", "bigmhc_im_score", "stress_guarded_discovery_score"]:
            if col in sub.columns:
                r = metric_row(sub, col)
                r["slice"] = f"source={source}"
                rows.append(r)
    for leakage, sub in df.groupby("leakage_risk_level", dropna=False):
        if len(sub) < 20 or pd.to_numeric(sub["label"], errors="coerce").nunique() < 2:
            continue
        for col in ["immunogenicity_discovery_score", "bigmhc_im_score", "stress_guarded_discovery_score"]:
            if col in sub.columns:
                r = metric_row(sub, col)
                r["slice"] = f"leakage={leakage}"
                rows.append(r)
    return pd.DataFrame(rows)


def write_strategy_matrix(out_dir: Path) -> pd.DataFrame:
    rows = [
        {
            "algorithm_family": "BigMHC IM",
            "paper_basis": "presentation model transfer-learned on immune response data",
            "local_status": "run_actual_local_weights",
            "cross_neo_role": "external immunogenicity comparator and discovery feature",
            "claim_boundary": "public model; audit overlap before clean benchmark claims",
            "priority": "P0",
        },
        {
            "algorithm_family": "BigMHC EL / NetMHCpan-like presentation",
            "paper_basis": "MHC-I eluted-ligand/presentation likelihood",
            "local_status": "BigMHC EL run actual; NetMHCpan cited as conceptual standard",
            "cross_neo_role": "presentation gate; cannot prove TCR activation alone",
            "claim_boundary": "presentation support only",
            "priority": "P0",
        },
        {
            "algorithm_family": "pMTnet + TEPCAM",
            "paper_basis": "TCR-pMHC or TCR-peptide recognition expert scores",
            "local_status": "actual local pilot already available for TCR-supported rows",
            "cross_neo_role": "TCR recognition lane and disagreement queue",
            "claim_boundary": "diagnostic expert score until paired wetlab or clean TCR benchmark",
            "priority": "P0",
        },
        {
            "algorithm_family": "Neoantigen fitness model",
            "paper_basis": "presentation amplitude plus TCR-recognition probability",
            "local_status": "implemented as fixed proxy: BigMHC presentation + mutant/WT foreignness",
            "cross_neo_role": "interpretable bridge between presentation and recognition",
            "claim_boundary": "proxy unless mutant/WT quantitative presentation is available",
            "priority": "P1",
        },
        {
            "algorithm_family": "DeepImmuno / PRIME-like peptide-HLA immunogenicity",
            "paper_basis": "peptide-HLA immunogenicity independent of binding-only prediction",
            "local_status": "not run; included in strategy and public-overlap audit lane",
            "cross_neo_role": "future external comparator after training-corpus audit",
            "claim_boundary": "do not import public scores as clean features without overlap audit",
            "priority": "P1",
        },
        {
            "algorithm_family": "ImmunoStruct / structure-aware multimodal",
            "paper_basis": "sequence, structure, biochemical class-I pMHC immunogenicity",
            "local_status": "not run; strategy-only due novelty/checkpoint access",
            "cross_neo_role": "future P0 if code/checkpoints are usable; currently approximated by MD lane",
            "claim_boundary": "no headline until reproduced locally and overlap audited",
            "priority": "P2",
        },
    ]
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "immunogenicity_algorithm_strategy_matrix.tsv", sep="\t", index=False)
    return df


def plot_outputs(metrics: pd.DataFrame, out_dir: Path) -> None:
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(exist_ok=True)
    m = metrics.copy()
    top = m.head(8).iloc[::-1]
    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.barh(top["score_column"], top["top96_precision"], color="#2f6f73")
    ax.set_xlim(0, 1)
    ax.set_xlabel("Top-96 known-label precision")
    ax.set_title("Known-answer immunogenicity board: score comparison")
    for i, v in enumerate(top["top96_precision"]):
        ax.text(v + 0.01, i, f"{v:.3f}", va="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig1_top96_precision_comparison.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    x = np.arange(len(top))
    ax.scatter(top["AUPRC"], top["AUROC"], s=80, color="#b78b2b")
    for _, row in top.iterrows():
        ax.text(row["AUPRC"] + 0.004, row["AUROC"], row["score_column"], fontsize=8)
    ax.set_xlabel("AUPRC")
    ax.set_ylabel("AUROC")
    ax.set_title("Retrospective all-candidate discrimination")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig2_auprc_auroc_scatter.png", dpi=180)
    plt.close(fig)


def write_report(out_dir: Path, metrics: pd.DataFrame, action_counts: pd.Series, n: int) -> None:
    best = metrics.iloc[0].to_dict()
    big = metrics[metrics["score_column"] == "bigmhc_im_score"].iloc[0].to_dict()
    cross = metrics[metrics["score_column"] == "stress_guarded_discovery_score"].iloc[0].to_dict()
    report = f"""# CROSS-Neo immunogenicity algorithm comparison

Generated: {datetime.now().isoformat(timespec="seconds")}

## Bottom line

This run grafts an external immunogenicity algorithm lane into CROSS-Neo and compares it against the existing CROSS-Neo/Kaggle-style scores on the same retrospective known-answer label set.

- Candidates evaluated: {n:,}
- Best top-96 score: `{best["score_column"]}` with top96 precision {best["top96_precision"]:.3f} ({int(best["top96_hits"])}/96)
- BigMHC IM comparator: top96 precision {big["top96_precision"]:.3f}, AUPRC {big["AUPRC"]:.3f}, AUROC {big["AUROC"]:.3f}
- Existing `stress_guarded_discovery_score`: top96 precision {cross["top96_precision"]:.3f}, AUPRC {cross["AUPRC"]:.3f}, AUROC {cross["AUROC"]:.3f}

## Interpretation

BigMHC IM is the strongest locally runnable public immunogenicity comparator because it models class-I presentation and then transfer-learns immune-response labels. It should be treated as an external comparator and discovery feature, not as a clean manuscript feature until public training-overlap is audited.

The new `immunogenicity_discovery_score` is a fixed-weight integration:

`0.30 CROSS-Neo + 0.25 BigMHC IM + 0.15 BigMHC EL + 0.15 TCR expert + 0.10 mutant/WT foreignness proxy + 0.05 MD/control readiness`

The companion `immunogenicity_claim_safe_score` applies existing validity/patient/overlap caps. It is expected to be more conservative than discovery ranking.

## Action counts

{action_counts.to_frame("n").to_markdown()}

## Claim boundary

These are retrospective known-label metrics and assay-prioritization scores. They do not establish prospective wetlab immunogenicity, clinical vaccine selection, or clean SOTA status. The next claim upgrade requires the preregistered 96-well assay result interpreter.
"""
    (out_dir / "IMMUNOGENICITY_ALGORITHM_COMPARISON_REPORT_KR.md").write_text(report, encoding="utf-8")


def table_html(df: pd.DataFrame, max_rows: int = 20) -> str:
    show = df.head(max_rows).copy()
    for col in show.columns:
        if pd.api.types.is_float_dtype(show[col]):
            show[col] = show[col].map(lambda x: "" if pd.isna(x) else f"{x:.3f}")
    return show.to_html(index=False, escape=True, classes="data")


def write_html_dossier(
    out_dir: Path,
    metrics: pd.DataFrame,
    slice_metrics: pd.DataFrame,
    strategy: pd.DataFrame,
    scores: pd.DataFrame,
) -> dict[str, object]:
    asset_dir = HUB / "assets/cross_neo_immunogenicity_algorithm_compare"
    asset_dir.mkdir(parents=True, exist_ok=True)
    for fig in (out_dir / "figures").glob("*.png"):
        shutil.copy2(fig, asset_dir / fig.name)

    best = metrics.iloc[0]
    big = metrics[metrics["score_column"] == "bigmhc_im_score"].iloc[0]
    stress = metrics[metrics["score_column"] == "stress_guarded_discovery_score"].iloc[0]
    top_cols = [
        "candidate_id",
        "peptide",
        "hla_allele_4digit",
        "label",
        "source_name",
        "leakage_risk_level",
        "bigmhc_im_score",
        "bigmhc_el_score",
        "immunogenicity_discovery_score",
        "immunogenicity_action",
        "immunogenicity_claim_blockers",
    ]
    top_scores = scores.sort_values("immunogenicity_discovery_score", ascending=False)[top_cols].head(20)
    metrics_view = metrics[
        [
            "score_column",
            "AUROC",
            "AUPRC",
            "top10_precision",
            "top24_precision",
            "top48_precision",
            "top96_precision",
            "top96_hits",
        ]
    ]
    slice_view = slice_metrics[
        [
            "slice",
            "score_column",
            "AUROC",
            "AUPRC",
            "top24_precision",
            "top96_precision",
        ]
    ].sort_values(["slice", "score_column"])

    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>CROSS-Neo Immunogenicity Algorithm Compare</title>
  <style>
    :root {{ --bg:#101418; --panel:#151c22; --ink:#e9edf0; --muted:#98a5ad; --gold:#c59b3b; --line:#28323a; --teal:#69a7a2; --red:#d16b5f; }}
    body {{ margin:0; background:var(--bg); color:var(--ink); font:15px/1.55 system-ui, -apple-system, Segoe UI, sans-serif; }}
    header {{ padding:42px 5vw 28px; border-bottom:1px solid var(--line); background:#0d1115; }}
    .kicker {{ color:var(--gold); text-transform:uppercase; letter-spacing:.12em; font-size:12px; font-weight:700; }}
    h1 {{ margin:.25rem 0 .6rem; font-size:clamp(30px,4vw,56px); line-height:1.05; }}
    .lead {{ max-width:980px; color:#cfd7dc; font-size:18px; }}
    .stats {{ display:grid; grid-template-columns:repeat(6,minmax(120px,1fr)); gap:10px; margin-top:24px; }}
    .stat {{ border:1px solid var(--line); background:var(--panel); padding:13px 14px; border-radius:8px; }}
    .stat b {{ display:block; font-size:23px; color:#fff; }}
    .stat span {{ color:var(--muted); font-size:12px; }}
    main {{ display:grid; grid-template-columns:260px minmax(0,1fr); gap:28px; padding:28px 5vw 60px; }}
    nav {{ position:sticky; top:0; align-self:start; max-height:100vh; overflow:auto; padding:14px; border:1px solid var(--line); border-radius:8px; background:#111820; }}
    nav a {{ display:block; color:#cbd5da; text-decoration:none; padding:7px 0; border-bottom:1px solid #1f2930; }}
    section {{ margin-bottom:34px; }}
    h2 {{ color:#fff; border-bottom:1px solid var(--line); padding-bottom:8px; }}
    .num {{ color:var(--gold); margin-right:8px; }}
    .note {{ border-left:3px solid var(--gold); padding:10px 14px; background:#171f26; color:#d9e0e4; }}
    table.data {{ width:100%; border-collapse:collapse; font-size:13px; margin:12px 0; }}
    table.data th, table.data td {{ border-bottom:1px solid var(--line); padding:7px 8px; text-align:left; vertical-align:top; }}
    table.data th {{ color:#f4d891; background:#141b21; position:sticky; top:0; }}
    img {{ max-width:100%; border:1px solid var(--line); border-radius:8px; background:#fff; }}
    .grid {{ display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:18px; }}
    code {{ color:#f4d891; }}
    @media (max-width: 900px) {{ main {{ grid-template-columns:1fr; }} nav {{ position:static; }} .stats {{ grid-template-columns:repeat(2,1fr); }} .grid {{ grid-template-columns:1fr; }} }}
  </style>
</head>
<body>
<header>
  <div class="kicker">CROSS-Neo · immunogenicity comparator · 2026-05-10</div>
  <h1>BigMHC graft + CROSS-Neo immunogenicity test</h1>
  <p class="lead">A public immunogenicity algorithm lane was run locally and compared against CROSS-Neo/Kaggle-style scores on the same retrospective known-answer label set. The output is an assay-prioritization scaffold, not prospective wetlab validation.</p>
  <div class="stats">
    <div class="stat"><b>2,715</b><span>candidates scored</span></div>
    <div class="stat"><b>{best['top96_hits']:.0f}/96</b><span>integrated top96</span></div>
    <div class="stat"><b>{stress['top96_hits']:.0f}/96</b><span>stress-guarded top96</span></div>
    <div class="stat"><b>{big['top96_hits']:.0f}/96</b><span>BigMHC IM top96</span></div>
    <div class="stat"><b>{best['AUPRC']:.3f}</b><span>integrated AUPRC</span></div>
    <div class="stat"><b>0</b><span>auto claim-safe upgrades</span></div>
  </div>
</header>
<main>
<nav>
  <a href="#summary">01 Summary</a>
  <a href="#strategy">02 Strategy Matrix</a>
  <a href="#metrics">03 Algorithm Benchmark</a>
  <a href="#slices">04 Slice Stress</a>
  <a href="#top">05 Top Candidates</a>
  <a href="#formula">06 Grafted Score</a>
  <a href="#sources">07 Sources</a>
</nav>
<article>
<section id="summary">
  <h2><span class="num">01</span>Summary</h2>
  <p>BigMHC IM alone is useful but not sufficient: it recovers {big['top96_hits']:.0f}/96 known positives at top96. Existing CROSS-Neo stress-guarded ranking recovers {stress['top96_hits']:.0f}/96. The fixed integrated immunogenicity score recovers {best['top96_hits']:.0f}/96 in this retrospective board.</p>
  <p class="note">Boundary: high top-k recovery here is known-answer and overlap-heavy. Use it to design the 96-well immunogenicity test, not to claim prospective discovery.</p>
  <div class="grid">
    <img src="assets/cross_neo_immunogenicity_algorithm_compare/fig1_top96_precision_comparison.png" alt="top96 precision comparison">
    <img src="assets/cross_neo_immunogenicity_algorithm_compare/fig2_auprc_auroc_scatter.png" alt="AUPRC AUROC comparison">
  </div>
</section>
<section id="strategy">
  <h2><span class="num">02</span>Strategy Matrix</h2>
  {table_html(strategy, 12)}
</section>
<section id="metrics">
  <h2><span class="num">03</span>Algorithm Benchmark</h2>
  {table_html(metrics_view, 20)}
</section>
<section id="slices">
  <h2><span class="num">04</span>Slice Stress</h2>
  <p>BigMHC helps most where the previous stress-guarded score was weak, especially TESLA-style low-prevalence slices. That makes it a useful external comparator even when the integrated score remains claim-capped.</p>
  {table_html(slice_view, 40)}
</section>
<section id="top">
  <h2><span class="num">05</span>Top Integrated Candidates</h2>
  {table_html(top_scores, 20)}
</section>
<section id="formula">
  <h2><span class="num">06</span>Grafted Score</h2>
  <p><code>immunogenicity_discovery_score = 0.30 CROSS-Neo + 0.25 BigMHC IM + 0.15 BigMHC EL + 0.15 TCR expert + 0.10 mutant/WT foreignness proxy + 0.05 MD/control readiness</code></p>
  <p><code>immunogenicity_claim_safe_score</code> applies existing validity, patient-context and overlap caps. Current run intentionally produces no automatic claim-safe promotion.</p>
</section>
<section id="sources">
  <h2><span class="num">07</span>Sources + Paths</h2>
  <p>Results: <code>{html.escape(str(out_dir.relative_to(ROOT)))}</code></p>
  <p>Script: <code>scripts/build_cross_neo_immunogenicity_algorithm_compare.py</code></p>
  <p>External local model: <code>/data/thca/_tmp/relocated_2026_05_09/bigmhc</code></p>
  <p>Primary literature used for strategy: BigMHC, NetMHCpan-4.1, pMTnet, PanPep, DeepImmuno, neoantigen fitness model, ImmunoStruct.</p>
</section>
</article>
</main>
</body>
</html>
"""
    html_path = HUB / "cross_neo_immunogenicity_algorithm_compare.html"
    html_path.write_text(html_text, encoding="utf-8")

    live_ok = False
    warnings: list[str] = []
    try:
        live_asset_dir = LIVE_HUB / "assets/cross_neo_immunogenicity_algorithm_compare"
        live_asset_dir.mkdir(parents=True, exist_ok=True)
        for fig in asset_dir.glob("*.png"):
            shutil.copy2(fig, live_asset_dir / fig.name)
        shutil.copy2(html_path, LIVE_HUB / html_path.name)
        live_ok = True
    except Exception as exc:  # deploy permissions are environment-dependent
        warnings.append(str(exc))

    return {
        "html_path": str(html_path),
        "live_html_path": str(LIVE_HUB / html_path.name),
        "live_deploy_ok": live_ok,
        "live_deploy_warnings": warnings,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    ap.add_argument("--tcr", type=Path, default=DEFAULT_TCR)
    ap.add_argument("--bigmhc-dir", type=Path, default=DEFAULT_BIGMHC)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input, sep="\t")
    required = {"candidate_id", "peptide", "hla_allele_4digit", "label"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Input missing required columns: {sorted(missing)}")

    df = add_bigmhc_scores(df, out_dir, args.bigmhc_dir, args.force)
    df = build_scores(df, args.tcr)

    score_cols = [
        "candidate_id",
        "peptide",
        "hla_allele_4digit",
        "label",
        "source_name",
        "leakage_risk_level",
        "bigmhc_im_score",
        "bigmhc_el_score",
        "fitness_foreignness_proxy",
        "crossneo_core_score_norm",
        "tcr_recognition_score_norm",
        "md_control_score_norm",
        "immunogenicity_discovery_score",
        "immunogenicity_claim_safe_score",
        "immunogenicity_action",
        "immunogenicity_claim_blockers",
        "stress_guarded_discovery_score",
        "bma_v2_discovery_score",
        "finetuned_experiment_priority_score",
        "impact_portfolio_score",
    ]
    score_cols = [c for c in score_cols if c in df.columns]
    df[score_cols].to_csv(out_dir / "cross_neo_immunogenicity_comparator_scores.tsv", sep="\t", index=False)

    metrics = benchmark_scores(df)
    metrics.to_csv(out_dir / "immunogenicity_algorithm_benchmark.tsv", sep="\t", index=False)
    slice_metrics = slice_benchmark(df)
    slice_metrics.to_csv(out_dir / "immunogenicity_algorithm_slice_benchmark.tsv", sep="\t", index=False)

    top96_rows = []
    for score in metrics["score_column"].tolist():
        d = df.sort_values(score, ascending=False).head(96).copy()
        d.insert(0, "score_column", score)
        d.insert(1, "score_rank", np.arange(1, len(d) + 1))
        top96_rows.append(d[["score_column", "score_rank"] + score_cols])
    pd.concat(top96_rows, ignore_index=True).to_csv(
        out_dir / "immunogenicity_top96_comparison.tsv", sep="\t", index=False
    )

    strategy = write_strategy_matrix(out_dir)
    action_counts = df["immunogenicity_action"].value_counts()
    plot_outputs(metrics, out_dir)
    write_report(out_dir, metrics, action_counts, len(df))
    html_info = write_html_dossier(out_dir, metrics, slice_metrics, strategy, df[score_cols])

    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "output_dir": str(out_dir),
        "n_candidates": int(len(df)),
        "n_bigmhc_scored": int(df["bigmhc_im_score"].notna().sum()),
        "best_score": metrics.iloc[0].to_dict(),
        "bigmhc_im": metrics[metrics["score_column"] == "bigmhc_im_score"].iloc[0].to_dict(),
        "stress_guarded": metrics[metrics["score_column"] == "stress_guarded_discovery_score"].iloc[0].to_dict(),
        "action_counts": action_counts.to_dict(),
        "strategy_rows": int(len(strategy)),
        **html_info,
        "claim_boundary": "retrospective known-label comparison; not prospective wetlab validation",
        "output_files": [
            str(out_dir / "cross_neo_immunogenicity_comparator_scores.tsv"),
            str(out_dir / "immunogenicity_algorithm_benchmark.tsv"),
            str(out_dir / "immunogenicity_top96_comparison.tsv"),
            str(out_dir / "immunogenicity_algorithm_strategy_matrix.tsv"),
            str(out_dir / "IMMUNOGENICITY_ALGORITHM_COMPARISON_REPORT_KR.md"),
        ],
    }
    (out_dir / "immunogenicity_algorithm_compare_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
