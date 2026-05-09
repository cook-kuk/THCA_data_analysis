#!/usr/bin/env python3
"""Write a BAR-Neo-X claim-safe interpretability dossier page."""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path
from typing import Any

import pandas as pd


PAGE_NAME = "cancer_vaccine_barneo_x_interpretability_2026_05_09.html"
ASSET_DIR_NAME = "barneo_x"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root; accepted for pipeline compatibility")
    parser.add_argument("--output-root", required=True, help="CLEAN-NeoBench BAR-Neo output directory")
    parser.add_argument("--hub-root", default="project/papers_hub_2026_05_04", help="HTML hub directory")
    parser.add_argument("--page-name", default=PAGE_NAME)
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except Exception:
        return {}


def fmt_num(value: object, digits: int = 3) -> str:
    try:
        if pd.isna(value):
            return "NA"
        x = float(value)
        if abs(x) >= 1000:
            return f"{x:,.0f}"
        return f"{x:.{digits}f}"
    except Exception:
        return html.escape(str(value))


def safe(value: object) -> str:
    return html.escape("" if pd.isna(value) else str(value))


def stat(value: str, label: str, note: str = "") -> str:
    return (
        "<div class=\"stat\">"
        f"<div class=\"stat-value\">{html.escape(value)}</div>"
        f"<div class=\"stat-label\">{html.escape(label)}</div>"
        f"<div class=\"stat-note\">{html.escape(note)}</div>"
        "</div>"
    )


def table_html(df: pd.DataFrame, columns: list[str], limit: int = 20, class_map: dict[str, str] | None = None) -> str:
    if df.empty:
        return "<p class=\"muted\">No rows available.</p>"
    cols = [c for c in columns if c in df.columns]
    if not cols:
        return "<p class=\"muted\">Requested columns unavailable.</p>"
    view = df.loc[:, cols].head(limit).copy()
    head = "".join(f"<th>{html.escape(c)}</th>" for c in cols)
    rows = []
    class_map = class_map or {}
    for _, row in view.iterrows():
        cells = []
        for col in cols:
            value = row[col]
            if isinstance(value, float):
                value = fmt_num(value, 4 if abs(value) < 1 else 2)
            cls = f" class=\"{class_map.get(col, '')}\"" if col in class_map else ""
            cells.append(f"<td{cls}>{safe(value)}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    return f"<div class=\"table-wrap\"><table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table></div>"


def metric_value(metrics: pd.DataFrame, score: str, field: str) -> float:
    if metrics.empty or "score" not in metrics.columns or field not in metrics.columns:
        return float("nan")
    sub = metrics.loc[metrics["score"].eq(score), field]
    if sub.empty:
        return float("nan")
    try:
        return float(sub.iloc[0])
    except Exception:
        return float("nan")


def topk_value(topk: pd.DataFrame, score: str, k: int, field: str) -> float:
    if topk.empty or "score" not in topk.columns or "top_k" not in topk.columns or field not in topk.columns:
        return float("nan")
    sub = topk.loc[topk["score"].eq(score) & topk["top_k"].eq(k), field]
    if sub.empty:
        return float("nan")
    try:
        return float(sub.iloc[0])
    except Exception:
        return float("nan")


def write_tradeoff_plot(metrics: pd.DataFrame, topk: pd.DataFrame, asset_dir: Path) -> str:
    asset_dir.mkdir(parents=True, exist_ok=True)
    out = asset_dir / "barneo_x_tradeoff.png"
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return ""

    labels = ["BAR-Neo", "BMA", "X-discovery", "X-claim-safe"]
    score_cols = ["barneo_score", "barneo_bma_score", "barneo_x_discovery_score", "barneo_x_claim_safe_score"]
    auprc = [metric_value(metrics, s, "apparent_auprc") for s in score_cols]
    leakage = [metric_value(metrics, s, "top10_high_leakage_fraction") for s in score_cols]
    precision = [
        topk_value(topk, "barneo_x_discovery_score", 10, "precision"),
        topk_value(topk, "barneo_x_claim_safe_score", 10, "precision"),
    ]
    top_labels = ["X-discovery", "X-claim-safe"]

    plt.style.use("dark_background")
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.8), dpi=180)
    fig.patch.set_facecolor("#0d1117")
    for ax in axes:
        ax.set_facecolor("#101820")
        ax.grid(axis="y", color="#2a3441", alpha=0.55, linewidth=0.8)
        ax.tick_params(colors="#c8d1dc", labelsize=8)
        for spine in ax.spines.values():
            spine.set_color("#2a3441")

    axes[0].bar(labels, auprc, color=["#e3b341", "#b58cff", "#5eead4", "#86efac"])
    axes[0].set_ylim(0, 1.05)
    axes[0].set_title("Apparent AUPRC", color="#e6edf3", fontsize=10)
    axes[0].tick_params(axis="x", rotation=22)

    axes[1].bar(labels, leakage, color=["#f87171", "#f87171", "#e3b341", "#86efac"])
    axes[1].set_ylim(0, 1.05)
    axes[1].set_title("Top10 high-leakage fraction", color="#e6edf3", fontsize=10)
    axes[1].tick_params(axis="x", rotation=22)

    axes[2].bar(top_labels, precision, color=["#5eead4", "#86efac"])
    axes[2].set_ylim(0, 1.05)
    axes[2].set_title("X-layer top10 precision", color="#e6edf3", fontsize=10)
    axes[2].tick_params(axis="x", rotation=22)

    fig.suptitle("BAR-Neo-X keeps precision while removing high-leakage top10 rows", color="#f8fafc", fontsize=12)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return f"assets/{ASSET_DIR_NAME}/{out.name}"


def component_frequency(explain: pd.DataFrame, column: str, limit: int = 12) -> pd.DataFrame:
    if explain.empty or column not in explain.columns:
        return pd.DataFrame()
    counts: dict[str, int] = {}
    for raw in explain[column].fillna("").astype(str):
        for piece in raw.split(";"):
            name = piece.strip().split("=")[0].strip()
            if name:
                counts[name] = counts.get(name, 0) + 1
    return (
        pd.DataFrame([{"factor": k, "n_top50_rows": v} for k, v in counts.items()])
        .sort_values(["n_top50_rows", "factor"], ascending=[False, True])
        .head(limit)
    )


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root)
    hub_root = Path(args.hub_root)
    asset_dir = hub_root / "assets" / ASSET_DIR_NAME
    hub_root.mkdir(parents=True, exist_ok=True)
    asset_dir.mkdir(parents=True, exist_ok=True)

    metrics = read_tsv(output_root / "barneo_x_metric_audit.tsv")
    topk = read_tsv(output_root / "barneo_x_topk_safety_audit.tsv")
    ablation = read_tsv(output_root / "barneo_x_ablation_audit.tsv")
    explain = read_tsv(output_root / "barneo_x_candidate_explanations.tsv")
    why_not = read_tsv(output_root / "barneo_x_why_not_audit.tsv")
    why_reasons = read_tsv(output_root / "barneo_x_why_not_reason_summary.tsv")
    why_blockers = read_tsv(output_root / "barneo_x_primary_blocker_summary.tsv")
    why_lanes = read_tsv(output_root / "barneo_x_rescue_lane_summary.tsv")
    summary = read_json(output_root / "BAR_NEO_X_INTERPRETABILITY_BOOST_SUMMARY.json")
    why_summary = read_json(output_root / "BAR_NEO_X_WHY_NOT_SUMMARY.json")
    manifest = read_json(output_root / "run_manifest.json")
    plot_rel = write_tradeoff_plot(metrics, topk, asset_dir)

    n_candidates = int(summary.get("n_candidates", len(explain)))
    priority = int(summary.get("n_priority_review_candidates", 0))
    top_claim = float(summary.get("top_claim_safe_score", metric_value(explain, "NA", "NA")))
    raw_auprc = metric_value(metrics, "barneo_score", "apparent_auprc")
    x_auprc = metric_value(metrics, "barneo_x_claim_safe_score", "apparent_auprc")
    raw_precision = metric_value(metrics, "barneo_score", "top10_precision")
    x_precision = metric_value(metrics, "barneo_x_claim_safe_score", "top10_precision")
    raw_leak = metric_value(metrics, "barneo_score", "top10_high_leakage_fraction")
    x_leak = metric_value(metrics, "barneo_x_claim_safe_score", "top10_high_leakage_fraction")
    x_top20_precision = metric_value(metrics, "barneo_x_claim_safe_score", "top20_precision")
    x_top20_leak = topk_value(topk, "barneo_x_claim_safe_score", 20, "high_leakage_fraction")
    hard_blocked = int(why_summary.get("n_hard_claim_blocked", 0))
    metadata_rescuable = int(why_summary.get("n_metadata_rescuable", 0))
    priority_now = int(why_summary.get("n_priority_now", 0))
    priority_if_metadata = int(why_summary.get("n_would_be_priority_if_metadata_complete", 0))

    top_claim_rows = explain.sort_values("barneo_x_claim_safe_rank").head(20) if "barneo_x_claim_safe_rank" in explain.columns else explain.head(20)
    top50 = explain.sort_values("barneo_x_claim_safe_rank").head(50) if "barneo_x_claim_safe_rank" in explain.columns else explain.head(50)
    why_top = why_not.sort_values("barneo_x_claim_safe_rank").head(30) if "barneo_x_claim_safe_rank" in why_not.columns else why_not.head(30)
    positive_frequency = component_frequency(top50, "barneo_x_top_positive_factors")
    negative_frequency = component_frequency(top50, "barneo_x_top_negative_factors")

    decision = pd.DataFrame(
        [
            {
                "surface": "raw BAR-Neo score",
                "best use": "internal discovery signal",
                "evidence": f"AUPRC {fmt_num(raw_auprc)}, top10 precision {fmt_num(raw_precision)}",
                "blocker": f"top10 high-leakage {fmt_num(raw_leak)}",
                "disposition": "do not use as clean reviewer-facing top10",
            },
            {
                "surface": "BAR-Neo-X discovery score",
                "best use": "hypothesis queue before manual audit",
                "evidence": f"top10 precision {fmt_num(metric_value(metrics, 'barneo_x_discovery_score', 'top10_precision'))}",
                "blocker": f"top20 high-leakage {fmt_num(topk_value(topk, 'barneo_x_discovery_score', 20, 'high_leakage_fraction'))}",
                "disposition": "internal only",
            },
            {
                "surface": "BAR-Neo-X claim-safe score",
                "best use": "reviewer-facing candidate triage",
                "evidence": f"top10 precision {fmt_num(x_precision)}, top20 precision {fmt_num(x_top20_precision)}",
                "blocker": "lower apparent AUPRC by design",
                "disposition": "primary dossier ranking",
            },
            {
                "surface": "public SOTA claim",
                "best use": "none in current evidence package",
                "evidence": "candidate-level apparent audit only",
                "blocker": "overlap and leakage controls dominate the clean claim",
                "disposition": "forbidden",
            },
            {
                "surface": "quantum advantage claim",
                "best use": "none",
                "evidence": "QUBO/QK pieces are bounded engineering heuristics",
                "blocker": "no quantum runtime or advantage test",
                "disposition": "forbidden",
            },
        ]
    )

    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>BAR-Neo-X Interpretability Boost</title>
  <style>
    :root {{ --bg:#0d1117; --panel:#151b23; --panel2:#101820; --ink:#e6edf3; --muted:#9aa7b4; --line:#2a3441; --gold:#e3b341; --cyan:#5eead4; --red:#f87171; --green:#86efac; --violet:#b58cff; }}
    * {{ box-sizing:border-box; }}
    body {{ margin:0; background:var(--bg); color:var(--ink); font-family:"JetBrains Mono", ui-monospace, SFMono-Regular, Menlo, monospace; line-height:1.55; }}
    a {{ color:var(--cyan); text-decoration:none; }}
    .hero {{ min-height:76vh; display:flex; align-items:flex-end; padding:58px 48px 40px; border-bottom:1px solid var(--line); background:linear-gradient(145deg, rgba(94,234,212,.10), transparent 40%), radial-gradient(circle at 78% 14%, rgba(134,239,172,.14), transparent 30%), #0d1117; }}
    .hero-inner {{ max-width:1220px; width:100%; }}
    .kicker {{ color:var(--gold); text-transform:uppercase; letter-spacing:.08em; font-size:13px; font-weight:700; }}
    h1 {{ font-family:Georgia, "Times New Roman", serif; font-size:clamp(42px, 7vw, 84px); line-height:.97; margin:14px 0 18px; letter-spacing:0; max-width:1080px; }}
    .lead {{ max-width:1000px; font-size:18px; color:#c8d1dc; }}
    .hero-links {{ display:flex; flex-wrap:wrap; gap:10px; margin-top:18px; }}
    .hero-links a {{ border:1px solid var(--line); padding:8px 11px; background:#101820; color:#d7fff6; }}
    .stats {{ display:grid; grid-template-columns:repeat(7, minmax(120px,1fr)); gap:12px; margin-top:28px; }}
    .stat {{ border:1px solid var(--line); background:rgba(21,27,35,.72); padding:14px; min-height:104px; }}
    .stat-value {{ font-size:24px; color:white; font-weight:800; }}
    .stat-label {{ color:var(--cyan); font-size:12px; margin-top:4px; }}
    .stat-note {{ color:var(--muted); font-size:11px; margin-top:8px; }}
    .layout {{ display:grid; grid-template-columns:260px minmax(0,1fr); gap:30px; max-width:1420px; margin:0 auto; padding:34px 28px 80px; }}
    nav {{ position:sticky; top:0; align-self:start; max-height:100vh; overflow:auto; padding:18px 0; }}
    nav a {{ display:block; padding:8px 10px; color:#b9c4d0; border-left:2px solid transparent; font-size:13px; }}
    nav a:hover {{ border-left-color:var(--cyan); color:white; }}
    section {{ padding:24px 0 36px; border-bottom:1px solid var(--line); }}
    h2 {{ font-family:Georgia, "Times New Roman", serif; font-size:34px; margin:0 0 14px; letter-spacing:0; }}
    h3 {{ font-size:18px; margin:24px 0 10px; color:var(--gold); }}
    .num {{ color:var(--gold); margin-right:10px; font-size:18px; }}
    .grid {{ display:grid; grid-template-columns:repeat(2, minmax(0,1fr)); gap:14px; }}
    .cards {{ display:grid; grid-template-columns:repeat(3, minmax(0,1fr)); gap:14px; }}
    .box,.card {{ border:1px solid var(--line); background:var(--panel); padding:16px; }}
    .box strong,.card strong {{ color:white; }}
    .muted {{ color:var(--muted); }}
    .warn {{ color:var(--red); font-weight:700; }}
    .ok {{ color:var(--green); font-weight:700; }}
    .gold {{ color:var(--gold); font-weight:700; }}
    .cyan {{ color:var(--cyan); font-weight:700; }}
    .table-wrap {{ overflow:auto; border:1px solid var(--line); background:var(--panel2); margin:14px 0; }}
    table {{ width:100%; border-collapse:collapse; font-size:12px; }}
    th, td {{ padding:9px 10px; border-bottom:1px solid var(--line); text-align:left; vertical-align:top; white-space:nowrap; }}
    th {{ color:var(--gold); background:#10151c; position:sticky; top:0; }}
    td {{ color:#d6dee7; }}
    td.good {{ color:var(--green); font-weight:700; }}
    td.warn {{ color:var(--red); font-weight:700; }}
    .figure {{ border:1px solid var(--line); background:#101820; padding:12px; margin:16px 0; }}
    .figure img {{ display:block; max-width:100%; height:auto; margin:0 auto; }}
    .path {{ color:#b9c4d0; font-size:12px; overflow-wrap:anywhere; }}
    .kakao {{ white-space:pre-wrap; background:#0b1320; border:1px solid var(--line); padding:16px; color:#edf7ff; }}
    @media (max-width: 920px) {{ .layout {{ grid-template-columns:1fr; padding:24px 16px 64px; }} nav {{ position:relative; max-height:none; }} .stats,.cards,.grid {{ grid-template-columns:1fr; }} .hero {{ padding:42px 20px 28px; }} th,td {{ white-space:normal; }} }}
  </style>
</head>
<body>
  <header class="hero">
    <div class="hero-inner">
      <div class="kicker">Cancer vaccine benchmark · claim-safe XAI layer · 2026-05-09</div>
      <h1>BAR-Neo-X Interpretability Boost</h1>
      <p class="lead">A post-hoc explanation and claim-safe reranking layer for CLEAN-NeoBench/BAR-Neo. It does not claim public SOTA, clinical vaccine selection, or quantum advantage. Its value is narrower and stronger: keep candidate-level precision while making the top reviewer-facing list low-leakage and explainable.</p>
      <div class="hero-links">
        <a href="index.html">Hub index</a>
        <a href="clean_neobench_barneo_dossier_2026_05_09.html">CLEAN-NeoBench/BAR-Neo dossier</a>
        <a href="cancer_vaccine_full_dossier.html">Cancer vaccine full dossier</a>
        <a href="assets/barneo_x/barneo_x_reviewer_packet_2026_05_09.zip">Reviewer packet ZIP</a>
      </div>
      <div class="stats">
        {stat(f"{n_candidates:,}", "Candidates", "CLEAN-NeoBench rows")}
        {stat(fmt_num(raw_auprc), "Raw BAR-Neo AUPRC", "apparent candidate-level audit")}
        {stat(fmt_num(x_auprc), "X claim-safe AUPRC", "lower by design")}
        {stat(f"{fmt_num(raw_leak)} -> {fmt_num(x_leak)}", "Top10 high-leakage", "BAR-Neo to BAR-Neo-X claim-safe")}
        {stat(fmt_num(x_precision), "Claim-safe top10 precision", f"top20 {fmt_num(x_top20_precision)}, top20 leak {fmt_num(x_top20_leak)}")}
        {stat(str(priority), "Priority rows", f"top claim-safe score {fmt_num(top_claim)}")}
        {stat(f"{hard_blocked:,}", "Hard blocked", f"{metadata_rescuable:,} metadata-rescuable")}
      </div>
    </div>
  </header>

  <div class="layout">
    <nav>
      <a href="#tldr">01 TL;DR</a>
      <a href="#task">02 Data, task, training</a>
      <a href="#tradeoff">03 Performance tradeoff</a>
      <a href="#ablation">04 Ablation audit</a>
      <a href="#algorithm">05 Algorithm</a>
      <a href="#candidates">06 Top candidates</a>
      <a href="#factors">07 Factor audit</a>
      <a href="#why-not">08 Why-not audit</a>
      <a href="#decision">09 Decision matrix</a>
      <a href="#kakao">10 Kakao payload</a>
      <a href="#paths">11 Sources + paths</a>
    </nav>

    <main>
      <section id="tldr">
        <h2><span class="num">01</span>TL;DR</h2>
        <div class="grid">
          <div class="box"><strong>What improved:</strong> top10 high-leakage fraction drops from <span class="warn">{fmt_num(raw_leak)}</span> to <span class="ok">{fmt_num(x_leak)}</span> while claim-safe top10 precision stays <span class="ok">{fmt_num(x_precision)}</span>.</div>
          <div class="box"><strong>What did not improve:</strong> apparent AUPRC is lower than raw BAR-Neo, from <span class="gold">{fmt_num(raw_auprc)}</span> to <span class="gold">{fmt_num(x_auprc)}</span>, because leakage-heavy rows are deliberately downranked.</div>
          <div class="box"><strong>What became interpretable:</strong> each candidate receives positive components and negative penalties: BAR-Neo evidence, BMA consensus, internal predictor support, confidence, leakage, overlap, patient gate incompleteness, uncertainty, and disagreement.</div>
          <div class="box"><strong>Why-not layer:</strong> {hard_blocked:,} rows are hard claim-blocked, {metadata_rescuable:,} are metadata-rescuable, {priority_if_metadata:,} would cross priority if metadata were complete, and {priority_now:,} is priority now.</div>
        </div>
      </section>

      <section id="task">
        <h2><span class="num">02</span>Data, Task, Training</h2>
        <div class="cards">
          <div class="card"><strong>Data used</strong><br>CLEAN-NeoBench master rows, method scores, BAR-Neo candidate scores, BAR-Neo-BMA scores, and row-level leakage/overlap flags.</div>
          <div class="card"><strong>Task</strong><br>Candidate-level neoantigen/immunogenicity triage with a separate reviewer-facing clean-claim ranking.</div>
          <div class="card"><strong>Training</strong><br>No new base predictor is trained in BAR-Neo-X. It is a deterministic post-hoc reranking and explanation layer on existing BAR-Neo/BMA outputs.</div>
          <div class="card"><strong>Model family</strong><br>Additive component attribution plus leakage-aware claim gate and metadata cap.</div>
          <div class="card"><strong>Quantum status</strong><br>None. Existing QUBO/QK terms remain bounded heuristics in the broader pipeline, not a quantum advantage claim.</div>
          <div class="card"><strong>Clinical status</strong><br>Research triage only. Candidate review requires patient metadata, disease context, presentation status, expression, clonality, and wet-lab validation.</div>
        </div>
      </section>

      <section id="tradeoff">
        <h2><span class="num">03</span>Performance Tradeoff</h2>
        <p class="muted">The right comparison is not raw metric maximization. The reviewer-facing objective is precision under leakage control.</p>
        {'<div class="figure"><img src="' + html.escape(plot_rel) + '" alt="BAR-Neo-X metric tradeoff plot"></div>' if plot_rel else ''}
        {table_html(metrics, ["score", "n", "apparent_auprc", "apparent_auroc", "top10_precision", "top20_precision", "top10_high_leakage_fraction"], 10)}
        <h3>Top-k safety audit</h3>
        {table_html(topk, ["score", "top_k", "precision", "high_leakage_fraction", "median_claim_safe_score", "priority_review_rows"], 12)}
      </section>

      <section id="ablation">
        <h2><span class="num">04</span>Ablation Audit</h2>
        <p class="muted">This table shows what breaks when penalties, leakage gating, or metadata capping are removed. The goal is to justify the conservative score, not maximize apparent AUPRC.</p>
        {table_html(ablation, ["score", "n", "apparent_auprc", "apparent_auroc", "top10_precision", "top20_precision", "top10_high_leakage_fraction", "top20_high_leakage_fraction", "interpretation"], 10)}
      </section>

      <section id="algorithm">
        <h2><span class="num">05</span>Algorithm</h2>
        <div class="grid">
          <div class="box"><strong>Positive evidence:</strong> BAR-Neo model evidence, BMA expert consensus, clean internal predictor support, anchor support, confidence support.</div>
          <div class="box"><strong>Negative evidence:</strong> posterior disagreement, posterior uncertainty, patient gate incompleteness, high/medium leakage, exact or near peptide overlap, study/patient overlap, public pretrained overlap, low-prevalence source, sparse expert support.</div>
          <div class="box"><strong>Discovery score:</strong> positive component sum minus a penalty-weighted risk sum. Use for internal hypothesis generation.</div>
          <div class="box"><strong>Claim-safe score:</strong> discovery score multiplied by leakage gate and metadata cap. Use for reviewer-facing candidate review.</div>
        </div>
      </section>

      <section id="candidates">
        <h2><span class="num">06</span>Top Claim-Safe Candidates</h2>
        <p class="muted">Rows are sorted by BAR-Neo-X claim-safe rank. Labels are benchmark labels, not clinical validation.</p>
        {table_html(top_claim_rows, ["candidate_id", "peptide", "hla_allele_4digit", "label", "leakage_risk_level", "barneo_score", "barneo_bma_score", "barneo_x_claim_safe_score", "barneo_x_primary_action", "barneo_x_top_positive_factors", "barneo_x_top_negative_factors"], 20, {"barneo_x_primary_action": "good"})}
      </section>

      <section id="factors">
        <h2><span class="num">07</span>Factor Audit</h2>
        <div class="grid">
          <div>
            <h3>Top positive factors in claim-safe top50</h3>
            {table_html(positive_frequency, ["factor", "n_top50_rows"], 12)}
          </div>
          <div>
            <h3>Top negative factors in claim-safe top50</h3>
            {table_html(negative_frequency, ["factor", "n_top50_rows"], 12)}
          </div>
        </div>
      </section>

      <section id="why-not">
        <h2><span class="num">08</span>Why-Not / Rescue Audit</h2>
        <p class="muted">This layer answers the failure case question directly: whether a row is blocked by leakage/identity overlap, missing patient context, model uncertainty, sparse support, or low-priority benchmark status.</p>
        <div class="grid">
          <div class="box"><strong>Hard claim-blocked:</strong> <span class="warn">{hard_blocked:,}</span> candidates need a new holdout, row-level overlap audit, or independent external cohort before clean benchmark claims.</div>
          <div class="box"><strong>Metadata-rescuable:</strong> <span class="gold">{metadata_rescuable:,}</span> candidates can be reprioritized after cancer type, disease context, stage, expression, VAF/clonality, HLA-LOH/B2M, and immune context are completed.</div>
          <div class="box"><strong>Priority now:</strong> <span class="ok">{priority_now:,}</span> candidate remains reviewer-triage priority under the claim-safe score.</div>
          <div class="box"><strong>Upside if metadata completes:</strong> <span class="cyan">{priority_if_metadata:,}</span> candidates would cross the priority threshold without removing the leakage gate.</div>
        </div>
        <h3>Primary blocker summary</h3>
        {table_html(why_blockers, ["primary_blocker", "n_candidates", "n_positive_label", "median_claim_safe_score", "n_would_be_priority_if_metadata_complete", "top_rescue_lane", "example_required_next_data"], 12)}
        <h3>Rescue lane summary</h3>
        {table_html(why_lanes, ["rescue_lane", "n_candidates", "n_positive_label", "median_claim_safe_score", "example_required_next_data", "claim_boundary"], 12)}
        <h3>Top why-not reasons</h3>
        {table_html(why_reasons, ["reason", "n_candidates", "n_positive_label", "median_claim_safe_score", "top_rescue_lane"], 16)}
        <h3>Top rows with blockers</h3>
        {table_html(why_top, ["candidate_id", "peptide", "hla_allele_4digit", "label", "leakage_risk_level", "barneo_x_claim_safe_score", "primary_blocker", "rescue_lane", "would_be_priority_if_metadata_complete", "required_next_data"], 30)}
      </section>

      <section id="decision">
        <h2><span class="num">09</span>Decision Matrix</h2>
        {table_html(decision, ["surface", "best use", "evidence", "blocker", "disposition"], 10)}
      </section>

      <section id="kakao">
        <h2><span class="num">10</span>Kakao Payload</h2>
        <div class="kakao">Cancer vaccine neoantigen benchmark에서 BAR-Neo 기반 BAR-Neo-X를 만들었고, 후보별 BAR-Neo evidence/BMA consensus/internal predictor/confidence와 leakage-overlap-uncertainty penalty를 분해해 설명 가능하게 만들었습니다. raw BAR-Neo는 apparent AUPRC {fmt_num(raw_auprc)}지만 top10이 전부 high-leakage라 SOTA 주장엔 위험하고, BAR-Neo-X claim-safe score는 AUPRC {fmt_num(x_auprc)}으로 보수화되는 대신 top10 high-leakage를 {fmt_num(raw_leak)}->{fmt_num(x_leak)}로 제거하면서 top10 precision {fmt_num(x_precision)}를 유지합니다. 안 될 때도 이유가 분해됩니다: {hard_blocked:,}개는 leakage/identity overlap 때문에 clean claim 불가, {metadata_rescuable:,}개는 patient/disease metadata 보강 시 구제 가능, {priority_if_metadata:,}개는 metadata 완성 시 priority threshold 통과 가능, 현재 priority-now는 {priority_now:,}개입니다. 즉 무조건 SOTA/clinical/quantum advantage가 아니라, 리뷰어 방어 가능한 해석형-누수차단 neoantigen triage layer입니다.</div>
      </section>

      <section id="paths">
        <h2><span class="num">11</span>Sources + Paths</h2>
        <p class="path">Output root: {safe(output_root)}</p>
        <p class="path">Metric audit: {safe(output_root / "barneo_x_metric_audit.tsv")}</p>
        <p class="path">Top-k safety audit: {safe(output_root / "barneo_x_topk_safety_audit.tsv")}</p>
        <p class="path">Candidate explanations: {safe(output_root / "barneo_x_candidate_explanations.tsv")}</p>
        <p class="path">Candidate scores: {safe(output_root / "barneo_x_candidate_scores.tsv")}</p>
        <p class="path">Why-not audit: {safe(output_root / "barneo_x_why_not_audit.tsv")}</p>
        <p class="path">Primary blocker summary: {safe(output_root / "barneo_x_primary_blocker_summary.tsv")}</p>
        <p class="path">Summary JSON: {safe(output_root / "BAR_NEO_X_INTERPRETABILITY_BOOST_SUMMARY.json")}</p>
        <p class="path">Reviewer packet ZIP: assets/barneo_x/barneo_x_reviewer_packet_2026_05_09.zip</p>
        <p class="path">Manifest stage: {safe('barneo_x_interpretability_boost' if 'barneo_x_interpretability_boost' in manifest.get('stages', {}) else 'recorded in run_manifest.json')}</p>
      </section>
    </main>
  </div>
</body>
</html>
"""
    out_path = hub_root / args.page_name
    out_path.write_text(html_text)
    print(json.dumps({"html": str(out_path), "plot": plot_rel, "n_candidates": n_candidates}, indent=2))


if __name__ == "__main__":
    main()
