#!/usr/bin/env python3
"""Write a CLEAN-NeoBench challenge-pack dossier and portable packet."""

from __future__ import annotations

import argparse
import csv
import html
import json
import shutil
import zipfile
from pathlib import Path
from typing import Any

import pandas as pd


PAGE_NAME = "clean_neobench_challenge_pack_2026_05_10.html"
ASSET_DIR_NAME = "clean_neobench_challenge"
PACKET_NAME = "clean_neobench_challenge_pack_2026_05_10"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", required=True, help="CLEAN-NeoBench BAR-Neo output directory")
    parser.add_argument("--hub-root", default="project/papers_hub_2026_05_04", help="HTML hub directory")
    parser.add_argument("--page-name", default=PAGE_NAME)
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def safe(value: object) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    return html.escape(str(value))


def fmt(value: object, digits: int = 3) -> str:
    try:
        if pd.isna(value):
            return "NA"
        x = float(value)
        if abs(x) >= 1000:
            return f"{x:,.0f}"
        return f"{x:.{digits}f}"
    except Exception:
        return safe(value)


def stat(value: str, label: str, note: str = "") -> str:
    return (
        "<div class=\"stat\">"
        f"<div class=\"stat-value\">{safe(value)}</div>"
        f"<div class=\"stat-label\">{safe(label)}</div>"
        f"<div class=\"stat-note\">{safe(note)}</div>"
        "</div>"
    )


def table_html(df: pd.DataFrame, columns: list[str], limit: int = 20) -> str:
    if df.empty:
        return "<p class=\"muted\">No rows available.</p>"
    cols = [c for c in columns if c in df.columns]
    if not cols:
        return "<p class=\"muted\">Requested columns unavailable.</p>"
    view = df.loc[:, cols].head(limit).copy()
    head = "".join(f"<th>{safe(c)}</th>" for c in cols)
    rows = []
    for _, row in view.iterrows():
        cells = []
        for col in cols:
            value = row[col]
            if isinstance(value, float):
                value = fmt(value, 4 if abs(value) < 1 else 2)
            cells.append(f"<td>{safe(value)}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    return f"<div class=\"table-wrap\"><table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table></div>"


def axis_count(summary: pd.DataFrame, axis: str, field: str = "n_unique_candidates") -> int:
    if summary.empty or "challenge_axis" not in summary.columns or field not in summary.columns:
        return 0
    sub = summary.loc[summary["challenge_axis"].eq(axis), field]
    if sub.empty:
        return 0
    try:
        return int(float(sub.iloc[0]))
    except Exception:
        return 0


def write_axis_plot(summary: pd.DataFrame, asset_dir: Path) -> str:
    if summary.empty or "challenge_axis" not in summary.columns:
        return ""
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return ""
    asset_dir.mkdir(parents=True, exist_ok=True)
    out = asset_dir / "clean_neobench_challenge_axis_summary.png"
    view = summary.sort_values("n_unique_candidates", ascending=True).copy()
    labels = view["challenge_axis"].astype(str).str.replace("_", "\n")
    counts = pd.to_numeric(view["n_unique_candidates"], errors="coerce").fillna(0)
    prevalence = pd.to_numeric(view.get("positive_prevalence"), errors="coerce").fillna(0)
    colors = ["#86efac" if p >= 0.75 else "#e3b341" if p > 0 else "#f87171" for p in prevalence]

    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(9.5, 5.7), dpi=180)
    fig.patch.set_facecolor("#0d1117")
    ax.set_facecolor("#101820")
    ax.barh(labels, counts, color=colors)
    ax.grid(axis="x", color="#2a3441", alpha=0.55)
    ax.tick_params(colors="#c8d1dc", labelsize=8)
    for spine in ax.spines.values():
        spine.set_color("#2a3441")
    ax.set_xlabel("unique candidates", color="#c8d1dc")
    ax.set_title("CLEAN-NeoBench challenge axes", color="#f8fafc", fontsize=12)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return f"assets/{ASSET_DIR_NAME}/{out.name}"


def parse_available_count(value: object) -> int:
    total = 0
    for piece in str(value).split(";"):
        if "=" not in piece:
            continue
        try:
            total += int(float(piece.rsplit("=", 1)[1].strip()))
        except Exception:
            continue
    return total


def write_experiment_plot(plan: pd.DataFrame, asset_dir: Path) -> str:
    if plan.empty or "experiment" not in plan.columns:
        return ""
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return ""
    asset_dir.mkdir(parents=True, exist_ok=True)
    out = asset_dir / "clean_neobench_next_experiment_board.png"
    view = plan.copy()
    priority_order = {"P0": 0, "P1": 1, "P2": 2}
    view["_priority_order"] = view["priority"].map(priority_order).fillna(9)
    available_raw = view["challenge_rows_available"] if "challenge_rows_available" in view.columns else pd.Series("", index=view.index)
    view["_available"] = available_raw.map(parse_available_count).fillna(0).astype(int)
    view = view.sort_values(["_priority_order", "_available"], ascending=[True, False])
    display_labels = {
        "public_training_corpus_row_overlap_audit": "public training\noverlap audit",
        "exact_near_overlap_lockdown": "exact/near\noverlap lockdown",
        "low_prevalence_topk_stress": "low-prevalence\ntop-k stress",
        "rare_hla_and_korean_hla_calibration": "rare/Korean HLA\ncalibration",
        "contextual_bma_clean_internal_ablation": "clean-internal\nBMA ablation",
        "PAAD_THCA_patient_gate_live_demo": "PAAD/THCA\npatient gate",
        "MHC_II_separate_benchmark": "MHC-II\nseparate board",
    }
    labels = view["experiment"].map(display_labels).fillna(view["experiment"].astype(str).str.replace("_", "\n"))
    values = view["_available"].clip(lower=3)
    colors = view["priority"].map({"P0": "#f87171", "P1": "#e3b341", "P2": "#5eead4"}).fillna("#9aa7b4")

    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(10.5, 5.6), dpi=180)
    fig.patch.set_facecolor("#0d1117")
    ax.set_facecolor("#101820")
    bars = ax.barh(labels, values, color=colors)
    ax.invert_yaxis()
    ax.grid(axis="x", color="#2a3441", alpha=0.55)
    ax.tick_params(colors="#c8d1dc", labelsize=8)
    for spine in ax.spines.values():
        spine.set_color("#2a3441")
    ax.set_xlabel("challenge rows available", color="#c8d1dc")
    ax.set_title("Next-experiment priority board", color="#f8fafc", fontsize=12)
    for bar, priority, available, experiment in zip(bars, view["priority"], view["_available"], view["experiment"], strict=False):
        value_label = "corpus audit" if available == 0 and experiment == "public_training_corpus_row_overlap_audit" else f"{available}"
        label = f"{priority} | {value_label}"
        ax.text(bar.get_width() + 1.8, bar.get_y() + bar.get_height() / 2, label, va="center", color="#e6edf3", fontsize=8)
    fig.subplots_adjust(left=0.24, right=0.94, top=0.90, bottom=0.14)
    fig.savefig(out, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return f"assets/{ASSET_DIR_NAME}/{out.name}"


def write_contract_matrix_plot(summary: pd.DataFrame, asset_dir: Path) -> str:
    if summary.empty or "recommended_split_contract" not in summary.columns:
        return ""
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception:
        return ""
    asset_dir.mkdir(parents=True, exist_ok=True)
    out = asset_dir / "clean_neobench_split_contract_matrix.png"
    tokens = [
        ("exact\npHLA", "exact_phla"),
        ("near\npeptide", "near_peptide"),
        ("source\nheldout", "source_heldout"),
        ("HLA\nheldout", "HLA_heldout"),
        ("low\nprev", "low_prevalence"),
        ("Korean\nHLA", "korean_hla"),
        ("public\noverlap", "public_overlap"),
        ("clean\ninternal", "clean_internal"),
        ("patient\ngate", "patient_gated"),
    ]
    view = summary.sort_values(["n_unique_candidates", "mean_priority_score"], ascending=[False, False]).copy()
    contracts = view["recommended_split_contract"].fillna("").astype(str).str.lower()
    matrix = []
    for contract in contracts:
        matrix.append([1 if token.lower() in contract else 0 for _, token in tokens])
    mat = np.array(matrix)
    labels = view["challenge_axis"].astype(str).str.replace("_", "\n")

    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(10.5, 6.2), dpi=180)
    fig.patch.set_facecolor("#0d1117")
    ax.set_facecolor("#101820")
    ax.imshow(mat, cmap="viridis", vmin=0, vmax=1, aspect="auto")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, color="#c8d1dc", fontsize=8)
    ax.set_xticks(range(len(tokens)))
    ax.set_xticklabels([label for label, _ in tokens], color="#c8d1dc", fontsize=8)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_color("#2a3441")
    for y in range(mat.shape[0]):
        for x in range(mat.shape[1]):
            ax.text(x, y, "1" if mat[y, x] else "", ha="center", va="center", color="#f8fafc", fontsize=8)
    ax.set_title("Split-contract matrix for claim eligibility", color="#f8fafc", fontsize=12)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return f"assets/{ASSET_DIR_NAME}/{out.name}"


def copy_if_exists(src: Path, dst: Path) -> bool:
    if not src.exists():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return True


def write_file_manifest(packet_dir: Path) -> None:
    rows = []
    for path in sorted(p for p in packet_dir.rglob("*") if p.is_file()):
        rows.append({"relative_path": str(path.relative_to(packet_dir)), "bytes": str(path.stat().st_size)})
    out = packet_dir / "FILE_MANIFEST.tsv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["relative_path", "bytes"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def zip_dir(packet_dir: Path, zip_path: Path) -> None:
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(p for p in packet_dir.rglob("*") if p.is_file()):
            zf.write(path, path.relative_to(packet_dir.parent))


def write_packet(
    output_root: Path,
    hub_root: Path,
    page_path: Path,
    plot_rels: list[str],
    summary_payload: dict[str, Any],
    kakao_text: str,
) -> Path:
    packet_dir = output_root / PACKET_NAME
    if packet_dir.exists():
        shutil.rmtree(packet_dir)
    packet_dir.mkdir(parents=True, exist_ok=True)
    (packet_dir / "CLEAN_NEOBENCH_CHALLENGE_CLAIM_CARD.json").write_text(
        json.dumps(summary_payload, indent=2, sort_keys=True) + "\n"
    )
    (packet_dir / "CLEAN_NEOBENCH_CHALLENGE_KAKAO_ONE_SHOT.txt").write_text(kakao_text + "\n")
    (packet_dir / "README.md").write_text(
        "\n".join(
            [
                "# CLEAN-NeoBench Challenge Pack",
                "",
                "Portable evidence bundle for converting benchmark failure modes into next experiments.",
                "",
                "Allowed claim: challenge benchmark design, manual-review prioritization, distribution stress testing, and reliability refinement.",
                "",
                "Forbidden claims: clinical vaccine selection, public SOTA, quantum advantage, and clean public baseline claims without row-level overlap audit.",
                "",
            ]
        )
    )
    for src in [
        output_root / "CLEAN_NEOBENCH_CHALLENGE_PACK.md",
        output_root / "CLEAN_NEOBENCH_CHALLENGE_PACK_KR.md",
        output_root / "clean_neobench_challenge_pack.tsv",
        output_root / "clean_neobench_challenge_axis_summary.tsv",
        output_root / "clean_neobench_next_experiment_plan.tsv",
        page_path,
    ]:
        copy_if_exists(src, packet_dir / src.name)
    for plot_rel in plot_rels:
        copy_if_exists(hub_root / plot_rel, packet_dir / Path(plot_rel).name)
    write_file_manifest(packet_dir)
    zip_path = output_root / f"{PACKET_NAME}.zip"
    zip_dir(packet_dir, zip_path)
    copy_if_exists(zip_path, hub_root / "assets" / ASSET_DIR_NAME / zip_path.name)
    return zip_path


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root)
    hub_root = Path(args.hub_root)
    asset_dir = hub_root / "assets" / ASSET_DIR_NAME
    hub_root.mkdir(parents=True, exist_ok=True)
    asset_dir.mkdir(parents=True, exist_ok=True)

    pack = read_tsv(output_root / "clean_neobench_challenge_pack.tsv")
    summary = read_tsv(output_root / "clean_neobench_challenge_axis_summary.tsv")
    plan = read_tsv(output_root / "clean_neobench_next_experiment_plan.tsv")
    if pack.empty or summary.empty or plan.empty:
        raise SystemExit("missing challenge-pack inputs")

    axis_plot_rel = write_axis_plot(summary, asset_dir)
    experiment_plot_rel = write_experiment_plot(plan, asset_dir)
    contract_plot_rel = write_contract_matrix_plot(summary, asset_dir)
    plot_rels = [p for p in [axis_plot_rel, experiment_plot_rel, contract_plot_rel] if p]
    n_rows = len(pack)
    n_axes = int(summary["challenge_axis"].nunique()) if "challenge_axis" in summary.columns else 0
    n_p0 = int(plan["priority"].astype(str).eq("P0").sum()) if "priority" in plan.columns else 0
    priority_now = axis_count(summary, "claim_safe_priority_review")
    high_blocked = axis_count(summary, "high_score_claim_blocked")
    low_prev = axis_count(summary, "low_prevalence_false_positive_stress")
    rare_korean = axis_count(summary, "rare_hla_support_gap") + axis_count(summary, "korean_hla_focus_stress")

    top_global = pack.sort_values("global_challenge_rank").head(80) if "global_challenge_rank" in pack.columns else pack.head(80)
    axis_table = summary.sort_values(["n_unique_candidates", "mean_priority_score"], ascending=[False, False])

    summary_payload = {
        "date": "2026-05-10",
        "n_challenge_rows": int(n_rows),
        "n_challenge_axes": int(n_axes),
        "n_p0_experiments": int(n_p0),
        "claim_safe_priority_review_candidates": int(priority_now),
        "high_score_claim_blocked_candidates": int(high_blocked),
        "low_prevalence_false_positive_stress_candidates": int(low_prev),
        "rare_or_korean_hla_candidates": int(rare_korean),
        "visualizations": [
            "clean_neobench_challenge_axis_summary.png",
            "clean_neobench_next_experiment_board.png",
            "clean_neobench_split_contract_matrix.png",
        ],
        "allowed_claim": "Challenge benchmark design and next-experiment prioritization.",
        "forbidden_claims": [
            "clinical vaccine selection",
            "public SOTA",
            "quantum advantage",
            "clean public baseline without row-level overlap audit",
        ],
    }
    kakao_text = (
        f"시각화까지 붙였습니다. CLEAN-NeoBench Challenge Pack은 총 {n_rows:,}개 challenge rows / {n_axes}개 axes를 "
        f"axis summary, next-experiment board, split-contract matrix 3장으로 보여줍니다. P0는 public-training row overlap audit와 "
        f"exact/near overlap lockdown이고, 핵심 큐는 high-score claim-blocked {high_blocked}개, "
        f"low-prevalence false-positive stress {low_prev}개, rare/Korean-HLA calibration {rare_korean}개, "
        f"claim-safe priority {priority_now}개입니다. 결론은 SOTA/clinical claim이 아니라, 어떤 후보가 어떤 split contract와 "
        "success metric을 통과해야 claim이 되는지 고정한 reviewer 방어용 next-experiment board입니다."
    )

    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>CLEAN-NeoBench Challenge Pack</title>
  <style>
    :root {{ --bg:#0d1117; --panel:#151b23; --panel2:#101820; --ink:#e6edf3; --muted:#9aa7b4; --line:#2a3441; --gold:#e3b341; --cyan:#5eead4; --red:#f87171; --green:#86efac; --violet:#b58cff; }}
    * {{ box-sizing:border-box; }}
    body {{ margin:0; background:var(--bg); color:var(--ink); font-family:"JetBrains Mono", ui-monospace, SFMono-Regular, Menlo, monospace; line-height:1.55; }}
    a {{ color:var(--cyan); text-decoration:none; }}
    .hero {{ min-height:70vh; display:flex; align-items:flex-end; padding:54px 48px 38px; border-bottom:1px solid var(--line); background:#0d1117; }}
    .hero-inner {{ max-width:1240px; width:100%; }}
    .kicker {{ color:var(--gold); text-transform:uppercase; letter-spacing:.08em; font-size:13px; font-weight:700; }}
    h1 {{ font-family:Georgia, "Times New Roman", serif; font-size:clamp(42px, 7vw, 82px); line-height:.98; margin:14px 0 18px; letter-spacing:0; max-width:1100px; }}
    .lead {{ max-width:1040px; font-size:18px; color:#c8d1dc; }}
    .hero-links {{ display:flex; flex-wrap:wrap; gap:10px; margin-top:18px; }}
    .hero-links a {{ border:1px solid var(--line); padding:8px 11px; background:#101820; color:#d7fff6; }}
    .stats {{ display:grid; grid-template-columns:repeat(6, minmax(130px,1fr)); gap:12px; margin-top:28px; }}
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
    .figure {{ border:1px solid var(--line); background:#101820; padding:12px; margin:16px 0; }}
    .figure img {{ display:block; max-width:100%; height:auto; margin:0 auto; }}
    .path {{ color:#b9c4d0; font-size:12px; overflow-wrap:anywhere; }}
    .kakao {{ white-space:pre-wrap; background:#0b1320; border:1px solid var(--line); padding:16px; color:#edf7ff; }}
    @media (max-width: 980px) {{ .layout {{ grid-template-columns:1fr; padding:24px 16px 64px; }} nav {{ position:relative; max-height:none; }} .stats,.cards,.grid {{ grid-template-columns:1fr; }} .hero {{ padding:42px 20px 28px; }} th,td {{ white-space:normal; }} }}
  </style>
</head>
<body>
  <header class="hero">
    <div class="hero-inner">
      <div class="kicker">CLEAN-NeoBench · failure-to-experiment conversion · 2026-05-10</div>
      <h1>CLEAN-NeoBench Challenge Pack</h1>
      <p class="lead">A reviewer-facing queue that turns benchmark failure modes into fixed challenge axes, split contracts, and success metrics. It complements BAR-Neo-X: after explaining why a row is blocked or rescuable, this pack says what experiment should happen next.</p>
      <div class="hero-links">
        <a href="index.html">Hub index</a>
        <a href="cancer_vaccine_barneo_x_interpretability_2026_05_09.html">BAR-Neo-X dossier</a>
        <a href="clean_neobench_barneo_dossier_2026_05_09.html">CLEAN-NeoBench dossier</a>
        <a href="assets/{ASSET_DIR_NAME}/{PACKET_NAME}.zip">Reviewer packet ZIP</a>
      </div>
      <div class="stats">
        {stat(f"{n_rows:,}", "Challenge rows", "fixed stress/manual-review queue")}
        {stat(str(n_axes), "Challenge axes", "distribution and claim blockers")}
        {stat(str(n_p0), "P0 experiments", "overlap audit + lockdown")}
        {stat(str(high_blocked), "High-score blocked", "exact/near overlap contract")}
        {stat(str(low_prev), "Low-prevalence stress", "negative-heavy sources")}
        {stat(str(rare_korean), "Rare/Korean HLA", "allele calibration queue")}
      </div>
    </div>
  </header>

  <div class="layout">
    <nav>
      <a href="#tldr">01 TL;DR</a>
      <a href="#axes">02 Challenge axes</a>
      <a href="#experiments">03 Next experiments</a>
      <a href="#rows">04 Top rows</a>
      <a href="#contracts">05 Claim contracts</a>
      <a href="#kakao">06 Kakao payload</a>
      <a href="#paths">07 Sources + paths</a>
    </nav>
    <main>
      <section id="tldr">
        <h2><span class="num">01</span>TL;DR</h2>
        <div class="grid">
          <div class="box"><strong>What this adds:</strong> {n_rows:,} challenge rows across {n_axes} axes, each with a recommended action, split contract, success metric, and claim boundary.</div>
          <div class="box"><strong>Immediate P0:</strong> public-training row-overlap audit and exact/near overlap lockdown before clean public comparator or SOTA language.</div>
          <div class="box"><strong>Practical high value:</strong> {low_prev} low-prevalence false-positive stress rows, {rare_korean} rare/Korean-HLA calibration rows, and {priority_now} claim-safe priority row.</div>
          <div class="box"><strong>Allowed claim:</strong> challenge benchmark design and next-experiment prioritization. <span class="warn">Forbidden:</span> clinical selection, public SOTA, or quantum advantage.</div>
        </div>
      </section>

      <section id="axes">
        <h2><span class="num">02</span>Challenge Axes</h2>
        <p class="muted">Axes are fixed failure contexts, not a validation leaderboard. They are designed to make future comparisons harder and more honest.</p>
        {'<div class="figure"><img src="' + safe(axis_plot_rel) + '" alt="CLEAN-NeoBench challenge axis summary"></div>' if axis_plot_rel else ''}
        {table_html(axis_table, ["challenge_axis", "n_unique_candidates", "n_pos", "positive_prevalence", "mean_priority_score", "mean_contextual_bma_score", "mean_confidence", "recommended_next_action", "recommended_split_contract", "success_metric"], 20)}
      </section>

      <section id="experiments">
        <h2><span class="num">03</span>Next Experiment Plan</h2>
        {'<div class="figure"><img src="' + safe(experiment_plot_rel) + '" alt="CLEAN-NeoBench next experiment priority board"></div>' if experiment_plot_rel else ''}
        {table_html(plan, ["priority", "experiment", "input_rows", "why", "success_metric", "expected_output", "challenge_rows_available"], 20)}
      </section>

      <section id="rows">
        <h2><span class="num">04</span>Top Challenge Rows</h2>
        <p class="muted">These are review prompts. A high row means the candidate is informative for a stress test, not clinically selected.</p>
        {table_html(top_global, ["global_challenge_rank", "challenge_axis", "candidate_id", "label", "source_name", "source_positive_prevalence", "hla_allele_4digit", "hla_allele_support_count", "contextual_bma_score", "clean_contextual_bma_score", "contextual_confidence_score", "recommended_next_action", "recommended_split_contract"], 80)}
      </section>

      <section id="contracts">
        <h2><span class="num">05</span>Claim Contracts</h2>
        {'<div class="figure"><img src="' + safe(contract_plot_rel) + '" alt="CLEAN-NeoBench split contract matrix"></div>' if contract_plot_rel else ''}
        <div class="cards">
          <div class="card"><strong>Overlap lockdown</strong><br>Exact peptide-HLA, near peptide cluster, source-heldout, and HLA-heldout contracts decide whether high-score rows can become clean evidence.</div>
          <div class="card"><strong>Low-prevalence stress</strong><br>TESLA-like negative-heavy sources require top-k false-positive pressure, not only aggregate AUPRC.</div>
          <div class="card"><strong>Rare/Korean HLA</strong><br>Allele-specific support and calibration are required before using rare or Korean-focus rows as evidence.</div>
          <div class="card"><strong>Public tools</strong><br>Public pretrained agreement stays caveated until row-level training-corpus overlap is audited.</div>
          <div class="card"><strong>Patient gate</strong><br>PAAD/THCA patient-gated triage is demo-only until disease, presentation, expression, clonality, HLA-LOH/B2M, and immune-context metadata are complete.</div>
          <div class="card"><strong>Class boundary</strong><br>MHC-I and MHC-II need separate benchmark contracts; pooled predictor claims remain disallowed.</div>
        </div>
      </section>

      <section id="kakao">
        <h2><span class="num">06</span>Kakao Payload</h2>
        <div class="kakao">더 붙였습니다. BAR-Neo-X가 “왜 안 되는지”를 설명했다면, CLEAN-NeoBench Challenge Pack은 그 실패 이유를 바로 다음 실험 큐로 바꿉니다. 총 {n_rows:,}개 challenge rows / {n_axes}개 axes를 만들었고, P0는 public-training row overlap audit와 exact/near overlap lockdown입니다. 핵심 큐는 high-score claim-blocked {high_blocked}개, low-prevalence false-positive stress {low_prev}개, rare/Korean-HLA calibration {rare_korean}개, claim-safe priority {priority_now}개입니다. 즉 “좋아 보이는 후보”를 그냥 미는 게 아니라, 어떤 split contract와 success metric을 통과해야 claim이 되는지 고정한 reviewer 방어용 next-experiment board입니다.</div>
      </section>

      <section id="paths">
        <h2><span class="num">07</span>Sources + Paths</h2>
        <p class="path">Output root: {safe(output_root)}</p>
        <p class="path">Challenge pack TSV: {safe(output_root / "clean_neobench_challenge_pack.tsv")}</p>
        <p class="path">Axis summary: {safe(output_root / "clean_neobench_challenge_axis_summary.tsv")}</p>
        <p class="path">Next experiment plan: {safe(output_root / "clean_neobench_next_experiment_plan.tsv")}</p>
        <p class="path">KR report: {safe(output_root / "CLEAN_NEOBENCH_CHALLENGE_PACK_KR.md")}</p>
        <p class="path">Reviewer packet ZIP: assets/{ASSET_DIR_NAME}/{PACKET_NAME}.zip</p>
      </section>
    </main>
  </div>
</body>
</html>
"""
    page_path = hub_root / args.page_name
    page_path.write_text(html_text)
    zip_path = write_packet(output_root, hub_root, page_path, plot_rels, summary_payload, kakao_text)
    print(
        json.dumps(
            {
                "html": str(page_path),
                "plots": plot_rels,
                "zip": str(zip_path),
                "n_rows": n_rows,
                "n_axes": n_axes,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
