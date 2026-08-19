#!/usr/bin/env python3
"""Compare v6 assay hit/fail calls against BigMHC, CROSS-Neo and integrated scores.

The default input is the existing confirmatory-threshold smoke sheet. When real
wetlab results are available, point --candidate-calls at the interpreter output
directory's candidate_calls_v6.tsv and set --result-source accordingly.
"""

from __future__ import annotations

import argparse
import html
import json
import shutil
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
EXEC = ROOT / "project/results/cross_neo_kaggle_winner_playbook_2026_05_10/execution_packet_v6_2026_05_10"
DEFAULT_CALLS = EXEC / "interpreter_smoke_confirmatory_threshold_calls_v6/candidate_calls_v6.tsv"
DEFAULT_ENDPOINT_CALLS = EXEC / "interpreter_smoke_confirmatory_threshold_calls_v6/endpoint_calls_v6.tsv"
DEFAULT_MANIFEST = EXEC / "execution_candidate_manifest_v6.tsv"
DEFAULT_SCORES = ROOT / "project/results/cross_neo_immunogenicity_algorithm_compare_2026_05_10/cross_neo_immunogenicity_comparator_scores.tsv"
DEFAULT_OUT = ROOT / "project/results/cross_neo_wetlab_algorithm_comparison_v6_2026_05_10"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")

SCORE_SPECS = [
    ("BigMHC_IM", "bigmhc_im_score"),
    ("BigMHC_EL", "bigmhc_el_score"),
    ("CROSS_Neo_stress_guarded", "stress_guarded_discovery_score"),
    ("CROSS_Neo_BMA_v2", "bma_v2_discovery_score"),
    ("CROSS_Neo_finetuned_priority", "finetuned_experiment_priority_score"),
    ("CROSS_Neo_integrated_immunogenicity", "immunogenicity_discovery_score"),
    ("CROSS_Neo_claim_safe_integrated", "immunogenicity_claim_safe_score"),
]


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def numeric(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def metric_row(df: pd.DataFrame, label_col: str, score_col: str, name: str) -> dict[str, object]:
    d = df[[label_col, score_col]].copy()
    d[label_col] = numeric(d[label_col])
    d[score_col] = numeric(d[score_col])
    d = d.dropna()
    row: dict[str, object] = {
        "algorithm": name,
        "score_column": score_col,
        "n": int(len(d)),
        "hits": int(d[label_col].sum()) if len(d) else 0,
        "hit_rate": float(d[label_col].mean()) if len(d) else np.nan,
    }
    if len(d) and d[label_col].nunique() == 2:
        row["AUROC"] = float(roc_auc_score(d[label_col], d[score_col]))
        row["AUPRC"] = float(average_precision_score(d[label_col], d[score_col]))
    else:
        row["AUROC"] = np.nan
        row["AUPRC"] = np.nan
    ordered = d.sort_values(score_col, ascending=False)
    for k in [1, 2, 3, 5, 10, 12, 24]:
        top = ordered.head(k)
        row[f"top{k}_hits"] = int(top[label_col].sum()) if len(top) else 0
        row[f"top{k}_precision"] = float(top[label_col].mean()) if len(top) else np.nan
    return row


def classify_quartile_breakdown(df: pd.DataFrame, score_col: str, label_col: str) -> pd.DataFrame:
    d = df[["candidate_id", score_col, label_col]].copy()
    d[score_col] = numeric(d[score_col])
    d[label_col] = numeric(d[label_col])
    d = d.dropna()
    if d.empty:
        return pd.DataFrame()
    ranks = d[score_col].rank(method="first", ascending=True)
    d["score_bin"] = pd.qcut(ranks, q=min(4, len(d)), labels=False, duplicates="drop") + 1
    d["score_bin"] = d["score_bin"].map(lambda x: f"Q{int(x)}_low_to_high")
    out = (
        d.groupby("score_bin", observed=False)
        .agg(n=("candidate_id", "size"), hits=(label_col, "sum"), mean_score=(score_col, "mean"))
        .reset_index()
    )
    out["hit_rate"] = out["hits"] / out["n"]
    out["score_column"] = score_col
    return out[["score_column", "score_bin", "n", "hits", "hit_rate", "mean_score"]]


def join_inputs(calls: pd.DataFrame, manifest: pd.DataFrame, scores: pd.DataFrame, result_source: str) -> pd.DataFrame:
    keep_call_cols = [
        "candidate_id",
        "plate_v4_slot",
        "plate_v4_arm",
        "execution_blind_id",
        "peptide",
        "hla_allele_4digit",
        "impact_lane",
        "normalized_candidate_success",
        "candidate_call_source",
        "endpoint_countable",
        "endpoint_success",
    ]
    keep_call_cols = [c for c in keep_call_cols if c in calls.columns]
    out = calls[keep_call_cols].copy()
    out["wetlab_result_source"] = result_source
    out["wetlab_hit_binary"] = out.get("normalized_candidate_success", "").astype(str).eq("PASS").astype(int)
    out["wetlab_fail_binary"] = out.get("normalized_candidate_success", "").astype(str).eq("FAIL").astype(int)
    out["wetlab_countable"] = out.get("normalized_candidate_success", "").astype(str).isin(["PASS", "FAIL"])

    manifest_keep = [
        "candidate_id",
        "source_name",
        "leakage_risk_level",
        "claim_layer_after_unlock",
        "main_text_claim_tier",
        "order_priority_score",
        "order_batch",
        "claim_boundary_execution",
        "wildtype_peptide",
        "wt_status",
        "anchor_preserved_decoy_control",
    ]
    manifest_keep = [c for c in manifest_keep if c in manifest.columns]
    out = out.merge(manifest[manifest_keep].drop_duplicates("candidate_id"), on="candidate_id", how="left")

    score_cols = ["candidate_id"] + [col for _, col in SCORE_SPECS if col in scores.columns]
    extra = [
        "bigmhc_el_score",
        "fitness_foreignness_proxy",
        "tcr_recognition_score_norm",
        "md_control_score_norm",
        "immunogenicity_action",
        "immunogenicity_claim_blockers",
    ]
    score_cols += [c for c in extra if c in scores.columns and c not in score_cols]
    out = out.merge(scores[score_cols].drop_duplicates("candidate_id"), on="candidate_id", how="left")

    # Failure-mode decomposition. These are operational labels for assay triage.
    labels = []
    for _, row in out.iterrows():
        call = str(row.get("normalized_candidate_success", ""))
        arm = str(row.get("plate_v4_arm", ""))
        if call == "PASS":
            if arm == "A_clean_discovery":
                labels.append("hit_clean_antigen")
            elif arm == "B_mechanism_TCR_MD":
                labels.append("hit_tcr_md_mechanism")
            elif arm == "C_label_rescue":
                labels.append("hit_label_rescue")
            elif arm == "D_specificity_moat":
                labels.append("hit_specificity_control")
            elif arm == "E_positive_QC_control":
                labels.append("hit_positive_qc")
            elif arm == "F_model_boundary":
                labels.append("hit_model_boundary_resolution")
            else:
                labels.append("hit_other")
        elif call == "FAIL":
            if arm == "A_clean_discovery":
                labels.append("fail_presentation_or_specificity")
            elif arm == "B_mechanism_TCR_MD":
                labels.append("fail_tcr_or_pmhc_mechanism")
            elif arm == "C_label_rescue":
                labels.append("fail_label_rescue")
            elif arm == "D_specificity_moat":
                labels.append("fail_specificity_control")
            elif arm == "E_positive_QC_control":
                labels.append("fail_assay_qc")
            elif arm == "F_model_boundary":
                labels.append("fail_model_boundary")
            else:
                labels.append("fail_other")
        else:
            labels.append("pending_or_excluded")
    out["hit_fail_decomposition"] = labels
    return out


def build_metrics(joined: pd.DataFrame) -> pd.DataFrame:
    countable = joined[joined["wetlab_countable"]].copy()
    rows = []
    for name, col in SCORE_SPECS:
        if col in countable.columns:
            rows.append(metric_row(countable, "wetlab_hit_binary", col, name))
    return pd.DataFrame(rows).sort_values(["top5_precision", "AUPRC"], ascending=False)


def build_breakdown(joined: pd.DataFrame) -> pd.DataFrame:
    countable = joined[joined["wetlab_countable"]].copy()
    rows = []
    for name, col in SCORE_SPECS:
        if col in countable.columns:
            b = classify_quartile_breakdown(countable, col, "wetlab_hit_binary")
            if not b.empty:
                b.insert(0, "algorithm", name)
                rows.append(b)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def plot_outputs(out_dir: Path, joined: pd.DataFrame, metrics: pd.DataFrame) -> None:
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    top = metrics.head(7).iloc[::-1]
    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.barh(top["algorithm"], top["AUPRC"], color="#2f6f73", label="AUPRC")
    ax.scatter(top["AUROC"], top["algorithm"], color="#c59b3b", s=70, label="AUROC", zorder=3)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Score")
    ax.set_title("v6 assay hit/fail comparison by algorithm")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(fig_dir / "fig1_v6_algorithm_auprc_auroc.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    ks = [1, 2, 3, 5, 10, 12, 24]
    for _, row in metrics.head(5).iterrows():
        vals = [row.get(f"top{k}_precision", np.nan) for k in ks]
        ax.plot(ks, vals, marker="o", label=row["algorithm"])
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Top-k by score")
    ax.set_ylabel("Wetlab hit precision")
    ax.set_title("Top-k wetlab hit recovery")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig2_v6_topk_hit_curve.png", dpi=180)
    plt.close(fig)

    countable = joined[joined["wetlab_countable"]].copy()
    if not countable.empty:
        ordered = countable.sort_values("immunogenicity_discovery_score", ascending=False)
        colors = ordered["wetlab_hit_binary"].map({1: "#2f6f73", 0: "#d16b5f"})
        fig, ax = plt.subplots(figsize=(10, 4.8))
        ax.bar(np.arange(len(ordered)), ordered["immunogenicity_discovery_score"], color=colors)
        ax.set_xticks(np.arange(len(ordered)))
        ax.set_xticklabels(ordered["candidate_id"], rotation=70, ha="right", fontsize=8)
        ax.set_ylabel("Integrated immunogenicity score")
        ax.set_title("Candidate-level hit/fail decomposition")
        fig.tight_layout()
        fig.savefig(fig_dir / "fig3_v6_candidate_hit_fail_rank.png", dpi=180)
        plt.close(fig)


def table_html(df: pd.DataFrame, max_rows: int = 30) -> str:
    show = df.head(max_rows).copy()
    for col in show.columns:
        if pd.api.types.is_float_dtype(show[col]):
            show[col] = show[col].map(lambda x: "" if pd.isna(x) else f"{x:.3f}")
    return show.to_html(index=False, escape=True, classes="data")


def write_report(out_dir: Path, joined: pd.DataFrame, metrics: pd.DataFrame, endpoint_calls: pd.DataFrame, result_source: str) -> None:
    best = metrics.iloc[0] if len(metrics) else None
    hit_counts = joined["hit_fail_decomposition"].value_counts().rename_axis("call").reset_index(name="n")
    endpoint_view = endpoint_calls.copy() if endpoint_calls is not None and not endpoint_calls.empty else pd.DataFrame()
    lines = [
        "# CROSS-Neo v6 wetlab algorithm comparison",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        f"Result source: `{result_source}`",
        "",
        "## Bottom line",
        "",
    ]
    if best is not None:
        lines.append(
            f"- Best comparator in this result sheet: `{best['algorithm']}` "
            f"(AUPRC {best['AUPRC']:.3f}, AUROC {best['AUROC']:.3f}, top5 precision {best['top5_precision']:.3f})."
        )
    lines.extend(
        [
            f"- Countable candidates: {int(joined['wetlab_countable'].sum())}/{len(joined)}",
            f"- PASS calls: {int(joined['wetlab_hit_binary'].sum())}",
            "",
            "## Algorithm metrics",
            "",
            metrics.to_markdown(index=False),
            "",
            "## Hit/fail decomposition",
            "",
            hit_counts.to_markdown(index=False),
        ]
    )
    if not endpoint_view.empty:
        cols = [
            "plate_v4_arm",
            "success_count",
            "threshold_successes",
            "confirmatory_threshold_successes",
            "exploratory_unlock_status",
            "confirmatory_unlock_status",
            "claim_layer_after_unlock",
        ]
        cols = [c for c in cols if c in endpoint_view.columns]
        lines.extend(["", "## Endpoint unlocks", "", endpoint_view[cols].to_markdown(index=False)])
    lines.extend(
        [
            "",
            "## Claim boundary",
            "",
            "The default run uses a smoke-confirmatory result sheet, not real prospective wetlab data. When actual assay calls are entered, rerun this same script with `--candidate-calls <interpreter output>/candidate_calls_v6.tsv --result-source actual_wetlab_v6`.",
            "",
        ]
    )
    (out_dir / "V6_WETLAB_ALGORITHM_COMPARISON_REPORT_KR.md").write_text("\n".join(lines), encoding="utf-8")


def write_html(out_dir: Path, joined: pd.DataFrame, metrics: pd.DataFrame, breakdown: pd.DataFrame, endpoint_calls: pd.DataFrame, result_source: str) -> dict[str, object]:
    asset_dir = HUB / "assets/cross_neo_v6_wetlab_algorithm_comparison"
    asset_dir.mkdir(parents=True, exist_ok=True)
    for fig in (out_dir / "figures").glob("*.png"):
        shutil.copy2(fig, asset_dir / fig.name)

    best = metrics.iloc[0]
    top_candidates = joined.sort_values("immunogenicity_discovery_score", ascending=False)[
        [
            "candidate_id",
            "peptide",
            "hla_allele_4digit",
            "plate_v4_arm",
            "normalized_candidate_success",
            "hit_fail_decomposition",
            "bigmhc_im_score",
            "stress_guarded_discovery_score",
            "immunogenicity_discovery_score",
            "claim_boundary_execution",
        ]
    ]
    endpoint_cols = [
        "plate_v4_arm",
        "success_count",
        "threshold_successes",
        "confirmatory_threshold_successes",
        "exploratory_unlock_status",
        "confirmatory_unlock_status",
        "claim_layer_after_unlock",
    ]
    endpoint_cols = [c for c in endpoint_cols if c in endpoint_calls.columns]
    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>CROSS-Neo v6 Wetlab Algorithm Comparison</title>
  <style>
    :root {{ --bg:#101418; --panel:#151c22; --ink:#e9edf0; --muted:#9aa7af; --gold:#c59b3b; --line:#29343c; --teal:#69a7a2; }}
    body {{ margin:0; background:var(--bg); color:var(--ink); font:15px/1.55 system-ui, -apple-system, Segoe UI, sans-serif; }}
    header {{ padding:42px 5vw 30px; background:#0c1115; border-bottom:1px solid var(--line); }}
    .kicker {{ color:var(--gold); text-transform:uppercase; letter-spacing:.12em; font-weight:700; font-size:12px; }}
    h1 {{ margin:.3rem 0 .65rem; font-size:clamp(30px,4vw,54px); line-height:1.05; }}
    .lead {{ color:#d0d8dd; max-width:980px; font-size:18px; }}
    .stats {{ display:grid; grid-template-columns:repeat(5,minmax(120px,1fr)); gap:10px; margin-top:22px; }}
    .stat {{ border:1px solid var(--line); background:var(--panel); padding:13px 14px; border-radius:8px; }}
    .stat b {{ display:block; font-size:23px; }}
    .stat span {{ color:var(--muted); font-size:12px; }}
    main {{ padding:30px 5vw 60px; display:grid; grid-template-columns:260px minmax(0,1fr); gap:28px; }}
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
  <div class="kicker">CROSS-Neo · v6 assay interpreter · algorithm comparison</div>
  <h1>96-well hit/fail comparison scaffold</h1>
  <p class="lead">The v6 interpreter output is joined to BigMHC, CROSS-Neo and integrated immunogenicity scores. This default view uses the confirmatory-threshold smoke sheet; the same script accepts actual wetlab calls when available.</p>
  <div class="stats">
    <div class="stat"><b>{int(joined['wetlab_countable'].sum())}/{len(joined)}</b><span>countable candidates</span></div>
    <div class="stat"><b>{int(joined['wetlab_hit_binary'].sum())}</b><span>PASS calls</span></div>
    <div class="stat"><b>{html.escape(str(best['algorithm']))}</b><span>best comparator</span></div>
    <div class="stat"><b>{best['AUPRC']:.3f}</b><span>best AUPRC</span></div>
    <div class="stat"><b>{html.escape(result_source)}</b><span>result source</span></div>
  </div>
</header>
<main>
<nav>
  <a href="#summary">01 Summary</a>
  <a href="#metrics">02 Algorithm Metrics</a>
  <a href="#endpoint">03 Endpoint Unlocks</a>
  <a href="#breakdown">04 Hit/Fail Breakdown</a>
  <a href="#candidates">05 Candidate Table</a>
  <a href="#handoff">06 Actual Wetlab Handoff</a>
</nav>
<article>
<section id="summary">
  <h2><span class="num">01</span>Summary</h2>
  <p class="note">This is the exact product surface needed after the assay: input candidate calls, output endpoint unlocks, algorithm comparison, and failure-mode labels. Current result source is <code>{html.escape(result_source)}</code>.</p>
  <div class="grid">
    <img src="assets/cross_neo_v6_wetlab_algorithm_comparison/fig1_v6_algorithm_auprc_auroc.png" alt="algorithm AUPRC AUROC">
    <img src="assets/cross_neo_v6_wetlab_algorithm_comparison/fig2_v6_topk_hit_curve.png" alt="top-k hit curve">
  </div>
  <img src="assets/cross_neo_v6_wetlab_algorithm_comparison/fig3_v6_candidate_hit_fail_rank.png" alt="candidate hit fail rank">
</section>
<section id="metrics">
  <h2><span class="num">02</span>Algorithm Metrics</h2>
  {table_html(metrics, 20)}
</section>
<section id="endpoint">
  <h2><span class="num">03</span>Endpoint Unlocks</h2>
  {table_html(endpoint_calls[endpoint_cols], 10) if endpoint_cols else '<p>No endpoint calls provided.</p>'}
</section>
<section id="breakdown">
  <h2><span class="num">04</span>Hit/Fail Breakdown</h2>
  {table_html(breakdown, 40)}
</section>
<section id="candidates">
  <h2><span class="num">05</span>Candidate-Level Calls</h2>
  {table_html(top_candidates, 30)}
</section>
<section id="handoff">
  <h2><span class="num">06</span>Actual Wetlab Handoff</h2>
  <p>When real results are entered, run:</p>
  <p><code>python scripts/interpret_cross_neo_assay_results_v6.py --results .../candidate_result_entry_v6.tsv --endpoint-plan .../preregistered_endpoint_plan_v5.tsv --out-dir .../interpreted_results_v6</code></p>
  <p>Then run:</p>
  <p><code>python scripts/build_cross_neo_wetlab_algorithm_comparison_v6.py --candidate-calls .../interpreted_results_v6/candidate_calls_v6.tsv --endpoint-calls .../interpreted_results_v6/endpoint_calls_v6.tsv --result-source actual_wetlab_v6</code></p>
</section>
</article>
</main>
</body>
</html>
"""
    html_path = HUB / "cross_neo_v6_wetlab_algorithm_comparison.html"
    html_path.write_text(html_text, encoding="utf-8")

    live_ok = False
    warnings: list[str] = []
    try:
        live_asset_dir = LIVE_HUB / "assets/cross_neo_v6_wetlab_algorithm_comparison"
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
    parser.add_argument("--candidate-calls", type=Path, default=DEFAULT_CALLS)
    parser.add_argument("--endpoint-calls", type=Path, default=DEFAULT_ENDPOINT_CALLS)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--scores", type=Path, default=DEFAULT_SCORES)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--result-source", default="smoke_confirmatory_threshold_v6")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    calls = read_tsv(args.candidate_calls)
    endpoint_calls = read_tsv(args.endpoint_calls) if args.endpoint_calls.exists() else pd.DataFrame()
    manifest = read_tsv(args.manifest)
    scores = read_tsv(args.scores)

    joined = join_inputs(calls, manifest, scores, args.result_source)
    metrics = build_metrics(joined)
    breakdown = build_breakdown(joined)

    joined.to_csv(args.out_dir / "v6_wetlab_algorithm_joined_calls.tsv", sep="\t", index=False)
    metrics.to_csv(args.out_dir / "v6_wetlab_algorithm_comparison_metrics.tsv", sep="\t", index=False)
    breakdown.to_csv(args.out_dir / "v6_wetlab_algorithm_score_bin_breakdown.tsv", sep="\t", index=False)
    endpoint_calls.to_csv(args.out_dir / "v6_endpoint_calls_used.tsv", sep="\t", index=False)

    plot_outputs(args.out_dir, joined, metrics)
    write_report(args.out_dir, joined, metrics, endpoint_calls, args.result_source)
    html_info = write_html(args.out_dir, joined, metrics, breakdown, endpoint_calls, args.result_source)

    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "result_source": args.result_source,
        "candidate_calls_input": str(args.candidate_calls),
        "endpoint_calls_input": str(args.endpoint_calls),
        "score_input": str(args.scores),
        "n_candidates": int(len(joined)),
        "n_countable": int(joined["wetlab_countable"].sum()),
        "n_pass": int(joined["wetlab_hit_binary"].sum()),
        "best_algorithm": metrics.iloc[0].to_dict() if len(metrics) else {},
        "hit_fail_decomposition_counts": joined["hit_fail_decomposition"].value_counts().to_dict(),
        "output_dir": str(args.out_dir),
        **html_info,
        "claim_boundary": "default is smoke/dry-run unless result_source is actual_wetlab_v6 and inputs are filled wetlab calls",
        "output_files": [
            str(args.out_dir / "v6_wetlab_algorithm_joined_calls.tsv"),
            str(args.out_dir / "v6_wetlab_algorithm_comparison_metrics.tsv"),
            str(args.out_dir / "v6_wetlab_algorithm_score_bin_breakdown.tsv"),
            str(args.out_dir / "V6_WETLAB_ALGORITHM_COMPARISON_REPORT_KR.md"),
        ],
    }
    (args.out_dir / "v6_wetlab_algorithm_comparison_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
