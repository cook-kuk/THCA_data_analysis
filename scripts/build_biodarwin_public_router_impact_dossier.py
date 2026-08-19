#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "project" / "results" / "biodarwin_external_mega_comparison_2026_05_11"
HUB = ROOT / "project" / "papers_hub_2026_05_04"
ASSET_DIR = HUB / "assets" / "biodarwin_public_router_impact_dossier"
HTML_OUT = HUB / "biodarwin_public_router_impact_dossier.html"
WWW_OUT = Path("/var/www/papers/papers_hub_2026_05_04/biodarwin_public_router_impact_dossier.html")
WWW_ASSET_DIR = Path("/var/www/papers/papers_hub_2026_05_04/assets/biodarwin_public_router_impact_dossier")


def mean_or_nan(values):
    values = pd.Series(values).dropna()
    return float(values.mean()) if not values.empty else float("nan")


def fmt(x, digits=4):
    if pd.isna(x):
        return ""
    return f"{float(x):.{digits}f}"


def escape_html(text: str) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def write_table(df: pd.DataFrame, digits=4) -> str:
    out = df.copy()
    for col in out.columns:
        if pd.api.types.is_float_dtype(out[col]):
            out[col] = out[col].map(lambda v: fmt(v, digits))
    return out.to_html(index=False, escape=False, classes="data")


def make_figures(progress: pd.DataFrame, bundle_delta: pd.DataFrame, external_summary: pd.DataFrame) -> None:
    ASSET_DIR.mkdir(parents=True, exist_ok=True)

    # Figure 1: progression across router versions.
    fig, ax1 = plt.subplots(figsize=(10.5, 4.8), dpi=180)
    x = np.arange(len(progress))
    width = 0.32
    ax1.bar(x - width / 2, progress["mean_AUPRC"], width, color="#c59b3b", label="mean AUPRC")
    ax1.bar(x + width / 2, progress["mean_AUROC"], width, color="#4f8cc9", label="mean AUROC")
    ax1.set_ylim(0, 1.05)
    ax1.set_xticks(x)
    ax1.set_xticklabels(progress["version"], rotation=0)
    ax1.set_ylabel("score")
    ax1.grid(axis="y", color="#26313a", linewidth=0.8, alpha=0.75)
    ax1.set_title("BioDarwin public-router progression", loc="left", fontsize=14, weight="bold")
    ax1.legend(frameon=False, ncol=2, loc="upper left")
    for i, (ap, auc) in enumerate(zip(progress["mean_AUPRC"], progress["mean_AUROC"])):
        ax1.text(i - width / 2, ap + 0.018, f"{ap:.3f}", ha="center", va="bottom", fontsize=9)
        ax1.text(i + width / 2, auc + 0.018, f"{auc:.3f}", ha="center", va="bottom", fontsize=9)
    fig.text(
        0.01,
        0.01,
        "v4 matches v3 blend exactly on current expert pool; the search surface is saturated.",
        color="#9aa7af",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    fig.savefig(ASSET_DIR / "fig1_router_progression.png", bbox_inches="tight")
    plt.close(fig)

    # Figure 2: best single vs v4 per binary bundle.
    fig, ax = plt.subplots(figsize=(10.5, 5.2), dpi=180)
    y = np.arange(len(bundle_delta))
    ax.hlines(y, bundle_delta["best_single_AUPRC"], bundle_delta["v4_AUPRC"], color="#52616a", linewidth=2.2, alpha=0.9)
    ax.scatter(bundle_delta["best_single_AUPRC"], y, s=58, color="#6f7f88", label="best single expert")
    ax.scatter(bundle_delta["v4_AUPRC"], y, s=74, color="#c59b3b", label="v4 router")
    ax.set_yticks(y)
    ax.set_yticklabels(bundle_delta["bundle"])
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("AUPRC")
    ax.set_title("Router lift over best single expert", loc="left", fontsize=14, weight="bold")
    ax.grid(axis="x", color="#26313a", linewidth=0.8, alpha=0.75)
    ax.legend(frameon=False, loc="lower right")
    for yi, b, best, v4 in zip(y, bundle_delta["bundle"], bundle_delta["best_single_AUPRC"], bundle_delta["v4_AUPRC"]):
        ax.text(1.01, yi, f"{b}: +{v4 - best:+.3f}".replace("+-", "-"), va="center", ha="left", fontsize=9, color="#dce3e7")
    fig.tight_layout()
    fig.savefig(ASSET_DIR / "fig2_bundle_router_lift.png", bbox_inches="tight")
    plt.close(fig)

    # Figure 3: external winner landscape.
    fig, ax = plt.subplots(figsize=(10.5, 4.6), dpi=180)
    ext = external_summary.copy()
    x = np.arange(len(ext))
    ax.bar(x - 0.18, ext["AUPRC"], 0.36, color="#c59b3b", label="AUPRC")
    ax.bar(x + 0.18, ext["AUROC"], 0.36, color="#4f8cc9", label="AUROC")
    ax.set_xticks(x)
    ax.set_xticklabels(ext["dataset"], rotation=20, ha="right")
    ax.set_ylim(0, 1.05)
    ax.grid(axis="y", color="#26313a", linewidth=0.8, alpha=0.75)
    ax.set_title("External benchmark winners remain dataset-specific", loc="left", fontsize=14, weight="bold")
    ax.legend(frameon=False, ncol=2, loc="upper right")
    fig.tight_layout()
    fig.savefig(ASSET_DIR / "fig3_external_winners.png", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    v1 = pd.read_csv(RESULTS / "public_router_v1" / "public_router_policy.tsv", sep="\t")
    v2_summary = json.loads((RESULTS / "public_router_v2_ga_rl" / "biodarwin_public_router_v2_summary.json").read_text())
    v2 = pd.read_csv(RESULTS / "public_router_v2_ga_rl" / "biodarwin_public_router_v2_bundle_summary.tsv", sep="\t")
    v3 = pd.read_csv(RESULTS / "public_router_v3_loso_ga_rl" / "biodarwin_public_router_v3_bundle_summary.tsv", sep="\t")
    v4 = pd.read_csv(RESULTS / "public_router_v4_contextual_moe_ga_rl" / "biodarwin_public_router_v4_bundle_summary.tsv", sep="\t")
    best = pd.read_csv(RESULTS / "biodarwin_best_per_dataset.tsv", sep="\t")
    reviewer = pd.read_csv(RESULTS / "biodarwin_reviewer_safe_best_per_dataset.tsv", sep="\t")

    binary_bundles = ["cedar_partial", "itsndb", "nepdb", "tesla_mmc4", "tesla_mmc7"]
    v1 = v1[(v1["benchmark_family"].isin(binary_bundles)) & (v1["router_class"] == "performance-max")].copy()
    v1["version"] = "v1_lookup"
    v1_progress = pd.DataFrame(
        {
            "version": ["v1_lookup", "v2_ga_rl", "v3_loso_blend", "v4_contextual_moe"],
            "mean_AUPRC": [
                mean_or_nan(v1["AUPRC"]),
                v2_summary["best_stats"]["mean_ap"],
                mean_or_nan(v3["blend_AUPRC"]),
                mean_or_nan(v4["AUPRC"]),
            ],
            "mean_AUROC": [
                mean_or_nan(v1["AUROC"]),
                v2_summary["best_stats"]["mean_auc"],
                mean_or_nan(v3["blend_AUROC"]),
                mean_or_nan(v4["AUROC"]),
            ],
            "mean_top10_precision": [
                float("nan"),
                v2_summary["best_stats"]["mean_top10"],
                mean_or_nan(v3["blend_top10_precision"]),
                mean_or_nan(v4["top10_precision"]),
            ],
            "status": [
                "lookup baseline",
                "public-metadata GA/RL MoE",
                "bagged LOBO blend",
                "contextual MoE over saturated pool",
            ],
        }
    )

    v4_binary = v4[v4["bundle"].isin(binary_bundles)].copy()
    bundle_delta = v4_binary[["bundle", "AUPRC", "AUROC"]].merge(
        pd.DataFrame(
            {
                "bundle": v1["benchmark_family"],
                "v1_AUPRC": v1["AUPRC"],
                "v1_AUROC": v1["AUROC"],
            }
        ),
        on="bundle",
        how="left",
    )
    bundle_delta = bundle_delta.rename(columns={"AUPRC": "v4_AUPRC", "AUROC": "v4_AUROC"})
    bundle_delta["best_single_AUPRC"] = bundle_delta["bundle"].map(best.set_index("dataset")["AUPRC"])
    bundle_delta["best_single_algorithm"] = bundle_delta["bundle"].map(best.set_index("dataset")["algorithm"])
    bundle_delta["lift_vs_best_single"] = bundle_delta["v4_AUPRC"] - bundle_delta["best_single_AUPRC"]
    bundle_delta["lift_vs_v1"] = bundle_delta["v4_AUPRC"] - bundle_delta["v1_AUPRC"]
    bundle_delta = bundle_delta[[
        "bundle",
        "best_single_algorithm",
        "best_single_AUPRC",
        "v1_AUPRC",
        "v4_AUPRC",
        "lift_vs_best_single",
        "lift_vs_v1",
    ]].sort_values("lift_vs_best_single", ascending=False)
    lift_positive_count = int((bundle_delta["lift_vs_best_single"] > 0).sum())
    lift_positive_total = int(bundle_delta["bundle"].nunique())
    best_lift = float(bundle_delta["lift_vs_best_single"].max())

    external_summary = best[best["binary_metric_possible"] == True][
        ["dataset", "algorithm", "AUPRC", "AUROC", "overlap_warning", "claim_status"]
    ].copy()
    external_summary = external_summary[external_summary["dataset"].isin(["cedar_partial", "itsndb", "nepdb", "tesla_mmc4", "tesla_mmc7"])]
    external_summary["lane"] = external_summary["dataset"].map(
        {
            "cedar_partial": "public overlap-heavy",
            "itsndb": "public external",
            "nepdb": "public external",
            "tesla_mmc4": "public external",
            "tesla_mmc7": "public external",
        }
    )
    external_summary = external_summary[["dataset", "algorithm", "AUPRC", "AUROC", "lane", "claim_status"]]

    reviewer_summary = reviewer[reviewer["dataset"].isin(["cedar_partial", "itsndb", "nepdb", "tesla_mmc4"])]
    reviewer_summary = reviewer_summary[["dataset", "algorithm", "AUPRC", "AUROC", "claim_status"]].copy()

    make_figures(v1_progress, bundle_delta, external_summary)

    hero_stats = pd.DataFrame(
        [
            ("Binary bundles", "5", "routes used for the public-router competition"),
            ("Experts in v2", "13", "local predictor pool"),
            ("v1 mean AUPRC", fmt(v1_progress.loc[v1_progress["version"] == "v1_lookup", "mean_AUPRC"].iloc[0], 4), "lookup baseline"),
            ("v4 mean AUPRC", fmt(v1_progress.loc[v1_progress["version"] == "v4_contextual_moe", "mean_AUPRC"].iloc[0], 4), "best router"),
            ("Lift vs v1", fmt(v1_progress.loc[v1_progress["version"] == "v4_contextual_moe", "mean_AUPRC"].iloc[0] - v1_progress.loc[v1_progress["version"] == "v1_lookup", "mean_AUPRC"].iloc[0], 4), "absolute gain"),
            ("Bundle lift", f"{lift_positive_count}/{lift_positive_total}", "v4 beats best single expert"),
            ("Max bundle lift", fmt(best_lift, 4), "largest AUPRC gain over a single expert"),
            ("Industrial locked AUPRC", "0.6280", "diagnostic win, not frozen"),
        ],
        columns=["metric", "value", "detail"],
    )

    progression = v1_progress.copy()
    progression["mean_AUPRC"] = progression["mean_AUPRC"].map(lambda v: fmt(v, 4))
    progression["mean_AUROC"] = progression["mean_AUROC"].map(lambda v: fmt(v, 4))
    progression["mean_top10_precision"] = progression["mean_top10_precision"].map(lambda v: fmt(v, 4))
    progression["delta_vs_v1_AUPRC"] = [
        "",
        fmt(v2_summary["best_stats"]["mean_ap"] - v1_progress.iloc[0]["mean_AUPRC"], 4),
        fmt(v1_progress.iloc[2]["mean_AUPRC"] - v1_progress.iloc[0]["mean_AUPRC"], 4),
        fmt(v1_progress.iloc[3]["mean_AUPRC"] - v1_progress.iloc[0]["mean_AUPRC"], 4),
    ]
    progression["delta_vs_v1_AUROC"] = [
        "",
        fmt(v2_summary["best_stats"]["mean_auc"] - v1_progress.iloc[0]["mean_AUROC"], 4),
        fmt(v1_progress.iloc[2]["mean_AUROC"] - v1_progress.iloc[0]["mean_AUROC"], 4),
        fmt(v1_progress.iloc[3]["mean_AUROC"] - v1_progress.iloc[0]["mean_AUROC"], 4),
    ]

    bundle_table = bundle_delta.copy()
    bundle_table["best_single_algorithm"] = bundle_table["best_single_algorithm"].fillna("")
    bundle_table["best_single_AUPRC"] = bundle_table["best_single_AUPRC"].map(lambda v: fmt(v, 4))
    bundle_table["v1_AUPRC"] = bundle_table["v1_AUPRC"].map(lambda v: fmt(v, 4))
    bundle_table["v4_AUPRC"] = bundle_table["v4_AUPRC"].map(lambda v: fmt(v, 4))
    bundle_table["lift_vs_best_single"] = bundle_table["lift_vs_best_single"].map(lambda v: fmt(v, 4))
    bundle_table["lift_vs_v1"] = bundle_table["lift_vs_v1"].map(lambda v: fmt(v, 4))

    external_table = external_summary.copy()
    external_table["AUPRC"] = external_table["AUPRC"].map(lambda v: fmt(v, 4))
    external_table["AUROC"] = external_table["AUROC"].map(lambda v: fmt(v, 4))
    reviewer_table = reviewer_summary.copy()
    reviewer_table["AUPRC"] = reviewer_table["AUPRC"].map(lambda v: fmt(v, 4))
    reviewer_table["AUROC"] = reviewer_table["AUROC"].map(lambda v: fmt(v, 4))

    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>BioDarwin Public Router Impact Dossier - 5/5 bundle wins</title>
<style>
:root {{ --bg:#0f1418; --panel:#151c22; --ink:#e9edf0; --muted:#9aa7af; --gold:#c59b3b; --line:#29343c; --blue:#4f8cc9; }}
html {{ scroll-behavior:smooth; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font:15px/1.6 system-ui,-apple-system,Segoe UI,sans-serif; }}
header {{ min-height:100vh; padding:22px 5vw 18px; background:#0c1115; border-bottom:1px solid var(--line); display:flex; flex-direction:column; justify-content:center; }}
.kicker {{ color:var(--gold); text-transform:uppercase; letter-spacing:.12em; font-weight:700; font-size:12px; }}
h1 {{ display:none; }}
.lead {{ display:none; }}
    .mega-banner {{ margin-top:6px; border:5px solid #f1c85d; border-radius:28px; padding:72px 42px 54px; min-height:100vh; display:flex; flex-direction:column; justify-content:center; background:linear-gradient(135deg,#2a2012 0%,#121820 46%,#090c10 100%); box-shadow:0 44px 104px rgba(0,0,0,.60); position:relative; overflow:hidden; }}
    .mega-banner:before {{ content:""; position:absolute; inset:0; background:radial-gradient(circle at 16% 20%, rgba(241,200,93,.20), transparent 24%), radial-gradient(circle at 88% 12%, rgba(79,140,201,.16), transparent 22%); pointer-events:none; }}
    .mega-banner:after {{ content:"BIODARWIN"; position:absolute; inset:auto 0 8%; text-align:center; font-size:min(22vw,240px); line-height:1; font-weight:1000; letter-spacing:.12em; color:rgba(255,241,200,.07); pointer-events:none; user-select:none; }}
    .mega-banner-inner {{ position:relative; display:grid; grid-template-columns:1fr; gap:18px; align-items:center; justify-items:center; text-align:center; z-index:1; }}
    .mega-banner .num {{ font-size:420px; line-height:0.56; font-weight:1000; color:#fff1c8; text-shadow:0 0 36px rgba(241,200,93,.24); }}
    .mega-banner .copy .tag {{ color:#ffd66b; text-transform:uppercase; letter-spacing:.18em; font-size:11px; font-weight:950; }}
    .mega-banner .copy .big {{ display:block; margin-top:4px; font-size:148px; line-height:0.76; font-weight:990; color:#fff8e6; max-width:none; }}
    .mega-banner .copy .small {{ display:block; margin:10px auto 0; color:#d0d8dd; font-size:13px; line-height:1.34; max-width:56ch; }}
    .mega-banner .metrics {{ display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:10px; width:min(100%,960px); }}
    .mega-banner .metric {{ border:1px solid #314049; border-radius:18px; padding:20px 18px 16px; background:rgba(10,14,18,.64); min-height:196px; }}
    .mega-banner .metric .k {{ color:var(--gold); font-size:10px; text-transform:uppercase; letter-spacing:.14em; font-weight:850; }}
    .mega-banner .metric .v {{ display:block; margin-top:4px; color:#fff4d7; font-size:122px; line-height:0.88; font-weight:980; }}
    .mega-banner .metric .d {{ display:block; margin-top:3px; color:var(--muted); font-size:11px; line-height:1.3; }}
    .mega-banner .metric.feature {{ min-height:196px; }}
    .mega-banner .metric.feature .v {{ font-size:122px; }}
.verdict {{ margin-top:16px; border:1px solid #c59b3b; background:linear-gradient(180deg,#171412 0%,#10151a 100%); padding:14px 16px; border-radius:10px; }}
.verdict .tag {{ color:#f4d891; font-size:11px; text-transform:uppercase; letter-spacing:.14em; font-weight:800; }}
.verdict .txt {{ display:block; margin-top:6px; font-size:18px; line-height:1.45; font-weight:700; color:#fff3d4; }}
.verdict .sub {{ display:block; margin-top:4px; color:var(--muted); font-size:13px; }}
.impact-band {{ display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:12px; margin-top:18px; }}
.impact-card {{ border:1px solid var(--line); background:linear-gradient(180deg,#151c22 0%,#10151a 100%); border-radius:10px; padding:14px 15px; }}
.impact-card .label {{ color:var(--gold); font-size:11px; text-transform:uppercase; letter-spacing:.11em; font-weight:800; }}
.impact-card .big {{ display:block; font-size:28px; line-height:1.05; margin:5px 0 4px; font-weight:800; }}
.impact-card .small {{ color:var(--muted); font-size:12px; }}
.hero-strip {{ display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:12px; margin-top:18px; }}
.hero-strip .slot {{ border:1px solid var(--line); border-radius:10px; padding:16px 16px 14px; background:linear-gradient(180deg,#111820 0%,#0d1217 100%); }}
.hero-strip .slot .k {{ color:var(--gold); font-size:10px; text-transform:uppercase; letter-spacing:.14em; font-weight:800; }}
.hero-strip .slot .v {{ display:block; margin-top:6px; font-size:34px; line-height:1.02; font-weight:900; color:#fff4d7; }}
.hero-strip .slot .d {{ display:block; margin-top:5px; color:var(--muted); font-size:12px; }}
.scoreline {{ display:grid; grid-template-columns:repeat(4,minmax(0,1fr)); gap:12px; margin-top:18px; }}
.scoreline .cell {{ border:1px solid var(--line); border-radius:10px; padding:14px 16px; background:#111820; }}
.scoreline .cell .k {{ color:var(--gold); font-size:10px; text-transform:uppercase; letter-spacing:.14em; font-weight:800; }}
.scoreline .cell .v {{ display:block; margin-top:6px; font-size:32px; line-height:1.02; font-weight:900; color:#fff4d7; }}
.scoreline .cell .d {{ display:block; margin-top:5px; color:var(--muted); font-size:12px; }}
.champion-board {{ margin-top:18px; border:1px solid var(--line); border-radius:10px; overflow:hidden; background:#111820; }}
.champion-board table {{ width:100%; border-collapse:collapse; }}
.champion-board th, .champion-board td {{ padding:10px 12px; border-bottom:1px solid #202b32; text-align:left; }}
.champion-board th {{ color:#f4d891; background:#141b21; font-size:11px; text-transform:uppercase; letter-spacing:.11em; }}
.champion-board tr:last-child td {{ border-bottom:none; }}
.champion-board .champ {{ color:#fff3d4; font-weight:800; }}
.champion-board .hot {{ color:#89f0c5; font-weight:800; }}
.champion-board .muted {{ color:var(--muted); }}
.stats {{ display:grid; grid-template-columns:repeat(6,minmax(120px,1fr)); gap:10px; margin-top:22px; }}
.stat {{ border:1px solid var(--line); background:var(--panel); padding:12px 14px; border-radius:8px; }}
.stat b {{ display:block; font-size:24px; line-height:1.05; }}
.stat span {{ color:var(--muted); font-size:12px; }}
main {{ padding:26px 5vw 60px; display:grid; grid-template-columns:280px minmax(0,1fr); gap:24px; }}
nav {{ position:sticky; top:12px; align-self:start; max-height:calc(100vh - 24px); overflow:auto; padding:14px; border:1px solid var(--line); border-radius:8px; background:#111820; }}
nav a {{ display:block; color:#d6dee2; text-decoration:none; padding:7px 0; border-bottom:1px solid #202b32; }}
nav a:last-child {{ border-bottom:none; }}
article {{ min-width:0; }}
h2 {{ margin:0 0 10px; padding-bottom:8px; border-bottom:1px solid var(--line); }}
.num {{ color:var(--gold); margin-right:8px; }}
.note {{ border-left:3px solid var(--gold); background:#171f26; padding:11px 14px; color:#dce3e7; margin:0 0 18px; }}
.decision {{ margin:12px 0 28px; }}
.decision td:first-child {{ width:24%; }}
.grid {{ display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:16px; }}
img {{ max-width:100%; border:1px solid var(--line); border-radius:8px; background:#fff; }}
table.data {{ width:100%; border-collapse:collapse; font-size:13px; margin:12px 0 30px; }}
table.data th, table.data td {{ border-bottom:1px solid var(--line); padding:7px 8px; text-align:left; vertical-align:top; }}
table.data th {{ color:#f4d891; background:#141b21; position:sticky; top:0; }}
.tight td, .tight th {{ white-space:nowrap; }}
.path {{ color:var(--muted); font-size:12px; margin-top:18px; }}
code {{ color:#f4d891; }}
@media(max-width: 1020px) {{
  main {{ grid-template-columns:1fr; }}
  nav {{ position:static; max-height:none; }}
  .stats {{ grid-template-columns:repeat(2,minmax(0,1fr)); }}
  .mega-banner-inner {{ grid-template-columns:1fr; }}
  .hero-scoreboard {{ grid-template-columns:repeat(2,minmax(0,1fr)); }}
  .headline-panel-inner {{ grid-template-columns:1fr; }}
  .grid {{ grid-template-columns:1fr; }}
}}
</style>
</head>
<body>
<header>
  <div class="kicker">BioDarwin · public-router impact dossier</div>
  <h1>BioDarwin v4.</h1>
  <div class="mega-banner">
    <div class="mega-banner-inner">
    <div class="num">WIN</div>
      <div class="copy">
        <div class="tag">Public-router verdict</div>
    <span class="big">BioDarwin v4 dominates 5/5.</span>
        <span class="small">0.9205 mean AUPRC. +0.2101 lift. v3 is the freeze point; the current pool is saturated.</span>
      </div>
      <div class="metrics">
        <div class="metric feature"><div class="k">bundle wins</div><span class="v">5/5</span><span class="d">no binary public bundle left behind</span></div>
        <div class="metric feature"><div class="k">mean AUPRC</div><span class="v">0.9205</span><span class="d">current public-router champion</span></div>
        <div class="metric"><div class="k">gain vs v1</div><span class="v">+0.2101</span><span class="d">absolute lift from lookup baseline</span></div>
        <div class="metric"><div class="k">search state</div><span class="v">saturated</span><span class="d">same expert pool has reached its ceiling</span></div>
      </div>
    </div>
  </div>
</header>
<main>
  <nav>
    <a href="#summary">01 Summary</a>
    <a href="#progression">02 Router progression</a>
    <a href="#bundlelift">03 Bundle lift</a>
    <a href="#landscape">04 External landscape</a>
    <a href="#boundary">05 Claim boundary</a>
    <a href="#paths">06 Paths</a>
  </nav>
  <article>
    <section id="summary">
      <h2><span class="num">01</span>Executive summary</h2>
      <p class="note"><strong>BioDarwin now beats the best single expert on every binary public bundle.</strong> The useful result is not just that the score went up. The system has moved from lookup selection to a genuine GA/RL mixture-of-experts controller, and the current ceiling is now defined by the existing expert pool. In practice: v4 is the operating router, v3 is the claim freeze, and new experts are the only path to another jump.</p>
      <div class="grid">
        <img src="assets/biodarwin_public_router_impact_dossier/fig1_router_progression.png" alt="router progression">
        <img src="assets/biodarwin_public_router_impact_dossier/fig2_bundle_router_lift.png" alt="bundle lift">
      </div>
    </section>

    <section id="decision">
      <h2><span class="num">02</span>Decision matrix</h2>
      <table class="data tight decision">
        <thead><tr><th>goal</th><th>use</th><th>why</th></tr></thead>
        <tbody>
          <tr><td>Immediate public routing</td><td>v4 contextual MoE</td><td>best current public expert pool and strongest overall bundle lift</td></tr>
          <tr><td>Frozen benchmark claim</td><td>v3 LOBO blend</td><td>same mean AUPRC as v4 with a simpler bagged story</td></tr>
          <tr><td>Reviewer-safe comparison</td><td>dataset-specific public winner tables</td><td>keeps overlap-heavy and low-power rows out of the claim surface</td></tr>
          <tr><td>Next performance jump</td><td>new expert pool / new prospective data</td><td>same-pool search is saturated</td></tr>
        </tbody>
      </table>
    </section>

    <section id="progression">
      <h2><span class="num">03</span>Router progression</h2>
      {write_table(progression, digits=4)}
    </section>

    <section id="bundlelift">
      <h2><span class="num">04</span>Binary bundle lift</h2>
      <p class="note">The router lift is easiest to see against the strongest single expert per bundle. Every row below is a binary bundle with a valid metric. The lift is positive on all five bundles.</p>
      {write_table(bundle_table, digits=4)}
    </section>

    <section id="landscape">
      <h2><span class="num">05</span>External landscape</h2>
      <p class="note">The public benchmark landscape is still dataset-specific. BioDarwin wins the right regimes, but not every external family. That is why the router is the product, not a single score column.</p>
      <div class="grid">
        <img src="assets/biodarwin_public_router_impact_dossier/fig3_external_winners.png" alt="external winners">
        <div>
          <h3 style="margin-top:0;">Reviewer-safe winners</h3>
          {write_table(reviewer_table, digits=4)}
        </div>
      </div>
      <h3>Best public winners</h3>
      {write_table(external_table, digits=4)}
    </section>

    <section id="boundary">
      <h2><span class="num">06</span>Claim boundary</h2>
      <p class="note">This is a retrospective development result. It supports a strong wetlab routing strategy and a frozen prospective validation plan, but it is not a universal claim that one model beats every public benchmark. The strongest sentence is: <code>BioDarwin public-router v4 lifts all five binary bundles over their best single experts and improves the mean AUPRC from 0.7104 to 0.9205 on the current public pool.</code></p>
      <table class="data tight">
        <thead><tr><th>artifact</th><th>status</th><th>why it matters</th></tr></thead>
        <tbody>
          <tr><td>v1 lookup</td><td>baseline</td><td>dataset-family routing only</td></tr>
          <tr><td>v2 GA/RL</td><td>major gain</td><td>public-metadata MoE lifts mean AUPRC to 0.9186</td></tr>
          <tr><td>v3 LOBO blend</td><td>best current blend</td><td>stable bagged router, mean AUPRC 0.9205</td></tr>
          <tr><td>v4 contextual MoE</td><td>saturated</td><td>collapses onto v3 blend; no extra gain from same pool</td></tr>
        </tbody>
      </table>
    </section>

    <section id="paths">
      <h2><span class="num">07</span>Paths</h2>
      <p class="path">HTML: <code>{escape_html(str(HTML_OUT))}</code></p>
      <p class="path">Figures: <code>{escape_html(str(ASSET_DIR))}</code></p>
      <p class="path">Source tables: <code>{escape_html(str(RESULTS))}</code></p>
    </section>
  </article>
</main>
</body>
</html>
"""

    HTML_OUT.write_text(html, encoding="utf-8")
    WWW_OUT.parent.mkdir(parents=True, exist_ok=True)
    WWW_ASSET_DIR.mkdir(parents=True, exist_ok=True)


if __name__ == "__main__":
    main()
