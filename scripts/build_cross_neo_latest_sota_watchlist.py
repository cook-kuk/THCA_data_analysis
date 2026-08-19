#!/usr/bin/env python3
"""Build a latest-SOTA watchlist for neoantigen/immunogenicity models.

This compares current local CROSS-Neo / BAR-Neo results against recent
published SOTA methods. It keeps apples-to-apples local comparisons separate
from literature-only methods whose datasets and endpoints differ.
"""

from __future__ import annotations

import html
import json
import shutil
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
GA = ROOT / "project/results/cross_neo_ga_rl_full_comparison_2026_05_10"
MOD = ROOT / "project/results/cross_neo_moderna_bd_package_2026_05_11"
OUT = ROOT / "project/results/cross_neo_latest_sota_watchlist_2026_05_11"
FIG = OUT / "figures"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET = HUB / "assets/cross_neo_latest_sota_watchlist"
LIVE_ASSET = LIVE / "assets/cross_neo_latest_sota_watchlist"
PAGE = HUB / "cross_neo_latest_sota_watchlist.html"
LIVE_PAGE = LIVE / "cross_neo_latest_sota_watchlist.html"


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


def fmt(value: object, digits: int = 3) -> str:
    try:
        if pd.isna(value):
            return "NA"
        return f"{float(value):.{digits}f}"
    except Exception:
        return esc(value)


def safe_copy(src: Path, dst: Path) -> None:
    try:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    except PermissionError:
        print(f"[sota-watchlist] skip permission-denied copy {dst}")


def table_html(df: pd.DataFrame, cols: list[str], max_rows: int = 50) -> str:
    if df.empty:
        return "<p class='muted'>No rows.</p>"
    keep = [c for c in cols if c in df.columns]
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


def build_local_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    master = read_tsv(GA / "ga_rl_master_algorithm_comparison.tsv")
    public = read_tsv(GA / "ga_rl_public_competitor_context_summary.tsv")
    strict = master[master["algorithm"].isin([
        "BAR-Neo_confidence", "BAR-Neo", "BAR-Neo_patient_gated", "BAR-Neo-X",
        "BAR-Neo-X_claim_safe", "RF_biophys", "MHCflurry_2.0_presentation", "KG_GA_evolved_controller",
    ])].copy()
    return master, public, strict


def build_literature_watchlist() -> pd.DataFrame:
    rows = [
        {
            "method": "NetMHCpan-4.2",
            "task_scope": "MHC-I antigen presentation / neoepitope mode",
            "latest_signal": "Official DTU server updated 2026-03-24; transfer learning + structural features; neoepitope output mode available.",
            "reported_metric": "No single headline AUPRC in abstract; performance improved modestly vs 4.1-family baselines.",
            "source": "https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/",
        },
        {
            "method": "ImmugenX",
            "task_scope": "pMHC immunogenicity",
            "latest_signal": "Modular protein language model with optional TCR context.",
            "reported_metric": "AUROC 0.619, AP 0.514; +7% AP vs next best model on pMHC immunogenicity.",
            "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11581412/",
        },
        {
            "method": "PepFore / NeoaPred",
            "task_scope": "neoantigen foreignness / immunogenicity",
            "latest_signal": "Structure- and surface-aware foreignness score; compared against MixMHCpred, NetMHCpan, MHCflurry, PRIME, BigMHC.",
            "reported_metric": "Best comparison in paper: AUROC 0.81, AUPRC 0.54; BigMHC 0.70/0.30 on same test set.",
            "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11419954/",
        },
        {
            "method": "NeoGuider",
            "task_scope": "neoepitope prioritization / immunogenicity",
            "latest_signal": "Advanced feature engineering with KDE + centered isotonic regression; benchmarked on 7 cohorts, 113 patients, 635 immunogenic candidates.",
            "reported_metric": "Outperformed existing methods; abstract does not expose a single numeric headline metric.",
            "source": "https://pubmed.ncbi.nlm.nih.gov/41437280/",
        },
        {
            "method": "TrambaHLApan",
            "task_scope": "presentation + immunogenicity joint predictor",
            "latest_signal": "Transformer + Mamba hybrid, explicit EL and IM heads.",
            "reported_metric": "Abstract says it outperforms SOTA on independent datasets; no headline metric exposed in abstract.",
            "source": "https://pubmed.ncbi.nlm.nih.gov/41087632/",
        },
        {
            "method": "NeoTImmuML",
            "task_scope": "tumor neoantigen immunogenicity",
            "latest_signal": "Weighted ensemble over LightGBM/XGBoost/RF; independent test set assembled from 2024-2025 literature.",
            "reported_metric": "Best model AUROC 0.8707 on independent test set; weighted ensemble slightly above voting ensemble.",
            "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12585993/",
        },
        {
            "method": "NetMHCpan-4.1 / BigMHC / MHCflurry / PRIME / TransPHLA",
            "task_scope": "baseline comparators",
            "latest_signal": "Common baselines in current papers and local repo context.",
            "reported_metric": "Used as baselines; local repo context shows RF_biophys/PRIME/MHCflurry/TransPHLA/BigMHC ranges, but cross-paper metrics are not directly comparable.",
            "source": "local-repo context + cited papers above",
        },
    ]
    return pd.DataFrame(rows)


def build_html(local: pd.DataFrame, public: pd.DataFrame, strict: pd.DataFrame, lit: pd.DataFrame) -> str:
    kg = local[local["algorithm"].eq("KG_GA_evolved")].iloc[0]
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Cross-Neo Latest SOTA Watchlist</title>
<style>
:root {{ --bg:#07111f; --panel:#101d2d; --ink:#edf5ff; --muted:#9fb0c5; --gold:#f2c46d; --line:#27384f; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font-family:Inter,Arial,sans-serif; line-height:1.5; }}
header {{ padding:48px 5vw 28px; border-bottom:1px solid var(--line); background:#0b1626; }}
h1 {{ font-size:clamp(30px,5vw,60px); margin:8px 0 10px; line-height:1.0; }}
.lead {{ color:var(--muted); max-width:1120px; font-size:18px; }}
.stats {{ display:grid; grid-template-columns:repeat(4,minmax(140px,1fr)); gap:12px; margin-top:24px; }}
.stat {{ background:var(--panel); border:1px solid var(--line); padding:14px; border-radius:8px; }}
.stat b {{ display:block; font-size:26px; color:var(--gold); }}
.stat span {{ color:var(--muted); font-size:12px; }}
main {{ padding:28px 5vw 60px; max-width:1500px; margin:0 auto; }}
section {{ margin:34px 0; }}
h2 {{ font-size:24px; margin-bottom:12px; }}
.note {{ border-left:4px solid var(--gold); background:#101d2d; padding:14px 16px; color:var(--muted); }}
table {{ width:100%; border-collapse:collapse; font-size:13px; background:var(--panel); }}
th, td {{ border:1px solid var(--line); padding:8px; vertical-align:top; }}
th {{ color:var(--gold); text-align:left; position:sticky; top:0; background:#101d2d; }}
.tablewrap {{ overflow:auto; border:1px solid var(--line); }}
.grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(320px,1fr)); gap:16px; }}
.card {{ border:1px solid var(--line); background:var(--panel); border-radius:8px; padding:12px; }}
.muted {{ color:var(--muted); }}
a {{ color:var(--gold); }}
@media (max-width:900px) {{ .stats {{ grid-template-columns:repeat(2,1fr); }} }}
</style></head><body>
<header>
<h1>Latest SOTA watchlist, with local comparison kept separate.</h1>
<p class="lead">This page does not pretend all neoantigen models solve the same task. It separates same-board local ranking, strict no-overlap reliability, and literature-only latest SOTA methods whose datasets and metrics are not directly comparable.</p>
<div class="stats">
<div class="stat"><b>{kg['all_AUPRC']:.3f}</b><span>KG-GA same-board AUPRC</span></div>
<div class="stat"><b>{kg['frozen_validation_AUPRC']:.3f}</b><span>frozen validation-like AUPRC</span></div>
<div class="stat"><b>{float(strict[strict['algorithm'].eq('BAR-Neo_confidence')]['neo_strict_AUPRC'].iloc[0]):.3f}</b><span>strict no-overlap winner</span></div>
<div class="stat"><b>{len(lit)}</b><span>latest literature method groups</span></div>
</div>
</header><main>
<section><h2>Rule</h2><p class="note"><b>Same-board AUPRC beats paper-by-paper bragging.</b> Literature SOTA only means the paper says it outperformed its baselines on its own cohorts. That is useful, but not a substitute for a locked same-board benchmark.</p></section>
<section><h2>Local Same-Board Comparison</h2><div class="tablewrap">{table_html(local[local['algorithm'].isin(['KG_GA_evolved','CROSS_integrated','CROSS_claimsafe','CROSS_stress','CROSS_BMA','CROSS_finetuned','BigMHC_IM','BigMHC_EL'])], ['algorithm','all_AUPRC','all_AUROC','frozen_validation_AUPRC','low_medium_leakage_AUPRC','v6_smoke_AUPRC','paperability_reward'], 20)}</div></section>
<section><h2>Strict No-Overlap Reliability</h2><div class="tablewrap">{table_html(strict, ['algorithm','neo_strict_AUPRC','neo_strict_AUROC'], 20)}</div></section>
<section><h2>Broader Public Context</h2><div class="tablewrap">{table_html(public, ['method','n_context_rows','mean_context_AUROC','best_context_AUROC','claim_use'], 20)}</div></section>
<section><h2>Latest Literature SOTA Watchlist</h2><div class="tablewrap">{table_html(lit, ['method','task_scope','latest_signal','reported_metric','source'], 20)}</div></section>
<section><h2>What Is Actually Strong Right Now</h2>
<div class="card">
<p><b>Selection / prioritization:</b> KG_GA_evolved is the local champion on the 2,715-row board.</p>
<p><b>Reliability / strict no-overlap:</b> BAR-Neo_confidence is still stronger on the strict lane.</p>
<p><b>Presentation modeling:</b> NetMHCpan-4.2 is the newest official MHC-I presentation server update in the literature set reviewed here.</p>
<p><b>Direct immunogenicity:</b> ImmugenX, PepFore/NeoaPred, NeoGuider, TrambaHLApan, and NeoTImmuML are the newest literature models in this scan, but their metrics live on different datasets and should not be merged into one leaderboard without reruns.</p>
</div></section>
</main></body></html>"""


def build_figures(local: pd.DataFrame, strict: pd.DataFrame, lit: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    plt.style.use("default")
    same = local[local["algorithm"].isin(["KG_GA_evolved", "CROSS_integrated", "CROSS_claimsafe", "CROSS_stress", "CROSS_BMA", "CROSS_finetuned", "BigMHC_IM", "BigMHC_EL"])].copy()
    fig, ax = plt.subplots(figsize=(10.5, 5))
    same = same.sort_values("all_AUPRC")
    ax.barh(same["algorithm"], same["all_AUPRC"], color=["#c7922b" if a == "KG_GA_evolved" else "#315f7d" for a in same["algorithm"]])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Same-board AUPRC")
    ax.set_title("Local same-board comparison")
    fig.tight_layout()
    fig.savefig(FIG / "fig1_same_board.png", dpi=220)
    fig.savefig(FIG / "fig1_same_board.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    strict = strict.sort_values("neo_strict_AUPRC")
    ax.barh(strict["algorithm"], strict["neo_strict_AUPRC"], color=["#b24a62" if "KG_GA" in a else "#315f7d" for a in strict["algorithm"]])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Strict no-overlap AUPRC")
    ax.set_title("Strict no-overlap reliability lane")
    fig.tight_layout()
    fig.savefig(FIG / "fig2_strict.png", dpi=220)
    fig.savefig(FIG / "fig2_strict.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    y = np.arange(len(lit))
    colors = ["#c7922b" if "NetMHCpan-4.2" in m else "#315f7d" for m in lit["method"]]
    ax.barh(y, np.arange(len(lit)), color=colors)
    ax.set_yticks(y, lit["method"])
    ax.set_xlabel("Method index only")
    ax.set_title("Latest literature SOTA watchlist")
    fig.tight_layout()
    fig.savefig(FIG / "fig3_lit_watchlist.png", dpi=220)
    fig.savefig(FIG / "fig3_lit_watchlist.pdf")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    ASSET.mkdir(parents=True, exist_ok=True)
    LIVE_ASSET.mkdir(parents=True, exist_ok=True)

    local, public, strict = build_local_tables()
    lit = build_literature_watchlist()

    local.to_csv(OUT / "local_same_board_comparison.tsv", sep="\t", index=False)
    public.to_csv(OUT / "local_public_context_summary.tsv", sep="\t", index=False)
    strict.to_csv(OUT / "local_strict_reliability.tsv", sep="\t", index=False)
    lit.to_csv(OUT / "latest_literature_sota_watchlist.tsv", sep="\t", index=False)

    build_figures(local, strict, lit)
    PAGE.write_text(build_html(local, public, strict, lit), encoding="utf-8")

    for src in FIG.glob("*"):
        if src.is_file():
            safe_copy(src, ASSET / src.name)
            safe_copy(src, LIVE_ASSET / src.name)
    for src in [
        OUT / "local_same_board_comparison.tsv",
        OUT / "local_public_context_summary.tsv",
        OUT / "local_strict_reliability.tsv",
        OUT / "latest_literature_sota_watchlist.tsv",
    ]:
        safe_copy(src, ASSET / src.name)
        safe_copy(src, LIVE_ASSET / src.name)
    safe_copy(PAGE, LIVE_PAGE)

    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "out_dir": str(OUT),
        "html_path": str(PAGE),
        "live_html_path": str(LIVE_PAGE),
        "local_champion": "KG_GA_evolved",
        "strict_winner": "BAR-Neo_confidence",
        "literature_methods": int(len(lit)),
    }
    (OUT / "latest_sota_watchlist_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
