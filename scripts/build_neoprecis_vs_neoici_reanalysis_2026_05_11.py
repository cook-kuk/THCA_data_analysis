#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "project" / "results" / "p_neoprecis_vs_neoici_2026_05_11"
HUB = ROOT / "project" / "papers_hub_2026_05_04"
WWW = Path("/var/www/papers/papers_hub_2026_05_04")


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def load_neoici_metrics() -> pd.DataFrame:
    path = ROOT / "project" / "results" / "p_neo_ici_combo_2026_05_11" / "neoici_algorithm_gauntlet_metrics.tsv"
    df = pd.read_csv(path, sep="\t")
    keep = df[df["split"].isin(["academic_validation_like", "industrial_locked_v0", "academic_all"])].copy()
    return keep


def build_rows() -> pd.DataFrame:
    rows = [
        {
            "system": "NeoPrecis",
            "layer": "immunogenicity",
            "benchmark": "CEDAR / NCI / T-cell and landscape evaluation",
            "published_contract": "NeoPrecis-Immuno and NeoPrecis landscape are benchmarked on CEDAR, NCI GI, melanoma, and NSCLC cohorts; MHC-II AUROC is stated as best, and NP-Integrated / clonality-aware landscape improves ICI-response modeling.",
            "comparable_now": "partial",
            "current_status": "code available in GitHub repo; no same-benchmark local re-run completed here",
            "source": "Nature Communications article + GitHub repo",
        },
        {
            "system": "NeoICI_GA_RL_combo_v1",
            "layer": "combo prioritization",
            "benchmark": "our validation-like local gauntlet",
            "published_contract": "experiment-priority score combining HLA-I, helper/dual-class, TCR/MD, BioDarwin mode-bank, and GA/RL",
            "comparable_now": "yes",
            "current_status": "AUPRC 0.863 / AUROC 0.983 on academic_validation_like; 0.563 / 0.597 on industrial_locked_v0",
            "source": "local TSV gauntlet",
        },
        {
            "system": "KG_GA_evolved",
            "layer": "prior GA/RL champion",
            "benchmark": "our validation-like local gauntlet",
            "published_contract": "retrospective + validation-like champion, frozen before prospective wetlab",
            "comparable_now": "yes",
            "current_status": "AUPRC 0.978 / AUROC 0.999 on academic_validation_like; 0.933 / 0.943 on academic_all",
            "source": "local TSV gauntlet",
        },
        {
            "system": "BioDarwin_mode_bank_v4",
            "layer": "industrial anchor",
            "benchmark": "our locked public mini-set",
            "published_contract": "mode-bank / public anchor comparator",
            "comparable_now": "yes",
            "current_status": "AUPRC 0.628 / AUROC 0.639 on industrial_locked_v0",
            "source": "local TSV gauntlet",
        },
    ]
    return pd.DataFrame(rows)


def make_plot(df: pd.DataFrame) -> Path:
    score_rows = []
    for _, row in df.iterrows():
        if row["system"] in {"NeoICI_GA_RL_combo_v1", "KG_GA_evolved", "BioDarwin_mode_bank_v4"}:
            metrics = row["current_status"]
            vals = {}
            for part in metrics.split(";"):
                part = part.strip()
                if "AUPRC" in part and "/" in part:
                    vals["AUPRC"] = float(part.split("AUPRC")[1].split("/")[0].strip())
                    vals["AUROC"] = float(part.split("AUROC")[1].split("on")[0].strip())
            if vals:
                score_rows.append({"system": row["system"], **vals})
    score_df = pd.DataFrame(score_rows)
    fig, ax = plt.subplots(figsize=(6, 3.5), dpi=160)
    x = range(len(score_df))
    ax.bar([i - 0.15 for i in x], score_df["AUPRC"], width=0.3, label="AUPRC")
    ax.bar([i + 0.15 for i in x], score_df["AUROC"], width=0.3, label="AUROC")
    ax.set_xticks(list(x))
    ax.set_xticklabels(score_df["system"], rotation=25, ha="right")
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=False)
    ax.set_title("Local performance only")
    fig.tight_layout()
    path = OUT / "figures" / "Fig_NeoPrecis_vs_NeoICI_local_metrics.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def main() -> None:
    ensure_dir(OUT / "figures")
    ensure_dir(HUB)
    ensure_dir(WWW)

    neoici = load_neoici_metrics()
    rows = build_rows()
    rows.to_csv(OUT / "neo_precis_vs_neoici_comparison.tsv", sep="\t", index=False)

    summary = [
        "# NeoPrecis vs NeoICI reanalysis",
        "",
        "## What is directly comparable",
        "",
        "- NeoPrecis is a published immunogenicity + landscape framework with CEDAR/NCI/melanoma/NSCLC benchmarks.",
        "- NeoICI is our local combo-therapy prioritization layer tested on the current gauntlet tables.",
        "",
        "## What is not yet directly comparable",
        "",
        "- I do not have a same-dataset score table for NeoPrecis and NeoICI in this workspace.",
        "- So this is a performance-contract comparison, not a strict paired re-run.",
        "",
        "## Current local result",
        "",
        f"- `KG_GA_evolved` remains the strongest local scorer on the validation-like split: AUPRC {neoici.loc[neoici['algorithm'].eq('KG_GA_evolved') & neoici['split'].eq('academic_validation_like'), 'AUPRC'].iloc[0]:.3f}.",
        f"- `NeoICI_GA_RL_combo_v1` is strong but not the top local scorer: AUPRC {neoici.loc[neoici['algorithm'].eq('NeoICI_GA_RL_combo_v1') & neoici['split'].eq('academic_validation_like'), 'AUPRC'].iloc[0]:.3f}.",
        f"- `BioDarwin_mode_bank_v4` still anchors the locked industrial mini-set: AUPRC {neoici.loc[neoici['algorithm'].eq('BioDarwin_mode_bank_v4') & neoici['split'].eq('industrial_locked_v0'), 'AUPRC'].iloc[0]:.3f}.",
        "",
        "## Decision",
        "",
        "- No, we cannot say the new package is universally better than NeoPrecis yet.",
        "- Yes, we can say the package is broader in scope: it adds combo-therapy readiness, GA/RL search, and explicit claim boundaries.",
    ]
    (OUT / "NEOPRECIS_VS_NEOICI_SUMMARY.md").write_text("\n".join(summary) + "\n", encoding="utf-8")

    fig_path = make_plot(rows)

    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>NeoPrecis vs NeoICI Reanalysis</title>
  <style>
    body {{ margin:0; background:#0b1020; color:#e5e7eb; font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; }}
    .wrap {{ max-width: 1180px; margin: 0 auto; padding: 28px 22px 48px; }}
    h1, h2 {{ margin:0 0 12px; }}
    .lead {{ color:#cbd5e1; line-height:1.55; max-width: 80ch; }}
    table {{ width:100%; border-collapse:collapse; margin-top:16px; font-size:14px; }}
    th, td {{ border:1px solid #243047; padding:10px; vertical-align:top; }}
    th {{ background:#111827; text-align:left; }}
    tr:nth-child(even) td {{ background:#0f172a; }}
    .grid {{ display:grid; grid-template-columns: 1.1fr .9fr; gap:18px; margin-top:18px; }}
    .panel {{ border:1px solid #243047; background:#0f172a; padding:16px; }}
    img {{ max-width:100%; display:block; border:1px solid #243047; background:#0b1020; }}
    code {{ background:#111827; padding:1px 4px; border-radius:4px; }}
  </style>
</head>
<body>
<div class="wrap">
  <h1>NeoPrecis vs NeoICI</h1>
  <p class="lead">This is a performance-contract reanalysis. NeoPrecis has a published immunogenicity and clonality-aware landscape framework; NeoICI is our combo-therapy prioritization layer. The comparison is honest but not yet a strict paired rerun because the same benchmark table is not available locally for both systems.</p>
  <div class="grid">
    <div class="panel">
      <h2>Comparison table</h2>
      {rows.to_html(index=False, escape=False)}
    </div>
    <div class="panel">
      <h2>Local metrics</h2>
      <img src="assets/neo_precis_vs_neoici/Fig_NeoPrecis_vs_NeoICI_local_metrics.png" alt="Local performance metrics">
      <p class="lead">The strongest local validation-like scorer remains <code>KG_GA_evolved</code>; <code>NeoICI_GA_RL_combo_v1</code> is strong but not top. Industrial locked still favors <code>BioDarwin_mode_bank_v4</code>.</p>
    </div>
  </div>
</div>
</body>
</html>
"""
    html_path = HUB / "neoprecis_vs_neoici_reanalysis.html"
    html_path.write_text(html, encoding="utf-8")
    (WWW / "neoprecis_vs_neoici_reanalysis.html").write_text(html, encoding="utf-8")
    asset_dir = HUB / "assets" / "neo_precis_vs_neoici"
    asset_www = WWW / "assets" / "neo_precis_vs_neoici"
    ensure_dir(asset_dir)
    ensure_dir(asset_www)
    fig_dst = asset_dir / fig_path.name
    fig_www = asset_www / fig_path.name
    fig_dst.write_bytes(fig_path.read_bytes())
    fig_www.write_bytes(fig_path.read_bytes())

    meta = {
        "generated": "2026-05-11",
        "html": str(html_path),
        "summary": str(OUT / "NEOPRECIS_VS_NEOICI_SUMMARY.md"),
        "figure": str(fig_path),
    }
    (OUT / "neo_precis_vs_neoici_reanalysis.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
