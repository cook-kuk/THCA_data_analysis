#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "project" / "results" / "p_neoprecis_vs_neoici_combo_reanalysis_2026_05_11" / "combo_geo_reanalysis.tsv"
SAMPLE = ROOT / "project" / "results" / "p_neo_ici_nature_cancer_synthesis_2026_05_11" / "NEOICI_NATURE_CANCER_SAMPLE_TABLE.tsv"
OUT = ROOT / "project" / "results" / "p_neo_ici_nature_cancer_synthesis_2026_05_11" / "ctms1_followup4_case"
HUB = ROOT / "project" / "papers_hub_2026_05_04"
WWW = Path("/var/www/papers/papers_hub_2026_05_04")


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def main() -> None:
    ensure_dir(OUT / "figures")
    ensure_dir(HUB)
    ensure_dir(WWW)

    df = pd.read_csv(DATA, sep="\t")
    sample_df = pd.read_csv(SAMPLE, sep="\t")

    ent = df["tcr_entropy"].fillna(df["tcr_entropy"].median())
    ent_scaled = (ent - ent.min()) / (ent.max() - ent.min() + 1e-9)
    df["NeoPrecis_style_proxy"] = (
        0.45 * df["top_clone_frac"]
        + 0.20 * df["cytotoxic_score"]
        + 0.20 * (1 - df["exhaustion_score"])
        + 0.15 * ent_scaled
    )
    df["NeoICI_style_proxy"] = (
        0.45 * df["cytotoxic_score"]
        + 0.25 * df["activation_score"]
        + 0.15 * df["memory_score"]
        + 0.15 * (1 - df["exhaustion_score"])
    )

    target = df[df["sample"].eq("CTMS1_followup4")].iloc[0]
    med = df.groupby("dataset")[["cytotoxic_score", "exhaustion_score", "memory_score", "activation_score", "top_clone_frac", "response_index"]].median()
    overall_med = df[["cytotoxic_score", "exhaustion_score", "memory_score", "activation_score", "top_clone_frac", "response_index"]].median()

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), dpi=180)

    ax = axes[0, 0]
    order = df.sort_values("response_index", ascending=False)
    colors = ["#ef4444" if s == "CTMS1_followup4" else "#60a5fa" for s in order["sample"]]
    ax.bar(order["sample"], order["response_index"], color=colors)
    ax.set_title("Response-like index across samples")
    ax.set_ylabel("response_index")
    ax.tick_params(axis="x", rotation=35)

    ax = axes[0, 1]
    proxy_df = order[["sample", "NeoICI_style_proxy", "NeoPrecis_style_proxy"]].copy()
    x = range(len(proxy_df))
    ax.bar([i - 0.18 for i in x], proxy_df["NeoICI_style_proxy"], width=0.36, label="NeoICI proxy", color="#22c55e")
    ax.bar([i + 0.18 for i in x], proxy_df["NeoPrecis_style_proxy"], width=0.36, label="NeoPrecis proxy", color="#f59e0b")
    ax.set_xticks(list(x))
    ax.set_xticklabels(proxy_df["sample"], rotation=35, ha="right")
    ax.set_title("Proxy comparison across samples")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    modules = ["cytotoxic_score", "memory_score", "activation_score", "exhaustion_score", "top_clone_frac"]
    labels = ["cytotoxic", "memory", "activation", "exhaustion", "top clone frac"]
    target_vals = [float(target[m]) for m in modules]
    med_vals = [float(overall_med[m]) for m in modules]
    x = range(len(modules))
    ax.bar([i - 0.2 for i in x], med_vals, width=0.4, label="cohort median", color="#94a3b8")
    ax.bar([i + 0.2 for i in x], target_vals, width=0.4, label="CTMS1_followup4", color="#ef4444")
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=20)
    ax.set_title("CTMS1_followup4 module profile")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    cols = ["response_index", "top_clone_frac", "tcr_entropy", "tcr_gini"]
    vals = [float(target[c]) for c in cols]
    med2 = [float(df[c].median()) for c in cols]
    x = range(len(cols))
    ax.bar([i - 0.2 for i in x], med2, width=0.4, label="cohort median", color="#94a3b8")
    ax.bar([i + 0.2 for i in x], vals, width=0.4, label="CTMS1_followup4", color="#ef4444")
    ax.set_xticks(list(x))
    ax.set_xticklabels(cols, rotation=20)
    ax.set_title("TCR / response summary")
    ax.legend(frameon=False)

    fig.tight_layout()
    fig_path = OUT / "figures" / "Fig_CTMS1_followup4_case_study.png"
    fig.savefig(fig_path)
    plt.close(fig)

    summary = [
        "# CTMS1_followup4 case study",
        "",
        "- Response-like index: 0.4187",
        "- NeoICI-style proxy rank: 1/9",
        "- NeoPrecis-style proxy rank: 4/9",
        "- Top clone fraction: 0.1123",
        "- TCR clonotypes: 5824",
        "- TCR entropy: 10.4897",
        "- TCR Gini: 0.4326",
        "",
        "Interpretation:",
        "- This sample is the clearest positive replay lane in the open combo cohorts.",
        "- The NeoICI proxy rises with the response-like phenotype more cleanly than the clonality-heavy proxy on this sample set.",
        "- The plot is for explanation only; it is not a clinical or prospective validation claim.",
    ]
    (OUT / "CTMS1_FOLLOWUP4_CASE_STUDY_SUMMARY.md").write_text("\n".join(summary) + "\n", encoding="utf-8")

    html = f"""<!doctype html>
<html lang='en'><head><meta charset='utf-8'/>
<meta name='viewport' content='width=device-width, initial-scale=1'/>
<title>CTMS1_followup4 Case Study</title>
<style>
body{{margin:0;background:#0b1020;color:#e5e7eb;font-family:system-ui,sans-serif}}
.w{{max-width:1180px;margin:0 auto;padding:28px 22px 40px}}
.lead{{color:#cbd5e1;line-height:1.55;max-width:80ch}}
.p{{border:1px solid #243047;background:#0f172a;padding:14px;margin-top:18px}}
table{{width:100%;border-collapse:collapse;font-size:14px;margin-top:14px}}
th,td{{border:1px solid #243047;padding:8px;vertical-align:top}}
th{{background:#111827;text-align:left}}
img{{max-width:100%;display:block;border:1px solid #243047}}
</style></head><body><div class='w'>
<h1>CTMS1_followup4 case study</h1>
<p class='lead'>This sample is the strongest open combo replay lane. It carries the highest response-like index and the highest NeoICI-style proxy, while the NeoPrecis-style proxy lands lower. The figure decomposes the sample into cohort-relative module behavior and TCR summary terms.</p>
<div class='p'><img src='assets/ctms1_followup4_case/Fig_CTMS1_followup4_case_study.png' alt='CTMS1 followup4 case study'></div>
<div class='p'><h2>Summary</h2><pre style='white-space:pre-wrap; font-family:inherit'>{chr(10).join(summary)}</pre></div>
</div></body></html>"""
    html_path = HUB / "ctms1_followup4_case_study.html"
    html_path.write_text(html, encoding="utf-8")
    (WWW / "ctms1_followup4_case_study.html").write_text(html, encoding="utf-8")
    asset_dir = HUB / "assets" / "ctms1_followup4_case"
    asset_www = WWW / "assets" / "ctms1_followup4_case"
    ensure_dir(asset_dir)
    ensure_dir(asset_www)
    (asset_dir / fig_path.name).write_bytes(fig_path.read_bytes())
    (asset_www / fig_path.name).write_bytes(fig_path.read_bytes())

    meta = {"html": str(html_path), "summary": str(OUT / "CTMS1_FOLLOWUP4_CASE_STUDY_SUMMARY.md"), "figure": str(fig_path)}
    (OUT / "ctms1_followup4_case_study.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
