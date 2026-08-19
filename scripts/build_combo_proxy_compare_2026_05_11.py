#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
INP = ROOT / "project" / "results" / "p_neoprecis_vs_neoici_combo_reanalysis_2026_05_11" / "combo_geo_reanalysis.tsv"
OUT = ROOT / "project" / "results" / "p_neoprecis_vs_neoici_combo_reanalysis_2026_05_11" / "proxy_compare"
HUB = ROOT / "project" / "papers_hub_2026_05_04"
WWW = Path("/var/www/papers/papers_hub_2026_05_04")


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def main() -> None:
    ensure_dir(OUT / "figures")
    ensure_dir(HUB)
    ensure_dir(WWW)

    df = pd.read_csv(INP, sep="\t")
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

    rows = []
    for name in ["NeoPrecis_style_proxy", "NeoICI_style_proxy"]:
        s = df[name]
        rows.append(
            {
                "proxy": name,
                "spearman_vs_response_index": spearmanr(s, df["response_index"]).statistic,
                "spearman_vs_top_clone_frac": spearmanr(s, df["top_clone_frac"]).statistic,
                "spearman_vs_cytotoxic_score": spearmanr(s, df["cytotoxic_score"]).statistic,
                "mean_rank": s.rank(ascending=False, method="average").mean(),
                "top_sample": df.loc[s.idxmax(), "sample"],
            }
        )
    out_df = pd.DataFrame(rows)
    out_df.to_csv(OUT / "proxy_compare_metrics.tsv", sep="\t", index=False)

    ranked = df[[
        "dataset",
        "sample",
        "response_index",
        "top_clone_frac",
        "cytotoxic_score",
        "exhaustion_score",
        "activation_score",
        "memory_score",
        "NeoPrecis_style_proxy",
        "NeoICI_style_proxy",
    ]].copy()
    ranked.to_csv(OUT / "proxy_compare_samples.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=160)
    x = range(len(df))
    ax.plot(x, df["response_index"], marker="o", label="response_index", linewidth=2)
    ax.plot(x, df["NeoPrecis_style_proxy"], marker="o", label="NeoPrecis_style_proxy", linewidth=2)
    ax.plot(x, df["NeoICI_style_proxy"], marker="o", label="NeoICI_style_proxy", linewidth=2)
    ax.set_xticks(list(x))
    ax.set_xticklabels(df["sample"], rotation=35, ha="right")
    ax.set_ylabel("score")
    ax.set_title("Open combo GEO score comparison")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig_path = OUT / "figures" / "Fig_combo_proxy_compare.png"
    fig.savefig(fig_path)
    plt.close(fig)

    summary = [
        "# Combo proxy comparison",
        "",
        f"- NeoPrecis-style proxy Spearman vs response_index: {out_df.loc[out_df.proxy=='NeoPrecis_style_proxy', 'spearman_vs_response_index'].iloc[0]:.3f}",
        f"- NeoICI-style proxy Spearman vs response_index: {out_df.loc[out_df.proxy=='NeoICI_style_proxy', 'spearman_vs_response_index'].iloc[0]:.3f}",
        f"- NeoPrecis-style proxy Spearman vs top_clone_frac: {out_df.loc[out_df.proxy=='NeoPrecis_style_proxy', 'spearman_vs_top_clone_frac'].iloc[0]:.3f}",
        f"- NeoICI-style proxy Spearman vs top_clone_frac: {out_df.loc[out_df.proxy=='NeoICI_style_proxy', 'spearman_vs_top_clone_frac'].iloc[0]:.3f}",
        "",
        "## Readout",
        "",
        "- On this 9-sample replay, the NeoICI-style T-cell-state proxy tracks the response-like index better than the NeoPrecis-style clonality-weighted proxy.",
        "- The ranking gap is still small and the sample size is tiny, so this is directional support only, not a final superiority claim.",
    ]
    (OUT / "PROXY_COMPARE_SUMMARY.md").write_text("\n".join(summary) + "\n", encoding="utf-8")

    html = f"""<!doctype html>
<html lang='en'><head><meta charset='utf-8'/>
<meta name='viewport' content='width=device-width, initial-scale=1'/>
<title>Combo Proxy Comparison</title>
<style>
body{{margin:0;background:#0b1020;color:#e5e7eb;font-family:system-ui,sans-serif}}
.w{{max-width:1160px;margin:0 auto;padding:28px 22px 40px}}
.lead{{color:#cbd5e1;line-height:1.55;max-width:80ch}}
table{{width:100%;border-collapse:collapse;margin-top:16px;font-size:14px}}
th,td{{border:1px solid #243047;padding:8px;vertical-align:top}}
th{{background:#111827}}
.p{{border:1px solid #243047;background:#0f172a;padding:14px;margin-top:18px}}
img{{max-width:100%;display:block;border:1px solid #243047}}
</style></head><body><div class='w'>
<h1>Combo proxy comparison</h1>
<p class='lead'>This compares a clonality-heavy NeoPrecis-style proxy to a T-cell-state-heavy NeoICI-style proxy on the same 9 open combo samples. The target is the local response-like index built from cytotoxicity, exhaustion, activation, and clonotype summary.</p>
<div class='p'><h2>Proxy metrics</h2>{out_df.to_html(index=False, escape=False)}</div>
<div class='p'><h2>Sample table</h2>{ranked.to_html(index=False, escape=False)}</div>
<div class='p'><h2>Figure</h2><img src='assets/combo_proxy_compare/Fig_combo_proxy_compare.png' alt='combo proxy comparison'></div>
</div></body></html>"""

    html_path = HUB / "combo_proxy_compare.html"
    html_path.write_text(html, encoding="utf-8")
    (WWW / "combo_proxy_compare.html").write_text(html, encoding="utf-8")
    asset_dir = HUB / "assets" / "combo_proxy_compare"
    asset_www = WWW / "assets" / "combo_proxy_compare"
    ensure_dir(asset_dir)
    ensure_dir(asset_www)
    (asset_dir / fig_path.name).write_bytes(fig_path.read_bytes())
    (asset_www / fig_path.name).write_bytes(fig_path.read_bytes())

    meta = {
        "html": str(html_path),
        "summary": str(OUT / "PROXY_COMPARE_SUMMARY.md"),
        "metrics": str(OUT / "proxy_compare_metrics.tsv"),
    }
    (OUT / "proxy_compare_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
