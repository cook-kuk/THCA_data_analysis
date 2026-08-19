#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULT_A = ROOT / "project" / "results" / "p_neo_ici_combo_2026_05_11"
RESULT_B = ROOT / "project" / "results" / "p_neoprecis_vs_neoici_combo_reanalysis_2026_05_11"
HUB = ROOT / "project" / "papers_hub_2026_05_04"
WWW = Path("/var/www/papers/papers_hub_2026_05_04")
OUT = ROOT / "project" / "results" / "p_neo_ici_nature_cancer_synthesis_2026_05_11"


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def load_metrics() -> dict[str, float]:
    g = pd.read_csv(RESULT_A / "neoici_algorithm_gauntlet_metrics.tsv", sep="\t")
    c = pd.read_csv(RESULT_B / "combo_geo_reanalysis.tsv", sep="\t")
    p = pd.read_csv(RESULT_B / "proxy_compare" / "proxy_compare_metrics.tsv", sep="\t")
    s = c.copy()
    import numpy as np
    ent = s["tcr_entropy"].fillna(s["tcr_entropy"].median())
    ent_scaled = (ent - ent.min()) / (ent.max() - ent.min() + 1e-9)
    s["NeoPrecis_style_proxy"] = (
        0.45 * s["top_clone_frac"]
        + 0.20 * s["cytotoxic_score"]
        + 0.20 * (1 - s["exhaustion_score"])
        + 0.15 * ent_scaled
    )
    s["NeoICI_style_proxy"] = (
        0.45 * s["cytotoxic_score"]
        + 0.25 * s["activation_score"]
        + 0.15 * s["memory_score"]
        + 0.15 * (1 - s["exhaustion_score"])
    )
    return {
        "kg_ga_validation_auprc": float(g.query("split == 'academic_validation_like' and algorithm == 'KG_GA_evolved'")["AUPRC"].iloc[0]),
        "neoici_validation_auprc": float(g.query("split == 'academic_validation_like' and algorithm == 'NeoICI_GA_RL_combo_v1'")["AUPRC"].iloc[0]),
        "industrial_biodarwin_auprc": float(g.query("split == 'industrial_locked_v0' and algorithm == 'BioDarwin_mode_bank_v4'")["AUPRC"].iloc[0]),
        "combo_n": int(c.shape[0]),
        "combo_response_median": float(c["response_index"].median()),
        "combo_clone_median": float(c["top_clone_frac"].median()),
        "neoici_proxy_rho": float(p.query("proxy == 'NeoICI_style_proxy'")["spearman_vs_response_index"].iloc[0]),
        "neoprecis_proxy_rho": float(p.query("proxy == 'NeoPrecis_style_proxy'")["spearman_vs_response_index"].iloc[0]),
        "neoici_proxy_rho_clone": float(p.query("proxy == 'NeoICI_style_proxy'")["spearman_vs_top_clone_frac"].iloc[0]),
        "neoprecis_proxy_rho_clone": float(p.query("proxy == 'NeoPrecis_style_proxy'")["spearman_vs_top_clone_frac"].iloc[0]),
    }


def make_overview_figure() -> Path:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), dpi=170)
    panels = [
        (
            RESULT_A / "figures" / "Fig_NeoICI_algorithm_gauntlet_AUPRC.png",
            "Algorithm gauntlet",
        ),
        (
            RESULT_B / "figures" / "Fig_combo_response_index.png",
            "Open combo GEO replay",
        ),
        (
            RESULT_B / "proxy_compare" / "figures" / "Fig_combo_proxy_compare.png",
            "Proxy comparison",
        ),
    ]
    for ax, (path, title) in zip(axes, panels):
        img = plt.imread(path)
        ax.imshow(img)
        ax.set_title(title, fontsize=12, pad=8)
        ax.axis("off")
    fig.tight_layout()
    out = OUT / "figures" / "Fig_NeoICI_nature_cancer_synthesis_overview.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    ensure_dir(OUT / "figures")
    ensure_dir(HUB)
    ensure_dir(WWW)

    m = load_metrics()
    fig_path = make_overview_figure()

    summary = [
        "# NeoICI Nature Cancer synthesis",
        "",
        "## What this package now shows",
        "",
        f"- NeoICI is not the strongest universal scorer yet; `KG_GA_evolved` still leads the local validation-like split at AUPRC {m['kg_ga_validation_auprc']:.3f}.",
        f"- `NeoICI_GA_RL_combo_v1` is strong but lower at AUPRC {m['neoici_validation_auprc']:.3f}.",
        f"- Industrial locked still favors `BioDarwin_mode_bank_v4` at AUPRC {m['industrial_biodarwin_auprc']:.3f}.",
        "",
        "## Combo-replay evidence",
        "",
        f"- The open combo GEO replay retained {m['combo_n']} paired samples across GSE222011 and GSE255830.",
        f"- Median response-like index was {m['combo_response_median']:.3f}; median top-clone fraction was {m['combo_clone_median']:.3f}.",
        "",
        "## Proxy comparison",
        "",
        f"- NeoICI-style proxy Spearman vs response-like index: {m['neoici_proxy_rho']:.3f}.",
        f"- NeoPrecis-style proxy Spearman vs response-like index: {m['neoprecis_proxy_rho']:.3f}.",
        f"- NeoICI-style proxy Spearman vs top-clone fraction: {m['neoici_proxy_rho_clone']:.3f}.",
        f"- NeoPrecis-style proxy Spearman vs top-clone fraction: {m['neoprecis_proxy_rho_clone']:.3f}.",
        "",
        "## Write-up rule",
        "",
        "- Use NeoICI as the combo-therapy readiness layer.",
        "- Do not claim universal superiority over NeoPrecis or any published benchmark.",
        "- If a figure is needed for the manuscript, the synthesis overview should be the top-level entry point, with the proxy comparison as the supporting panel.",
    ]
    (OUT / "NEOICI_NATURE_CANCER_SYNTHESIS_SUMMARY.md").write_text("\n".join(summary) + "\n", encoding="utf-8")

    contract = pd.DataFrame(
        [
            {
                "layer": "algorithm gauntlet",
                "key_readout": f"KG_GA_evolved {m['kg_ga_validation_auprc']:.3f} vs NeoICI {m['neoici_validation_auprc']:.3f}",
                "interpretation": "NeoICI is strong, but not the top local scorer.",
            },
            {
                "layer": "open combo replay",
                "key_readout": f"n={m['combo_n']} paired samples; response median {m['combo_response_median']:.3f}",
                "interpretation": "Public combo cohorts show a measurable response-like phenotype.",
            },
            {
                "layer": "proxy comparison",
                "key_readout": f"NeoICI rho {m['neoici_proxy_rho']:.3f} vs NeoPrecis rho {m['neoprecis_proxy_rho']:.3f}",
                "interpretation": "NeoICI-style proxy tracks the response-like index better in this tiny replay.",
            },
        ]
    )
    contract.to_csv(OUT / "NEOICI_NATURE_CANCER_SYNTHESIS_TABLE.tsv", sep="\t", index=False)

    combo = pd.read_csv(RESULT_B / "combo_geo_reanalysis.tsv", sep="\t")
    combo = combo.sort_values("response_index", ascending=False).copy()
    combo["rank_response"] = range(1, len(combo) + 1)
    import numpy as np
    ent = combo["tcr_entropy"].fillna(combo["tcr_entropy"].median())
    ent_scaled = (ent - ent.min()) / (ent.max() - ent.min() + 1e-9)
    combo["NeoPrecis_style_proxy"] = (
        0.45 * combo["top_clone_frac"]
        + 0.20 * combo["cytotoxic_score"]
        + 0.20 * (1 - combo["exhaustion_score"])
        + 0.15 * ent_scaled
    )
    combo["NeoICI_style_proxy"] = (
        0.45 * combo["cytotoxic_score"]
        + 0.25 * combo["activation_score"]
        + 0.15 * combo["memory_score"]
        + 0.15 * (1 - combo["exhaustion_score"])
    )
    combo["rank_neoici"] = combo["NeoICI_style_proxy"].rank(ascending=False, method="average").astype(float)
    combo["rank_neoprecis"] = combo["NeoPrecis_style_proxy"].rank(ascending=False, method="average").astype(float)
    sample_table = combo[[
        "dataset",
        "sample",
        "response_index",
        "rank_response",
        "NeoICI_style_proxy",
        "rank_neoici",
        "NeoPrecis_style_proxy",
        "rank_neoprecis",
        "top_clone_frac",
        "cytotoxic_score",
        "activation_score",
        "exhaustion_score",
    ]].copy()
    sample_table.to_csv(OUT / "NEOICI_NATURE_CANCER_SAMPLE_TABLE.tsv", sep="\t", index=False)

    from scipy.stats import spearmanr
    loo_rows = []
    for name in ["NeoPrecis_style_proxy", "NeoICI_style_proxy"]:
        vals = []
        for i in range(len(combo)):
            sub = combo.drop(index=combo.index[i])
            vals.append(float(spearmanr(sub[name], sub["response_index"]).statistic))
        loo_rows.append(
            {
                "proxy": name,
                "loo_min": float(np.min(vals)),
                "loo_median": float(np.median(vals)),
                "loo_max": float(np.max(vals)),
                "loo_sd": float(np.std(vals, ddof=1)),
            }
        )
    loo_table = pd.DataFrame(loo_rows)
    loo_table.to_csv(OUT / "NEOICI_NATURE_CANCER_LOO_TABLE.tsv", sep="\t", index=False)

    claim = pd.DataFrame(
        [
            {
                "status": "Allowed now",
                "claim": "NeoICI is a combo-therapy readiness layer that can be benchmarked against open vaccine+ICI cohorts.",
                "evidence": "open combo GEO replay + proxy comparison",
            },
            {
                "status": "Allowed now",
                "claim": "The open combo replay shows a measurable response-like phenotype and clonotype expansion signal.",
                "evidence": f"n={m['combo_n']} paired samples; response median {m['combo_response_median']:.3f}",
            },
            {
                "status": "Allowed with caution",
                "claim": "NeoICI-style proxy tracks the response-like index better than the NeoPrecis-style proxy in this tiny replay.",
                "evidence": f"Spearman {m['neoici_proxy_rho']:.3f} vs {m['neoprecis_proxy_rho']:.3f}; leave-one-out median {loo_table.query('proxy == \"NeoICI_style_proxy\"')['loo_median'].iloc[0]:.3f}",
            },
            {
                "status": "Forbidden now",
                "claim": "Universal superiority over NeoPrecis or any published benchmark.",
                "evidence": "n=9 heuristic proxies only",
            },
            {
                "status": "Forbidden now",
                "claim": "Clinical-grade vaccine or ICI treatment selection.",
                "evidence": "post-hoc replay only",
            },
        ]
    )
    claim.to_csv(OUT / "NEOICI_NATURE_CANCER_CLAIM_TABLE.tsv", sep="\t", index=False)

    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>NeoICI Nature Cancer Synthesis</title>
  <style>
    body {{ margin:0; background:#0b1020; color:#e5e7eb; font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; }}
    .wrap {{ max-width: 1240px; margin: 0 auto; padding: 28px 22px 48px; }}
    h1, h2 {{ margin:0 0 12px; }}
    .lead {{ max-width: 84ch; color:#cbd5e1; line-height:1.6; margin: 0 0 18px; }}
    .grid {{ display:grid; grid-template-columns: 1.15fr .85fr; gap:18px; align-items:start; }}
    .panel {{ border:1px solid #243047; background:#0f172a; padding:16px; }}
    table {{ width:100%; border-collapse:collapse; font-size:14px; margin-top:14px; }}
    th, td {{ border:1px solid #243047; padding:10px; vertical-align:top; }}
    th {{ background:#111827; text-align:left; }}
    .hero {{ margin-top:16px; }}
    .hero img {{ width:100%; display:block; border:1px solid #243047; background:#0b1020; }}
    .small {{ color:#94a3b8; font-size:13px; line-height:1.5; }}
    code {{ background:#111827; padding:1px 4px; border-radius:4px; }}
  </style>
</head>
<body>
<div class="wrap">
  <h1>NeoICI Nature Cancer synthesis</h1>
  <p class="lead">This synthesis page merges the local algorithm gauntlet, the open combo GEO replay, and the proxy comparison into one reviewer-facing summary. The result is intentionally conservative: NeoICI is a strong combo-therapy readiness layer, not a universal winner over NeoPrecis or the existing local champions.</p>
  <div class="hero">
    <img src="assets/neoici_nature_cancer_synthesis/Fig_NeoICI_nature_cancer_synthesis_overview.png" alt="NeoICI synthesis overview" />
  </div>
  <div class="grid" style="margin-top:18px;">
    <div class="panel">
      <h2>Summary table</h2>
      {contract.to_html(index=False, escape=False)}
      <h2 style="margin-top:18px;">Claim boundary</h2>
      {claim.to_html(index=False, escape=False)}
    </div>
    <div class="panel">
      <h2>Decision boundary</h2>
      <p class="small">Use the NeoICI package as the combo-therapy readiness layer for Nature Cancer framing. Keep the published NeoPrecis contract as the external benchmark reference, and keep <code>KG_GA_evolved</code> / <code>BioDarwin_mode_bank_v4</code> as local upper-bound comparators.</p>
      <p class="small">The open combo replay is now a real phenotype support layer: 9 paired samples, measurable response-like index, and a proxy comparison that favors the NeoICI-style T-cell-state view on this tiny set.</p>
      <p class="small">Do not upgrade this to universal superiority. The correct claim is that NeoICI adds the combo-ICI axis that the older single-objective benchmark does not directly own.</p>
    </div>
  </div>
  <div class="panel" style="margin-top:18px;">
    <h2>Sample ranking</h2>
    {sample_table.to_html(index=False, escape=False)}
    <h2 style="margin-top:18px;">Leave-one-out robustness</h2>
    {loo_table.to_html(index=False, escape=False)}
  </div>
</div>
</body>
</html>
"""
    html_path = HUB / "neoici_nature_cancer_synthesis.html"
    html_path.write_text(html, encoding="utf-8")
    (WWW / "neoici_nature_cancer_synthesis.html").write_text(html, encoding="utf-8")

    asset_dir = HUB / "assets" / "neoici_nature_cancer_synthesis"
    asset_www = WWW / "assets" / "neoici_nature_cancer_synthesis"
    ensure_dir(asset_dir)
    ensure_dir(asset_www)
    (asset_dir / fig_path.name).write_bytes(fig_path.read_bytes())
    (asset_www / fig_path.name).write_bytes(fig_path.read_bytes())

    meta = {
        "html": str(html_path),
        "summary": str(OUT / "NEOICI_NATURE_CANCER_SYNTHESIS_SUMMARY.md"),
        "table": str(OUT / "NEOICI_NATURE_CANCER_SYNTHESIS_TABLE.tsv"),
        "sample_table": str(OUT / "NEOICI_NATURE_CANCER_SAMPLE_TABLE.tsv"),
        "loo_table": str(OUT / "NEOICI_NATURE_CANCER_LOO_TABLE.tsv"),
        "figure": str(fig_path),
    }
    (OUT / "neoici_nature_cancer_synthesis.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
