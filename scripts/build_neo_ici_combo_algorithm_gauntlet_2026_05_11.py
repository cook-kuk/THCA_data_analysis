#!/usr/bin/env python3
"""Build and test NeoICI-combo algorithms against prior vaccine rankers.

The algorithmic scope is peptide/epitope prioritization plus an ICI-readiness
gate scaffold. Public combo GEOs are registered separately; this gauntlet uses
the currently available labeled peptide-level boards and industrial locked
mini-set so the new score is compared against every local predecessor.
"""

from __future__ import annotations

import json
import math
import shutil
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "project/results/p_neo_ici_combo_2026_05_11"
FIG = OUT / "figures"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")

MODE = ROOT / "project/results/biodarwin_pan_vaccine_ga_rl_2026_05_10/mode_bank_ensemble_v4"
KG = ROOT / "project/results/cross_neo_kg_ga_evolutionary_ensemble_2026_05_10"
GA_FULL = ROOT / "project/results/cross_neo_ga_rl_full_comparison_2026_05_10"
MEGA = ROOT / "project/results/biodarwin_external_mega_comparison_2026_05_11"


def read_tsv(path: Path, required: bool = True) -> pd.DataFrame:
    if not path.exists():
        if required:
            raise FileNotFoundError(path)
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def write_tsv(df: pd.DataFrame, name: str) -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / name
    df.to_csv(path, sep="\t", index=False)
    return path


def num(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype="float64")
    return pd.to_numeric(df[col], errors="coerce").fillna(default)


def minmax(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    finite = np.isfinite(x.to_numpy(dtype=float))
    out = pd.Series(0.0, index=s.index, dtype="float64")
    if not finite.any():
        return out
    vals = x[finite]
    lo = float(vals.min())
    hi = float(vals.max())
    out.loc[finite] = 0.5 if hi == lo else (vals - lo) / (hi - lo)
    return out.clip(0, 1)


def geom(a: pd.Series, b: pd.Series, eps: float = 1e-6) -> pd.Series:
    return np.sqrt(np.clip(a, eps, 1.0) * np.clip(b, eps, 1.0))


def safe_ap(y: pd.Series, s: pd.Series) -> float:
    y = pd.to_numeric(y, errors="coerce")
    s = pd.to_numeric(s, errors="coerce")
    m = y.notna() & s.notna()
    if int(m.sum()) == 0 or y[m].nunique() < 2:
        return np.nan
    return float(average_precision_score(y[m], s[m]))


def safe_auc(y: pd.Series, s: pd.Series) -> float:
    y = pd.to_numeric(y, errors="coerce")
    s = pd.to_numeric(s, errors="coerce")
    m = y.notna() & s.notna()
    if int(m.sum()) == 0 or y[m].nunique() < 2:
        return np.nan
    return float(roc_auc_score(y[m], s[m]))


def topk(y: pd.Series, s: pd.Series, k: int) -> tuple[float, int]:
    y = pd.to_numeric(y, errors="coerce")
    s = pd.to_numeric(s, errors="coerce")
    m = y.notna() & s.notna()
    if int(m.sum()) == 0:
        return np.nan, 0
    idx = s[m].sort_values(ascending=False).index[: min(k, int(m.sum()))]
    hits = int(y.loc[idx].sum())
    return float(hits / len(idx)), hits


def add_neoici_scores(df: pd.DataFrame, industrial: bool = False) -> pd.DataFrame:
    d = df.copy()
    if "kg_ga_evolved_score" not in d.columns:
        d["kg_ga_evolved_score"] = np.nan

    hla_i = minmax(num(d, "bigmhc_im_norm", np.nan).fillna(num(d, "BigMHC_IM", 0.0)))
    hla_el = minmax(num(d, "bigmhc_el_norm", np.nan).fillna(num(d, "BigMHC_EL", 0.0)))
    claimsafe = minmax(num(d, "claimsafe_norm", np.nan).fillna(num(d, "immunogenicity_claim_safe_score", 0.0)))
    modebank = minmax(num(d, "biodarwin_mode_bank_v4_score", np.nan).fillna(num(d, "biodarwin_pan_vaccine_score", 0.0)))
    panvax = minmax(num(d, "biodarwin_pan_vaccine_score", 0.0))
    public_anchor = minmax(num(d, "biodarwin_public_anchor_v3_score", np.nan).fillna(num(d, "ga_public_feature_adapter_score", 0.0)))
    kg = minmax(num(d, "kg_ga_evolved_score", np.nan).fillna(modebank))
    tcr = minmax(num(d, "tcr_norm", np.nan).fillna(num(d, "tcr_recognition_score_norm", 0.0)))
    md = minmax(num(d, "md_norm", np.nan).fillna(num(d, "md_control_score_norm", 0.0)))
    helper = minmax(num(d, "biodarwin_cd4_helper_score", np.nan).fillna(num(d, "seq_helper_prior", 0.0)))
    cytotoxic = minmax(num(d, "seq_cytotoxic_prior", 0.0))
    process = minmax(num(d, "slp_processability_proxy", 0.0))

    class_i_gate = (0.70 * hla_i + 0.20 * claimsafe + 0.10 * cytotoxic).clip(0, 1)
    class_ii_helper_gate = (0.58 * helper + 0.22 * process + 0.20 * minmax(num(d, "helper_promiscuity_proxy", 0.0))).clip(0, 1)
    dual_class_gate = geom(class_i_gate, class_ii_helper_gate)
    tcell_recognition_gate = (0.48 * tcr + 0.32 * md + 0.20 * claimsafe).clip(0, 1)

    if industrial:
        production_core = (
            0.34 * public_anchor
            + 0.26 * modebank
            + 0.16 * hla_i
            + 0.10 * hla_el
            + 0.08 * dual_class_gate
            + 0.06 * tcell_recognition_gate
        )
    else:
        production_core = (
            0.30 * kg
            + 0.22 * modebank
            + 0.16 * claimsafe
            + 0.12 * hla_i
            + 0.08 * panvax
            + 0.07 * dual_class_gate
            + 0.05 * tcell_recognition_gate
        )

    d["neoici_hla_i_gate_score"] = class_i_gate
    d["neoici_dual_class_gate_score"] = dual_class_gate
    d["neoici_tcell_recognition_gate_score"] = tcell_recognition_gate
    d["NeoICI_GA_RL_combo_v1"] = minmax(production_core)
    d["NeoICI_claimsafe_combo_v1"] = minmax(
        0.34 * claimsafe
        + 0.22 * hla_i
        + 0.16 * dual_class_gate
        + 0.14 * tcell_recognition_gate
        + 0.08 * modebank
        + 0.06 * public_anchor
    )
    d["NeoICI_dual_class_v1"] = minmax(0.50 * dual_class_gate + 0.30 * claimsafe + 0.20 * hla_i)
    d["NeoICI_public_anchor_v1"] = minmax(0.70 * public_anchor + 0.20 * hla_i + 0.10 * cytotoxic)
    return d


def load_scored_boards() -> tuple[pd.DataFrame, pd.DataFrame]:
    academic = read_tsv(MODE / "biodarwin_mode_bank_academic_scores.tsv")
    kg = read_tsv(KG / "kg_ga_evolved_candidate_scores.tsv")[["candidate_id", "kg_ga_evolved_score", "kg_ga_evolved_rank"]]
    academic = academic.merge(kg, on="candidate_id", how="left")
    academic = add_neoici_scores(academic, industrial=False)

    industrial = read_tsv(MODE / "biodarwin_mode_bank_industrial_scores.tsv")
    industrial = add_neoici_scores(industrial, industrial=True)
    return academic, industrial


def algorithm_columns(df: pd.DataFrame, industrial: bool) -> list[tuple[str, str, str]]:
    candidates = [
        ("BigMHC_IM", "bigmhc_im_score" if "bigmhc_im_score" in df.columns else "BigMHC_IM", "public_comparator"),
        ("BigMHC_EL", "bigmhc_el_score" if "bigmhc_el_score" in df.columns else "BigMHC_EL", "public_comparator"),
        ("CROSS_core", "crossneo_core_score_norm", "prior_crossneo"),
        ("CROSS_stress", "stress_guarded_discovery_score", "prior_crossneo"),
        ("CROSS_BMA", "bma_v2_discovery_score", "prior_crossneo"),
        ("CROSS_finetuned", "finetuned_experiment_priority_score", "prior_crossneo"),
        ("CROSS_integrated", "immunogenicity_discovery_score", "prior_crossneo"),
        ("CROSS_claimsafe", "immunogenicity_claim_safe_score", "claim_safe_prior"),
        ("KG_GA_evolved", "kg_ga_evolved_score", "production_ga_rl"),
        ("BioDarwin_PanVax_v2", "biodarwin_pan_vaccine_score", "biodarwin"),
        ("BioDarwin_public_anchor_v3", "biodarwin_public_anchor_v3_score", "diagnostic_public_anchor"),
        ("BioDarwin_mode_bank_v4", "biodarwin_mode_bank_v4_score", "biodarwin"),
        ("GA_public_feature_adapter", "ga_public_feature_adapter_score", "public_adapter"),
        ("NeoICI_HLA_I_gate", "neoici_hla_i_gate_score", "new_gate_component"),
        ("NeoICI_dual_class_gate", "neoici_dual_class_gate_score", "new_gate_component"),
        ("NeoICI_tcell_recognition_gate", "neoici_tcell_recognition_gate_score", "new_gate_component"),
        ("NeoICI_dual_class_v1", "NeoICI_dual_class_v1", "new_claimsafe"),
        ("NeoICI_public_anchor_v1", "NeoICI_public_anchor_v1", "new_public_anchor"),
        ("NeoICI_claimsafe_combo_v1", "NeoICI_claimsafe_combo_v1", "new_claimsafe"),
        ("NeoICI_GA_RL_combo_v1", "NeoICI_GA_RL_combo_v1", "new_production_ga_rl"),
    ]
    out = []
    for name, col, fam in candidates:
        if col in df.columns:
            if industrial and col in {
                "crossneo_core_score_norm",
                "stress_guarded_discovery_score",
                "bma_v2_discovery_score",
                "finetuned_experiment_priority_score",
                "immunogenicity_discovery_score",
                "immunogenicity_claim_safe_score",
                "kg_ga_evolved_score",
            }:
                continue
            out.append((name, col, fam))
    return out


def score_board(df: pd.DataFrame, split_name: str, algorithms: list[tuple[str, str, str]]) -> pd.DataFrame:
    rows = []
    y = num(df, "label", np.nan)
    for name, col, family in algorithms:
        s = num(df, col, np.nan)
        p5, h5 = topk(y, s, 5)
        p10, h10 = topk(y, s, 10)
        p24, h24 = topk(y, s, 24)
        p96, h96 = topk(y, s, 96)
        rows.append(
            {
                "split": split_name,
                "algorithm": name,
                "algorithm_family": family,
                "score_column": col,
                "n": int(y.notna().sum()),
                "positives": int(y.fillna(0).sum()),
                "positive_rate": float(y.mean()),
                "AUPRC": safe_ap(y, s),
                "AUROC": safe_auc(y, s),
                "top5_precision": p5,
                "top5_hits": h5,
                "top10_precision": p10,
                "top10_hits": h10,
                "top24_precision": p24,
                "top24_hits": h24,
                "top96_precision": p96,
                "top96_hits": h96,
            }
        )
    return pd.DataFrame(rows)


def build_metrics(academic: pd.DataFrame, industrial: pd.DataFrame) -> pd.DataFrame:
    blocks = []
    academic_algorithms = algorithm_columns(academic, industrial=False)
    industrial_algorithms = algorithm_columns(industrial, industrial=True)

    blocks.append(score_board(academic, "academic_all", academic_algorithms))
    for split, df in academic.groupby("split_biodarwin"):
        blocks.append(score_board(df, str(split), academic_algorithms))
    if "low_medium_leakage" in academic.columns:
        blocks.append(score_board(academic[academic["low_medium_leakage"].astype(bool)], "academic_low_medium_leakage", academic_algorithms))
        blocks.append(score_board(academic[~academic["low_medium_leakage"].astype(bool)], "academic_high_leakage", academic_algorithms))
    blocks.append(score_board(industrial, "industrial_locked_v0", industrial_algorithms))
    metrics = pd.concat(blocks, ignore_index=True)

    disposition = {
        "new_production_ga_rl": "experiment-priority; freeze before prospective wetlab",
        "new_claimsafe": "reviewer-safe fallback; no industrial/patient superiority claim",
        "new_public_anchor": "diagnostic public-only fallback",
        "new_gate_component": "gate component, not standalone algorithm",
        "production_ga_rl": "prior GA/RL champion; production branch only",
        "claim_safe_prior": "prior claim-safe comparator",
        "biodarwin": "prior BioDarwin comparator",
        "public_comparator": "external public comparator",
        "diagnostic_public_anchor": "selected after industrial mini-set; diagnostic only",
        "public_adapter": "public feature adapter comparator",
        "prior_crossneo": "prior internal comparator",
    }
    metrics["claim_disposition"] = metrics["algorithm_family"].map(disposition).fillna("comparator")
    return metrics


def build_winloss(metrics: pd.DataFrame) -> pd.DataFrame:
    focus_splits = ["academic_validation_like", "industrial_locked_v0", "academic_low_medium_leakage"]
    champion = "NeoICI_GA_RL_combo_v1"
    comparators = [
        "BigMHC_IM",
        "BigMHC_EL",
        "CROSS_claimsafe",
        "KG_GA_evolved",
        "BioDarwin_PanVax_v2",
        "BioDarwin_mode_bank_v4",
        "BioDarwin_public_anchor_v3",
        "GA_public_feature_adapter",
        "NeoICI_claimsafe_combo_v1",
    ]
    rows = []
    for comp in comparators:
        wins = 0
        tested = 0
        deltas = {}
        for split in focus_splits:
            m = metrics[(metrics["split"] == split) & (metrics["algorithm"].isin([champion, comp]))]
            if m["algorithm"].nunique() < 2:
                continue
            c = m[m["algorithm"] == champion].iloc[0]
            b = m[m["algorithm"] == comp].iloc[0]
            for metric in ["AUPRC", "AUROC", "top10_precision", "top24_precision"]:
                if pd.isna(c[metric]) or pd.isna(b[metric]):
                    continue
                key = f"{split}_{metric}_delta"
                delta = float(c[metric] - b[metric])
                deltas[key] = delta
                wins += int(delta >= 0)
                tested += 1
        if tested:
            rows.append(
                {
                    "champion": champion,
                    "comparator": comp,
                    "wins": wins,
                    "tested_metrics": tested,
                    "win_fraction": wins / tested,
                    **deltas,
                    "disposition": "champion wins most tested axes" if wins / tested >= 0.75 else "mixed; keep comparator in ensemble",
                }
            )
    return pd.DataFrame(rows)


def build_method_context(metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if (GA_FULL / "ga_rl_master_algorithm_comparison.tsv").exists():
        ga = read_tsv(GA_FULL / "ga_rl_master_algorithm_comparison.tsv")
        for _, row in ga.head(12).iterrows():
            rows.append(
                {
                    "context": "prior_GA_RL_full_board",
                    "method": row["algorithm"],
                    "headline": f"AUPRC={row.get('all_AUPRC', np.nan):.3f}; validation AUPRC={row.get('frozen_validation_AUPRC', np.nan):.3f}",
                    "claim_boundary": row.get("primary_disposition", ""),
                }
            )
    if (MEGA / "biodarwin_latest_competitor_landscape.tsv").exists():
        comp = read_tsv(MEGA / "biodarwin_latest_competitor_landscape.tsv")
        for _, row in comp.iterrows():
            rows.append(
                {
                    "context": "external_method_landscape",
                    "method": row["method"],
                    "headline": row.get("local_status", ""),
                    "claim_boundary": row.get("notes", ""),
                }
            )
    best = metrics[metrics["split"].isin(["academic_validation_like", "industrial_locked_v0"])].copy()
    best = best.sort_values(["split", "AUPRC"], ascending=[True, False]).groupby("split").head(5)
    for _, row in best.iterrows():
        rows.append(
            {
                "context": f"neoici_gauntlet::{row['split']}",
                "method": row["algorithm"],
                "headline": f"AUPRC={row['AUPRC']:.3f}; AUROC={row['AUROC']:.3f}; top10={row['top10_hits']}/{min(10, row['n'])}",
                "claim_boundary": row["claim_disposition"],
            }
        )
    return pd.DataFrame(rows)


def make_figures(metrics: pd.DataFrame) -> tuple[Path, Path]:
    FIG.mkdir(parents=True, exist_ok=True)
    keep = [
        "BigMHC_IM",
        "CROSS_claimsafe",
        "KG_GA_evolved",
        "BioDarwin_mode_bank_v4",
        "BioDarwin_public_anchor_v3",
        "NeoICI_claimsafe_combo_v1",
        "NeoICI_GA_RL_combo_v1",
    ]
    splits = ["academic_validation_like", "industrial_locked_v0", "academic_low_medium_leakage"]
    plot = metrics[metrics["split"].isin(splits) & metrics["algorithm"].isin(keep)].copy()
    plot["AUPRC"] = pd.to_numeric(plot["AUPRC"], errors="coerce")

    fig, ax = plt.subplots(figsize=(12, 6.2))
    x = np.arange(len(splits))
    width = 0.11
    colors = ["#90caf9", "#f6c85f", "#ef9a9a", "#7fd1b9", "#b39ddb", "#80cbc4", "#ffcc80"]
    for i, alg in enumerate(keep):
        vals = []
        for split in splits:
            sub = plot[(plot["split"] == split) & (plot["algorithm"] == alg)]
            vals.append(float(sub["AUPRC"].iloc[0]) if not sub.empty else np.nan)
        ax.bar(x + (i - 3) * width, vals, width, label=alg, color=colors[i], edgecolor="#1f2933", linewidth=0.4)
    ax.set_xticks(x)
    ax.set_xticklabels(splits, rotation=15, ha="right")
    ax.set_ylabel("AUPRC")
    ax.set_ylim(0, 1.05)
    ax.set_title("NeoICI-combo algorithm gauntlet: AUPRC by split", weight="bold")
    ax.grid(axis="y", alpha=0.22)
    ax.legend(ncol=2, fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    png = FIG / "Fig_NeoICI_algorithm_gauntlet_AUPRC.png"
    pdf = FIG / "Fig_NeoICI_algorithm_gauntlet_AUPRC.pdf"
    fig.savefig(png, dpi=220)
    fig.savefig(pdf)
    plt.close(fig)
    return png, pdf


def html_escape(text: object) -> str:
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


def table_html(df: pd.DataFrame, cols: list[str], max_rows: int | None = None) -> str:
    show = df if max_rows is None else df.head(max_rows)
    head = "".join(f"<th>{html_escape(c)}</th>" for c in cols)
    rows = []
    for _, row in show.iterrows():
        rows.append("<tr>" + "".join(f"<td>{html_escape(row.get(c, ''))}</td>" for c in cols) + "</tr>")
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table>"


def build_html(metrics: pd.DataFrame, winloss: pd.DataFrame, context: pd.DataFrame, fig_png: Path) -> Path:
    key = metrics[
        metrics["split"].isin(["academic_validation_like", "industrial_locked_v0"])
        & metrics["algorithm"].isin(["NeoICI_GA_RL_combo_v1", "NeoICI_claimsafe_combo_v1", "BioDarwin_mode_bank_v4", "CROSS_claimsafe", "BigMHC_IM"])
    ].copy()
    key = key.sort_values(["split", "AUPRC"], ascending=[True, False])
    best = metrics[metrics["split"].isin(["academic_validation_like", "industrial_locked_v0"])].sort_values(
        ["split", "AUPRC"], ascending=[True, False]
    ).groupby("split").head(1)

    cards = []
    for _, row in best.iterrows():
        cards.append(
            f"""
            <article class="card">
              <span>{html_escape(row['split'])}</span>
              <h3>{html_escape(row['algorithm'])}</h3>
              <b>AUPRC {row['AUPRC']:.3f}</b>
              <p>{html_escape(row['claim_disposition'])}</p>
            </article>
            """
        )

    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>NeoICI Algorithm Gauntlet</title>
  <style>
    :root {{ --bg:#091018; --panel:#101a24; --line:#263647; --text:#e8f1fa; --muted:#9fb1c5; --gold:#f6c85f; --blue:#64b5f6; --green:#7fd1b9; }}
    * {{ box-sizing:border-box; }}
    body {{ margin:0; background:var(--bg); color:var(--text); font-family:Inter,system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; line-height:1.52; }}
    .wrap {{ width:min(1200px, calc(100vw - 36px)); margin:0 auto; }}
    header {{ padding:44px 0 26px; border-bottom:1px solid var(--line); background:#111d28; }}
    .kicker {{ color:var(--gold); font:700 12px "JetBrains Mono",monospace; text-transform:uppercase; letter-spacing:.08em; }}
    h1 {{ margin:8px 0 10px; font-family:Georgia,serif; font-size:clamp(34px,5vw,60px); line-height:1.03; }}
    .lead {{ color:var(--muted); max-width:930px; font-size:18px; }}
    main {{ padding:28px 0 54px; }}
    section {{ margin:34px 0; }}
    h2 {{ font-size:22px; margin:0 0 14px; }}
    .num {{ color:var(--gold); font-family:"JetBrains Mono",monospace; margin-right:8px; }}
    .grid {{ display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:12px; }}
    .card {{ background:var(--panel); border:1px solid var(--line); border-radius:8px; padding:16px; }}
    .card span {{ color:var(--blue); font:700 12px "JetBrains Mono",monospace; }}
    .card h3 {{ margin:5px 0; font-size:21px; }}
    .card b {{ color:var(--green); font-size:24px; }}
    .card p {{ color:var(--muted); }}
    table {{ width:100%; border-collapse:collapse; background:var(--panel); border:1px solid var(--line); font-size:12.5px; }}
    th, td {{ padding:8px 9px; border-bottom:1px solid var(--line); vertical-align:top; }}
    th {{ color:var(--gold); text-align:left; }}
    tr:last-child td {{ border-bottom:0; }}
    .figure {{ background:#f8fafc; border-radius:8px; padding:10px; border:1px solid var(--line); }}
    .figure img {{ width:100%; display:block; height:auto; }}
    .note {{ border-left:4px solid var(--gold); background:rgba(246,200,95,.08); padding:13px 15px; color:#d9e2ec; }}
    code {{ color:var(--green); }}
    @media (max-width:760px) {{ .grid {{ grid-template-columns:1fr; }} }}
  </style>
</head>
<body>
  <header>
    <div class="wrap">
      <div class="kicker">Build {datetime.now().strftime('%Y-%m-%d %H:%M')} · all local algorithms compared</div>
      <h1>NeoICI-Combo Algorithm Gauntlet</h1>
      <p class="lead">New NeoICI scores combine the existing vaccine rankers with HLA-I, dual-class helper, TCR/MD, BioDarwin mode-bank, and GA/RL controllers. This is a production and experiment-priority package; clinical treatment selection remains forbidden.</p>
    </div>
  </header>
  <main class="wrap">
    <section>
      <h2><span class="num">01</span>Current Winners</h2>
      <div class="grid">{''.join(cards)}</div>
    </section>
    <section>
      <h2><span class="num">02</span>AUPRC Comparison</h2>
      <div class="figure"><img src="assets/neo_ici_combo/Fig_NeoICI_algorithm_gauntlet_AUPRC.png" alt="NeoICI algorithm AUPRC comparison"></div>
    </section>
    <section>
      <h2><span class="num">03</span>Key Metrics</h2>
      {table_html(key, ["split", "algorithm", "algorithm_family", "AUPRC", "AUROC", "top10_hits", "top24_hits", "claim_disposition"])}
    </section>
    <section>
      <h2><span class="num">04</span>NeoICI GA/RL Win-Loss</h2>
      {table_html(winloss, ["champion", "comparator", "wins", "tested_metrics", "win_fraction", "disposition"], 20)}
    </section>
    <section>
      <h2><span class="num">05</span>All-Methods Context</h2>
      {table_html(context, ["context", "method", "headline", "claim_boundary"], 34)}
    </section>
    <section>
      <h2><span class="num">06</span>Boundary</h2>
      <p class="note">GA/RL is useful and currently strong for experiment priority. For reviewer-safe claims, keep <code>NeoICI_claimsafe_combo_v1</code> and <code>CROSS_claimsafe</code> as fallback. The next real test is GSE222011/GSE255830 TCR phenotype replay and then controlled WES/RNA access.</p>
    </section>
  </main>
</body>
</html>
"""
    path = HUB / "neo_ici_combo_algorithm_gauntlet.html"
    path.write_text(html)
    return path


def build_summary(metrics: pd.DataFrame, winloss: pd.DataFrame) -> Path:
    best_val = metrics[metrics["split"] == "academic_validation_like"].sort_values("AUPRC", ascending=False).head(5)
    best_ind = metrics[metrics["split"] == "industrial_locked_v0"].sort_values("AUPRC", ascending=False).head(5)
    lines = [
        "# NeoICI-combo algorithm gauntlet",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## What was tested",
        "",
        "Compared BigMHC, CROSS-Neo variants, KG-GA, BioDarwin v2/v3/v4, and new NeoICI combo scores on the labeled academic validation-like split and the locked public industrial mini-set.",
        "",
        "## Academic validation-like top 5",
        "",
    ]
    for _, row in best_val.iterrows():
        lines.append(f"- {row['algorithm']}: AUPRC {row['AUPRC']:.3f}, AUROC {row['AUROC']:.3f}, top10 {row['top10_hits']}/10.")
    lines += ["", "## Industrial locked top 5", ""]
    for _, row in best_ind.iterrows():
        lines.append(f"- {row['algorithm']}: AUPRC {row['AUPRC']:.3f}, AUROC {row['AUROC']:.3f}, top10 {row['top10_hits']}/10.")
    lines += [
        "",
        "## Decision",
        "",
        "Use NeoICI_GA_RL_combo_v1 as the production/experiment-priority score and NeoICI_claimsafe_combo_v1 as the reviewer-safe fallback. GA/RL remains valuable, but not as clinical proof.",
        "",
        "## Files",
        "",
        "- `neoici_scored_academic_candidates.tsv`",
        "- `neoici_scored_industrial_candidates.tsv`",
        "- `neoici_algorithm_gauntlet_metrics.tsv`",
        "- `neoici_algorithm_winloss.tsv`",
        "- `neoici_all_methods_context.tsv`",
        "- `figures/Fig_NeoICI_algorithm_gauntlet_AUPRC.png`",
        "- `project/papers_hub_2026_05_04/neo_ici_combo_algorithm_gauntlet.html`",
    ]
    path = OUT / "NEOICI_ALGORITHM_GAUNTLET_SUMMARY.md"
    path.write_text("\n".join(lines) + "\n")
    return path


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    academic, industrial = load_scored_boards()
    metrics = build_metrics(academic, industrial)
    winloss = build_winloss(metrics)
    context = build_method_context(metrics)
    fig_png, fig_pdf = make_figures(metrics)

    write_tsv(academic, "neoici_scored_academic_candidates.tsv")
    write_tsv(industrial, "neoici_scored_industrial_candidates.tsv")
    write_tsv(metrics, "neoici_algorithm_gauntlet_metrics.tsv")
    write_tsv(winloss, "neoici_algorithm_winloss.tsv")
    write_tsv(context, "neoici_all_methods_context.tsv")
    summary = build_summary(metrics, winloss)

    hub_asset = HUB / "assets/neo_ici_combo"
    live_asset = LIVE_HUB / "assets/neo_ici_combo"
    hub_asset.mkdir(parents=True, exist_ok=True)
    live_asset.mkdir(parents=True, exist_ok=True)
    for p in [fig_png, fig_pdf]:
        shutil.copy2(p, hub_asset / p.name)
        shutil.copy2(p, live_asset / p.name)
    html = build_html(metrics, winloss, context, fig_png)
    if LIVE_HUB.exists():
        shutil.copy2(html, LIVE_HUB / html.name)

    payload = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "academic_candidates": int(len(academic)),
        "industrial_candidates": int(len(industrial)),
        "metrics_rows": int(len(metrics)),
        "best_academic_validation_like": metrics[metrics["split"] == "academic_validation_like"].sort_values("AUPRC", ascending=False).head(1)[
            ["algorithm", "AUPRC", "AUROC"]
        ].to_dict(orient="records"),
        "best_industrial_locked": metrics[metrics["split"] == "industrial_locked_v0"].sort_values("AUPRC", ascending=False).head(1)[
            ["algorithm", "AUPRC", "AUROC"]
        ].to_dict(orient="records"),
        "summary": str(summary),
        "html": str(html),
    }
    (OUT / "neoici_algorithm_gauntlet_summary.json").write_text(json.dumps(payload, indent=2))
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
