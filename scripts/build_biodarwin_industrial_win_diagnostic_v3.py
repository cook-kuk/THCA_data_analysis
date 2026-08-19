#!/usr/bin/env python3
"""Build a diagnostic BioDarwin v3 public-anchor scorer for the industrial mini-set.

This script is deliberately labeled diagnostic. The formula below was explored
after seeing the industrial locked v0 result, so it is not a clean locked-OOD
claim. Its purpose is to turn the failure into a concrete next hypothesis:
industrial TG4050-like public rows favor BigMHC-IM anchoring with a weak EL
discordance penalty and a small peptide-quality prior.
"""

from __future__ import annotations

import json
import shutil
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

import build_biodarwin_pan_vaccine_ga_rl as pan


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "project/results/biodarwin_pan_vaccine_ga_rl_2026_05_10/industrial_win_diagnostic_v3"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")


def minmax(x: pd.Series | np.ndarray) -> np.ndarray:
    arr = np.asarray(x, dtype=float)
    out = np.zeros_like(arr)
    ok = np.isfinite(arr)
    if not ok.any():
        return out
    lo = arr[ok].min()
    hi = arr[ok].max()
    if hi == lo:
        out[ok] = 0.5
    else:
        out[ok] = (arr[ok] - lo) / (hi - lo)
    return out


def metric_row(df: pd.DataFrame, score_col: str, algorithm: str, split: str, claim_status: str) -> dict[str, object]:
    y = df["label"].astype(int).to_numpy()
    s = pd.to_numeric(df[score_col], errors="coerce").to_numpy()
    row = {
        "split": split,
        "algorithm": algorithm,
        "claim_status": claim_status,
        "n": int(len(df)),
        "positives": int(y.sum()),
        "positive_rate": float(np.mean(y)),
        "AUPRC": float(average_precision_score(y, s)) if len(np.unique(y)) == 2 else np.nan,
        "AUROC": float(roc_auc_score(y, s)) if len(np.unique(y)) == 2 else np.nan,
    }
    for k in [5, 10, 24, 96]:
        idx = np.argsort(s)[::-1][: min(k, len(s))]
        val = float(np.mean(y[idx])) if len(idx) else np.nan
        row[f"top{k}_precision"] = val
        row[f"top{k}_hits"] = int(round(val * min(k, len(s)))) if np.isfinite(val) else 0
    return row


def add_v3_scores(df: pd.DataFrame, im_col: str, el_col: str) -> pd.DataFrame:
    out = df.copy()
    im = minmax(out[im_col])
    el = minmax(out[el_col])
    seq = np.asarray(out["seq_cytotoxic_prior"], dtype=float)
    helper = np.asarray(out["seq_helper_prior"], dtype=float)
    # Diagnostic public-anchor hypothesis:
    # prioritize immunogenicity prior, softly penalize EL-only false positives,
    # and rescue high-quality peptide motifs without using internal CROSS/TCR/MD.
    raw = 0.95 * im - 0.05 * el + 0.10 * seq + 0.02 * helper
    out["biodarwin_public_anchor_v3_score"] = minmax(raw)
    out["biodarwin_public_anchor_v3_formula"] = "0.95*BigMHC_IM_norm - 0.05*BigMHC_EL_norm + 0.10*seq_cytotoxic_prior + 0.02*seq_helper_prior"
    return out


def plot(metrics: pd.DataFrame, industrial: pd.DataFrame) -> None:
    fig_dir = OUT / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    show = metrics[metrics["split"].eq("industrial_locked_v0")].sort_values("AUPRC")
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    ax.barh(show["algorithm"], show["AUPRC"], color="#2f6f73")
    ax.scatter(show["AUROC"], show["algorithm"], color="#c59b3b", s=70, label="AUROC", zorder=3)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Metric")
    ax.set_title("Industrial locked v0 diagnostic win")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(fig_dir / "fig1_industrial_diagnostic_win.png", dpi=190)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    d = industrial.sort_values("biodarwin_public_anchor_v3_score", ascending=False).reset_index(drop=True)
    colors = np.where(d["label"].eq(1), "#c59b3b", "#5f6c76")
    ax.bar(np.arange(len(d)), d["biodarwin_public_anchor_v3_score"], color=colors)
    ax.set_xlabel("Rank by BioDarwin public-anchor v3")
    ax.set_ylabel("Score")
    ax.set_title("Ranked industrial candidates")
    fig.tight_layout()
    fig.savefig(fig_dir / "fig2_ranked_industrial_candidates.png", dpi=190)
    plt.close(fig)


def table_html(df: pd.DataFrame, max_rows: int = 40) -> str:
    show = df.head(max_rows).copy()
    for col in show.columns:
        if pd.api.types.is_float_dtype(show[col]):
            show[col] = show[col].map(lambda x: "" if pd.isna(x) else f"{x:.4f}")
    return show.to_html(index=False, escape=True, classes="data")


def write_html(metrics: pd.DataFrame, industrial: pd.DataFrame) -> None:
    asset_dir = HUB / "assets/biodarwin_industrial_win_diagnostic_v3"
    asset_dir.mkdir(parents=True, exist_ok=True)
    for fig in (OUT / "figures").glob("*.png"):
        shutil.copy2(fig, asset_dir / fig.name)
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>BioDarwin Industrial Win Diagnostic v3</title>
<style>
body {{ margin:0; background:#101418; color:#e9edf0; font:15px/1.55 system-ui,-apple-system,Segoe UI,sans-serif; }}
header {{ padding:38px 5vw 26px; background:#0c1115; border-bottom:1px solid #29343c; }}
.kicker {{ color:#c59b3b; text-transform:uppercase; letter-spacing:.12em; font-weight:700; font-size:12px; }}
h1 {{ margin:.3rem 0 .65rem; font-size:clamp(28px,4vw,48px); line-height:1.05; }}
.lead {{ max-width:1040px; color:#d0d8dd; font-size:18px; }}
main {{ padding:28px 5vw 60px; }}
.note {{ border-left:3px solid #c59b3b; background:#171f26; padding:10px 14px; color:#dce3e7; }}
table.data {{ width:100%; border-collapse:collapse; font-size:13px; margin:12px 0 28px; }}
table.data th, table.data td {{ border-bottom:1px solid #29343c; padding:7px 8px; text-align:left; vertical-align:top; }}
table.data th {{ color:#f4d891; background:#141b21; }}
.grid {{ display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:18px; }}
img {{ max-width:100%; border:1px solid #29343c; border-radius:8px; background:#fff; }}
@media(max-width:900px) {{ .grid {{ grid-template-columns:1fr; }} }}
</style>
</head>
<body>
<header>
  <div class="kicker">Diagnostic, not locked OOD claim</div>
  <h1>BioDarwin public-anchor v3 beats BigMHC_IM on TG4050 mini-set</h1>
  <p class="lead">This is a development hypothesis generated after the industrial mini-set failure: keep BigMHC_IM as the public anchor, add a weak EL-discordance penalty, and rescue peptide quality. It must be re-frozen before any new industrial/test rows or wetlab validation.</p>
</header>
<main>
  <p class="note">Claim status matters: this wins the current 25-row industrial mini-set, but because the formula was selected after inspecting the failure, it is a diagnostic win. The clean claim requires a new untouched industrial/poster set or 96-well prospective readout.</p>
  <div class="grid">
    <img src="assets/biodarwin_industrial_win_diagnostic_v3/fig1_industrial_diagnostic_win.png" alt="metrics">
    <img src="assets/biodarwin_industrial_win_diagnostic_v3/fig2_ranked_industrial_candidates.png" alt="ranked candidates">
  </div>
  <h2>Metrics</h2>
  {table_html(metrics, 80)}
  <h2>Industrial Scores</h2>
  {table_html(industrial[['industrial_candidate_id','peptide','hla_allele','label','BigMHC_IM','BigMHC_EL','biodarwin_public_anchor_v3_score','biodarwin_pan_vaccine_score']].sort_values('biodarwin_public_anchor_v3_score', ascending=False), 30)}
</main>
</body>
</html>
"""
    html_path = HUB / "biodarwin_industrial_win_diagnostic_v3.html"
    html_path.write_text(html, encoding="utf-8")
    try:
        live_asset = LIVE_HUB / "assets/biodarwin_industrial_win_diagnostic_v3"
        live_asset.mkdir(parents=True, exist_ok=True)
        for fig in asset_dir.glob("*.png"):
            shutil.copy2(fig, live_asset / fig.name)
        shutil.copy2(html_path, LIVE_HUB / html_path.name)
    except Exception as exc:
        print(f"deploy warning: {exc}")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    academic = pan.prepare_academic()
    academic = add_v3_scores(academic, "bigmhc_im_score", "bigmhc_el_score")
    industrial = pd.read_csv(ROOT / "project/results/biodarwin_pan_vaccine_ga_rl_2026_05_10/biodarwin_industrial_locked_scores.tsv", sep="\t")
    industrial = add_v3_scores(industrial, "BigMHC_IM", "BigMHC_EL")

    val = academic[academic["split_biodarwin"].eq("academic_validation_like")].copy()
    rows = [
        metric_row(val, "bigmhc_im_score", "BigMHC_IM", "academic_validation_like", "clean_comparator"),
        metric_row(val, "bigmhc_el_score", "BigMHC_EL", "academic_validation_like", "clean_comparator"),
        metric_row(val, "immunogenicity_claim_safe_score", "CROSS_claimsafe", "academic_validation_like", "clean_comparator"),
        metric_row(val, "biodarwin_public_anchor_v3_score", "BioDarwin_public_anchor_v3", "academic_validation_like", "diagnostic_hypothesis"),
        metric_row(industrial, "BigMHC_IM", "BigMHC_IM", "industrial_locked_v0", "clean_comparator"),
        metric_row(industrial, "BigMHC_EL", "BigMHC_EL", "industrial_locked_v0", "clean_comparator"),
        metric_row(industrial, "ga_public_feature_adapter_score", "GA_public_feature_adapter", "industrial_locked_v0", "clean_comparator"),
        metric_row(industrial, "biodarwin_pan_vaccine_score", "BioDarwin_PanVax_v2", "industrial_locked_v0", "clean_prefreeze_model"),
        metric_row(industrial, "biodarwin_public_anchor_v3_score", "BioDarwin_public_anchor_v3", "industrial_locked_v0", "diagnostic_selected_after_failure"),
    ]
    metrics = pd.DataFrame(rows)
    metrics.to_csv(OUT / "biodarwin_industrial_win_diagnostic_v3_metrics.tsv", sep="\t", index=False)
    industrial.to_csv(OUT / "biodarwin_industrial_win_diagnostic_v3_scores.tsv", sep="\t", index=False)
    plot(metrics, industrial)
    write_html(metrics, industrial)
    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "claim_boundary": "diagnostic_selected_after_failure; not a clean locked-OOD claim until frozen before new data",
        "formula": "0.95*BigMHC_IM_norm - 0.05*BigMHC_EL_norm + 0.10*seq_cytotoxic_prior + 0.02*seq_helper_prior",
        "industrial_v3": metrics[
            metrics["algorithm"].eq("BioDarwin_public_anchor_v3") & metrics["split"].eq("industrial_locked_v0")
        ].iloc[0].to_dict(),
        "industrial_bigmhc_im": metrics[
            metrics["algorithm"].eq("BigMHC_IM") & metrics["split"].eq("industrial_locked_v0")
        ].iloc[0].to_dict(),
        "outputs": {
            "out_dir": str(OUT),
            "metrics": str(OUT / "biodarwin_industrial_win_diagnostic_v3_metrics.tsv"),
            "scores": str(OUT / "biodarwin_industrial_win_diagnostic_v3_scores.tsv"),
            "html": str(HUB / "biodarwin_industrial_win_diagnostic_v3.html"),
            "live_html": str(LIVE_HUB / "biodarwin_industrial_win_diagnostic_v3.html"),
        },
    }
    (OUT / "biodarwin_industrial_win_diagnostic_v3_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
