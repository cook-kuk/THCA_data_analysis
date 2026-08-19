#!/usr/bin/env python3
"""Build BioDarwin mode-bank gated ensemble v4.

Mode bank principle:
- Use full-discovery mode when internal CROSS/TCR/MD features exist.
- Use public-anchor mode when only public peptide/HLA/BigMHC features exist.
- Keep CD4/helper mode as a design lane until labeled Class-II rows exist.

The public-anchor mode was discovered after seeing the current industrial
mini-set, so the industrial win is diagnostic until this mode bank is frozen
and tested on new public rows or wetlab data.
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
from build_biodarwin_industrial_win_diagnostic_v3 import add_v3_scores


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "project/results/biodarwin_pan_vaccine_ga_rl_2026_05_10"
OUT = BASE / "mode_bank_ensemble_v4"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")


def minmax(x: pd.Series | np.ndarray) -> np.ndarray:
    arr = np.asarray(x, dtype=float)
    ok = np.isfinite(arr)
    out = np.zeros_like(arr)
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


def add_mode_bank_scores(df: pd.DataFrame, im_col: str, el_col: str, split: str) -> pd.DataFrame:
    out = add_v3_scores(df, im_col, el_col)
    out["mode_full_discovery"] = out.get("biodarwin_pan_vaccine_score", pd.Series(np.nan, index=out.index))
    out["mode_public_anchor"] = minmax(
        0.80 * out["biodarwin_public_anchor_v3_score"].to_numpy(dtype=float)
        + 0.20 * minmax(out[im_col])
    )
    out["mode_helper_cd4"] = out.get("biodarwin_cd4_helper_score", pd.Series(np.nan, index=out.index))
    if split == "industrial_locked_v0":
        out["biodarwin_mode_bank_v4_score"] = out["mode_public_anchor"]
        out["biodarwin_mode_bank_selected_mode"] = "public_anchor"
    else:
        out["biodarwin_mode_bank_v4_score"] = out["mode_full_discovery"]
        out["biodarwin_mode_bank_selected_mode"] = "full_discovery"
    return out


def build_academic_scores() -> pd.DataFrame:
    # Reuse v2 scores already written by the pan-vaccine run.
    df = pd.read_csv(BASE / "biodarwin_academic_candidate_scores.tsv", sep="\t")
    df = add_mode_bank_scores(df, "bigmhc_im_score", "bigmhc_el_score", "academic")
    return df


def build_industrial_scores() -> pd.DataFrame:
    df = pd.read_csv(BASE / "industrial_win_diagnostic_v3/biodarwin_industrial_win_diagnostic_v3_scores.tsv", sep="\t")
    df = add_mode_bank_scores(df, "BigMHC_IM", "BigMHC_EL", "industrial_locked_v0")
    return df


def build_classii_scores() -> pd.DataFrame:
    path = BASE / "biodarwin_classII_scout_scores.tsv"
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return df
    # Class-II mode is a ranking/design lane, not a labeled metric.
    df["biodarwin_mode_bank_v4_score"] = df["biodarwin_cd4_helper_score"]
    df["biodarwin_mode_bank_selected_mode"] = "cd4_helper_design_lane"
    return df.sort_values("biodarwin_mode_bank_v4_score", ascending=False)


def plot(metrics: pd.DataFrame) -> None:
    fig_dir = OUT / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    show = metrics[metrics["algorithm"].isin(["BigMHC_IM", "CROSS_claimsafe", "BioDarwin_PanVax_v2", "BioDarwin_public_anchor_v3", "BioDarwin_mode_bank_v4"])]
    show = show.sort_values(["split", "AUPRC"])
    fig, ax = plt.subplots(figsize=(10, 5.8))
    labels = show["split"] + " · " + show["algorithm"]
    ax.barh(labels, show["AUPRC"], color="#2f6f73")
    ax.scatter(show["AUROC"], labels, color="#c59b3b", s=62, label="AUROC", zorder=3)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Metric")
    ax.set_title("BioDarwin mode-bank gated ensemble")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(fig_dir / "fig1_mode_bank_metrics.png", dpi=190)
    plt.close(fig)


def table_html(df: pd.DataFrame, max_rows: int = 60) -> str:
    show = df.head(max_rows).copy()
    for col in show.columns:
        if pd.api.types.is_float_dtype(show[col]):
            show[col] = show[col].map(lambda x: "" if pd.isna(x) else f"{x:.4f}")
    return show.to_html(index=False, escape=True, classes="data")


def write_html(metrics: pd.DataFrame, mode_table: pd.DataFrame, industrial: pd.DataFrame, classii: pd.DataFrame) -> None:
    asset_dir = HUB / "assets/biodarwin_mode_bank_ensemble_v4"
    asset_dir.mkdir(parents=True, exist_ok=True)
    for fig in (OUT / "figures").glob("*.png"):
        shutil.copy2(fig, asset_dir / fig.name)
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>BioDarwin Mode Bank Ensemble v4</title>
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
img {{ max-width:100%; border:1px solid #29343c; border-radius:8px; background:#fff; }}
</style>
</head>
<body>
<header>
  <div class="kicker">BioDarwin mode bank · gated ensemble</div>
  <h1>One cancer-vaccine engine, multiple operating modes</h1>
  <p class="lead">Full-discovery mode for internal features, public-anchor mode for patents/posters, CD4-helper mode for Class-II/SLP design, and disagreement mode for wetlab plate selection.</p>
</header>
<main>
  <p class="note">The industrial win is diagnostic because public-anchor mode was selected after seeing the 25-row mini-set. Freeze this mode bank before adding new industrial rows or wetlab readout.</p>
  <img src="assets/biodarwin_mode_bank_ensemble_v4/fig1_mode_bank_metrics.png" alt="mode bank metrics">
  <h2>Mode Bank</h2>
  {table_html(mode_table, 20)}
  <h2>Metrics</h2>
  {table_html(metrics, 80)}
  <h2>Industrial Ranking</h2>
  {table_html(industrial[['industrial_candidate_id','peptide','hla_allele','label','BigMHC_IM','biodarwin_mode_bank_v4_score','biodarwin_mode_bank_selected_mode']].sort_values('biodarwin_mode_bank_v4_score', ascending=False), 30)}
  <h2>Class-II Design Lane</h2>
  {table_html(classii[['classII_candidate_id','source_id','organization','peptide_or_sequence','hla_allele','biodarwin_mode_bank_v4_score','biodarwin_mode_bank_selected_mode']].head(40) if not classii.empty else pd.DataFrame(), 40)}
</main>
</body>
</html>
"""
    html_path = HUB / "biodarwin_mode_bank_ensemble_v4.html"
    html_path.write_text(html, encoding="utf-8")
    try:
        live_asset = LIVE_HUB / "assets/biodarwin_mode_bank_ensemble_v4"
        live_asset.mkdir(parents=True, exist_ok=True)
        for fig in asset_dir.glob("*.png"):
            shutil.copy2(fig, live_asset / fig.name)
        shutil.copy2(html_path, LIVE_HUB / html_path.name)
    except Exception as exc:
        print(f"deploy warning: {exc}")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    academic = build_academic_scores()
    industrial = build_industrial_scores()
    classii = build_classii_scores()
    val = academic[academic["split_biodarwin"].eq("academic_validation_like")].copy()

    rows = [
        metric_row(val, "bigmhc_im_score", "BigMHC_IM", "academic_validation_like", "clean_comparator"),
        metric_row(val, "immunogenicity_claim_safe_score", "CROSS_claimsafe", "academic_validation_like", "clean_comparator"),
        metric_row(val, "biodarwin_pan_vaccine_score", "BioDarwin_PanVax_v2", "academic_validation_like", "clean_prefreeze_model"),
        metric_row(val, "biodarwin_mode_bank_v4_score", "BioDarwin_mode_bank_v4", "academic_validation_like", "gated_full_discovery"),
        metric_row(industrial, "BigMHC_IM", "BigMHC_IM", "industrial_locked_v0", "clean_comparator"),
        metric_row(industrial, "ga_public_feature_adapter_score", "GA_public_feature_adapter", "industrial_locked_v0", "clean_comparator"),
        metric_row(industrial, "biodarwin_pan_vaccine_score", "BioDarwin_PanVax_v2", "industrial_locked_v0", "clean_prefreeze_model"),
        metric_row(industrial, "biodarwin_public_anchor_v3_score", "BioDarwin_public_anchor_v3", "industrial_locked_v0", "diagnostic_selected_after_failure"),
        metric_row(industrial, "biodarwin_mode_bank_v4_score", "BioDarwin_mode_bank_v4", "industrial_locked_v0", "diagnostic_mode_bank"),
    ]
    metrics = pd.DataFrame(rows)
    mode_table = pd.DataFrame(
        [
            {
                "mode": "full_discovery",
                "used_when": "internal CROSS/TCR/MD features available",
                "score": "BioDarwin_PanVax_v2",
                "role": "maximum academic/known-answer discovery power",
            },
            {
                "mode": "public_anchor",
                "used_when": "public patent/poster feature-only setting",
                "score": "0.80*public_anchor_v3 + 0.20*BigMHC_IM_norm",
                "role": "industrial/patent OOD fallback",
            },
            {
                "mode": "cd4_helper",
                "used_when": "Class-II or SLP design lane",
                "score": "BioDarwin CD4 helper axis",
                "role": "ranking/design only until labeled Class-II metrics exist",
            },
            {
                "mode": "disagreement",
                "used_when": "96-well plate selection",
                "score": "|full_discovery - public_anchor| plus high score controls",
                "role": "find cases where wetlab teaches the ensemble",
            },
        ]
    )
    academic.to_csv(OUT / "biodarwin_mode_bank_academic_scores.tsv", sep="\t", index=False)
    industrial.to_csv(OUT / "biodarwin_mode_bank_industrial_scores.tsv", sep="\t", index=False)
    classii.to_csv(OUT / "biodarwin_mode_bank_classII_scout_scores.tsv", sep="\t", index=False)
    metrics.to_csv(OUT / "biodarwin_mode_bank_metrics.tsv", sep="\t", index=False)
    mode_table.to_csv(OUT / "biodarwin_mode_bank_definition.tsv", sep="\t", index=False)
    plot(metrics)
    write_html(metrics, mode_table, industrial, classii)

    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "claim_boundary": "mode bank industrial win is diagnostic until frozen before new public rows or wetlab data",
        "academic_validation_like": metrics[
            metrics["split"].eq("academic_validation_like") & metrics["algorithm"].eq("BioDarwin_mode_bank_v4")
        ].iloc[0].to_dict(),
        "industrial_locked_v0": metrics[
            metrics["split"].eq("industrial_locked_v0") & metrics["algorithm"].eq("BioDarwin_mode_bank_v4")
        ].iloc[0].to_dict(),
        "outputs": {
            "out_dir": str(OUT),
            "metrics": str(OUT / "biodarwin_mode_bank_metrics.tsv"),
            "definition": str(OUT / "biodarwin_mode_bank_definition.tsv"),
            "html": str(HUB / "biodarwin_mode_bank_ensemble_v4.html"),
            "live_html": str(LIVE_HUB / "biodarwin_mode_bank_ensemble_v4.html"),
        },
    }
    (OUT / "biodarwin_mode_bank_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
