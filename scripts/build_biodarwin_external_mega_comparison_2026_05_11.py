#!/usr/bin/env python3
"""Build a unified external benchmark comparison for BioDarwin.

This script consolidates:
- local public external benchmark bundles from p_neo_bayesian wave11
- BioDarwin academic/industrial mode-bank results
- CleanNeo/BarNeo internal leaderboard summaries
- a current competitor landscape table from primary/public sources

The output is a comparison scaffold. Rows are explicitly flagged for claim
status, low power, overlap, and diagnostic-after-failure status.
"""

from __future__ import annotations

import json
import re
import shutil
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

import build_biodarwin_pan_vaccine_ga_rl as pan


ROOT = Path(__file__).resolve().parents[1]
WAVE11 = ROOT / "project/results/p_neo_bayesian_2026_05_09/wave11"
PRED_DIR = WAVE11 / "wave11_predictions"
BIODARWIN = ROOT / "project/results/biodarwin_pan_vaccine_ga_rl_2026_05_10"
MODE_BANK = BIODARWIN / "mode_bank_ensemble_v4"
CLEAN_NEO = ROOT / "project/results/clean_neobench_barneo_2026_05_09"
OUT = ROOT / "project/results/biodarwin_external_mega_comparison_2026_05_11"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")


def minmax(x: pd.Series | np.ndarray) -> np.ndarray:
    arr = np.asarray(x, dtype=float)
    ok = np.isfinite(arr)
    out = np.zeros_like(arr, dtype=float)
    if not ok.any():
        return out
    lo = arr[ok].min()
    hi = arr[ok].max()
    if hi == lo:
        out[ok] = 0.5
    else:
        out[ok] = (arr[ok] - lo) / (hi - lo)
    return out


def metric_dict(y: np.ndarray, s: np.ndarray) -> dict[str, float | int | bool]:
    y = np.asarray(y, dtype=int)
    s = np.asarray(s, dtype=float)
    ok = np.isfinite(s) & np.isfinite(y)
    y = y[ok]
    s = s[ok]
    n = int(len(y))
    n_pos = int(y.sum()) if n else 0
    n_neg = int(n - n_pos)
    row: dict[str, float | int | bool] = {
        "n": n,
        "n_pos": n_pos,
        "n_neg": n_neg,
        "positive_rate": float(np.mean(y)) if n else np.nan,
        "binary_metric_possible": bool(n_pos > 0 and n_neg > 0),
        "low_power": bool(n_pos < 10 or n_neg < 10),
    }
    if n_pos > 0 and n_neg > 0:
        row["AUPRC"] = float(average_precision_score(y, s))
        row["AUROC"] = float(roc_auc_score(y, s))
    else:
        row["AUPRC"] = np.nan
        row["AUROC"] = np.nan
    for k in [5, 10, 20, 50, 96]:
        if n:
            idx = np.argsort(s)[::-1][: min(k, n)]
            val = float(np.mean(y[idx])) if len(idx) else np.nan
        else:
            val = np.nan
        row[f"top{k}_precision"] = val
        row[f"top{k}_hits"] = int(round(val * min(k, n))) if np.isfinite(val) else 0
    return row


def score_column(df: pd.DataFrame) -> str | None:
    candidates = [c for c in df.columns if c.startswith("score_") and not c.endswith("_rank")]
    if not candidates:
        return None
    preferred = [c for c in candidates if c not in {"score_mhcflurry_aff"}]
    return preferred[0] if preferred else candidates[0]


def read_overlap_audit() -> pd.DataFrame:
    path = WAVE11 / "wave11_overlap_audit.tsv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def wave11_prediction_metrics() -> pd.DataFrame:
    rows = []
    overlap = read_overlap_audit()
    overlap_key = {}
    if not overlap.empty:
        for _, r in overlap.iterrows():
            overlap_key[(str(r["algorithm"]), str(r["test_bundle"]))] = float(r.get("frac_in_training", np.nan))
    for path in sorted(PRED_DIR.glob("*.tsv")):
        m = re.match(r"(.+)__(.+)\.tsv$", path.name)
        if not m:
            continue
        algorithm, bundle = m.group(1), m.group(2)
        df = pd.read_csv(path, sep="\t")
        col = score_column(df)
        if col is None or "label" not in df:
            continue
        row = {
            "benchmark_family": "wave11_public_external",
            "dataset": bundle,
            "algorithm": algorithm,
            "score_col": col,
            "claim_status": "local_prediction_file",
            "locally_evaluated": True,
            "diagnostic_after_failure": False,
            "notes": "",
        }
        row.update(metric_dict(df["label"].to_numpy(), df[col].to_numpy()))
        frac = overlap_key.get((algorithm, bundle), np.nan)
        row["training_overlap_fraction"] = frac
        row["overlap_warning"] = bool(np.isfinite(frac) and frac >= 0.20)
        rows.append(row)
    return pd.DataFrame(rows)


def peptide_public_features(peptides: pd.Series) -> pd.DataFrame:
    return pd.DataFrame([pan.peptide_features(x) for x in peptides], index=peptides.index)


def biodarwin_anchor_wave11_metrics() -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    scores = []
    for path in sorted(PRED_DIR.glob("BigMHC_IM__*.tsv")):
        bundle = path.name.split("__", 1)[1].replace(".tsv", "")
        df = pd.read_csv(path, sep="\t")
        if "score_bigmhc_im" not in df or "label" not in df:
            continue
        pf = peptide_public_features(df["peptide"])
        im = minmax(df["score_bigmhc_im"])
        # No BigMHC_EL is available across wave11, so this is the no-EL public
        # anchor variant. It is diagnostic until frozen against new data.
        raw = 0.95 * im + 0.10 * pf["seq_cytotoxic_prior"].to_numpy(float) + 0.02 * pf["seq_helper_prior"].to_numpy(float)
        df["biodarwin_public_anchor_noEL_v4_score"] = minmax(raw)
        df["test_bundle"] = bundle
        scores.append(df)
        row = {
            "benchmark_family": "wave11_public_external",
            "dataset": bundle,
            "algorithm": "BioDarwin_public_anchor_noEL_v4",
            "score_col": "biodarwin_public_anchor_noEL_v4_score",
            "claim_status": "diagnostic_public_anchor_noEL",
            "locally_evaluated": True,
            "diagnostic_after_failure": True,
            "training_overlap_fraction": np.nan,
            "overlap_warning": False,
            "notes": "Derived from mode-bank public-anchor concept; no industrial labels used here, but mode selected after TG4050 failure.",
        }
        row.update(metric_dict(df["label"].to_numpy(), df["biodarwin_public_anchor_noEL_v4_score"].to_numpy()))
        rows.append(row)
    score_df = pd.concat(scores, ignore_index=True) if scores else pd.DataFrame()
    return pd.DataFrame(rows), score_df


def biodarwin_mode_bank_metrics() -> pd.DataFrame:
    path = MODE_BANK / "biodarwin_mode_bank_metrics.tsv"
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t")
    out = df.rename(columns={"split": "dataset"}).copy()
    out.insert(0, "benchmark_family", "biodarwin_mode_bank")
    out["score_col"] = "reported"
    out["locally_evaluated"] = True
    out["diagnostic_after_failure"] = out["claim_status"].astype(str).str.contains("diagnostic", case=False)
    out["training_overlap_fraction"] = np.nan
    out["overlap_warning"] = False
    out["binary_metric_possible"] = True
    out["low_power"] = out["positives"].lt(10) | (out["n"] - out["positives"]).lt(10)
    out = out.rename(columns={"positives": "n_pos"})
    out["n_neg"] = out["n"] - out["n_pos"]
    out["notes"] = "Mode-bank v4; industrial win diagnostic until frozen before new data."
    return out[
        [
            "benchmark_family",
            "dataset",
            "algorithm",
            "score_col",
            "claim_status",
            "locally_evaluated",
            "diagnostic_after_failure",
            "training_overlap_fraction",
            "overlap_warning",
            "binary_metric_possible",
            "low_power",
            "n",
            "n_pos",
            "n_neg",
            "positive_rate",
            "AUPRC",
            "AUROC",
            "top5_precision",
            "top5_hits",
            "top10_precision",
            "top10_hits",
            "top24_precision",
            "top24_hits",
            "top96_precision",
            "top96_hits",
            "notes",
        ]
    ]


def clean_neobench_summary() -> pd.DataFrame:
    path = CLEAN_NEO / "clean_neobench_leaderboard.tsv"
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t").head(40).copy()
    df.insert(0, "benchmark_family", "clean_neobench_leaderboard")
    return df


def source_inventory() -> pd.DataFrame:
    rows = []
    manifest = WAVE11 / "test_bundle_manifest.tsv"
    if manifest.exists():
        df = pd.read_csv(manifest, sep="\t")
        for _, r in df.iterrows():
            rows.append(
                {
                    "benchmark_family": "wave11_public_external",
                    "dataset": r["bundle"],
                    "source": r["source"],
                    "validation_type": r["validation"],
                    "n": int(r["n"]),
                    "n_pos": int(r["n_pos"]),
                    "n_with_hla": int(r["n_with_hla"]),
                    "n_8_15": int(r["n_8_15"]),
                    "status": r["status"],
                }
            )
    industrial_summary = ROOT / "project/results/public_industrial_neoantigen_benchmark_scout_2026_05_10/public_industrial_benchmark_scout_summary.json"
    if industrial_summary.exists():
        data = json.loads(industrial_summary.read_text(encoding="utf-8"))
        rows.append(
            {
                "benchmark_family": "public_industrial_scout",
                "dataset": "industrial_patent_company_poster_scout",
                "source": "59 public sources / 19 org families",
                "validation_type": "LOCKED_PUBLIC_INDUSTRIAL_SCOUT",
                "n": int(data["n_candidate_pairs"]),
                "n_pos": np.nan,
                "n_with_hla": int(data["n_direct_classI_score_ready"] + data["n_direct_classII_score_ready"]),
                "n_8_15": int(data["n_direct_classI_score_ready"]),
                "status": f"locked_metric_classI={data['n_locked_classI_testset_v0']}; locked_metric_classII={data['n_locked_classII_testset_v0']}",
            }
        )
    return pd.DataFrame(rows)


def latest_competitor_landscape() -> pd.DataFrame:
    rows = [
        {
            "method": "BigMHC",
            "year": 2023,
            "type": "MHC-I presentation + immunogenicity transfer learning",
            "local_status": "evaluated_locally",
            "local_score_files": "wave11 BigMHC_IM; BioDarwin industrial comparator",
            "primary_source": "https://www.nature.com/articles/s42256-023-00694-6",
            "code_or_server": "https://github.com/KarchinLab/bigmhc",
            "notes": "Strong public comparator; industrial TG4050 mini-set baseline to beat.",
        },
        {
            "method": "NetMHCpan-4.1 / NetMHCIIpan",
            "year": 2020,
            "type": "MHC-I/MHC-II binding and eluted ligand prediction",
            "local_status": "partial_local_ITSNdb_only",
            "local_score_files": "wave11 NetMHCpan_4.1__itsndb; prior Class-II runs",
            "primary_source": "https://academic.oup.com/nar/article/48/W1/W449/5837056",
            "code_or_server": "https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/",
            "notes": "Class-II server/predictor is important for next CD4-helper benchmark.",
        },
        {
            "method": "MHCflurry 2.0",
            "year": 2020,
            "type": "MHC-I binding, antigen processing, presentation",
            "local_status": "evaluated_locally",
            "local_score_files": "wave11 MHCflurry across CEDAR/dbPepNeo2/ITSNdb/McPAS/NeoDB/NEPdb/TESLA",
            "primary_source": "https://doi.org/10.1016/j.cels.2020.06.010",
            "code_or_server": "https://openvax.github.io/mhcflurry/",
            "notes": "Best clean external AUROC on several wave11 bundles.",
        },
        {
            "method": "PRIME / PRIME2.0",
            "year": 2021,
            "type": "presentation + peptide-intrinsic TCR recognition propensity",
            "local_status": "evaluated_locally",
            "local_score_files": "wave11 PRIME across multiple bundles",
            "primary_source": "https://www.sciencedirect.com/science/article/pii/S2666379121000057",
            "code_or_server": "http://prime.gfellerlab.org/",
            "notes": "Useful immunogenicity comparator; local rows available.",
        },
        {
            "method": "TransPHLA",
            "year": 2022,
            "type": "Transformer peptide-HLA binding predictor",
            "local_status": "evaluated_locally",
            "local_score_files": "wave11 TransPHLA across multiple bundles",
            "primary_source": "https://www.nature.com/articles/s42256-022-00459-7",
            "code_or_server": "https://github.com/a96123155/TransPHLA-AOMP",
            "notes": "Binding-oriented; some local scores saturate near 1.0 and require calibration caution.",
        },
        {
            "method": "DeepImmuno",
            "year": 2021,
            "type": "deep immunogenicity predictor",
            "local_status": "ITSNdb_local",
            "local_score_files": "wave11 DeepImmuno__itsndb",
            "primary_source": "https://github.com/frankligy/DeepImmuno",
            "code_or_server": "https://github.com/frankligy/DeepImmuno",
            "notes": "Included in local and literature comparisons.",
        },
        {
            "method": "IMPROVE",
            "year": 2024,
            "type": "feature model for neoepitope immunogenicity validation",
            "local_status": "dataset_manifest_only",
            "local_score_files": "wave11 improve_cedar bundle has no positives in current manifest",
            "primary_source": "https://pubmed.ncbi.nlm.nih.gov/38633261/",
            "code_or_server": "",
            "notes": "Useful benchmark/source; current local extract needs label reconstruction before metrics.",
        },
        {
            "method": "MUNIS",
            "year": 2024,
            "type": "deep learning HLA-I presented CD8 epitope model",
            "local_status": "not_yet_local",
            "local_score_files": "",
            "primary_source": "https://www.nature.com/articles/s42256-024-00971-y",
            "code_or_server": "",
            "notes": "Presentation/epitope competitor to add if code/model access is available.",
        },
        {
            "method": "NeoTImmuML",
            "year": 2025,
            "type": "weighted ensemble immunogenicity ML",
            "local_status": "not_yet_local",
            "local_score_files": "",
            "primary_source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12585993/",
            "code_or_server": "",
            "notes": "Recent 2025 immunogenicity model; compare as literature/implementation target.",
        },
        {
            "method": "PHLA",
            "year": 2025,
            "type": "pretrained model incorporating HLA-peptide binding and immunogenicity",
            "local_status": "not_yet_local",
            "local_score_files": "",
            "primary_source": "https://www.sciencedirect.com/science/article/abs/pii/S0925231225019642",
            "code_or_server": "",
            "notes": "Recent 2025 method; needs install/source availability check.",
        },
        {
            "method": "weakly supervised peptide-TCR model",
            "year": 2025,
            "type": "patient/TCR-aware neoantigen identification",
            "local_status": "not_yet_local",
            "local_score_files": "",
            "primary_source": "https://www.sciencedirect.com/science/article/pii/S2405471225002364",
            "code_or_server": "",
            "notes": "TCR-profile mode is relevant to future BioDarwin patient-specific mode.",
        },
        {
            "method": "NeoPrecis-Immuno",
            "year": 2026,
            "type": "qualified immunogenicity + clonality-aware neoantigen landscape",
            "local_status": "not_yet_local",
            "local_score_files": "",
            "primary_source": "https://www.nature.com/articles/s41467-026-68651-6",
            "code_or_server": "",
            "notes": "Current high-level comparator concept; may be response-prediction rather than peptide-only benchmark.",
        },
    ]
    return pd.DataFrame(rows)


def algorithm_comparison_table(all_metrics: pd.DataFrame, landscape: pd.DataFrame) -> pd.DataFrame:
    """Create a compact 11-row table: our BioDarwin row plus 10 comparator algorithms."""

    def best_row_for_algorithm(alg: str) -> dict:
        sub = all_metrics[
            all_metrics["algorithm"].astype(str).str.fullmatch(re.escape(alg), case=False, na=False)
        ].copy()
        if sub.empty:
            sub = all_metrics[all_metrics["algorithm"].astype(str).str.contains(re.escape(alg), case=False, na=False)].copy()
        sub = sub[sub["AUPRC"].notna()]
        if sub.empty:
            return {
                "algorithm": alg,
                "dataset": "",
                "benchmark_family": "",
                "claim_status": "not_locally_evaluated",
                "AUPRC": np.nan,
                "AUROC": np.nan,
            }
        row = sub.sort_values(["AUPRC", "AUROC"], ascending=False).iloc[0].to_dict()
        return row

    ours_alg = "BioDarwin_mode_bank_v4"
    ours_row = best_row_for_algorithm(ours_alg)
    rows = [
        {
            "rank": 1,
            "algorithm": f"ours · {ours_alg}",
            "source_type": str(ours_row.get("benchmark_family", "biodarwin")),
            "dataset_or_bundle": str(ours_row.get("dataset", "")),
            "status": str(ours_row.get("claim_status", "")),
            "AUPRC": ours_row.get("AUPRC", np.nan),
            "AUROC": ours_row.get("AUROC", np.nan),
            "note": "BioDarwin representative row; label intentionally marked ours.",
        }
    ]

    comparator_order = [
        "RF_biophys",
        "MHCflurry",
        "PRIME",
        "TransPHLA",
        "BigMHC_IM",
        "ESM2_Bayesian",
        "Structure_LR",
        "DeepImmuno",
        "GP_quantum",
        "GA_public_feature_adapter",
    ]
    for rank, alg in enumerate(comparator_order, start=2):
        row = best_row_for_algorithm(alg)
        rows.append(
            {
                "rank": rank,
                "algorithm": alg,
                "source_type": str(row.get("benchmark_family", "")),
                "dataset_or_bundle": str(row.get("dataset", "")),
                "status": str(row.get("claim_status", "")),
                "AUPRC": row.get("AUPRC", np.nan),
                "AUROC": row.get("AUROC", np.nan),
                "note": "Best local row for this algorithm.",
            }
        )

    return pd.DataFrame(rows)[
        ["rank", "algorithm", "source_type", "dataset_or_bundle", "status", "AUPRC", "AUROC", "note"]
    ]


def plot_outputs(all_metrics: pd.DataFrame, inventory: pd.DataFrame) -> None:
    fig_dir = OUT / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    metricable = all_metrics[all_metrics["binary_metric_possible"].eq(True)].copy()
    metricable = metricable[metricable["AUPRC"].notna()]
    best = (
        metricable.sort_values("AUPRC", ascending=False)
        .groupby(["benchmark_family", "dataset"], as_index=False)
        .head(1)
        .sort_values("AUPRC")
    )
    fig, ax = plt.subplots(figsize=(10.5, max(4.8, 0.34 * len(best))))
    labels = best["benchmark_family"] + " · " + best["dataset"] + " · " + best["algorithm"]
    colors = np.where(best["algorithm"].astype(str).str.contains("BioDarwin"), "#c59b3b", "#2f6f73")
    ax.barh(labels, best["AUPRC"], color=colors)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Best AUPRC in available local/diagnostic rows")
    ax.set_title("Best local algorithm per external dataset")
    fig.tight_layout()
    fig.savefig(fig_dir / "fig1_best_per_dataset_auprc.png", dpi=190)
    plt.close(fig)

    wave = metricable[metricable["benchmark_family"].eq("wave11_public_external")].copy()
    if not wave.empty:
        pivot = wave.pivot_table(index="algorithm", columns="dataset", values="AUROC", aggfunc="max")
        fig, ax = plt.subplots(figsize=(11, max(5.5, 0.28 * len(pivot))))
        im = ax.imshow(pivot.fillna(np.nan).to_numpy(), aspect="auto", cmap="viridis", vmin=0.45, vmax=0.9)
        ax.set_xticks(range(len(pivot.columns)))
        ax.set_xticklabels(pivot.columns, rotation=35, ha="right")
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels(pivot.index)
        ax.set_title("Wave11 local AUROC matrix")
        fig.colorbar(im, ax=ax, label="AUROC")
        fig.tight_layout()
        fig.savefig(fig_dir / "fig2_wave11_auroc_matrix.png", dpi=190)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    inv = inventory.copy()
    inv["metricable_binary"] = (inv["n_pos"].fillna(0) > 0) & ((inv["n"] - inv["n_pos"].fillna(0)) > 0)
    counts = inv.groupby("benchmark_family")["dataset"].count().sort_values()
    ax.barh(counts.index, counts.values, color="#2f6f73")
    ax.set_xlabel("Datasets/source groups")
    ax.set_title("External benchmark inventory")
    fig.tight_layout()
    fig.savefig(fig_dir / "fig3_inventory_counts.png", dpi=190)
    plt.close(fig)


def table_html(df: pd.DataFrame, max_rows: int = 60) -> str:
    if df.empty:
        return "<p>No rows.</p>"
    show = df.head(max_rows).copy()
    for col in show.columns:
        if pd.api.types.is_float_dtype(show[col]):
            show[col] = show[col].map(lambda x: "" if pd.isna(x) else f"{x:.4f}")
    return show.to_html(index=False, escape=True, classes="data")


def write_report_and_html(all_metrics: pd.DataFrame, inventory: pd.DataFrame, clean: pd.DataFrame, landscape: pd.DataFrame) -> None:
    metricable = all_metrics[all_metrics["binary_metric_possible"].eq(True) & all_metrics["AUPRC"].notna()].copy()
    best = metricable.sort_values("AUPRC", ascending=False).groupby(["benchmark_family", "dataset"], as_index=False).head(1)
    safe = metricable[
        ~metricable["diagnostic_after_failure"].astype(bool)
        & ~metricable["overlap_warning"].astype(bool)
        & ~metricable["low_power"].astype(bool)
    ].copy()
    safe_best = safe.sort_values("AUPRC", ascending=False).groupby(["benchmark_family", "dataset"], as_index=False).head(1)
    algo_table = algorithm_comparison_table(all_metrics, landscape)

    report = [
        "# BioDarwin external mega comparison 2026-05-11",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Bottom line",
        "",
        f"- External/source inventory rows: {len(inventory)}.",
        f"- Local metric rows: {len(all_metrics)}.",
        f"- Binary metric-capable rows: {len(metricable)}.",
        f"- Reviewer-safer non-diagnostic/non-overlap/non-low-power rows: {len(safe)}.",
        "- BioDarwin mode-bank v4 wins the current industrial 25-row diagnostic mini-set, but that win is not a locked claim until frozen before new data.",
        "- Wave11 still shows MHCflurry/BigMHC/PRIME as strong public comparators depending on dataset; several bundles are all-positive and cannot support AUROC/AUPRC.",
        "",
        "## 11-algorithm comparison table",
        "",
        algo_table.to_markdown(index=False),
        "",
        "## Best per dataset",
        "",
        best.to_markdown(index=False),
        "",
        "## Reviewer-safer best per dataset",
        "",
        safe_best.to_markdown(index=False) if not safe_best.empty else "No reviewer-safer rows.",
        "",
    ]
    (OUT / "BIODARWIN_EXTERNAL_MEGA_COMPARISON_REPORT.md").write_text("\n".join(report), encoding="utf-8")

    asset_dir = HUB / "assets/biodarwin_external_mega_comparison_2026_05_11"
    asset_dir.mkdir(parents=True, exist_ok=True)
    for fig in (OUT / "figures").glob("*.png"):
        shutil.copy2(fig, asset_dir / fig.name)

    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>BioDarwin External Mega Comparison</title>
<style>
body {{ margin:0; background:#101418; color:#e9edf0; font:15px/1.55 system-ui,-apple-system,Segoe UI,sans-serif; }}
header {{ padding:38px 5vw 26px; background:#0c1115; border-bottom:1px solid #29343c; }}
.kicker {{ color:#c59b3b; text-transform:uppercase; letter-spacing:.12em; font-weight:700; font-size:12px; }}
h1 {{ margin:.3rem 0 .65rem; font-size:clamp(28px,4vw,50px); line-height:1.05; }}
.lead {{ max-width:1100px; color:#d0d8dd; font-size:18px; }}
.stats {{ display:grid; grid-template-columns:repeat(5,minmax(130px,1fr)); gap:10px; margin-top:20px; }}
.stat {{ border:1px solid #29343c; background:#151c22; border-radius:8px; padding:12px 14px; }}
.stat b {{ display:block; font-size:23px; }}
.stat span {{ color:#9aa7af; font-size:12px; }}
main {{ padding:28px 5vw 60px; }}
.note {{ border-left:3px solid #c59b3b; background:#171f26; padding:10px 14px; color:#dce3e7; }}
table.data {{ width:100%; border-collapse:collapse; font-size:13px; margin:12px 0 28px; }}
table.data th, table.data td {{ border-bottom:1px solid #29343c; padding:7px 8px; text-align:left; vertical-align:top; }}
table.data th {{ color:#f4d891; background:#141b21; }}
.grid {{ display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:18px; }}
img {{ max-width:100%; border:1px solid #29343c; border-radius:8px; background:#fff; }}
@media(max-width:900px) {{ .stats {{ grid-template-columns:repeat(2,1fr); }} .grid {{ grid-template-columns:1fr; }} }}
</style>
</head>
<body>
<header>
  <div class="kicker">BioDarwin · external mega comparison</div>
  <h1>Public external datasets, industrial scout and competitor algorithms</h1>
  <p class="lead">Unified comparison across local wave11 public external bundles, BioDarwin mode-bank results, CleanNeo/BarNeo leaderboards and current public competitor landscape. Claim flags are kept beside the metrics.</p>
  <div class="stats">
    <div class="stat"><b>{len(inventory)}</b><span>dataset/source inventory rows</span></div>
    <div class="stat"><b>{len(all_metrics)}</b><span>local metric rows</span></div>
    <div class="stat"><b>{len(metricable)}</b><span>binary metric-capable rows</span></div>
    <div class="stat"><b>{len(landscape)}</b><span>competitor methods tracked</span></div>
    <div class="stat"><b>{int(all_metrics['algorithm'].astype(str).str.contains('BioDarwin').sum())}</b><span>BioDarwin metric rows</span></div>
  </div>
</header>
<main>
  <h2>11-Algorithm Comparison</h2>
  <p class="note">AUPRC/AUROC are the best locally evaluated row per algorithm, so every comparator has a numeric value where one exists.</p>
  {table_html(algo_table, 20)}
  <p class="note">Do not collapse these into one marketing number. Some rows are diagnostic, low-power, all-positive, or overlap-flagged. The strongest manuscript claim must use the reviewer-safer subset or new frozen prospective rows.</p>
  <div class="grid">
    <img src="assets/biodarwin_external_mega_comparison_2026_05_11/fig1_best_per_dataset_auprc.png" alt="best per dataset">
    <img src="assets/biodarwin_external_mega_comparison_2026_05_11/fig2_wave11_auroc_matrix.png" alt="wave11 matrix">
  </div>
  <img src="assets/biodarwin_external_mega_comparison_2026_05_11/fig3_inventory_counts.png" alt="inventory">
  <h2>Best Per Dataset</h2>
  {table_html(best.sort_values(['benchmark_family','dataset']), 80)}
  <h2>Reviewer-Safer Best Per Dataset</h2>
  {table_html(safe_best.sort_values(['benchmark_family','dataset']), 80)}
  <h2>All Local Metrics</h2>
  {table_html(all_metrics.sort_values(['benchmark_family','dataset','AUPRC'], ascending=[True, True, False]), 160)}
  <h2>Source Inventory</h2>
  {table_html(inventory, 80)}
  <h2>Current Competitor Landscape</h2>
  {table_html(landscape, 80)}
  <h2>CleanNeo/BarNeo Top Local Leaderboard</h2>
  {table_html(clean, 50)}
</main>
</body>
</html>
"""
    html_path = HUB / "biodarwin_external_mega_comparison_2026_05_11.html"
    html_path.write_text(html, encoding="utf-8")
    try:
        live_asset = LIVE_HUB / "assets/biodarwin_external_mega_comparison_2026_05_11"
        live_asset.mkdir(parents=True, exist_ok=True)
        for fig in asset_dir.glob("*.png"):
            shutil.copy2(fig, live_asset / fig.name)
        shutil.copy2(html_path, LIVE_HUB / html_path.name)
    except Exception as exc:
        print(f"deploy warning: {exc}")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    wave = wave11_prediction_metrics()
    bio_wave, bio_wave_scores = biodarwin_anchor_wave11_metrics()
    mode = biodarwin_mode_bank_metrics()
    # Normalize columns across wave/mode metric tables.
    all_metrics = pd.concat([wave, bio_wave, mode], ignore_index=True, sort=False)
    clean = clean_neobench_summary()
    inventory = source_inventory()
    landscape = latest_competitor_landscape()

    all_metrics.to_csv(OUT / "biodarwin_external_all_local_metrics.tsv", sep="\t", index=False)
    bio_wave_scores.to_csv(OUT / "biodarwin_wave11_public_anchor_scores.tsv", sep="\t", index=False)
    clean.to_csv(OUT / "biodarwin_clean_neobench_top40.tsv", sep="\t", index=False)
    inventory.to_csv(OUT / "biodarwin_external_dataset_inventory.tsv", sep="\t", index=False)
    landscape.to_csv(OUT / "biodarwin_latest_competitor_landscape.tsv", sep="\t", index=False)

    metricable = all_metrics[all_metrics["binary_metric_possible"].eq(True) & all_metrics["AUPRC"].notna()].copy()
    best = metricable.sort_values("AUPRC", ascending=False).groupby(["benchmark_family", "dataset"], as_index=False).head(1)
    best.to_csv(OUT / "biodarwin_best_per_dataset.tsv", sep="\t", index=False)
    safe = metricable[
        ~metricable["diagnostic_after_failure"].astype(bool)
        & ~metricable["overlap_warning"].astype(bool)
        & ~metricable["low_power"].astype(bool)
    ].copy()
    safe_best = safe.sort_values("AUPRC", ascending=False).groupby(["benchmark_family", "dataset"], as_index=False).head(1)
    safe_best.to_csv(OUT / "biodarwin_reviewer_safe_best_per_dataset.tsv", sep="\t", index=False)

    plot_outputs(all_metrics, inventory)
    write_report_and_html(all_metrics, inventory, clean, landscape)

    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "out_dir": str(OUT),
        "n_inventory_rows": int(len(inventory)),
        "n_local_metric_rows": int(len(all_metrics)),
        "n_binary_metric_rows": int(len(metricable)),
        "n_reviewer_safe_metric_rows": int(len(safe)),
        "n_competitor_landscape_rows": int(len(landscape)),
        "best_per_dataset": best[["benchmark_family", "dataset", "algorithm", "AUPRC", "AUROC", "claim_status"]].to_dict("records"),
        "claim_boundary": "Mega comparison scaffold; diagnostic/overlap/low-power rows are not locked claims.",
        "outputs": {
            "metrics": str(OUT / "biodarwin_external_all_local_metrics.tsv"),
            "best": str(OUT / "biodarwin_best_per_dataset.tsv"),
            "safe_best": str(OUT / "biodarwin_reviewer_safe_best_per_dataset.tsv"),
            "landscape": str(OUT / "biodarwin_latest_competitor_landscape.tsv"),
            "html": str(HUB / "biodarwin_external_mega_comparison_2026_05_11.html"),
            "live_html": str(LIVE_HUB / "biodarwin_external_mega_comparison_2026_05_11.html"),
        },
    }
    (OUT / "biodarwin_external_mega_comparison_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
