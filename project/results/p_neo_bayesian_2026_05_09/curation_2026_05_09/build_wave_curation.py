#!/usr/bin/env python3
"""Build a curated registry and fair-comparison pack for all neo waves.

This is intentionally conservative:
  - raw wave outputs are left untouched;
  - per-method metrics are recomputed from row-level predictions where possible;
  - Wave 8's large TCR/SelfSim signal is stress-tested for exact-reference artifacts.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
OUT = ROOT / "curation_2026_05_09"
PRED11 = ROOT / "wave11" / "wave11_predictions"
WAVE9 = ROOT / "wave9"
WAVE8 = ROOT / "wave8"
RNG_SEED = 42


def as_bool(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s
    return s.astype(str).str.lower().isin(["true", "1", "yes"])


def key_cols(df: pd.DataFrame, hla_col: str) -> pd.DataFrame:
    out = df.copy()
    out["peptide"] = out["peptide"].astype(str)
    out["hla"] = out[hla_col].astype(str)
    return out


def ece_score(y: np.ndarray, p: np.ndarray, bins: int = 10) -> float:
    y = np.asarray(y, dtype=float)
    p = np.asarray(p, dtype=float)
    p01 = np.clip(p, 0, 1)
    edges = np.linspace(0, 1, bins + 1)
    ece = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (p01 >= lo) & (p01 < hi if hi < 1 else p01 <= hi)
        if not mask.any():
            continue
        ece += mask.mean() * abs(y[mask].mean() - p01[mask].mean())
    return float(ece)


def boot_auc_ci(y: np.ndarray, p: np.ndarray, n_boot: int = 300, seed: int = RNG_SEED) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    y = np.asarray(y, dtype=int)
    p = np.asarray(p, dtype=float)
    vals: list[float] = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y), len(y))
        yy = y[idx]
        if yy.min() == yy.max():
            continue
        vals.append(roc_auc_score(yy, p[idx]))
    if len(vals) < 50:
        return math.nan, math.nan
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def paired_auc_delta(
    y: np.ndarray,
    p_method: np.ndarray,
    p_ref: np.ndarray,
    n_boot: int = 500,
    seed: int = RNG_SEED,
) -> tuple[float, float, float, float]:
    """Return observed delta, CI, and a two-sided paired-bootstrap sign p-value."""
    y = np.asarray(y, dtype=int)
    p_method = np.asarray(p_method, dtype=float)
    p_ref = np.asarray(p_ref, dtype=float)
    obs = roc_auc_score(y, p_method) - roc_auc_score(y, p_ref)
    rng = np.random.default_rng(seed)
    vals: list[float] = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y), len(y))
        yy = y[idx]
        if yy.min() == yy.max():
            continue
        vals.append(roc_auc_score(yy, p_method[idx]) - roc_auc_score(yy, p_ref[idx]))
    if len(vals) < 50:
        return float(obs), math.nan, math.nan, math.nan
    arr = np.asarray(vals)
    p_two = 2 * min(float((arr <= 0).mean()), float((arr >= 0).mean()))
    return float(obs), float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5)), min(1.0, p_two)


def load_bundle() -> pd.DataFrame:
    df = pd.read_csv(ROOT / "bundle.tsv", sep="\t", low_memory=False)
    df = key_cols(df, "HLA_norm")
    df["in_master"] = as_bool(df["in_master"])
    return df


def load_score_file(path: Path, method: str, score_col: str, hla_col: str = "hla") -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df = key_cols(df, hla_col)
    keep = ["peptide", "hla", "label", score_col]
    if "in_master" in df.columns:
        keep.append("in_master")
    if "split" in df.columns:
        keep.append("split")
    out = df[keep].copy()
    out = out.rename(columns={score_col: "score"})
    out["method"] = method
    out["score"] = pd.to_numeric(out["score"], errors="coerce")
    return out.dropna(subset=["score"])


def add_wave8_predictions(score_frames: list[pd.DataFrame]) -> dict[str, list[str]]:
    w8 = pd.read_csv(WAVE8 / "wave8_combined_features.tsv", sep="\t", low_memory=False)
    w8 = key_cols(w8, "hla")
    w8["in_master"] = as_bool(w8["in_master"])
    feature_blocks = {
        "Wave8_TCR_only": ["tcr_motif_score", "tcr_motif_count_log", "tcr_class_score"],
        "Wave8_TCR_motif_only": ["tcr_motif_score", "tcr_motif_count_log"],
        "Wave8_SelfSim_full": [
            "self_exact_match",
            "self_hamming1_count_log",
            "self_hamming2_count_log",
            "self_blosum_max",
        ],
        "Wave8_SelfSim_no_exact": [
            "self_hamming1_count_log",
            "self_hamming2_count_log",
            "self_blosum_max",
        ],
        "Wave8_TCR_SelfSim_full": [
            "tcr_motif_score",
            "tcr_motif_count_log",
            "tcr_class_score",
            "self_exact_match",
            "self_hamming1_count_log",
            "self_hamming2_count_log",
            "self_blosum_max",
        ],
        "Wave8_TCR_SelfSim_no_exact": [
            "tcr_motif_score",
            "tcr_motif_count_log",
            "tcr_class_score",
            "self_hamming1_count_log",
            "self_hamming2_count_log",
            "self_blosum_max",
        ],
    }
    train = w8[w8["split"] == "train"].copy()
    test = w8[w8["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])].copy()
    ytr = train["label"].astype(int).to_numpy()
    for method, cols in feature_blocks.items():
        clf = Pipeline(
            [
                ("scaler", StandardScaler()),
                ("lr", LogisticRegression(max_iter=2000, C=1.0, solver="lbfgs")),
            ]
        )
        clf.fit(train[cols].to_numpy(dtype=float), ytr)
        pred = test[["peptide", "hla", "label", "in_master", "split"]].copy()
        pred["score"] = clf.predict_proba(test[cols].to_numpy(dtype=float))[:, 1]
        pred["method"] = method
        score_frames.append(pred)
    return feature_blocks


def build_prediction_long() -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    wave11_specs = [
        ("MHCflurry", "MHCflurry__itsndb.tsv", "score_mhcflurry"),
        ("BigMHC_IM", "BigMHC_IM__itsndb.tsv", "score_bigmhc_im"),
        ("DeepImmuno", "DeepImmuno__itsndb.tsv", "score_deepimmuno"),
        ("PRIME", "PRIME__itsndb.tsv", "score_prime"),
        ("NetMHCpan_4.1", "NetMHCpan_4.1__itsndb.tsv", "score_netmhcpan"),
        ("TransPHLA", "TransPHLA__itsndb.tsv", "score_transphla"),
        ("RF_biophys", "RF_biophys__itsndb.tsv", "score_rf_biophys"),
        ("Structure_LR", "Structure_LR__itsndb.tsv", "score_structure_lr"),
        ("ESM2_Bayesian", "ESM2_Bayesian__itsndb.tsv", "score_esm2_bayesian"),
        ("VQC", "VQC__itsndb.tsv", "score_vqc"),
        ("GP_quantum", "GP_quantum__itsndb.tsv", "score_gp_quantum"),
    ]
    for method, fn, score_col in wave11_specs:
        path = PRED11 / fn
        if path.exists():
            frames.append(load_score_file(path, method, score_col))

    wave9_specs = [
        ("TSCAPE_TITANiAN", WAVE9 / "predictions_tscape.tsv", "score_tscape"),
        ("MHCnuggets_2", WAVE9 / "predictions_mhcnuggets.tsv", "score_mhcnuggets"),
        ("NetMHCstabpan", WAVE9 / "predictions_netmhcstabpan.tsv", "score_stabpan"),
    ]
    for method, path, score_col in wave9_specs:
        if path.exists():
            frames.append(load_score_file(path, method, score_col))

    w4c = ROOT / "wave4c" / "ensemble_predictions.tsv"
    if w4c.exists():
        for method, score_col in [
            ("Stack_mean_E1", "score_E1"),
            ("Stack_median_E2", "score_E2"),
            ("Stack_LR_inmaster_E3a", "score_E3a"),
            ("Stack_LR_OOF_E3b", "score_E3b"),
        ]:
            frames.append(load_score_file(w4c, method, score_col))

    w5b_specs = [
        ("MultiTask_A_only", ROOT / "wave5b" / "predictions_A_only.tsv"),
        ("MultiTask_AB", ROOT / "wave5b" / "predictions_AB.tsv"),
        ("MultiTask_AC", ROOT / "wave5b" / "predictions_AC.tsv"),
        ("MultiTask_AD", ROOT / "wave5b" / "predictions_AD.tsv"),
        ("MultiTask_full", ROOT / "wave5b" / "predictions_full.tsv"),
    ]
    for method, path in w5b_specs:
        if path.exists():
            frames.append(load_score_file(path, method, "pred_A_mean"))

    extra_specs = [
        ("kNN", ROOT / "wave4b" / "predictions_knn.tsv", "pred_knn", "HLA_norm"),
        ("TENT", ROOT / "wave4b" / "predictions_tent.tsv", "pred_tent", "HLA_norm"),
        ("W7A_full", ROOT / "wave7" / "predictions_combined.tsv", "p_w7a", "HLA_norm"),
        ("W7A_QK_only", ROOT / "wave7" / "predictions_combined.tsv", "p_w7a_qk", "HLA_norm"),
        ("W7B_stacked", ROOT / "wave7" / "predictions_combined.tsv", "p_w7b", "HLA_norm"),
    ]
    for method, path, score_col, hla_col in extra_specs:
        if path.exists():
            frames.append(load_score_file(path, method, score_col, hla_col=hla_col))

    feature_blocks = add_wave8_predictions(frames)
    pred = pd.concat(frames, ignore_index=True)

    bundle = load_bundle()
    its = bundle[bundle["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])][
        ["peptide", "hla", "label", "split", "in_master"]
    ].copy()
    pred = pred.drop(columns=[c for c in ["split", "in_master"] if c in pred.columns], errors="ignore")
    pred = pred.merge(its, on=["peptide", "hla", "label"], how="inner")
    pred["row_id"] = pred["peptide"] + "|" + pred["hla"]
    pred["score"] = pd.to_numeric(pred["score"], errors="coerce")
    pred = pred.dropna(subset=["score"]).drop_duplicates(["method", "row_id"])
    pred.to_csv(OUT / "curated_predictions_itsndb_long.tsv", sep="\t", index=False)

    wide = pred.pivot_table(
        index=["row_id", "peptide", "hla", "label", "split", "in_master"],
        columns="method",
        values="score",
        aggfunc="first",
    ).reset_index()
    wide.to_csv(OUT / "curated_predictions_itsndb_wide.tsv", sep="\t", index=False)

    (OUT / "wave8_feature_blocks.json").write_text(json.dumps(feature_blocks, indent=2))
    return pred, wide


FAMILY = {
    "MHCflurry": "off_the_shelf",
    "BigMHC_IM": "off_the_shelf",
    "DeepImmuno": "off_the_shelf",
    "PRIME": "off_the_shelf",
    "NetMHCpan_4.1": "off_the_shelf",
    "TransPHLA": "off_the_shelf",
    "TSCAPE_TITANiAN": "off_the_shelf_2025",
    "MHCnuggets_2": "off_the_shelf",
    "NetMHCstabpan": "off_the_shelf",
    "Structure_LR": "custom_biophys",
    "RF_biophys": "custom_biophys",
    "ESM2_Bayesian": "custom_deep",
    "MultiTask_A_only": "custom_deep",
    "MultiTask_AB": "custom_deep",
    "MultiTask_AC": "custom_deep",
    "MultiTask_AD": "custom_deep",
    "MultiTask_full": "custom_deep",
    "kNN": "custom_retrieval",
    "TENT": "custom_adaptation",
    "VQC": "quantum",
    "GP_quantum": "quantum",
    "W7A_full": "synthesis",
    "W7A_QK_only": "synthesis",
    "W7B_stacked": "synthesis",
    "Stack_mean_E1": "ensemble",
    "Stack_median_E2": "ensemble",
    "Stack_LR_inmaster_E3a": "ensemble",
    "Stack_LR_OOF_E3b": "ensemble",
}


def metric_rows(pred: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for subset_name, mask in [
        ("ITSNdb_no_overlap", ~pred["in_master"]),
        ("ITSNdb_in_master", pred["in_master"]),
        ("ITSNdb_combined", pd.Series(True, index=pred.index)),
    ]:
        sub_all = pred[mask].copy()
        for method, sub in sub_all.groupby("method"):
            y = sub["label"].astype(int).to_numpy()
            p = sub["score"].astype(float).to_numpy()
            n_pos = int(y.sum())
            n_neg = len(y) - n_pos
            if len(y) < 10 or n_pos < 5 or n_neg < 5:
                continue
            auc = roc_auc_score(y, p)
            lo, hi = boot_auc_ci(y, p)
            rows.append(
                {
                    "method": method,
                    "family": FAMILY.get(method, "wave8_candidate" if method.startswith("Wave8") else "custom"),
                    "subset": subset_name,
                    "n": len(y),
                    "n_pos": n_pos,
                    "AUROC": auc,
                    "AUROC_lo95": lo,
                    "AUROC_hi95": hi,
                    "AUPRC": average_precision_score(y, p),
                    "Brier": brier_score_loss(y, np.clip(p, 0, 1)),
                    "ECE": ece_score(y, p),
                }
            )
    metrics = pd.DataFrame(rows)

    # Pairwise paired-bootstrap deltas on no-overlap against MHCflurry and Wave8 full.
    wide = pred[~pred["in_master"]].pivot_table(
        index=["row_id", "label"], columns="method", values="score", aggfunc="first"
    ).reset_index()
    delta_rows = []
    delta_methods = {
        "MHCflurry",
        "Structure_LR",
        "BigMHC_IM",
        "TSCAPE_TITANiAN",
        "Wave8_TCR_only",
        "Wave8_TCR_motif_only",
        "Wave8_TCR_SelfSim_full",
        "Wave8_TCR_SelfSim_no_exact",
        "Wave8_SelfSim_no_exact",
    }
    for method in sorted(m for m in pred["method"].unique() if m in delta_methods):
        row = {"method": method}
        for ref in ["MHCflurry", "Wave8_TCR_SelfSim_full"]:
            if method == ref or method not in wide.columns or ref not in wide.columns:
                row[f"paired_n_vs_{ref}"] = math.nan
                row[f"delta_AUROC_vs_{ref}"] = math.nan
                row[f"delta_lo95_vs_{ref}"] = math.nan
                row[f"delta_hi95_vs_{ref}"] = math.nan
                row[f"paired_p_vs_{ref}"] = math.nan
                continue
            comp = wide.dropna(subset=[method, ref]).copy()
            y = comp["label"].astype(int).to_numpy()
            if len(comp) < 10 or y.min() == y.max():
                continue
            delta, lo, hi, p_two = paired_auc_delta(y, comp[method].to_numpy(), comp[ref].to_numpy())
            row[f"paired_n_vs_{ref}"] = len(comp)
            row[f"delta_AUROC_vs_{ref}"] = delta
            row[f"delta_lo95_vs_{ref}"] = lo
            row[f"delta_hi95_vs_{ref}"] = hi
            row[f"paired_p_vs_{ref}"] = p_two
        delta_rows.append(row)
    deltas = pd.DataFrame(delta_rows)
    metrics = metrics.merge(deltas, on="method", how="left")
    metrics = metrics.sort_values(["subset", "AUROC"], ascending=[True, False])
    metrics.to_csv(OUT / "curated_method_metrics.tsv", sep="\t", index=False)

    no_overlap = metrics[metrics["subset"] == "ITSNdb_no_overlap"].copy()
    no_overlap = no_overlap.sort_values("AUROC", ascending=False)
    no_overlap.to_csv(OUT / "curated_no_overlap_ranking.tsv", sep="\t", index=False)
    return metrics


def wave8_artifact_audit(wide: pd.DataFrame) -> pd.DataFrame:
    w8 = pd.read_csv(WAVE8 / "wave8_combined_features.tsv", sep="\t", low_memory=False)
    w8 = key_cols(w8, "hla")
    w8["row_id"] = w8["peptide"] + "|" + w8["hla"]
    its = w8[w8["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])].copy()
    its["in_master"] = as_bool(its["in_master"])
    score_cols = [
        "MHCflurry",
        "Structure_LR",
        "BigMHC_IM",
        "Wave8_TCR_only",
        "Wave8_TCR_motif_only",
        "Wave8_SelfSim_full",
        "Wave8_SelfSim_no_exact",
        "Wave8_TCR_SelfSim_full",
        "Wave8_TCR_SelfSim_no_exact",
    ]
    joined = its.merge(
        wide[["row_id"] + [c for c in score_cols if c in wide.columns]],
        on="row_id",
        how="left",
    )
    subsets = {
        "no_overlap_all": (~joined["in_master"]),
        "no_overlap_no_same_hla_tcr_exact": (~joined["in_master"]) & (joined["tcr_motif_score"] < 0.999),
        "no_overlap_no_global_or_same_hla_tcr_exact": (~joined["in_master"]) & (joined["tcr_motif_score"] < 0.499),
        "no_overlap_no_self_exact": (~joined["in_master"]) & (joined["self_exact_match"] < 0.5),
        "no_overlap_no_tcr_or_self_exact": (~joined["in_master"])
        & (joined["tcr_motif_score"] < 0.999)
        & (joined["self_exact_match"] < 0.5),
    }
    rows: list[dict] = []
    for subset, mask in subsets.items():
        dat = joined[mask].copy()
        base = {
            "subset": subset,
            "n_rows": len(dat),
            "n_pos": int(dat["label"].sum()) if len(dat) else 0,
            "tcr_same_hla_exact_n": int((dat["tcr_motif_score"] >= 0.999).sum()) if len(dat) else 0,
            "tcr_global_or_same_hla_exact_n": int((dat["tcr_motif_score"] >= 0.499).sum()) if len(dat) else 0,
            "self_exact_n": int((dat["self_exact_match"] >= 0.5).sum()) if len(dat) else 0,
        }
        for method in score_cols:
            if method not in dat.columns:
                continue
            sub = dat.dropna(subset=[method])
            y = sub["label"].astype(int).to_numpy()
            if len(sub) < 10 or y.sum() < 5 or len(y) - y.sum() < 5:
                auc = math.nan
            else:
                auc = roc_auc_score(y, sub[method].to_numpy(dtype=float))
            base[f"{method}__AUROC"] = auc
            base[f"{method}__n"] = len(sub)
        rows.append(base)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "wave8_artifact_stress.tsv", sep="\t", index=False)

    summary = {
        "interpretation": (
            "Wave8 is the only new block that beats MHCflurry on all no-overlap rows, "
            "but the signal weakens after removing exact TCR/self overlaps; treat as candidate, not locked."
        ),
        "subsets": out.to_dict(orient="records"),
    }
    (OUT / "wave8_artifact_stress.json").write_text(json.dumps(summary, indent=2))
    return out


def write_registry(metrics: pd.DataFrame, artifact: pd.DataFrame) -> None:
    no = metrics[metrics["subset"] == "ITSNdb_no_overlap"].sort_values("AUROC", ascending=False)
    top_lines = []
    for _, r in no.head(12).iterrows():
        top_lines.append(
            f"| {r['method']} | {r['family']} | {r['AUROC']:.3f} | "
            f"[{r['AUROC_lo95']:.3f}, {r['AUROC_hi95']:.3f}] | {int(r['n'])} |"
        )

    wave_rows = [
        ("wave1 / root", "ESM2-Bayesian + DANN", "ARCHIVE_NEGATIVE", "No-overlap AUROC ~0.41; useful as leakage-failure baseline."),
        ("wave2", "Calibration / OOD / Structure_LR / quantum", "KEEP_PARTIAL", "Keep Structure_LR and calibration; quantum is negative/no advantage."),
        ("wave2.5", "PWM feature source", "KEEP_AS_FEATURE", "Used by stacking and Wave7; not a standalone headline."),
        ("wave3_algorithm_sweep + wave3_*", "Off-the-shelf algorithm sweep", "KEEP", "Paper-grade baseline; MHCflurry remains clean single-method leader."),
        ("wave4a", "GroupDRO / MIRO / MoLE", "ARCHIVE_NEGATIVE", "DG variants inflate in-master and fail no-overlap."),
        ("wave4a_lora_fix", "LoRA fix skeleton", "NO_CLAIM", "No completed outputs on disk."),
        ("wave4b", "kNN / TENT / conformal", "KEEP_AS_NEGATIVE_CONTROL", "TENT modest but below chance; conformal shows undercoverage on shifted rows."),
        ("wave4c", "Stacked ensemble", "ARCHIVE_NEGATIVE", "No ensemble beats MHCflurry on no-overlap."),
        ("wave5a", "MHCflurry distillation", "STOPPED_INCOMPLETE", "Stopped during seed 2; no completed result table; no claim."),
        ("wave5b", "Multi-task auxiliary heads", "KEEP_AS_ABLATION", "Best full model ~0.48 no-overlap; regularizes but not enough."),
        ("wave5c", "Modern PEFT skeleton", "NO_CLAIM", "Scripts only; no results."),
        ("wave6", "Encoder hierarchy / MHCflurry features", "STOPPED_PARTIAL", "Partial negative results; useful as notes only."),
        ("wave7", "Synthesis classifier", "ARCHIVE_NEGATIVE", "Combining known strengths does not beat MHCflurry."),
        ("wave8", "TCR + self-similarity", "CANDIDATE_RETEST", "Only block exceeding MHCflurry; artifact stress test required."),
        ("wave9", "3 added tools", "KEEP", "T-SCAPE/MHCnuggets/stabpan extend method panel; none beats MHCflurry."),
        ("wave10", "Mega comparison", "KEEP_CANONICAL", "Primary existing fair-comparison table and figures."),
        ("wave11", "Expanded test sets", "KEEP_SUPPLEMENT", "Useful external-set support; many all-positive/low-power bundles flagged."),
        ("wave12 / wave15 / wave16", "Claude UI labels only", "ABSENT_ON_DISK", "No local output directories found after token-limit interruption."),
    ]

    lines = [
        "# Neo Wave Curation Registry — 2026-05-09",
        "",
        "## Immediate Process Cleanup",
        "",
        "Stopped stale or low-value background processes before curation:",
        "- `wave3_algorithm_sweep` PRIME tail jobs for A0217/B4102 already documented as timeout/skipped.",
        "- `wave6/track6b_simple_reps.py` and associated log watchers.",
        "- `wave5a/train_student_distill.py` and associated log watchers; output remains incomplete and is not claimable.",
        "- Remote `neo_wave4a` log tail watcher.",
        "",
        "No `p_neo_bayesian_2026_05_09` jobs were running after cleanup.",
        "",
        "## Canonical No-Overlap Ranking",
        "",
        "Recomputed from row-level predictions where available. Primary subset is `ITSNdb_no_overlap`.",
        "",
        "| Method | Family | AUROC | 95% CI | n |",
        "|---|---|---:|---:|---:|",
        *top_lines,
        "",
        "## Wave Disposition",
        "",
        "| Wave | Content | Disposition | Reason |",
        "|---|---|---|---|",
    ]
    for wave, content, disp, reason in wave_rows:
        lines.append(f"| {wave} | {content} | {disp} | {reason} |")

    art = artifact.set_index("subset")
    if "no_overlap_all" in art.index:
        all_row = art.loc["no_overlap_all"]
        strict_row = art.loc["no_overlap_no_tcr_or_self_exact"] if "no_overlap_no_tcr_or_self_exact" in art.index else None
        lines.extend(
            [
                "",
                "## Wave 8 Audit Read",
                "",
                f"All no-overlap rows: n={int(all_row['n_rows'])}, "
                f"TCR same-HLA exact hits={int(all_row['tcr_same_hla_exact_n'])}, "
                f"TCR same/global exact-or-near flag (`score>=0.499`)={int(all_row['tcr_global_or_same_hla_exact_n'])}, "
                f"self exact hits={int(all_row['self_exact_n'])}.",
            ]
        )
        if strict_row is not None:
            lines.append(
                f"After removing same-HLA TCR exact hits and self exact hits: n={int(strict_row['n_rows'])}; "
                f"Wave8 full AUROC={strict_row.get('Wave8_TCR_SelfSim_full__AUROC', math.nan):.3f}; "
                f"Wave8 no-exact AUROC={strict_row.get('Wave8_TCR_SelfSim_no_exact__AUROC', math.nan):.3f}; "
                f"MHCflurry AUROC={strict_row.get('MHCflurry__AUROC', math.nan):.3f}."
            )

    lines.extend(
        [
            "",
            "## Files",
            "",
            "- `curated_predictions_itsndb_long.tsv`: row-level predictions in long form.",
            "- `curated_predictions_itsndb_wide.tsv`: row-level predictions in wide form.",
            "- `curated_method_metrics.tsv`: combined / in-master / no-overlap metrics with paired bootstrap deltas.",
            "- `curated_no_overlap_ranking.tsv`: sorted no-overlap ranking.",
            "- `wave8_artifact_stress.tsv`: exact-overlap stress subsets for Wave 8.",
            "- `wave8_feature_blocks.json`: exact feature definitions used for Wave 8 refits.",
            "",
            "## Bottom Line",
            "",
            "For a paper-grade method comparison, use `wave10` plus this curation pack as the source of truth. "
            "MHCflurry / Structure_LR / BigMHC are the stable baselines. "
            "Wave8 is the only possible new paper-making signal, but it must be treated as a candidate until the exact-reference audit is resolved on additional external sets.",
        ]
    )
    (OUT / "WAVE_REGISTRY.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pred, wide = build_prediction_long()
    metrics = metric_rows(pred)
    artifact = wave8_artifact_audit(wide)
    write_registry(metrics, artifact)
    print(f"[done] wrote curation pack to {OUT}")


if __name__ == "__main__":
    main()
