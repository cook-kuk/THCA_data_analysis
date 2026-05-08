#!/usr/bin/env python3
"""Fast pathology bridge for ST-derived CNV territories.

Ground truth is the spatial transcriptome-derived T00 CNV territory produced by
`spatial_cnv_territory_gate_2026_05_08.py`. This tests whether pathology image
tiles can recover that territory in leave-one-slide-out evaluation.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import average_precision_score, balanced_accuracy_score, roc_auc_score
from sklearn.preprocessing import OneHotEncoder, StandardScaler


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/pathology_to_st_cnv_territory_fast_2026_05_08"
REPORT = ROOT / "project/reports/2026_05_08_pathology_to_st_cnv_territory_fast.md"

TERR = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_territory_gate_2026_05_08/spatial_cnv_territory_per_spot.tsv.gz"
UNI_META = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase1_gse250521/uni_embed_metadata_size224.tsv"
UNI_NPZ = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase1_gse250521/uni_embeddings_size224.npz"

MORPH_AXES = [
    "B1_stromal_encased_thyrocyte_idx",
    "B3_nuclear_eccentricity_var",
    "B4_tumor_stroma_ratio",
    "C2_spatial_celltype_entropy",
    "B5_thyrocyte_cluster_med",
    "n_thy",
    "n_lym",
    "n_fib",
]
CONT_TARGETS = [
    "t00_spatial_cnv_signature",
    "TACSTD2_z",
    "Macrophage_TAM_z",
]
BIN_TARGETS = ["t00_territory_high"]


def spearman(x, y) -> tuple[float, float, int]:
    tmp = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    if tmp.shape[0] < 5 or tmp["x"].nunique() < 2 or tmp["y"].nunique() < 2:
        return np.nan, np.nan, int(tmp.shape[0])
    r, p = stats.spearmanr(tmp["x"], tmp["y"])
    return float(r), float(p), int(tmp.shape[0])


def build_designs(dat: pd.DataFrame, emb: np.ndarray, train: np.ndarray, test: np.ndarray) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    uni_scaler = StandardScaler().fit(emb[train])
    xtr_uni0 = uni_scaler.transform(emb[train])
    xte_uni0 = uni_scaler.transform(emb[test])
    n_comp = min(64, xtr_uni0.shape[1], train.sum() - 1)
    pca = PCA(n_components=n_comp, random_state=8).fit(xtr_uni0)
    xtr_uni = pca.transform(xtr_uni0)
    xte_uni = pca.transform(xte_uni0)

    morph = dat[MORPH_AXES].replace([np.inf, -np.inf], np.nan).fillna(dat[MORPH_AXES].median(numeric_only=True)).to_numpy(float)
    morph_scaler = StandardScaler().fit(morph[train])
    xtr_morph = morph_scaler.transform(morph[train])
    xte_morph = morph_scaler.transform(morph[test])

    enc = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
    text_all = enc.fit_transform(dat[["condition"]])
    xtr_text = text_all[train]
    xte_text = text_all[test]

    return {
        "UNI_image_only": (xtr_uni, xte_uni),
        "cell_morph_only": (xtr_morph, xte_morph),
        "stage_text_only": (xtr_text, xte_text),
        "UNI_plus_cell_morph": (np.hstack([xtr_uni, xtr_morph]), np.hstack([xte_uni, xte_morph])),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    terr = pd.read_csv(TERR, sep="\t")
    terr["t00_territory_high"] = terr["cnv_territory"].astype(str).str.contains("high").astype(int)
    terr["t00_norm_q95_high"] = terr["t00_signature_norm_q95_high"].astype(int)

    meta = pd.read_csv(UNI_META, sep="\t").rename(columns={"slide": "sample_id"})
    emb = np.load(UNI_NPZ)["embeddings"]
    dat = meta.merge(
        terr[["sample_id", "spot_id", "condition", "cnv_territory"] + MORPH_AXES + CONT_TARGETS + BIN_TARGETS],
        on=["sample_id", "spot_id"],
        how="inner",
    )
    emb = emb[dat.index.to_numpy()]

    pred_frames = []
    for held in sorted(dat["sample_id"].unique()):
        train = (dat["sample_id"] != held).to_numpy()
        test = ~train
        if train.sum() < 100 or test.sum() < 20:
            continue
        designs = build_designs(dat, emb, train, test)
        test_df = dat.loc[test].reset_index(drop=True)
        for mode, (xtr, xte) in designs.items():
            for target in CONT_TARGETS:
                ytr = dat.loc[train, target].to_numpy(float)
                yte = dat.loc[test, target].to_numpy(float)
                model = Ridge(alpha=10.0)
                model.fit(xtr, ytr)
                pred = model.predict(xte)
                pred_frames.append(
                    pd.DataFrame(
                        {
                            "sample_id": test_df["sample_id"].to_numpy(),
                            "spot_id": test_df["spot_id"].to_numpy(),
                            "condition": test_df["condition"].to_numpy(),
                            "cnv_territory": test_df["cnv_territory"].to_numpy(),
                            "mode": mode,
                            "target": target,
                            "observed": yte,
                            "predicted": pred,
                        }
                    )
                )
            for target in BIN_TARGETS:
                ytr = dat.loc[train, target].astype(int).to_numpy()
                yte = dat.loc[test, target].astype(int).to_numpy()
                if len(np.unique(ytr)) < 2:
                    continue
                model = LogisticRegression(max_iter=500, class_weight="balanced", C=0.2, solver="liblinear")
                model.fit(xtr, ytr)
                prob = model.predict_proba(xte)[:, 1]
                pred_frames.append(
                    pd.DataFrame(
                        {
                            "sample_id": test_df["sample_id"].to_numpy(),
                            "spot_id": test_df["spot_id"].to_numpy(),
                            "condition": test_df["condition"].to_numpy(),
                            "cnv_territory": test_df["cnv_territory"].to_numpy(),
                            "mode": mode,
                            "target": target,
                            "observed": yte,
                            "predicted": prob,
                        }
                    )
                )

    pred = pd.concat(pred_frames, ignore_index=True)
    pred.to_csv(OUT / "pathology_to_st_cnv_predictions.tsv.gz", sep="\t", index=False, compression="gzip")

    metric_rows = []
    for (mode, target), sub in pred.groupby(["mode", "target"]):
        if target in CONT_TARGETS:
            r, p, n = spearman(sub["observed"], sub["predicted"])
            metric_rows.append({"mode": mode, "target": target, "metric": "spearman", "value": r, "p": p, "n": n})
        else:
            y = sub["observed"].astype(int)
            if y.nunique() < 2:
                auc = ap = bacc = np.nan
            else:
                auc = roc_auc_score(y, sub["predicted"])
                ap = average_precision_score(y, sub["predicted"])
                bacc = balanced_accuracy_score(y, sub["predicted"] >= 0.5)
            metric_rows.append({"mode": mode, "target": target, "metric": "AUROC", "value": auc, "p": np.nan, "n": sub.shape[0]})
            metric_rows.append({"mode": mode, "target": target, "metric": "average_precision", "value": ap, "p": np.nan, "n": sub.shape[0]})
            metric_rows.append({"mode": mode, "target": target, "metric": "balanced_accuracy@0.5", "value": bacc, "p": np.nan, "n": sub.shape[0]})
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(OUT / "pathology_to_st_cnv_metrics.tsv", sep="\t", index=False)

    per_slide = []
    for (mode, target, sample_id), sub in pred.groupby(["mode", "target", "sample_id"]):
        if target in CONT_TARGETS:
            r, p, n = spearman(sub["observed"], sub["predicted"])
            per_slide.append({"mode": mode, "target": target, "sample_id": sample_id, "metric": "spearman", "value": r, "n": n})
        else:
            y = sub["observed"].astype(int)
            auc = roc_auc_score(y, sub["predicted"]) if y.nunique() == 2 else np.nan
            per_slide.append({"mode": mode, "target": target, "sample_id": sample_id, "metric": "AUROC", "value": auc, "n": sub.shape[0]})
    pd.DataFrame(per_slide).to_csv(OUT / "pathology_to_st_cnv_per_slide_metrics.tsv", sep="\t", index=False)

    write_report(dat, metrics)
    print(f"[write] {OUT}")
    print(f"[write] {REPORT}")


def fmt(x) -> str:
    if pd.isna(x):
        return "NA"
    if isinstance(x, (float, np.floating)):
        if abs(x) > 0 and abs(x) < 1e-4:
            return f"{x:.2e}"
        return f"{x:.3g}"
    return str(x)


def best(metrics: pd.DataFrame, target: str, metric: str, mode: str | None = None):
    sub = metrics[(metrics["target"] == target) & (metrics["metric"] == metric)].copy()
    if mode:
        sub = sub[sub["mode"] == mode]
    if sub.empty:
        return None
    return sub.sort_values("value", ascending=False).iloc[0]


def write_report(dat: pd.DataFrame, metrics: pd.DataFrame) -> None:
    lines = []
    lines.append("# Pathology to ST-CNV territory bridge - 2026-05-08\n")
    lines.append("Ground truth: spatial transcriptome-derived T00 CNV territory from the territory gate. Prediction: leave-one-slide-out from UNI tile embeddings and image-derived cell morphology features.\n")

    uni_auc = best(metrics, "t00_territory_high", "AUROC", "UNI_image_only")
    uni_sig = best(metrics, "t00_spatial_cnv_signature", "spearman", "UNI_image_only")
    cell_auc = best(metrics, "t00_territory_high", "AUROC", "cell_morph_only")
    text_auc = best(metrics, "t00_territory_high", "AUROC", "stage_text_only")
    best_auc = best(metrics, "t00_territory_high", "AUROC")

    lines.append("## Verdict\n")
    lines.append(f"- UNI image-only to ST-CNV high territory: AUROC **{fmt(uni_auc['value'] if uni_auc is not None else np.nan)}**.")
    lines.append(f"- UNI image-only to continuous ST-CNV signature: rho **{fmt(uni_sig['value'] if uni_sig is not None else np.nan)}**.")
    lines.append(f"- Cell morphology-only AUROC: **{fmt(cell_auc['value'] if cell_auc is not None else np.nan)}**.")
    lines.append(f"- Stage/text-only AUROC: **{fmt(text_auc['value'] if text_auc is not None else np.nan)}**.")
    lines.append(f"- Best overall mode: **{best_auc['mode'] if best_auc is not None else 'NA'}**, AUROC **{fmt(best_auc['value'] if best_auc is not None else np.nan)}**.")
    if uni_auc is not None and float(uni_auc["value"]) >= 0.70:
        lines.append("- Interpretation: pathology image bridge is real enough to become a major paper layer.")
    elif uni_auc is not None and float(uni_auc["value"]) >= 0.60:
        lines.append("- Interpretation: pathology has a weak/moderate signal, useful as a companion figure but not the main claim.")
    else:
        lines.append("- Interpretation: pathology image alone does not recover ST-CNV territory strongly; this is an important negative boundary.")

    lines.append("\n## Data")
    lines.append(f"- Tile embeddings evaluated: **{dat.shape[0]}** spots.")
    lines.append(f"- Slides: **{dat['sample_id'].nunique()}**.")
    lines.append(f"- T00-high territory prevalence in evaluated tiles: **{dat['t00_territory_high'].mean():.1%}**.\n")

    lines.append("## Main Metrics\n")
    lines.append("| Mode | Target | Metric | Value | n |")
    lines.append("|---|---|---|---:|---:|")
    show_targets = ["t00_territory_high", "t00_norm_q95_high", "t00_spatial_cnv_signature", "TACSTD2_z", "Macrophage_TAM_z", "TLS_B_z"]
    show_modes = ["UNI_image_only", "cell_morph_only", "stage_text_only", "UNI_plus_cell_morph", "UNI_plus_stage_text", "UNI_cell_stage"]
    for target in show_targets:
        for mode in show_modes:
            sub = metrics[(metrics["target"] == target) & (metrics["mode"] == mode)]
            if target in BIN_TARGETS:
                sub = sub[sub["metric"] == "AUROC"]
            else:
                sub = sub[sub["metric"] == "spearman"]
            if sub.empty:
                continue
            r = sub.iloc[0]
            lines.append(f"| {mode} | {target} | {r['metric']} | {fmt(r['value'])} | {int(r['n'])} |")

    lines.append("\n## Paper Meaning\n")
    lines.append("If image-only stays modest, the paper should not claim H&E can replace ST-CNV mapping. The stronger high-IF framing is: spatial transcriptomics reveals CNV-like territories; pathology can partly read ecosystem programs such as TROP2/TAM/TLS, but direct ST/CNV measurement remains needed for the CNV state itself.\n")

    lines.append("## Output files\n")
    for p in sorted(OUT.glob("*")):
        lines.append(f"- `{p}`")
    REPORT.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
