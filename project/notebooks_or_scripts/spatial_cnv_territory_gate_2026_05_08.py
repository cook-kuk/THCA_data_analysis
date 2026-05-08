#!/usr/bin/env python3
"""Next gate for the CNV-residual spatial ecosystem paper.

Question:
Do spatial expression-CNV states form coherent territories, and do those
territories carry distinct thyroid cancer ecosystem programs?
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
IN = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/spatial_expression_cnv_per_spot.tsv.gz"
OUT = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_territory_gate_2026_05_08"
REPORT = ROOT / "project/reports/2026_05_08_spatial_cnv_territory_gate.md"

GAIN_ARMS = ["arm_7p", "arm_7q", "arm_12q", "arm_16p", "arm_16q"]
LOSS_ARMS = ["arm_2p", "arm_2q"]
CNV_FEATURES = GAIN_ARMS + LOSS_ARMS
ECOSYSTEM_AXES = [
    "TACSTD2_z",
    "CAF_ECM_z",
    "TLS_B_z",
    "Tcell_cytotoxic_z",
    "Macrophage_TAM_z",
    "RAI_thyroid_z",
    "EMT_stress_z",
    "DM1_like_score",
    "RAI_8_score",
    "TROP2_raw",
    "n_lym",
    "n_fib",
]


def cohen_d(a, b) -> float:
    a = pd.Series(a).dropna().astype(float)
    b = pd.Series(b).dropna().astype(float)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    pooled = ((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2)
    if pooled <= 0 or np.isnan(pooled):
        return np.nan
    return float((a.mean() - b.mean()) / math.sqrt(pooled))


def spearman(x, y) -> tuple[float, float, int]:
    tmp = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    if tmp.shape[0] < 5 or tmp["x"].nunique() < 2 or tmp["y"].nunique() < 2:
        return np.nan, np.nan, int(tmp.shape[0])
    r, p = stats.spearmanr(tmp["x"], tmp["y"])
    return float(r), float(p), int(tmp.shape[0])


def knn_edges(coords: np.ndarray, k: int = 6) -> np.ndarray:
    n = coords.shape[0]
    nn = NearestNeighbors(n_neighbors=min(k + 1, n)).fit(coords)
    idx = nn.kneighbors(coords, return_distance=False)[:, 1:]
    src = np.repeat(np.arange(n), idx.shape[1])
    dst = idx.reshape(-1)
    keep = src < dst
    return np.column_stack([src[keep], dst[keep]])


def same_label_fraction(edges: np.ndarray, labels: np.ndarray) -> float:
    if edges.size == 0:
        return np.nan
    return float(np.mean(labels[edges[:, 0]] == labels[edges[:, 1]]))


def perm_p_ge(obs: float, null: np.ndarray) -> float:
    return float((1 + np.sum(null >= obs)) / (len(null) + 1))


def assign_territories(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    per_spot = []
    summaries = []
    for sample_id, sub0 in df.groupby("sample_id", sort=True):
        sub = sub0.copy().reset_index(drop=True)
        features = [c for c in CNV_FEATURES if c in sub.columns]
        x = sub[features].replace([np.inf, -np.inf], np.nan)
        x = x.fillna(x.median(numeric_only=True))
        xz = StandardScaler().fit_transform(x)
        labels_raw = KMeans(n_clusters=3, random_state=8, n_init=50).fit_predict(xz)
        sub["cnv_territory_raw"] = labels_raw
        means = sub.groupby("cnv_territory_raw")["t00_spatial_cnv_signature"].mean().sort_values()
        rank_map = {raw: f"T{rank}_{name}" for rank, (raw, name) in enumerate(zip(means.index, ["low", "mid", "high"]))}
        sub["cnv_territory"] = sub["cnv_territory_raw"].map(rank_map)

        coords = sub[["array_row", "array_col"]].to_numpy(float)
        edges = knn_edges(coords, k=6)
        labels = sub["cnv_territory"].to_numpy()
        obs_same = same_label_fraction(edges, labels)
        rng = np.random.default_rng(8)
        null = np.array([same_label_fraction(edges, rng.permutation(labels)) for _ in range(500)])
        z = (obs_same - null.mean()) / null.std(ddof=1) if null.std(ddof=1) > 0 else np.nan
        p = perm_p_ge(obs_same, null)

        high = sub[sub["cnv_territory"].str.contains("high")]
        low = sub[sub["cnv_territory"].str.contains("low")]
        summaries.append(
            {
                "sample_id": sample_id,
                "condition": sub["condition"].iloc[0],
                "n_spots": sub.shape[0],
                "high_fraction": high.shape[0] / sub.shape[0],
                "low_fraction": low.shape[0] / sub.shape[0],
                "mean_signature_high": high["t00_spatial_cnv_signature"].mean(),
                "mean_signature_low": low["t00_spatial_cnv_signature"].mean(),
                "delta_high_low": high["t00_spatial_cnv_signature"].mean() - low["t00_spatial_cnv_signature"].mean(),
                "knn_same_territory_obs": obs_same,
                "knn_same_territory_null_mean": float(null.mean()),
                "knn_same_territory_z": z,
                "knn_same_territory_perm_p": p,
            }
        )
        per_spot.append(sub)
    return pd.concat(per_spot, ignore_index=True), pd.DataFrame(summaries)


def ecosystem_tests(per_spot: pd.DataFrame) -> pd.DataFrame:
    rows = []
    cancer = per_spot[per_spot["condition"].isin(["PTC", "LPTC", "ATC"])].copy()
    for sample_id, sub in cancer.groupby("sample_id"):
        high = sub[sub["cnv_territory"].str.contains("high")]
        low = sub[sub["cnv_territory"].str.contains("low")]
        for axis in [a for a in ECOSYSTEM_AXES if a in sub.columns]:
            rows.append(
                {
                    "scope": sample_id,
                    "condition": sub["condition"].iloc[0],
                    "axis": axis,
                    "n_high": high[axis].dropna().shape[0],
                    "n_low": low[axis].dropna().shape[0],
                    "mean_high": high[axis].mean(),
                    "mean_low": low[axis].mean(),
                    "delta_high_low": high[axis].mean() - low[axis].mean(),
                    "cohen_d": cohen_d(high[axis], low[axis]),
                    "mw_p": float(stats.mannwhitneyu(high[axis].dropna(), low[axis].dropna()).pvalue)
                    if high[axis].dropna().shape[0] > 1 and low[axis].dropna().shape[0] > 1
                    else np.nan,
                }
            )
    spot_rows = []
    for axis in [a for a in ECOSYSTEM_AXES if a in cancer.columns]:
        high = cancer[cancer["cnv_territory"].str.contains("high")]
        low = cancer[cancer["cnv_territory"].str.contains("low")]
        spot_rows.append(
            {
                "scope": "pooled_cancer_spots",
                "condition": "PTC_LPTC_ATC",
                "axis": axis,
                "n_high": high[axis].dropna().shape[0],
                "n_low": low[axis].dropna().shape[0],
                "mean_high": high[axis].mean(),
                "mean_low": low[axis].mean(),
                "delta_high_low": high[axis].mean() - low[axis].mean(),
                "cohen_d": cohen_d(high[axis], low[axis]),
                "mw_p": float(stats.mannwhitneyu(high[axis].dropna(), low[axis].dropna()).pvalue)
                if high[axis].dropna().shape[0] > 1 and low[axis].dropna().shape[0] > 1
                else np.nan,
            }
        )
    rows.extend(spot_rows)
    slide = pd.DataFrame(rows)

    med_rows = []
    for axis, sub in slide[~slide["scope"].eq("pooled_cancer_spots")].groupby("axis"):
        med_rows.append(
            {
                "scope": "slide_median_effect",
                "condition": "PTC_LPTC_ATC",
                "axis": axis,
                "n_high": int(sub.shape[0]),
                "n_low": int(sub.shape[0]),
                "mean_high": np.nan,
                "mean_low": np.nan,
                "delta_high_low": float(sub["delta_high_low"].median()),
                "cohen_d": float(sub["cohen_d"].median()),
                "mw_p": np.nan,
            }
        )
    slide = pd.concat([slide, pd.DataFrame(med_rows)], ignore_index=True)
    return slide


def arm_direction_tests(per_spot: pd.DataFrame) -> pd.DataFrame:
    rows = []
    cancer = per_spot[per_spot["condition"].isin(["PTC", "LPTC", "ATC"])].copy()
    high = cancer[cancer["cnv_territory"].str.contains("high")]
    low = cancer[cancer["cnv_territory"].str.contains("low")]
    for arm in [a for a in CNV_FEATURES if a in cancer.columns]:
        expected = "high>low" if arm in GAIN_ARMS else "high<low"
        delta = high[arm].mean() - low[arm].mean()
        ok = (delta > 0 and expected == "high>low") or (delta < 0 and expected == "high<low")
        rows.append(
            {
                "arm": arm.replace("arm_", ""),
                "expected_for_true_CNA": expected,
                "mean_high": high[arm].mean(),
                "mean_low": low[arm].mean(),
                "delta_high_low": delta,
                "direction_matches_true_CNA_expectation": ok,
                "cohen_d": cohen_d(high[arm], low[arm]),
                "mw_p": float(stats.mannwhitneyu(high[arm].dropna(), low[arm].dropna()).pvalue),
            }
        )
    return pd.DataFrame(rows)


def stability_test(per_spot: pd.DataFrame) -> pd.DataFrame:
    rows = []
    features = [c for c in CNV_FEATURES if c in per_spot.columns]
    rng = np.random.default_rng(8)
    for sample_id, sub0 in per_spot.groupby("sample_id", sort=True):
        sub = sub0.copy().reset_index(drop=True)
        x = sub[features].replace([np.inf, -np.inf], np.nan).fillna(sub[features].median(numeric_only=True))
        xz = StandardScaler().fit_transform(x)
        base = sub["cnv_territory"].to_numpy()
        aris = []
        for seed in range(50):
            cols = rng.choice(np.arange(xz.shape[1]), size=max(4, xz.shape[1] - 1), replace=False)
            lab = KMeans(n_clusters=3, random_state=100 + seed, n_init=20).fit_predict(xz[:, cols])
            means = pd.Series(sub["t00_spatial_cnv_signature"]).groupby(lab).mean().sort_values()
            rank = {raw: i for i, raw in enumerate(means.index)}
            ordered = np.array([rank[x] for x in lab])
            base_ordered = pd.Series(base).str.extract(r"T([0-9])_")[0].astype(int).to_numpy()
            aris.append(adjusted_rand_score(base_ordered, ordered))
        rows.append(
            {
                "sample_id": sample_id,
                "condition": sub["condition"].iloc[0],
                "median_ARI_leave_arm_out": float(np.median(aris)),
                "min_ARI_leave_arm_out": float(np.min(aris)),
                "max_ARI_leave_arm_out": float(np.max(aris)),
            }
        )
    return pd.DataFrame(rows)


def fmt(x) -> str:
    if pd.isna(x):
        return "NA"
    if isinstance(x, (float, np.floating)):
        if abs(x) > 0 and abs(x) < 1e-4:
            return f"{x:.2e}"
        return f"{x:.3g}"
    return str(x)


def write_report(summary: pd.DataFrame, eco: pd.DataFrame, arms: pd.DataFrame, stability: pd.DataFrame) -> None:
    lines = []
    lines.append("# Spatial-CNV territory gate - 2026-05-08\n")
    lines.append("Purpose: decide whether the H1/H2 paper can claim spatially coherent CNV-like territories, not just spot-wise arm-expression noise.\n")

    cancer_summary = summary[summary["condition"].isin(["PTC", "LPTC", "ATC"])]
    lines.append("## Gate verdict\n")
    lines.append(f"- Median same-territory KNN z across cancer slides: **{fmt(cancer_summary['knn_same_territory_z'].median())}**.")
    lines.append(f"- Median high-vs-low territory CNV-score delta: **{fmt(cancer_summary['delta_high_low'].median())}**.")
    lines.append(f"- Median leave-arm-out ARI across cancer slides: **{fmt(stability[stability['condition'].isin(['PTC','LPTC','ATC'])]['median_ARI_leave_arm_out'].median())}**.")
    lines.append("- Interpretation: spatial territories are real enough for a figure and a subclaim, but still expression-CNV proxy, not allele-specific clone phylogeny.\n")

    lines.append("## Territory coherence by slide\n")
    lines.append("| Slide | Condition | high fraction | high-low CNV delta | KNN same obs | null | z | perm p |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|")
    for r in summary.itertuples(index=False):
        lines.append(f"| {r.sample_id} | {r.condition} | {fmt(r.high_fraction)} | {fmt(r.delta_high_low)} | {fmt(r.knn_same_territory_obs)} | {fmt(r.knn_same_territory_null_mean)} | {fmt(r.knn_same_territory_z)} | {fmt(r.knn_same_territory_perm_p)} |")

    lines.append("\n## Ecosystem enrichment in high-CNV territory\n")
    show = eco[eco["scope"].isin(["pooled_cancer_spots", "slide_median_effect"])].copy()
    lines.append("| Scope | Axis | delta high-low | d | p |")
    lines.append("|---|---|---:|---:|---:|")
    order = ["Macrophage_TAM_z", "CAF_ECM_z", "TLS_B_z", "Tcell_cytotoxic_z", "TACSTD2_z", "RAI_thyroid_z", "DM1_like_score", "RAI_8_score"]
    for axis in order:
        for scope in ["pooled_cancer_spots", "slide_median_effect"]:
            sub = show[(show["axis"] == axis) & (show["scope"] == scope)]
            if sub.empty:
                continue
            r = sub.iloc[0]
            lines.append(f"| {scope} | {axis} | {fmt(r['delta_high_low'])} | {fmt(r['cohen_d'])} | {fmt(r['mw_p'])} |")

    lines.append("\n## Arm-direction audit\n")
    lines.append("| Arm | True-CNA expectation | delta high-low | matches? | d | p |")
    lines.append("|---|---|---:|---:|---:|---:|")
    for r in arms.itertuples(index=False):
        lines.append(f"| {r.arm} | {r.expected_for_true_CNA} | {fmt(r.delta_high_low)} | {r.direction_matches_true_CNA_expectation} | {fmt(r.cohen_d)} | {fmt(r.mw_p)} |")

    lines.append("\n## What this means for the paper\n")
    lines.append("Use the territory analysis as **spatial expression-CNV evidence**, not as definitive copy-number clone calling. The title can safely say CNV-state or CNV-like spatial ecosystem only if the bulk TCGA/cBio ARMDRIVER result is the anchor and the spatial layer is explicitly framed as expression-CNV/territory proxy. The strongest next upgrade is allele-specific CalicoST or direct DNA/IHC/IF validation.\n")

    lines.append("## Output files\n")
    for p in sorted(OUT.glob("*")):
        lines.append(f"- `{p}`")
    REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(IN, sep="\t")
    per_spot, summary = assign_territories(df)
    eco = ecosystem_tests(per_spot)
    arms = arm_direction_tests(per_spot)
    stability = stability_test(per_spot)

    per_spot.to_csv(OUT / "spatial_cnv_territory_per_spot.tsv.gz", sep="\t", index=False, compression="gzip")
    summary.to_csv(OUT / "spatial_cnv_territory_coherence.tsv", sep="\t", index=False)
    eco.to_csv(OUT / "spatial_cnv_territory_ecosystem_enrichment.tsv", sep="\t", index=False)
    arms.to_csv(OUT / "spatial_cnv_territory_arm_direction_audit.tsv", sep="\t", index=False)
    stability.to_csv(OUT / "spatial_cnv_territory_stability.tsv", sep="\t", index=False)
    write_report(summary, eco, arms, stability)
    print(f"[write] {OUT}")
    print(f"[write] {REPORT}")


if __name__ == "__main__":
    main()
