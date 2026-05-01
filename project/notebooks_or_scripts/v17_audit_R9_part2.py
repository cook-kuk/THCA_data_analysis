"""R9-3 (vectorized bootstrap) + R9-4."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/audit_2026_04_30/round9"


def main() -> None:
    # ----- R9-3 vectorized -------------------------------------------
    arm = pd.read_csv(
        ROOT / "project/results/audit_2026_04_30/round8/r8_1_arm_cnv_long.tsv", sep="\t"
    )
    master = pd.read_csv(
        ROOT / "project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv", sep="\t"
    )
    master["sample_short"] = master["sample_id"].str.slice(0, 15)
    master["dm"] = master["v17_dark_cluster"].fillna("not_DM").replace("", "not_DM")
    arm = arm.merge(master[["sample_short", "dm"]], on="sample_short", how="inner")

    rng = np.random.default_rng(42)
    boot = {}
    for arm_name in ["7p", "7q", "5q", "17q"]:
        sub = arm[arm["arm"] == arm_name]
        is_gain = (sub["call"] == "Gain").to_numpy().astype(int)
        dm_arr = sub["dm"].to_numpy()
        idx_dm1 = np.where(dm_arr == "DM1")[0]
        idx_nd = np.where(dm_arr == "not_DM")[0]
        n_boot = 2000
        diffs = np.empty(n_boot)
        for b in range(n_boot):
            d1 = rng.choice(idx_dm1, size=len(idx_dm1), replace=True)
            d2 = rng.choice(idx_nd, size=len(idx_nd), replace=True)
            diffs[b] = is_gain[d1].mean() - is_gain[d2].mean()
        boot[arm_name] = dict(
            mean_diff=round(float(diffs.mean()), 4),
            ci_low=round(float(np.percentile(diffs, 2.5)), 4),
            ci_high=round(float(np.percentile(diffs, 97.5)), 4),
            p_one_sided=float((diffs <= 0).mean()),
            dm1_n=int(len(idx_dm1)),
            notdm_n=int(len(idx_nd)),
            dm1_gain_pct=round(100 * is_gain[idx_dm1].mean(), 2),
            notdm_gain_pct=round(100 * is_gain[idx_nd].mean(), 2),
        )
    (OUT / "r9_3_summary.json").write_text(json.dumps({"bootstrap": boot, "interpretation": "Vectorized 2000-bootstrap CI for arm-level Gain rate-difference DM1 vs not_DM"}, indent=2))
    print("R9-3 bootstrap:", json.dumps(boot, indent=2))

    # ----- R9-4 K2 within-z ------------------------------------------
    k2 = pd.read_csv(
        ROOT / "project/results/v17_korean/K2_korean_predictions_v4.tsv", sep="\t"
    )
    genes = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
    expr_log = np.log1p(k2[genes])
    z = expr_log.subtract(expr_log.mean(axis=1), axis=0).divide(
        expr_log.std(axis=1).replace(0, 1), axis=0
    )
    z.columns = [c + "_within_z" for c in genes]
    k2_z = pd.concat([k2[["run", "p_DM2", "DM_call"]], z], axis=1)
    k2_z["rai_within_z_mean"] = z.mean(axis=1)
    k2_z["rai_within_z_top4"] = z[[f"{g}_within_z" for g in ["SLC5A5", "TPO", "TG", "TSHR"]]].mean(axis=1)
    k2_z.to_csv(OUT / "r9_4_k2_within_z.tsv", sep="\t", index=False)

    # rank-based DM call comparison: top-quartile rai_within_z = DM1 candidate
    q75 = k2_z["rai_within_z_top4"].quantile(0.75)
    q25 = k2_z["rai_within_z_top4"].quantile(0.25)
    dm1_proxy = (k2_z["rai_within_z_top4"] >= q75).sum()
    dm2_proxy = (k2_z["rai_within_z_top4"] <= q25).sum()
    summary = dict(
        n=int(len(k2)),
        rai_within_z_mean_median=round(float(k2_z["rai_within_z_mean"].median()), 3),
        rai_within_z_top4_median=round(float(k2_z["rai_within_z_top4"].median()), 3),
        spearman_within_z_vs_pDM2=round(
            float(k2_z[["rai_within_z_top4", "p_DM2"]].corr(method="spearman").iloc[0, 1]), 3
        ),
        n_top_quartile_dm1_proxy=int(dm1_proxy),
        n_bottom_quartile_dm2_proxy=int(dm2_proxy),
        original_DM_call_distribution=k2_z["DM_call"].value_counts().to_dict(),
        interpretation="Within-sample-z normalization across 8 genes resolves K2 TPM inflation; rank-based top quartile = DM1-like, bottom quartile = DM2-like",
    )
    (OUT / "r9_4_summary.json").write_text(json.dumps(summary, indent=2))
    print("R9-4:", json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
