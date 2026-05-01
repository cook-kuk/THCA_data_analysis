"""Honest pivot: 8-gene paper's REAL endpoint is RAI-responsiveness, not survival.
Test: DM2 vs DM1 RAI score, DICER1/EIF1AX × RAI score, histology distribution within DM.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project")
OUT = ROOT / "results" / "dark_matter_phase1"

# dark_matter_cohort.tsv has rai_score_v17 + tds16_score_v17 + histology
dm_cohort = pd.read_csv(ROOT / "results/v17/tables/dark_matter_cohort.tsv", sep="\t")
dm_cohort["tcga_short"] = dm_cohort["sample_id"].str[:12]

dm_master = pd.read_csv(OUT / "tcga_dm_master_with_pfi.tsv", sep="\t")
df = dm_master.merge(
    dm_cohort[["tcga_short", "rai_score_v17", "tds16_score_v17", "histology_subtype"]],
    on="tcga_short", how="left", suffixes=("", "_m")
)
print(f"Rows: {len(df)}")
print(f"DM patients with RAI score: {df.loc[df['dm_status'], 'rai_score_v17'].notna().sum()}")

# === 1. DM2 vs DM1 on RAI score within DM cohort ===
print("\n=== RAI score: DM2 vs DM1 within Dark Matter ===")
dm = df[df["dm_status"] & df["v17_dark_cluster"].isin(["DM1", "DM2"])].copy()
g1 = dm.loc[dm["v17_dark_cluster"] == "DM1", "rai_score_v17"].dropna()
g2 = dm.loc[dm["v17_dark_cluster"] == "DM2", "rai_score_v17"].dropna()
t = stats.ttest_ind(g1, g2, equal_var=False)
mw = stats.mannwhitneyu(g1, g2, alternative="two-sided")
print(f"  DM1 (n={len(g1)}): RAI mean={g1.mean():.3f}, median={g1.median():.3f}")
print(f"  DM2 (n={len(g2)}): RAI mean={g2.mean():.3f}, median={g2.median():.3f}")
print(f"  Welch t={t.statistic:.3f}, p={t.pvalue:.2e}")
print(f"  Mann-Whitney U p={mw.pvalue:.2e}")
print(f"  Cohen's d ≈ {(g1.mean() - g2.mean()) / np.sqrt((g1.std()**2 + g2.std()**2)/2):.2f}")

# === 2. RAI by alt-driver status within DM ===
print("\n=== RAI score: alt-driver+ vs alt-driver- within DM ===")
alt = ["DICER1_EIF1AX_PPM1D", "RET", "NTRK_fusion", "RET_fusion", "ALK_fusion", "TP53"]
dm["alt_pos"] = dm["driver_anchor_v17"].isin(alt)
gp = dm.loc[dm["alt_pos"], "rai_score_v17"].dropna()
gn = dm.loc[~dm["alt_pos"], "rai_score_v17"].dropna()
print(f"  Alt+ (n={len(gp)}): RAI mean={gp.mean():.3f}")
print(f"  Alt- (n={len(gn)}): RAI mean={gn.mean():.3f}")
if len(gp) >= 3:
    t2 = stats.ttest_ind(gp, gn, equal_var=False)
    print(f"  Welch t={t2.statistic:.3f}, p={t2.pvalue:.4f}")

# DICER1/EIF1AX/PPM1D specifically
dpp = dm.loc[dm["driver_anchor_v17"] == "DICER1_EIF1AX_PPM1D", "rai_score_v17"].dropna()
others = dm.loc[~(dm["driver_anchor_v17"] == "DICER1_EIF1AX_PPM1D"), "rai_score_v17"].dropna()
print(f"\n  DICER1/EIF1AX/PPM1D (n={len(dpp)}): RAI mean={dpp.mean():.3f}")
print(f"  Other DM (n={len(others)}): RAI mean={others.mean():.3f}")
if len(dpp) >= 3:
    t3 = stats.ttest_ind(dpp, others, equal_var=False)
    print(f"  Welch t={t3.statistic:.3f}, p={t3.pvalue:.4f}")

# === 3. Whole-cohort RAI: DM2-in-DM vs all rest ===
print("\n=== Whole-cohort RAI: DM2 (in DM) vs all-rest TCGA ===")
all_tcga = df.copy()
all_tcga["dm2_in_dm"] = all_tcga["dm_status"] & (all_tcga["v17_dark_cluster"] == "DM2")
g_dm2 = all_tcga.loc[all_tcga["dm2_in_dm"], "rai_score_v17"].dropna()
g_rest = all_tcga.loc[~all_tcga["dm2_in_dm"], "rai_score_v17"].dropna()
print(f"  DM2-in-DM (n={len(g_dm2)}): RAI mean={g_dm2.mean():.3f}")
print(f"  Rest TCGA (n={len(g_rest)}): RAI mean={g_rest.mean():.3f}")
t4 = stats.ttest_ind(g_dm2, g_rest, equal_var=False)
print(f"  Welch t={t4.statistic:.3f}, p={t4.pvalue:.2e}")
print(f"  Cohen's d ≈ {(g_dm2.mean() - g_rest.mean()) / np.sqrt((g_dm2.std()**2 + g_rest.std()**2)/2):.2f}")

# === 4. Histology distribution within DM ===
print("\n=== Histology within DM, by cluster ===")
hist = pd.crosstab(dm["v17_dark_cluster"], dm["histology_subtype"], dropna=False)
print(hist.to_string())

# === 5. TDS score (differentiation) ===
print("\n=== TDS16 score: DM2 vs DM1 within DM ===")
g1t = dm.loc[dm["v17_dark_cluster"] == "DM1", "tds16_score_v17"].dropna()
g2t = dm.loc[dm["v17_dark_cluster"] == "DM2", "tds16_score_v17"].dropna()
print(f"  DM1 (n={len(g1t)}): TDS16 mean={g1t.mean():.3f}")
print(f"  DM2 (n={len(g2t)}): TDS16 mean={g2t.mean():.3f}")
t5 = stats.ttest_ind(g1t, g2t, equal_var=False)
print(f"  Welch t={t5.statistic:.3f}, p={t5.pvalue:.2e}")

# === Summary save ===
result = {
    "rai_dm2_vs_dm1_within_DM": {
        "DM1_mean": float(g1.mean()), "DM1_n": len(g1),
        "DM2_mean": float(g2.mean()), "DM2_n": len(g2),
        "welch_p": float(t.pvalue),
        "cohen_d": float((g1.mean() - g2.mean()) / np.sqrt((g1.std()**2 + g2.std()**2)/2)),
    },
    "rai_dm2_vs_rest_whole_cohort": {
        "DM2_mean": float(g_dm2.mean()), "DM2_n": len(g_dm2),
        "rest_mean": float(g_rest.mean()), "rest_n": len(g_rest),
        "welch_p": float(t4.pvalue),
        "cohen_d": float((g_dm2.mean() - g_rest.mean()) / np.sqrt((g_dm2.std()**2 + g_rest.std()**2)/2)),
    },
    "rai_dicer1_eif1ax": {
        "alt_n": len(dpp), "alt_mean": float(dpp.mean()) if len(dpp) else None,
        "other_n": len(others), "other_mean": float(others.mean()),
        "welch_p": float(t3.pvalue) if len(dpp) >= 3 else None,
    },
    "tds_dm2_vs_dm1": {
        "DM1_mean": float(g1t.mean()), "DM2_mean": float(g2t.mean()), "welch_p": float(t5.pvalue),
    },
}
(OUT / "step2_rai_pivot.json").write_text(json.dumps(result, indent=2))
hist.to_csv(OUT / "step2_dm_histology.tsv", sep="\t")
print(f"\nSaved.")
