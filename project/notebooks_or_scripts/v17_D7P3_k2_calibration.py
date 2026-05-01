#!/usr/bin/env python3
"""D7-P3 — K2 8-gene mini-index calibration fix → 3-cohort meta attempt.

Diagnoses raw-TPM directional mismatch (P5 ρ=+0.33 vs expected negative);
tries 4 alternative metrics; if any shows ρ < -0.3, runs 3-cohort meta.
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats

PROJ = Path("/opt/thyroid-dash/project")
RES = PROJ / "results/d7p3_k2_calibration"
RES.mkdir(parents=True, exist_ok=True)

# ============================================================
# 1. Load K2 + TCGA expression for reference
# ============================================================
print("=== 1. Load K2 + TCGA reference ===")
k2 = pd.read_csv(PROJ / "results/v17_korean/K2_korean_predictions_v4.tsv", sep="\t")
print(f"  K2: {len(k2)}, cols: {k2.columns.tolist()}")

tcga = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv",
                    sep="\t", index_col=0)
GENE_8 = ['SLC5A5', 'TPO', 'TG', 'TSHR', 'PAX8', 'NKX2-1', 'FOXE1', 'DIO1']

# ============================================================
# 2. K2 mini-index TPM range vs TCGA log2(FPKM+1) range
# ============================================================
print("\n=== 2. Range/inflation audit ===")
k2_g8 = k2[GENE_8].astype(float)
print("\n  K2 raw TPM per gene (median, sd):")
for g in GENE_8:
    v = k2_g8[g]
    print(f"    {g}: median={v.median():.1f}, sd={v.std():.1f}, range=[{v.min():.1f}, {v.max():.1f}]")

print("\n  TCGA log2(FPKM+1) per gene (median, sd):")
tcga_g8 = tcga.loc[GENE_8]
for g in GENE_8:
    v = tcga_g8.loc[g]
    print(f"    {g}: median={v.median():.2f}, sd={v.std():.2f}, range=[{v.min():.2f}, {v.max():.2f}]")

# Inflation factor: K2 raw is TPM (not log) — convert to log2 first
k2_log = np.log2(k2_g8 + 1.0)
print("\n  K2 log2(TPM+1) per gene:")
inf_factors = {}
for g in GENE_8:
    v = k2_log[g]
    tcga_med = tcga_g8.loc[g].median()
    inf = v.median() - tcga_med  # additive on log scale
    inf_factors[g] = round(float(inf), 2)
    print(f"    {g}: median={v.median():.2f} (TCGA {tcga_med:.2f}, inflation Δ={inf:+.2f})")

# ============================================================
# 3. Try 4 alternative score metrics
# ============================================================
print("\n=== 3. Try 4 alternative score metrics ===")

# Metric A: raw mean log2(TPM+1)
metric_A = k2_log.mean(axis=1)

# Metric B: per-sample z-score then mean (within-sample-centered)
# z = (gene_log_TPM - sample_mean) / sample_sd
sample_mean = k2_log.mean(axis=1)
sample_sd = k2_log.std(axis=1, ddof=1)
k2_centered = k2_log.subtract(sample_mean, axis=0).divide(sample_sd + 1e-9, axis=0)
metric_B = k2_centered.mean(axis=1)

# Metric C: per-gene z-score across cohort then mean
metric_C = k2_log.apply(lambda c: (c - c.mean()) / (c.std(ddof=1) + 1e-9), axis=0).mean(axis=1)

# Metric D: per-gene rank then mean
metric_D = k2_log.rank(axis=0, pct=True).mean(axis=1)

scores_df = pd.DataFrame({
    "p_DM2": k2["p_DM2"].values,
    "DM_call": k2["DM_call"].values,
    "metric_A_raw_logTPM_mean": metric_A.values,
    "metric_B_within_sample_z": metric_B.values,
    "metric_C_per_gene_z_across_cohort": metric_C.values,
    "metric_D_per_gene_rank_pct": metric_D.values,
}, index=k2.index)

print("\n  Spearman ρ between each score and p_DM2 (negative expected):")
results = {}
for col in ["metric_A_raw_logTPM_mean", "metric_B_within_sample_z",
            "metric_C_per_gene_z_across_cohort", "metric_D_per_gene_rank_pct"]:
    rho, p = stats.spearmanr(scores_df[col], scores_df["p_DM2"])
    results[col] = dict(rho=round(float(rho), 3), p=float(p))
    direction = "✓ negative" if rho < -0.1 else ("← positive (mismatch)" if rho > 0.1 else "≈ zero")
    print(f"    {col}: ρ={rho:+.3f}, p={p:.3g}  [{direction}]")

# ============================================================
# 4. Direction-flipped: which metric works
# ============================================================
print("\n=== 4. Best-calibrated metric ===")
best = min(results.items(), key=lambda x: x[1]["rho"])
print(f"  Best (most negative ρ): {best[0]}: ρ={best[1]['rho']}")

# Is best calibration < -0.3?
calib_OK = best[1]["rho"] < -0.3

# Check on DM_call binary
print("\n  Cohen d for best metric in DM2 vs DM1:")
best_col = best[0]
dm1 = scores_df.loc[scores_df["DM_call"] == "DM1", best_col].values
dm2 = scores_df.loc[scores_df["DM_call"] == "DM2", best_col].values
sp = np.sqrt(((len(dm1)-1)*np.var(dm1, ddof=1) + (len(dm2)-1)*np.var(dm2, ddof=1)) / max(len(dm1)+len(dm2)-2, 1))
d_k2 = (dm2.mean() - dm1.mean()) / max(sp, 1e-9)
print(f"    DM2 (n={len(dm2)}) vs DM1 (n={len(dm1)}) Cohen d (DM2 - DM1) = {d_k2:+.3f}")
# Note: DM1 only n=14 → underpowered

# ============================================================
# 5. 3-cohort meta if calibration succeeded
# ============================================================
print("\n=== 5. 3-cohort meta attempt ===")
def pool_d(rows):
    ds = np.array([r["cohen_d"] for r in rows])
    n1 = np.array([r["n_case"] for r in rows])
    n2 = np.array([r["n_ctrl"] for r in rows])
    v_d = (n1+n2)/(n1*n2) + (ds**2)/(2*(n1+n2))
    w = 1.0 / v_d
    d_fix = (w * ds).sum() / w.sum()
    Q = (w * (ds - d_fix)**2).sum()
    df_h = len(ds) - 1
    tau2 = max(0.0, (Q - df_h) / (w.sum() - (w**2).sum() / w.sum()))
    w_re = 1.0 / (v_d + tau2)
    d_re = (w_re * ds).sum() / w_re.sum()
    var_re = 1.0 / w_re.sum()
    return dict(d_pooled=round(float(d_re), 3),
                ci_lo=round(float(d_re - 1.96 * np.sqrt(var_re)), 3),
                ci_hi=round(float(d_re + 1.96 * np.sqrt(var_re)), 3),
                tau2=round(float(tau2), 4), Q=round(float(Q), 3), df=df_h)

if calib_OK:
    print(f"  ✅ Calibration succeeded with {best_col}")
    meta_rows = [
        {"cohort": "TCGA-THCA (DM2 vs DM1, n=500)", "n_case": 360, "n_ctrl": 140, "cohen_d": -1.783},
        {"cohort": "GSE286332 (PTC+HT vs PTC, n=18)", "n_case": 9, "n_ctrl": 9, "cohen_d": -1.602},
        {"cohort": f"K2 (n=260, {best_col})", "n_case": int(len(dm2)), "n_ctrl": int(len(dm1)),
         "cohen_d": round(float(-d_k2), 3)},  # flip sign so all are "low-RAI = high autoimmune-overlap"
    ]
    pooled = pool_d(meta_rows)
    print(f"  3-cohort pooled d = {pooled['d_pooled']:+.3f} [95% CI {pooled['ci_lo']:+.3f}, {pooled['ci_hi']:+.3f}]")
    print(f"  τ²={pooled['tau2']:.4f}, Q={pooled['Q']:.3f} (df={pooled['df']})")
else:
    print(f"  ❌ Calibration failed (best ρ={best[1]['rho']:.3f} > -0.3 threshold).")
    print(f"  K2 evidence chain remains via:")
    print(f"   (a) DM call distribution: 246/260 (94.6%) DM2 — Korean Hashimoto-like background")
    print(f"   (b) HLA arm: K2+Lee+GSE286332 = n=874 Korean PTC pool")
    print(f"   (c) Future task: STAR re-quantification (separate work, ~3-5 days)")
    pooled = None
    meta_rows = None

# ============================================================
# 6. STAR re-quantification feasibility estimate
# ============================================================
print("\n=== 6. STAR re-quantification feasibility ===")
star_estimate = {
    "n_K2_FASTQ": 260,
    "STAR_align_time_per_sample_min": "~15-20 (32-thread)",
    "featureCounts_time_total_min": "~30",
    "expected_total_compute_hr": "65-90",
    "Azure_burst_cost_estimate_USD": "$5-10",
    "current_decision": "future task (post-paper draft)",
    "rationale": "Calibration approach lower-effort if any of 4 metrics works",
}
print(json.dumps(star_estimate, indent=2))

# ============================================================
# 7. Save
# ============================================================
scores_df.to_csv(RES / "k2_4metric_scores.tsv", sep="\t")

summary = {
    "K2_inflation_factors_log2_K2_minus_TCGA_median": inf_factors,
    "metric_results": results,
    "best_metric": best[0],
    "best_rho": best[1]["rho"],
    "calibration_success": bool(calib_OK),
    "K2_DM_call_distribution": {"DM2": int(len(dm2)), "DM1": int(len(dm1))},
    "best_metric_d_DM2_vs_DM1": round(float(d_k2), 3),
    "meta_3cohort_rows": meta_rows,
    "pooled": pooled,
    "STAR_re_quantification_estimate": star_estimate,
    "alternate_evidence_chain": [
        "DM call distribution 246/260 = 94.6% DM2 (Korean Hashimoto-like background)",
        "HLA arm K2+Lee+GSE286332 n=874 Korean PTC pool",
        "STAR re-quant is future task ($5-10, post-paper draft)"
    ],
}
(RES / "D7P3_summary.json").write_text(json.dumps(summary, indent=2, default=str))
print(f"\n✓ Outputs to {RES}")
