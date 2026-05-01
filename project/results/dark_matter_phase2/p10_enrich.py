"""Phase 9 ENRICH: 5 simple analyses for paper rigor."""
import warnings
warnings.filterwarnings("ignore")
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project")
OUT = ROOT / "results/dark_matter_phase2"
FIG = OUT / "web/figures"
DATA_OUT = OUT / "web/data"

expr = pd.read_csv(ROOT / "data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv", sep="\t", index_col=0)
master = pd.read_csv(OUT / "../dark_matter_phase1/tcga_dm_master_with_pfi.tsv", sep="\t")
master_v17 = pd.read_csv(ROOT / "results/v17/tables/sample_master_v17_full.tsv", sep="\t", low_memory=False)
master_v17["tcga_short"] = master_v17["sample_id"].str[:12]
master_v17_tcga = master_v17[(master_v17["dataset"] == "TCGA-THCA") & (master_v17["normal_vs_tumor"] == "tumor")]
master = master.merge(master_v17_tcga[["tcga_short", "histology_subtype", "sample_id"]].drop_duplicates("tcga_short"),
                      on="tcga_short", how="left")
master["expr_col"] = master["sample_id"]
master = master[master["expr_col"].isin(expr.columns)]

PANEL = ["DIO1", "FOXE1", "NKX2-1", "PAX8", "SLC5A5", "TG", "TPO", "TSHR"]

# ============================================================================
# 1. Brown-Forsythe / Levene differential variability
# ============================================================================
print("=== 1. Brown-Forsythe variability test (DM2 more homogeneous?) ===")
dm = master[master["dm_status"] & master["v17_dark_cluster"].isin(["DM1", "DM2"])].copy()

# Per-gene Brown-Forsythe
bf_results = []
for gene in PANEL:
    if gene not in expr.index: continue
    v_dm1 = expr.loc[gene, dm.loc[dm["v17_dark_cluster"] == "DM1", "expr_col"]].values
    v_dm2 = expr.loc[gene, dm.loc[dm["v17_dark_cluster"] == "DM2", "expr_col"]].values
    bf_stat, bf_p = stats.levene(v_dm1, v_dm2, center="median")  # Brown-Forsythe = Levene with median
    bf_results.append({"gene": gene, "var_DM1": float(np.var(v_dm1)), "var_DM2": float(np.var(v_dm2)),
                       "var_ratio_DM2/DM1": float(np.var(v_dm2) / np.var(v_dm1)),
                       "BF_stat": float(bf_stat), "BF_p": float(bf_p)})
bf_df = pd.DataFrame(bf_results)
print(bf_df.round(4).to_string(index=False))
bf_df.to_csv(DATA_OUT / "brown_forsythe_variability.tsv", sep="\t", index=False)

# Pathway-level variability
PATHWAYS = {
    "Thyroid_diff": ["TG", "TPO", "TSHR", "DIO1", "DIO2", "SLC5A5", "PAX8", "FOXE1", "NKX2-1"],
    "Immune_T": ["CD3D", "CD3E", "CD8A", "CD4", "GZMB"],
    "EMT": ["VIM", "FN1", "SNAI1", "SNAI2", "ZEB1", "ZEB2"],
    "let7_targets": ["HMGA2", "LIN28B", "MYC"],
}
pw_var = []
for pw, genes in PATHWAYS.items():
    avail = [g for g in genes if g in expr.index]
    if len(avail) < 3: continue
    score = expr.loc[avail].mean(axis=0)
    s_dm1 = score.loc[dm.loc[dm["v17_dark_cluster"] == "DM1", "expr_col"]].values
    s_dm2 = score.loc[dm.loc[dm["v17_dark_cluster"] == "DM2", "expr_col"]].values
    bf_stat, bf_p = stats.levene(s_dm1, s_dm2, center="median")
    pw_var.append({"pathway": pw, "var_DM1": float(np.var(s_dm1)), "var_DM2": float(np.var(s_dm2)),
                   "ratio_DM2/DM1": float(np.var(s_dm2) / np.var(s_dm1)),
                   "BF_p": float(bf_p)})
pw_var_df = pd.DataFrame(pw_var)
print(f"\nPathway variability:\n{pw_var_df.round(4).to_string(index=False)}")

# Plot
fig, ax = plt.subplots(figsize=(10, 5))
all_var = pd.concat([bf_df.assign(level="gene").rename(columns={"gene": "label", "var_ratio_DM2/DM1": "ratio"}),
                      pw_var_df.assign(level="pathway").rename(columns={"pathway": "label", "ratio_DM2/DM1": "ratio"})], ignore_index=True)
all_var = all_var.sort_values("ratio")
colors = ["#2a9d8f" if r < 1 else "#e76f51" for r in all_var["ratio"]]
ax.barh(all_var["label"], np.log2(all_var["ratio"]), color=colors, alpha=0.85)
for i, (_, r) in enumerate(all_var.iterrows()):
    sig = " *" if r["BF_p"] < 0.05 else ""
    ax.text(np.log2(r["ratio"]) + (0.02 if r["ratio"] >= 1 else -0.02), i,
            f"{r['ratio']:.2f}{sig}", va="center", fontsize=9,
            ha="left" if r["ratio"] >= 1 else "right")
ax.axvline(0, color="black", lw=0.5)
ax.set_xlabel("log2 (var_DM2 / var_DM1) — negative = DM2 more homogeneous")
ax.set_title("Figure E30 — Differential variability DM1 vs DM2 (Brown-Forsythe)\n* p<0.05")
ax.grid(alpha=0.3, axis="x")
fig.tight_layout()
fig.savefig(FIG / "figE30_variability.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved figE30")

# ============================================================================
# 2. 4-way histology × cluster cross-stratification
# ============================================================================
print("\n=== 2. 4-way histology × cluster cross-stratification ===")
dm["combo"] = dm["histology_subtype"].astype(str) + "-" + dm["v17_dark_cluster"]
combo_counts = dm["combo"].value_counts()
print(f"Combo counts:\n{combo_counts}")

# PFI by combo
combo_pfi = []
for combo in dm["combo"].unique():
    if "unknown" in combo: continue
    sub = dm[dm["combo"] == combo]
    n = int(sub["PFI"].notna().sum())
    e = int(sub["PFI"].sum())
    combo_pfi.append({"combo": combo, "n": n, "events": e,
                      "rate": float(e / n * 100) if n else 0})
cpfi = pd.DataFrame(combo_pfi).sort_values("rate", ascending=False)
print(f"\nPFI rate by combo:\n{cpfi.to_string(index=False)}")
cpfi.to_csv(DATA_OUT / "histology_cluster_combo_pfi.tsv", sep="\t", index=False)

# Plot
fig, ax = plt.subplots(figsize=(10, 5))
ax.bar(cpfi["combo"], cpfi["rate"], color=["#e76f51", "#f4a261", "#4361ee", "#2a9d8f"][:len(cpfi)], alpha=0.85)
for i, (_, r) in enumerate(cpfi.iterrows()):
    ax.text(i, r["rate"] + 0.3, f"{r['events']}/{r['n']}\n({r['rate']:.1f}%)",
            ha="center", va="bottom", fontsize=9, fontweight="bold")
ax.set_ylabel("PFI event rate (%)")
ax.set_title("Figure E31 — 4-way histology × cluster PFI rate (TCGA Liu 2018 CDR)")
ax.grid(alpha=0.3, axis="y")
plt.setp(ax.get_xticklabels(), rotation=15, ha="right")
fig.tight_layout()
fig.savefig(FIG / "figE31_histology_cluster_combo.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved figE31")

# ============================================================================
# 3. Harrell C-index for 8-gene continuous score
# ============================================================================
print("\n=== 3. Harrell C-index ===")
# Score 8-gene via mean
panel_present = [g for g in PANEL if g in expr.index]
gene_score = expr.loc[panel_present].mean(axis=0)
master["panel_score"] = master["expr_col"].map(gene_score)

# Within DM, on PFI
dm_c = master[master["dm_status"] & master["PFI"].notna()].copy()
print(f"DM with PFI: n={len(dm_c)}, events={int(dm_c['PFI'].sum())}")

# C-index for cluster (binary)
dm_c["cluster_dm2"] = (dm_c["v17_dark_cluster"] == "DM2").astype(int)
dm_c["age_n"] = pd.to_numeric(dm_c["age"], errors="coerce")
dm_c_clean = dm_c.dropna(subset=["panel_score", "age_n", "PFI", "PFI.time"])
print(f"After dropna: {len(dm_c_clean)}")
c_cluster = concordance_index(dm_c_clean["PFI.time"], dm_c_clean["cluster_dm2"], dm_c_clean["PFI"])
c_score = concordance_index(dm_c_clean["PFI.time"], dm_c_clean["panel_score"], dm_c_clean["PFI"])
c_age = concordance_index(dm_c_clean["PFI.time"], dm_c_clean["age_n"], dm_c_clean["PFI"])
print(f"  C-index cluster (DM2 binary): {c_cluster:.3f}")
print(f"  C-index 8-gene continuous score: {c_score:.3f}")
print(f"  C-index age: {c_age:.3f}")

# Combined
cox = CoxPHFitter(penalizer=0.05)
sub = dm_c_clean[["panel_score", "age_n", "PFI", "PFI.time"]].rename(columns={"PFI": "event", "PFI.time": "time"})
cox.fit(sub, duration_col="time", event_col="event")
risk = cox.predict_partial_hazard(sub.drop(columns=["time", "event"]))
c_combined = concordance_index(sub["time"], -risk, sub["event"])  # negate for higher=earlier
print(f"  C-index combined (Cox panel+age): {c_combined:.3f}")

# Quartile-based PFI
master["panel_q"] = pd.qcut(master["panel_score"], q=4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
q_results = []
for q in ["Q1", "Q2", "Q3", "Q4"]:
    sub = master[(master["panel_q"] == q) & master["dm_status"]]
    n = int(sub["PFI"].notna().sum()); e = int(sub["PFI"].sum())
    q_results.append({"quartile": q, "n": n, "events": e, "rate_pct": float(e/n*100) if n else 0})
qdf = pd.DataFrame(q_results)
print(f"\nQuartile PFI rate:\n{qdf.to_string(index=False)}")

# ============================================================================
# 4. Drop-one cohort sensitivity meta-analysis
# ============================================================================
print("\n=== 4. Drop-one sensitivity for meta-analysis ===")
cohorts = [
    {"cohort": "GSE241184", "n": 2427, "r": 0.905},
    {"cohort": "GSE193581 PTC", "n": 8590, "r": 0.687},
    {"cohort": "GSE193581 PTC+ATC", "n": 14624, "r": 0.893},
    {"cohort": "GSE184362", "n": 21821, "r": 0.889},
    {"cohort": "GSE184362 multi-site", "n": 13005, "r": 0.914},
]
mdf = pd.DataFrame(cohorts)
mdf["z"] = np.arctanh(mdf["r"])
mdf["w"] = mdf["n"] - 3

drop_results = []
for i in range(len(mdf)):
    leave_out = mdf.iloc[i]["cohort"]
    rest = mdf.drop(i)
    pooled_z = (rest["z"] * rest["w"]).sum() / rest["w"].sum()
    pooled_r = np.tanh(pooled_z)
    se = 1 / np.sqrt(rest["w"].sum())
    ci_lo = np.tanh(pooled_z - 1.96 * se); ci_hi = np.tanh(pooled_z + 1.96 * se)
    drop_results.append({"dropped": leave_out, "pooled_r": float(pooled_r),
                         "ci_lo": float(ci_lo), "ci_hi": float(ci_hi)})
ddf = pd.DataFrame(drop_results)
# All-in baseline
pooled_z_all = (mdf["z"] * mdf["w"]).sum() / mdf["w"].sum()
pooled_r_all = np.tanh(pooled_z_all)
se_all = 1 / np.sqrt(mdf["w"].sum())
ci_all = (np.tanh(pooled_z_all - 1.96 * se_all), np.tanh(pooled_z_all + 1.96 * se_all))
print(f"All cohorts pooled r = {pooled_r_all:.3f} [{ci_all[0]:.3f}-{ci_all[1]:.3f}]")
print(f"\nDrop-one sensitivity:\n{ddf.round(4).to_string(index=False)}")
ddf.to_csv(DATA_OUT / "drop_one_sensitivity.tsv", sep="\t", index=False)

# Plot
fig, ax = plt.subplots(figsize=(11, 5))
y = np.arange(len(ddf))[::-1]
ax.errorbar(ddf["pooled_r"], y, xerr=[ddf["pooled_r"]-ddf["ci_lo"], ddf["ci_hi"]-ddf["pooled_r"]],
            fmt="s", color="#1864ab", markersize=10, capsize=4)
for i, (_, r) in enumerate(ddf.iterrows()):
    ax.text(r["ci_hi"] + 0.005, y[i], f"  {r['dropped']} 제외 시 r={r['pooled_r']:.3f}", va="center", fontsize=9)
ax.axvline(pooled_r_all, color="red", linestyle="--", lw=1.5, label=f"All-in r={pooled_r_all:.3f}")
ax.fill_betweenx([min(y)-0.5, max(y)+0.5], ci_all[0], ci_all[1], color="red", alpha=0.1, label="All-in 95% CI")
ax.set_yticks(list(y))
ax.set_yticklabels([f"Drop {d}" for d in ddf["dropped"]])
ax.set_xlabel("Pooled r (5−1=4 cohorts)")
ax.set_title("Figure E32 — Drop-one cohort sensitivity meta-analysis")
ax.set_xlim(0.7, 1.0)
ax.legend(loc="lower left")
ax.grid(alpha=0.3, axis="x")
fig.tight_layout()
fig.savefig(FIG / "figE32_drop_one_sensitivity.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved figE32")

# ============================================================================
# 5. Power analysis for HR detection
# ============================================================================
print("\n=== 5. Power analysis ===")
# How many events needed for 80% power to detect HR=2.0, HR=1.5, HR=3.0?
# Schoenfeld formula: events = (z_alpha + z_beta)^2 / (log(HR)^2 * p1 * (1-p1))
# where p1 = proportion in DM2 group
p1 = 0.4  # 40% DM2 in DM cohort
z_alpha = 1.96  # two-sided alpha 0.05
z_beta_80 = 0.84  # 80% power
z_beta_90 = 1.28

power_data = []
for hr in [1.3, 1.5, 1.8, 2.0, 2.5, 3.0, 4.0]:
    events_80 = (z_alpha + z_beta_80) ** 2 / (np.log(hr) ** 2 * p1 * (1 - p1))
    events_90 = (z_alpha + z_beta_90) ** 2 / (np.log(hr) ** 2 * p1 * (1 - p1))
    power_data.append({"target_HR": hr, "events_for_80pct_power": int(np.ceil(events_80)),
                       "events_for_90pct_power": int(np.ceil(events_90))})
pdf = pd.DataFrame(power_data)
print(pdf.to_string(index=False))
pdf.to_csv(DATA_OUT / "power_analysis.tsv", sep="\t", index=False)

# Plot
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(pdf["target_HR"], pdf["events_for_80pct_power"], "o-", lw=2, markersize=10, label="80% power", color="#4361ee")
ax.plot(pdf["target_HR"], pdf["events_for_90pct_power"], "s--", lw=2, markersize=10, label="90% power", color="#e63946")
ax.axhline(9, color="orange", lw=2, alpha=0.6, label="Current TCGA-THCA DM PFI events (n=9)")
ax.set_xlabel("Target Hazard Ratio (DM2/DM1)")
ax.set_ylabel("Events needed (Schoenfeld formula)")
ax.set_yscale("log")
ax.set_title(f"Figure E33 — Power analysis for HR detection within DM\n(p1=DM2 fraction=0.4, two-sided α=0.05)")
ax.legend()
ax.grid(alpha=0.3, which="both")
fig.tight_layout()
fig.savefig(FIG / "figE33_power_analysis.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved figE33")

# Save final
final = {
    "brown_forsythe_per_gene": bf_df.to_dict(orient="records"),
    "brown_forsythe_per_pathway": pw_var_df.to_dict(orient="records"),
    "histology_cluster_combo_pfi": cpfi.to_dict(orient="records"),
    "harrell_c": {"cluster": float(c_cluster), "panel_continuous": float(c_score), "age": float(c_age),
                  "combined_panel_age": float(c_combined)},
    "quartile_pfi": qdf.to_dict(orient="records"),
    "all_in_pooled_r": float(pooled_r_all), "all_in_CI": list(ci_all),
    "drop_one_sensitivity": ddf.to_dict(orient="records"),
    "power_analysis": pdf.to_dict(orient="records"),
}
(DATA_OUT / "p10_enrich_summary.json").write_text(json.dumps(final, indent=2, default=str))
print(f"\n=== ALL P10 ENRICH DONE ===")
