#!/usr/bin/env python3
"""Paper 9 SL first-pass — local-data-only, honest scope.

Inputs (all local, no download):
  - CCLE thyroid: project/results/v8p1_rigor/f_ccle_validation/{ccle_thyroid_expression.tsv, ccle_brs_labels.tsv, ccle_thyroid_metadata.tsv}
  - DM1 vs DM2 DGE: project/results/dark_matter_phase2/web/data/dge_dm1_vs_dm2.tsv
  - TCGA-THCA scored: project_external_st/results/extra/s_tcga_thca_scored.tsv  (DM1_like + os_event + os_days + tert_pos)
  - TCGA pancan expression: project/data/raw/TCGA_pancan/pancan_geneExp.gz  (HiSeqV2 log2(norm_count+1))
  - TCGA pancan phenotype: project/data/raw/TCGA_pancan/phenotype.tsv.gz

Outputs:
  project/results/paper9_sl_first_pass/
    candidate_targets.tsv             — Paper 9 candidate gene panel + class
    ccle_dm1_target_corr.tsv          — CCLE n=13 DM1-score × target-gene Pearson r
    dge_candidate_ranking.tsv         — DM1-vs-DM2 effect sizes for candidate genes
    tcga_thca_target_cox.tsv          — TCGA-THCA Cox HR per target gene (univariate, OS)
    tcga_thca_dm1_target_cox.tsv      — TCGA-THCA Cox HR for target gene WITHIN DM1-high subset
    figures/
      F2_ccle_dm1_target_corr.png
      F3_dge_candidate_ranking.png
      F6_tcga_target_survival.png
"""
from __future__ import annotations

import gzip
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/paper9_sl_first_pass"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

# Paper 9 candidate target panel — class-tagged
CANDIDATES = {
    "NAMPT":   "NAD_salvage",
    "NAPRT":   "NAD_de_novo_SR_partner",
    "JAK1":    "JAK_STAT",
    "JAK2":    "JAK_STAT",
    "TYK2":    "JAK_STAT",
    "STAT3":   "JAK_STAT",
    "IL6R":    "JAK_STAT_upstream",
    "OSMR":    "JAK_STAT_upstream",
    "DNMT1":   "Epigenetic",
    "DNMT3A":  "Epigenetic",
    "DNMT3B":  "Epigenetic",
    "HDAC1":   "Epigenetic",
    "HDAC2":   "Epigenetic",
    "HDAC6":   "Epigenetic",
    "KDM1A":   "Epigenetic",
    "EZH2":    "Epigenetic",
    "LYN":     "SFK",
    "FYN":     "SFK",
    "SRC":     "SFK",
    "YES1":    "SFK",
    "KCNN4":   "Ion_channel",
    "PARP1":   "DDR",
    "PARP2":   "DDR",
    "ATR":     "DDR",
    "CHEK1":   "DDR",
    "CHEK2":   "DDR",
    "GLS":     "Metabolism",
    "LDHA":    "Metabolism_glycolysis",
    "HK2":     "Metabolism_glycolysis",
    "IDH2":    "Metabolism",
    "MYC":     "Metabolism_master",
    "TACSTD2": "TROP2_ADC_non_SL",
    # TF collapse anchor (for cross-reference)
    "FOXE1":   "Lineage_TF_anchor",
    "NKX2-1":  "Lineage_TF_anchor",
    "PAX8":    "Lineage_TF_anchor",
    "HHEX":    "Lineage_TF_anchor",
}
TF_COLLAPSE = ["FOXE1", "NKX2-1", "PAX8", "HHEX"]
RAI_8 = ["TPO", "DIO1", "TSHR", "PAX8", "TG", "FOXE1", "NKX2-1", "SLC5A5"]

candidates_df = pd.DataFrame(
    [{"gene": g, "class": c} for g, c in CANDIDATES.items()]
)
candidates_df.to_csv(OUT / "candidate_targets.tsv", sep="\t", index=False)
print(f"Wrote candidate_targets.tsv  ({len(candidates_df)} genes)")


# =====================================================================
# B2 — CCLE n=13 DM1 score + target-gene correlation
# =====================================================================
print("\n=== B2: CCLE thyroid DM1 score + target correlation ===")
ccle_dir = ROOT / "project/results/v8p1_rigor/f_ccle_validation"
expr = pd.read_csv(ccle_dir / "ccle_thyroid_expression.tsv", sep="\t")
expr = expr.set_index("entrezGeneId")  # Despite the name, this column holds gene SYMBOLS
brs = pd.read_csv(ccle_dir / "ccle_brs_labels.tsv", sep="\t")
print(f"CCLE thyroid: {expr.shape[1]} lines × {expr.shape[0]} genes")

# Within-sample z per gene then mean across the RAI_8 set → DM1_like = -RAI_8
# Within-sample z = (gene_expr_in_line - mean_over_genes_in_line) / std_over_genes_in_line.
# This is the same logic as score_8gene_dm1.py.
def within_sample_z_score(expr_df: pd.DataFrame, geneset: list[str]) -> pd.Series:
    """Mean of within-sample z-scores across `geneset` genes, per cell line."""
    avail = [g for g in geneset if g in expr_df.index]
    missing = [g for g in geneset if g not in expr_df.index]
    if missing:
        print(f"  missing in CCLE: {missing}")
    sub = expr_df.loc[avail]
    # Within-sample z: per-line, z-score over the line's full transcriptome.
    line_mu = expr_df.mean(axis=0)
    line_sd = expr_df.std(axis=0).replace(0, np.nan)
    z = (sub.sub(line_mu, axis=1)).div(line_sd, axis=1)
    return z.mean(axis=0)


rai_score = within_sample_z_score(expr, RAI_8)
tf_collapse_score = within_sample_z_score(expr, TF_COLLAPSE)
dm1_score = -rai_score  # DM1_like = -RAI_8 (lineage-collapsed = high DM1)

ccle_state = pd.DataFrame({
    "line": dm1_score.index,
    "RAI_8_score_z": rai_score.values,
    "TF_collapse_score_z": tf_collapse_score.values,
    "DM1_like_score": dm1_score.values,
}).merge(
    brs[["sample", "name", "cancer_type_detailed", "oncotree", "hist_subtype",
         "mutation_label"]],
    left_on="line", right_on="sample", how="left",
)
ccle_state["DM1_high"] = (ccle_state["DM1_like_score"] >
                          ccle_state["DM1_like_score"].median()).astype(int)
ccle_state.to_csv(OUT / "ccle_dm1_state.tsv", sep="\t", index=False)
print(f"CCLE DM1 state: {ccle_state['DM1_high'].sum()} DM1-high / "
      f"{(1-ccle_state['DM1_high']).sum()} DM1-low (median split, n={len(ccle_state)})")
print(ccle_state[["line", "name", "oncotree", "DM1_like_score",
                  "TF_collapse_score_z", "DM1_high"]].sort_values(
                      "DM1_like_score", ascending=False).to_string(index=False))

# Target-gene correlation
corr_rows = []
for g, cls in CANDIDATES.items():
    if g not in expr.index:
        corr_rows.append({"gene": g, "class": cls, "in_ccle": False,
                          "pearson_r": np.nan, "p": np.nan,
                          "mean_dm1_high": np.nan, "mean_dm1_low": np.nan,
                          "delta": np.nan})
        continue
    g_expr = expr.loc[g]
    r, p = stats.pearsonr(dm1_score.values, g_expr.values)
    high = ccle_state["DM1_high"] == 1
    mh = g_expr[ccle_state.loc[high, "line"]].mean()
    ml = g_expr[ccle_state.loc[~high, "line"]].mean()
    corr_rows.append({
        "gene": g, "class": cls, "in_ccle": True,
        "pearson_r": r, "p": p,
        "mean_dm1_high": mh, "mean_dm1_low": ml,
        "delta": mh - ml,
    })
corr_df = pd.DataFrame(corr_rows).sort_values("pearson_r",
                                              ascending=False, na_position="last")
corr_df.to_csv(OUT / "ccle_dm1_target_corr.tsv", sep="\t", index=False)
print(f"\nWrote ccle_dm1_target_corr.tsv  ({len(corr_df)} rows)")
print("\nTop +5 / Bottom -5 (Pearson r vs DM1_like_score):")
print(corr_df.head(5)[["gene", "class", "pearson_r", "p", "delta"]].to_string(index=False))
print(corr_df.tail(5)[["gene", "class", "pearson_r", "p", "delta"]].to_string(index=False))


# =====================================================================
# B3 — DGE-based candidate ranking (DM1 vs DM2 in TCGA-THCA bulk)
# =====================================================================
print("\n=== B3: DM1 vs DM2 DGE-based candidate ranking ===")
dge = pd.read_csv(ROOT / "project/results/dark_matter_phase2/web/data/dge_dm1_vs_dm2.tsv", sep="\t")
dge_cand = dge[dge["gene"].isin(CANDIDATES.keys())].copy()
dge_cand["class"] = dge_cand["gene"].map(CANDIDATES)
# log2FC_DM2vDM1 < 0 means DM1 > DM2 (target up in DM1).
dge_cand["dm1_up"] = dge_cand["log2FC_DM2vDM1"] < 0
dge_cand = dge_cand.sort_values("log2FC_DM2vDM1")
dge_cand.to_csv(OUT / "dge_candidate_ranking.tsv", sep="\t", index=False)
print(f"Wrote dge_candidate_ranking.tsv  ({len(dge_cand)} genes found in DGE)")
missing = sorted(set(CANDIDATES) - set(dge_cand["gene"]))
print(f"Missing in DGE table: {missing}")
print("\nTop 8 DM1-up / Top 8 DM2-up:")
print(dge_cand.head(8)[["gene", "class", "log2FC_DM2vDM1", "t",
                         "q_BH"]].to_string(index=False))
print(dge_cand.tail(8)[["gene", "class", "log2FC_DM2vDM1", "t",
                         "q_BH"]].to_string(index=False))


# =====================================================================
# B4 — TCGA pancan extract THCA + per-target Cox HR (Wald)
# =====================================================================
print("\n=== B4: TCGA-THCA per-target gene survival ===")
scored = pd.read_csv(ROOT / "project_external_st/results/extra/s_tcga_thca_scored.tsv", sep="\t")
scored = scored.dropna(subset=["os_event", "os_days", "DM1_like"]).copy()
scored["os_event"] = scored["os_event"].astype(int)
scored["os_days"] = scored["os_days"].astype(float)
print(f"TCGA-THCA scored (after OS dropna): {len(scored)} samples; "
      f"OS event rate {scored['os_event'].mean():.3f}")
scored["DM1_high"] = (scored["DM1_like"] > scored["DM1_like"].median()).astype(int)

# Stream-parse the pancan expression matrix; rows = gene symbols, cols = sample IDs.
pancan_path = ROOT / "project/data/raw/TCGA_pancan/pancan_geneExp.gz"
print(f"Streaming {pancan_path.name} for candidate genes …")
candidate_set = set(CANDIDATES.keys())

with gzip.open(pancan_path, "rt") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    sample_cols = header[1:]
    rows = []
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        gene = parts[0]
        if gene in candidate_set:
            vals = parts[1:]
            rows.append([gene] + vals)
expr_pancan = pd.DataFrame(rows, columns=header)
expr_pancan = expr_pancan.set_index(header[0])
expr_pancan = expr_pancan.replace(["", "NA", "na", "NaN"], np.nan).astype(float)
print(f"Pancan candidate expression: {expr_pancan.shape[0]} genes × "
      f"{expr_pancan.shape[1]} samples")

# Match THCA samples — TCGA scored 'patient' is TCGA-XX-XXXX; pancan cols are TCGA-XX-XXXX-NNX.
# Build a map from patient → list of sample columns; preferentially keep tumor (-01).
pat_to_cols = {}
for c in expr_pancan.columns:
    pat = "-".join(c.split("-")[:3])
    pat_to_cols.setdefault(pat, []).append(c)

def pick_tumor(cols: list[str]) -> str | None:
    tumor = [c for c in cols if c.split("-")[3].startswith("01")]
    if tumor:
        return tumor[0]
    return cols[0] if cols else None

scored["sample_col"] = scored["patient"].map(
    lambda p: pick_tumor(pat_to_cols.get(p, []))
)
matched = scored.dropna(subset=["sample_col"])
print(f"TCGA scored matched to pancan expression: {len(matched)} / {len(scored)}")

# Per-gene univariate Cox HR via z = -ln(HR_lo/HR_hi) ... use lifelines if available.
try:
    from lifelines import CoxPHFitter
    HAS_LIFELINES = True
except ImportError:
    HAS_LIFELINES = False
    print("WARN lifelines not available — falling back to log-rank.")

cox_rows = []
cox_dm1_rows = []

dm1_high_mask = matched["DM1_high"] == 1

for g in CANDIDATES:
    if g not in expr_pancan.index:
        continue
    g_expr = expr_pancan.loc[g, matched["sample_col"].values].values.astype(float)
    df = pd.DataFrame({
        "gene_expr": g_expr,
        "os_event": matched["os_event"].values,
        "os_days": matched["os_days"].values,
        "DM1_high": matched["DM1_high"].values,
    }).dropna()
    df["os_event"] = df["os_event"].astype(int)
    df["os_days"] = df["os_days"].astype(float)
    df["DM1_high"] = df["DM1_high"].astype(int)
    if len(df) < 30:
        continue

    # Whole cohort
    if HAS_LIFELINES:
        try:
            cph = CoxPHFitter()
            cph.fit(df[["gene_expr", "os_event", "os_days"]],
                    duration_col="os_days", event_col="os_event")
            s = cph.summary.loc["gene_expr"]
            cox_rows.append({
                "gene": g, "class": CANDIDATES[g], "n": len(df),
                "hr": s["exp(coef)"], "hr_lo": s["exp(coef) lower 95%"],
                "hr_hi": s["exp(coef) upper 95%"], "p": s["p"],
            })
        except Exception as e:
            cox_rows.append({"gene": g, "class": CANDIDATES[g], "n": len(df),
                             "hr": np.nan, "hr_lo": np.nan, "hr_hi": np.nan,
                             "p": np.nan, "error": str(e)})
    # Within DM1-high
    sub = df[df["DM1_high"] == 1]
    if len(sub) >= 30 and HAS_LIFELINES:
        try:
            cph2 = CoxPHFitter()
            cph2.fit(sub[["gene_expr", "os_event", "os_days"]],
                     duration_col="os_days", event_col="os_event")
            s2 = cph2.summary.loc["gene_expr"]
            cox_dm1_rows.append({
                "gene": g, "class": CANDIDATES[g], "n": len(sub),
                "hr": s2["exp(coef)"], "hr_lo": s2["exp(coef) lower 95%"],
                "hr_hi": s2["exp(coef) upper 95%"], "p": s2["p"],
            })
        except Exception:
            pass

cox_df = pd.DataFrame(cox_rows).sort_values("p")
cox_dm1_df = pd.DataFrame(cox_dm1_rows).sort_values("p")
cox_df.to_csv(OUT / "tcga_thca_target_cox.tsv", sep="\t", index=False)
cox_dm1_df.to_csv(OUT / "tcga_thca_dm1_target_cox.tsv", sep="\t", index=False)
print(f"Wrote tcga_thca_target_cox.tsv  ({len(cox_df)} rows, whole TCGA-THCA)")
print(f"Wrote tcga_thca_dm1_target_cox.tsv  ({len(cox_dm1_df)} rows, DM1-high subset)")
if len(cox_df):
    print("\nTop 10 by p (whole cohort):")
    print(cox_df.head(10)[["gene", "class", "n", "hr", "p"]].to_string(index=False))
if len(cox_dm1_df):
    print("\nTop 10 by p (DM1-high subset):")
    print(cox_dm1_df.head(10)[["gene", "class", "n", "hr", "p"]].to_string(index=False))


# =====================================================================
# Figures
# =====================================================================
print("\n=== Figures ===")

# F2 — CCLE DM1 × target Pearson r barplot
fig, ax = plt.subplots(figsize=(8, 9))
plot_df = corr_df.dropna(subset=["pearson_r"]).copy()
plot_df = plot_df.sort_values("pearson_r")
classes = plot_df["class"].unique()
palette = {c: plt.cm.tab20(i / max(1, len(classes)))
           for i, c in enumerate(classes)}
colors = plot_df["class"].map(palette).values
ax.barh(plot_df["gene"], plot_df["pearson_r"], color=colors)
ax.axvline(0, color="black", linewidth=0.5)
ax.set_xlabel("Pearson r (target gene expr × DM1_like_score, CCLE n=13)")
ax.set_title("F2 — CCLE thyroid: DM1 × candidate-target correlation\n"
             "(positive r = target up in DM1-like)")
plt.tight_layout()
plt.savefig(FIG / "F2_ccle_dm1_target_corr.png", dpi=140)
plt.close()
print("Wrote F2_ccle_dm1_target_corr.png")

# F3 — DGE log2FC ranking (TCGA-THCA bulk, DM1 vs DM2)
fig, ax = plt.subplots(figsize=(8, 9))
plot2 = dge_cand.dropna(subset=["log2FC_DM2vDM1"]).copy()
plot2 = plot2.sort_values("log2FC_DM2vDM1")
plot2["dm1_minus_dm2"] = -plot2["log2FC_DM2vDM1"]  # so positive = DM1 up
colors2 = plot2["class"].map(palette).fillna("#888").values
ax.barh(plot2["gene"], plot2["dm1_minus_dm2"], color=colors2)
ax.axvline(0, color="black", linewidth=0.5)
ax.set_xlabel("log2(DM1 / DM2) — DGE TCGA-THCA bulk")
ax.set_title("F3 — DGE candidate ranking: DM1 vs DM2\n"
             "(positive = target up in DM1)")
# Significance markers
for i, (idx, row) in enumerate(plot2.iterrows()):
    if pd.notna(row["q_BH"]) and row["q_BH"] < 0.05:
        ax.text(row["dm1_minus_dm2"], i, "  *", va="center")
plt.tight_layout()
plt.savefig(FIG / "F3_dge_candidate_ranking.png", dpi=140)
plt.close()
print("Wrote F3_dge_candidate_ranking.png")

# F6 — TCGA-THCA per-target Cox HR forest (whole cohort + DM1-high subset)
if len(cox_df):
    fig, axes = plt.subplots(1, 2, figsize=(13, 9), sharey=True)
    for ax, df, title in zip(axes, [cox_df, cox_dm1_df],
                              ["Whole TCGA-THCA", "DM1-high subset"]):
        d = df.dropna(subset=["hr"]).copy()
        d = d.sort_values("hr")
        if len(d) == 0:
            ax.set_title(f"{title}\n(no Cox fits)")
            continue
        d["log_hr"] = np.log(d["hr"])
        d["log_hr_lo"] = np.log(d["hr_lo"])
        d["log_hr_hi"] = np.log(d["hr_hi"])
        y = np.arange(len(d))
        ax.errorbar(d["log_hr"], y,
                    xerr=[d["log_hr"] - d["log_hr_lo"],
                          d["log_hr_hi"] - d["log_hr"]],
                    fmt="o", color="#34547A", ecolor="#7A8AA8")
        ax.axvline(0, color="black", linewidth=0.5)
        ax.set_yticks(y)
        ax.set_yticklabels(d["gene"])
        ax.set_xlabel("log HR (per +1 expression unit)")
        ax.set_title(f"{title}\n(n={d['n'].iloc[0] if len(d) else 0})")
    fig.suptitle("F6 — TCGA-THCA OS Cox HR per Paper 9 candidate target", y=1.0)
    plt.tight_layout()
    plt.savefig(FIG / "F6_tcga_target_survival.png", dpi=140)
    plt.close()
    print("Wrote F6_tcga_target_survival.png")

# Summary header
summary = {
    "ccle_n": int(expr.shape[1]),
    "ccle_dm1_high_n": int(ccle_state["DM1_high"].sum()),
    "ccle_dm1_low_n": int((1 - ccle_state["DM1_high"]).sum()),
    "candidates_in_ccle": int(corr_df["in_ccle"].sum()),
    "candidates_in_dge": int(len(dge_cand)),
    "tcga_thca_n_scored": int(len(scored)),
    "tcga_thca_n_matched_pancan": int(len(matched)),
    "tcga_thca_event_rate": float(matched["os_event"].mean()),
    "cox_rows_whole": int(len(cox_df)),
    "cox_rows_dm1_high": int(len(cox_dm1_df)),
}
import json
(OUT / "summary.json").write_text(json.dumps(summary, indent=2))
print("\n=== summary ===")
print(json.dumps(summary, indent=2))
print("\nDone.")
