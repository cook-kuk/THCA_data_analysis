#!/usr/bin/env python3
"""v17 HLA analysis — does HLA expression/genotype reveal NEW biology in THCA?

H1: HLA Class I suppressed in BRAF V600E (TCGA-THCA)
H2: HLA Class I higher in DM1 vs DM2 (TCGA-THCA)
H3: HLA Class II correlates with lymphocytic infiltration (autoimmune-PTC axis)
H4: HLA decreases with dedifferentiation Normal -> PTC -> PDFP -> ATC (Korean GSE213647)
H5: HLA-A 4-digit genotype distribution by DM/BRAF (deferred — requires Synapse Thorsson)
"""
from __future__ import annotations
import json, sys, traceback
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy import stats

RNG = 42
PROJECT = Path("/opt/thyroid-dash/project")
TAB = PROJECT / "results/v17_hla"
FIG = PROJECT / "reports/html/figs_interactive/v17"
LOG = PROJECT / "logs/v17_hla.log"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

BG = "#0b0e12"; INK = "#F2F2F2"
DARK = dict(template="plotly_dark", paper_bgcolor=BG, plot_bgcolor=BG, font=dict(color=INK, size=14))

CLASS_I = ["HLA-A", "HLA-B", "HLA-C", "B2M"]
CLASS_I_REG = ["TAP1", "TAP2", "NLRC5"]
CLASS_II = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1"]
CLASS_II_REG = ["CIITA"]

ENSG_HLA = {
    "ENSG00000206503": "HLA-A",
    "ENSG00000234745": "HLA-B",
    "ENSG00000204525": "HLA-C",
    "ENSG00000166710": "B2M",
    "ENSG00000168394": "TAP1",
    "ENSG00000204267": "TAP2",
    "ENSG00000140853": "NLRC5",
    "ENSG00000204287": "HLA-DRA",
    "ENSG00000196126": "HLA-DRB1",
    "ENSG00000204252": "HLA-DPA1",
    "ENSG00000223865": "HLA-DPB1",
    "ENSG00000196735": "HLA-DQA1",
    "ENSG00000179344": "HLA-DQB1",
    "ENSG00000179583": "CIITA",
}

LOG_FH = open(LOG, "w")
def log(msg=""):
    print(msg); LOG_FH.write(str(msg) + "\n"); LOG_FH.flush()


def cohens_d(a, b):
    a = np.asarray(a); b = np.asarray(b)
    s = np.sqrt(((len(a)-1)*a.var(ddof=1) + (len(b)-1)*b.var(ddof=1)) / (len(a)+len(b)-2))
    if s == 0: return 0.0
    return float((a.mean() - b.mean()) / s)


def zscore(df):
    """Row-wise z-score (genes x samples)."""
    return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1) + 1e-9, axis=0)


# =============================================================
# 1. Load TCGA-THCA expression + metadata
# =============================================================
log("="*70)
log("v17 HLA analysis — START")
log("="*70)

expr_path = "/data/thca/data_processed/bulk_rnaseq_v3/TCGA-THCA_v3_log2.tsv"
log(f"\n[1] Loading TCGA-THCA expression: {expr_path}")
expr = pd.read_csv(expr_path, sep="\t", index_col=0)
log(f"  expr shape: {expr.shape} (genes x samples)")

hla_genes = CLASS_I + CLASS_I_REG + CLASS_II + CLASS_II_REG
present = [g for g in hla_genes if g in expr.index]
missing = [g for g in hla_genes if g not in expr.index]
log(f"  HLA genes present: {len(present)}/{len(hla_genes)}")
if missing: log(f"  missing: {missing}")

hla_expr = expr.loc[present].copy()
log(f"  hla_expr range: {hla_expr.values.min():.2f} - {hla_expr.values.max():.2f}")

# z-score within TCGA cohort
hla_z = zscore(hla_expr)

# Per-sample composite scores
ci_genes = [g for g in CLASS_I + CLASS_I_REG if g in hla_z.index]
cii_genes = [g for g in CLASS_II + CLASS_II_REG if g in hla_z.index]
score_i = hla_z.loc[ci_genes].mean(axis=0)
score_ii = hla_z.loc[cii_genes].mean(axis=0)
log(f"  Class I composite: {len(ci_genes)} genes, samples={len(score_i)}")
log(f"  Class II composite: {len(cii_genes)} genes, samples={len(score_ii)}")

# Sample IDs in TCGA: TCGA-XX-XXXX-01A. Patient ID = first 12 chars
def pid(sid): return sid[:12]

per_sample = pd.DataFrame({
    "sample_id": score_i.index,
    "patient_id": [pid(s) for s in score_i.index],
    "hla_class_I_score": score_i.values,
    "hla_class_II_score": score_ii.values,
})
log(f"  per-sample df: {per_sample.shape}")

# =============================================================
# 2. Merge metadata: DM1/DM2 + BRAF V600E + RAS + TERT
# =============================================================
log("\n[2] Merging metadata...")

# DM cluster from v17p3 (full cohort, 514 TCGA samples)
dm_path = PROJECT / "results/v17p3/tables/A2_dm_score_full_cohort.tsv"
dm = pd.read_csv(dm_path, sep="\t")
dm["patient_id"] = dm["sample_id"].str[:12]
dm_tcga = dm[dm["dataset"] == "TCGA-THCA"][["sample_id", "patient_id", "dm_like", "prob_dm1", "prob_dm2"]]
log(f"  DM clusters TCGA: {dm_tcga['dm_like'].value_counts().to_dict()}")

# BRAF/TERT external combined
ext = pd.read_csv(PROJECT / "results/v17_tert_recovery/v3/v3_external_combined_BRAF_TERT.tsv", sep="\t")
ext_tcga = ext[ext["cohort"] == "TCGA-THCA"][["patient_id", "braf_class", "tert", "histology"]]
log(f"  TCGA-THCA BRAF/TERT calls: {len(ext_tcga)}")
log(f"  BRAF distribution: {ext_tcga['braf_class'].value_counts().to_dict()}")
log(f"  TERT distribution: {ext_tcga['tert'].value_counts().to_dict()}")

# RAS from sample master
master = pd.read_csv(PROJECT / "results/v17/tables/sample_master_v17_full.tsv", sep="\t")
master["patient_id"] = master["sample_id"].str[:12]
master_tcga = master[master["dataset"] == "TCGA-THCA"][["sample_id", "patient_id", "ras_hotspot"]]
log(f"  RAS hotspot calls: {master_tcga['ras_hotspot'].value_counts().to_dict()}")

# Merge
merged = per_sample.merge(dm_tcga[["sample_id", "dm_like"]], on="sample_id", how="left")
merged = merged.merge(ext_tcga, on="patient_id", how="left")
merged = merged.merge(master_tcga[["sample_id", "ras_hotspot"]], on="sample_id", how="left")
log(f"  merged shape: {merged.shape}")
log(f"  merged with DM cluster: {merged['dm_like'].notna().sum()}")
log(f"  merged with BRAF: {merged['braf_class'].notna().sum()}")

# Save per-sample table
merged.to_csv(TAB / "tcga_thca_hla_per_sample.tsv", sep="\t", index=False)
log(f"  wrote {TAB}/tcga_thca_hla_per_sample.tsv")

# =============================================================
# 3. H1 — HLA Class I in BRAF V600E vs WT
# =============================================================
log("\n" + "="*70)
log("[H1] HLA Class I score: BRAF V600E vs non-V600E (TCGA-THCA)")
log("="*70)

m_braf = merged[merged["braf_class"].notna()].copy()
m_braf["braf_v600e"] = (m_braf["braf_class"] == "V600E").astype(int)
braf_pos = m_braf[m_braf["braf_v600e"] == 1]["hla_class_I_score"].values
braf_neg = m_braf[m_braf["braf_v600e"] == 0]["hla_class_I_score"].values
log(f"  n BRAF V600E+: {len(braf_pos)}, V600E-: {len(braf_neg)}")
log(f"  median I score: V600E+ {np.median(braf_pos):.3f}, V600E- {np.median(braf_neg):.3f}")
u_h1, p_h1 = stats.mannwhitneyu(braf_pos, braf_neg, alternative="two-sided")
d_h1 = cohens_d(braf_pos, braf_neg)
log(f"  Mann-Whitney U={u_h1:.0f}, p={p_h1:.3e}, Cohen's d={d_h1:.3f}")

# Also Class II by BRAF
braf_pos_ii = m_braf[m_braf["braf_v600e"] == 1]["hla_class_II_score"].values
braf_neg_ii = m_braf[m_braf["braf_v600e"] == 0]["hla_class_II_score"].values
u_h1ii, p_h1ii = stats.mannwhitneyu(braf_pos_ii, braf_neg_ii, alternative="two-sided")
d_h1ii = cohens_d(braf_pos_ii, braf_neg_ii)
log(f"  Class II by BRAF: medians {np.median(braf_pos_ii):.3f} vs {np.median(braf_neg_ii):.3f}, "
    f"U={u_h1ii:.0f}, p={p_h1ii:.3e}, d={d_h1ii:.3f}")

# Figure: violin BRAF V600E vs WT for Class I
fig = go.Figure()
for label, vals, color in [("BRAF V600E", braf_pos, "#c24c4c"), ("BRAF WT", braf_neg, "#7ccfcd")]:
    fig.add_trace(go.Violin(y=vals, name=f"{label} (n={len(vals)})", box_visible=True, meanline_visible=True,
                            line_color=color, fillcolor=color, opacity=0.6, points="all", jitter=0.3, marker_size=3))
fig.update_layout(title=dict(text=f"<b>TCGA-THCA — HLA Class I composite z-score by BRAF V600E status</b><br>" +
                                  f"<sub style='color:#7ccfcd'>HLA-A,B,C + B2M + TAP1/2 + NLRC5 mean z. " +
                                  f"Mann-Whitney p={p_h1:.2e}, Cohen's d={d_h1:.3f}</sub>",
                             font=dict(size=14, color=INK)),
                  yaxis=dict(title=dict(text="HLA Class I composite z", font=dict(color=INK)),
                             gridcolor="rgba(255,255,255,0.08)", tickfont=dict(color=INK), zeroline=True),
                  xaxis=dict(tickfont=dict(color=INK)),
                  showlegend=False, height=480, margin=dict(l=80, r=20, t=110, b=60), **DARK)
fig.write_html(FIG / "v17_hla_braf_split.html", include_plotlyjs="cdn", full_html=True)
log(f"  wrote {FIG}/v17_hla_braf_split.html")

# =============================================================
# 4. H2 — HLA scores in DM1 vs DM2
# =============================================================
log("\n" + "="*70)
log("[H2] HLA Class I & II: DM1 vs DM2 (TCGA-THCA)")
log("="*70)

m_dm = merged[merged["dm_like"].notna()].copy()
log(f"  n with DM label: {len(m_dm)}, distribution: {m_dm['dm_like'].value_counts().to_dict()}")
dm1_i = m_dm[m_dm["dm_like"] == "DM1_like"]["hla_class_I_score"].values
dm2_i = m_dm[m_dm["dm_like"] == "DM2_like"]["hla_class_I_score"].values
u_h2, p_h2 = stats.mannwhitneyu(dm1_i, dm2_i, alternative="two-sided")
d_h2 = cohens_d(dm1_i, dm2_i)
log(f"  Class I: DM1 median={np.median(dm1_i):.3f} (n={len(dm1_i)}) vs DM2 median={np.median(dm2_i):.3f} (n={len(dm2_i)})")
log(f"  Mann-Whitney U={u_h2:.0f}, p={p_h2:.3e}, Cohen's d={d_h2:.3f}")

dm1_ii = m_dm[m_dm["dm_like"] == "DM1_like"]["hla_class_II_score"].values
dm2_ii = m_dm[m_dm["dm_like"] == "DM2_like"]["hla_class_II_score"].values
u_h2b, p_h2b = stats.mannwhitneyu(dm1_ii, dm2_ii, alternative="two-sided")
d_h2b = cohens_d(dm1_ii, dm2_ii)
log(f"  Class II: DM1 median={np.median(dm1_ii):.3f} vs DM2 median={np.median(dm2_ii):.3f}")
log(f"  Mann-Whitney U={u_h2b:.0f}, p={p_h2b:.3e}, Cohen's d={d_h2b:.3f}")

# Correlation between HLA-I and DM probabilities (orthogonality check)
dm_probs = dm[dm["dataset"] == "TCGA-THCA"][["sample_id", "prob_dm1"]].merge(
    per_sample[["sample_id", "hla_class_I_score", "hla_class_II_score"]], on="sample_id")
rho_i, p_rho_i = stats.spearmanr(dm_probs["prob_dm1"], dm_probs["hla_class_I_score"])
rho_ii, p_rho_ii = stats.spearmanr(dm_probs["prob_dm1"], dm_probs["hla_class_II_score"])
log(f"  Spearman prob_dm1 vs HLA-I: rho={rho_i:.3f}, p={p_rho_i:.2e}")
log(f"  Spearman prob_dm1 vs HLA-II: rho={rho_ii:.3f}, p={p_rho_ii:.2e}")

# Figure: violin DM1 vs DM2 for Class I + Class II
fig = go.Figure()
colors = {"DM1_like": "#7ccfcd", "DM2_like": "#c24c4c"}
for cls, score_label, vals_dm1, vals_dm2 in [
    ("Class I", "I", dm1_i, dm2_i),
    ("Class II", "II", dm1_ii, dm2_ii),
]:
    for grp, vals, dx in [("DM1", vals_dm1, -0.2), ("DM2", vals_dm2, 0.2)]:
        col = colors["DM1_like" if grp == "DM1" else "DM2_like"]
        fig.add_trace(go.Violin(
            y=vals, x=[cls]*len(vals), name=f"{grp} {cls} (n={len(vals)})",
            side="negative" if grp == "DM1" else "positive",
            line_color=col, fillcolor=col, opacity=0.6,
            box_visible=True, meanline_visible=True, points="all", jitter=0.3, marker_size=3,
            legendgroup=grp, scalegroup=cls))
fig.update_layout(title=dict(
    text=f"<b>TCGA-THCA — HLA composite z by Dark-Matter cluster</b><br>" +
         f"<sub style='color:#7ccfcd'>Class I p={p_h2:.2e} (d={d_h2:.2f}); Class II p={p_h2b:.2e} (d={d_h2b:.2f}). " +
         f"Spearman prob_DM1 vs HLA-I rho={rho_i:.2f}, vs HLA-II rho={rho_ii:.2f}</sub>",
    font=dict(size=14, color=INK)),
    yaxis=dict(title=dict(text="HLA composite z", font=dict(color=INK)),
               gridcolor="rgba(255,255,255,0.08)", tickfont=dict(color=INK), zeroline=True),
    xaxis=dict(tickfont=dict(color=INK)),
    violinmode="overlay", height=520, margin=dict(l=80, r=20, t=110, b=60), **DARK)
fig.write_html(FIG / "v17_hla_dm1dm2.html", include_plotlyjs="cdn", full_html=True)
log(f"  wrote {FIG}/v17_hla_dm1dm2.html")

# =============================================================
# 5. H3 — HLA Class II vs lymphocyte signature (proxy for autoimmune/Hashimoto axis)
# =============================================================
log("\n" + "="*70)
log("[H3] HLA Class II vs lymphocyte/CD8 T-cell signature score")
log("="*70)

# Build lymphocyte signature: classic CD8/T-cell markers
LYMPHO_GENES = ["CD8A", "CD8B", "CD3D", "CD3E", "CD3G", "GZMA", "GZMB", "PRF1", "IFNG", "CXCL9", "CXCL10"]
lymph_present = [g for g in LYMPHO_GENES if g in expr.index]
log(f"  lymphocyte sig genes present: {len(lymph_present)}/{len(LYMPHO_GENES)}")
lymph_z = zscore(expr.loc[lymph_present])
lymph_score = lymph_z.mean(axis=0)

merged_lymph = per_sample.copy()
merged_lymph["lymph_score"] = lymph_score.reindex(merged_lymph["sample_id"]).values

# Spearman: HLA-II vs lymphocyte score
rho_h3, p_h3 = stats.spearmanr(merged_lymph["hla_class_II_score"], merged_lymph["lymph_score"])
rho_h3i, p_h3i = stats.spearmanr(merged_lymph["hla_class_I_score"], merged_lymph["lymph_score"])
log(f"  Spearman HLA-II vs lymph score: rho={rho_h3:.3f}, p={p_h3:.2e}")
log(f"  Spearman HLA-I  vs lymph score: rho={rho_h3i:.3f}, p={p_h3i:.2e}")

# Define Hashimoto-like top quartile of lymph + DM1
lymph_q3 = merged_lymph["lymph_score"].quantile(0.75)
hashi_like = merged_lymph[merged_lymph["lymph_score"] >= lymph_q3].copy()
log(f"  'Hashimoto-like' (lymph top quartile): n={len(hashi_like)}, "
    f"median HLA-II={hashi_like['hla_class_II_score'].median():.3f}, "
    f"median HLA-I={hashi_like['hla_class_I_score'].median():.3f}")

# Save lymph + HLA crosswalk
merged_lymph.to_csv(TAB / "tcga_thca_hla_lymph.tsv", sep="\t", index=False)

# =============================================================
# 6. H4 — Korean GSE213647: HLA decreases with dedifferentiation
# =============================================================
log("\n" + "="*70)
log("[H4] HLA Class I & II by histology (Korean GSE213647)")
log("="*70)

clin = pd.read_csv(PROJECT / "results/v17_korean/GSE213647_panel_score.tsv", sep="\t")
clin = clin.set_index("gsm")
log(f"  clin: {clin.shape}, histology: {clin['histology'].value_counts().to_dict()}")

kexpr = pd.read_csv("/data/thca/v17_korean/GSE213647/expression_counts.tsv.gz", sep="\t", index_col=0)
kexpr.index = kexpr.index.str.split(".").str[0]
log(f"  Korean expr: {kexpr.shape}")

ensg_present = [e for e in ENSG_HLA if e in kexpr.index]
log(f"  HLA ENSG present: {len(ensg_present)}/{len(ENSG_HLA)}")
khla = kexpr.loc[ensg_present].copy()
khla.index = [ENSG_HLA[e] for e in khla.index]

# CPM normalize using full library size, then log2(CPM+1)
lib = kexpr.sum(axis=0)
khla_cpm = khla.div(lib, axis=1) * 1e6
khla_log = np.log2(khla_cpm + 1)
log(f"  log2-CPM range: {khla_log.values.min():.2f} - {khla_log.values.max():.2f}")
khla_z = zscore(khla_log)

ki_genes = [g for g in CLASS_I + CLASS_I_REG if g in khla_z.index]
kii_genes = [g for g in CLASS_II + CLASS_II_REG if g in khla_z.index]
k_score_i = khla_z.loc[ki_genes].mean(axis=0)
k_score_ii = khla_z.loc[kii_genes].mean(axis=0)

clin["hla_class_I_score"] = k_score_i.reindex(clin.index)
clin["hla_class_II_score"] = k_score_ii.reindex(clin.index)

order = ["Normal", "PTC", "PDFP", "UTC/ATC"]
clin_ord = clin[clin["histology"].isin(order)].copy()
clin_ord["histology"] = pd.Categorical(clin_ord["histology"], categories=order, ordered=True)

groups_i = [clin_ord[clin_ord["histology"] == h]["hla_class_I_score"].dropna().values for h in order]
groups_ii = [clin_ord[clin_ord["histology"] == h]["hla_class_II_score"].dropna().values for h in order]
H_i, p_kw_i = stats.kruskal(*groups_i)
H_ii, p_kw_ii = stats.kruskal(*groups_ii)
log(f"  Class I  Kruskal-Wallis H={H_i:.2f}, p={p_kw_i:.2e}")
log(f"  Class II Kruskal-Wallis H={H_ii:.2f}, p={p_kw_ii:.2e}")
for h, gi, gii in zip(order, groups_i, groups_ii):
    log(f"    {h}: n={len(gi)}, median I={np.median(gi):.3f}, median II={np.median(gii):.3f}")

# Pairwise vs Normal
for i, h in enumerate(order):
    if i == 0: continue
    u_i, p_i = stats.mannwhitneyu(groups_i[0], groups_i[i], alternative="two-sided")
    u_ii, p_ii = stats.mannwhitneyu(groups_ii[0], groups_ii[i], alternative="two-sided")
    log(f"    Normal vs {h}: I p={p_i:.2e}, II p={p_ii:.2e}")

clin_ord.to_csv(TAB / "korean_GSE213647_hla_per_sample.tsv", sep="\t")
log(f"  wrote {TAB}/korean_GSE213647_hla_per_sample.tsv")

# Figure: box HLA-I and HLA-II by histology
fig = go.Figure()
hcolors = {"Normal": "#2ECC71", "PTC": "#F5A623", "PDFP": "#E67E22", "UTC/ATC": "#c24c4c"}
for cls, key, p_kw in [("Class I", "hla_class_I_score", p_kw_i), ("Class II", "hla_class_II_score", p_kw_ii)]:
    for h in order:
        sub = clin_ord[clin_ord["histology"] == h]
        fig.add_trace(go.Box(
            y=sub[key], x=[f"{h}<br>{cls}"]*len(sub), name=f"{h} {cls} (n={len(sub)})",
            marker=dict(color=hcolors[h]), boxpoints="outliers", line=dict(width=2),
            showlegend=False))
fig.update_layout(title=dict(
    text=f"<b>GSE213647 (Korean, Lee 2024) — HLA composite z by histology</b><br>" +
         f"<sub style='color:#7ccfcd'>Class I Kruskal-Wallis p={p_kw_i:.2e}; Class II p={p_kw_ii:.2e}. " +
         f"Order: Normal -> PTC -> PDFP -> UTC/ATC</sub>",
    font=dict(size=14, color=INK)),
    yaxis=dict(title=dict(text="HLA composite z", font=dict(color=INK)),
               gridcolor="rgba(255,255,255,0.08)", tickfont=dict(color=INK), zeroline=True),
    xaxis=dict(tickfont=dict(color=INK, size=10)),
    height=520, margin=dict(l=80, r=20, t=110, b=80), **DARK)
fig.write_html(FIG / "v17_hla_korean_histology.html", include_plotlyjs="cdn", full_html=True)
log(f"  wrote {FIG}/v17_hla_korean_histology.html")

# =============================================================
# 7. H5 — HLA 4-digit genotype distribution (Thorsson 2018)
# =============================================================
log("\n" + "="*70)
log("[H5] HLA 4-digit genotype (Thorsson 2018 PanCanAtlas)")
log("="*70)
log("  STATUS: Thorsson immune subtypes / OptiType HLA calls require")
log("  PanCanAtlas Synapse access (syn4602499). Public cBioPortal subset")
log("  (thca_tcga_pan_can_atlas_2018) does not export 4-digit HLA-A/B/C.")
log("  Skipping — figure v17_hla_genotype.html NOT generated.")
log("  Future: download Thorsson 2018 Cell Rep TCGA.HLA.4digit.byPatient.txt")
log("  via authenticated Synapse client and re-run Fisher exact for")
log("  HLA-A*02:01 vs DM1/DM2 vs BRAF V600E.")
h5_status = "deferred_thorsson_synapse_required"

# =============================================================
# 8. Summary JSON
# =============================================================
log("\n" + "="*70)
log("Writing summary JSON")
log("="*70)

summary = {
    "random_state": RNG,
    "tcga_thca": {
        "n_samples": int(len(per_sample)),
        "n_with_dm_label": int(merged["dm_like"].notna().sum()),
        "n_with_braf_call": int(merged["braf_class"].notna().sum()),
        "hla_genes_used_class_I": ci_genes,
        "hla_genes_used_class_II": cii_genes,
        "hla_genes_missing": missing,
    },
    "H1_braf_v600e_class_I": {
        "hypothesis": "BRAF V600E suppresses HLA Class I (immune evasion)",
        "n_braf_pos": int(len(braf_pos)),
        "n_braf_neg": int(len(braf_neg)),
        "median_braf_pos": float(np.median(braf_pos)),
        "median_braf_neg": float(np.median(braf_neg)),
        "mannwhitney_p": float(p_h1),
        "cohens_d": float(d_h1),
        "verdict": "supported" if (p_h1 < 0.05 and np.median(braf_pos) < np.median(braf_neg)) else "not_supported",
    },
    "H1b_braf_v600e_class_II": {
        "median_braf_pos": float(np.median(braf_pos_ii)),
        "median_braf_neg": float(np.median(braf_neg_ii)),
        "mannwhitney_p": float(p_h1ii),
        "cohens_d": float(d_h1ii),
    },
    "H2_dm1_vs_dm2": {
        "hypothesis": "DM1 (immune-hot) > DM2 (cold) for HLA",
        "class_I": {
            "median_dm1": float(np.median(dm1_i)),
            "median_dm2": float(np.median(dm2_i)),
            "mannwhitney_p": float(p_h2),
            "cohens_d": float(d_h2),
        },
        "class_II": {
            "median_dm1": float(np.median(dm1_ii)),
            "median_dm2": float(np.median(dm2_ii)),
            "mannwhitney_p": float(p_h2b),
            "cohens_d": float(d_h2b),
        },
        "spearman_prob_dm1_vs_hla_I":  {"rho": float(rho_i),  "p": float(p_rho_i)},
        "spearman_prob_dm1_vs_hla_II": {"rho": float(rho_ii), "p": float(p_rho_ii)},
    },
    "H3_class_II_vs_lymph": {
        "hypothesis": "HLA Class II tracks lymphocyte infiltration (autoimmune-PTC axis)",
        "spearman_class_II_vs_lymph": {"rho": float(rho_h3),  "p": float(p_h3)},
        "spearman_class_I_vs_lymph":  {"rho": float(rho_h3i), "p": float(p_h3i)},
    },
    "H4_korean_dedifferentiation": {
        "hypothesis": "HLA decreases Normal -> PTC -> PDFP -> ATC",
        "n_by_histology": {h: int(len(g)) for h, g in zip(order, groups_i)},
        "median_class_I": {h: float(np.median(g)) for h, g in zip(order, groups_i)},
        "median_class_II": {h: float(np.median(g)) for h, g in zip(order, groups_ii)},
        "kruskal_p_class_I": float(p_kw_i),
        "kruskal_p_class_II": float(p_kw_ii),
    },
    "H5_genotype": {
        "status": h5_status,
        "note": "Thorsson 2018 OptiType HLA-A/B/C 4-digit requires Synapse syn4602499; deferred.",
    },
    "figures": {
        "h1_braf_split": "reports/html/figs_interactive/v17/v17_hla_braf_split.html",
        "h2_dm1dm2": "reports/html/figs_interactive/v17/v17_hla_dm1dm2.html",
        "h4_korean_histology": "reports/html/figs_interactive/v17/v17_hla_korean_histology.html",
    },
}
(TAB / "v17_hla_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
log(f"  wrote {TAB}/v17_hla_summary.json")

log("\n" + "="*70)
log("v17 HLA analysis — DONE")
log("="*70)
LOG_FH.close()
