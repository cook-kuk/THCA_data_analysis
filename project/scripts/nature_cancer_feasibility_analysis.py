#!/usr/bin/env python3
from pathlib import Path
import json
import math

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "project/results/nature_cancer_feasibility"
OUT.mkdir(parents=True, exist_ok=True)


def cohen_d(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return np.nan
    pooled = math.sqrt(((nx - 1) * x.var(ddof=1) + (ny - 1) * y.var(ddof=1)) / (nx + ny - 2))
    if pooled == 0:
        return np.nan
    return (x.mean() - y.mean()) / pooled


def auc_binary(y, score):
    y = np.asarray(y, dtype=int)
    score = np.asarray(score, dtype=float)
    ok = np.isfinite(score)
    y = y[ok]
    score = score[ok]
    pos = score[y == 1]
    neg = score[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    u = stats.mannwhitneyu(pos, neg, alternative="two-sided").statistic
    return float(u / (len(pos) * len(neg)))


summary = []

# 1) External expression validation
ext_meta = pd.read_csv(ROOT / "project/results/p_external_expression_validation/sample_metadata.tsv", sep="\t")
ext_tests = pd.read_csv(ROOT / "project/results/p_external_expression_validation/external_score_tests.tsv", sep="\t")
ext_spear = pd.read_csv(ROOT / "project/results/p_external_expression_validation/external_spearman.tsv", sep="\t")

lineage_axes = [
    ("RAI_8_score", -1),
    ("DM1_like_score", 1),
    ("THYROID_NONOVERLAP_score", -1),
    ("TDS_like_score", -1),
    ("TF_collapse_score", -1),
    ("STAT3_AP1_DNMT_score", 1),
]
headline = ext_tests[ext_tests["contrast"].isin(["ATC_vs_PTC", "ATC_vs_normal", "PDTC_vs_PTC", "advanced_vs_DTC", "advanced_vs_normal"])]
direction_cells = []
for axis, sign in lineage_axes:
    col = f"{axis}_d"
    if col not in headline:
        continue
    vals = headline[["dataset", "contrast", col]].dropna()
    for _, row in vals.iterrows():
        d = row[col]
        direction_cells.append({
            "axis": axis,
            "dataset": row["dataset"],
            "contrast": row["contrast"],
            "cohens_d": d,
            "expected_sign": sign,
            "direction_consistent": np.sign(d) == sign,
        })
direction_df = pd.DataFrame(direction_cells)
direction_df.to_csv(OUT / "external_direction_cells.tsv", sep="\t", index=False)
dm1_non = ext_spear[(ext_spear["x"] == "DM1_like_score") & (ext_spear["y"] == "THYROID_NONOVERLAP_score")]
summary.append({
    "layer": "External expression replication",
    "n": int(len(ext_meta)),
    "datasets": int(ext_meta["dataset"].nunique()),
    "best_number": f"{int(direction_df['direction_consistent'].sum())}/{len(direction_df)} direction-consistent cells",
    "support": "Strong",
    "gap_closed": "Replicates lineage-state direction across independent processed expression cohorts.",
    "remaining_gap": "Still expression-only; no clinical-endpoint or prospective validation.",
})
summary.append({
    "layer": "Zero-overlap anti-circularity",
    "n": int(dm1_non["n"].sum()),
    "datasets": int(dm1_non["dataset"].nunique()),
    "best_number": f"rho range {dm1_non['rho'].min():.2f} to {dm1_non['rho'].max():.2f}",
    "support": "Strong",
    "gap_closed": "Shows DM1_like is not just the original RAI_8 gene panel recapitulating itself.",
    "remaining_gap": "Still correlative module agreement, not causality.",
})

# 2) RAI-labeled clinical-endpoint proxy dataset
rai_path = ROOT / "project/results/v17_ultimate/U1B_RAI_clinical_validation.tsv"
rai = pd.read_csv(rai_path, sep="\t")
for label, sub in [("all", rai), ("primary_only", rai[rai["is_primary_tumor"].astype(bool)])]:
    avid = sub[sub["y"] == 0]["p_DM1_score"]
    refr = sub[sub["y"] == 1]["p_DM1_score"]
    mw = stats.mannwhitneyu(refr, avid, alternative="two-sided") if len(avid) and len(refr) else None
    pd.DataFrame([{
        "subset": label,
        "n": len(sub),
        "n_avid": len(avid),
        "n_refractory": len(refr),
        "mean_refractory": refr.mean(),
        "mean_avid": avid.mean(),
        "cohens_d_refractory_vs_avid": cohen_d(refr, avid),
        "mw_p": mw.pvalue if mw else np.nan,
        "auc": auc_binary(sub["y"], sub["p_DM1_score"]),
    }]).to_csv(OUT / f"rai_refractory_signal_{label}.tsv", sep="\t", index=False)
rai_all = pd.read_csv(OUT / "rai_refractory_signal_all.tsv", sep="\t").iloc[0]
summary.append({
    "layer": "RAI-refractory labelled public cohort",
    "n": int(rai_all["n"]),
    "datasets": 1,
    "best_number": f"AUC {rai_all['auc']:.2f}; d {rai_all['cohens_d_refractory_vs_avid']:.2f}; p {rai_all['mw_p']:.2g}",
    "support": "Weak-to-moderate",
    "gap_closed": "Adds a labelled RAI-refractory/avid endpoint proxy.",
    "remaining_gap": "Small n, mixed primary/metastatic context, not an in-house cohort.",
})

# 3) Methylation state support
meth = pd.read_csv(ROOT / "project/results/audit_2026_04_30/round5/r5_2_per_gene_methylation_DM.tsv", sep="\t")
meth_sig = meth[meth["mw_p"] < 0.05]
meth_strong = meth[meth["cohens_d"] >= 0.8]
summary.append({
    "layer": "Promoter methylation support",
    "n": int(meth["n_DM1"].max() + meth["n_DM2"].max()),
    "datasets": 1,
    "best_number": f"{len(meth_sig)}/8 genes p<0.05; {len(meth_strong)}/8 genes d>=0.8; TPO d={meth.loc[meth['gene']=='TPO','cohens_d'].iloc[0]:.2f}",
    "support": "Moderate",
    "gap_closed": "Supports epigenetic silencing as a plausible lineage-collapse mechanism.",
    "remaining_gap": "Does not prove methylation causes expression loss without paired causal/functional validation.",
})

# 4) Existing lineage-state master
lineage = pd.read_csv(ROOT / "project/results/p3_p9_full_execution/shared/lineage_state_master.tsv", sep="\t")
hist_counts = lineage.groupby(["dataset", "histology"]).size().reset_index(name="n")
hist_counts.to_csv(OUT / "lineage_state_sample_inventory.tsv", sep="\t", index=False)
summary.append({
    "layer": "Integrated lineage-state sample inventory",
    "n": int(len(lineage)),
    "datasets": int(lineage["dataset"].nunique()),
    "best_number": f"{len(lineage)} samples with RAI/DM1/non-overlap/TF/STAT3 scores",
    "support": "Moderate",
    "gap_closed": "Creates a unified state-score table across already-local cohorts.",
    "remaining_gap": "Needs careful cohort separation; cannot be pooled naively.",
})

# 5) Cell-line / DepMap translational foothold
cell = pd.read_csv(ROOT / "project/results/paper9_sl_first_pass/ccle_dm1_state.tsv", sep="\t")
thyroid_lines = len(cell)
dm1_high = int(cell["DM1_high"].sum())
atc = cell[cell["hist_subtype"].astype(str).str.contains("anaplastic", case=False, na=False)]
atc_high = int(atc["DM1_high"].sum())
summary.append({
    "layer": "Thyroid cell-line foothold",
    "n": thyroid_lines,
    "datasets": 1,
    "best_number": f"{dm1_high}/{thyroid_lines} thyroid lines DM1-high; {atc_high}/{len(atc)} ATC lines DM1-high",
    "support": "Weak",
    "gap_closed": "Shows possible model systems for follow-up perturbation/drug work.",
    "remaining_gap": "Tiny cell-line set; descriptive only unless dependency/drug-response validation is added.",
})

# 6) RAI data landscape
land = pd.read_csv(ROOT / "project/results/v17_ultimate/U3A_rai_landscape.tsv", sep="\t")
open_high = land[(land["access"].astype(str).str.contains("OPEN", na=False)) & (land["score"] >= 4)]
controlled = land[land["access"].astype(str).str.contains("CONTROLLED", na=False)]
summary.append({
    "layer": "Additional public / controlled RAI-data opportunity",
    "n": int(land["n"].sum()),
    "datasets": int(len(land)),
    "best_number": f"{len(open_high)} high-value open entries; controlled gold-standard n={int(controlled['n'].sum()) if len(controlled) else 0}",
    "support": "Opportunity",
    "gap_closed": "Identifies where to seek stronger clinical-endpoint evidence.",
    "remaining_gap": "Some strongest labels require DAC/controlled access and time.",
})

scorecard = pd.DataFrame(summary)
scorecard.to_csv(OUT / "nature_cancer_gap_scorecard.tsv", sep="\t", index=False)

opportunities = pd.DataFrame([
    {
        "analysis": "Clinical-endpoint replication using GSE151181/GSE151179",
        "can_do_now": "Partly already done",
        "expected_gain": "Moderate",
        "why_it_matters": "Moves from expression-only validation toward RAI-refractory/avid relevance.",
        "risk": "Small n and mixed primary/LNM context.",
    },
    {
        "analysis": "Formal DAC request for Mu 2024 / HRA004166 n=214",
        "can_do_now": "Prepare request, not analyze immediately",
        "expected_gain": "High",
        "why_it_matters": "Best labelled I-RAIA vs I-RAIR clinical cohort in current registry.",
        "risk": "Controlled access, timeline uncertainty.",
    },
    {
        "analysis": "Paired methylation-expression coupling for RAI-lineage genes",
        "can_do_now": "Likely possible if paired TCGA matrices are local",
        "expected_gain": "Moderate-to-high",
        "why_it_matters": "Strengthens DNMT/methylation mechanism beyond promoter beta differences.",
        "risk": "Still correlative without perturbation.",
    },
    {
        "analysis": "IHC proxy mining from pathology annotations or public tissue images",
        "can_do_now": "Only if curated labels/images exist locally",
        "expected_gain": "Low-to-moderate",
        "why_it_matters": "Could create a bridge toward own-cohort validation.",
        "risk": "Previous H&E-to-DM1 battery was negative; avoid image prediction claims.",
    },
    {
        "analysis": "DepMap dependency/drug-response correlation restricted to thyroid lines",
        "can_do_now": "Possible with local DepMap tables if already downloaded",
        "expected_gain": "Low-to-moderate",
        "why_it_matters": "Suggests model systems and functional hypotheses.",
        "risk": "n=13 thyroid lines is too small for strong claims.",
    },
    {
        "analysis": "Minimal in-house validation: qPCR/IHC on archived FFPE",
        "can_do_now": "Requires sample access and wet-lab",
        "expected_gain": "High",
        "why_it_matters": "This is the cleanest way to make the Nature Cancer story credible.",
        "risk": "Not computational; needs approvals, tissue, budget, time.",
    },
])
opportunities.to_csv(OUT / "next_analysis_opportunities.tsv", sep="\t", index=False)

# Probability scenarios as an explicit, honest table.
scenarios = pd.DataFrame([
    {"scenario": "Current computational package", "probability_low": 5, "probability_high": 10, "note": "Strong story pressure-test, but no own data or causal validation."},
    {"scenario": "Lineage-axis rewrite + current scorecard", "probability_low": 10, "probability_high": 20, "note": "Best no-new-data route; still a reach."},
    {"scenario": "Add stronger external clinical-endpoint replication", "probability_low": 15, "probability_high": 25, "note": "Requires robust RAI-refractory/avid data, ideally independent and larger."},
    {"scenario": "Add in-house IHC/qPCR or functional validation", "probability_low": 25, "probability_high": 40, "note": "First point where Nature Cancer becomes a serious reach rather than symbolic reach."},
])
scenarios.to_csv(OUT / "venue_probability_scenarios.tsv", sep="\t", index=False)

# Gap ladder figure.
pillars = [
    ("Novel state\nbiology", 4.0),
    ("Anti-circular\nvalidation", 4.5),
    ("External\nexpression", 4.0),
    ("Clinical\nendpoint", 2.0),
    ("Mechanism\ncausality", 2.0),
    ("Own cohort /\nwet-lab", 0.5),
]
labels, values = zip(*pillars)
colors = ["#2f6f73", "#2f6f73", "#2e5b87", "#a97922", "#a97922", "#8b2635"]
fig, ax = plt.subplots(figsize=(9, 4.8))
ax.bar(labels, values, color=colors)
ax.set_ylim(0, 5)
ax.set_ylabel("Current evidence strength (0-5)")
ax.set_title("Nature Cancer Gap Ladder for Paper 1")
ax.axhline(3.5, color="#111827", linestyle="--", linewidth=1, alpha=0.7)
ax.text(5.35, 3.58, "reach bar", ha="right", va="bottom", fontsize=9, color="#111827")
for i, v in enumerate(values):
    ax.text(i, v + 0.08, f"{v:.1f}", ha="center", va="bottom", fontsize=10)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(OUT / "nature_cancer_gap_ladder.png", dpi=220)
plt.close(fig)

report = {
    "out_dir": str(OUT),
    "scorecard_rows": len(scorecard),
    "external_direction_consistent": int(direction_df["direction_consistent"].sum()),
    "external_direction_cells": int(len(direction_df)),
    "rai_refractory_auc_all": float(rai_all["auc"]),
    "methylation_sig_genes": int(len(meth_sig)),
    "cellline_thyroid_n": thyroid_lines,
}
with open(OUT / "summary.json", "w") as f:
    json.dump(report, f, indent=2)

print(json.dumps(report, indent=2))
