#!/usr/bin/env python3
"""Sprint M — DoRothEA / CollecTRI TF activity inference on TCGA-THCA.
Identifies TFs differentially active in DM1-high vs DM1-low patients."""
from pathlib import Path
import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import decoupler as dc
from scipy.stats import mannwhitneyu

OUT = Path("project_external_st/results/extra")

# Get CollecTRI regulons (decoupler v2)
print("Loading CollecTRI regulons...")
try:
    net = dc.op.collectri(organism="human")
    print(f"  CollecTRI: {net['source'].nunique()} TFs, {net['target'].nunique()} targets, {len(net)} interactions")
except Exception as e:
    print(f"  collectri failed: {e}, trying dorothea...")
    net = dc.op.dorothea(organism="human", levels=["A","B","C"])
    print(f"  DoRothEA A/B/C: {net['source'].nunique()} TFs, {net['target'].nunique()} targets, {len(net)} interactions")

needed_genes = set(net["target"].unique()) | set(net["source"].unique())
print(f"  needed genes: {len(needed_genes)}")

# Stream pancan, keep needed genes
print("Streaming pancan EBPP (TF + target genes only)...")
rows = {}
with gzip.open("project/data/raw/TCGA_pancan/pancan_geneExp.gz","rt") as f:
    header = next(f).strip().split("\t")
    samples = header[1:]
    for line in f:
        parts = line.rstrip("\n").split("\t")
        sym = parts[0].strip()
        if sym in needed_genes:
            rows[sym] = [float(x) if x not in ("","NA","NaN") else np.nan for x in parts[1:]]
expr = pd.DataFrame(rows, index=samples).T
print(f"  expression: {expr.shape[0]} genes × {expr.shape[1]} samples")

# Filter to TCGA-THCA samples (match by patient prefix)
score = pd.read_csv(OUT / "s_tcga_thca_scored.tsv", sep="\t")
thca_patients = set(score["patient"].tolist())
sample_to_patient = {s: "-".join(s.split("-")[:3]) for s in expr.columns}
thca_samples_in_expr = [s for s, p in sample_to_patient.items() if p in thca_patients]
expr_thca = expr[thca_samples_in_expr].copy()
print(f"  THCA samples in expr: {expr_thca.shape[1]}")

# Drop genes with too many NaN or all-zero
expr_thca = expr_thca.dropna(thresh=int(0.8 * expr_thca.shape[1]))
expr_thca = expr_thca.fillna(0)
print(f"  after NaN filter: {expr_thca.shape[0]} genes")

# decoupler expects samples × genes
mat = expr_thca.T  # samples × genes

# Run univariate linear model
print("Running ULM TF activity inference...")
acts = dc.mt.ulm(data=mat, net=net)
# decoupler v2 returns (estimate_df, pvals_df) tuple
if isinstance(acts, tuple):
    tf_act = acts[0]
elif isinstance(acts, dict):
    tf_act = acts.get("ulm_estimate", list(acts.values())[0])
else:
    tf_act = acts
if hasattr(tf_act, "to_df"):
    tf_act = tf_act.to_df()
print(f"  tf_act shape: {tf_act.shape}")
print(f"  sample row preview: {tf_act.iloc[0, :5].to_dict()}")

# Match samples to DM1 (dedupe patient first to avoid Series-valued .get())
patient_to_dm1 = score.drop_duplicates("patient").set_index("patient")["DM1_like"]
sample_to_dm1 = pd.Series({s: float(patient_to_dm1[p])
                           for s, p in sample_to_patient.items()
                           if s in tf_act.index and p in patient_to_dm1.index})
sample_to_dm1 = sample_to_dm1.dropna()
common = sample_to_dm1.index.intersection(tf_act.index)
score_idx = sample_to_dm1.loc[common].to_frame("DM1_like")
tf_act = tf_act.loc[common]

# DM1-high vs DM1-low: top 25% vs bottom 25%
vals = score_idx["DM1_like"].astype(float).to_numpy()
q25, q75 = np.percentile(vals, [25, 75])
hi_samples = score_idx[score_idx["DM1_like"] >= q75].index
lo_samples = score_idx[score_idx["DM1_like"] <= q25].index
print(f"  DM1-high n={len(hi_samples)}, DM1-low n={len(lo_samples)}")

# Mann-Whitney per TF
diff_rows = []
for tf in tf_act.columns:
    a = tf_act.loc[hi_samples, tf].dropna()
    b = tf_act.loc[lo_samples, tf].dropna()
    if len(a) < 5 or len(b) < 5: continue
    try:
        U, p = mannwhitneyu(a, b, alternative="two-sided")
    except ValueError:
        continue
    diff_rows.append({"TF": tf, "n_high": len(a), "n_low": len(b),
                      "mean_high": a.mean(), "mean_low": b.mean(),
                      "delta": a.mean() - b.mean(), "p": p})
diff = pd.DataFrame(diff_rows)
diff["fdr"] = diff["p"].rank() / len(diff) * diff["p"].shape[0]
from statsmodels.stats.multitest import multipletests
_, fdr_corr, _, _ = multipletests(diff["p"], method="fdr_bh")
diff["fdr"] = fdr_corr
diff = diff.sort_values("p")
diff.to_csv(OUT / "m_tf_activity_diff.tsv", sep="\t", index=False)

print("\n=== Top 20 TFs differential DM1-high vs DM1-low (by p) ===")
print(diff.head(20).to_string(index=False))
print(f"\n  TFs with FDR < 0.05: {(diff['fdr'] < 0.05).sum()}")
print(f"  TFs with FDR < 0.05 AND |delta| > 0.5: {((diff['fdr'] < 0.05) & (diff['delta'].abs() > 0.5)).sum()}")

# Figure: volcano
fig, ax = plt.subplots(figsize=(10, 7))
diff["nlog10p"] = -np.log10(diff["p"].clip(lower=1e-300))
diff["sig"] = (diff["fdr"] < 0.05) & (diff["delta"].abs() > 0.5)
sig = diff[diff["sig"]]; nonsig = diff[~diff["sig"]]
ax.scatter(nonsig["delta"], nonsig["nlog10p"], color="#999999", s=10, alpha=0.5)
ax.scatter(sig["delta"], sig["nlog10p"], color="#962E2E", s=20, alpha=0.8)
for _, r in diff.head(15).iterrows():
    ax.annotate(r["TF"], (r["delta"], r["nlog10p"]), fontsize=9,
                ha="left" if r["delta"] > 0 else "right",
                va="bottom", color="black")
ax.axhline(-np.log10(0.05), color="grey", lw=0.5, ls=":")
ax.axvline(0, color="grey", lw=0.5, ls=":")
ax.set_xlabel("Δ TF activity (DM1-high − DM1-low)")
ax.set_ylabel("-log10(p)")
ax.set_title(f"M — TF master regulators in DM1-high vs DM1-low PTC (TCGA-THCA, n={len(common)})\n"
             f"{(diff['fdr']<0.05).sum()} TFs FDR<0.05, {((diff['fdr']<0.05)&(diff['delta'].abs()>0.5)).sum()} with |Δ|>0.5",
             fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "m_master_regulator_volcano.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"→ {OUT / 'm_master_regulator_volcano.png'}")
