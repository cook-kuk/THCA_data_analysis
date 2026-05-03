#!/usr/bin/env python3
"""Sprint T (simplified) — Druggable target overlap.
DM1-high signature top genes → DGIdb / Open Targets API → druggable target candidates.
This is NOT direct drug-sensitivity (would need DepMap full); it identifies TARGETS in the DM1-high
gene signature that have approved or investigational drugs."""
from pathlib import Path
import gzip
import json
import time
import numpy as np
import pandas as pd
import requests
import matplotlib.pyplot as plt

OUT = Path("project_external_st/results/extra")

# ===== Step 1: differential expression DM1-high vs DM1-low (TCGA-THCA) =====
# Read scored TCGA — score TSV is indexed by patient with DM1_like
score = pd.read_csv(OUT / "s_tcga_thca_scored.tsv", sep="\t")
thca_patients = set(score["patient"].tolist())
print(f"THCA patients: {len(thca_patients)}")

# Stream pancan, load all genes for THCA only (~561 columns × 20k+ genes)
print("Streaming pancan for THCA samples (all genes)...")
header = None
data_rows = {}
target_indices = None
with gzip.open("project/data/raw/TCGA_pancan/pancan_geneExp.gz","rt") as f:
    header = next(f).strip().split("\t")
    samples = header[1:]
    # match by patient prefix (TCGA-XX-XXXX from TCGA-XX-XXXX-01)
    target_indices = [i for i, s in enumerate(samples)
                      if "-".join(s.split("-")[:3]) in thca_patients]
    target_samples = [samples[i] for i in target_indices]
    print(f"  THCA samples found in pancan: {len(target_samples)}")
    n_genes_loaded = 0
    for line in f:
        parts = line.rstrip("\n").split("\t")
        sym = parts[0].strip()
        if sym.startswith("?"): continue
        try:
            vals = [float(parts[i+1]) if parts[i+1] not in ("","NA","NaN") else np.nan for i in target_indices]
        except (ValueError, IndexError):
            continue
        data_rows[sym] = vals
        n_genes_loaded += 1
expr = pd.DataFrame(data_rows, index=target_samples).T
print(f"  loaded {expr.shape[0]} genes × {expr.shape[1]} samples")

# DM1 high/low groups — match pancan sample → score patient
sample_to_patient = {s: "-".join(s.split("-")[:3]) for s in expr.columns}
# dedupe: first DM1 per patient
patient_to_dm1 = score.drop_duplicates("patient").set_index("patient")["DM1_like"]
score_idx = pd.Series({s: float(patient_to_dm1[p]) for s, p in sample_to_patient.items()
                        if p in patient_to_dm1.index})
score_idx = score_idx.dropna().astype(float)
common = score_idx.index.intersection(expr.columns)
score_idx = score_idx.loc[common]
expr = expr[common]
vals = score_idx.to_numpy()
q25, q75 = np.percentile(vals, [25, 75])
hi = score_idx[score_idx >= q75].index
lo = score_idx[score_idx <= q25].index
print(f"  DM1 hi n={len(hi)}, lo n={len(lo)}")

# Differential mean (Welch t-test)
from scipy.stats import ttest_ind
de_rows = []
for g in expr.index:
    a = expr.loc[g, hi].dropna()
    b = expr.loc[g, lo].dropna()
    if len(a) < 10 or len(b) < 10: continue
    try:
        t, p = ttest_ind(a, b, equal_var=False)
    except Exception: continue
    de_rows.append({"gene": g, "mean_hi": a.mean(), "mean_lo": b.mean(),
                    "log2FC": np.log2((a.mean() + 1e-3) / (b.mean() + 1e-3))
                                if (a.mean() > 0 and b.mean() > 0) else (a.mean() - b.mean()),
                    "delta_mean": a.mean() - b.mean(), "t": t, "p": p})
de = pd.DataFrame(de_rows)
from statsmodels.stats.multitest import multipletests
_, fdr, _, _ = multipletests(de["p"].fillna(1.0), method="fdr_bh")
de["fdr"] = fdr
de_up = de.sort_values("delta_mean", ascending=False).query("fdr < 0.01 and delta_mean > 0").head(100)
de_dn = de.sort_values("delta_mean").query("fdr < 0.01 and delta_mean < 0").head(100)
de.to_csv(OUT / "t_DM1_DE_genes.tsv", sep="\t", index=False)
print(f"  DM1-high UP genes (FDR<0.01, top 100): {len(de_up)}")
print(f"  DM1-low UP (DM1-high DOWN, FDR<0.01, top 100): {len(de_dn)}")
print("\n=== Top 20 DM1-high UP-regulated genes ===")
print(de_up.head(20)[["gene","mean_hi","mean_lo","delta_mean","p","fdr"]].to_string(index=False))

# ===== Step 2: DGIdb interactions for top DM1-high UP genes =====
print("\nQuerying DGIdb for druggable targets...")
genes = de_up["gene"].head(50).tolist()
batches = [genes[i:i+20] for i in range(0, len(genes), 20)]
all_interactions = []
for batch in batches:
    try:
        r = requests.get("https://dgidb.org/api/v2/interactions.json",
                         params={"genes": ",".join(batch)}, timeout=20)
        if r.status_code == 200:
            data = r.json()
            for match in data.get("matchedTerms", []):
                gene = match.get("geneName", "")
                for interaction in match.get("interactions", []):
                    all_interactions.append({
                        "gene": gene,
                        "drug": interaction.get("drugName", ""),
                        "drug_chembl": interaction.get("drugConceptId", ""),
                        "interaction_types": ",".join(interaction.get("interactionTypes", [])),
                        "sources": ",".join(interaction.get("sources", [])),
                        "score": interaction.get("score", 0),
                    })
        time.sleep(1)
    except Exception as e:
        print(f"  batch failed: {e}")

interactions = pd.DataFrame(all_interactions)
if not interactions.empty:
    print(f"  {len(interactions)} drug-gene interactions for {interactions['gene'].nunique()} genes")
    print("\n=== Top druggable DM1-high targets (by interaction count) ===")
    by_gene = interactions.groupby("gene").agg(
        n_drugs=("drug", "nunique"),
        sample_drugs=("drug", lambda x: ",".join(x.unique()[:5]))
    ).sort_values("n_drugs", ascending=False).head(15)
    print(by_gene.to_string())
    interactions.to_csv(OUT / "t_dgidb_interactions.tsv", sep="\t", index=False)
else:
    print("  No interactions returned (DGIdb API may be down or rate-limited)")

# ===== Step 3: Figure =====
fig, ax = plt.subplots(figsize=(11, 7))
de["nlog10p"] = -np.log10(de["p"].fillna(1.0).clip(lower=1e-300))
de["category"] = "ns"
de.loc[(de["fdr"] < 0.01) & (de["delta_mean"] > 0), "category"] = "up_DM1high"
de.loc[(de["fdr"] < 0.01) & (de["delta_mean"] < 0), "category"] = "up_DM1low"
colors = {"ns":"#cccccc", "up_DM1high":"#962E2E", "up_DM1low":"#3C6B4F"}
for cat, c in colors.items():
    sub = de[de["category"] == cat]
    ax.scatter(sub["delta_mean"], sub["nlog10p"], color=c, s=8, alpha=0.5)
# annotate top 15 druggable
if not interactions.empty:
    top_drug_genes = by_gene.head(15).index.tolist()
    for g in top_drug_genes:
        if g in de["gene"].values:
            r = de[de["gene"] == g].iloc[0]
            ax.annotate(g, (r["delta_mean"], r["nlog10p"]), fontsize=9,
                        ha="left" if r["delta_mean"] > 0 else "right",
                        color="black", fontweight="bold")
ax.axhline(-np.log10(0.01), color="grey", lw=0.5, ls=":")
ax.axvline(0, color="grey", lw=0.5, ls=":")
ax.set_xlabel("Δ mean (DM1-high − DM1-low)"); ax.set_ylabel("-log10(p)")
ax.set_title(f"T — Druggable targets in DM1-high TCGA-THCA signature\n"
             f"DE {len(de_up)} UP + {len(de_dn)} DOWN at FDR<0.01; "
             f"{len(interactions)} DGIdb drug-gene interactions",
             fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "t_drug_target_volcano.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"→ {OUT / 't_drug_target_volcano.png'}")
