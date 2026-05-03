#!/usr/bin/env python3
"""N3 — CCLE thyroid cell line panel + TROP2 expression + DM1_like score.
Identifies thyroid cell lines that recapitulate DM1-high phenotype = candidates
for downstream sacituzumab govitecan (anti-TROP2 ADC) testing."""
from pathlib import Path
import sys
import urllib.request
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUT = Path("project_external_st/results/extra")

# CCLE expression: use Xena CCLE TPM (smaller, faster than full DepMap)
# Try multiple CCLE sources
CCLE_URLS = [
    "https://depmap.org/portal/download/api/download/external?file_name=expression_protein_coding_genes_expected_count.csv&release_name=DepMap+Public+22Q4",
    "https://figshare.com/ndownloader/files/35020924",  # CCLE expression 22Q4
    "https://depmap.org/api/downloads/external?file_name=CCLE_expression.csv&release_name=DepMap+Public+22Q4",
]
LOCAL = Path("project/data/raw/CCLE_expression.csv")
LOCAL.parent.mkdir(parents=True, exist_ok=True)
if not LOCAL.exists():
    print(f"Downloading CCLE expression...")
    for url in CCLE_URLS:
        try:
            print(f"  trying {url[:80]}...")
            req = urllib.request.Request(url, headers={"User-Agent":"Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                with open(LOCAL, "wb") as out:
                    out.write(resp.read())
            if LOCAL.stat().st_size > 1000:
                print(f"  ✓ {LOCAL.stat().st_size / 1e6:.0f} MB")
                break
            else:
                LOCAL.unlink(missing_ok=True)
        except Exception as e:
            print(f"  failed: {e}")
            LOCAL.unlink(missing_ok=True)
    if not LOCAL.exists():
        print("All CCLE sources failed. Falling back to literature-based thyroid TROP2 panel.")
        sys.exit(2)

# Parse GCT format: header line, dimensions, then expression matrix
import gzip
print("Parsing CCLE GCT...")
RAI_8 = ["TPO","DIO1","TSHR","PAX8","TG","FOXE1","NKX2-1","SLC5A5","NKX2_1","TITF1"]
NONOVERLAP = ["SLC26A4","IYD","DUOX1","DUOX2","TFF3","HHEX","GLIS3","DIO2"]
DM1_TARGETS = ["TACSTD2","FN1","KCNN4","NAMPT","LYN","CYP1B1","CREB5","ELF3","MUC21","CEACAM6"]
ALL_GENES = set(RAI_8 + NONOVERLAP + DM1_TARGETS)

rows = {}
samples = None
with gzip.open(LOCAL, "rt") as f:
    line1 = next(f).strip()  # #1.2
    line2 = next(f).strip()  # n_rows n_cols
    header = next(f).strip().split("\t")
    samples = header[2:]  # first two cols are Name, Description
    print(f"  CCLE: {len(samples)} cell lines")
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3: continue
        # Try Description (gene symbol) at col 1
        sym = parts[1].strip()
        if sym in ALL_GENES:
            try:
                vals = [float(x) if x not in ("","NA","NaN") else np.nan for x in parts[2:]]
                rows[sym] = vals
            except ValueError:
                continue
expr = pd.DataFrame(rows, index=samples).T
if "NKX2_1" in expr.index and "NKX2-1" not in expr.index:
    expr = expr.rename(index={"NKX2_1":"NKX2-1"})
elif "TITF1" in expr.index and "NKX2-1" not in expr.index:
    expr = expr.rename(index={"TITF1":"NKX2-1"})
expr = expr[~expr.index.duplicated(keep="first")]
print(f"  expr: {expr.shape[0]} genes × {expr.shape[1]} cell lines")
print(f"  found: {sorted(expr.index)}")

# Filter to thyroid cell lines (CCLE cell line names contain '_THYROID')
thyroid_cells = [c for c in expr.columns if "THYROID" in c.upper()]
print(f"\n  thyroid cell lines: {len(thyroid_cells)}")
print(f"    {thyroid_cells[:20]}")

if len(thyroid_cells) < 5:
    # fallback: cell line name matching
    thyroid_keywords = ["TPC1","KAT","BCPAP","FRO","BHT","C643","SW1736","8505C","TT2609","CAL62","HTH","TT-","ARO","FB-1"]
    thyroid_cells = [c for c in expr.columns if any(k in c.upper() for k in thyroid_keywords)]
    print(f"  by keyword match: {len(thyroid_cells)} thyroid lines")
    print(f"    {thyroid_cells}")

if len(thyroid_cells) == 0:
    print("WARNING: 0 thyroid cell lines, dumping summary on all CCLE...")
    thyroid_cells = expr.columns.tolist()[:50]

# Score
def score(genes, ref_cols):
    avail = [g for g in genes if g in expr.index]
    if not avail: return None
    z = expr.loc[avail].sub(expr.loc[avail, ref_cols].mean(axis=1), axis=0).div(
        expr.loc[avail, ref_cols].std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)

# Use thyroid cells as reference for z-scoring (within-tissue ranking)
ref = thyroid_cells
df = pd.DataFrame(index=thyroid_cells)
df["RAI_8"] = score([g for g in RAI_8 if g not in ("NKX2_1","TITF1")], ref)
df["NONOVERLAP"] = score(NONOVERLAP, ref)
df["DM1_like"] = -df["RAI_8"]
df["TACSTD2_raw_RPKM"] = expr.loc["TACSTD2", thyroid_cells] if "TACSTD2" in expr.index else np.nan
df["FN1_raw_RPKM"] = expr.loc["FN1", thyroid_cells] if "FN1" in expr.index else np.nan
df = df.sort_values("DM1_like", ascending=False)

print("\n=== Thyroid cell lines ranked by DM1_like ===")
print(df.head(15).to_string())
print("\n=== Cell lines in DM1-high (top 5) ===")
top5 = df.head(5)
print(top5.to_string())
print("\n=== TACSTD2 (TROP2) in thyroid cell lines ===")
print(f"  Top expressed (RPKM): {df.nlargest(5, 'TACSTD2_raw_RPKM')[['TACSTD2_raw_RPKM','DM1_like']].to_string()}")

# DM1 vs TROP2 correlation in cell lines
ok = df[["DM1_like","TACSTD2_raw_RPKM"]].dropna()
if len(ok) >= 5:
    from scipy.stats import spearmanr, pearsonr
    rho, p = spearmanr(ok["DM1_like"], ok["TACSTD2_raw_RPKM"])
    r, pp = pearsonr(ok["DM1_like"], ok["TACSTD2_raw_RPKM"])
    print(f"\n  DM1_like vs TACSTD2 in {len(ok)} thyroid cell lines: Spearman ρ = {rho:.3f}, p = {p:.2e}; Pearson r = {r:.3f}")

df.to_csv(OUT / "n3_ccle_thyroid_dm1.tsv", sep="\t")

# Figure
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
ax = axes[0]
df_sorted = df.dropna(subset=["DM1_like"]).sort_values("DM1_like")
colors = ["#962E2E" if x > 0 else "#3C6B4F" for x in df_sorted["DM1_like"]]
ax.barh(df_sorted.index, df_sorted["DM1_like"], color=colors, edgecolor="black")
ax.axvline(0, color="grey", lw=0.5, ls=":")
ax.set_xlabel("DM1_like score (within thyroid cell lines)")
ax.set_title(f"N3.A — Thyroid CCLE cell lines (n={len(df)}) by DM1_like\n"
             f"Top DM1-high lines = candidates for sacituzumab govitecan testing",
             fontsize=11)
ax.tick_params(axis="y", labelsize=8)

ax = axes[1]
ax.scatter(df["DM1_like"], df["TACSTD2_raw_RPKM"], s=80, c="#962E2E",
           edgecolor="black", alpha=0.8)
for i, name in enumerate(df.index):
    if pd.notna(df.loc[name, "TACSTD2_raw_RPKM"]) and pd.notna(df.loc[name, "DM1_like"]):
        ax.annotate(name.split("_")[0], (df.loc[name, "DM1_like"], df.loc[name, "TACSTD2_raw_RPKM"]),
                    fontsize=8, alpha=0.7)
ax.set_xlabel("DM1_like score (cell line)")
ax.set_ylabel("TACSTD2 / TROP2 (RPKM)")
title = f"N3.B — DM1_like vs TROP2 in thyroid cell lines (n={len(df)})"
if len(ok) >= 5:
    title += f"\nSpearman ρ = {rho:.3f}, p = {p:.2e}"
ax.set_title(title, fontsize=11)
ax.axhline(df["TACSTD2_raw_RPKM"].median(), color="grey", lw=0.5, ls=":", label="median")

fig.tight_layout()
fig.savefig(OUT / "n3_ccle_thyroid.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"→ {OUT / 'n3_ccle_thyroid.png'}")
