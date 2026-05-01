#!/usr/bin/env python3
"""D5-P6 — GSE286332 BCR repertoire diversity + clonality + TLS + AICDA.

Tests whether PTC+HT immunoglobulin signature reflects antigen-driven clonal
expansion vs polyclonal infiltration. Mechanism layer for autoimmune-PTC paper.
"""
from __future__ import annotations
from pathlib import Path
import json, re
import numpy as np
import pandas as pd
from scipy import stats

PROJ = Path("/opt/thyroid-dash/project")
RES = PROJ / "results/d5p6_bcr_repertoire"
RES.mkdir(parents=True, exist_ok=True)

# ============================================================
# 1. Load GSE286332 (FPKM)
# ============================================================
DATA = "/data/thca/v17_korean/GSE286332/GSE286332_all_sample_rawdata.txt.gz"
print("=== 1. Load GSE286332 ===")
df = pd.read_csv(DATA, sep="\t", low_memory=False)
fpkm_cols = [c for c in df.columns if c.endswith("_FPKM")]
fpkm = df.groupby("Gene_Symbol", as_index=True)[fpkm_cols].sum()
fpkm.columns = [c.replace("_FPKM", "") for c in fpkm.columns]
log_fpkm = np.log2(fpkm + 1.0)
print(f"  log_fpkm shape: {log_fpkm.shape}")

# Group labels
samples = log_fpkm.columns.tolist()
group = pd.Series(["PTC" if s.startswith("NG_") else "PTC_HT" for s in samples], index=samples, name="group")

# ============================================================
# 2. Extract IG gene expression (V/J/C all chains)
# ============================================================
print("\n=== 2. IG repertoire ===")
ig_pattern = re.compile(r"^(IGH|IGK|IGL)[VJCD]")
ig_genes = [g for g in fpkm.index if isinstance(g, str) and ig_pattern.match(g)]
print(f"  IG gene-symbols found: {len(ig_genes)}")

ighv = [g for g in ig_genes if g.startswith("IGHV")]
ighj = [g for g in ig_genes if g.startswith("IGHJ")]
ighc_d = [g for g in ig_genes if g.startswith("IGH") and re.match(r"^IGH[CDMJAEG]", g) and not g.startswith("IGHV") and not g.startswith("IGHJ")]
igkv = [g for g in ig_genes if g.startswith("IGKV")]
iglv = [g for g in ig_genes if g.startswith("IGLV")]
print(f"  IGHV: {len(ighv)}, IGHJ: {len(ighj)}, IGHC*: {len(ighc_d)}, IGKV: {len(igkv)}, IGLV: {len(iglv)}")

# ============================================================
# 3. Per-sample diversity metrics
# ============================================================
print("\n=== 3. Per-sample diversity (Shannon entropy + clonality) ===")

def shannon_entropy(p):
    p = p[p > 0]
    return float(-np.sum(p * np.log(p)))

def clonality(p):
    """1 - normalized Shannon = Pielou's J subtracted from 1."""
    p = p[p > 0]
    if len(p) <= 1:
        return 1.0
    H = -np.sum(p * np.log(p))
    Hmax = np.log(len(p))
    return float(1 - H / Hmax) if Hmax > 0 else 0.0

def family_bias(values, names, prefix_len=4):
    """Count fraction in dominant V family (e.g., IGHV3)."""
    fam = {}
    for n, v in zip(names, values):
        f = n[:prefix_len]
        fam[f] = fam.get(f, 0) + v
    total = sum(fam.values())
    if total == 0:
        return {}
    return {k: v / total for k, v in sorted(fam.items(), key=lambda x: -x[1])}

div_rows = []
for s in samples:
    row = {"sample": s, "group": group[s]}
    for chain_name, chain_genes in [("IGHV", ighv), ("IGKV", igkv), ("IGLV", iglv)]:
        if not chain_genes:
            continue
        vals = fpkm.loc[chain_genes, s].values  # raw FPKM
        total = vals.sum()
        if total == 0:
            row[f"{chain_name}_total"] = 0.0
            row[f"{chain_name}_entropy"] = np.nan
            row[f"{chain_name}_clonality"] = np.nan
            continue
        p = vals / total
        row[f"{chain_name}_total"] = float(total)
        row[f"{chain_name}_entropy"] = shannon_entropy(p)
        row[f"{chain_name}_clonality"] = clonality(p)
        # Top 1 dominance
        sorted_p = np.sort(p)[::-1]
        row[f"{chain_name}_top1_pct"] = float(sorted_p[0] * 100)
        row[f"{chain_name}_top3_pct"] = float(sorted_p[:3].sum() * 100)
    # Family bias for IGHV (e.g., IGHV1/2/3/4 dominance)
    if ighv:
        vals = fpkm.loc[ighv, s].values
        fam = family_bias(vals, ighv, prefix_len=5)  # IGHV1, IGHV2, IGHV3, IGHV4
        for k, v in fam.items():
            row[f"fam_{k}"] = float(v)
    div_rows.append(row)
div_df = pd.DataFrame(div_rows).set_index("sample")
print(div_df[["group", "IGHV_total", "IGHV_entropy", "IGHV_clonality", "IGHV_top1_pct"]].to_string())
div_df.to_csv(RES / "per_sample_diversity.tsv", sep="\t")

# ============================================================
# 4. Group comparison
# ============================================================
print("\n=== 4. PTC vs PTC+HT comparison ===")
def cohen_d(a, b):
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan, np.nan
    sp = np.sqrt(((len(a)-1)*np.var(a, ddof=1) + (len(b)-1)*np.var(b, ddof=1)) / (len(a)+len(b)-2))
    d = (np.mean(a) - np.mean(b)) / max(sp, 1e-9)
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return d, p, u

cmp_rows = []
metrics = [c for c in div_df.columns if c not in ("group",)]
for m in metrics:
    a = div_df.loc[div_df["group"] == "PTC_HT", m].values.astype(float)
    b = div_df.loc[div_df["group"] == "PTC", m].values.astype(float)
    if np.all(np.isnan(a)) or np.all(np.isnan(b)):
        continue
    d, p, u = cohen_d(a, b)
    if np.isnan(d):
        continue
    cmp_rows.append(dict(metric=m, mean_PTC_HT=round(np.nanmean(a), 3),
                          mean_PTC=round(np.nanmean(b), 3),
                          cohen_d=round(d, 3), mw_p=p))
cmp_df = pd.DataFrame(cmp_rows).sort_values("cohen_d", key=abs, ascending=False)
print(cmp_df.to_string(index=False))
cmp_df.to_csv(RES / "group_comparison.tsv", sep="\t", index=False)

# ============================================================
# 5. AICDA (SHM proxy)
# ============================================================
print("\n=== 5. AICDA (SHM enzyme) expression ===")
if "AICDA" in fpkm.index:
    aicda_ptc = log_fpkm.loc["AICDA", group=="PTC"].values
    aicda_th = log_fpkm.loc["AICDA", group=="PTC_HT"].values
    d, p, _ = cohen_d(aicda_th, aicda_ptc)
    print(f"  AICDA log2 FPKM: PTC mean={aicda_ptc.mean():.3f}, PTC+HT mean={aicda_th.mean():.3f}, d={d:+.3f}, p={p:.3g}")
    aicda_result = dict(d=round(float(d), 3), p=float(p), mean_PTC=round(float(aicda_ptc.mean()), 3),
                         mean_PTC_HT=round(float(aicda_th.mean()), 3))
else:
    aicda_result = dict(error="AICDA not in expression matrix")
    print("  AICDA not found")

# ============================================================
# 6. TLS 12-gene signature (Cabrita 2020)
# ============================================================
print("\n=== 6. TLS 12-gene (Cabrita 2020) ===")
TLS = ["CCL19", "CCL21", "CXCL13", "CCR7", "CXCR5", "SELL", "LAMP3",
        "MS4A1", "CD79A", "CD79B", "PTGDS", "TRBC2"]
tls_in = [g for g in TLS if g in log_fpkm.index]
print(f"  TLS genes available: {len(tls_in)}/{len(TLS)}: {tls_in}")
tls_z = log_fpkm.loc[tls_in].apply(lambda c: (c - c.mean()) / (c.std(ddof=1) + 1e-9), axis=1)
tls_score = tls_z.mean(axis=0)
tls_ptc = tls_score[group == "PTC"].values
tls_th = tls_score[group == "PTC_HT"].values
d_tls, p_tls, _ = cohen_d(tls_th, tls_ptc)
print(f"  TLS score: PTC mean={tls_ptc.mean():.3f}, PTC+HT mean={tls_th.mean():.3f}, d={d_tls:+.3f}, p={p_tls:.3g}")
tls_df = pd.DataFrame({"TLS_score": tls_score, "group": group})
tls_df.to_csv(RES / "tls_score_per_sample.tsv", sep="\t")

# ============================================================
# 7. Clonality × TF backbone (P3 link)
# ============================================================
print("\n=== 7. IGHV clonality vs 8-gene RAI ===")
panel = pd.read_csv(PROJ / "results/p3_gse286332/8gene_panel_per_sample.tsv", sep="\t", index_col=0)
g8_score = panel["RAI_score_8gene"]
clo = div_df["IGHV_clonality"].astype(float)
joined = pd.concat([clo, g8_score, tls_score, group], axis=1)
joined.columns = ["IGHV_clonality", "g8_RAI", "TLS_score", "group"]
print(joined.to_string())

if joined["IGHV_clonality"].dropna().std() > 1e-3:
    rho_clo_g8, p_clo_g8 = stats.spearmanr(joined["IGHV_clonality"], joined["g8_RAI"])
    rho_tls_g8, p_tls_g8 = stats.spearmanr(joined["TLS_score"], joined["g8_RAI"])
    rho_clo_tls, p_clo_tls = stats.spearmanr(joined["IGHV_clonality"], joined["TLS_score"])
    print(f"\n  ρ(IGHV_clonality, g8_RAI) = {rho_clo_g8:+.3f}, p={p_clo_g8:.3g}")
    print(f"  ρ(TLS_score, g8_RAI) = {rho_tls_g8:+.3f}, p={p_tls_g8:.3g}")
    print(f"  ρ(IGHV_clonality, TLS_score) = {rho_clo_tls:+.3f}, p={p_clo_tls:.3g}")
    spearman_results = dict(
        clo_vs_g8=dict(rho=round(float(rho_clo_g8), 3), p=float(p_clo_g8)),
        tls_vs_g8=dict(rho=round(float(rho_tls_g8), 3), p=float(p_tls_g8)),
        clo_vs_tls=dict(rho=round(float(rho_clo_tls), 3), p=float(p_clo_tls)),
    )
else:
    spearman_results = dict(error="clonality variation too low")

# ============================================================
# 8. Decision
# ============================================================
ighv_ent_d = next((r["cohen_d"] for r in cmp_rows if r["metric"] == "IGHV_entropy"), None)
ighv_clo_d = next((r["cohen_d"] for r in cmp_rows if r["metric"] == "IGHV_clonality"), None)

if ighv_ent_d is not None and ighv_clo_d is not None:
    aicda_d = aicda_result.get("d", 0) if isinstance(aicda_result.get("d"), (int, float)) else 0
    if (ighv_clo_d > 0.5) and (d_tls > 1.0) and (aicda_d > 0.5):
        decision = "STRONG"
        msg = ("PTC+HT shows clonal IGHV (clonality d > 0.5), elevated TLS (d > 1.0), and "
               "AICDA up-regulation — antigen-driven B cell response confirmed.")
    elif d_tls > 1.0:
        decision = "MODERATE"
        msg = ("PTC+HT shows TLS signature elevated but clonality not strongly differentiated "
               "— polyclonal infiltration with germinal center activity.")
    else:
        decision = "WEAK"
        msg = "Insufficient clonal signal at gene-level resolution (BCR-seq would be definitive)."
else:
    decision = "INDETERMINATE"
    msg = "IGHV diversity metrics could not be computed reliably."

print(f"\n=== DECISION: {decision} ===\n  {msg}")

summary = {
    "n_IG_genes": len(ig_genes),
    "n_IGHV": len(ighv), "n_IGKV": len(igkv), "n_IGLV": len(iglv),
    "group_comparison": cmp_rows,
    "AICDA": aicda_result,
    "TLS": dict(d=round(float(d_tls), 3), p=float(p_tls),
                 mean_PTC=round(float(tls_ptc.mean()), 3),
                 mean_PTC_HT=round(float(tls_th.mean()), 3),
                 genes_used=tls_in),
    "spearman_clo_g8_tls": spearman_results,
    "decision": decision,
    "message": msg,
}
(RES / "D5P6_summary.json").write_text(json.dumps(summary, indent=2, default=str))
print(f"\n✓ Outputs to {RES}")
