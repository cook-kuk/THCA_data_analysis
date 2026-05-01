"""Round 6 — R6-1 fusion- DM1 immune deep, R6-2 CNV, R6-3 RPPA, R6-4 methylation-expression correlation.
cBioPortal API parallel fetches."""
from __future__ import annotations
import json
import requests
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, fisher_exact, spearmanr, pearsonr
from concurrent.futures import ThreadPoolExecutor
import warnings
warnings.filterwarnings("ignore")

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/audit_2026_04_30/round6"
OUT.mkdir(parents=True, exist_ok=True)
api = "https://www.cbioportal.org/api"

# Common load
sm = pd.read_csv(ROOT / "project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv",
                  sep="\t", low_memory=False)
sm = sm[(sm["dataset"] == "TCGA-THCA") & (sm["normal_vs_tumor"] == "tumor")].copy()
sm["tcga_short"] = sm["sample_id"].str.slice(0, 12)
sm["age_n"] = pd.to_numeric(sm["age"], errors="coerce")

xing = pd.read_csv(ROOT / "project/results/dark_matter_phase1/step6_xing_rescue.tsv", sep="\t")
xing["dm_cluster"] = xing["v17_dark_cluster"].fillna("not_DM")
if "v17_dark_cluster" in sm.columns:
    sm = sm.drop(columns=["v17_dark_cluster"])
sm = sm.merge(xing[["tcga_short", "dm_cluster"]], on="tcga_short", how="left")
sm["dm_cluster"] = sm["dm_cluster"].fillna("not_DM")

sv = pd.read_csv(ROOT / "project/results/audit_2026_04_30/round3/cbio_sv_thca.tsv", sep="\t", low_memory=False)
sv["sample_short"] = sv["sampleId"].str.slice(0, 12)
sm["any_fusion"] = sm["tcga_short"].isin(sv["sample_short"]).astype(int)

immune = pd.read_csv(ROOT / "project/results/tables/immune_signatures_tcga.tsv", sep="\t")

# ============================================================
# R6-1 — Fusion- DM1 immune-hot deep dive (TLS + BCR)
# ============================================================
print("\n=== R6-1 — Fusion- DM1 TLS/BCR signature ===\n")

# Cabrita 2020 12-gene TLS signature: CCL19, CCL21, CXCL13, CCR6, CCR7, CXCR5, SELL, LAMP3, FCRL5, CD79B, MS4A1, IGHA1
# (Need to fetch via cBioPortal mRNA expression)
TLS_genes = ["CCL19", "CCL21", "CXCL13", "CCR6", "CCR7", "CXCR5",
              "SELL", "LAMP3", "FCRL5", "CD79B", "MS4A1", "IGHA1"]
BCR_genes = ["IGHV3-66", "IGHV3-13", "IGHV2-70", "IGKV1-8", "IGLV6-57",
              "IGHG1", "IGHM", "IGKC"]

# Get gene IDs via cBioPortal /genes endpoint
def fetch_gene_ids(symbols):
    """Get Entrez IDs for gene symbols."""
    ids = {}
    for sym in symbols:
        try:
            r = requests.get(f"{api}/genes/{sym}", timeout=10, headers={"Accept":"application/json"})
            if r.status_code == 200:
                d = r.json()
                ids[sym] = d.get("entrezGeneId")
        except Exception:
            pass
    return ids

print("Fetching gene IDs (TLS + BCR)...")
all_genes = TLS_genes + BCR_genes
gene_ids = fetch_gene_ids(all_genes)
print(f"  Got {len(gene_ids)}/{len(all_genes)} gene IDs")

# Fetch mRNA expression via cBioPortal
study_id = "thca_tcga_pan_can_atlas_2018"
# Get mRNA profile id
r = requests.get(f"{api}/studies/{study_id}/molecular-profiles", timeout=20, headers={"Accept":"application/json"})
mrna_profile = None
if r.status_code == 200:
    profiles = r.json()
    for p in profiles:
        if "mrna_seq" in p.get("molecularProfileId", "").lower() and "v2" in p.get("molecularProfileId", "").lower():
            mrna_profile = p["molecularProfileId"]
            break
    if not mrna_profile:
        for p in profiles:
            if "MRNA" in p.get("molecularAlterationType", "") and "z_scores" in p.get("molecularProfileId", "").lower():
                mrna_profile = p["molecularProfileId"]
                break
print(f"  mRNA profile: {mrna_profile}")

def fetch_gene_expr(gene_eid):
    gene, eid = gene_eid
    payload = {"entrezGeneIds": [eid], "sampleListId": f"{study_id}_all"}
    try:
        r = requests.post(f"{api}/molecular-profiles/{mrna_profile}/molecular-data/fetch",
                           json=payload, timeout=30)
        if r.status_code == 200:
            data = r.json()
            if data:
                df = pd.DataFrame(data)
                df["gene"] = gene
                return gene, df
    except Exception:
        pass
    return gene, None

if mrna_profile and gene_ids:
    valid_pairs = [(g, eid) for g, eid in gene_ids.items() if eid is not None]
    print(f"  Parallel fetch {len(valid_pairs)} genes...")
    with ThreadPoolExecutor(max_workers=4) as ex:
        results = dict(ex.map(fetch_gene_expr, valid_pairs))
    expr_dfs = [df for g, df in results.items() if df is not None and not df.empty]
    if expr_dfs:
        all_expr = pd.concat(expr_dfs, ignore_index=True)
        all_expr["sample_short"] = all_expr["sampleId"].astype(str).str.slice(0, 12)
        all_expr["expr_z"] = pd.to_numeric(all_expr["value"], errors="coerce")

        # Per-sample TLS score (mean z across TLS genes)
        tls_data = all_expr[all_expr["gene"].isin(TLS_genes)]
        tls_score = tls_data.groupby("sample_short")["expr_z"].mean().rename("TLS_score").reset_index()

        bcr_data = all_expr[all_expr["gene"].isin(BCR_genes)]
        bcr_score = bcr_data.groupby("sample_short")["expr_z"].mean().rename("BCR_score").reset_index()

        # Merge with DM cluster
        dm_imm = sm[["tcga_short", "dm_cluster", "any_fusion"]].merge(tls_score, left_on="tcga_short", right_on="sample_short", how="left")
        dm_imm = dm_imm.merge(bcr_score, on="sample_short", how="left")

        # DM1 fusion+ vs fusion-
        dm1_dm_imm = dm_imm[dm_imm["dm_cluster"] == "DM1"]
        fp = dm1_dm_imm[dm1_dm_imm["any_fusion"] == 1]
        fn = dm1_dm_imm[dm1_dm_imm["any_fusion"] == 0]

        r6_1_results = []
        for col in ["TLS_score", "BCR_score"]:
            a = fp[col].dropna()
            b = fn[col].dropna()
            if len(a) > 5 and len(b) > 3:
                mw = mannwhitneyu(a, b, alternative="two-sided")
                cd = (a.mean() - b.mean()) / np.sqrt((a.var() + b.var()) / 2) if (a.var() + b.var()) > 0 else 0
                r6_1_results.append({"score": col, "n_fusion+": int(len(a)), "n_fusion-": int(len(b)),
                                      "median_+": round(float(a.median()), 3),
                                      "median_-": round(float(b.median()), 3),
                                      "cohens_d": round(float(cd), 3),
                                      "mw_p": float(mw.pvalue)})

        print("\nDM1 fusion+ vs fusion- TLS/BCR:")
        print(pd.DataFrame(r6_1_results).to_string(index=False))
        pd.DataFrame(r6_1_results).to_csv(OUT / "r6_1_dm1_fusion_TLS_BCR.tsv", sep="\t", index=False)

        # Also compare DM1 vs DM2
        for col in ["TLS_score", "BCR_score"]:
            a = dm_imm[dm_imm["dm_cluster"] == "DM1"][col].dropna()
            b = dm_imm[dm_imm["dm_cluster"] == "DM2"][col].dropna()
            if len(a) > 5 and len(b) > 5:
                mw = mannwhitneyu(a, b, alternative="two-sided")
                cd = (a.mean() - b.mean()) / np.sqrt((a.var() + b.var()) / 2)
                print(f"\n  DM1 vs DM2 {col}: d={cd:.3f}, p={mw.pvalue:.4e}")
                r6_1_results.append({"comparison": "DM1 vs DM2", "score": col,
                                      "median_DM1": round(float(a.median()), 3),
                                      "median_DM2": round(float(b.median()), 3),
                                      "cohens_d": round(float(cd), 3),
                                      "mw_p": float(mw.pvalue)})

        r6_1 = {
            "tls_genes": TLS_genes,
            "bcr_genes": BCR_genes,
            "dm1_fusion_pos_neg_TLS_BCR": [r for r in r6_1_results if "n_fusion+" in r],
            "DM1_vs_DM2_TLS_BCR": [r for r in r6_1_results if "comparison" in r],
        }
        (OUT / "r6_1_dm1_fusion_immune.json").write_text(json.dumps(r6_1, indent=2, default=str))

# ============================================================
# R6-2 — TCGA CNV via cBioPortal API
# ============================================================
print("\n\n=== R6-2 — TCGA CNV via cBioPortal API ===\n")

# Discrete CNV (GISTIC2 levels: -2 deep del, -1 shallow del, 0 diploid, 1 gain, 2 amp)
cnv_profile = None
r = requests.get(f"{api}/studies/{study_id}/molecular-profiles", timeout=20, headers={"Accept":"application/json"})
if r.status_code == 200:
    profiles = r.json()
    for p in profiles:
        mt = p.get("molecularAlterationType", "")
        if "COPY_NUMBER" in mt and "DISCRETE" in mt:
            cnv_profile = p["molecularProfileId"]
            break
print(f"  CNV profile: {cnv_profile}")

# Key thyroid cancer CNV target genes
cnv_genes = {
    "CDKN2A": 1029, "CDKN2B": 1030, "TP53": 7157, "MYC": 4609,
    "PTEN": 5728, "RB1": 5925, "EGFR": 1956, "MET": 4233,
    "PIK3CA": 5290, "AKT1": 207, "CCND1": 595, "MDM2": 4193,
}

if cnv_profile:
    def fetch_cnv(gene_eid):
        gene, eid = gene_eid
        payload = {"entrezGeneIds": [eid], "sampleListId": f"{study_id}_all"}
        try:
            r = requests.post(f"{api}/molecular-profiles/{cnv_profile}/molecular-data/fetch",
                               json=payload, timeout=30)
            if r.status_code == 200:
                data = r.json()
                if data:
                    df = pd.DataFrame(data)
                    df["gene"] = gene
                    return gene, df
        except Exception:
            pass
        return gene, None

    print("  Parallel CNV fetch...")
    with ThreadPoolExecutor(max_workers=4) as ex:
        cnv_results = dict(ex.map(fetch_cnv, cnv_genes.items()))
    cnv_dfs = [df for g, df in cnv_results.items() if df is not None and not df.empty]

    if cnv_dfs:
        cnv_all = pd.concat(cnv_dfs, ignore_index=True)
        cnv_all["sample_short"] = cnv_all["sampleId"].astype(str).str.slice(0, 12)
        cnv_all["cnv_value"] = pd.to_numeric(cnv_all["value"], errors="coerce")

        # CNV burden per sample = number of altered genes (|CNV| >= 1)
        cnv_all["altered"] = (cnv_all["cnv_value"].abs() >= 1).astype(int)
        cnv_burden = cnv_all.groupby("sample_short")["altered"].sum().rename("cnv_burden").reset_index()

        # Per-gene CNV by DM cluster
        per_gene_cnv = []
        for g in cnv_genes.keys():
            gd = cnv_all[cnv_all["gene"] == g].merge(sm[["tcga_short", "dm_cluster"]],
                                                       left_on="sample_short", right_on="tcga_short", how="left")
            for dm in ["DM1", "DM2"]:
                sub = gd[gd["dm_cluster"] == dm]
                if len(sub) > 5:
                    altered_pct = float(sub["altered"].mean() * 100)
                    deep_del_pct = float((sub["cnv_value"] == -2).mean() * 100)
                    amp_pct = float((sub["cnv_value"] == 2).mean() * 100)
                    per_gene_cnv.append({"gene": g, "dm_cluster": dm, "n": int(len(sub)),
                                          "altered_pct": round(altered_pct, 1),
                                          "deep_del_pct": round(deep_del_pct, 1),
                                          "amp_pct": round(amp_pct, 1)})

        # CNV burden DM1 vs DM2
        cnv_burden_dm = cnv_burden.merge(sm[["tcga_short", "dm_cluster", "any_fusion"]],
                                          left_on="sample_short", right_on="tcga_short", how="left")
        a = cnv_burden_dm[cnv_burden_dm["dm_cluster"] == "DM1"]["cnv_burden"].dropna()
        b = cnv_burden_dm[cnv_burden_dm["dm_cluster"] == "DM2"]["cnv_burden"].dropna()
        cnv_burden_test = None
        if len(a) > 5 and len(b) > 5:
            mw = mannwhitneyu(a, b, alternative="two-sided")
            cd = (a.mean() - b.mean()) / np.sqrt((a.var() + b.var()) / 2)
            cnv_burden_test = {"DM1_mean_burden": round(float(a.mean()), 2),
                                "DM2_mean_burden": round(float(b.mean()), 2),
                                "cohens_d": round(float(cd), 3),
                                "mw_p": float(mw.pvalue)}
            print(f"\n  CNV burden DM1 vs DM2: d={cd:.3f}, p={mw.pvalue:.4e}")

        # Per-gene table
        per_gene_df = pd.DataFrame(per_gene_cnv).pivot_table(
            index="gene", columns="dm_cluster", values=["altered_pct", "deep_del_pct", "amp_pct"]
        )
        print("\n  Per-gene CNV by DM cluster:")
        print(per_gene_df)
        per_gene_df.to_csv(OUT / "r6_2_cnv_per_gene_DM.tsv", sep="\t")

        r6_2 = {
            "cnv_profile": cnv_profile,
            "n_samples_with_cnv": int(cnv_all["sample_short"].nunique()),
            "genes_fetched": list(cnv_genes.keys()),
            "cnv_burden_DM1_vs_DM2": cnv_burden_test,
            "per_gene_summary": per_gene_cnv,
        }
        (OUT / "r6_2_cnv_DM.json").write_text(json.dumps(r6_2, indent=2, default=str))

# ============================================================
# R6-3 — TCGA RPPA (protein) data for DM cluster
# ============================================================
print("\n\n=== R6-3 — TCGA RPPA via cBioPortal API ===\n")

rppa_profile = None
if r.status_code == 200:
    for p in profiles:
        if "PROTEIN" in p.get("molecularAlterationType", ""):
            rppa_profile = p["molecularProfileId"]
            break
print(f"  RPPA profile: {rppa_profile}")

# Key thyroid signaling proteins
rppa_genes = {
    "MAPK1": 5594, "MAPK3": 5595, "AKT1": 207, "MTOR": 2475,
    "PTEN": 5728, "TP53": 7157, "EGFR": 1956, "RAF1": 5894,
    "PAX8": 7849, "NKX2-1": 7080,
}

if rppa_profile:
    def fetch_rppa(gene_eid):
        gene, eid = gene_eid
        payload = {"entrezGeneIds": [eid], "sampleListId": f"{study_id}_all"}
        try:
            r = requests.post(f"{api}/molecular-profiles/{rppa_profile}/molecular-data/fetch",
                               json=payload, timeout=30)
            if r.status_code == 200:
                data = r.json()
                if data:
                    df = pd.DataFrame(data)
                    df["gene"] = gene
                    return gene, df
        except Exception:
            pass
        return gene, None

    print("  Parallel RPPA fetch...")
    with ThreadPoolExecutor(max_workers=4) as ex:
        rppa_results = dict(ex.map(fetch_rppa, rppa_genes.items()))
    rppa_dfs = [df for g, df in rppa_results.items() if df is not None and not df.empty]

    if rppa_dfs:
        rppa_all = pd.concat(rppa_dfs, ignore_index=True)
        rppa_all["sample_short"] = rppa_all["sampleId"].astype(str).str.slice(0, 12)
        rppa_all["rppa_value"] = pd.to_numeric(rppa_all["value"], errors="coerce")

        rppa_with_dm = rppa_all.merge(sm[["tcga_short", "dm_cluster"]],
                                       left_on="sample_short", right_on="tcga_short", how="left")
        per_gene_rppa = []
        for g in rppa_genes.keys():
            gd = rppa_with_dm[rppa_with_dm["gene"] == g]
            a = gd[gd["dm_cluster"] == "DM1"]["rppa_value"].dropna()
            b = gd[gd["dm_cluster"] == "DM2"]["rppa_value"].dropna()
            if len(a) > 5 and len(b) > 5:
                mw = mannwhitneyu(a, b, alternative="two-sided")
                cd = (a.mean() - b.mean()) / np.sqrt((a.var() + b.var()) / 2) if (a.var() + b.var()) > 0 else 0
                per_gene_rppa.append({"protein": g, "n_DM1": int(len(a)), "n_DM2": int(len(b)),
                                       "mean_DM1": round(float(a.mean()), 3),
                                       "mean_DM2": round(float(b.mean()), 3),
                                       "cohens_d": round(float(cd), 3),
                                       "mw_p": float(mw.pvalue)})
        print("\n  RPPA DM1 vs DM2:")
        print(pd.DataFrame(per_gene_rppa).to_string(index=False))
        pd.DataFrame(per_gene_rppa).to_csv(OUT / "r6_3_rppa_DM.tsv", sep="\t", index=False)
        r6_3 = {"rppa_profile": rppa_profile, "per_protein_DM1_vs_DM2": per_gene_rppa,
                 "n_samples_with_rppa": int(rppa_all["sample_short"].nunique())}
        (OUT / "r6_3_rppa_DM.json").write_text(json.dumps(r6_3, indent=2, default=str))

# ============================================================
# R6-4 — Methylation-expression correlation (R5-2 follow-up)
# ============================================================
print("\n\n=== R6-4 — Methylation-expression correlation ===\n")

# Load R5-2 methylation
meth_r5 = pd.read_csv(ROOT / "project/results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv", sep="\t")

# Get 8-gene mRNA expression via cBioPortal
P8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
P8_ids = {"SLC5A5": 6528, "TPO": 7173, "TG": 7038, "TSHR": 7253,
           "PAX8": 7849, "NKX2-1": 7080, "FOXE1": 2304, "DIO1": 1733}

if mrna_profile:
    def fetch_p8_expr(gene_eid):
        gene, eid = gene_eid
        payload = {"entrezGeneIds": [eid], "sampleListId": f"{study_id}_all"}
        try:
            r = requests.post(f"{api}/molecular-profiles/{mrna_profile}/molecular-data/fetch",
                               json=payload, timeout=30)
            if r.status_code == 200:
                data = r.json()
                if data:
                    df = pd.DataFrame(data)
                    df["gene"] = gene
                    return gene, df
        except Exception:
            pass
        return gene, None

    print("  Parallel P8 mRNA fetch...")
    with ThreadPoolExecutor(max_workers=4) as ex:
        p8_results = dict(ex.map(fetch_p8_expr, P8_ids.items()))
    p8_dfs = [df for g, df in p8_results.items() if df is not None and not df.empty]
    if p8_dfs:
        p8_expr = pd.concat(p8_dfs, ignore_index=True)
        p8_expr["sample_short"] = p8_expr["sampleId"].astype(str).str.slice(0, 12)
        p8_expr["expr_z"] = pd.to_numeric(p8_expr["value"], errors="coerce")
        p8_pivot = p8_expr.pivot_table(index="sample_short", columns="gene", values="expr_z", aggfunc="mean")

        # Compute methylation-expression correlation per gene
        meth_pivot = meth_r5.set_index("sample_short")
        # rename columns if needed
        if "mean_8g_beta" in meth_pivot.columns:
            meth_pivot = meth_pivot.drop(columns=["mean_8g_beta"])

        common = meth_pivot.index.intersection(p8_pivot.index)
        print(f"  Samples with both meth + expr: {len(common)}")

        corr_results = []
        for g in P8:
            if g in meth_pivot.columns and g in p8_pivot.columns:
                m = meth_pivot.loc[common, g]
                e = p8_pivot.loc[common, g]
                valid = m.notna() & e.notna()
                if valid.sum() > 30:
                    rho, p = spearmanr(m[valid], e[valid])
                    corr_results.append({"gene": g, "n": int(valid.sum()),
                                          "spearman_rho": round(float(rho), 3),
                                          "p": float(p),
                                          "interpretation": "negative ρ = hypermethylation silences expression"})
        print("\n  Methylation × Expression correlation per gene:")
        print(pd.DataFrame(corr_results).to_string(index=False))
        pd.DataFrame(corr_results).to_csv(OUT / "r6_4_meth_expr_correlation.tsv", sep="\t", index=False)
        r6_4 = {"correlation_per_gene": corr_results,
                 "n_samples_paired": int(len(common)),
                 "interpretation": "Strong negative correlation (ρ < -0.4) supports causal hypermethylation → silencing model. Weak/positive correlation suggests methylation is marker not driver."}
        (OUT / "r6_4_methylation_expression.json").write_text(json.dumps(r6_4, indent=2, default=str))

print("\n\n=== R6 done ===")
