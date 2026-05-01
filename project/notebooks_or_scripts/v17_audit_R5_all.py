"""Round 5 — Next-step analyses: R5-1 MSK fusion, R5-2 methylation, R5-3 DM1 fusion survival, R5-4 Korean fusion proxy."""
from __future__ import annotations
import json
import requests
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, fisher_exact, spearmanr
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import warnings
warnings.filterwarnings("ignore")

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/audit_2026_04_30/round5"
OUT.mkdir(parents=True, exist_ok=True)

# Common
sm = pd.read_csv(ROOT / "project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv",
                  sep="\t", low_memory=False)
sm = sm[(sm["dataset"] == "TCGA-THCA") & (sm["normal_vs_tumor"] == "tumor")].copy()
sm["tcga_short"] = sm["sample_id"].str.slice(0, 12)
sm["age_n"] = pd.to_numeric(sm["age"], errors="coerce")
sm["os_event"] = pd.to_numeric(sm["os_event"], errors="coerce")
sm["os_days"] = pd.to_numeric(sm["os_days"], errors="coerce")

xing = pd.read_csv(ROOT / "project/results/dark_matter_phase1/step6_xing_rescue.tsv", sep="\t")
xing["dm_cluster"] = xing["v17_dark_cluster"].fillna("not_DM")
if "v17_dark_cluster" in sm.columns:
    sm = sm.drop(columns=["v17_dark_cluster"])
sm = sm.merge(xing[["tcga_short", "xing_group", "dm_cluster"]], on="tcga_short", how="left")
sm["dm_cluster"] = sm["dm_cluster"].fillna("not_DM")

# Load TCGA SV from R3
sv = pd.read_csv(ROOT / "project/results/audit_2026_04_30/round3/cbio_sv_thca.tsv",
                  sep="\t", low_memory=False)
sv["sample_short"] = sv["sampleId"].str.slice(0, 12)
sm["any_fusion"] = sm["tcga_short"].isin(sv["sample_short"]).astype(int)

# =================================================================
# R5-1 — MSK-IMPACT fusion validation via cBioPortal
# =================================================================
print("\n=== R5-1 — MSK-IMPACT fusion via cBioPortal ===\n")

api_base = "https://www.cbioportal.org/api"
msk_study = "thyroid_mskcc_2016"

# Check MSK study profiles
r = requests.get(f"{api_base}/studies/{msk_study}/molecular-profiles", timeout=20,
                  headers={"Accept":"application/json"})
print(f"MSK profiles status: {r.status_code}")
msk_sv_profile = None
if r.status_code == 200:
    profiles = r.json()
    for p in profiles:
        if "STRUCTURAL" in p.get("molecularAlterationType", "") or "FUSION" in p.get("molecularAlterationType", ""):
            msk_sv_profile = p["molecularProfileId"]
            print(f"  MSK SV profile: {msk_sv_profile}")
            break
    if not msk_sv_profile:
        print(f"  All profiles: {[p['molecularProfileId'] for p in profiles]}")

# Fetch MSK structural variants
msk_sv = None
if msk_sv_profile:
    payload = {"molecularProfileIds": [msk_sv_profile]}
    r2 = requests.post(f"{api_base}/structural-variant/fetch", json=payload, timeout=30)
    print(f"  MSK SV fetch status: {r2.status_code}")
    if r2.status_code == 200:
        msk_sv = pd.DataFrame(r2.json())
        if not msk_sv.empty:
            print(f"  MSK SV records: {len(msk_sv)}, samples with SV: {msk_sv['sampleId'].nunique()}")
            msk_sv.to_csv(OUT / "msk_sv_thca.tsv", sep="\t", index=False)

            # Top MSK fusion partners
            if "site1HugoSymbol" in msk_sv.columns and "site2HugoSymbol" in msk_sv.columns:
                msk_sv["fusion_pair"] = msk_sv["site1HugoSymbol"].astype(str) + "-" + msk_sv["site2HugoSymbol"].astype(str)
                print(f"\n  Top MSK fusion pairs:")
                print(msk_sv["fusion_pair"].value_counts().head(15))

                def ftype(row):
                    s1 = str(row.get("site1HugoSymbol", ""))
                    s2 = str(row.get("site2HugoSymbol", ""))
                    pair = f"{s1}-{s2}"
                    if "RET" in pair: return "RET"
                    if "NTRK" in pair: return "NTRK"
                    if "ALK" in pair: return "ALK"
                    if "BRAF" in pair: return "BRAF"
                    if "PAX8" in pair and "PPARG" in pair: return "PAX8-PPARG"
                    if "THADA" in pair: return "THADA"
                    return "Other"
                msk_sv["fusion_class"] = msk_sv.apply(ftype, axis=1)
                print(f"\n  MSK fusion class breakdown:")
                print(msk_sv["fusion_class"].value_counts())

# MSK clinical sample-level
msk_clin_s = pd.read_csv(ROOT / "project/results/v17_tert_recovery/v3/v3_thyroid_mskcc_2016_clinical_sample.tsv",
                          sep="\t", comment="#")
print(f"\nMSK clinical samples: {len(msk_clin_s)}")

r5_1 = {
    "msk_sv_profile": msk_sv_profile,
    "msk_sv_records": int(len(msk_sv)) if msk_sv is not None and not msk_sv.empty else 0,
    "msk_samples_with_sv": int(msk_sv['sampleId'].nunique()) if msk_sv is not None and not msk_sv.empty else 0,
    "msk_clinical_n": int(len(msk_clin_s)),
    "fusion_class_breakdown": msk_sv['fusion_class'].value_counts().to_dict() if msk_sv is not None and not msk_sv.empty and 'fusion_class' in msk_sv.columns else None,
    "msk_cohort_context": "MSK-IMPACT thyroid 2016 (Landa 2016 Cell) is enriched for advanced disease (PDTC + ATC). Higher fusion+ rate expected vs TCGA primary PTC.",
}
(OUT / "r5_1_msk_fusion.json").write_text(json.dumps(r5_1, indent=2, default=str), encoding="utf-8")

# =================================================================
# R5-2 — TCGA methylation via cBioPortal
# =================================================================
print("\n\n=== R5-2 — TCGA methylation via cBioPortal API ===\n")

study_id = "thca_tcga_pan_can_atlas_2018"
# Get methylation profile id
r = requests.get(f"{api_base}/studies/{study_id}/molecular-profiles", timeout=20,
                  headers={"Accept":"application/json"})
meth_profile = None
if r.status_code == 200:
    profiles = r.json()
    for p in profiles:
        mt = p.get("molecularAlterationType", "")
        if "METHYLATION" in mt:
            meth_profile = p["molecularProfileId"]
            break
print(f"  Methylation profile: {meth_profile}")

# Try fetch methylation for 8-gene promoters
P8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
gene_ids = {"SLC5A5": 6528, "TPO": 7173, "TG": 7038, "TSHR": 7253,
             "PAX8": 7849, "NKX2-1": 7080, "FOXE1": 2304, "DIO1": 1733}

meth_data = []
sample_list_id = f"{study_id}_all"
for gene_name, eid in gene_ids.items():
    payload = {
        "entrezGeneIds": [eid],
        "sampleListId": sample_list_id,
    }
    try:
        r3 = requests.post(f"{api_base}/molecular-profiles/{meth_profile}/molecular-data/fetch",
                            json=payload, timeout=60)
        if r3.status_code == 200:
            data = r3.json()
            if data:
                df = pd.DataFrame(data)
                df["gene"] = gene_name
                meth_data.append(df)
                print(f"  {gene_name}: {len(df)} sample β-values fetched")
    except Exception as e:
        print(f"  {gene_name}: {e}")

if meth_data:
    meth_df = pd.concat(meth_data, ignore_index=True)
    meth_df["sample_short"] = meth_df["sampleId"].astype(str).str.slice(0, 12)
    meth_df["beta"] = pd.to_numeric(meth_df["value"], errors="coerce")

    # Average β per sample (across HM450 probes for each gene)
    sample_gene_beta = meth_df.groupby(["sample_short", "gene"])["beta"].mean().unstack()
    print(f"\n  Sample × gene β matrix shape: {sample_gene_beta.shape}")

    # Per-sample mean β for the 8 genes (= "promoter methylation index")
    sample_gene_beta["mean_8g_beta"] = sample_gene_beta.mean(axis=1)
    sample_gene_beta = sample_gene_beta.reset_index()

    # Merge with DM cluster
    dm_meth = sample_gene_beta.merge(sm[["tcga_short", "dm_cluster", "any_fusion"]],
                                       left_on="sample_short", right_on="tcga_short", how="left")

    # 8-gene mean β by DM cluster
    print(f"\n  Mean 8-gene β by DM cluster:")
    print(dm_meth.groupby("dm_cluster")["mean_8g_beta"].agg(["mean", "median", "std", "count"]))

    # MW DM1 vs DM2
    dm1_b = dm_meth[dm_meth["dm_cluster"] == "DM1"]["mean_8g_beta"].dropna()
    dm2_b = dm_meth[dm_meth["dm_cluster"] == "DM2"]["mean_8g_beta"].dropna()
    if len(dm1_b) > 5 and len(dm2_b) > 5:
        mw = mannwhitneyu(dm1_b, dm2_b, alternative="two-sided")
        cd = (dm1_b.mean() - dm2_b.mean()) / np.sqrt((dm1_b.var() + dm2_b.var()) / 2)
        print(f"  DM1 vs DM2 mean β: d={cd:.3f}, p={mw.pvalue:.4e}")

    # Per-gene β by DM cluster
    perg_results = []
    for g in P8:
        if g in dm_meth.columns:
            a = dm_meth[dm_meth["dm_cluster"] == "DM1"][g].dropna()
            b = dm_meth[dm_meth["dm_cluster"] == "DM2"][g].dropna()
            if len(a) > 5 and len(b) > 5:
                mw = mannwhitneyu(a, b, alternative="two-sided")
                cd = (a.mean() - b.mean()) / np.sqrt((a.var() + b.var()) / 2)
                perg_results.append({"gene": g, "n_DM1": int(len(a)), "n_DM2": int(len(b)),
                                      "mean_β_DM1": round(float(a.mean()), 3),
                                      "mean_β_DM2": round(float(b.mean()), 3),
                                      "cohens_d": round(float(cd), 3),
                                      "mw_p": float(mw.pvalue)})
    perg_df = pd.DataFrame(perg_results)
    print(f"\n  Per-gene β DM1 vs DM2:")
    print(perg_df.to_string(index=False))
    perg_df.to_csv(OUT / "r5_2_per_gene_methylation_DM.tsv", sep="\t", index=False)

    # Fusion+ vs fusion- methylation in DM1
    dm1_data = dm_meth[dm_meth["dm_cluster"] == "DM1"].dropna(subset=["mean_8g_beta", "any_fusion"])
    fp_b = dm1_data[dm1_data["any_fusion"] == 1]["mean_8g_beta"]
    fn_b = dm1_data[dm1_data["any_fusion"] == 0]["mean_8g_beta"]
    fusion_meth_result = None
    if len(fp_b) > 5 and len(fn_b) > 5:
        mw = mannwhitneyu(fp_b, fn_b, alternative="two-sided")
        cd = (fp_b.mean() - fn_b.mean()) / np.sqrt((fp_b.var() + fn_b.var()) / 2)
        fusion_meth_result = {"n_fusion+": int(len(fp_b)), "n_fusion-": int(len(fn_b)),
                               "mean_β_fusion+": round(float(fp_b.mean()), 3),
                               "mean_β_fusion-": round(float(fn_b.mean()), 3),
                               "cohens_d": round(float(cd), 3),
                               "mw_p": float(mw.pvalue)}
        print(f"\n  Within DM1, β fusion+ vs fusion-: d={cd:.3f}, p={mw.pvalue:.4f}")

    r5_2 = {
        "n_samples_with_meth": int(sample_gene_beta["sample_short"].nunique()) if "sample_short" in sample_gene_beta.columns else 0,
        "genes_fetched": list(gene_ids.keys()),
        "per_gene_DM1_vs_DM2": perg_results,
        "DM1_fusion_status_methylation": fusion_meth_result,
    }
    (OUT / "r5_2_methylation_DM.json").write_text(json.dumps(r5_2, indent=2, default=str))
    sample_gene_beta.to_csv(OUT / "r5_2_sample_methylation_8gene.tsv", sep="\t", index=False)
else:
    r5_2 = {"status": "no_methylation_data_fetched"}
    print("  No methylation data retrieved")

# =================================================================
# R5-3 — DM1 fusion+ vs fusion- survival Cox HR
# =================================================================
print("\n\n=== R5-3 — DM1 fusion+ vs fusion- survival ===\n")

dm1 = sm[sm["dm_cluster"] == "DM1"].copy()
dm1_surv = dm1.dropna(subset=["os_event", "os_days"])
print(f"DM1 with survival: {len(dm1_surv)}")
print(f"DM1 events: {int(dm1_surv['os_event'].sum())}")

dm1_fp = dm1_surv[dm1_surv["any_fusion"] == 1]
dm1_fn = dm1_surv[dm1_surv["any_fusion"] == 0]
print(f"DM1 fusion+ n={len(dm1_fp)}, events={int(dm1_fp['os_event'].sum())}")
print(f"DM1 fusion- n={len(dm1_fn)}, events={int(dm1_fn['os_event'].sum())}")

# Cox HR fusion+ vs fusion-
sub = dm1_surv.copy()
sub["G"] = sub["any_fusion"]
cox_result = None
if sub["G"].nunique() == 2 and sub["os_event"].sum() >= 2:
    try:
        c = CoxPHFitter(penalizer=0.05)
        c.fit(sub[["os_days", "os_event", "G"]].rename(columns={"os_days":"T", "os_event":"E"}),
              duration_col="T", event_col="E", show_progress=False)
        coef = c.summary.loc["G", "coef"]
        ci_lo = c.summary.loc["G", "coef lower 95%"]
        ci_hi = c.summary.loc["G", "coef upper 95%"]
        cox_result = {"hr": round(float(np.exp(coef)), 3),
                       "ci_lo": round(float(np.exp(ci_lo)), 3),
                       "ci_hi": round(float(np.exp(ci_hi)), 3),
                       "p": round(float(c.summary.loc["G", "p"]), 4),
                       "interpretation": "HR < 1 → fusion+ better survival; HR > 1 → fusion- better"}
        print(f"\nCox HR (fusion+ vs fusion- in DM1):")
        print(json.dumps(cox_result, indent=2))
    except Exception as e:
        cox_result = {"error": str(e)}

# Logrank
try:
    lr = logrank_test(dm1_fp["os_days"], dm1_fn["os_days"],
                       dm1_fp["os_event"], dm1_fn["os_event"])
    logrank_p = float(lr.p_value)
except Exception:
    logrank_p = None

r5_3 = {
    "n_dm1_with_survival": int(len(dm1_surv)),
    "dm1_total_events": int(dm1_surv["os_event"].sum()),
    "n_fusion+": int(len(dm1_fp)),
    "n_fusion-": int(len(dm1_fn)),
    "events_fusion+": int(dm1_fp["os_event"].sum()),
    "events_fusion-": int(dm1_fn["os_event"].sum()),
    "cox_HR_fusion+_vs_fusion-": cox_result,
    "logrank_p": logrank_p,
    "underpowered_warning": "events < 5 — exploratory only" if dm1_surv["os_event"].sum() < 5 else "moderate",
}
(OUT / "r5_3_dm1_fusion_survival.json").write_text(json.dumps(r5_3, indent=2, default=str))

# =================================================================
# R5-4 — Korean cohort fusion proxy (K2 RET expression)
# =================================================================
print("\n\n=== R5-4 — K2 RET expression aberrant pattern ===\n")

# K2 8-gene TPM
k2_tpm = pd.read_csv(ROOT / "project/results/v17_korean/K2_8gene_tpm_matrix_v4.tsv", sep="\t")
print(f"K2 samples: {len(k2_tpm)}")
print(f"K2 cols: {list(k2_tpm.columns)}")

# Whole transcriptome k2 TPM (검색)
k2_full_paths = list((ROOT / "project/results/v17_korean").glob("*tpm*.tsv")) + \
                 list((ROOT / "project/results/v17_korean").glob("*.h5ad"))
print(f"K2 expression files: {[str(p.name) for p in k2_full_paths]}")

# K2 panel score breakdown
if "p_DM2" in k2_tpm.columns:
    print("K2 has p_DM2 column")

# Check predictions file
k2_pred = pd.read_csv(ROOT / "project/results/v17_korean/K2_korean_predictions_v4.tsv", sep="\t")
print(f"\nK2 predictions: {len(k2_pred)} samples")
print(f"  Columns: {list(k2_pred.columns)}")

# K2에 RET expression이 직접 없으니 panel score 분포 + DM_call로 fusion proxy 추정
# DM_call = "DM1" 환자가 fusion+ 가능성 높다는 hypothesis from R3 F4
if "DM_call" in k2_pred.columns:
    print(f"\nK2 DM_call distribution: {k2_pred['DM_call'].value_counts().to_dict()}")

# Korean GSE213647 — has full panel_z for predictions
gse_score = pd.read_csv(ROOT / "project/results/v17_korean/GSE213647_panel_score.tsv", sep="\t")
print(f"\nGSE213647 samples: {len(gse_score)}")
print(f"  Columns: {list(gse_score.columns)[:10]}")
if "panel_z" in gse_score.columns:
    # qcut into 3 categories
    gse_score["panel_z"] = pd.to_numeric(gse_score["panel_z"], errors="coerce")
    gse_score["dm_score_tertile"] = pd.qcut(gse_score["panel_z"].dropna(), 3, labels=["low_DM1like", "mid", "high_DM2like"])
    print(f"\nGSE213647 panel_z tertiles: {gse_score['dm_score_tertile'].value_counts().to_dict()}")
    # By histology
    if "histology" in gse_score.columns:
        ct = pd.crosstab(gse_score["histology"], gse_score["dm_score_tertile"])
        print(f"\nHistology × DM tertile:")
        print(ct)
        ct.to_csv(OUT / "r5_4_gse213647_histology_dm_tertile.tsv", sep="\t")

r5_4 = {
    "k2_n": int(len(k2_tpm)),
    "k2_DM_call_distribution": k2_pred["DM_call"].value_counts().to_dict() if "DM_call" in k2_pred.columns else None,
    "gse213647_n": int(len(gse_score)),
    "limitation": (
        "K2 / GSE213647 raw RNA-seq fusion calling not performed. "
        "Direct RET fusion identification requires STAR-Fusion or Arriba on raw fastq. "
        "K2_8gene_tpm_matrix_v4 contains only 8 panel gene TPMs (mini-index quantification), "
        "not whole transcriptome — direct fusion detection not possible from this file. "
        "Full transcriptome processing required for cross-cohort fusion validation."
    ),
    "next_step": (
        "Pipeline: K2 raw fastq → STAR-Fusion → DM1+ samples 의 RET / NTRK / ALK fusion 빈도 측정. "
        "K2 28명이 DM1-classified (mini-index calibration adjustment 후) 에서 fusion+ 비율 expected ~20-30%. "
        "이는 별도 작업 (1-2주)."
    ),
}
(OUT / "r5_4_korean_fusion_proxy.json").write_text(json.dumps(r5_4, indent=2, default=str))

print("\n\n=== R5 done ===")
