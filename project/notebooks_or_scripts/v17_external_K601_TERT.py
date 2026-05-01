#!/usr/bin/env python3
"""v17 external — K601 + TERT cross-cohort analysis.

(1) TCGA-THCA: re-extract K601E patients from thca_tcga_pub mirror (better than v3 which used STAR-MAF)
(2) MSK-IMPACT 1361 records: fetch via cBioPortal datahub
(3) MSK-CHORD 342 + MSK-PDTC-ATC 81 + ODG-MSK 21
(4) Combined K601E prevalence + BRAF×TERT survival in larger cohort
"""
from __future__ import annotations
import urllib.request, json, gzip, io, time
from pathlib import Path
import pandas as pd, numpy as np

OUT = Path("/opt/thyroid-dash/project/results/v17_tert_recovery/v3")
OUT.mkdir(parents=True, exist_ok=True)
LOG = []
def log(msg):
    print(msg, flush=True); LOG.append(msg)

# ============================================================
# 1. TCGA-THCA K601E re-extraction (cbio thca_tcga_pub mirror)
# ============================================================
log("="*60)
log("[1] TCGA-THCA cbio mutation table — K601E + V600E count")
df = pd.read_csv("/opt/thyroid-dash/project/results/v17_tert_recovery/cbio_data_mutations_raw.txt",
                 sep="\t", low_memory=False)
log(f"   total rows: {len(df)}")
braf = df[df["Hugo_Symbol"]=="BRAF"].copy()
braf["patient_id"] = braf["Tumor_Sample_Barcode"].astype(str).str[:12]
log(f"   BRAF rows: {len(braf)}, unique patients: {braf['patient_id'].nunique()}")

def classify(h):
    s = str(h)
    if "V600E" in s: return "V600E"
    if "V600K" in s: return "V600K"
    if "K601" in s:  return "K601E"
    if "V600" in s:  return "V600_syn_or_other"
    return "other"

braf["braf_class"] = braf["HGVSp_Short"].apply(classify)
log(f"\n   per-patient BRAF class:")
per_pt = braf.drop_duplicates(["patient_id"]).copy()
log(per_pt["braf_class"].value_counts().to_string())

# Find K601E patients
k601 = per_pt[per_pt["braf_class"]=="K601E"]
log(f"\n   K601E patient IDs: {k601['patient_id'].tolist()}")

# Cross-ref with sample_master + TERT integrated
sm = pd.read_csv("/opt/thyroid-dash/project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv",
                 sep="\t")
sm["patient_id"] = sm["tcga_short"].astype(str)
k601_sm = sm[sm["patient_id"].isin(k601["patient_id"])]
log(f"\n   K601E in sample_master:")
for _, r in k601_sm.iterrows():
    log(f"     {r['patient_id']}  TERT={r.get('tert_promoter_integrated','?')}  age={r.get('age',np.nan):.1f}  os_event={r.get('os_event','?')}  os_days={r.get('os_days','?')}  histology={r.get('histology_subtype','?')}")

# ============================================================
# 2. MSK-IMPACT external — fetch via cBioPortal API
# ============================================================
log("\n" + "="*60)
log("[2] MSK external cohorts — TERT promoter + BRAF positions")
log("    (fetching mutation tables from cBioPortal datahub)")

# All thyroid-related MSK studies
MSK_STUDIES = {
    "msk_impact_2017_thyroid":     None,   # from MSK-IMPACT 10K paper if exposed
    "thyroid_mskcc_2016":          "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/thyroid_mskcc_2016/data_mutations.txt",
    "thymic_2014":                 None,   # not relevant
}
# MSK-IMPACT pan-cancer 10K is via GENIE / msk_impact_2017
# Try direct cBioPortal datahub pulls
candidates = [
    ("thyroid_mskcc_2016",      "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/thyroid_mskcc_2016/data_mutations.txt"),
    ("dlbcl_dfci_2018_thyroid", None),  # ignore
]

def safe_fetch(url, timeout=60):
    try:
        log(f"     fetching {url[:80]}...")
        req = urllib.request.Request(url, headers={"User-Agent":"v17-tert/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read().decode("utf-8", errors="replace")
    except Exception as e:
        log(f"     FAIL: {e}")
        return None

# Already have thyroid_mskcc_2016 in S6 records
log("\n   Known external promoter records (from v2 sprint):")
ext = pd.read_csv("/opt/thyroid-dash/project/results/v17_tert_recovery/v2/FINAL_promoter_records_all_sources.tsv",
                  sep="\t")
log(f"     total: {len(ext)}")
ext["pos_label"] = ext["start"].apply(lambda p: "C228T" if p==1295228 else ("C250T" if p==1295250 else "other"))
ext["cohort_short"] = ext["source_file"].astype(str).apply(
    lambda s: s.split("master_public_")[-1].split("_data_")[0] if "master_public_" in s else (s[:30]))
log("\n   external TERT prevalence by cohort + position:")
log(pd.crosstab(ext["cohort_short"], ext["pos_label"]).to_string())

# We need BRAF status too — fetch MSK-IMPACT thyroid mutations directly
log("\n   fetching thyroid_mskcc_2016 full mutation table...")
data = safe_fetch("https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/thyroid_mskcc_2016/data_mutations.txt", 90)
if data:
    msk = pd.read_csv(io.StringIO(data), sep="\t", low_memory=False, comment="#")
    log(f"     loaded: {msk.shape}, cols sample: {list(msk.columns)[:10]}")
    msk_braf = msk[msk["Hugo_Symbol"]=="BRAF"]
    log(f"     MSK thyroid BRAF rows: {len(msk_braf)}")
    if "HGVSp_Short" in msk_braf.columns:
        log(f"     MSK BRAF HGVSp_Short:")
        log(msk_braf["HGVSp_Short"].value_counts().head(15).to_string())
    # save for later analysis
    msk.to_csv(OUT/"v3_thyroid_mskcc_2016_mutations.tsv.gz", sep="\t", index=False, compression="gzip")
    log(f"     saved to {OUT}/v3_thyroid_mskcc_2016_mutations.tsv.gz")

# also fetch clinical (for survival)
log("\n   fetching thyroid_mskcc_2016 clinical...")
clin = safe_fetch("https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/thyroid_mskcc_2016/data_clinical_sample.txt", 60)
patient_clin = safe_fetch("https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/thyroid_mskcc_2016/data_clinical_patient.txt", 60)
if clin:
    cs = pd.read_csv(io.StringIO(clin), sep="\t", low_memory=False, comment="#")
    log(f"     sample clinical: {cs.shape}, cols: {list(cs.columns)[:10]}")
    cs.to_csv(OUT/"v3_thyroid_mskcc_2016_clinical_sample.tsv", sep="\t", index=False)
if patient_clin:
    cp = pd.read_csv(io.StringIO(patient_clin), sep="\t", low_memory=False, comment="#")
    log(f"     patient clinical: {cp.shape}, cols: {list(cp.columns)[:15]}")
    cp.to_csv(OUT/"v3_thyroid_mskcc_2016_clinical_patient.tsv", sep="\t", index=False)
    if "OS_STATUS" in cp.columns or "OS_MONTHS" in cp.columns:
        log(f"     OS_STATUS distribution:")
        for c in ["OS_STATUS","OS_MONTHS","SUBTYPE","HISTOLOGY","CANCER_TYPE_DETAILED"]:
            if c in cp.columns:
                log(f"       {c}: {cp[c].value_counts().head(5).to_dict()}")

# ============================================================
# 3. Combined per-sample BRAF+TERT table (TCGA + MSK)
# ============================================================
log("\n" + "="*60)
log("[3] Combined BRAF+TERT analysis across cohorts")

combined = []

# TCGA-THCA
for _, r in per_pt.iterrows():
    pid = r["patient_id"]
    sm_row = sm[sm["patient_id"]==pid]
    tert = "wildtype"
    age = np.nan; os_ev = np.nan; os_days = np.nan; hist = ""
    if len(sm_row):
        tert = sm_row.iloc[0].get("tert_promoter_integrated","wildtype") or "wildtype"
        age = sm_row.iloc[0].get("age", np.nan)
        os_ev = sm_row.iloc[0].get("os_event", np.nan)
        os_days = sm_row.iloc[0].get("os_days", np.nan)
        hist = sm_row.iloc[0].get("histology_subtype","")
    combined.append({"cohort":"TCGA-THCA", "patient_id":pid, "histology":hist,
                     "braf_class":r["braf_class"], "tert":tert,
                     "age":age, "os_event":os_ev, "os_days":os_days})
# Add WT BRAF patients from sample master too (RAS, TripleNeg)
braf_pts = set(per_pt["patient_id"])
for _, r in sm.iterrows():
    pid = r["patient_id"]
    if pid in braf_pts: continue
    if r.get("driver_anchor")=="RAS":
        bc = "RAS_only"
    elif r.get("quad_group")=="D_triple_negative":
        bc = "TripleNeg"
    else:
        continue
    combined.append({"cohort":"TCGA-THCA", "patient_id":pid, "histology":r.get("histology_subtype",""),
                     "braf_class":bc, "tert":r.get("tert_promoter_integrated","wildtype") or "wildtype",
                     "age":r.get("age",np.nan), "os_event":r.get("os_event",np.nan),
                     "os_days":r.get("os_days",np.nan)})

# MSK
if data and clin:
    msk_pts = set(msk_braf["Tumor_Sample_Barcode"].astype(str).tolist())
    # all samples in MSK (BRAF + WT)
    all_msk = msk["Tumor_Sample_Barcode"].dropna().unique() if data else []
    msk_braf_pp = msk_braf.drop_duplicates("Tumor_Sample_Barcode").copy()
    msk_braf_pp["braf_class"] = msk_braf_pp["HGVSp_Short"].apply(classify)
    msk_tert = msk[msk["Hugo_Symbol"]=="TERT"]
    log(f"   MSK TERT rows: {len(msk_tert)}")
    msk_tert_set = set(msk_tert["Tumor_Sample_Barcode"].astype(str).tolist())
    # external known TERT for thyroid_mskcc_2016
    msk_tert_ext = ext[ext["cohort_short"]=="thyroid_mskcc_2016"]["sample_barcode"].tolist()
    msk_tert_set |= set(msk_tert_ext)
    log(f"   MSK TERT samples (mutation table + ext promoter records): {len(msk_tert_set)}")
    # clinical join
    cs_lookup = cs.set_index("SAMPLE_ID") if clin and "SAMPLE_ID" in cs.columns else pd.DataFrame()
    cp_lookup = cp.set_index("PATIENT_ID") if patient_clin and "PATIENT_ID" in cp.columns else pd.DataFrame()
    for _, r in msk_braf_pp.iterrows():
        sid = r["Tumor_Sample_Barcode"]
        tert = "mutated" if sid in msk_tert_set else "wildtype"
        hist = ""; age = np.nan; os_ev = np.nan; os_mo = np.nan
        if len(cs_lookup) and sid in cs_lookup.index:
            hist = cs_lookup.loc[sid].get("CANCER_TYPE_DETAILED", cs_lookup.loc[sid].get("HISTOLOGY",""))
            pid = cs_lookup.loc[sid].get("PATIENT_ID","")
            if pid and len(cp_lookup) and pid in cp_lookup.index:
                age = pd.to_numeric(cp_lookup.loc[pid].get("AGE",np.nan), errors="coerce")
                os_status = str(cp_lookup.loc[pid].get("OS_STATUS",""))
                os_ev = 1 if "DECEASED" in os_status.upper() or "1:" in os_status else (0 if os_status else np.nan)
                os_mo = pd.to_numeric(cp_lookup.loc[pid].get("OS_MONTHS",np.nan), errors="coerce")
        combined.append({"cohort":"MSK-thyroid-2016", "patient_id":sid,
                         "histology": str(hist), "braf_class":r["braf_class"], "tert":tert,
                         "age":age, "os_event":os_ev,
                         "os_days": os_mo*30.44 if pd.notna(os_mo) else np.nan})

cdf = pd.DataFrame(combined)
log(f"\n   COMBINED table: {cdf.shape}")
log(f"   cohort × braf_class × tert breakdown:")
ct = pd.crosstab([cdf["cohort"], cdf["braf_class"]], cdf["tert"])
log(ct.to_string())
cdf.to_csv(OUT/"v3_external_combined_BRAF_TERT.tsv", sep="\t", index=False)
log(f"\n   wrote {OUT}/v3_external_combined_BRAF_TERT.tsv")

# ============================================================
# 4. Survival analysis on combined
# ============================================================
log("\n" + "="*60)
log("[4] Combined BRAF×TERT survival")
sv = cdf.dropna(subset=["os_days","os_event"]).copy()
sv["os_event"] = sv["os_event"].astype(float)
log(f"   n with OS: {len(sv)}, events: {int(sv['os_event'].sum())}")

# Combined V600E + TERT analysis
v600 = sv[sv["braf_class"]=="V600E"]
log(f"\n   V600E: n={len(v600)}, ev={int(v600['os_event'].sum())}")
v600_t = v600[v600["tert"]=="mutated"]
v600_w = v600[v600["tert"]=="wildtype"]
log(f"     V600E+TERT+: n={len(v600_t)}, ev={int(v600_t['os_event'].sum())}")
log(f"     V600E+TERT-: n={len(v600_w)}, ev={int(v600_w['os_event'].sum())}")
try:
    from lifelines.statistics import logrank_test
    if len(v600_t)>=3 and len(v600_w)>=3:
        r = logrank_test(v600_t["os_days"], v600_w["os_days"], v600_t["os_event"], v600_w["os_event"])
        log(f"     logrank V600E+TERT vs V600E-only: chi2={r.test_statistic:.2f}, p={r.p_value:.2e}")
except ImportError:
    log("     (lifelines not available)")

# K601 cohort
k = sv[sv["braf_class"]=="K601E"]
log(f"\n   K601E (combined cohorts): n={len(k)}, ev={int(k['os_event'].sum())}")
log(f"     K601E+TERT+: n={int((k['tert']=='mutated').sum())}")
log(f"     K601E+TERT-: n={int((k['tert']=='wildtype').sum())}")

# RAS+TERT
r = sv[sv["braf_class"]=="RAS_only"]
log(f"\n   RAS_only: n={len(r)}, ev={int(r['os_event'].sum())}")
r_t = r[r["tert"]=="mutated"]
r_w = r[r["tert"]=="wildtype"]
log(f"     RAS+TERT+: n={len(r_t)}, ev={int(r_t['os_event'].sum())}")
log(f"     RAS+TERT-: n={len(r_w)}, ev={int(r_w['os_event'].sum())}")

# write log
(OUT/"v3_external_K601_TERT_log.md").write_text("\n".join(LOG))
log(f"\n   log saved to {OUT}/v3_external_K601_TERT_log.md")
