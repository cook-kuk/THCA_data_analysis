"""Build the joined TCGA master table for the NC 6-figure builds.

Joins:
  d4p2 master (DM, BRAF/RAS, age, sig_*, hashi_*)
  tcga_thca_clinical_extended (vital, OS days, stage, age, gender)
  r5_2 HM450 mean β per sample (8 genes)
  cbio_sv_thca structural variants (fusion partners)
  dm_calls_thpa_tcga_gdc (alt DM call provenance)

Outputs: master_tcga.tsv with sample_short, patient_id, DM, drivers, fusion class, OS, mean β.
"""
import pandas as pd, numpy as np, json, sys, os

ROOT = "/home/seungho/personal/THCA_data_analysis/project/results"
OUT = ROOT + "/manuscript_v8_nc_main"
os.makedirs(OUT, exist_ok=True)

# ---- master DM + clinical + drivers ----
d4 = pd.read_csv(f"{ROOT}/d4p2_tcga_hashimoto_signature/tcga_with_clinical_mutations.tsv", sep="\t", index_col=0)
d4.index.name = "sample_id"
d4["patient_id"] = d4.index.str[:12]
d4["sample_short"] = d4.index.str[:12]
print(f"d4p2 master: {len(d4)} rows, DM counts: {d4['DM'].value_counts().to_dict()}")

# ---- clinical extended ----
cli = pd.read_csv(f"{ROOT}/tables/tcga_thca_clinical_extended.tsv", sep="\t")
cli["patient_id"] = cli["case_id"]
print(f"clinical: {len(cli)} rows; columns: {list(cli.columns)[:10]}")

# ---- HM450 mean β ----
hm = pd.read_csv(f"{ROOT}/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv", sep="\t")
hm = hm.rename(columns={"sample_short": "patient_id"})
print(f"HM450: {len(hm)} samples, mean β range [{hm['mean_8g_beta'].min():.3f}, {hm['mean_8g_beta'].max():.3f}]")

# ---- Fusion (cbio_sv) ----
sv = pd.read_csv(f"{ROOT}/audit_2026_04_30/round3/cbio_sv_thca.tsv", sep="\t")
sv["patient_id"] = sv["patientId"]
sv["partner_a"] = sv["site1HugoSymbol"].astype(str)
sv["partner_b"] = sv["site2HugoSymbol"].astype(str)
sv["fusion_pair"] = sv["partner_a"] + "-" + sv["partner_b"]

def classify(row):
    a, b = row["partner_a"], row["partner_b"]
    for kinase in ["RET", "NTRK1", "NTRK2", "NTRK3", "ALK"]:
        if kinase in (a, b): return kinase if kinase.startswith("NTRK") and kinase != "NTRK" else ("NTRK" if kinase.startswith("NTRK") else kinase)
    if "BRAF" in (a, b): return "BRAF"
    if "PAX8" in (a, b) and "PPARG" in (a, b): return "PAX8-PPARG"
    return "other"
sv["fusion_class"] = sv.apply(classify, axis=1)

# Per-patient: any fusion? top fusion class
fus_per_pt = sv.groupby("patient_id").agg(
    fusion_any=("fusion_class", lambda x: True),
    fusion_class=("fusion_class", lambda x: x.value_counts().index[0] if len(x) else "none"),
    fusion_partners=("fusion_pair", lambda x: ";".join(sorted(set(x))[:5]))
).reset_index()
print(f"SV: {len(sv)} events, {fus_per_pt['patient_id'].nunique()} unique patients")
print(f"Fusion class counts: {fus_per_pt['fusion_class'].value_counts().to_dict()}")

# ---- Merge ----
master = d4.reset_index().merge(cli[["case_id", "vital_status", "os_event", "os_days", "stage", "tumor_size_mm", "gender"]], left_on="patient_id", right_on="case_id", how="left")
master = master.merge(hm, on="patient_id", how="left")
master = master.merge(fus_per_pt, on="patient_id", how="left")
master["fusion_any"] = master["fusion_any"].fillna(False)
master["fusion_class"] = master["fusion_class"].fillna("none")

# Driver class derived
def driver_class(r):
    if r.get("has_braf_v600e", 0) > 0: return "BRAF_V600E"
    if r.get("has_ras_mut", 0) > 0: return "RAS"
    if r["fusion_class"] in ("RET", "NTRK1", "NTRK2", "NTRK3", "ALK", "BRAF"): return f"FUSION_{r['fusion_class']}"
    return "DRIVER_NEG"
master["driver_class"] = master.apply(driver_class, axis=1)

# ATA tier proxy (stage-based; intermediate = T3 or N1 proxy via stage III/IV)
def ata_tier(s):
    if pd.isna(s): return "unknown"
    s = str(s)
    if "I" == s.strip() or "Stage I" == s.strip(): return "low"
    if "Stage II" in s and "IV" not in s and "III" not in s: return "low"
    if "Stage III" in s: return "intermediate"
    if "Stage IV" in s: return "high"
    return "intermediate"
master["ata_tier_proxy"] = master["stage"].apply(ata_tier)

# Save
out_tsv = f"{OUT}/master_tcga.tsv"
master.to_csv(out_tsv, sep="\t", index=False)
print(f"\nMASTER written: {out_tsv}  shape={master.shape}")
print("\nDM × driver_class crosstab:")
print(pd.crosstab(master["DM"], master["driver_class"]))
print("\nDM × fusion_any crosstab:")
print(pd.crosstab(master["DM"], master["fusion_any"]))
print(f"\nOS available: {master['os_days'].notna().sum()} samples; events: {master['os_event'].sum() if master['os_event'].notna().any() else 'n/a'}")
print(f"HM450 mean β available: {master['mean_8g_beta'].notna().sum()} samples")
