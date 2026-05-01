"""Re-run Step 2 (within-DM cluster survival HR) with PFI from TCGA-CDR (Liu 2018)."""
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from scipy import stats

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project")
OUT = ROOT / "results" / "dark_matter_phase1"

# Load CDR — Liu 2018 has all 33 TCGA cohorts in one sheet
cdr_xlsx = OUT / "tcga_cdr.xlsx"
xl = pd.ExcelFile(cdr_xlsx)
print("Sheets:", xl.sheet_names)
cdr = pd.read_excel(cdr_xlsx, sheet_name=0)
print("Columns:", list(cdr.columns)[:30])
print("Shape:", cdr.shape)

# Filter THCA
thca = cdr[cdr["type"] == "THCA"].copy()
print(f"\nTHCA rows: {len(thca)}")
endpoints = ["OS", "OS.time", "DSS", "DSS.time", "DFI", "DFI.time", "PFI", "PFI.time"]
for e in endpoints:
    if e in thca.columns:
        nn = thca[e].notna().sum()
        if e.endswith(".time"):
            print(f"  {e}: non-null {nn}")
        else:
            evt = (thca[e] == 1).sum()
            print(f"  {e}: events {evt} / non-null {nn}")

# bcr_patient_barcode is the 12-char short form
thca["tcga_short"] = thca["bcr_patient_barcode"]

# Load our DM master built earlier
dm_master = pd.read_csv(OUT / "tcga_dark_matter_master.tsv", sep="\t")

# Merge in PFI / DFI / DSS
keep_cdr = ["tcga_short", "OS", "OS.time", "PFI", "PFI.time", "DFI", "DFI.time", "DSS", "DSS.time", "ajcc_pathologic_tumor_stage"]
keep_cdr = [c for c in keep_cdr if c in thca.columns]
df = dm_master.merge(thca[keep_cdr], on="tcga_short", how="left")
print(f"\nMerged DM-master with CDR: {df.shape}")
print(f"DM patients with PFI: {df.loc[df['dm_status'], 'PFI'].notna().sum()}")
print(f"DM PFI events: {(df.loc[df['dm_status'], 'PFI'] == 1).sum()}")
print(f"All TCGA-THCA PFI events (any driver): {(df['PFI'] == 1).sum()} / {df['PFI'].notna().sum()}")

# DM-only Cox HR (DM2 vs DM1) on PFI
dm = df[df["dm_status"] & df["v17_dark_cluster"].isin(["DM1", "DM2"])].copy()
results = {}

for endpoint, time_col in [("PFI", "PFI.time"), ("DFI", "DFI.time"), ("DSS", "DSS.time"), ("OS", "OS.time")]:
    if endpoint not in dm.columns:
        continue
    sub = dm[[endpoint, time_col, "v17_dark_cluster"]].dropna()
    sub["cluster_dm2"] = (sub["v17_dark_cluster"] == "DM2").astype(int)
    n = len(sub)
    events = int(sub[endpoint].sum())
    print(f"\n[{endpoint}] DM-only n={n}, events={events}")
    if events < 3 or sub["cluster_dm2"].nunique() < 2:
        results[endpoint] = {"n": n, "events": events, "skipped": "insufficient"}
        continue
    cph = CoxPHFitter()
    try:
        cph.fit(sub.rename(columns={endpoint: "event", time_col: "time"})[["time", "event", "cluster_dm2"]],
                duration_col="time", event_col="event")
        s = cph.summary.iloc[0]
        results[endpoint] = {
            "n": n, "events": events,
            "HR": float(s["exp(coef)"]),
            "CI_low": float(s["exp(coef) lower 95%"]),
            "CI_high": float(s["exp(coef) upper 95%"]),
            "p": float(s["p"]),
        }
        print(f"  HR(DM2/DM1) = {results[endpoint]['HR']:.2f} "
              f"[{results[endpoint]['CI_low']:.2f}-{results[endpoint]['CI_high']:.2f}], p={results[endpoint]['p']:.4f}")
    except Exception as e:
        results[endpoint] = {"n": n, "events": events, "error": str(e)}
        print(f"  ERROR: {e}")

# Also: DM2 vs everything-else (BRAF+/RAS+/DM1) — broader survival claim
print("\n\n=== Whole-cohort comparison: DM2 vs rest ===")
all_with_cluster = df[df["v17_dark_cluster"].isin(["DM1", "DM2"])].copy()
all_with_cluster["dm2_in_dm"] = (all_with_cluster["dm_status"] & (all_with_cluster["v17_dark_cluster"] == "DM2")).astype(int)
for endpoint, time_col in [("PFI", "PFI.time"), ("DFI", "DFI.time"), ("DSS", "DSS.time"), ("OS", "OS.time")]:
    if endpoint not in all_with_cluster.columns:
        continue
    sub = all_with_cluster[[endpoint, time_col, "dm2_in_dm"]].dropna()
    n = len(sub); events = int(sub[endpoint].sum())
    print(f"[{endpoint}] n={n}, events={events}")
    if events < 3 or sub["dm2_in_dm"].nunique() < 2:
        continue
    cph = CoxPHFitter()
    cph.fit(sub.rename(columns={endpoint: "event", time_col: "time"})[["time", "event", "dm2_in_dm"]],
            duration_col="time", event_col="event")
    s = cph.summary.iloc[0]
    print(f"  HR(DM2-in-DM vs rest) = {s['exp(coef)']:.2f} "
          f"[{s['exp(coef) lower 95%']:.2f}-{s['exp(coef) upper 95%']:.2f}], p={s['p']:.4f}")

import json
(OUT / "step2_pfi_results.json").write_text(json.dumps(results, indent=2))
df.to_csv(OUT / "tcga_dm_master_with_pfi.tsv", sep="\t", index=False)
print(f"\nWrote step2_pfi_results.json and tcga_dm_master_with_pfi.tsv")
