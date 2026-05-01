#!/usr/bin/env python3
"""v17 TERT recovery v3 — extended analyses on the 36 TERT+ TCGA-THCA cohort.

Adds five things the v2 sprint did not produce:
  1. BRAF V600E vs K601E split (re-parse MAFs)
  2. BRAF + TERT double-hit survival interaction
  3. Per-quad-group Kaplan-Meier
  4. TERT x histology / stage / age interaction tests
  5. Interactive plotly KM + 2x2 figure
"""
from __future__ import annotations
import json, gzip, os, sys
from pathlib import Path
import pandas as pd, numpy as np

OUT = Path("/opt/thyroid-dash/project/results/v17_tert_recovery/v3")
FIG = Path("/opt/thyroid-dash/project/reports/html/figs_interactive/v17")
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

# ============================================================
# 1. BRAF V600 vs K601 split
# ============================================================
MAF_DIR = Path("/data/thca/data_raw/gdc/TCGA-THCA/mutation")
def parse_maf_for_braf():
    rows = []
    for f in sorted(MAF_DIR.glob("*.maf.gz")):
        try:
            df = pd.read_csv(f, sep="\t", compression="gzip", comment="#", low_memory=False)
            sub = df[df["Hugo_Symbol"]=="BRAF"][["Tumor_Sample_Barcode","HGVSp_Short","Variant_Classification"]]
            if len(sub):
                sub = sub.copy()
                sub["maf_file"] = f.name
                rows.append(sub)
        except Exception as e:
            print(f"skip {f.name}: {e}", file=sys.stderr)
    if not rows:
        return pd.DataFrame()
    out = pd.concat(rows, ignore_index=True)
    out["patient_id"] = out["Tumor_Sample_Barcode"].astype(str).str[:12]
    return out

print("[1] parsing MAFs for BRAF (this can take a minute) ...")
braf_mut = parse_maf_for_braf()
print(f"    BRAF mutation rows: {len(braf_mut)}, patients: {braf_mut['patient_id'].nunique()}")

# Classify BRAF position per patient.
# GDC MAFs use the newer BRAF transcript NM_001374258 (BRAF-201) which shifts
# all amino acid numbers by +40 vs the legacy clinical numbering (NM_004333).
# So clinical V600 = GDC V640, K601 = K641, G469 = G509, etc. Map both.
def position(hgvsp):
    s = str(hgvsp)
    # canonical hotspots in either numbering
    if "V600E" in s or "V640E" in s: return "V600E"
    if "V600K" in s or "V640K" in s: return "V600K"
    if "V600" in s or "V640" in s:   return "V600_other"
    if "K601" in s or "K641" in s:   return "K601E_RAS_like"
    # G-loop RAS-like (less common, behaves like RAS)
    if any(t in s for t in ["G469","G464","G466","G509","G504","G506"]):
        return "G_loop_RAS_like"
    return "other_BRAF"

braf_mut["braf_class"] = braf_mut["HGVSp_Short"].apply(position)
print(f"    BRAF class breakdown:")
print(braf_mut.drop_duplicates(["patient_id","braf_class"])["braf_class"].value_counts().to_string())

per_patient_braf = (braf_mut.sort_values("braf_class")
                            .drop_duplicates("patient_id", keep="first")
                            [["patient_id","braf_class","HGVSp_Short"]])
per_patient_braf.to_csv(OUT/"v3_braf_position_per_patient.tsv", sep="\t", index=False)

# ============================================================
# 2. Merge with sample master + TERT integrated
# ============================================================
sm = pd.read_csv("/opt/thyroid-dash/project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv", sep="\t")
sm["patient_id"] = sm["tcga_short"].astype(str)

merged = sm.merge(per_patient_braf, on="patient_id", how="left")
merged["braf_class_filled"] = merged["braf_class"].fillna("WT")

# Tert integrated already in sm: tert_promoter_integrated
merged["tert"] = merged["tert_promoter_integrated"].fillna("wildtype")
merged["tert_bin"] = (merged["tert"]=="mutated").astype(int)

print(f"\n[2] merged shape: {merged.shape}")

# Quint group: BRAF_V600 / BRAF_K601 / RAS / TripleNeg, each ± TERT
def quint(row):
    bc = row["braf_class_filled"]
    qg_lower = str(row["quad_group"]).lower()
    t = "TERT+" if row["tert_bin"]==1 else "TERT-"
    if "braf" in qg_lower:
        if bc == "V600E": return f"BRAF_V600E·{t}"
        if bc == "K601E_RAS_like":  return f"BRAF_K601E·{t}"
        if bc == "V600K": return f"BRAF_V600K·{t}"
        if bc.startswith("V600") or bc == "G_loop_RAS_like": return f"BRAF_other·{t}"
        return f"BRAF_unspecified·{t}"
    if "ras" in qg_lower: return f"RAS·{t}"
    if "triple" in qg_lower: return f"TripleNeg·{t}"
    return f"unknown·{t}"

merged["quint_group"] = merged.apply(quint, axis=1)
print("\n[3] quint_group distribution:")
print(merged["quint_group"].value_counts().to_string())

# ============================================================
# 3. Per-quad-group survival (KM curves) — lifelines if available, else manual
# ============================================================
print("\n[4] survival analysis ...")
have_surv = merged.dropna(subset=["os_days","os_event"])
print(f"    n with OS: {len(have_surv)}, events: {int(have_surv['os_event'].sum())}")

try:
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import multivariate_logrank_test, logrank_test
    HAS_LL = True
except ImportError:
    HAS_LL = False
    print("    lifelines not installed - using manual logrank")

# 4-group logrank
groups_4 = have_surv["quad_group"].fillna("unknown")
if HAS_LL:
    res = multivariate_logrank_test(have_surv["os_days"], groups_4, have_surv["os_event"])
    print(f"    4-group logrank: chi2={res.test_statistic:.2f}, p={res.p_value:.2e}")
    surv4_p = float(res.p_value)
else:
    surv4_p = float("nan")

# BRAF+TERT vs BRAF alone (within BRAF-only patients)
braf_only = have_surv[have_surv["quad_group"]=="A_braf_only"]
braf_only_tert_pos = braf_only[braf_only["tert_bin"]==1]
braf_only_tert_neg = braf_only[braf_only["tert_bin"]==0]
print(f"\n    BRAF-only: TERT+ n={len(braf_only_tert_pos)} ev={int(braf_only_tert_pos['os_event'].sum())} | TERT- n={len(braf_only_tert_neg)} ev={int(braf_only_tert_neg['os_event'].sum())}")

if HAS_LL and len(braf_only_tert_pos)>=3 and len(braf_only_tert_neg)>=3:
    r = logrank_test(braf_only_tert_pos["os_days"], braf_only_tert_neg["os_days"],
                     braf_only_tert_pos["os_event"], braf_only_tert_neg["os_event"])
    print(f"    BRAF+TERT vs BRAF-only logrank: p={r.p_value:.2e}, chi2={r.test_statistic:.2f}")

# ============================================================
# 4. TERT x histology / stage / age tests
# ============================================================
print("\n[5] TERT × clinical interactions:")
from scipy import stats

# histology
ct_hist = pd.crosstab(merged["histology_subtype"], merged["tert"])
chi2, ph, _, _ = stats.chi2_contingency(ct_hist)
print(f"    histology x TERT: chi2={chi2:.2f}, p={ph:.2e}")
print(ct_hist)

# stage
if "stage" in merged.columns:
    stage_simple = merged["stage"].astype(str).str.replace("Stage ","").str[0]
    stage_simple = stage_simple.replace("n","unknown").replace("u","unknown")
    ct_st = pd.crosstab(stage_simple, merged["tert"])
    chi2s, ps, _, _ = stats.chi2_contingency(ct_st)
    print(f"\n    stage x TERT: chi2={chi2s:.2f}, p={ps:.2e}")
    print(ct_st)

# age
if "age" in merged.columns:
    age_pos = merged[merged["tert_bin"]==1]["age"].dropna()
    age_neg = merged[merged["tert_bin"]==0]["age"].dropna()
    u, pa = stats.mannwhitneyu(age_pos, age_neg, alternative="two-sided")
    print(f"\n    age TERT+ ({age_pos.median():.1f}) vs TERT- ({age_neg.median():.1f}): U={u:.0f}, p={pa:.2e}")

# tds16 / rai_score
for col in ["tds16_score_v17","rai_score_v17","tds_score","dedifferentiation_proxy_score"]:
    if col in merged.columns:
        v_pos = merged[merged["tert_bin"]==1][col].dropna()
        v_neg = merged[merged["tert_bin"]==0][col].dropna()
        if len(v_pos)>=5 and len(v_neg)>=5:
            u, p = stats.mannwhitneyu(v_pos, v_neg, alternative="two-sided")
            print(f"    {col}: TERT+ {v_pos.median():.3f} vs TERT- {v_neg.median():.3f}, p={p:.2e}")

# Save merged
merged_out = merged[["patient_id","quad_group","quint_group","braf_class_filled","tert","tert_bin",
                     "histology_subtype","stage","age","os_event","os_days",
                     "tds16_score_v17","rai_score_v17"]]
merged_out.to_csv(OUT/"v3_merged_quint.tsv", sep="\t", index=False)

# ============================================================
# 5. Plotly figures
# ============================================================
print("\n[6] generating figures ...")
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Figure 1: KM curves per quad group + per quint group
def km_curve(times, events):
    """Return (t, S(t)) Kaplan-Meier estimate."""
    df = pd.DataFrame({"t":times, "e":events}).sort_values("t").reset_index(drop=True)
    n = len(df)
    s = 1.0
    ts, ss = [0.0], [1.0]
    at_risk = n
    for _, r in df.iterrows():
        if r["e"]==1 and at_risk>0:
            s *= (1 - 1/at_risk)
        ts.append(r["t"]); ss.append(s)
        at_risk -= 1
    return ts, ss

fig = make_subplots(rows=1, cols=2, subplot_titles=("Quad-group KM", "BRAF×TERT 2×2 KM"),
                    specs=[[{}, {}]])

palette4 = {"A_braf_only":"#F5A623","B_ras_only":"#7ccfcd","D_triple_negative":"#9b59b6"}
for g, color in palette4.items():
    sub = have_surv[have_surv["quad_group"]==g]
    if len(sub)<3: continue
    t, s = km_curve(sub["os_days"].values, sub["os_event"].values)
    fig.add_trace(go.Scatter(x=t, y=s, mode="lines", name=f"{g} (n={len(sub)})",
                             line=dict(color=color, width=2), legendgroup="quad"), row=1, col=1)

palette_quint = {
    "BRAF_V600E·TERT-":"#F5A623",
    "BRAF_V600E·TERT+":"#c24c4c",
    "BRAF_K601·TERT-":"#e8b86b",
    "RAS·TERT-":"#7ccfcd",
    "RAS·TERT+":"#3498db",
    "TripleNeg·TERT-":"#9b59b6",
    "TripleNeg·TERT+":"#6c2ba0",
}
for g, color in palette_quint.items():
    sub = have_surv[have_surv["quint_group"]==g]
    if len(sub)<3: continue
    t, s = km_curve(sub["os_days"].values, sub["os_event"].values)
    dash = "dash" if "TERT+" in g else "solid"
    fig.add_trace(go.Scatter(x=t, y=s, mode="lines",
                             name=f"{g} (n={len(sub)}, ev={int(sub['os_event'].sum())})",
                             line=dict(color=color, width=2, dash=dash), legendgroup="quint"),
                  row=1, col=2)

fig.update_xaxes(title_text="Days", row=1, col=1, gridcolor="rgba(255,255,255,0.05)")
fig.update_xaxes(title_text="Days", row=1, col=2, gridcolor="rgba(255,255,255,0.05)")
fig.update_yaxes(title_text="OS probability", row=1, col=1, range=[0.5,1.02], gridcolor="rgba(255,255,255,0.05)")
fig.update_yaxes(range=[0.5,1.02], row=1, col=2, gridcolor="rgba(255,255,255,0.05)")
fig.update_layout(template="plotly_dark",
                  title=f"v17 TERT v3 — Quad/Quint group survival (TCGA-THCA, n={len(have_surv)}, ev={int(have_surv['os_event'].sum())})",
                  paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                  font=dict(color="#EAEAEA"), height=520, margin=dict(t=80,l=60,r=20,b=60))

fig.write_html(FIG/"v17_tert_v3_KM.html", include_plotlyjs="cdn", full_html=True)
print(f"    wrote {FIG}/v17_tert_v3_KM.html")

# Figure 2: BRAF V600 vs K601 vs other position breakdown
braf_dist = (per_patient_braf["braf_class"].value_counts()
             .reindex(["V600E","V600K","K601E_RAS_like","V600_other","G_loop_RAS_like","other_BRAF"]).fillna(0))
fig2 = go.Figure(go.Bar(x=braf_dist.index, y=braf_dist.values,
                        text=braf_dist.values.astype(int), textposition="outside",
                        marker=dict(color=["#F5A623","#c24c4c","#e8b86b","#7ccfcd","#9b59b6"])))
fig2.update_layout(template="plotly_dark",
                   title=f"BRAF mutation position breakdown (TCGA-THCA, n_patients={int(braf_dist.sum())})",
                   yaxis=dict(title="patients", gridcolor="rgba(255,255,255,0.05)"),
                   paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                   font=dict(color="#EAEAEA"), height=420, margin=dict(t=80,l=60,r=20,b=60))
fig2.write_html(FIG/"v17_tert_v3_braf_positions.html", include_plotlyjs="cdn", full_html=True)
print(f"    wrote {FIG}/v17_tert_v3_braf_positions.html")

# ============================================================
# Summary
# ============================================================
summary = {
    "n_TCGA_THCA": int(len(merged)),
    "n_TERT_pos": int(merged["tert_bin"].sum()),
    "n_OS_avail": int(len(have_surv)),
    "n_OS_events": int(have_surv["os_event"].sum()),
    "BRAF_position_breakdown": dict(zip(braf_dist.index.tolist(), braf_dist.values.astype(int).tolist())),
    "quint_group_counts": merged["quint_group"].value_counts().to_dict(),
}
(OUT/"v3_summary.json").write_text(json.dumps(summary, indent=2))
print(f"\n[7] DONE — summary at {OUT}/v3_summary.json")
print(json.dumps(summary, indent=2))
