"""
Dark Matter Phase 1 viability test — steps 1, 2, 3, 6 (TCGA-THCA only).

Step 1: Dark Matter cohort = BRAF V600E negative AND NRAS/HRAS/KRAS hotspot negative.
        Compute %.
Step 2: Within DM, examine 8-gene cluster (v17_dark_cluster) sizes + OS Cox HR.
Step 3: Within DM, count DICER1/EIF1AX/fusion/TERT per cluster, Fisher exact.
Step 6: Xing-axis (BRAF+TERT) 4-group → DM rescue rate.

NOT executing GSE33630 / K2 / Bundang / sc — mutation calls absent (see data_inventory).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from scipy import stats

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project")
OUT = ROOT / "results" / "dark_matter_phase1"
OUT.mkdir(parents=True, exist_ok=True)


def short12(s):
    return s[:12] if isinstance(s, str) else s


# ============================================================================
# Load source tables
# ============================================================================
mut = pd.read_csv(ROOT / "results/tables/tcga_thca_mutation_groups.tsv", sep="\t")
clin = pd.read_csv(ROOT / "results/tables/tcga_thca_clinical_extended.tsv", sep="\t")
master = pd.read_csv(ROOT / "results/v17/tables/sample_master_v17_full.tsv", sep="\t", low_memory=False)
tert = pd.read_csv(ROOT / "results/v17_tert_recovery/v2/FINAL_tert_status_integrated.tsv", sep="\t")

# Filter master to TCGA-THCA primary tumors
master_tcga = master[(master["dataset"] == "TCGA-THCA") & (master["normal_vs_tumor"] == "tumor")].copy()
master_tcga["tcga_short"] = master_tcga["sample_id"].apply(short12)
mut["tcga_short"] = mut["sample_id"].apply(short12)
clin["tcga_short"] = clin["sample_id"].apply(short12)

# TERT calls (12-char short barcode)
tert_set = set(tert.loc[tert["status"] == "mutated", "tcga_short"])

# Build merged DM-analysis dataframe — INNER on mutation_groups (definitional)
df = mut[["tcga_short", "has_braf_v600e", "has_ras_mut", "mutation_group"]].copy()
_keep = [c for c in ["tcga_short", "v17_dark_cluster", "driver_anchor_v17", "mutation_genes",
                     "tds16_score_v17", "rai_score_v17", "age", "sex"] if c in master_tcga.columns]
df = df.merge(master_tcga[_keep], on="tcga_short", how="left")
df = df.merge(
    clin[["tcga_short", "os_event", "os_days", "stage"]],
    on="tcga_short", how="left"
)
df["tert_pos"] = df["tcga_short"].isin(tert_set)
df["dm_status"] = (~df["has_braf_v600e"]) & (~df["has_ras_mut"])

df.to_csv(OUT / "tcga_dark_matter_master.tsv", sep="\t", index=False)

# ============================================================================
# STEP 1 — Dark Matter %
# ============================================================================
n_total = len(df)
n_dm = int(df["dm_status"].sum())
n_braf = int(df["has_braf_v600e"].sum())
n_ras = int(df["has_ras_mut"].sum())
pct_dm = n_dm / n_total * 100

step1 = {
    "n_total_tcga_thca_mutation_groups": n_total,
    "n_braf_v600e": n_braf,
    "n_ras_hotspot": n_ras,
    "n_dark_matter": n_dm,
    "pct_dark_matter": round(pct_dm, 2),
    "pass_threshold_5pct": pct_dm >= 5,
    "pass_expected_15to25pct": 15 <= pct_dm <= 35,
}
print(f"\n[STEP 1] Dark Matter % = {pct_dm:.2f}% ({n_dm}/{n_total})  PASS={step1['pass_threshold_5pct']}")

# ============================================================================
# STEP 2 — Within-DM cluster + survival HR
# ============================================================================
dm = df[df["dm_status"]].copy()
n_dm_with_cluster = dm["v17_dark_cluster"].notna().sum()
cluster_counts = dm["v17_dark_cluster"].value_counts(dropna=False).to_dict()

# Survival data — filter to samples with both cluster and survival
surv = dm[dm["v17_dark_cluster"].isin(["DM1", "DM2"]) & dm["os_days"].notna() & dm["os_event"].notna()].copy()
surv["cluster_dm2"] = (surv["v17_dark_cluster"] == "DM2").astype(int)
n_surv = len(surv)
n_events = int(surv["os_event"].sum())

cox_result = {}
if n_events >= 3 and surv["cluster_dm2"].nunique() == 2:
    cph = CoxPHFitter()
    try:
        cph.fit(surv[["os_days", "os_event", "cluster_dm2"]], duration_col="os_days", event_col="os_event")
        s = cph.summary.iloc[0]
        cox_result = {
            "HR": float(s["exp(coef)"]),
            "HR_lower": float(s["exp(coef) lower 95%"]),
            "HR_upper": float(s["exp(coef) upper 95%"]),
            "p_value": float(s["p"]),
            "n": n_surv,
            "events": n_events,
            "ref": "DM1",
            "tested": "DM2 vs DM1",
        }
    except Exception as e:
        cox_result = {"error": str(e), "n": n_surv, "events": n_events}
else:
    cox_result = {"error": "insufficient events or single cluster", "n": n_surv, "events": n_events}

step2 = {
    "n_dm_total": int(len(dm)),
    "n_dm_with_cluster": int(n_dm_with_cluster),
    "cluster_counts_in_dm": {str(k): int(v) for k, v in cluster_counts.items()},
    "cox_dm2_vs_dm1_OS": cox_result,
    "decision_rule": "HR>2 with p<0.05 → publishable",
    "pass": cox_result.get("HR", 0) > 2 and cox_result.get("p_value", 1) < 0.05,
}
print(f"[STEP 2] DM-cluster: DM1 n={cluster_counts.get('DM1', 0)}, DM2 n={cluster_counts.get('DM2', 0)}; "
      f"OS events={n_events}; HR(DM2/DM1)={cox_result.get('HR', 'NA')}")

# ============================================================================
# STEP 3 — Alt-driver / fusion / TERT enrichment per cluster within DM
# ============================================================================
dm_c = dm[dm["v17_dark_cluster"].isin(["DM1", "DM2"])].copy()

def enrich_test(gene_label):
    """Fisher exact: count gene_label-positive samples per cluster among DM."""
    if gene_label == "TERT":
        dm_c["pos"] = dm_c["tert_pos"]
    elif gene_label == "anchor_alt":
        # any non-BRAF/RAS/unknown driver_anchor_v17
        alt = ["DICER1_EIF1AX_PPM1D", "RET", "NTRK_fusion", "RET_fusion", "ALK_fusion", "TP53", "PAX8PPARG"]
        dm_c["pos"] = dm_c["driver_anchor_v17"].isin(alt)
    else:
        dm_c["pos"] = dm_c["driver_anchor_v17"] == gene_label
    a = ((dm_c["v17_dark_cluster"] == "DM1") & dm_c["pos"]).sum()
    b = ((dm_c["v17_dark_cluster"] == "DM1") & ~dm_c["pos"]).sum()
    c = ((dm_c["v17_dark_cluster"] == "DM2") & dm_c["pos"]).sum()
    d = ((dm_c["v17_dark_cluster"] == "DM2") & ~dm_c["pos"]).sum()
    table = [[int(a), int(b)], [int(c), int(d)]]
    if a + c == 0:
        return {"label": gene_label, "DM1_pos": 0, "DM1_neg": int(b), "DM2_pos": 0, "DM2_neg": int(d), "p": None, "OR": None}
    odds, p = stats.fisher_exact(table)
    return {
        "label": gene_label,
        "DM1_pos": int(a), "DM1_neg": int(b),
        "DM2_pos": int(c), "DM2_neg": int(d),
        "DM1_freq": round(a / (a + b) if (a + b) else 0, 3),
        "DM2_freq": round(c / (c + d) if (c + d) else 0, 3),
        "OR": round(float(odds), 3),
        "p": round(float(p), 4),
    }

step3_rows = [
    enrich_test("DICER1_EIF1AX_PPM1D"),
    enrich_test("RET"),
    enrich_test("RET_fusion"),
    enrich_test("NTRK_fusion"),
    enrich_test("ALK_fusion"),
    enrich_test("TP53"),
    enrich_test("PAX8PPARG"),
    enrich_test("anchor_alt"),
    enrich_test("TERT"),
]
step3_df = pd.DataFrame(step3_rows)
step3_df.to_csv(OUT / "step3_dm_cluster_alt_driver_enrichment.tsv", sep="\t", index=False)
print(f"[STEP 3] Alt-driver enrichment table:\n{step3_df.to_string(index=False)}")

# ============================================================================
# STEP 6 — Xing axis (BRAF + TERT) 4-group rescue rate
# ============================================================================
# Xing 4-group: BRAF-/TERT-, BRAF+/TERT-, BRAF-/TERT+, BRAF+/TERT+
xing_full = df.copy()
xing_full["xing_group"] = "BRAF-/TERT-"
xing_full.loc[xing_full["has_braf_v600e"] & ~xing_full["tert_pos"], "xing_group"] = "BRAF+/TERT-"
xing_full.loc[~xing_full["has_braf_v600e"] & xing_full["tert_pos"], "xing_group"] = "BRAF-/TERT+"
xing_full.loc[xing_full["has_braf_v600e"] & xing_full["tert_pos"], "xing_group"] = "BRAF+/TERT+"

xing_full["xing_low_risk"] = xing_full["xing_group"].isin(["BRAF-/TERT-", "BRAF+/TERT-"])
# DM cluster as risk: assume DM2 = high-risk per existing v17 work convention if HR>1 in step 2
high_risk_cluster = "DM2" if cox_result.get("HR", 1) > 1 else "DM1"
xing_full["dm_high_risk"] = xing_full["v17_dark_cluster"] == high_risk_cluster

# Rescue: Xing low-risk but DM high-risk
xing_full["rescued"] = xing_full["xing_low_risk"] & xing_full["dm_high_risk"]
n_xing_low = int(xing_full["xing_low_risk"].sum())
n_rescued = int(xing_full["rescued"].sum())
rescue_rate = n_rescued / n_xing_low * 100 if n_xing_low else 0

xing_dist = xing_full["xing_group"].value_counts().to_dict()
xing_full[["tcga_short", "xing_group", "v17_dark_cluster", "rescued"]].to_csv(
    OUT / "step6_xing_rescue.tsv", sep="\t", index=False
)
step6 = {
    "xing_group_distribution": {k: int(v) for k, v in xing_dist.items()},
    "high_risk_cluster_assigned_to": high_risk_cluster,
    "n_xing_low_risk": n_xing_low,
    "n_rescued_by_dm_cluster": n_rescued,
    "rescue_rate_pct": round(rescue_rate, 2),
    "decision_rule": "rescue >5% → clinical impact claim",
    "pass": rescue_rate > 5,
}
print(f"[STEP 6] Xing low-risk={n_xing_low}, rescued by DM cluster={n_rescued} ({rescue_rate:.2f}%)")

# ============================================================================
# Final go/nogo
# ============================================================================
overall_pass = step1["pass_threshold_5pct"] and step2.get("pass", False) and step6.get("pass", False)
summary = {
    "step1_DM_pct": step1,
    "step2_within_DM_cluster_HR": step2,
    "step3_alt_driver_enrichment": step3_rows,
    "step6_xing_rescue": step6,
    "GO_NOGO": "GO" if overall_pass else "PARTIAL_OR_NOGO",
    "notes": "Steps 4 (sc) and 5 (trajectory) require GSE241184 reanalysis — deferred to next sub-task.",
}
(OUT / "step1_to_6_summary.json").write_text(json.dumps(summary, indent=2))
print(f"\n=== OVERALL: {summary['GO_NOGO']} ===")
print(f"Step1 pass: {step1['pass_threshold_5pct']}")
print(f"Step2 pass: {step2.get('pass', False)}")
print(f"Step3 (informational, no hard pass): see table")
print(f"Step6 pass: {step6.get('pass', False)}")
