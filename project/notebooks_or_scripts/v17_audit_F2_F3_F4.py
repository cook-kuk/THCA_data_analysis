"""F2 (TCGA Hashimoto-like × DM cluster) + F3 (DM1 sub mechanism) + F4 (cBioPortal API fusion)."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, fisher_exact, chi2_contingency
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings("ignore")

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/audit_2026_04_30/round3"
OUT.mkdir(parents=True, exist_ok=True)

# Load core data
sm = pd.read_csv(ROOT / "project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv",
                  sep="\t", low_memory=False)
sm = sm[(sm["dataset"] == "TCGA-THCA") & (sm["normal_vs_tumor"] == "tumor")].copy()
sm["tcga_short"] = sm["sample_id"].str.slice(0, 12)
sm["age_n"] = pd.to_numeric(sm["age"], errors="coerce")

xing = pd.read_csv(ROOT / "project/results/dark_matter_phase1/step6_xing_rescue.tsv", sep="\t")
xing["dm_cluster"] = xing["v17_dark_cluster"].fillna("not_DM")
if "v17_dark_cluster" in sm.columns:
    sm = sm.drop(columns=["v17_dark_cluster"])
sm = sm.merge(xing[["tcga_short", "xing_group", "dm_cluster"]], on="tcga_short", how="left")
sm["dm_cluster"] = sm["dm_cluster"].fillna("not_DM")

immune = pd.read_csv(ROOT / "project/results/tables/immune_signatures_tcga.tsv", sep="\t")
sm = sm.merge(immune, on="sample_id", how="left", suffixes=("", "_imm"))

hla_tc = pd.read_csv(ROOT / "project/results/v17_hla/tcga_thca_hla_per_sample.tsv", sep="\t")
hla_tc["tcga_short_h"] = hla_tc["sample_id"].str.slice(0, 12)
sm = sm.merge(hla_tc[["tcga_short_h", "hla_class_I_score", "hla_class_II_score"]],
              left_on="tcga_short", right_on="tcga_short_h", how="left")

# =================================================================
# F2 — TCGA Hashimoto-like × DM cluster
# =================================================================
print("\n=== F2 — TCGA Hashimoto-like × DM cluster ===\n")

# Hashimoto-like proxy: B_cell + IFN_gamma_response score
sm["hashi_proxy"] = pd.to_numeric(sm["B_cell"], errors="coerce") * 0.5 + \
                     pd.to_numeric(sm["IFN_gamma_response"], errors="coerce") * 0.5
proxy_valid = sm.dropna(subset=["hashi_proxy", "dm_cluster"])
print(f"TCGA samples with hashi_proxy + DM cluster: {len(proxy_valid)}")

# Bimodality via GMM-like: simple percentile
hashi_threshold = proxy_valid["hashi_proxy"].quantile(0.80)  # top 20% as "Hashimoto-like"
sm["hashi_like_bin"] = (sm["hashi_proxy"] >= hashi_threshold).astype("Int64")
proxy_valid = sm.dropna(subset=["hashi_proxy", "dm_cluster"]).copy()
proxy_valid["hashi_like_bin"] = (proxy_valid["hashi_proxy"] >= hashi_threshold).astype(int)

n_hashi_like = int((proxy_valid["hashi_like_bin"] == 1).sum())
print(f"Hashimoto-like (top 20% B_cell+IFN-γ): {n_hashi_like}")

# Cross-tab with DM cluster
ct = pd.crosstab(proxy_valid["dm_cluster"], proxy_valid["hashi_like_bin"])
ct = ct.reindex(["DM1", "DM2", "not_DM"]).fillna(0).astype(int)
print("DM cluster × Hashimoto-like crosstab:")
print(ct)

# DM2 vs DM1 enrichment (Fisher) — handle 0 cells with continuity correction
if "DM1" in ct.index and "DM2" in ct.index:
    ct_2x2 = ct.loc[["DM1", "DM2"]]
    if ct_2x2.shape == (2, 2):
        # Add 0.5 continuity correction if any zero
        ct_corrected = ct_2x2.values.astype(float) + 0.5
        odds, p = fisher_exact(ct_corrected)
        # Original (no correction) for OR display
        try:
            odds_raw, p_raw = fisher_exact(ct_2x2.values)
        except Exception:
            odds_raw, p_raw = float("inf"), 0.0
        f2_dm12 = {
            "n_DM1": int(ct_2x2.loc["DM1"].sum()),
            "n_DM2": int(ct_2x2.loc["DM2"].sum()),
            "hashi_yes_DM1": int(ct_2x2.loc["DM1", 1]),
            "hashi_yes_DM2": int(ct_2x2.loc["DM2", 1]),
            "hashi_pct_DM1": round(100 * ct_2x2.loc["DM1", 1] / ct_2x2.loc["DM1"].sum(), 1),
            "hashi_pct_DM2": round(100 * ct_2x2.loc["DM2", 1] / ct_2x2.loc["DM2"].sum(), 1),
            "odds_ratio_corrected": round(float(odds), 3),
            "fisher_p_corrected": round(float(p), 6),
            "fisher_p_raw": round(float(p_raw), 6),
            "DIRECTION": "DM1 has MORE Hashimoto-like" if ct_2x2.loc["DM1", 1] > ct_2x2.loc["DM2", 1] else "DM2 has MORE",
        }
    else:
        f2_dm12 = {"error": "shape mismatch"}
else:
    f2_dm12 = None
print(f"\nDM2 vs DM1 enrichment: {json.dumps(f2_dm12, indent=2)}")

# HLA-II residualization: remove Hashimoto-like, recompute Cohen's d
hla_data = sm.dropna(subset=["dm_cluster", "hla_class_II_score", "hashi_proxy"]).copy()
hla_data["hashi_like_bin"] = (hla_data["hashi_proxy"] >= hashi_threshold).astype(int)

# Original (with all)
dm1_all = hla_data[hla_data["dm_cluster"] == "DM1"]["hla_class_II_score"]
dm2_all = hla_data[hla_data["dm_cluster"] == "DM2"]["hla_class_II_score"]
d_all = (dm1_all.mean() - dm2_all.mean()) / np.sqrt((dm1_all.var() + dm2_all.var()) / 2) if len(dm1_all) >= 5 and len(dm2_all) >= 5 else None

# Without Hashimoto-like
no_hashi = hla_data[hla_data["hashi_like_bin"] == 0]
dm1_nh = no_hashi[no_hashi["dm_cluster"] == "DM1"]["hla_class_II_score"]
dm2_nh = no_hashi[no_hashi["dm_cluster"] == "DM2"]["hla_class_II_score"]
d_nh = (dm1_nh.mean() - dm2_nh.mean()) / np.sqrt((dm1_nh.var() + dm2_nh.var()) / 2) if len(dm1_nh) >= 5 and len(dm2_nh) >= 5 else None

print(f"\nHLA-II Cohen's d:")
print(f"  All samples: d={d_all:.3f} (n_DM1={len(dm1_all)}, n_DM2={len(dm2_all)})")
if d_nh is not None:
    print(f"  Excluding Hashimoto-like: d={d_nh:.3f} (n_DM1={len(dm1_nh)}, n_DM2={len(dm2_nh)})")
    print(f"  Δd = {d_all - d_nh:.3f}")
else:
    print(f"  Excluding Hashimoto-like: NA (insufficient n)")

# Same for HLA-I
dm1_all_I = hla_data.dropna(subset=["hla_class_I_score"])
dm1_I = dm1_all_I[dm1_all_I["dm_cluster"] == "DM1"]["hla_class_I_score"]
dm2_I = dm1_all_I[dm1_all_I["dm_cluster"] == "DM2"]["hla_class_I_score"]
d_I_all = (dm1_I.mean() - dm2_I.mean()) / np.sqrt((dm1_I.var() + dm2_I.var()) / 2) if len(dm1_I) >= 5 and len(dm2_I) >= 5 else None
no_hashi_I = no_hashi.dropna(subset=["hla_class_I_score"])
dm1_I_nh = no_hashi_I[no_hashi_I["dm_cluster"] == "DM1"]["hla_class_I_score"]
dm2_I_nh = no_hashi_I[no_hashi_I["dm_cluster"] == "DM2"]["hla_class_I_score"]
d_I_nh = (dm1_I_nh.mean() - dm2_I_nh.mean()) / np.sqrt((dm1_I_nh.var() + dm2_I_nh.var()) / 2) if len(dm1_I_nh) >= 5 and len(dm2_I_nh) >= 5 else None

# Hashimoto-like × Xing rescue
sm_xing = sm.dropna(subset=["xing_group", "dm_cluster", "hashi_proxy"]).copy()
sm_xing["hashi_like_bin"] = (sm_xing["hashi_proxy"] >= hashi_threshold).astype(int)
xing_btb = sm_xing[sm_xing["xing_group"] == "BRAF-/TERT-"]
xing_3way = pd.crosstab([xing_btb["dm_cluster"], xing_btb["hashi_like_bin"]], xing_btb["hashi_like_bin"])
print(f"\nBRAF-/TERT- × DM × Hashi-like 3-way:")
ct_xing = pd.crosstab(xing_btb["dm_cluster"], xing_btb["hashi_like_bin"])
print(ct_xing)

# Hashimoto-like × age
hashi_yes = sm[sm["hashi_proxy"] >= hashi_threshold]["age_n"].dropna()
hashi_no = sm[sm["hashi_proxy"] < hashi_threshold]["age_n"].dropna()
mw_age = mannwhitneyu(hashi_yes, hashi_no, alternative="two-sided")
cd_age = (hashi_yes.mean() - hashi_no.mean()) / np.sqrt((hashi_yes.var() + hashi_no.var()) / 2)
print(f"\nHashi-like vs not Hashi age: median {hashi_yes.median():.1f} vs {hashi_no.median():.1f} (Cohen's d={cd_age:.2f}, MW p={mw_age.pvalue:.2e})")

f2_summary = {
    "method": "TCGA-THCA Hashimoto-like proxy = (B_cell score + IFN_gamma_response) / 2, top-20% threshold",
    "threshold": float(hashi_threshold),
    "n_total": int(len(proxy_valid)),
    "n_hashi_like": int(n_hashi_like),
    "hashi_pct": round(100 * n_hashi_like / len(proxy_valid), 1),
    "DM1_vs_DM2_enrichment": f2_dm12,
    "hla_II_cohens_d": {
        "all_samples": round(float(d_all), 3) if d_all is not None else None,
        "excluding_hashi": round(float(d_nh), 3) if d_nh is not None else None,
        "delta": round(float(d_all - d_nh), 3) if (d_all and d_nh) else None,
        "interpretation": (
            "If delta < 0.3 → HLA-II differential is independent of Hashimoto-overlap. "
            "If delta > 0.5 → significant fraction mediated by autoimmune signature."
        ),
    },
    "hla_I_cohens_d": {
        "all_samples": round(float(d_I_all), 3) if d_I_all is not None else None,
        "excluding_hashi": round(float(d_I_nh), 3) if d_I_nh is not None else None,
        "delta": round(float(d_I_all - d_I_nh), 3) if (d_I_all and d_I_nh) else None,
    },
    "hashi_like_age": {
        "hashi_yes_median": round(float(hashi_yes.median()), 1),
        "hashi_no_median": round(float(hashi_no.median()), 1),
        "cohens_d": round(float(cd_age), 3),
        "mw_p": float(mw_age.pvalue),
    },
    "BRAF-/TERT-_dm_x_hashi": ct_xing.to_dict() if not ct_xing.empty else None,
}
(OUT / "f2_tcga_hashimoto_generalization.json").write_text(json.dumps(f2_summary, indent=2, default=str), encoding="utf-8")

# =================================================================
# F3 — DM1 sub-A vs sub-B mechanism (via available signatures)
# =================================================================
print("\n\n=== F3 — DM1 sub-A vs sub-B mechanism ===\n")

dm1 = sm[sm["dm_cluster"] == "DM1"].copy()
score_cols = ["tds_score", "rai_score_v17", "tds16_score_v17", "dedifferentiation_proxy_score"]
score_cols = [c for c in score_cols if c in dm1.columns]
dm1_scores = dm1[score_cols].apply(pd.to_numeric, errors="coerce").dropna()

scaler = StandardScaler()
sc_data = scaler.fit_transform(dm1_scores)
km = KMeans(n_clusters=2, random_state=42, n_init=10)
sub_labs = km.fit_predict(sc_data)
dm1_with_sub = dm1.loc[dm1_scores.index].copy()
dm1_with_sub["dm1_sub"] = pd.Series(sub_labs, index=dm1_scores.index).map({0: "subA", 1: "subB"})

# Confirm sub-A vs sub-B (sub-A = high RAI, n=72; sub-B = low RAI, n=19)
print(f"Sub-A: {(dm1_with_sub['dm1_sub']=='subA').sum()}, Sub-B: {(dm1_with_sub['dm1_sub']=='subB').sum()}")

f3_results = []

# Immune signatures × sub-cluster
imm_cols = ["T_cell_total", "CD8_cytotoxic", "Treg", "Macrophage_M1", "Macrophage_M2",
             "NK_cell", "B_cell", "Checkpoint_exhaustion", "IFN_gamma_response"]
for col in imm_cols:
    if col not in dm1_with_sub.columns: continue
    a = pd.to_numeric(dm1_with_sub[dm1_with_sub["dm1_sub"] == "subA"][col], errors="coerce").dropna()
    b = pd.to_numeric(dm1_with_sub[dm1_with_sub["dm1_sub"] == "subB"][col], errors="coerce").dropna()
    if len(a) >= 5 and len(b) >= 3:
        mw = mannwhitneyu(a, b, alternative="two-sided")
        cd = (a.mean() - b.mean()) / np.sqrt((a.var() + b.var()) / 2) if (a.var() + b.var()) > 0 else 0
        f3_results.append({"angle": f"immune:{col}", "n_A": int(len(a)), "n_B": int(len(b)),
                            "median_A": round(float(a.median()), 3),
                            "median_B": round(float(b.median()), 3),
                            "cohens_d": round(float(cd), 3),
                            "mw_p": float(mw.pvalue)})

# Hashimoto-like × sub-cluster
dm1_with_sub["hashi_proxy"] = pd.to_numeric(dm1_with_sub["B_cell"], errors="coerce") * 0.5 + \
                                pd.to_numeric(dm1_with_sub["IFN_gamma_response"], errors="coerce") * 0.5
dm1_with_sub["hashi_like_bin"] = (dm1_with_sub["hashi_proxy"] >= hashi_threshold).astype(int)
ct_sub_hashi = pd.crosstab(dm1_with_sub["dm1_sub"], dm1_with_sub["hashi_like_bin"])
if ct_sub_hashi.shape == (2, 2):
    odds, p = fisher_exact(ct_sub_hashi.values + 0.5)
    f3_results.append({"angle": "Hashimoto-like × sub", "n_A": int(ct_sub_hashi.loc["subA"].sum()),
                        "n_B": int(ct_sub_hashi.loc["subB"].sum()),
                        "hashi_pct_A": round(100 * ct_sub_hashi.loc["subA", 1] / ct_sub_hashi.loc["subA"].sum(), 1),
                        "hashi_pct_B": round(100 * ct_sub_hashi.loc["subB", 1] / ct_sub_hashi.loc["subB"].sum(), 1),
                        "odds_ratio": round(float(odds), 3),
                        "p": round(float(p), 4),
                        "test": "Fisher exact"})

# Age × sub-cluster (P7 connection)
a_age = dm1_with_sub[dm1_with_sub["dm1_sub"] == "subA"]["age_n"].dropna()
b_age = dm1_with_sub[dm1_with_sub["dm1_sub"] == "subB"]["age_n"].dropna()
if len(a_age) >= 5 and len(b_age) >= 3:
    mw = mannwhitneyu(a_age, b_age, alternative="two-sided")
    cd = (a_age.mean() - b_age.mean()) / np.sqrt((a_age.var() + b_age.var()) / 2)
    f3_results.append({"angle": "Age × sub", "n_A": int(len(a_age)), "n_B": int(len(b_age)),
                        "median_A": round(float(a_age.median()), 1),
                        "median_B": round(float(b_age.median()), 1),
                        "cohens_d": round(float(cd), 3),
                        "mw_p": float(mw.pvalue)})

# Young-onset (<45) × sub
young_a = (dm1_with_sub[dm1_with_sub["dm1_sub"] == "subA"]["age_n"] < 45).sum()
young_b = (dm1_with_sub[dm1_with_sub["dm1_sub"] == "subB"]["age_n"] < 45).sum()
n_a_age = dm1_with_sub[dm1_with_sub["dm1_sub"] == "subA"]["age_n"].notna().sum()
n_b_age = dm1_with_sub[dm1_with_sub["dm1_sub"] == "subB"]["age_n"].notna().sum()
try:
    odds, p = fisher_exact([[young_a, n_a_age - young_a], [young_b, n_b_age - young_b]])
    f3_results.append({"angle": "Young-onset (<45) × sub", "n_A": int(n_a_age), "n_B": int(n_b_age),
                        "young_pct_A": round(100 * young_a / n_a_age, 1) if n_a_age else None,
                        "young_pct_B": round(100 * young_b / n_b_age, 1) if n_b_age else None,
                        "odds_ratio": round(float(odds), 3),
                        "p": round(float(p), 4)})
except Exception:
    pass

f3_df = pd.DataFrame(f3_results)
f3_df.to_csv(OUT / "f3_dm1_subcluster_mechanism.tsv", sep="\t", index=False)
print(f3_df.to_string(index=False))

# =================================================================
# F4 — cBioPortal API attempt
# =================================================================
print("\n\n=== F4 — cBioPortal API attempt ===\n")
import requests

api_base = "https://www.cbioportal.org/api"
study_id = "thca_tcga_pan_can_atlas_2018"

# Try structural variants endpoint
def try_get(url, params=None, timeout=15):
    try:
        r = requests.get(url, params=params, timeout=timeout, headers={"Accept": "application/json"})
        return r.status_code, (r.json() if r.status_code == 200 else r.text[:200])
    except Exception as e:
        return None, str(e)[:200]

# Test API
code, response = try_get(f"{api_base}/studies/{study_id}")
print(f"Test ping (study info): {code}")
if code == 200:
    print(f"Study: {response.get('name', '')[:80]}")

# Try fetching structural variants — endpoints vary
sv_url = f"{api_base}/structural-variants/fetch"
code2, resp2 = try_get(sv_url, params={
    "molecularProfileId": f"{study_id}_structural_variants",
    "sampleListId": f"{study_id}_all",
    "structuralVariantQueries": "[]",
})

# Try alternative: study profiles
prof_code, prof_resp = try_get(f"{api_base}/studies/{study_id}/molecular-profiles")
print(f"Profiles fetch: {prof_code}")
profiles_with_sv = []
if prof_code == 200:
    for p in prof_resp:
        mt = p.get("molecularAlterationType", "")
        if "STRUCTURAL" in mt or "FUSION" in mt or "MUTATION_EXTENDED" in mt:
            profiles_with_sv.append({"id": p.get("molecularProfileId"),
                                     "name": p.get("name"), "type": mt})
print(f"Mutation/SV profiles: {profiles_with_sv}")

f4_summary = {
    "api_status": code,
    "structural_variant_endpoint_status": code2 if 'code2' in dir() else None,
    "profiles_query_status": prof_code,
    "profiles_with_mut_or_sv": profiles_with_sv,
    "verdict": "API accessible" if code == 200 else "API blocked or down",
    "next_steps": (
        "Try /structural-variants/fetch POST with molecularProfileId + sampleIds. "
        "If Pan-Cancer Atlas study has fusion data, retrieve. "
        "Otherwise pivot to gdc-client manifest download."
    ),
}
(OUT / "f4_cbioportal_api_attempt.json").write_text(json.dumps(f4_summary, indent=2, default=str), encoding="utf-8")
print(json.dumps(f4_summary, indent=2)[:500])

# Try POST structural-variants/fetch
print("\nAttempting POST /structural-variants/fetch...")
try:
    payload = {
        "molecularProfileIds": [f"{study_id}_structural_variants"],
        "sampleIdentifiers": [],
    }
    r = requests.post(f"{api_base}/structural-variants/fetch",
                       json=payload, timeout=30,
                       headers={"Content-Type": "application/json", "Accept": "application/json"})
    print(f"  POST status: {r.status_code}")
    if r.status_code == 200:
        sv_data = r.json()
        print(f"  Returned {len(sv_data)} structural variant records")
        if len(sv_data) > 0:
            sv_df = pd.DataFrame(sv_data)
            sv_df.to_csv(OUT / "f4_tcga_thca_structural_variants.tsv", sep="\t", index=False)
            print(f"  Saved to f4_tcga_thca_structural_variants.tsv")
            print(f"  Columns: {list(sv_df.columns)[:10]}")
    else:
        print(f"  Body: {r.text[:300]}")
except Exception as e:
    print(f"  POST error: {e}")

print("\n\n=== Done ===")
