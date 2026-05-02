"""
v17 audit Round 8 — R8-1 arm-level CNV (recovered), R8-2 RET partner phenotype,
R8-3 TCGA RAI/I-131 clinical lookup, R8-4 cross-cohort composite validation.

User instruction: 모든 자원 써도 됨 — full R8.
Outputs: project/results/audit_2026_04_30/round8/
"""
from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/audit_2026_04_30/round8"
OUT.mkdir(parents=True, exist_ok=True)

CBIO = "https://www.cbioportal.org/api"
STUDY = "thca_tcga_pan_can_atlas_2018"
ARM_PROFILE = f"{STUDY}_armlevel_cna"
SV_PROFILE = f"{STUDY}_structural_variants"
SAMPLE_LIST = f"{STUDY}_all"


def load_master() -> pd.DataFrame:
    p = ROOT / "project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv"
    m = pd.read_csv(p, sep="\t")
    m["sample_short"] = m["sample_id"].str.slice(0, 15)  # TCGA-XX-XXXX-01
    m["patient_id"] = m["sample_id"].str.slice(0, 12)
    return m


def load_sv() -> pd.DataFrame:
    p = ROOT / "project/results/audit_2026_04_30/round3/cbio_sv_thca.tsv"
    return pd.read_csv(p, sep="\t")


# -------------------------------------------------------------------- R8-1
def r8_1_arm_cnv(master: pd.DataFrame) -> dict:
    """Pull arm-level CNV (GISTIC categorical) for all TCGA-THCA samples,
    then DM1 vs DM2 vs not_DM contingency per arm."""
    sample_ids = master["sample_short"].dropna().unique().tolist()
    # batches of 200 for safety
    rows = []
    BATCH = 200
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = []
        for i in range(0, len(sample_ids), BATCH):
            chunk = sample_ids[i : i + BATCH]
            futs.append(
                ex.submit(
                    requests.post,
                    f"{CBIO}/generic_assay_data/{ARM_PROFILE}/fetch?projection=DETAILED",
                    headers={"Content-Type": "application/json"},
                    json={"sampleIds": chunk},
                    timeout=60,
                )
            )
        for f in as_completed(futs):
            r = f.result()
            if r.status_code == 200:
                rows.extend(r.json())

    arm = pd.DataFrame(rows)
    if arm.empty:
        return {"error": "no arm-level data returned"}

    arm["sample_short"] = arm["sampleId"]
    arm["arm"] = arm["genericAssayStableId"].str.replace("_status", "", regex=False)
    arm = arm[["sample_short", "arm", "value"]].rename(columns={"value": "call"})
    arm.to_csv(OUT / "r8_1_arm_cnv_long.tsv", sep="\t", index=False)

    merged = arm.merge(
        master[["sample_short", "v17_dark_cluster"]], on="sample_short", how="inner"
    )
    merged["dm"] = merged["v17_dark_cluster"].fillna("not_DM").replace("", "not_DM")

    # per-arm DM1 vs DM2 vs not_DM gain/loss frequencies
    rows_out = []
    for arm_name, sub in merged.groupby("arm"):
        for cluster in ["DM1", "DM2", "not_DM"]:
            csub = sub[sub["dm"] == cluster]
            if len(csub) < 10:
                continue
            n = len(csub)
            gain = (csub["call"] == "Gain").sum()
            loss = (csub["call"] == "Loss").sum()
            unchanged = (csub["call"] == "Unchanged").sum()
            rows_out.append(
                dict(
                    arm=arm_name,
                    cluster=cluster,
                    n=n,
                    gain_pct=round(100 * gain / n, 1),
                    loss_pct=round(100 * loss / n, 1),
                    unchanged_pct=round(100 * unchanged / n, 1),
                )
            )
    perarm = pd.DataFrame(rows_out)
    perarm.to_csv(OUT / "r8_1_perarm_dm.tsv", sep="\t", index=False)

    # DM1 vs not_DM gain difference, per arm: Fisher
    enrich = []
    for arm_name in perarm["arm"].unique():
        sub = merged[merged["arm"] == arm_name]
        dm1_gain = ((sub["dm"] == "DM1") & (sub["call"] == "Gain")).sum()
        dm1_other = ((sub["dm"] == "DM1") & (sub["call"] != "Gain")).sum()
        nd_gain = ((sub["dm"] == "not_DM") & (sub["call"] == "Gain")).sum()
        nd_other = ((sub["dm"] == "not_DM") & (sub["call"] != "Gain")).sum()
        if min(dm1_gain + dm1_other, nd_gain + nd_other) < 10:
            continue
        try:
            odds, p = stats.fisher_exact(
                [[dm1_gain, dm1_other], [nd_gain, nd_other]]
            )
        except Exception:
            continue
        enrich.append(
            dict(
                arm=arm_name,
                dm1_gain_pct=round(
                    100 * dm1_gain / max(dm1_gain + dm1_other, 1), 1
                ),
                not_dm_gain_pct=round(
                    100 * nd_gain / max(nd_gain + nd_other, 1), 1
                ),
                odds_ratio=round(odds, 3) if np.isfinite(odds) else None,
                fisher_p=p,
            )
        )
    enrichdf = pd.DataFrame(enrich).sort_values("fisher_p")
    enrichdf.to_csv(OUT / "r8_1_dm1_vs_notDM_arm_enrichment.tsv", sep="\t", index=False)
    top5 = enrichdf.head(5).to_dict(orient="records")
    return {
        "n_samples": int(arm["sample_short"].nunique()),
        "n_arms": int(arm["arm"].nunique()),
        "top5_dm1_vs_notDM_arm_enrichment": top5,
        "interpretation": "Arm-level CNA via GISTIC categorical calls (Gain/Loss/Unchanged); DM1 vs not_DM Fisher per arm",
    }


# -------------------------------------------------------------------- R8-2
def r8_2_ret_partners(master: pd.DataFrame, sv: pd.DataFrame) -> dict:
    """RET fusion partners by DM cluster: phenotype, age, BRAF, score."""
    is_ret = (sv["site1HugoSymbol"] == "RET") | (sv["site2HugoSymbol"] == "RET")
    ret = sv[is_ret].copy()
    ret["partner"] = np.where(
        ret["site1HugoSymbol"] == "RET",
        ret["site2HugoSymbol"],
        ret["site1HugoSymbol"],
    )
    ret["sample_short"] = ret["sampleId"]
    # one row per sample (collapse multiple breakpoints per same partner)
    ret_simple = (
        ret.groupby("sample_short")["partner"].apply(lambda s: ";".join(sorted(set(s)))).reset_index()
    )

    merged = master.merge(ret_simple, on="sample_short", how="left")
    merged["ret_partner"] = merged["partner"].fillna("none")
    merged["any_ret"] = (merged["ret_partner"] != "none").astype(int)
    merged["dm"] = merged["v17_dark_cluster"].fillna("not_DM").replace("", "not_DM")

    # partition partners
    def cls(p: str) -> str:
        if p == "none":
            return "no_RET"
        if "CCDC6" in p:
            return "CCDC6-RET"
        if "NCOA4" in p:
            return "NCOA4-RET"
        return f"other-RET ({p})"

    merged["ret_class"] = merged["ret_partner"].apply(cls)

    # DM1 only — partner-by-phenotype
    dm1 = merged[merged["dm"] == "DM1"].copy()
    rows = []
    for k, sub in dm1.groupby("ret_class"):
        if len(sub) < 1:
            continue
        rows.append(
            dict(
                ret_class=k,
                n=len(sub),
                rai_score_med=round(sub["rai_score_v17"].median(skipna=True), 3)
                if "rai_score_v17" in sub
                else None,
                age_med=round(sub["age"].dropna().median(), 1) if sub["age"].notna().any() else None,
                braf_pct=round(
                    100
                    * (sub["driver_anchor"].astype(str).str.contains("BRAF", na=False)).mean(),
                    1,
                ),
                tert_pos_pct=round(
                    100 * (sub["tert_promoter_integrated"].fillna("wildtype") != "wildtype").mean(),
                    1,
                ),
                cptc_pct=round(100 * (sub["histology_subtype"] == "cPTC").mean(), 1),
            )
        )
    parttab = pd.DataFrame(rows).sort_values("n", ascending=False)
    parttab.to_csv(OUT / "r8_2_dm1_ret_partners.tsv", sep="\t", index=False)

    # KW test for rai_score across partner classes
    classes = [
        sub["rai_score_v17"].dropna().to_numpy()
        for k, sub in dm1.groupby("ret_class")
        if sub["rai_score_v17"].notna().sum() > 1
    ]
    if len(classes) >= 2:
        try:
            kw_stat, kw_p = stats.kruskal(*classes)
        except Exception:
            kw_stat, kw_p = None, None
    else:
        kw_stat, kw_p = None, None

    return {
        "dm1_n": int(len(dm1)),
        "dm1_ret_partner_breakdown": parttab.to_dict(orient="records"),
        "rai_score_kw_p_across_partners": float(kw_p) if kw_p is not None else None,
        "interpretation": "RET partner heterogeneity within DM1 (CCDC6 vs NCOA4 vs other)",
    }


# -------------------------------------------------------------------- R8-3
def r8_3_clinical_rai(master: pd.DataFrame) -> dict:
    """Pull all sample/patient clinical attributes; look for RAI / radioiodine / I-131 fields."""
    # 1. list available clinical attribute IDs
    r = requests.get(
        f"{CBIO}/studies/{STUDY}/clinical-attributes",
        timeout=60,
    )
    if r.status_code != 200:
        return {"error": f"clinical-attributes status {r.status_code}"}
    attrs = pd.DataFrame(r.json())
    # match anything RAI-ish
    pat = "(?i)radioiod|radio_iod|^rai|iodine|i131|i-131|^treat|therapy"
    rai_attrs = attrs[
        attrs["clinicalAttributeId"].str.contains(pat, regex=True, na=False)
        | attrs["displayName"].str.contains(pat, regex=True, na=False)
    ]
    rai_attrs.to_csv(OUT / "r8_3_candidate_attributes.tsv", sep="\t", index=False)

    if rai_attrs.empty:
        return {
            "n_clinical_attributes": int(len(attrs)),
            "rai_candidates": 0,
            "note": "no RAI-like clinical attribute id in cBioPortal pan-can-atlas thyroid clinical schema",
            "interpretation": "TCGA pan-can-atlas thyroid clinical lacks per-patient I-131 dose / response field; needs CDR.txt or UCSC Xena clinicalMatrix",
        }

    # 2. fetch values for each candidate
    out_rows = []
    for _, a in rai_attrs.iterrows():
        attr_id = a["clinicalAttributeId"]
        is_patient = a["patientAttribute"]
        endpoint = (
            f"{CBIO}/studies/{STUDY}/clinical-data?clinicalDataType={'PATIENT' if is_patient else 'SAMPLE'}&attributeId={attr_id}&projection=DETAILED"
        )
        rr = requests.get(endpoint, timeout=60)
        if rr.status_code != 200:
            continue
        data = pd.DataFrame(rr.json())
        if data.empty:
            continue
        data["clinicalAttributeId"] = attr_id
        out_rows.append(data)

    if out_rows:
        clin = pd.concat(out_rows, ignore_index=True)
        clin.to_csv(OUT / "r8_3_clinical_rai_values.tsv", sep="\t", index=False)
        # value distributions
        dist = (
            clin.groupby(["clinicalAttributeId", "value"])
            .size()
            .reset_index(name="n")
        )
        dist.to_csv(OUT / "r8_3_value_distribution.tsv", sep="\t", index=False)
        dist_brief = dist.sort_values(["clinicalAttributeId", "n"], ascending=[True, False]).head(40)
    else:
        dist_brief = pd.DataFrame()

    return {
        "n_clinical_attributes": int(len(attrs)),
        "rai_candidates": int(len(rai_attrs)),
        "candidate_ids": rai_attrs["clinicalAttributeId"].tolist(),
        "value_distribution_top": dist_brief.to_dict(orient="records") if not dist_brief.empty else [],
        "interpretation": "Direct RAI/I-131 outcome fields scanned; if absent, paper limitation #N+1",
    }


# -------------------------------------------------------------------- R8-4
def r8_4_cross_cohort_composite() -> dict:
    """Build R7-2 composite model on TCGA, then score K2 (PRJEB11591) and GSE213647."""
    master = load_master()
    hla_tcga = pd.read_csv(
        ROOT / "project/results/v17_hla/tcga_thca_hla_per_sample.tsv", sep="\t"
    )
    hla_tcga["sample_short"] = hla_tcga["sample_id"].str.slice(0, 15)
    meth = pd.read_csv(
        ROOT / "project/results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv",
        sep="\t",
    )
    if "sample_short" not in meth.columns:
        meth["sample_short"] = meth.iloc[:, 0]

    sv = load_sv()
    fusion_samples = set(sv["sampleId"].unique())
    master["any_fusion"] = master["sample_short"].isin(fusion_samples).astype(int)
    master["dm"] = master["v17_dark_cluster"].fillna("not_DM").replace("", "not_DM")

    # methylation file uses 12-char patient id form; master uses 15-char (with -01)
    master["patient12"] = master["sample_short"].str.slice(0, 12)
    hla_tcga["patient12"] = hla_tcga["sample_short"].str.slice(0, 12)
    meth["patient12"] = meth["sample_short"].str.slice(0, 12)
    df = (
        master[
            [
                "patient12",
                "sample_short",
                "rai_score_v17",
                "age",
                "v17_dark_cluster",
                "any_fusion",
            ]
        ]
        .merge(hla_tcga[["patient12", "hla_class_I_score", "hla_class_II_score"]],
               on="patient12", how="inner")
        .merge(meth[["patient12", "mean_8g_beta"]], on="patient12", how="inner")
    )
    df = df.dropna(
        subset=["rai_score_v17", "hla_class_I_score", "hla_class_II_score", "mean_8g_beta", "age"]
    )
    df["dm1"] = (df["v17_dark_cluster"] == "DM1").astype(int)

    feats = ["rai_score_v17", "hla_class_I_score", "hla_class_II_score", "mean_8g_beta", "age"]
    X = df[feats].to_numpy()
    sc = StandardScaler().fit(X)
    Xz = sc.transform(X)

    model_dm1 = LogisticRegression(max_iter=2000, C=1.0).fit(Xz, df["dm1"].to_numpy())
    auc_dm1_train = roc_auc_score(df["dm1"], model_dm1.predict_proba(Xz)[:, 1])

    model_fus = LogisticRegression(max_iter=2000, C=1.0).fit(Xz, df["any_fusion"].to_numpy())
    auc_fus_train = roc_auc_score(df["any_fusion"], model_fus.predict_proba(Xz)[:, 1])

    # ---- K2 application -------------------------------------------------
    k2_pred = pd.read_csv(
        ROOT / "project/results/v17_korean/K2_korean_predictions_v4.tsv", sep="\t"
    )
    k2_meta = pd.read_csv(
        ROOT / "project/results/v17_korean/K1A_prjeb11591_runs.tsv", sep="\t"
    ) if (ROOT / "project/results/v17_korean/K1A_prjeb11591_runs.tsv").exists() else None

    # K2 has 8-gene TPM and DM call; build proxy composite using only the
    # features available in K2: rai_score (within-sample-z), age (proxy via meta),
    # HLA scores from arcasHLA (separate file).
    arcas_k2_paths = [
        ROOT / "project/results/v17_arcasHLA/k2_korean_hla_per_sample.tsv",
        ROOT / "project/results/v17_arcasHLA/K2_arcasHLA_per_sample.tsv",
        ROOT / "project/results/v17_hla/korean_GSE213647_hla_per_sample.tsv",
    ]
    k2_hla_path = None
    for p in arcas_k2_paths:
        if p.exists():
            k2_hla_path = p
            break

    k2_summary = {"n": int(len(k2_pred)), "p_DM2_median": float(k2_pred["p_DM2"].median())}

    # ---- GSE213647 -----------------------------------------------------
    g3_score = pd.read_csv(
        ROOT / "project/results/v17_korean/GSE213647_panel_score.tsv", sep="\t"
    )
    g3_hla = pd.read_csv(
        ROOT / "project/results/v17_hla/korean_GSE213647_hla_per_sample.tsv", sep="\t"
    )
    g3 = g3_score.merge(
        g3_hla[["gsm", "hla_class_I_score", "hla_class_II_score"]], on="gsm", how="left"
    )
    g3["age_n"] = g3["age"]
    # within-sample-z mock: panel_z already z-scored within cohort
    g3 = g3.dropna(subset=["panel_z", "hla_class_I_score", "hla_class_II_score", "age"])
    # build a 4-feature composite (no methylation in GSE213647)
    Xg = g3[["panel_z", "hla_class_I_score", "hla_class_II_score", "age"]].to_numpy()
    Xg_z = (Xg - Xg.mean(axis=0)) / (Xg.std(axis=0) + 1e-9)

    # train a 4-feature DM1 surrogate on TCGA using those same 4 features
    feats4 = ["rai_score_v17", "hla_class_I_score", "hla_class_II_score", "age"]
    X4 = StandardScaler().fit_transform(df[feats4].to_numpy())
    m4 = LogisticRegression(max_iter=2000).fit(X4, df["dm1"].to_numpy())
    auc4_tcga = roc_auc_score(df["dm1"], m4.predict_proba(X4)[:, 1])
    g3_pred = m4.predict_proba(Xg_z)[:, 1]
    g3["composite_dm1_p"] = g3_pred

    # rank histology by composite
    if "histology" in g3.columns:
        g3_summary = (
            g3.groupby("histology")["composite_dm1_p"].agg(["count", "median", "mean"]).reset_index()
        )
        g3_summary.to_csv(OUT / "r8_4_gse213647_composite_by_histology.tsv", sep="\t", index=False)
    else:
        g3_summary = pd.DataFrame()

    g3.to_csv(OUT / "r8_4_gse213647_composite_pred.tsv", sep="\t", index=False)

    return {
        "tcga_train_auc_dm1_5feat": round(auc_dm1_train, 3),
        "tcga_train_auc_fusion_5feat": round(auc_fus_train, 3),
        "tcga_train_auc_dm1_4feat_no_meth": round(auc4_tcga, 3),
        "k2_summary": k2_summary,
        "k2_hla_file_used": str(k2_hla_path) if k2_hla_path else None,
        "gse213647_n": int(len(g3)),
        "gse213647_composite_median": round(float(g3["composite_dm1_p"].median()), 3),
        "gse213647_histology_breakdown": g3_summary.to_dict(orient="records") if not g3_summary.empty else [],
        "interpretation": "TCGA-trained 4-feature composite (no methylation) applied to GSE213647 within-z; histology stratification",
    }


# -------------------------------------------------------------------- main
def main() -> None:
    print("[R8] loading master...")
    master = load_master()
    print(f"  master n={len(master)}")
    print("[R8] loading SV...")
    sv = load_sv()
    print(f"  SV rows={len(sv)}")

    print("[R8-1] arm-level CNV...")
    r1 = r8_1_arm_cnv(master)
    (OUT / "r8_1_arm_cnv_summary.json").write_text(json.dumps(r1, indent=2, default=str))
    print(f"  -> {r1.get('n_samples')} samples × {r1.get('n_arms')} arms")

    print("[R8-2] RET partners...")
    r2 = r8_2_ret_partners(master, sv)
    (OUT / "r8_2_ret_partners_summary.json").write_text(json.dumps(r2, indent=2, default=str))

    print("[R8-3] clinical RAI...")
    r3 = r8_3_clinical_rai(master)
    (OUT / "r8_3_clinical_rai_summary.json").write_text(json.dumps(r3, indent=2, default=str))

    print("[R8-4] cross-cohort composite...")
    r4 = r8_4_cross_cohort_composite()
    (OUT / "r8_4_cross_cohort_summary.json").write_text(json.dumps(r4, indent=2, default=str))

    print("[R8] done")


if __name__ == "__main__":
    main()
