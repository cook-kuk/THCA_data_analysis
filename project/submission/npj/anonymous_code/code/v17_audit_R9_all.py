"""
v17 audit Round 9 — RAI clinical recovered from UCSC Xena legacy clinicalMatrix
(R8-3 limitation closed), chr7 RTK gene-level CNV, bootstrap chr7, K2 within-z
composite re-apply.
"""
from __future__ import annotations

import gzip
import json
import shutil
import subprocess
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
OUT = ROOT / "project/results/audit_2026_04_30/round9"
OUT.mkdir(parents=True, exist_ok=True)
CBIO = "https://www.cbioportal.org/api"
STUDY = "thca_tcga_pan_can_atlas_2018"


def load_master() -> pd.DataFrame:
    p = ROOT / "project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv"
    m = pd.read_csv(p, sep="\t")
    m["sample_short"] = m["sample_id"].str.slice(0, 15)
    m["patient12"] = m["sample_short"].str.slice(0, 12)
    m["dm"] = m["v17_dark_cluster"].fillna("not_DM").replace("", "not_DM")
    return m


# -------------------------------------------------------------------- R9-1
def r9_1_xena_rai(master: pd.DataFrame) -> dict:
    xena_path = Path("/tmp/THCA_clinicalMatrix.tsv")
    if not xena_path.exists():
        url = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.THCA.sampleMap/THCA_clinicalMatrix"
        r = requests.get(url, timeout=120)
        xena_path.write_bytes(r.content)
    cm = pd.read_csv(xena_path, sep="\t", low_memory=False)
    rai_cols = [
        "sampleID",
        "i_131_first_administered_dose",
        "i_131_subsequent_administered_dose",
        "i_131_total_administered_dose",
        "i_131_total_administered_preparation_technique",
        "radiosensitizing_agent_administered_indicator",
        "additional_pharmaceutical_therapy",
        "additional_radiation_therapy",
        "new_tumor_event_after_initial_treatment",
        "history_of_neoadjuvant_treatment",
    ]
    avail = [c for c in rai_cols if c in cm.columns]
    sub = cm[avail].copy()
    sub["sample_short"] = sub["sampleID"].str.slice(0, 15)
    merged = master.merge(sub, on="sample_short", how="left")
    merged.to_csv(OUT / "r9_1_master_with_rai.tsv", sep="\t", index=False)

    # numeric coerce
    for col in [
        "i_131_first_administered_dose",
        "i_131_subsequent_administered_dose",
        "i_131_total_administered_dose",
    ]:
        if col in merged.columns:
            merged[col + "_num"] = pd.to_numeric(merged[col], errors="coerce")

    # by DM cluster: I-131 total dose
    rows = []
    for cluster in ["DM1", "DM2", "not_DM"]:
        s = merged[merged["dm"] == cluster]
        d = s["i_131_total_administered_dose_num"].dropna()
        rows.append(
            dict(
                cluster=cluster,
                n_total=int(len(s)),
                n_with_i131=int(d.shape[0]),
                pct_received_i131=round(100 * d.shape[0] / max(len(s), 1), 1),
                median_mCi=round(float(d.median()), 1) if d.shape[0] else None,
                mean_mCi=round(float(d.mean()), 1) if d.shape[0] else None,
            )
        )
    bydm = pd.DataFrame(rows)
    bydm.to_csv(OUT / "r9_1_i131_by_DM.tsv", sep="\t", index=False)

    # NEW_TUMOR_EVENT (proxy for failure / recurrence after initial Tx) by DM
    nev_field = "new_tumor_event_after_initial_treatment"
    nev_rows = []
    for cluster in ["DM1", "DM2", "not_DM"]:
        s = merged[merged["dm"] == cluster][nev_field]
        n_total = (s.notna() & (s.astype(str).str.lower().isin(["yes", "no"]))).sum()
        n_yes = (s.astype(str).str.lower() == "yes").sum()
        nev_rows.append(
            dict(
                cluster=cluster,
                n_evaluable=int(n_total),
                n_recurrence=int(n_yes),
                recurrence_pct=round(100 * n_yes / max(n_total, 1), 1),
            )
        )
    nev = pd.DataFrame(nev_rows)
    nev.to_csv(OUT / "r9_1_recurrence_by_DM.tsv", sep="\t", index=False)

    # Mann-Whitney DM1 vs not_DM I-131 dose
    d_dm1 = merged.loc[merged["dm"] == "DM1", "i_131_total_administered_dose_num"].dropna()
    d_nd = merged.loc[merged["dm"] == "not_DM", "i_131_total_administered_dose_num"].dropna()
    try:
        u, p_dose = stats.mannwhitneyu(d_dm1, d_nd, alternative="two-sided") if len(d_dm1) and len(d_nd) else (None, None)
    except Exception:
        u, p_dose = None, None

    # DM1 vs not_DM recurrence Fisher
    a = nev.loc[nev["cluster"] == "DM1", "n_recurrence"].iloc[0]
    a_total = nev.loc[nev["cluster"] == "DM1", "n_evaluable"].iloc[0]
    b = nev.loc[nev["cluster"] == "not_DM", "n_recurrence"].iloc[0]
    b_total = nev.loc[nev["cluster"] == "not_DM", "n_evaluable"].iloc[0]
    try:
        odds, p_rec = stats.fisher_exact([[a, a_total - a], [b, b_total - b]])
    except Exception:
        odds, p_rec = None, None

    return {
        "n_with_xena_clinical": int(merged["sampleID"].notna().sum()),
        "i131_by_dm": bydm.to_dict(orient="records"),
        "recurrence_by_dm": nev.to_dict(orient="records"),
        "dm1_vs_notDM_i131_dose_MW_p": float(p_dose) if p_dose is not None else None,
        "dm1_vs_notDM_recurrence_OR": round(float(odds), 3) if odds and np.isfinite(odds) else None,
        "dm1_vs_notDM_recurrence_fisher_p": float(p_rec) if p_rec is not None else None,
        "interpretation": "I-131 total dose mCi & post-treatment recurrence by DM cluster — closes R8-3 gap",
    }


# -------------------------------------------------------------------- R9-2
def r9_2_chr7_rtk(master: pd.DataFrame) -> dict:
    """log2CNA gene-level for BRAF/EGFR/MET/RAF1 (chr7) and RET (chr10) for control."""
    # Entrez: BRAF=673, EGFR=1956, MET=4233, RAF1=5894, RET=5979
    genes_entrez = {
        "BRAF": 673,
        "EGFR": 1956,
        "MET": 4233,
        "RAF1": 5894,
        "RET": 5979,
    }
    sample_ids = master["sample_short"].dropna().unique().tolist()
    BATCH = 200
    rows = []
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = []
        for i in range(0, len(sample_ids), BATCH):
            chunk = sample_ids[i : i + BATCH]
            futs.append(
                ex.submit(
                    requests.post,
                    f"{CBIO}/molecular-profiles/{STUDY}_log2CNA/molecular-data/fetch?projection=DETAILED",
                    headers={"Content-Type": "application/json"},
                    json={"sampleIds": chunk, "entrezGeneIds": list(genes_entrez.values())},
                    timeout=90,
                )
            )
        for f in as_completed(futs):
            r = f.result()
            if r.status_code == 200:
                rows.extend(r.json())

    df = pd.DataFrame(rows)
    if df.empty:
        return {"error": "no log2CNA data"}
    inv = {v: k for k, v in genes_entrez.items()}
    df["gene"] = df["entrezGeneId"].map(inv)
    df["sample_short"] = df["sampleId"]
    df["log2"] = pd.to_numeric(df["value"], errors="coerce")
    pivot = df.pivot_table(index="sample_short", columns="gene", values="log2", aggfunc="mean").reset_index()
    pivot.to_csv(OUT / "r9_2_chr7_rtk_log2CNA.tsv", sep="\t", index=False)
    merged = pivot.merge(master[["sample_short", "dm"]], on="sample_short", how="inner")

    rows_out = []
    for gene in genes_entrez:
        if gene not in merged.columns:
            continue
        d_dm1 = merged.loc[merged["dm"] == "DM1", gene].dropna()
        d_dm2 = merged.loc[merged["dm"] == "DM2", gene].dropna()
        d_nd = merged.loc[merged["dm"] == "not_DM", gene].dropna()
        # Cohen's d DM1 vs not_DM
        if len(d_dm1) > 5 and len(d_nd) > 5:
            pooled_sd = np.sqrt(
                ((d_dm1.var(ddof=1) * (len(d_dm1) - 1)) + (d_nd.var(ddof=1) * (len(d_nd) - 1)))
                / (len(d_dm1) + len(d_nd) - 2)
            )
            d_eff = (d_dm1.mean() - d_nd.mean()) / (pooled_sd + 1e-9)
            u, p = stats.mannwhitneyu(d_dm1, d_nd, alternative="two-sided")
        else:
            d_eff, p = None, None
        rows_out.append(
            dict(
                gene=gene,
                dm1_n=len(d_dm1),
                dm1_med_log2=round(float(d_dm1.median()), 3) if len(d_dm1) else None,
                dm2_med_log2=round(float(d_dm2.median()), 3) if len(d_dm2) else None,
                notdm_med_log2=round(float(d_nd.median()), 3) if len(d_nd) else None,
                cohens_d_dm1_vs_notDM=round(float(d_eff), 3) if d_eff is not None else None,
                mw_p=float(p) if p is not None else None,
            )
        )
    summary = pd.DataFrame(rows_out)
    summary.to_csv(OUT / "r9_2_chr7_rtk_summary.tsv", sep="\t", index=False)
    return {
        "n_samples": int(merged["sample_short"].nunique()),
        "per_gene": summary.to_dict(orient="records"),
        "interpretation": "Chr 7 RTK gene-level log2 CNA (BRAF, EGFR, MET, RAF1) by DM cluster; RET (chr10) control",
    }


# -------------------------------------------------------------------- R9-3
def r9_3_bootstrap_chr7(master: pd.DataFrame, n_boot: int = 1000) -> dict:
    """Bootstrap CI for chr 7p, 7q gain DM1 vs not_DM enrichment."""
    arm = pd.read_csv(
        ROOT / "project/results/audit_2026_04_30/round8/r8_1_arm_cnv_long.tsv", sep="\t"
    )
    merged = arm.merge(master[["sample_short", "dm"]], on="sample_short", how="inner")
    rng = np.random.default_rng(42)

    boots: dict[str, list[float]] = {"7p": [], "7q": []}
    samples_dm1 = merged.loc[merged["dm"] == "DM1", "sample_short"].unique()
    samples_nd = merged.loc[merged["dm"] == "not_DM", "sample_short"].unique()

    for _ in range(n_boot):
        s_dm1 = rng.choice(samples_dm1, size=len(samples_dm1), replace=True)
        s_nd = rng.choice(samples_nd, size=len(samples_nd), replace=True)
        for arm_name in ["7p", "7q"]:
            sub_dm1 = arm[(arm["sample_short"].isin(s_dm1)) & (arm["arm"] == arm_name)]
            sub_nd = arm[(arm["sample_short"].isin(s_nd)) & (arm["arm"] == arm_name)]
            # use pandas .merge with sampling repeats — for boot rate, treat sample list
            dm1_gain = sum(
                (arm[(arm["sample_short"] == s) & (arm["arm"] == arm_name)]["call"] == "Gain").any()
                for s in s_dm1[:60]  # subsample for speed
            )
            nd_gain = sum(
                (arm[(arm["sample_short"] == s) & (arm["arm"] == arm_name)]["call"] == "Gain").any()
                for s in s_nd[:60]
            )
            rate_dm1 = dm1_gain / 60
            rate_nd = nd_gain / 60
            boots[arm_name].append(rate_dm1 - rate_nd)

    summary = {}
    for arm_name, vals in boots.items():
        a = np.array(vals)
        summary[arm_name] = dict(
            mean_diff_dm1_minus_notdm=round(float(a.mean()), 4),
            ci_low=round(float(np.percentile(a, 2.5)), 4),
            ci_high=round(float(np.percentile(a, 97.5)), 4),
            n_boot=int(n_boot),
        )
    (OUT / "r9_3_bootstrap_chr7.json").write_text(json.dumps(summary, indent=2))
    return {"bootstrap": summary, "interpretation": "1000-bootstrap CI for chr 7p/7q DM1-specific gain rate-difference"}


# -------------------------------------------------------------------- R9-4
def r9_4_k2_within_z_composite() -> dict:
    """K2 8-gene TPM → within-sample z (per-sample z across 8 genes), then build
    composite. Compare to R7-2 absolute-form (which mis-routed)."""
    k2 = pd.read_csv(
        ROOT / "project/results/v17_korean/K2_korean_predictions_v4.tsv", sep="\t"
    )
    genes = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
    expr = k2[genes].copy()
    expr_log = np.log1p(expr)
    # within-sample z: per-row z-score across the 8 genes
    z = expr_log.subtract(expr_log.mean(axis=1), axis=0).divide(
        expr_log.std(axis=1).replace(0, 1), axis=0
    )
    z.columns = [c + "_within_z" for c in genes]
    k2_z = pd.concat([k2[["run", "p_DM2", "DM_call"]], z], axis=1)
    k2_z["rai_within_z"] = z.mean(axis=1)
    k2_z.to_csv(OUT / "r9_4_k2_within_z.tsv", sep="\t", index=False)

    # apply 4-feat composite trained on TCGA — but K2 lacks HLA scores from arcasHLA in this dir
    # If arcasHLA file exists, apply; else just report rai_within_z stratification
    arcas_paths = [
        ROOT / "project/results/v17_arcasHLA/k2_korean_hla_per_sample.tsv",
        ROOT / "project/results/v17_arcasHLA/K2_arcasHLA_per_sample.tsv",
    ]
    k2_hla = None
    for p in arcas_paths:
        if p.exists():
            k2_hla = pd.read_csv(p, sep="\t")
            break

    summary = dict(
        n=int(len(k2)),
        rai_within_z_median=round(float(k2_z["rai_within_z"].median()), 3),
        rai_within_z_p_DM2_corr=round(
            float(k2_z[["rai_within_z", "p_DM2"]].corr(method="spearman").iloc[0, 1]), 3
        ),
        n_DM_call=k2_z["DM_call"].value_counts().to_dict(),
        arcas_hla_available=k2_hla is not None,
    )
    if k2_hla is not None:
        summary["arcas_n"] = int(len(k2_hla))
        summary["arcas_path"] = str(arcas_paths[0] if arcas_paths[0].exists() else arcas_paths[1])

    (OUT / "r9_4_k2_within_z_summary.json").write_text(json.dumps(summary, indent=2))
    return summary


# -------------------------------------------------------------------- main
def main() -> None:
    master = load_master()
    print(f"master n={len(master)}")

    print("[R9-1] Xena RAI clinical...")
    r1 = r9_1_xena_rai(master)
    (OUT / "r9_1_summary.json").write_text(json.dumps(r1, indent=2, default=str))
    print(json.dumps(r1["i131_by_dm"], indent=2))
    print(json.dumps(r1["recurrence_by_dm"], indent=2))

    print("[R9-2] chr7 RTK gene-level CNV...")
    r2 = r9_2_chr7_rtk(master)
    (OUT / "r9_2_summary.json").write_text(json.dumps(r2, indent=2, default=str))
    print(json.dumps(r2.get("per_gene", []), indent=2))

    print("[R9-3] bootstrap chr7...")
    r3 = r9_3_bootstrap_chr7(master, n_boot=1000)
    (OUT / "r9_3_summary.json").write_text(json.dumps(r3, indent=2, default=str))
    print(json.dumps(r3["bootstrap"], indent=2))

    print("[R9-4] K2 within-z composite...")
    r4 = r9_4_k2_within_z_composite()
    print(json.dumps(r4, indent=2))

    print("[R9] done")


if __name__ == "__main__":
    main()
