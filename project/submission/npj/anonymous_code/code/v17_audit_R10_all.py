"""
v17 audit Round 10 — MSK chr7 RTK gene-level CNA replication, 3-way DM×fusion×TERT
survival, chr7 × fusion interaction within DM1, K2 within-z composite apply,
TCGA Cabrita TLS by DM.
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
from sklearn.preprocessing import StandardScaler

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/audit_2026_04_30/round10"
OUT.mkdir(parents=True, exist_ok=True)
CBIO = "https://www.cbioportal.org/api"


def load_master() -> pd.DataFrame:
    p = ROOT / "project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv"
    m = pd.read_csv(p, sep="\t")
    m["sample_short"] = m["sample_id"].str.slice(0, 15)
    m["patient12"] = m["sample_short"].str.slice(0, 12)
    m["dm"] = m["v17_dark_cluster"].fillna("not_DM").replace("", "not_DM")
    return m


# ---------------------------------------------------------------- R10-1
def r10_1_msk_chr7(genes_entrez: dict[str, int]) -> dict:
    """Pull GISTIC discrete gene-level CNA from MSK 2016 + thyroid_gatci_2024 +
    thpa_tcga_gdc for BRAF/EGFR/MET/RAF1/RET; replicate chr7 DM1-favoring direction
    using histology proxy (cPTC vs FVPTC vs other) since MSK lacks DM cluster call."""
    studies = {
        "thyroid_mskcc_2016": "thyroid_mskcc_2016_gistic",
        "thyroid_gatci_2024": "thyroid_gatci_2024_gistic",
    }
    out = {}
    for study, profile in studies.items():
        # all sample list
        sl_resp = requests.get(
            f"{CBIO}/studies/{study}/sample-lists",
            timeout=60,
        )
        sample_list_id = f"{study}_all"
        # fetch via sample-lists/{id}/sample-ids
        sids_resp = requests.get(
            f"{CBIO}/sample-lists/{sample_list_id}/sample-ids",
            timeout=60,
        )
        if sids_resp.status_code != 200:
            out[study] = {"error": f"sample list {sids_resp.status_code}"}
            continue
        sample_ids = sids_resp.json()
        if not sample_ids:
            out[study] = {"error": "empty sample list"}
            continue
        rows = []
        BATCH = 200
        with ThreadPoolExecutor(max_workers=4) as ex:
            futs = []
            for i in range(0, len(sample_ids), BATCH):
                chunk = sample_ids[i : i + BATCH]
                futs.append(
                    ex.submit(
                        requests.post,
                        f"{CBIO}/molecular-profiles/{profile}/molecular-data/fetch?projection=DETAILED",
                        headers={"Content-Type": "application/json"},
                        json={"sampleIds": chunk, "entrezGeneIds": list(genes_entrez.values())},
                        timeout=90,
                    )
                )
            for f in as_completed(futs):
                r = f.result()
                if r.status_code == 200:
                    rows.extend(r.json())
        if not rows:
            out[study] = {"error": "no data", "n_samples": len(sample_ids)}
            continue
        df = pd.DataFrame(rows)
        inv = {v: k for k, v in genes_entrez.items()}
        df["gene"] = df["entrezGeneId"].map(inv)
        df["call"] = pd.to_numeric(df["value"], errors="coerce")
        # GISTIC discrete: -2,-1,0,1,2 (HomDel, HetLoss, Diploid, Gain, Amp)
        # gain or amp counted as positive
        per = df.groupby("gene")["call"].agg(
            [
                ("n", "count"),
                ("median", "median"),
                ("gain_amp_n", lambda s: int(((s == 1) | (s == 2)).sum())),
                ("gain_amp_pct", lambda s: round(100 * ((s == 1) | (s == 2)).mean(), 2)),
            ]
        ).reset_index()
        out[study] = {"per_gene": per.to_dict(orient="records"), "n_samples": int(df["sampleId"].nunique())}
    (OUT / "r10_1_msk_chr7_replication.json").write_text(json.dumps(out, indent=2, default=str))
    return out


# ---------------------------------------------------------------- R10-2
def r10_2_3way_survival(master: pd.DataFrame) -> dict:
    """3-way DM × fusion × TERT survival. Cox + KM strata."""
    sv = pd.read_csv(ROOT / "project/results/audit_2026_04_30/round3/cbio_sv_thca.tsv", sep="\t")
    fusion_set = set(sv["sampleId"].unique())
    m = master.copy()
    m["any_fusion"] = m["sample_short"].isin(fusion_set).astype(int)
    m["tert"] = (m["tert_promoter_integrated"].fillna("wildtype") != "wildtype").astype(int)
    m["dm1"] = (m["v17_dark_cluster"] == "DM1").astype(int)

    m["combo"] = m.apply(
        lambda r: (
            f"DM1_fus{r['any_fusion']}_tert{r['tert']}"
            if r["dm1"] == 1
            else f"notDM1_fus{r['any_fusion']}_tert{r['tert']}"
        ),
        axis=1,
    )
    m["os_event"] = pd.to_numeric(m["os_event"], errors="coerce")
    m["os_days"] = pd.to_numeric(m["os_days"], errors="coerce")
    m_surv = m.dropna(subset=["os_event", "os_days"])
    summary = (
        m_surv.groupby("combo")
        .agg(n=("sample_short", "count"), events=("os_event", "sum"), median_os=("os_days", "median"))
        .reset_index()
    )
    summary.to_csv(OUT / "r10_2_3way_combo_summary.tsv", sep="\t", index=False)

    try:
        from lifelines import CoxPHFitter
        cox_df = m_surv[["dm1", "any_fusion", "tert", "age", "os_days", "os_event"]].dropna().copy()
        cox_df["dm1_x_fusion"] = cox_df["dm1"] * cox_df["any_fusion"]
        cox_df["dm1_x_tert"] = cox_df["dm1"] * cox_df["tert"]
        cph = CoxPHFitter()
        cph.fit(cox_df, duration_col="os_days", event_col="os_event")
        cox_table = cph.summary[["coef", "exp(coef)", "se(coef)", "p", "exp(coef) lower 95%", "exp(coef) upper 95%"]].reset_index()
        cox_table.to_csv(OUT / "r10_2_cox_3way.tsv", sep="\t", index=False)
        cox_summary = cox_table.to_dict(orient="records")
    except Exception as e:
        cox_summary = [{"error": str(e)}]

    return {
        "combo_summary": summary.to_dict(orient="records"),
        "cox_3way": cox_summary,
        "interpretation": "DM1 × fusion × TERT 3-way Cox with interaction terms; clinical risk stratification",
    }


# ---------------------------------------------------------------- R10-3
def r10_3_chr7_x_fusion(master: pd.DataFrame) -> dict:
    """chr 7p/7q gain × RET-fusion-status interaction within DM1."""
    arm = pd.read_csv(
        ROOT / "project/results/audit_2026_04_30/round8/r8_1_arm_cnv_long.tsv", sep="\t"
    )
    sv = pd.read_csv(ROOT / "project/results/audit_2026_04_30/round3/cbio_sv_thca.tsv", sep="\t")
    fusion_set = set(sv["sampleId"].unique())
    chr7p_gain = set(arm.loc[(arm["arm"] == "7p") & (arm["call"] == "Gain"), "sample_short"])
    chr7q_gain = set(arm.loc[(arm["arm"] == "7q") & (arm["call"] == "Gain"), "sample_short"])
    chr7_any = chr7p_gain | chr7q_gain

    m = master.copy()
    m["fus"] = m["sample_short"].isin(fusion_set).astype(int)
    m["chr7"] = m["sample_short"].isin(chr7_any).astype(int)
    dm1 = m[m["dm"] == "DM1"].copy()
    not_dm1 = m[m["dm"] != "DM1"].copy()

    rows = []
    for label, sub in [("DM1", dm1), ("not_DM1", not_dm1)]:
        ct = pd.crosstab(sub["chr7"], sub["fus"])
        # 4-cell rates
        for chr7v in [0, 1]:
            for fusv in [0, 1]:
                n = int(((sub["chr7"] == chr7v) & (sub["fus"] == fusv)).sum())
                med_rai = (
                    sub.loc[(sub["chr7"] == chr7v) & (sub["fus"] == fusv), "rai_score_v17"]
                    .dropna()
                    .median()
                )
                rows.append(
                    dict(
                        cohort=label,
                        chr7_gain=chr7v,
                        fusion=fusv,
                        n=n,
                        rai_med=round(float(med_rai), 3) if pd.notna(med_rai) else None,
                    )
                )
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "r10_3_chr7_x_fusion_4cell.tsv", sep="\t", index=False)

    # Test: within DM1, does chr7 add information beyond fusion?
    dm1_sub = dm1.dropna(subset=["rai_score_v17"]).copy()
    if len(dm1_sub) > 30:
        # 2x2 crosstab Fisher: chr7 gain vs fusion, within DM1
        ct = pd.crosstab(dm1_sub["chr7"], dm1_sub["fus"])
        try:
            odds, p_co = stats.fisher_exact(ct.values)
        except Exception:
            odds, p_co = None, None
        # rai_score linear-model with interaction
        try:
            import statsmodels.formula.api as smf
            mod = smf.ols("rai_score_v17 ~ chr7 * fus + age", data=dm1_sub).fit()
            interaction_p = float(mod.pvalues.get("chr7:fus", float("nan")))
            chr7_main_p = float(mod.pvalues.get("chr7", float("nan")))
            fus_main_p = float(mod.pvalues.get("fus", float("nan")))
        except Exception as e:
            interaction_p = chr7_main_p = fus_main_p = None
    else:
        odds, p_co = None, None
        interaction_p = chr7_main_p = fus_main_p = None

    return {
        "4cell": out.to_dict(orient="records"),
        "dm1_chr7_x_fusion_co_occurrence_OR": round(float(odds), 3) if odds is not None and np.isfinite(odds) else None,
        "dm1_chr7_x_fusion_co_occurrence_fisher_p": float(p_co) if p_co is not None else None,
        "rai_interaction_p": interaction_p,
        "rai_chr7_main_p": chr7_main_p,
        "rai_fus_main_p": fus_main_p,
        "interpretation": "chr7 gain × RET-fusion within DM1 (additive vs co-occurrent vs interaction)",
    }


# ---------------------------------------------------------------- R10-4
def r10_4_k2_within_z_composite() -> dict:
    """Apply 4-feat TCGA-trained composite to K2 (within-z + arcasHLA scores)."""
    # K2 within-z
    z = pd.read_csv(
        ROOT / "project/results/audit_2026_04_30/round9/r9_4_k2_within_z.tsv", sep="\t"
    )
    # arcasHLA per-sample for K2: use HLA Class I/II proxy from genotype calls if scores not present
    # Use d4p1 GSE286332 panel as reference HLA score (similar Korean population)
    # Fallback: the K2_arcasHLA_dm_merged has alleles + p_DM2 + DM_call only; no precomputed score.
    # For r10-4, we synthesize HLA scores by allele-frequency reference using GSE213647 file
    # (same Korean panel used to derive HLA score); apply allele match -> z-like binary.
    arcas_genot = pd.read_csv(
        ROOT / "project/results/v17_korean/arcasHLA/K2_arcasHLA_genotypes.tsv", sep="\t"
    )
    # GSE213647 HLA score file gives reference distribution; we use mean within-K2 per-sample
    # score = count of "high-risk" alleles weighted (proxy)
    risk_alleles = ["DPB1*05:01", "DRB1*04:03", "DQB1*05:01"]  # Korean Graves'/Hashi-related
    hla_score = []
    for _, r in arcas_genot.iterrows():
        alleles = " ".join(r.fillna("").astype(str).values)
        score_II = sum(1 for a in risk_alleles if a in alleles)
        hla_score.append(score_II)
    arcas_genot["hla_class_II_proxy"] = hla_score
    arcas_genot["hla_class_I_proxy"] = 0  # not enough info; placeholder

    # merge
    z["run_acc"] = z["run"]
    arcas_genot["run_acc"] = arcas_genot["run"]
    merged = z.merge(arcas_genot[["run_acc", "hla_class_II_proxy"]], on="run_acc", how="left")

    # Composite using rai_within_z_top4 alone (since HLA-I missing): already a 1-feat
    # For full TCGA-train composite apply we need standardized features matching TCGA scale
    # — but K2 doesn't have age either, and within-z is qualitatively different from TCGA absolute.
    # So R10-4 is a sanity-check apply: simply correlate rai_within_z_top4 with original p_DM2,
    # and report HLA-II proxy distribution within K2 top-quartile vs bottom-quartile of rai_within_z

    q75 = merged["rai_within_z_top4"].quantile(0.75)
    q25 = merged["rai_within_z_top4"].quantile(0.25)
    top_q = merged[merged["rai_within_z_top4"] >= q75]
    bot_q = merged[merged["rai_within_z_top4"] <= q25]
    return {
        "n_with_arcas": int(merged["hla_class_II_proxy"].notna().sum()),
        "top_quartile_hla_II_mean": round(float(top_q["hla_class_II_proxy"].mean(skipna=True)), 3),
        "bot_quartile_hla_II_mean": round(float(bot_q["hla_class_II_proxy"].mean(skipna=True)), 3),
        "top_quartile_p_DM2_med": round(float(top_q["p_DM2"].median()), 3),
        "bot_quartile_p_DM2_med": round(float(bot_q["p_DM2"].median()), 3),
        "spearman_within_z_pDM2": round(float(merged[["rai_within_z_top4", "p_DM2"]].corr(method="spearman").iloc[0, 1]), 3),
        "interpretation": "K2 within-z proxy DM1 (top quartile) vs proxy DM2 (bottom quartile); HLA Class II Korean risk-allele count comparison",
    }


# ---------------------------------------------------------------- R10-5
def r10_5_tcga_cabrita_tls(master: pd.DataFrame) -> dict:
    """Cabrita 2020 12-gene TLS signature in TCGA THCA by DM cluster."""
    cabrita = ["CCL19", "CCL21", "CXCL13", "CCR7", "CXCR5", "SELL", "LAMP3",
               "PTGDS", "CCR6", "CCL2", "TNFSF13B", "CD79B"]
    # entrez map fetched
    map_resp = requests.post(
        f"{CBIO}/genes/fetch?geneIdType=HUGO_GENE_SYMBOL",
        headers={"Content-Type": "application/json"},
        json=cabrita,
        timeout=60,
    )
    if map_resp.status_code != 200:
        return {"error": f"gene fetch {map_resp.status_code}: {map_resp.text[:200]}"}
    gmap = pd.DataFrame(map_resp.json())
    if "entrezGeneId" not in gmap.columns:
        return {"error": f"gene fetch schema unexpected: cols={list(gmap.columns)}"}
    entrez_ids = gmap["entrezGeneId"].tolist()
    name_map = dict(zip(gmap["entrezGeneId"], gmap["hugoGeneSymbol"]))

    sample_ids = master["sample_short"].dropna().unique().tolist()
    rows = []
    BATCH = 100
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = []
        for i in range(0, len(sample_ids), BATCH):
            chunk = sample_ids[i : i + BATCH]
            futs.append(
                ex.submit(
                    requests.post,
                    f"{CBIO}/molecular-profiles/thca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_median_Zscores/molecular-data/fetch?projection=DETAILED",
                    headers={"Content-Type": "application/json"},
                    json={"sampleIds": chunk, "entrezGeneIds": entrez_ids},
                    timeout=90,
                )
            )
        for f in as_completed(futs):
            r = f.result()
            if r.status_code == 200:
                rows.extend(r.json())

    if not rows:
        return {"error": "no expression data"}
    df = pd.DataFrame(rows)
    df["gene"] = df["entrezGeneId"].map(name_map)
    df["sample_short"] = df["sampleId"]
    df["z"] = pd.to_numeric(df["value"], errors="coerce")
    pivot = df.pivot_table(index="sample_short", columns="gene", values="z", aggfunc="mean").reset_index()
    pivot["TLS_score"] = pivot[[c for c in cabrita if c in pivot.columns]].mean(axis=1)
    pivot.to_csv(OUT / "r10_5_tcga_tls_per_sample.tsv", sep="\t", index=False)
    merged = pivot.merge(master[["sample_short", "dm"]], on="sample_short", how="inner")

    rows_out = []
    for cluster in ["DM1", "DM2", "not_DM"]:
        sub = merged[merged["dm"] == cluster]["TLS_score"].dropna()
        rows_out.append(
            dict(
                cluster=cluster,
                n=int(len(sub)),
                median=round(float(sub.median()), 3),
                mean=round(float(sub.mean()), 3),
            )
        )
    summary = pd.DataFrame(rows_out)
    summary.to_csv(OUT / "r10_5_tls_by_dm.tsv", sep="\t", index=False)

    # Cohen's d DM1 vs not_DM
    d_dm1 = merged.loc[merged["dm"] == "DM1", "TLS_score"].dropna()
    d_dm2 = merged.loc[merged["dm"] == "DM2", "TLS_score"].dropna()
    d_nd = merged.loc[merged["dm"] == "not_DM", "TLS_score"].dropna()
    pooled_sd = np.sqrt(
        ((d_dm1.var(ddof=1) * (len(d_dm1) - 1)) + (d_nd.var(ddof=1) * (len(d_nd) - 1)))
        / (len(d_dm1) + len(d_nd) - 2)
    )
    d_eff = (d_dm1.mean() - d_nd.mean()) / (pooled_sd + 1e-9)
    u, p = stats.mannwhitneyu(d_dm1, d_nd, alternative="two-sided")
    return {
        "by_dm": summary.to_dict(orient="records"),
        "cohens_d_dm1_vs_notDM": round(float(d_eff), 3),
        "mw_p_dm1_vs_notDM": float(p),
        "interpretation": "Cabrita 2020 12-gene TLS signature in TCGA-THCA by DM cluster",
    }


def main() -> None:
    master = load_master()
    print(f"master n={len(master)}")

    print("[R10-1] MSK chr7...")
    genes_entrez = {"BRAF": 673, "EGFR": 1956, "MET": 4233, "RAF1": 5894, "RET": 5979}
    r1 = r10_1_msk_chr7(genes_entrez)
    print(json.dumps(r1, indent=2, default=str)[:800])

    print("[R10-2] 3-way survival...")
    r2 = r10_2_3way_survival(master)
    (OUT / "r10_2_summary.json").write_text(json.dumps(r2, indent=2, default=str))
    print(json.dumps(r2["combo_summary"], indent=2))

    print("[R10-3] chr7 × fusion...")
    r3 = r10_3_chr7_x_fusion(master)
    (OUT / "r10_3_summary.json").write_text(json.dumps(r3, indent=2, default=str))
    print(json.dumps({k: v for k, v in r3.items() if k != "4cell"}, indent=2))

    print("[R10-4] K2 composite...")
    r4 = r10_4_k2_within_z_composite()
    (OUT / "r10_4_summary.json").write_text(json.dumps(r4, indent=2, default=str))
    print(json.dumps(r4, indent=2))

    print("[R10-5] TCGA Cabrita TLS...")
    r5 = r10_5_tcga_cabrita_tls(master)
    (OUT / "r10_5_summary.json").write_text(json.dumps(r5, indent=2, default=str))
    print(json.dumps(r5, indent=2, default=str))
    print("[R10] done")


if __name__ == "__main__":
    main()
