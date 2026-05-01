"""R10-5 standalone — TCGA Cabrita TLS by DM."""
from __future__ import annotations
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy import stats

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/audit_2026_04_30/round10"
OUT.mkdir(parents=True, exist_ok=True)
CBIO = "https://www.cbioportal.org/api"


def main() -> None:
    master = pd.read_csv(
        ROOT / "project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv", sep="\t"
    )
    master["sample_short"] = master["sample_id"].str.slice(0, 15)
    master["dm"] = master["v17_dark_cluster"].fillna("not_DM").replace("", "not_DM")

    cabrita = ["CCL19", "CCL21", "CXCL13", "CCR7", "CXCR5", "SELL", "LAMP3",
               "PTGDS", "CCR6", "CCL2", "TNFSF13B", "CD79B"]
    map_resp = requests.post(
        f"{CBIO}/genes/fetch?geneIdType=HUGO_GENE_SYMBOL",
        headers={"Content-Type": "application/json"},
        json=cabrita,
        timeout=60,
    )
    gmap = pd.DataFrame(map_resp.json())
    entrez_ids = gmap["entrezGeneId"].tolist()
    name_map = dict(zip(gmap["entrezGeneId"], gmap["hugoGeneSymbol"]))
    print("genes mapped:", len(entrez_ids))

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

    df = pd.DataFrame(rows)
    print("rows:", len(df))
    df["gene"] = df["entrezGeneId"].map(name_map)
    df["sample_short"] = df["sampleId"]
    df["z"] = pd.to_numeric(df["value"], errors="coerce")
    pivot = df.pivot_table(index="sample_short", columns="gene", values="z", aggfunc="mean").reset_index()
    cabrita_in = [c for c in cabrita if c in pivot.columns]
    pivot["TLS_score"] = pivot[cabrita_in].mean(axis=1)
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

    d_dm1 = merged.loc[merged["dm"] == "DM1", "TLS_score"].dropna()
    d_nd = merged.loc[merged["dm"] == "not_DM", "TLS_score"].dropna()
    d_dm2 = merged.loc[merged["dm"] == "DM2", "TLS_score"].dropna()
    pooled_sd = np.sqrt(
        ((d_dm1.var(ddof=1) * (len(d_dm1) - 1)) + (d_nd.var(ddof=1) * (len(d_nd) - 1)))
        / (len(d_dm1) + len(d_nd) - 2)
    )
    d_eff_dm1_nd = (d_dm1.mean() - d_nd.mean()) / (pooled_sd + 1e-9)
    u, p_dm1_nd = stats.mannwhitneyu(d_dm1, d_nd, alternative="two-sided")

    pooled_sd2 = np.sqrt(
        ((d_dm1.var(ddof=1) * (len(d_dm1) - 1)) + (d_dm2.var(ddof=1) * (len(d_dm2) - 1)))
        / (len(d_dm1) + len(d_dm2) - 2)
    )
    d_eff_dm1_dm2 = (d_dm1.mean() - d_dm2.mean()) / (pooled_sd2 + 1e-9)
    u, p_dm1_dm2 = stats.mannwhitneyu(d_dm1, d_dm2, alternative="two-sided")

    summary_d = dict(
        genes_used=cabrita_in,
        n_genes_recovered=len(cabrita_in),
        by_dm=summary.to_dict(orient="records"),
        cohens_d_DM1_vs_notDM=round(float(d_eff_dm1_nd), 3),
        mw_p_DM1_vs_notDM=float(p_dm1_nd),
        cohens_d_DM1_vs_DM2=round(float(d_eff_dm1_dm2), 3),
        mw_p_DM1_vs_DM2=float(p_dm1_dm2),
        interpretation="Cabrita 2020 12-gene TLS in TCGA-THCA by DM cluster",
    )
    (OUT / "r10_5_summary.json").write_text(json.dumps(summary_d, indent=2))
    print(json.dumps(summary_d, indent=2))


if __name__ == "__main__":
    main()
