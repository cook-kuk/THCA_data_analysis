"""
TCGA-PAAD fetcher (cBioPortal REST, no auth).

Pulls:
  - clinical (patient + sample)
  - mutation calls (MAF) for the standard PDAC drivers + extended panel
  - mRNA z-scores for Moffitt classical/basal genes + immune modules
  - sample list

Outputs to /data/pdac_poc/raw/*.tsv and processed JSON manifest.
"""

import json
import time
from pathlib import Path
import pandas as pd
import requests

API = "https://www.cbioportal.org/api"
STUDY = "paad_tcga_pan_can_atlas_2018"
OUT = Path("/data/pdac_poc/raw")
OUT.mkdir(parents=True, exist_ok=True)

S = requests.Session()
S.headers.update({"Accept": "application/json"})


def get(url, **params):
    r = S.get(url, params=params, timeout=60)
    r.raise_for_status()
    return r.json()


def post(url, body, **params):
    r = S.post(url, json=body, params=params, timeout=120,
               headers={"Content-Type": "application/json"})
    r.raise_for_status()
    return r.json()


def fetch_clinical():
    samples = get(f"{API}/studies/{STUDY}/samples")
    sample_df = pd.DataFrame(samples)[["sampleId", "patientId", "sampleType"]]
    pat_clin = get(f"{API}/studies/{STUDY}/clinical-data", clinicalDataType="PATIENT")
    samp_clin = get(f"{API}/studies/{STUDY}/clinical-data", clinicalDataType="SAMPLE")
    pdf = (pd.DataFrame(pat_clin)
           .pivot_table(index="patientId", columns="clinicalAttributeId",
                        values="value", aggfunc="first").reset_index())
    sdf = (pd.DataFrame(samp_clin)
           .pivot_table(index="sampleId", columns="clinicalAttributeId",
                        values="value", aggfunc="first").reset_index())
    df = sample_df.merge(sdf, on="sampleId", how="left").merge(pdf, on="patientId", how="left")
    df.to_csv(OUT / "clinical.tsv", sep="\t", index=False)
    return df


def fetch_mutations(gene_symbols):
    genes = post(f"{API}/genes/fetch", body=gene_symbols, geneIdType="HUGO_GENE_SYMBOL")
    entrez = [g["entrezGeneId"] for g in genes]
    profile_id = f"{STUDY}_mutations"
    sample_list_id = f"{STUDY}_sequenced"
    body = {"entrezGeneIds": entrez, "sampleListId": sample_list_id}
    muts = post(f"{API}/molecular-profiles/{profile_id}/mutations/fetch",
                body=body, projection="DETAILED")
    if not muts:
        return pd.DataFrame()
    mdf = pd.DataFrame(muts)
    keep = [c for c in ["sampleId", "patientId", "entrezGeneId", "gene",
                        "proteinChange", "mutationType", "variantType",
                        "startPosition", "endPosition", "referenceAllele",
                        "variantAllele", "tumorAltCount", "tumorRefCount"]
            if c in mdf.columns]
    if "gene" in mdf.columns:
        mdf["hugo"] = mdf["gene"].apply(lambda g: g.get("hugoGeneSymbol") if isinstance(g, dict) else None)
    mdf[keep + (["hugo"] if "hugo" in mdf.columns else [])].to_csv(
        OUT / "mutations.tsv", sep="\t", index=False)
    return mdf


def fetch_mrna_zscores(gene_symbols, profile="paad_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_median_Zscores"):
    genes = post(f"{API}/genes/fetch", body=gene_symbols, geneIdType="HUGO_GENE_SYMBOL")
    entrez_to_sym = {g["entrezGeneId"]: g["hugoGeneSymbol"] for g in genes}
    body = {"entrezGeneIds": list(entrez_to_sym.keys()),
            "sampleListId": f"{STUDY}_rna_seq_v2_mrna"}
    data = post(f"{API}/molecular-profiles/{profile}/molecular-data/fetch", body=body)
    if not data:
        return pd.DataFrame()
    df = pd.DataFrame(data)
    df["gene"] = df["entrezGeneId"].map(entrez_to_sym)
    wide = df.pivot_table(index="sampleId", columns="gene", values="value", aggfunc="first")
    wide.to_csv(OUT / f"mrna_z_{profile.split('_')[-1]}.tsv", sep="\t")
    return wide


# --- gene panels --------------------------------------------------------
PDAC_DRIVERS = ["KRAS", "TP53", "CDKN2A", "SMAD4", "BRCA2", "BRCA1", "PALB2",
                "ATM", "MLH1", "MSH2", "MSH6", "PMS2", "STK11", "ARID1A",
                "RNF43", "GNAS", "PTEN", "PIK3CA", "TGFBR2", "ACVR1B",
                "MAP2K4", "RBM10", "EGFR", "ERBB2", "ERBB3", "BRAF", "NTRK1",
                "NTRK2", "NTRK3", "ROS1", "ALK", "FGFR2", "CDK4", "CDK6"]

MOFFITT_BASAL = ["VGLL1", "UCA1", "S100A2", "LY6D", "SPRR3", "SPRR1B", "LEMD1",
                 "LYPD3", "KRT15", "CTSL2", "DHRS9", "AREG", "CST6", "SERPINB3",
                 "KRT6A", "SERPINB4", "FAM83A", "SCEL", "FGFBP1", "KRT7",
                 "KRT17", "GPR87", "TNS4", "SLC2A1", "ANXA8L2"]

MOFFITT_CLASSICAL = ["BTNL8", "FAM3D", "ATAD4", "AGR3", "CTSE", "LOC400573",
                     "LYZ", "TFF2", "TFF1", "ANXA10", "LGALS4", "ECT2",
                     "CLRN3", "MYO1A", "CLDN18", "LRRC31", "TFF3", "BTNL3",
                     "CDX2", "SERPINA10", "VSIG2", "TSPAN8", "ST6GALNAC1",
                     "AGR2", "TOX3"]

# Bailey 2016 Nature subtypes use additional immune/squamous panels
BAILEY_IMMUNE = ["CD8A", "CD8B", "CXCL9", "CXCL10", "CXCL11", "GZMA", "GZMB",
                 "PRF1", "IFNG", "STAT1", "IRF1", "TBX21", "EOMES", "CD274",
                 "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "TOX",
                 "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "NLRC5",
                 "PSMB8", "PSMB9", "HLA-DRA", "HLA-DRB1", "HLA-DPA1",
                 "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "CIITA"]

PDAC_TME_MYELOID = ["CD68", "CD163", "MRC1", "CSF1R", "MARCO", "MSR1", "CD86",
                    "CD80", "ITGAM", "ITGAX", "ARG1", "IDO1", "TGFB1", "IL10",
                    "CCL2", "CCL5", "CXCL12", "VEGFA"]

PDAC_TME_CAF = ["ACTA2", "FAP", "POSTN", "COL1A1", "COL1A2", "COL3A1", "PDPN",
                "PDGFRA", "PDGFRB", "DES", "S100A4", "VIM", "DCN", "LUM",
                "MMP2", "MMP9", "TIMP1", "IL6", "CXCL12", "HAS2", "LRRC15"]

PDAC_TLS = ["CXCL13", "CCL19", "CCL21", "CCR7", "LAMP3", "CD79A", "MS4A1",
            "BCL6", "CXCR5", "AICDA", "JCHAIN"]


def main():
    print(f"[01] cBioPortal study={STUDY}")
    clin = fetch_clinical()
    print(f"  clinical rows={len(clin)} samples; cols={len(clin.columns)}")
    mut = fetch_mutations(PDAC_DRIVERS)
    print(f"  mutations rows={len(mut)} (drivers panel)")
    panel = sorted(set(MOFFITT_BASAL + MOFFITT_CLASSICAL + BAILEY_IMMUNE +
                       PDAC_TME_MYELOID + PDAC_TME_CAF + PDAC_TLS + PDAC_DRIVERS))
    z = fetch_mrna_zscores(panel)
    print(f"  mrna_z rows={len(z)} cols={z.shape[1] if len(z) else 0}")
    manifest = {
        "study": STUDY,
        "samples": int(len(clin)),
        "mutations": int(len(mut)),
        "mrna_samples": int(len(z)),
        "mrna_genes": int(z.shape[1]) if len(z) else 0,
        "panels": {
            "moffitt_basal": MOFFITT_BASAL,
            "moffitt_classical": MOFFITT_CLASSICAL,
            "bailey_immune": BAILEY_IMMUNE,
            "tme_myeloid": PDAC_TME_MYELOID,
            "tme_caf": PDAC_TME_CAF,
            "tls": PDAC_TLS,
            "drivers": PDAC_DRIVERS,
        },
        "fetched_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    Path("/data/pdac_poc/processed").mkdir(parents=True, exist_ok=True)
    with open("/data/pdac_poc/processed/manifest_01_fetch.json", "w") as f:
        json.dump(manifest, f, indent=2)
    print("[01] done")


if __name__ == "__main__":
    main()
