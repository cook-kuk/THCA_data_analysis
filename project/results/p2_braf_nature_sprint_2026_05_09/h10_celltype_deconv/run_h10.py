"""
Paper 1+2 BRAF Nature sprint H10 — cell-type deconvolution within BRAF-cPTC.

Question: Is the HT-13 panel's DM1 vs DM2 separation (AUC=0.93 within BRAF-cPTC)
explained by B cells alone, or by a coordinated immune compartment (T cells, DCs,
plasma cells, etc.)?

Approach (3 layers, all on the same BRAF-cPTC DM1/DM2 stratum, n=110):
 1. nu-SVR LM-class deconvolution (8 cell types, Lu 2023 reference) reused from
    `p_deconv_2026_05_08/fractions_nu_SVR.tsv` — ground truth at 8-class resolution.
 2. Direct gene-set scores at LM22-grade resolution (canonical markers per
    cell subtype). z-score is already centred per-gene across TCGA, so a
    simple mean(z) per gene-set per sample = relative module score.
 3. TLS Cabrita 2020 12-gene panel (CCL19, CCL21, CXCL13, CCR7, CXCR5, SELL,
    LAMP3, CD79B, CR2, FCRL2, MS4A1, JCHAIN/IGJ).

For each cell-type / module:
   Cohen's d (DM1 - DM2 sign convention), Mann-Whitney p, BH FDR.

Stratum
-------
master_TSV: tissue_type=='Primary Tumor', dm in {DM1, DM2},
            molecular_subtype=='BRAF_like', histology_subtype=='cPTC'
            -> n=110 (DM1=85, DM2=25)

Outputs
-------
- h10_celltype_d.tsv : rows = (layer, cell_type, cohort_cut), cols =
                       d, p, FDR, n_DM1, n_DM2, mean_DM1, mean_DM2.
- h10_per_sample_scores.tsv : rows = sample × layer-1+2+3 module scores.
- H10_REPORT.md : <400-word summary, lead with strongest finding.
"""
from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

OUT = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "p2_braf_nature_sprint_2026_05_09/h10_celltype_deconv"
)
OUT.mkdir(parents=True, exist_ok=True)

MASTER = "/data/thca/repo_results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv"
EXPR_Z = "/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_zscore.tsv"
NUSVR = (
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "p_deconv_2026_05_08/fractions_nu_SVR.tsv"
)
CLAM_MAN = (
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_clam/slide_manifest.tsv"
)


# ---------------------------------------------------------------------- Marker sets
# LM22-grade module markers (canonical literature)
GENESETS_LM22 = {
    "Plasma_cells":  ["MZB1", "JCHAIN", "IGJ", "IGHA1", "IGHG1", "XBP1", "CD38"],
    "Memory_B":      ["MS4A1", "CD27", "CD19", "BANK1"],
    "Naive_B":       ["TCL1A", "IL4R", "IGHD", "FCER2"],
    "CD8_T":         ["CD8A", "CD8B", "GZMA", "GZMB", "GZMK", "PRF1", "NKG7"],
    "CD4_Tfh":       ["CXCL13", "CXCR5", "BCL6", "ICOS", "PDCD1"],
    "Treg":          ["FOXP3", "IL2RA", "CTLA4", "TIGIT"],
    "DC":            ["ITGAX", "CLEC9A", "IRF8", "BATF3"],
    "Macrophage_M1": ["CXCL10", "CXCL11", "IDO1", "CD86"],
    "Macrophage_M2": ["MRC1", "CD163", "MERTK", "MARCO"],
}

# Cabrita 2020 12-gene TLS signature
TLS_CABRITA = ["CCL19", "CCL21", "CXCL13", "CCR7", "CXCR5", "SELL",
               "LAMP3", "CD79B", "CR2", "FCRL2", "MS4A1", "JCHAIN"]


# ---------------------------------------------------------------------- helpers
def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    sa, sb = a.var(ddof=1), b.var(ddof=1)
    pooled = np.sqrt(((len(a) - 1) * sa + (len(b) - 1) * sb) / (len(a) + len(b) - 2))
    if pooled == 0:
        return float("nan")
    return float((a.mean() - b.mean()) / pooled)


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    finite = ~np.isnan(p)
    out = np.full_like(p, np.nan)
    if finite.sum() == 0:
        return out
    pf = p[finite]
    n = len(pf)
    order = np.argsort(pf)
    ranks = np.empty(n)
    ranks[order] = np.arange(1, n + 1)
    q = pf * n / ranks
    # enforce monotonicity in sorted order
    q_sorted = q[order]
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q[order] = np.clip(q_sorted, 0, 1)
    out[finite] = q
    return out


def compare_dm(values: dict, dm1_idx: list, dm2_idx: list, layer: str, cohort_cut: str) -> list:
    """values: dict cell_type -> Series indexed by sample_id.  Returns rows."""
    rows = []
    pvals = []
    for ct, ser in values.items():
        a = ser.reindex(dm1_idx).dropna()
        b = ser.reindex(dm2_idx).dropna()
        d = cohens_d(a.values, b.values)
        try:
            _, p = stats.mannwhitneyu(a.values, b.values, alternative="two-sided")
        except ValueError:
            p = float("nan")
        rows.append({
            "layer": layer,
            "cohort_cut": cohort_cut,
            "cell_type": ct,
            "n_DM1": int(len(a)),
            "n_DM2": int(len(b)),
            "mean_DM1": float(a.mean()) if len(a) else float("nan"),
            "mean_DM2": float(b.mean()) if len(b) else float("nan"),
            "cohen_d": d,
            "MW_p": float(p),
        })
        pvals.append(p)
    fdr = bh_fdr(np.array(pvals))
    for r, q in zip(rows, fdr):
        r["BH_FDR"] = float(q) if not np.isnan(q) else float("nan")
    return rows


# ---------------------------------------------------------------------- load
print("[load] master")
master = pd.read_csv(MASTER, sep="\t", low_memory=False)
master = master[master["tissue_type"] == "Primary Tumor"].copy()
master = master[master["dm"].isin(["DM1", "DM2"])].copy()

# Master cohort cut: BRAF_like × cPTC
braf_cptc = master[
    (master["molecular_subtype"] == "BRAF_like")
    & (master["histology_subtype"] == "cPTC")
].copy()
print(
    f"  BRAF_like+cPTC DM1/DM2 master n={len(braf_cptc)} "
    f"(DM1={(braf_cptc['dm']=='DM1').sum()}, DM2={(braf_cptc['dm']=='DM2').sum()})"
)

# CLAM cohort cut (intersection with CLAM image manifest)
clam_man = pd.read_csv(CLAM_MAN, sep="\t")
# CLAM manifest: submitter_id is TCGA-XX-XXXX form; map to master sample_id
clam_pids = set(clam_man["submitter_id"].astype(str).tolist())
# master patient12 column may exist; otherwise derive from sample_id
def _pid12(sid):
    return "-".join(sid.split("-")[:3])
braf_cptc["_pid12"] = braf_cptc["sample_id"].astype(str).map(_pid12)
braf_cptc_clam = braf_cptc[braf_cptc["_pid12"].isin(clam_pids)].copy()
braf_cptc = braf_cptc.drop(columns=["_pid12"])
print(
    f"  BRAF_like+cPTC ∩ CLAM image cohort n={len(braf_cptc_clam)} "
    f"(DM1={(braf_cptc_clam['dm']=='DM1').sum()}, "
    f"DM2={(braf_cptc_clam['dm']=='DM2').sum()})"
)

# ---------------------------------------------------------------------- expression
print("[load] z-scored expression")
expr = pd.read_csv(EXPR_Z, sep="\t", index_col=0)
print(f"  expression shape: {expr.shape}")  # genes x samples


def gene_mean_z(gset: list, expr_df: pd.DataFrame) -> pd.Series:
    """Return per-sample mean z-score across genes that exist in expr_df."""
    present = [g for g in gset if g in expr_df.index]
    if not present:
        return pd.Series(dtype=float)
    return expr_df.loc[present].mean(axis=0)


# ---------------------------------------------------------------------- nu-SVR
print("[load] nu-SVR fractions (Lu 2023 ref, 8 cell types)")
nusvr = pd.read_csv(NUSVR, sep="\t", index_col=0)
print(f"  nu-SVR shape: {nusvr.shape}")


# ---------------------------------------------------------------------- compute
def run_cohort(label_df: pd.DataFrame, cohort_cut: str, all_rows: list,
               per_sample_records: dict):
    dm1_ids = label_df[label_df["dm"] == "DM1"]["sample_id"].tolist()
    dm2_ids = label_df[label_df["dm"] == "DM2"]["sample_id"].tolist()

    # Layer 1 — nu-SVR fractions (8 types)
    nu_dict = {}
    for ct in nusvr.columns:
        ser = nusvr[ct]
        ser = ser.reindex(label_df["sample_id"])
        nu_dict[ct] = ser
        for sid, val in ser.items():
            per_sample_records.setdefault(sid, {})[f"nuSVR_{ct}"] = val
    all_rows.extend(compare_dm(nu_dict, dm1_ids, dm2_ids,
                               layer="nuSVR_8class", cohort_cut=cohort_cut))

    # Layer 2 — gene-set scores at LM22 resolution
    gs_dict = {}
    for ct, gset in GENESETS_LM22.items():
        ser = gene_mean_z(gset, expr)
        ser = ser.reindex(label_df["sample_id"])
        gs_dict[ct] = ser
        for sid, val in ser.items():
            per_sample_records.setdefault(sid, {})[f"GS_{ct}"] = val
    all_rows.extend(compare_dm(gs_dict, dm1_ids, dm2_ids,
                               layer="geneset_LM22", cohort_cut=cohort_cut))

    # Layer 3 — TLS Cabrita 2020
    tls_present = [g for g in TLS_CABRITA if g in expr.index]
    print(f"  [{cohort_cut}] TLS Cabrita present: {len(tls_present)}/12 -> {tls_present}")
    tls_ser = gene_mean_z(TLS_CABRITA, expr).reindex(label_df["sample_id"])
    for sid, val in tls_ser.items():
        per_sample_records.setdefault(sid, {})["TLS_Cabrita12"] = val
    all_rows.extend(compare_dm({"TLS_Cabrita12": tls_ser}, dm1_ids, dm2_ids,
                               layer="TLS_signature", cohort_cut=cohort_cut))


all_rows: list = []
per_sample: dict = {}

run_cohort(braf_cptc, "BRAF_cPTC_n110", all_rows, per_sample)
run_cohort(braf_cptc_clam, f"BRAF_cPTC_CLAM_n{len(braf_cptc_clam)}",
           all_rows, per_sample)

df_out = pd.DataFrame(all_rows)
df_out = df_out.sort_values(["cohort_cut", "layer", "cohen_d"],
                            ascending=[True, True, False])
df_out.to_csv(OUT / "h10_celltype_d.tsv", sep="\t", index=False)
print(f"[write] {OUT / 'h10_celltype_d.tsv'}  rows={len(df_out)}")

ps = pd.DataFrame.from_dict(per_sample, orient="index")
ps.index.name = "sample_id"
# annotate dm/cohort
ann = master.set_index("sample_id")[["dm", "molecular_subtype", "histology_subtype"]]
ps = ps.join(ann, how="left")
ps.to_csv(OUT / "h10_per_sample_scores.tsv", sep="\t")
print(f"[write] {OUT / 'h10_per_sample_scores.tsv'}  shape={ps.shape}")

# ---------------------------------------------------------------------- print headline
print("\n=== TOP 5 by |d| within BRAF-cPTC master (n=110) ===")
sub = df_out[df_out["cohort_cut"] == "BRAF_cPTC_n110"].copy()
sub["abs_d"] = sub["cohen_d"].abs()
print(sub.sort_values("abs_d", ascending=False)
        .head(8)[["layer", "cell_type", "n_DM1", "n_DM2",
                  "cohen_d", "MW_p", "BH_FDR"]]
        .to_string(index=False))

print("\n=== TLS Cabrita 12-gene ===")
print(df_out[df_out["layer"] == "TLS_signature"]
       [["cohort_cut", "cell_type", "n_DM1", "n_DM2",
         "cohen_d", "MW_p", "BH_FDR", "mean_DM1", "mean_DM2"]]
       .to_string(index=False))
