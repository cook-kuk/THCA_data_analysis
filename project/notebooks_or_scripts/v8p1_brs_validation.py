#!/usr/bin/env python3
"""v8.1 — BRS-based cross-platform direction validation.

Goal: extend the GSE27155 sign-concordance test (83.4 %) to two additional
cohorts (GSE33630 PTC n=49, GSE29265 PTC n=20) that lack BRAF/RAS mutation
labels but have expression matrices on disk. We use the Chakravarty 2011
BRS52 gene-expression signature to classify each PTC sample as BRAF-like
or RAS-like, validate the signature on TCGA-THCA where mutation truth is
known, then for each BRS-labeled sub-cohort compute the per-gene
mean(BRAF-like) − mean(RAS-like) and check sign-agreement with the v8.1
raw-count DESeq2 log2FC.

Non-circularity guard: BRS labels are used ONLY for downstream direction
validation. They are NOT fed into DIAL or v5.1 LODO computations. The
Chakravarty BRS52 gene panel is also separate from the v5.1 3 000 top-
variance LODO feature pool (the panel is fixed, defined ex ante from
Cancer Cell 2011 — independent of our biomarker discovery).

Outputs (results/v8p1_rigor/e_brs_validation/):
  brs_tcga_validation.tsv          — BRS classification accuracy on TCGA
  brs_labels_GSE33630.tsv          — per-sample BRS label + score
  brs_labels_GSE29265.tsv          — same
  three_cohort_concordance.tsv     — sign agreement for each of 3 cohorts
                                      against TCGA raw-DESeq2 log2FC
  brs_validation_analysis.md       — narrative

Log: logs/v8p1_brs_validation.log
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/opt/thyroid-dash/project")
RES = ROOT / "results" / "v8p1_rigor"
OUT = RES / "e_brs_validation"
OUT.mkdir(parents=True, exist_ok=True)
LOG = ROOT / "logs" / "v8p1_brs_validation.log"

_log_f = open(LOG, "w")
def log(msg: str) -> None:
    stamp = time.strftime("[%Y-%m-%d %H:%M:%S]")
    print(f"{stamp} {msg}", flush=True)
    _log_f.write(f"{stamp} {msg}\n"); _log_f.flush()

log("v8.1 BRS-based 3-cohort cross-platform validation")
log("=" * 60)

# --------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------
BRS_FILE = ROOT / "results" / "v12_literature" / "crosscheck" / \
           "reference_gene_lists" / "Chakravarty_2011_BRS52.tsv"
TCGA_TPM = Path("/data/thca/v5_cross_cancer/raw/THCA/tcga_log2tpm.tsv.gz")
TCGA_LBL = Path("/data/thca/v5_cross_cancer/raw/THCA/tcga_samples_labels.tsv")
GEO33    = Path("/data/thca/data_processed/microarray_v3/GSE33630_v3_log2.tsv")
GEO29    = Path("/data/thca/data_processed/microarray_v3/GSE29265_v3_log2.tsv")
META33   = Path("/data/thca/data_raw/v3_ext/GSE33630/GSE33630_metadata.tsv")
META29   = Path("/data/thca/data_raw/v3_ext/GSE29265/GSE29265_metadata.tsv")
RAW_DE   = RES / "a_deseq2_raw" / "deseq2_tcga_raw.tsv"

brs = pd.read_csv(BRS_FILE, sep="\t")["gene"].tolist()
log(f"[brs] {len(brs)} BRS genes (Chakravarty 2011 BRS52)")

# --------------------------------------------------------------------
# 1. Load TCGA-THCA log2-TPM, restrict to BRS, split by mutation label
# --------------------------------------------------------------------
tpm = pd.read_csv(TCGA_TPM, sep="\t", index_col=0)
log(f"[tcga] log2-TPM shape: {tpm.shape}")
labels = pd.read_csv(TCGA_LBL, sep="\t").set_index("sample_id")["label"]

brs_in_tcga = [g for g in brs if g in tpm.index]
log(f"[tcga] BRS genes available: {len(brs_in_tcga)}/{len(brs)}")

# Use only BRS genes for BRS-score computation
tpm_brs = tpm.loc[brs_in_tcga]
samples_in_tpm = [s for s in tpm.columns if s in labels.index]
log(f"[tcga] samples with mutation label: {len(samples_in_tpm)}")

tpm_brs_lbl = tpm_brs[samples_in_tpm]
y_tcga = labels.loc[samples_in_tpm]
log(f"[tcga] subtype counts: {y_tcga.value_counts().to_dict()}")

# Compute per-gene centroids
centroid_braf = tpm_brs_lbl.loc[:, (y_tcga == "BRAF").values].mean(axis=1)
centroid_ras  = tpm_brs_lbl.loc[:, (y_tcga == "RAS").values].mean(axis=1)
log(f"[centroid] BRAF n={(y_tcga=='BRAF').sum()}  RAS n={(y_tcga=='RAS').sum()}")
log(f"[centroid] BRAF mean range: [{centroid_braf.min():.2f}, {centroid_braf.max():.2f}]")
log(f"[centroid] RAS  mean range: [{centroid_ras.min():.2f}, {centroid_ras.max():.2f}]")

def brs_score_from_matrix(expr: pd.DataFrame, c_braf: pd.Series,
                          c_ras: pd.Series) -> pd.DataFrame:
    """Per-sample BRS = corr(sample, c_braf) - corr(sample, c_ras).
       Positive → BRAF-like. Returns DataFrame with columns:
       brs_score, corr_braf, corr_ras, brs_label."""
    common = expr.index.intersection(c_braf.index)
    expr_m = expr.loc[common]
    cb = c_braf.loc[common]
    cr = c_ras.loc[common]
    # correlate each sample (column) with each centroid (1D)
    # Pearson via numpy
    Xc = expr_m.values - expr_m.values.mean(axis=0, keepdims=True)
    cb_c = cb.values - cb.values.mean()
    cr_c = cr.values - cr.values.mean()
    norm_x = np.sqrt((Xc**2).sum(axis=0))
    norm_b = np.sqrt((cb_c**2).sum())
    norm_r = np.sqrt((cr_c**2).sum())
    corr_b = (Xc * cb_c[:, None]).sum(axis=0) / (norm_x * norm_b + 1e-12)
    corr_r = (Xc * cr_c[:, None]).sum(axis=0) / (norm_x * norm_r + 1e-12)
    score = corr_b - corr_r
    label = np.where(score > 0, "BRAF_like", "RAS_like")
    return pd.DataFrame({
        "sample": expr_m.columns,
        "corr_braf": corr_b, "corr_ras": corr_r,
        "brs_score": score, "brs_label": label,
    }).set_index("sample")

# --------------------------------------------------------------------
# 2. Validate BRS classification on TCGA (where mutation labels are known)
# --------------------------------------------------------------------
brs_tcga = brs_score_from_matrix(tpm_brs_lbl, centroid_braf, centroid_ras)
brs_tcga["mutation_label"] = y_tcga.values

# Concordance: BRS BRAF_like ↔ mutation BRAF; BRS RAS_like ↔ mutation RAS
agree = ((brs_tcga["brs_label"] == "BRAF_like") &
         (brs_tcga["mutation_label"] == "BRAF")).sum() + \
        ((brs_tcga["brs_label"] == "RAS_like") &
         (brs_tcga["mutation_label"] == "RAS")).sum()
total = len(brs_tcga)
log(f"[tcga-validation] BRS-vs-mutation accuracy: {agree}/{total} ({100*agree/total:.1f}%)")
log(f"[tcga-validation] confusion:\n"
    f"{pd.crosstab(brs_tcga['mutation_label'], brs_tcga['brs_label']).to_string()}")

brs_tcga.to_csv(OUT / "brs_tcga_validation.tsv", sep="\t")

# --------------------------------------------------------------------
# 3. Apply BRS to GSE33630 + GSE29265 (PTC subset only)
# --------------------------------------------------------------------
def parse_ptc(meta_path: Path, expr_path: Path,
              tissue_col: str = "characteristics_ch1") -> pd.DataFrame:
    """Return expr matrix restricted to PTC samples."""
    meta = pd.read_csv(meta_path, sep="\t")
    char = meta[tissue_col].astype(str).str.lower()
    is_ptc = char.str.contains("papillary thyroid carcinoma", na=False)
    ptc_samples = meta.loc[is_ptc, "geo_accession"].tolist()
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    cols = [c for c in expr.columns if c in ptc_samples]
    return expr[cols]

ptc33 = parse_ptc(META33, GEO33)
ptc29 = parse_ptc(META29, GEO29)
log(f"[geo] GSE33630 PTC samples: {ptc33.shape[1]}")
log(f"[geo] GSE29265 PTC samples: {ptc29.shape[1]}")

brs33 = brs_score_from_matrix(ptc33, centroid_braf, centroid_ras)
brs29 = brs_score_from_matrix(ptc29, centroid_braf, centroid_ras)
brs33.to_csv(OUT / "brs_labels_GSE33630.tsv", sep="\t")
brs29.to_csv(OUT / "brs_labels_GSE29265.tsv", sep="\t")
log(f"[geo] GSE33630 BRS labels: {brs33['brs_label'].value_counts().to_dict()}")
log(f"[geo] GSE29265 BRS labels: {brs29['brs_label'].value_counts().to_dict()}")

# --------------------------------------------------------------------
# 4. Three-cohort direction concordance vs raw-DESeq2 log2FC
# --------------------------------------------------------------------
de = pd.read_csv(RAW_DE, sep="\t")
de = de.dropna(subset=["gene_name"]).drop_duplicates("gene_name")
de_sig = de[de["significant"]].set_index("gene_name")
log(f"[de] raw-count DESeq2 sig: {len(de_sig)}")

# Per-cohort: for each gene, compute mean(BRAF_like) − mean(RAS_like)
def cohort_delta(expr: pd.DataFrame, brs_lbl: pd.DataFrame) -> pd.Series:
    braf_samples = brs_lbl[brs_lbl["brs_label"] == "BRAF_like"].index
    ras_samples  = brs_lbl[brs_lbl["brs_label"] == "RAS_like"].index
    if len(braf_samples) < 2 or len(ras_samples) < 2:
        return pd.Series(dtype=float)
    return (expr.loc[:, braf_samples].mean(axis=1) -
            expr.loc[:, ras_samples].mean(axis=1))

delta33 = cohort_delta(ptc33, brs33)
delta29 = cohort_delta(ptc29, brs29)
log(f"[delta] GSE33630 delta computed: {len(delta33)} genes")
log(f"[delta] GSE29265 delta computed: {len(delta29)} genes")

# Also recompute GSE27155 mutation-truth delta for direct apples-to-apples
geo_pheno = pd.read_csv("/data/thca/v5_cross_cancer/raw/THCA/geo/GSE27155_pheno.tsv",
                        sep="\t").dropna(subset=["label"])
geo_expr = pd.read_csv("/data/thca/v5_cross_cancer/raw/THCA/geo/GSE27155_expr.tsv",
                       sep="\t", index_col=0)
ph = geo_pheno[geo_pheno["label"].isin(["BRAF","RAS"])]
braf27 = geo_expr[ph[ph["label"]=="BRAF"]["sample_id"].tolist()]
ras27  = geo_expr[ph[ph["label"]=="RAS"]["sample_id"].tolist()]
delta27 = braf27.mean(axis=1) - ras27.mean(axis=1)
log(f"[delta] GSE27155 (mutation-truth, n_BRAF=28 n_RAS=13): {len(delta27)} genes")

# Concordance per cohort
def concord(delta: pd.Series, label: str) -> tuple:
    common = delta.index.intersection(de_sig.index)
    if not len(common):
        return label, 0, 0, 0.0
    tcga_sign = (de_sig.loc[common, "log2FoldChange"] > 0).values
    geo_sign  = (delta.loc[common] > 0).values
    agree_n = int((tcga_sign == geo_sign).sum())
    return label, len(common), agree_n, agree_n/len(common)*100

rows = [
    concord(delta27, "GSE27155 (mutation truth)"),
    concord(delta33, "GSE33630 (BRS-inferred)"),
    concord(delta29, "GSE29265 (BRS-inferred)"),
]
conc_df = pd.DataFrame(rows, columns=[
    "cohort", "sig_overlap_n", "same_sign_n", "concordance_pct",
])
conc_df["concordance_pct"] = conc_df["concordance_pct"].round(2)
log(f"[concordance]\n{conc_df.to_string(index=False)}")

# Three-cohort consensus: gene level, in how many cohorts does the sign
# agree with TCGA?
common_all = (de_sig.index
              .intersection(delta27.index)
              .intersection(delta33.index)
              .intersection(delta29.index))
log(f"[consensus] genes in all 4 (TCGA+27155+33630+29265): {len(common_all)}")

if len(common_all):
    tcga_sgn = (de_sig.loc[common_all, "log2FoldChange"] > 0).values
    s27 = (delta27.loc[common_all] > 0).values
    s33 = (delta33.loc[common_all] > 0).values
    s29 = (delta29.loc[common_all] > 0).values
    n_match = (tcga_sgn==s27).astype(int) + (tcga_sgn==s33).astype(int) + \
              (tcga_sgn==s29).astype(int)
    consensus = pd.DataFrame({
        "gene": list(common_all),
        "match_count_3cohort": n_match,
    })
    log(f"[consensus] match-count distribution:\n"
        f"{consensus['match_count_3cohort'].value_counts().sort_index().to_string()}")
    log(f"[consensus] genes matching TCGA in all 3 GEO cohorts: "
        f"{(n_match==3).sum()}  ({100*(n_match==3).sum()/len(common_all):.1f}%)")
    log(f"[consensus] genes matching in ≥2 of 3:               "
        f"{(n_match>=2).sum()}  ({100*(n_match>=2).sum()/len(common_all):.1f}%)")
    consensus.to_csv(OUT / "three_cohort_per_gene_consensus.tsv",
                     sep="\t", index=False)

conc_df.to_csv(OUT / "three_cohort_concordance.tsv", sep="\t", index=False)
log(f"[write] {OUT / 'three_cohort_concordance.tsv'}")

log("v8.1 BRS validation complete")
_log_f.close()
