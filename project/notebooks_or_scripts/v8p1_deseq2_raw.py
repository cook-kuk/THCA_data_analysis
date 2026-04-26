#!/usr/bin/env python3
"""v8.1 Task A — DESeq2 on RAW STAR counts (no pseudo-count).

Fixes v8 T4 degradation: v8 used round((2^log2tpm - 1) * 50) pseudo-counts
because the script looked in /data/thca/v5_cross_cancer/raw/THCA/tcga_rnaseq
(doesn't exist). The files are actually at
/data/thca/data_raw/gdc/TCGA-THCA/counts/*.tsv (572 per-sample STAR outputs).

This script:
  1. Loads manifest to get file_id -> sample_submitter_id mapping
  2. Builds 60,660 gene × 351 sample count matrix (column = unstranded)
  3. Filters to 351 labeled THCA samples (BRAF=293 / RAS=58)
  4. Runs pydeseq2 with design ~subtype
  5. Compares to v8 pseudo-count DESeq2 output
  6. Validates direction on GSE27155 (microarray log2 intensities)

Outputs (results/v8p1_rigor/a_deseq2_raw/):
  deseq2_tcga_raw.tsv        — full gene-level DE table
  comparison_pseudo_vs_raw.tsv — v8 vs v8.1 concordance
  gse27155_validation.tsv    — direction validation on microarray
  summary.md                 — narrative

Log: logs/v8p1_deseq2_raw.log
"""
from __future__ import annotations
import sys, time, gzip
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/opt/thyroid-dash/project")
GDC  = Path("/data/thca/data_raw/gdc/TCGA-THCA")
V5   = Path("/data/thca/v5_cross_cancer/raw/THCA")
OUT  = ROOT / "results" / "v8p1_rigor" / "a_deseq2_raw"
LOG  = ROOT / "logs" / "v8p1_deseq2_raw.log"
OUT.mkdir(parents=True, exist_ok=True)
LOG.parent.mkdir(parents=True, exist_ok=True)

_log_f = open(LOG, "w")
def log(msg: str) -> None:
    stamp = time.strftime("[%Y-%m-%d %H:%M:%S]")
    line = f"{stamp} {msg}"
    print(line, flush=True)
    _log_f.write(line + "\n")
    _log_f.flush()

log("v8.1 Task A — raw-count DESeq2 on THCA")
log("=" * 60)

# ------------------------------------------------------------------
# Step 1: manifest -> file/sample mapping
# ------------------------------------------------------------------
manifest = pd.read_csv(GDC / "tcga_thca_star_counts_manifest.tsv", sep="\t")
log(f"[manifest] {len(manifest)} entries, cols={list(manifest.columns)}")

# file_id corresponds to the subdirectory; file_name is the actual TSV inside
# BUT the cached layout is GDC/counts/<file_id>.tsv (flat, renamed from
# the original bundle). Verify both.
count_files = sorted((GDC / "counts").glob("*.tsv"))
log(f"[counts] {len(count_files)} TSVs in counts/ dir")

# Map by file_id prefix (UUID)
fid2sample = dict(zip(manifest["file_id"], manifest["sample_submitter_id"]))
fid2case = dict(zip(manifest["file_id"], manifest["case_submitter_id"]))

# Each TSV is named <uuid>.tsv — extract UUID and map
mapped = 0
path2sample = {}
for f in count_files:
    uuid = f.stem
    if uuid in fid2sample:
        path2sample[f] = fid2sample[uuid]
        mapped += 1
log(f"[counts] mapped {mapped}/{len(count_files)} files to sample_submitter_id")

# ------------------------------------------------------------------
# Step 2: load labels, restrict to labeled samples
# ------------------------------------------------------------------
labels = pd.read_csv(V5 / "tcga_samples_labels.tsv", sep="\t")
log(f"[labels] {len(labels)} TCGA-THCA labels, subtype counts: "
    f"{dict(labels['label'].value_counts())}")

labeled_samples = set(labels["sample_id"])

# Filter count files to those with sample_submitter_id in labeled set
keep_files = [(p, s) for p, s in path2sample.items() if s in labeled_samples]
log(f"[filter] {len(keep_files)} count files match labeled samples")

# ------------------------------------------------------------------
# Step 3: build count matrix (genes × samples)
# ------------------------------------------------------------------
log("[build] reading per-sample TSVs, extracting unstranded counts")
first = pd.read_csv(keep_files[0][0], sep="\t", comment="#")
# The first 4 rows are N_unmapped/N_multimapping/N_noFeature/N_ambiguous
data_rows = first[first["gene_id"].str.startswith("ENSG")]
gene_ids = data_rows["gene_id"].values
gene_names = data_rows["gene_name"].values
log(f"[build] {len(gene_ids)} ENSG genes per file")

count_mat = np.zeros((len(gene_ids), len(keep_files)), dtype=np.int32)
sample_ids = []
for j, (p, sid) in enumerate(keep_files):
    df = pd.read_csv(p, sep="\t", comment="#")
    df = df[df["gene_id"].str.startswith("ENSG")].set_index("gene_id").loc[gene_ids]
    count_mat[:, j] = df["unstranded"].fillna(0).astype(np.int32).values
    sample_ids.append(sid)
    if (j + 1) % 50 == 0:
        log(f"  loaded {j+1}/{len(keep_files)}")

log(f"[build] count matrix: {count_mat.shape} (genes × samples)")
log(f"[build] non-zero fraction: {(count_mat > 0).mean():.3f}")

# Map to subtype in same order
lbl_map = dict(zip(labels["sample_id"], labels["label"]))
subtypes = np.array([lbl_map[s] for s in sample_ids])
log(f"[build] subtype distribution: "
    f"{dict(zip(*np.unique(subtypes, return_counts=True)))}")

# ------------------------------------------------------------------
# Step 4: pydeseq2
# ------------------------------------------------------------------
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

log("[deseq2] preparing DeseqDataSet (raw integer counts)")
# pydeseq2 expects counts shape (samples, genes) — transpose
counts_df = pd.DataFrame(
    count_mat.T, index=sample_ids,
    columns=gene_ids,
)
# Map ENSG to gene_name, but keep ENSG as row identifier to allow duplicates
meta_df = pd.DataFrame({"subtype": subtypes}, index=sample_ids)

t0 = time.time()
dds = DeseqDataSet(
    counts=counts_df,
    metadata=meta_df,
    design_factors="subtype",
    refit_cooks=True,
    inference=DefaultInference(n_cpus=8),
    quiet=False,
)
dds.deseq2()
log(f"[deseq2] dds.deseq2() done in {time.time()-t0:.1f}s")

t0 = time.time()
stat = DeseqStats(dds, contrast=["subtype", "BRAF", "RAS"],
                  inference=DefaultInference(n_cpus=8))
stat.summary()
log(f"[deseq2] summary() done in {time.time()-t0:.1f}s")
res = stat.results_df.copy()
res.index.name = "gene_id"
res["gene_name"] = pd.Series(dict(zip(gene_ids, gene_names)))

# Significance mask
res["significant"] = (res["padj"] < 0.05) & (res["log2FoldChange"].abs() > 1)
sig_n = int(res["significant"].sum())
log(f"[deseq2] RAW significant genes (FDR<0.05, |log2FC|>1): {sig_n}")

out_tsv = OUT / "deseq2_tcga_raw.tsv"
res.reset_index().to_csv(out_tsv, sep="\t", index=False,
                         float_format="%.6g")
log(f"[deseq2] wrote {out_tsv}")

# ------------------------------------------------------------------
# Step 5: compare to v8 pseudo-count DESeq2
# ------------------------------------------------------------------
V8 = ROOT / "results" / "v8_statgen" / "v8_biomarker_DE_recomputation.tsv"
v8_de = pd.read_csv(V8, sep="\t")
log(f"[compare] v8 DE table rows: {len(v8_de)}, cols: {list(v8_de.columns)}")

# The v8 table is keyed by HGNC 'gene'. Our raw table is keyed by ENSG.
# Map via gene_name.
raw_by_hgnc = res.dropna(subset=["gene_name"]).copy()
raw_by_hgnc = raw_by_hgnc[~raw_by_hgnc["gene_name"].duplicated()].set_index("gene_name")

v8_by_hgnc = v8_de.set_index("gene")
common = raw_by_hgnc.index.intersection(v8_by_hgnc.index)
log(f"[compare] gene-name intersection (raw vs v8 pseudo): {len(common)}")

# Direction agreement: same sign of log2FC
raw_lfc = raw_by_hgnc.loc[common, "log2FoldChange"].fillna(0)
v8_lfc = v8_by_hgnc.loc[common, "new_log2FC"].fillna(0)
same_sign = ((raw_lfc > 0) == (v8_lfc > 0)).sum()
direction_frac = same_sign / len(common)
log(f"[compare] same-sign log2FC: {same_sign}/{len(common)} ({100*direction_frac:.1f}%)")

# Gene overlap: significance co-occurrence
raw_sig = raw_by_hgnc.loc[common, "significant"].fillna(False)
v8_sig = v8_by_hgnc.loc[common, "new_significant"].fillna(False)
both_sig = int(((raw_sig) & (v8_sig)).sum())
raw_only = int(((raw_sig) & (~v8_sig)).sum())
v8_only = int(((~raw_sig) & (v8_sig)).sum())
log(f"[compare] both_sig={both_sig}  raw_only={raw_only}  v8_only={v8_only}")

comp_df = pd.DataFrame({
    "metric": ["gene_intersection", "same_sign_log2FC", "direction_fraction",
               "sig_in_both", "raw_only_sig", "v8_pseudo_only_sig",
               "raw_sig_total", "v8_pseudo_sig_total"],
    "value": [len(common), same_sign, f"{100*direction_frac:.2f}%",
              both_sig, raw_only, v8_only,
              int(raw_sig.sum()), int(v8_sig.sum())],
})
comp_df.to_csv(OUT / "comparison_pseudo_vs_raw.tsv", sep="\t", index=False)
log(f"[compare] wrote {OUT / 'comparison_pseudo_vs_raw.tsv'}")

# ------------------------------------------------------------------
# Step 6: GSE27155 direction validation
# ------------------------------------------------------------------
# Load GSE27155 log2 matrix + labels and compute mean(BRAF) - mean(RAS)
# For each DESeq2-significant gene, check sign agreement.
GEO = Path("/data/thca/v5_cross_cancer/raw/THCA/geo")
log(f"[gse] looking for GSE27155 under {GEO}")
geo_files = list(GEO.glob("GSE27155*"))
log(f"[gse] {len(geo_files)} files match: {[f.name for f in geo_files]}")

# Fallback — search
if not geo_files:
    geo_files = list(Path("/data/thca/data_processed/microarray").glob("GSE27155*"))
    log(f"[gse] fallback search found {len(geo_files)}: "
        f"{[f.name for f in geo_files]}")
if not geo_files:
    # Try microarray_v3
    geo_files = list(Path("/data/thca/data_processed/microarray_v3").glob("GSE27155*"))
    log(f"[gse] v3 search: {[f.name for f in geo_files]}")

try:
    # The expression matrix is typically a TSV with probes x samples; we want
    # gene-level. Find it and the matching label file.
    expr_file = next((f for f in geo_files if "log2" in f.name or "expr" in f.name), None)
    label_file = next((f for f in geo_files if "label" in f.name or "meta" in f.name), None)
    if expr_file is None or label_file is None:
        raise RuntimeError(f"could not identify expr/label files in {geo_files}")
    log(f"[gse] using expr={expr_file.name} label={label_file.name}")

    # Load
    if expr_file.suffix == ".gz":
        expr = pd.read_csv(expr_file, sep="\t", index_col=0, compression="gzip")
    else:
        expr = pd.read_csv(expr_file, sep="\t", index_col=0)
    log(f"[gse] expr shape: {expr.shape}")

    lbl = pd.read_csv(label_file, sep="\t")
    log(f"[gse] label cols: {list(lbl.columns)}")

    # Best-effort: find a label column and sample id column
    lbl_cols = [c for c in lbl.columns if "label" in c.lower() or
                "subtype" in c.lower() or "mutation" in c.lower()]
    sid_cols = [c for c in lbl.columns if "sample" in c.lower() or
                "gsm" in c.lower() or c.lower() in ("id", "name")]
    log(f"[gse] candidate label cols: {lbl_cols}")
    log(f"[gse] candidate sample-id cols: {sid_cols}")

    if lbl_cols and sid_cols:
        lbl_col = lbl_cols[0]
        sid_col = sid_cols[0]
        lbl = lbl[[sid_col, lbl_col]].dropna()
        lbl.columns = ["sample", "label"]
        # Normalise
        lbl["label"] = lbl["label"].str.upper().str.extract(r"(BRAF|RAS)", expand=False)
        lbl = lbl.dropna()
        log(f"[gse] usable labels: {lbl['label'].value_counts().to_dict()}")
    else:
        raise RuntimeError("could not identify label/sample columns")

    # Align
    cohort_samples = [s for s in lbl["sample"] if s in expr.columns]
    expr_sub = expr[cohort_samples]
    lbl_sub = lbl.set_index("sample").loc[cohort_samples, "label"]

    braf_mean = expr_sub.loc[:, (lbl_sub == "BRAF").values].mean(axis=1)
    ras_mean  = expr_sub.loc[:, (lbl_sub == "RAS").values].mean(axis=1)
    geo_delta = braf_mean - ras_mean

    # Align to raw DESeq2 significant genes via gene_name
    sig_genes = raw_by_hgnc[raw_by_hgnc["significant"]].index
    # Expression matrix probably has HGNC row names; check
    matching = sig_genes.intersection(expr.index)
    log(f"[gse] raw-sig ∩ GSE27155 gene universe: {len(matching)}")

    if len(matching) > 100:
        tcga_sign = (raw_by_hgnc.loc[matching, "log2FoldChange"] > 0).values
        geo_sign  = (geo_delta.loc[matching] > 0).values
        agree = (tcga_sign == geo_sign).sum()
        concordance = agree / len(matching)
        log(f"[gse] BRAF↑/RAS↑ direction concordance: "
            f"{agree}/{len(matching)} ({100*concordance:.1f}%)")
    else:
        log(f"[gse] too few overlapping genes ({len(matching)}); "
            "skipping concordance")
        concordance = None

    val_df = pd.DataFrame({
        "metric": ["sig_genes_raw_DESeq2", "gse27155_gene_universe",
                   "gse27155_overlap", "BRAF_up_direction_concordance"],
        "value":  [int(raw_by_hgnc["significant"].sum()),
                   expr.shape[0], len(matching),
                   f"{100*concordance:.1f}%" if concordance is not None else "NA"],
    })
except Exception as e:
    log(f"[gse] VALIDATION FAILED: {type(e).__name__}: {e}")
    val_df = pd.DataFrame({
        "metric": ["gse27155_validation_status"],
        "value":  [f"failed: {type(e).__name__}: {e}"],
    })

val_df.to_csv(OUT / "gse27155_validation.tsv", sep="\t", index=False)
log(f"[gse] wrote {OUT / 'gse27155_validation.tsv'}")

log("v8.1 Task A complete")
_log_f.close()
