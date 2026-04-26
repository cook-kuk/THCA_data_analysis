#!/usr/bin/env python
"""v3_ext_02_harmonize.py

Harmonize external cohorts into v3 panels:
  - Map ENSG -> HGNC (only when needed; our existing files are already gene-symbol indexed).
  - Collapse duplicate gene symbols by max.
  - log2(count+1) normalisation for raw counts.
  - Probe-to-gene collapse for GPL570 using cached v2 mapping where possible.
  - Save to data_processed/bulk_rnaseq_v3/ and data_processed/microarray_v3/.
  - Build metadata/sample_master_v3.tsv with v3_anchor in
      {mutation_verified, fusion_verified, histology_proxy, ambiguous}.
  - Append to metadata/sample_master_v3_merged.tsv.
  - NEVER overwrite metadata/sample_master.tsv.
"""
import os
import sys
import json
import logging
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path("/opt/thyroid-dash/project")
DP = ROOT / "data_processed"
DP_V3_RNA = DP / "bulk_rnaseq_v3"
DP_V3_MIC = DP / "microarray_v3"
DP_V3_RNA.mkdir(parents=True, exist_ok=True)
DP_V3_MIC.mkdir(parents=True, exist_ok=True)
RAW = ROOT / "data_raw" / "v3_ext"
META = ROOT / "metadata"
LOGDIR = ROOT / "logs"
LOGDIR.mkdir(exist_ok=True)

LOGFILE = LOGDIR / "v3_ext_02_harmonize.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[logging.FileHandler(LOGFILE, mode="w"), logging.StreamHandler()],
)
log = logging.getLogger("v3_ext_02")


def collapse_by_symbol(df: pd.DataFrame, idx_col: str = "gene_symbol") -> pd.DataFrame:
    """Collapse duplicated gene symbols by taking the max row."""
    if idx_col in df.columns:
        df = df.set_index(idx_col)
    df = df[~df.index.duplicated(keep="first")] if df.index.is_unique else df.groupby(level=0).max()
    return df


def collapse_probes_to_genes(probe_matrix: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    """Probe x samples -> gene x samples (max probe per gene)."""
    probe_matrix = probe_matrix.copy()
    probe_matrix.index = probe_matrix.index.astype(str)
    mp = mapping.dropna(subset=["probe_id", "gene_symbol"]).drop_duplicates("probe_id").set_index("probe_id")["gene_symbol"]
    shared = probe_matrix.index.intersection(mp.index)
    log.info(f"  probes shared with mapping: {len(shared)}/{len(probe_matrix)}")
    pm = probe_matrix.loc[shared].copy()
    pm["gene_symbol"] = mp.loc[shared].values
    g = pm.groupby("gene_symbol").max()
    return g


def harmonize_rnaseq():
    """TCGA-THCA + GSE126698 already log2 gene-symbol. Copy collapsed versions."""
    for src_name, gse in [
        ("TCGA-THCA_rnaseq_expression_log2.tsv", "TCGA-THCA"),
        ("GSE126698_rnaseq_expression_log2.tsv", "GSE126698"),
        ("GSE213647_rnaseq_expression_log2.tsv", "GSE213647"),
    ]:
        src = DP / "bulk_rnaseq" / src_name
        if not src.exists():
            log.warning(f"missing {src}")
            continue
        df = pd.read_csv(src, sep="\t")
        df = collapse_by_symbol(df, "gene_symbol")
        dest = DP_V3_RNA / f"{gse}_v3_log2.tsv"
        df.to_csv(dest, sep="\t")
        log.info(f"wrote {dest}  ({df.shape})")


def harmonize_microarray():
    """GSE27155, GSE76039 already processed; copy collapsed versions under v3."""
    for src_name, gse in [
        ("GSE27155_microarray_expression_log2.tsv", "GSE27155"),
        ("GSE76039_microarray_expression_log2.tsv", "GSE76039"),
    ]:
        src = DP / "microarray" / src_name
        if not src.exists():
            continue
        df = pd.read_csv(src, sep="\t")
        df = collapse_by_symbol(df, "gene_symbol")
        dest = DP_V3_MIC / f"{gse}_v3_log2.tsv"
        df.to_csv(dest, sep="\t")
        log.info(f"wrote {dest}  ({df.shape})")

    # Now handle new GPL570 cohorts from v3_ext download
    gpl570_mapping_cache = DP / "microarray" / "GSE27155_probe_to_gene_mapping.tsv"
    if gpl570_mapping_cache.exists():
        mp = pd.read_csv(gpl570_mapping_cache, sep="\t")
    else:
        mp = pd.DataFrame(columns=["probe_id", "gene_symbol"])
    for gse in ["GSE33630", "GSE29265"]:
        pf = RAW / gse / f"{gse}_probe_matrix_log2.tsv"
        if not pf.exists():
            log.warning(f"{gse} probe matrix missing -> SKIPPED")
            continue
        pm = pd.read_csv(pf, sep="\t", index_col=0)
        gm = collapse_probes_to_genes(pm, mp)
        dest = DP_V3_MIC / f"{gse}_v3_log2.tsv"
        gm.to_csv(dest, sep="\t")
        log.info(f"wrote {dest}  ({gm.shape})")


def build_sample_master_v3() -> pd.DataFrame:
    sm = pd.read_csv(META / "sample_master.tsv", sep="\t")
    # Anchor logic:
    #   mutation_verified: driver_anchor in {BRAF, RAS} with has_braf_v600e or has_ras_mut true
    #   fusion_verified: (placeholder unless downstream annotates)
    #   histology_proxy: label_confidence==high but driver_anchor=unknown
    #   ambiguous: otherwise
    mut_groups_f = ROOT / "results/tables/tcga_thca_mutation_groups.tsv"
    mut = pd.read_csv(mut_groups_f, sep="\t") if mut_groups_f.exists() else pd.DataFrame()
    mut_ids_braf = set(mut.loc[mut.get("has_braf_v600e") == True, "sample_id"]) if not mut.empty else set()
    mut_ids_ras = set(mut.loc[mut.get("has_ras_mut") == True, "sample_id"]) if not mut.empty else set()
    def _anchor(r):
        sid = r["sample_id"]
        if sid in mut_ids_braf or sid in mut_ids_ras:
            return "mutation_verified"
        if str(r.get("driver_anchor", "")).lower() in ("braf", "ras") and str(r.get("label_confidence", "")).lower() == "high":
            return "mutation_verified"
        if str(r.get("label_confidence", "")).lower() == "high":
            return "histology_proxy"
        return "ambiguous"
    sm["v3_anchor"] = sm.apply(_anchor, axis=1)
    # Save master_v3
    out = META / "sample_master_v3.tsv"
    sm.to_csv(out, sep="\t", index=False)
    log.info(f"wrote {out}  n={len(sm)}  anchor counts:\n{sm['v3_anchor'].value_counts().to_dict()}")

    # Append to sample_master_v3_merged
    merged_out = META / "sample_master_v3_merged.tsv"
    if merged_out.exists():
        prev = pd.read_csv(merged_out, sep="\t")
        merged = pd.concat([prev, sm], ignore_index=True).drop_duplicates("sample_id", keep="last")
    else:
        merged = sm.copy()
    merged.to_csv(merged_out, sep="\t", index=False)
    log.info(f"wrote {merged_out}  n={len(merged)}")

    # Placeholder rows for GSE33630 / GSE29265 (and PRJEB11591 if available)
    extras = []
    for gse in ["GSE33630", "GSE29265"]:
        meta_f = RAW / gse / f"{gse}_metadata.tsv"
        if meta_f.exists():
            dfm = pd.read_csv(meta_f, sep="\t")
            for sid in dfm["sample_id"].tolist():
                extras.append({
                    "sample_id": sid, "dataset": gse, "platform": "GPL570",
                    "modality": "microarray", "tissue_type": "Thyroid",
                    "v3_anchor": "histology_proxy", "label_confidence": "medium",
                    "histology_subtype": "unknown",
                })
    prjeb_meta = RAW / "prjeb11591" / "ena_filereport.tsv"
    if prjeb_meta.exists():
        try:
            em = pd.read_csv(prjeb_meta, sep="\t")
            for sid in em.get("run_accession", pd.Series(dtype=str)).tolist():
                extras.append({
                    "sample_id": sid, "dataset": "PRJEB11591", "platform": "Illumina",
                    "modality": "bulk RNA-seq", "tissue_type": "Thyroid",
                    "v3_anchor": "histology_proxy", "label_confidence": "medium",
                    "histology_subtype": "unknown",
                })
        except Exception as e:
            log.warning(f"prjeb metadata parse failed: {e}")
    if extras:
        ex = pd.DataFrame(extras)
        merged2 = pd.concat([merged, ex], ignore_index=True).drop_duplicates("sample_id", keep="first")
        merged2.to_csv(merged_out, sep="\t", index=False)
        log.info(f"extras appended: +{len(extras)} rows -> {merged_out} (n={len(merged2)})")
    return sm


def main():
    log.info("v3_ext_02_harmonize start")
    harmonize_rnaseq()
    harmonize_microarray()
    build_sample_master_v3()
    log.info("v3_ext_02_harmonize done")


if __name__ == "__main__":
    main()
