#!/usr/bin/env python3
"""Check immune-module gene coverage per downloaded dataset.

For each dataset, extract gene/probe identifiers from the series matrix
expression table (or from supplementary processed matrices when applicable),
then intersect with each module's gene list and report coverage.

For microarray (GPL570 = Affy HG U133 Plus 2): row IDs are probe IDs.
We do a probe→symbol mapping using the platform annotation in the series
matrix is typically not embedded; instead we report COVERAGE_NEEDS_PROBE_MAP.

For RNA-seq datasets (Hugo, Riaz, GSE126698 supplementary CSV/XLSX), we read
the supplementary processed expression and intersect symbols directly.

Output:
  paper3_ici_dataset_gene_coverage.tsv
"""
from __future__ import annotations
import gzip
import io
import os
import re
from collections import defaultdict
from pathlib import Path

DATA_DIR = Path("/home/seungho/personal/THCA_data_analysis/project/results/paper3_ici_public_data")
REG_DIR  = Path("/home/seungho/personal/THCA_data_analysis/project/reports/paper3_ici")
OUT_DIR  = Path("/home/seungho/personal/THCA_data_analysis/project/results/paper3_ici_data_registry")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Load module gene list
MODULES = defaultdict(list)
with open(REG_DIR / "paper3_ici_module_gene_list.tsv") as fh:
    next(fh)  # header
    for line in fh:
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 2: continue
        MODULES[cols[0]].append(cols[1])

ALL_MODULE_GENES = sorted({g for genes in MODULES.values() for g in genes})

def extract_genes_from_series_matrix(path: Path, max_rows: int = 200000) -> set[str]:
    """Read the data block after !series_matrix_table_begin. First column = ID."""
    genes = set()
    in_table = False
    n = 0
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if not in_table:
                if line.startswith("!series_matrix_table_begin"):
                    in_table = True
                    next(fh, None)  # header row
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            cols = line.rstrip("\n").split("\t")
            if cols and cols[0]:
                genes.add(cols[0].strip().strip('"'))
            n += 1
            if n > max_rows: break
    return genes

def extract_genes_from_csvgz(path: Path, max_rows: int = 200000) -> set[str]:
    genes = set()
    with gzip.open(path, "rt", errors="replace") as fh:
        next(fh, None)  # header
        for i, line in enumerate(fh):
            if i > max_rows: break
            cols = line.rstrip("\n").split(",")
            if cols and cols[0]:
                # may be quoted ensembl ID + symbol; take first
                tok = cols[0].strip().strip('"')
                # try second col if first is ensembl
                if re.match(r"^ENSG\d+", tok) and len(cols) > 1:
                    tok2 = cols[1].strip().strip('"')
                    if tok2 and not re.match(r"^ENSG", tok2):
                        genes.add(tok2)
                        continue
                genes.add(tok)
    return genes

# Datasets to scan: (id, file_path, file_kind)
TARGETS = [
    ("GSE76039",  DATA_DIR/"GSE76039/GSE76039_series_matrix.txt.gz",  "series_matrix"),
    ("GSE65144",  DATA_DIR/"GSE65144/GSE65144_series_matrix.txt.gz",  "series_matrix"),
    ("GSE29265",  DATA_DIR/"GSE29265/GSE29265_series_matrix.txt.gz",  "series_matrix"),
    ("GSE33630",  DATA_DIR/"GSE33630/GSE33630_series_matrix.txt.gz",  "series_matrix"),
    ("GSE60542",  DATA_DIR/"GSE60542/GSE60542_series_matrix.txt.gz",  "series_matrix"),
    ("GSE151179", DATA_DIR/"GSE151179/GSE151179_series_matrix.txt.gz","series_matrix"),
    ("GSE126698", DATA_DIR/"GSE126698/GSE126698_DE_Thyroid_totalRNA_all.csv.gz","csvgz"),
    ("GSE91061",  DATA_DIR/"GSE91061/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz", "csvgz"),
    ("GSE193581", DATA_DIR/"GSE193581/GSE193581_bulkcellline.gz", "csvgz"),
]

# Affy HG U133 Plus 2 platform — small list of representative probes per symbol of interest
# (probe→symbol full annotation deferred to Track B; here we mark microarray as needs-probe-map)
MICROARRAY_NEEDS_PROBE_MAP = True

def main():
    rows = []
    for ds, path, kind in TARGETS:
        if not path.exists():
            rows.append((ds, "MISSING", "", "", ""))
            continue
        if kind == "series_matrix":
            ids = extract_genes_from_series_matrix(path)
        else:
            ids = extract_genes_from_csvgz(path)
        # detect ID type
        sample_ids = list(ids)[:50]
        is_probe = any(re.match(r"^\d+_(at|s_at|x_at|i_at)$", x) for x in sample_ids)
        is_ensembl = any(re.match(r"^ENSG\d+", x) for x in sample_ids)
        is_symbol = any(re.match(r"^[A-Z][A-Z0-9-]+$", x) for x in sample_ids)
        id_type = "probe_affy" if is_probe else "ensembl" if is_ensembl else "symbol_or_other"
        # coverage per module
        module_cov = {}
        for mod, genes in MODULES.items():
            if id_type == "probe_affy":
                module_cov[mod] = "needs_probe_to_symbol_mapping_at_track_B"
            else:
                hit = sum(1 for g in genes if g in ids)
                module_cov[mod] = f"{hit}/{len(genes)}"
        rows.append((ds, id_type, len(ids), module_cov, kind))

    out = OUT_DIR / "paper3_ici_dataset_gene_coverage.tsv"
    cols = ["dataset_id","file_kind","id_type","n_features"] + list(MODULES.keys())
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for ds, idt, nf, mod_cov, kind in rows:
            cells = [ds, kind, idt, str(nf)] + [mod_cov[m] for m in MODULES]
            fh.write("\t".join(cells) + "\n")
    print(f"wrote {out}")
    print(f"n_datasets_scanned={len(rows)} module_count={len(MODULES)} module_genes_total={len(ALL_MODULE_GENES)}")

if __name__ == "__main__":
    main()
