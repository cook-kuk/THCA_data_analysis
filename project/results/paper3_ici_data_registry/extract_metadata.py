#!/usr/bin/env python3
"""Extract per-sample metadata + dataset QC from downloaded GEO series matrices.

Output:
  paper3_ici_public_sample_metadata.tsv
  paper3_ici_public_dataset_qc_summary.tsv
"""
from __future__ import annotations
import gzip
import os
import re
import sys
from pathlib import Path

DATA_DIR = Path("/home/seungho/personal/THCA_data_analysis/project/results/paper3_ici_public_data")
OUT_DIR  = Path("/home/seungho/personal/THCA_data_analysis/project/results/paper3_ici_data_registry")
OUT_DIR.mkdir(parents=True, exist_ok=True)

DATASETS = {
    # bulk / microarray (filename ends with _series_matrix.txt.gz)
    "GSE76039":   "GSE76039_series_matrix.txt.gz",
    "GSE65144":   "GSE65144_series_matrix.txt.gz",
    "GSE29265":   "GSE29265_series_matrix.txt.gz",
    "GSE33630":   "GSE33630_series_matrix.txt.gz",
    "GSE60542":   "GSE60542_series_matrix.txt.gz",
    "GSE151179":  "GSE151179_series_matrix.txt.gz",
    "GSE126698":  "GSE126698_series_matrix.txt.gz",
    "GSE78220":   "GSE78220_series_matrix.txt.gz",
    "GSE91061":   "GSE91061_series_matrix.txt.gz",
    # scRNA series matrices (named *_matrix.gz from earlier download)
    "GSE184362":  "GSE184362_matrix.gz",
    "GSE193581":  "GSE193581_matrix.gz",
    "GSE232237":  "GSE232237_matrix.gz",
    "GSE191288":  "GSE191288_matrix.gz",
    "GSE148673":  "GSE148673_matrix.gz",
}

def split_quoted(line: str) -> list[str]:
    parts = line.rstrip("\n").split("\t")
    return [p.strip().strip('"') for p in parts]

def parse_series_matrix(path: Path):
    series = {}
    samples = {}
    char_keys = []  # ordered keys observed in characteristics_ch1
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if not line.startswith("!"):
                continue
            cols = split_quoted(line)
            tag = cols[0]
            vals = cols[1:]
            if tag.startswith("!Series_"):
                key = tag[len("!Series_"):]
                # take first non-empty
                v = "; ".join([x for x in vals if x])
                series.setdefault(key, v)
                continue
            if tag.startswith("!Sample_"):
                key = tag[len("!Sample_"):]
                if key == "characteristics_ch1":
                    if "characteristics" not in samples:
                        samples["characteristics"] = []
                    samples["characteristics"].append(vals)
                else:
                    samples.setdefault(key, []).extend([])
                    samples[key] = vals
            if tag == "!series_matrix_table_begin":
                break
    return series, samples

def write_metadata(metadata_rows, qc_rows):
    smeta_path = OUT_DIR / "paper3_ici_public_sample_metadata.tsv"
    qc_path    = OUT_DIR / "paper3_ici_public_dataset_qc_summary.tsv"
    smeta_cols = [
        "dataset_id", "gsm_id", "sample_title", "source_name",
        "tumor_subtype", "tissue", "rai_status", "metastasis_site",
        "treatment", "response", "platform", "characteristics_full"
    ]
    with open(smeta_path, "w") as fh:
        fh.write("\t".join(smeta_cols) + "\n")
        for r in metadata_rows:
            fh.write("\t".join(r.get(c, "") for c in smeta_cols) + "\n")
    qc_cols = [
        "dataset_id", "n_samples", "platform_id", "platform_name_guess",
        "data_unit_guess", "histology_breakdown", "tumor_normal_breakdown",
        "rai_status_present", "treatment_present", "response_present",
        "immune_module_scoreable", "hla_typeable_arcasHLA",
        "neoantigen_calculable", "scrna_atlas_role", "notes"
    ]
    with open(qc_path, "w") as fh:
        fh.write("\t".join(qc_cols) + "\n")
        for r in qc_rows:
            fh.write("\t".join(str(r.get(c, "")) for c in qc_cols) + "\n")

def histology_breakdown(samples):
    chars = samples.get("characteristics", [])
    sub_idx = None
    for i, row in enumerate(chars):
        if any("subtype" in c.lower() or "histology" in c.lower() or "tumor type" in c.lower() for c in row if c):
            sub_idx = i
            break
    if sub_idx is None:
        # fallback: source_name
        sn = samples.get("source_name_ch1", [])
        from collections import Counter
        c = Counter(s for s in sn if s)
        return ";".join(f"{k}:{v}" for k, v in c.most_common())
    from collections import Counter
    vals = []
    for c in chars[sub_idx]:
        m = re.search(r":\s*(.+)$", c)
        vals.append(m.group(1).strip() if m else c.strip())
    cnt = Counter(vals)
    return ";".join(f"{k}:{v}" for k, v in cnt.most_common())

def detect_field(samples, keywords):
    chars = samples.get("characteristics", [])
    for row in chars:
        if any(any(kw.lower() in (c.lower()) for kw in keywords) for c in row if c):
            return True
    sn = " ".join(samples.get("source_name_ch1", []))
    return any(kw.lower() in sn.lower() for kw in keywords)

PLATFORM_GUESS = {
    "GPL570":   ("Affy_HG_U133_Plus_2", "microarray_intensity"),
    "GPL15456": ("Illumina_HiSeq_2500_RNAseq", "RNAseq_counts_or_FPKM_per_supp"),
    "GPL11154": ("Illumina_HiSeq_2000_RNAseq", "RNAseq_FPKM_or_counts"),
    "GPL9052":  ("Illumina_GA2_RNAseq",       "RNAseq_FPKM_or_counts"),
    "GPL23159": ("Affy_Clariom_D_Human",      "microarray_intensity"),
    "GPL24676": ("Illumina_NovaSeq_6000_scRNA","UMI_counts_via_10x"),
    "GPL20795": ("HiSeq_X_Ten_scRNA",         "UMI_counts_via_10x"),
}

def main():
    metadata_rows = []
    qc_rows = []
    for ds, fname in DATASETS.items():
        p = DATA_DIR / ds / fname
        if not p.exists():
            qc_rows.append({"dataset_id": ds, "notes": f"FILE_MISSING:{fname}"})
            continue
        try:
            series, samples = parse_series_matrix(p)
        except Exception as e:
            qc_rows.append({"dataset_id": ds, "notes": f"PARSE_ERROR:{e}"})
            continue
        gsm = samples.get("geo_accession", [])
        title = samples.get("title", [])
        sn    = samples.get("source_name_ch1", [])
        chars = samples.get("characteristics", [])
        plat  = series.get("platform_id", "")
        plat_name, unit = PLATFORM_GUESS.get(plat, (plat, "unknown"))

        # per-sample
        # transpose characteristics: list of rows with N values each → per-sample list of attributes
        n = len(gsm)
        per_sample_chars = []
        for j in range(n):
            row_attrs = []
            for row in chars:
                if j < len(row) and row[j]:
                    row_attrs.append(row[j])
            per_sample_chars.append(row_attrs)

        for j in range(n):
            attrs = per_sample_chars[j]
            attr_str = " | ".join(attrs)
            tumor_subtype = ""
            tissue = ""
            rai_status = ""
            metastasis = ""
            treatment = ""
            response = ""
            for a in attrs:
                low = a.lower()
                m = re.match(r"\s*([^:]+)\s*:\s*(.+)$", a)
                if not m:
                    continue
                k, v = m.group(1).strip().lower(), m.group(2).strip()
                if "subtype" in k or "tumor type" in k or "diagnosis" in k or "histology" in k:
                    tumor_subtype = v
                elif "tissue" in k:
                    tissue = v
                elif "iodine" in k or "rai" in k or "radioiodine" in k or "radioactive" in k or "refractory" in k or "avid" in k:
                    rai_status = v
                elif "metastas" in k or "site" in k:
                    metastasis = v
                elif "treatment" in k or "therapy" in k:
                    treatment = v
                elif "response" in k or "recist" in k or "best response" in k:
                    response = v
            metadata_rows.append({
                "dataset_id": ds,
                "gsm_id": gsm[j] if j < len(gsm) else "",
                "sample_title": title[j] if j < len(title) else "",
                "source_name": sn[j] if j < len(sn) else "",
                "tumor_subtype": tumor_subtype,
                "tissue": tissue,
                "rai_status": rai_status,
                "metastasis_site": metastasis,
                "treatment": treatment,
                "response": response,
                "platform": plat,
                "characteristics_full": attr_str,
            })

        # qc row
        hist = histology_breakdown(samples)
        from collections import Counter
        sn_counter = Counter(sn)
        tn = ";".join(f"{k}:{v}" for k, v in sn_counter.most_common())
        rai_present = detect_field(samples, ["iodine","rai","refractory","avid"])
        treat_present = detect_field(samples, ["treatment","therapy","pembro","nivo","atezo","ipi","chemo","tki","lenvatinib"])
        resp_present = detect_field(samples, ["response","recist","cr","pr","sd","pd","ORR"])
        # role guesses
        is_microarray = "microarray" in unit
        is_scrna = "scRNA" in plat_name or unit.startswith("UMI")
        immune_score = "yes_with_module_intersection" if not is_microarray else "yes_with_microarray_caveat_50pct_gene_coverage_required"
        hla = "yes_arcasHLA_if_FASTQ_available" if not is_microarray else "no_microarray"
        if is_scrna: hla = "limited_via_scHLA_or_arcasHLA_on_pseudobulk"
        neo = "yes_if_paired_WES_or_WGS_obtained_separately" if not is_microarray else "no_microarray"
        if is_scrna: neo = "no_no_paired_germline"
        scrna_role = "primary_scRNA_atlas_candidate" if is_scrna else "n_a"
        notes = []
        if ds == "GSE126698":
            notes.append("user_estimate_365_INCORRECT_actual_n=28_subtype_breakdown:ATC=10_FTC=6_NT=6_PTC=6_per_series_matrix")
        if ds == "GSE65144":
            notes.append("user_estimate_13_DIFFERENT_actual_n=25_per_series_matrix")
        if ds == "GSE60542":
            notes.append("user_unsure_actual_n=92_PTC_primary_and_nodal_metastases")
        if ds == "GSE193581":
            notes.append("processed_celltype_annotation_AND_cell_line_bulk_obtained_separately_in_suppl")
        if ds == "GSE184362":
            notes.append("11_patients_158577_cells_per_publication_23_GSM_libraries_RAW_tar_970MB_DEFER_full_to_track_B")

        qc_rows.append({
            "dataset_id": ds,
            "n_samples": n,
            "platform_id": plat,
            "platform_name_guess": plat_name,
            "data_unit_guess": unit,
            "histology_breakdown": hist,
            "tumor_normal_breakdown": tn,
            "rai_status_present": "yes" if rai_present else "no",
            "treatment_present": "yes" if treat_present else "no",
            "response_present": "yes" if resp_present else "no",
            "immune_module_scoreable": immune_score,
            "hla_typeable_arcasHLA": hla,
            "neoantigen_calculable": neo,
            "scrna_atlas_role": scrna_role,
            "notes": " | ".join(notes),
        })

    write_metadata(metadata_rows, qc_rows)
    print(f"sample_metadata rows: {len(metadata_rows)}")
    print(f"qc_summary rows:      {len(qc_rows)}")

if __name__ == "__main__":
    main()
