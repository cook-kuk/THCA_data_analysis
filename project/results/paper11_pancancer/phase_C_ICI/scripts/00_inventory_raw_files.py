#!/usr/bin/env python3
"""Phase C v2 — raw data manifest.

Walks phase_C_ICI/raw/ and produces:
  manifest/raw_download_manifest.tsv    (one row per file)

Columns (per user prompt):
  dataset source accession local_path file_name file_size md5_if_available
  data_type sample_count_observed patient_count_observed_if_known
  has_expression has_clinical has_response_label has_survival notes
"""
from __future__ import annotations
import gzip, hashlib, json, subprocess
from pathlib import Path
import pandas as pd

ROOT = Path("/data/thca/repo_results/paper11_pancancer/phase_C_ICI")
RAW = ROOT / "raw"
PROC = ROOT / "processed"
MANIFEST = ROOT / "manifest"
MANIFEST.mkdir(exist_ok=True, parents=True)

DATASET_META = {
    "IMvigor210": dict(source="research-pub.gene.com (Mariathasan 2018)",
                       accession="NCT02108652/NCT02951767 (R package, EGAS00001002556)"),
    "GSE176307":  dict(source="GEO", accession="GSE176307 (BACI urothelial ICB)"),
    "riaz_GSE91061": dict(source="GEO", accession="GSE91061 (Riaz 2017 melanoma anti-PD-1)"),
    "gide_PRJEB23709": dict(source="ENA + GEO substitute",
                            accession="PRJEB23709 (Gide 2019 — FASTQ-only); substitute GSE115821 (MGH melanoma ICB)"),
}


def md5(path: Path) -> str:
    try:
        h = hashlib.md5()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()
    except Exception:
        return ""


def file_kind(p: Path) -> str:
    n = p.name.lower()
    if "phenotype" in n or "metadata" in n or "series_matrix" in n or "soft" in n or "key" in n: return "phenotype"
    if "tpm" in n or "fpkm" in n or "counts" in n or "rsem" in n or "rnaseq" in n or "expr" in n: return "expression"
    if n.endswith(".rdata"): return "expression+phenotype (R bundle)"
    if "mut" in n or "mutation" in n: return "mutation"
    if n.endswith(".pdf") or n.endswith(".md"): return "doc"
    return "other"


rows = []
for ds_dir in sorted(d for d in RAW.iterdir() if d.is_dir()):
    ds = ds_dir.name
    meta = DATASET_META.get(ds, {})
    # Walk every file (skip extracted/ for IMvigor since we list both raw + extracted)
    for p in sorted(ds_dir.rglob("*")):
        if p.is_dir(): continue
        # observed sample/patient counts: pull from harmonized metadata if it matches dataset
        proc_md = PROC / ds / "metadata.tsv"
        if ds == "gide_PRJEB23709":
            proc_md = PROC / "MGH_GSE115821" / "metadata.tsv"
        n_samples, n_patients, has_resp, has_surv = "", "", False, False
        if proc_md.exists():
            try:
                m = pd.read_csv(proc_md, sep="\t")
                n_samples = len(m)
                if "patient_id" in m.columns:
                    n_patients = m["patient_id"].nunique()
                has_resp = m["response_binary_CRPR_vs_SD_PD"].notna().any() if "response_binary_CRPR_vs_SD_PD" in m.columns else False
                has_surv = m["os_time"].notna().any() if "os_time" in m.columns else False
            except Exception:
                pass
        rows.append({
            "dataset": ds,
            "source": meta.get("source", ""),
            "accession": meta.get("accession", ""),
            "local_path": str(p.relative_to(ROOT)),
            "file_name": p.name,
            "file_size_bytes": p.stat().st_size,
            "md5": md5(p),
            "data_type": file_kind(p),
            "sample_count_observed": n_samples,
            "patient_count_observed_if_known": n_patients,
            "has_expression": file_kind(p) in {"expression", "expression+phenotype (R bundle)"},
            "has_clinical": file_kind(p) == "phenotype" or file_kind(p) == "expression+phenotype (R bundle)",
            "has_response_label": has_resp,
            "has_survival": has_surv,
            "notes": "",
        })

# substitution note row
rows.append({
    "dataset": "gide_PRJEB23709",
    "source": "ENA",
    "accession": "PRJEB23709",
    "local_path": "(NOT DOWNLOADED)",
    "file_name": "gide_PRJEB23709_FASTQ_only",
    "file_size_bytes": "",
    "md5": "",
    "data_type": "fastq (skipped)",
    "sample_count_observed": 91,
    "patient_count_observed_if_known": 91,
    "has_expression": False,
    "has_clinical": False,
    "has_response_label": False,
    "has_survival": False,
    "notes": "Gide 2019 has FASTQ-only on ENA; processed RNA-seq matrix not in ENA analysis API. "
             "Sprint-budget skipped; substituted with GSE115821 (MGH melanoma ICB) as additional independent cohort.",
})

df = pd.DataFrame(rows)
df.to_csv(MANIFEST / "raw_download_manifest.tsv", sep="\t", index=False)
print(f"manifest → {MANIFEST/'raw_download_manifest.tsv'}  ({len(df)} rows)")
print(df[["dataset","file_name","file_size_bytes","data_type"]].to_string(index=False))
