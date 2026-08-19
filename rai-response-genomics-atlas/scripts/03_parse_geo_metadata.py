#!/usr/bin/env python3
"""03 — parse a GEO series matrix file and extract RAI-related metadata.

Usage:  python3 scripts/03_parse_geo_metadata.py <ACCESSION>

Outputs data/interim/<ACCESSION>_metadata.tsv with sample-level RAI labels.
"""
from __future__ import annotations
import gzip
import re
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"

# Search patterns for RAI / iodine / response in characteristics_ch1
RAI_TOKENS = [
    "iodine", "iodide", "i-131", "i131", "rai", "radioiodine", "radioactive",
    "avid", "refractory", "uptake", "remission", "persistence", "persistent",
    "response", "recurrence", "redifferentiation",
]


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def parse_series_matrix(path: Path) -> pd.DataFrame:
    """Read GEO series matrix; return a sample-level DataFrame of !Sample_* rows."""
    rows: dict[str, list[str]] = {}
    sample_ids: list[str] = []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            parts = line.split("\t")
            key = parts[0].lstrip("!")
            vals = [p.strip('"') for p in parts[1:]]
            rows.setdefault(key, []).append(vals)
    if not rows:
        raise ValueError(f"no !Sample_ rows parsed from {path}")
    n_samples = len(next(iter(rows.values()))[0])
    out: dict[str, list[str]] = {}
    for key, list_of_lists in rows.items():
        # GEO can repeat 'characteristics_ch1' multiple times; concatenate
        if len(list_of_lists) == 1:
            out[key] = list_of_lists[0]
        else:
            joined = []
            for i in range(n_samples):
                joined.append(" | ".join(ll[i] for ll in list_of_lists if i < len(ll)))
            out[f"{key}_joined"] = joined
    df = pd.DataFrame(out)
    if "Sample_geo_accession" in df.columns:
        df = df.set_index("Sample_geo_accession")
    return df


def extract_rai_labels(df: pd.DataFrame) -> pd.DataFrame:
    """Look for RAI / iodine / response tokens in any characteristics field.

    Also parses GEO 'key: value' pairs commonly used inside Sample_characteristics_ch1.
    """
    char_cols = [c for c in df.columns if "characteristics" in c.lower() or c.endswith("_joined")]
    label_blob = df[char_cols].astype(str).agg(" | ".join, axis=1) if char_cols else pd.Series("", index=df.index)
    df = df.copy()
    df["rai_label_blob"] = label_blob
    for tok in RAI_TOKENS:
        df[f"has_{tok}"] = label_blob.str.contains(tok, case=False, regex=False).astype(int)

    def kv_extract(s: str, key_regex: str) -> str:
        m = re.search(key_regex + r"\s*:\s*([^|]+?)(?:\s*\||$)", s, flags=re.I)
        return m.group(1).strip() if m else ""

    df["rai_response"]      = label_blob.map(lambda s: kv_extract(s, r"(?:patient[ _]?)?rai[ _]?respon[cs]e")).str.lower()
    df["rai_uptake_met"]    = label_blob.map(lambda s: kv_extract(s, r"rai[ _]?uptake(?:[ _]?at[ _]?the[ _]?metastatic[ _]?site)?")).str.lower()
    # disease status: field may be "disease", "disease status", "disease remission"
    df["disease_status"]    = label_blob.map(lambda s: kv_extract(s, r"disease(?:[ _]?(?:status|remission))?")).str.lower()
    # sample type: prefer "tissue type" / "tissue types" over generic "tissue"
    df["sample_type"]       = label_blob.map(lambda s: kv_extract(s, r"tissue[ _]?types?")).str.lower()
    # if the specific key wasn't found (e.g. only generic 'tissue: ...'), fall back
    fallback_tissue = label_blob.map(lambda s: kv_extract(s, r"tissue")).str.lower()
    df.loc[df["sample_type"].eq(""), "sample_type"] = fallback_tissue[df["sample_type"].eq("")]
    df["collection_timing"] = label_blob.map(lambda s: kv_extract(s, r"collection[ _]?before/after[ _]?rai")).str.lower()
    df["lesion_class"]      = label_blob.map(lambda s: kv_extract(s, r"lesion[ _]?class")).str.lower()
    df["lesion_driver"]     = label_blob.map(lambda s: kv_extract(s, r"lesion[ _]?by[ _]?ptc[-_]?ma")).str.lower()
    df["patient_id"]        = label_blob.map(lambda s: kv_extract(s, r"patient[ _]?id"))
    df["histology_variant"] = label_blob.map(lambda s: kv_extract(s, r"histolog(?:ical|y)[ _]?variant"))
    df["tumor_purity"]      = label_blob.map(lambda s: kv_extract(s, r"tumor[ _]?purity[ _]?class"))

    df["rai_avid_refractory"] = df["rai_response"].str.lower()
    df["sample_type_guess"] = df["sample_type"]
    return df


def main() -> None:
    if len(sys.argv) < 2:
        print("usage: 03_parse_geo_metadata.py <ACCESSION>", file=sys.stderr); sys.exit(1)
    acc = sys.argv[1].strip()
    src_dir = RAW / acc
    if not src_dir.exists():
        print(f"missing {src_dir}", file=sys.stderr); sys.exit(2)
    matrices = sorted(src_dir.glob(f"{acc}*series_matrix.txt*"))
    if not matrices:
        print(f"no series_matrix file in {src_dir}", file=sys.stderr); sys.exit(3)
    INTERIM.mkdir(parents=True, exist_ok=True)
    all_meta = []
    for m in matrices:
        print(f"# parsing {m.name}")
        df = parse_series_matrix(m)
        df = extract_rai_labels(df)
        df["_source"] = m.name
        all_meta.append(df)
    meta = pd.concat(all_meta, axis=0)
    out = INTERIM / f"{acc}_metadata.tsv"
    meta.to_csv(out, sep="\t")
    print(f"# wrote {out}  ({meta.shape[0]} samples × {meta.shape[1]} cols)")

    # Summary of label tokens
    token_cols = [c for c in meta.columns if c.startswith("has_")]
    summary = meta[token_cols].sum().sort_values(ascending=False)
    print("\n# RAI-related token counts:")
    print(summary[summary > 0].to_string())


if __name__ == "__main__":
    main()
