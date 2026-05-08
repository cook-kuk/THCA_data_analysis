#!/usr/bin/env python3
"""Repair phenotype labels by concatenating ALL Sample_characteristics_ch1
rows from each cached GEO series matrix, then re-classify."""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
GEO = ROOT / "geo_meta"
RAW = GEO / "_raw"

DATASETS = ["GSE65144", "GSE60542", "GSE82208", "GSE76039"]


def parse_full_meta(path: Path) -> pd.DataFrame:
    """Return sample × characteristic dataframe by concatenating ALL
    Sample_characteristics_ch1 lines and Sample_source_name_ch1."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        lines = []
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            lines.append(line.rstrip("\n"))

    sample_ids: list[str] = []
    rows_by_key: dict[str, list[list[str]]] = {}
    for ln in lines:
        if "\t" not in ln:
            continue
        head, *vals = ln.split("\t")
        head = head.lstrip("!").strip()
        vals = [v.strip().strip('"') for v in vals]
        if head == "Sample_geo_accession":
            sample_ids = vals
        if head in ("Sample_characteristics_ch1", "Sample_source_name_ch1",
                    "Sample_title", "Sample_description"):
            rows_by_key.setdefault(head, []).append(vals)

    if not sample_ids:
        return pd.DataFrame()

    # Concatenate all characteristic rows per sample
    by_sample: dict[str, list[str]] = {sid: [] for sid in sample_ids}
    for key, row_list in rows_by_key.items():
        for row in row_list:
            for sid, val in zip(sample_ids, row):
                if val and val.lower() not in ("--", "n/a", "na", "null"):
                    by_sample[sid].append(f"[{key}] {val}")

    return pd.DataFrame(
        {"phenotype_full": ["  ||  ".join(by_sample[s]) for s in sample_ids]},
        index=sample_ids,
    )


def classify_phenotype(text: str | float | None) -> str:
    if text is None or (isinstance(text, float) and np.isnan(text)):
        return "unknown"
    s = str(text).lower()
    # Order matters: ATC > PDTC > PTC > FTC/follicular > normal
    if any(k in s for k in ("anaplastic", " atc", "atc ", "atc:", "atc,",
                            "atc/", "undifferentiated thyroid")):
        return "ATC"
    if any(k in s for k in ("poorly differentiated", "pdtc", "pdc:",
                            "poorly diff", "poorly-differentiated")):
        return "PDTC"
    if any(k in s for k in ("papillary", " ptc", "ptc ", "ptc:", "ptc,", "ptc/")):
        return "PTC"
    if any(k in s for k in ("follicular thyroid", "ftc:", "ftc ", " ftc",
                            "follicular adenoma", "follicular cancer",
                            "follicular carcinoma", "fa:", "fc:")):
        return "FTC"
    if any(k in s for k in ("normal", "control", "non-tumor", "non tumor",
                            "non-cancer", "non cancer", "matched normal",
                            "healthy", "benign")):
        return "normal"
    return "unknown"


def main() -> None:
    for gse in DATASETS:
        matrix = RAW / f"{gse}_series_matrix.txt.gz"
        scores = GEO / f"scores_{gse}.tsv"
        if not (matrix.exists() and scores.exists()):
            print(f"skip {gse} (missing matrix or scores)")
            continue
        meta = parse_full_meta(matrix)
        df = pd.read_csv(scores, sep="\t", index_col=0)
        df = df.loc[~df.index.duplicated(keep="last")]
        df = df.join(meta, how="left")
        df["phenotype"] = df["phenotype_full"].map(classify_phenotype)
        df.to_csv(scores, sep="\t")
        counts = df["phenotype"].value_counts()
        print(f"{gse}: n={len(df)}; phenotype counts = {dict(counts)}")


if __name__ == "__main__":
    main()
