#!/usr/bin/env python3
"""Build a probe → gene_symbol map for GPL23159 (Clariom S) from the SOFT family file.

Output: data/raw/GSE151179/GPL23159_probe2gene.tsv  (columns: ID, gene_symbol)
"""
import re
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "data" / "raw" / "GSE151179" / "GPL23159_annot_extracted.tsv"
OUT = ROOT / "data" / "raw" / "GSE151179" / "GPL23159_probe2gene.tsv"


# Patterns ordered by preference. Each must match the gene symbol in group 1.
PATTERNS = [
    # RefSeq with optional transcript variant: "...(SYMBOL), [transcript variant X,] mRNA"
    re.compile(r"\(([A-Z][A-Z0-9\-]{0,18})\)\s*,\s*(?:transcript variant[^,]*,\s*)?mRNA"),
    # ENSEMBL: "ENSEMBL // SYMBOL [gene_biotype:..."
    re.compile(r"ENSEMBL\s*//\s*([A-Z][A-Z0-9\-]{0,18})\s*\["),
    # ccdsGene: "ccdsGene // <description> [Source:HGNC Symbol;Acc:HGNC:XXXXX]" — fallback only
    re.compile(r"Source:HGNC Symbol;Acc:HGNC:\d+[^/]*?//[^/]*?\(([A-Z][A-Z0-9\-]{0,18})\)"),
]


def extract_symbol(text):
    if not isinstance(text, str):
        return None
    for pat in PATTERNS:
        m = pat.search(text)
        if m:
            return m.group(1)
    return None


def main():
    df = pd.read_csv(SRC, sep="\t", header=0, dtype=str, low_memory=False)
    last_col = df.columns[-1]
    df["gene_symbol"] = df[last_col].apply(extract_symbol)
    mapped = df.dropna(subset=["gene_symbol"]).copy()
    print(f"probes with symbol: {len(mapped)} / {len(df)}")
    mapped[["ID", "gene_symbol"]].to_csv(OUT, sep="\t", index=False)
    print(f"wrote {OUT}")

    panel = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
    found = mapped[mapped["gene_symbol"].isin(panel)]
    print(f"panel-gene probes: {found['gene_symbol'].nunique()}/{len(panel)} unique")
    print(found.groupby("gene_symbol").size())


if __name__ == "__main__":
    main()
