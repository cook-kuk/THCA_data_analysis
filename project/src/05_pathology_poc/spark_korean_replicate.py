#!/usr/bin/env python3
"""
Korean Lee 2024 (GSE213647 n=632) external replicate of CCI signature.

Downloads GSE213647_RAW.tar (240 MB), extracts 632 STAR ReadsPerGene.out.tab.gz,
builds expression matrix for CCI module genes only, then replicates:
  CAF / M1_M2 / TLS / CD8 × RAI_lineage Spearman

If CAF × RAI = -0.31 (TCGA) replicates in Korean independent cohort,
that is 3-cohort cross-platform replication = paper-grade.
"""
from __future__ import annotations
import gzip, io, os, sys, tarfile
from pathlib import Path
from urllib.request import urlretrieve

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent.parent.parent
RAW_DIR = ROOT / "project/data/raw/GSE213647"
RES = ROOT / "project/results/03_pathology_poc"
TAR_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213647/suppl/GSE213647_RAW.tar"
TAR_PATH = RAW_DIR / "GSE213647_RAW.tar"
EXTRACT_DIR = RAW_DIR / "extracted"

MODULES = {
    "CAF":   ["FAP", "ACTA2", "PDGFRA", "PDGFRB", "COL1A1", "COL3A1", "DCN", "POSTN"],
    "M1_M2": ["CD68", "CD163", "MRC1", "MARCO", "C1QA", "C1QB", "C1QC", "TYROBP"],
    "CD8":   ["CD8A", "CD8B", "GZMB", "GZMK", "PRF1", "IFNG"],
    "TLS":   ["CXCL13", "CCL19", "CCR7", "LTB", "LTA"],
    "RAI":   ["TPO", "DIO1", "TSHR", "PAX8", "TG", "FOXE1", "NKX2-1", "SLC5A5"],
}


def download_tar():
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    if not TAR_PATH.exists() or TAR_PATH.stat().st_size < 200_000_000:
        print(f"downloading {TAR_URL} → {TAR_PATH}")
        urlretrieve(TAR_URL, TAR_PATH)
    print(f"tar size: {TAR_PATH.stat().st_size/1e6:.1f} MB")


def extract_tar():
    EXTRACT_DIR.mkdir(parents=True, exist_ok=True)
    if len(list(EXTRACT_DIR.glob("*ReadsPerGene*"))) < 600:
        print("extracting tar...")
        with tarfile.open(TAR_PATH) as tf:
            tf.extractall(EXTRACT_DIR, filter="data")
    print(f"extracted: {len(list(EXTRACT_DIR.glob('*ReadsPerGene*')))} files")


def parse_star_counts(path):
    """STAR ReadsPerGene.out.tab columns: ENSG, unstranded, forward, reverse.
       First 4 lines are summary (N_unmapped/multimapping/no_feature/ambiguous).
       Return dict of {ENSG_clean: count_col_4_reverse_stranded}."""
    counts = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4: continue
            gene = parts[0]
            if gene.startswith("N_"): continue
            # GENCODE/Ensembl ID like ENSG00000142973.13
            ensg = gene.split(".")[0]
            try:
                counts[ensg] = int(parts[3])  # reverse-stranded col
            except ValueError:
                continue
    return counts


def ensg_to_symbol_map():
    """Use cached GENCODE mapping, or grab a small one."""
    cache = ROOT / "project/data/raw/ensg_to_symbol.tsv"
    if cache.exists():
        return pd.read_csv(cache, sep="\t").set_index("ensg")["symbol"].to_dict()
    # Fallback: grab from a small public mapping
    print("downloading ENSG↔symbol map (small subset)...")
    url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.metadata.HGNC.gz"
    import urllib.request, gzip as gz
    cache.parent.mkdir(parents=True, exist_ok=True)
    try:
        urllib.request.urlretrieve(url, str(cache) + ".raw.gz")
        rows = []
        with gz.open(str(cache) + ".raw.gz", "rt") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    rows.append((parts[0].split(".")[0], parts[1]))
        df = pd.DataFrame(rows, columns=["txid", "symbol"]).drop_duplicates()
        # Build a fallback ensg→symbol via a different metadata file
        # Use GTF instead
    except Exception as e:
        print(f"GENCODE map fail: {e}")
    return {}


def main():
    download_tar()
    extract_tar()

    files = sorted(EXTRACT_DIR.glob("*ReadsPerGene*"))
    print(f"parsing {len(files)} STAR count files...")
    all_genes = set()
    sample_counts = {}
    for fp in files:
        sid = fp.name.split("_")[0]  # GSM6590734
        sample_counts[sid] = parse_star_counts(fp)
        all_genes.update(sample_counts[sid].keys())
    print(f"  collected {len(sample_counts)} samples × {len(all_genes)} ENSG")

    # Build minimal matrix: only need ENSGs that map to our module genes
    # We need an ENSG→symbol mapping. Try mygene.info first
    target_symbols = list({g for v in MODULES.values() for g in v})
    print(f"\ntarget symbols: {target_symbols}")

    # Use REST API for a quick mapping (no local cache)
    import urllib.request, json
    sym_to_ensg = {}
    for sym in target_symbols:
        try:
            url = f"https://mygene.info/v3/query?q={sym}&species=human&fields=ensembl.gene&size=1"
            with urllib.request.urlopen(url, timeout=15) as r:
                d = json.loads(r.read())
                hit = d.get("hits", [{}])[0]
                ens = hit.get("ensembl", {})
                if isinstance(ens, list): ens = ens[0]
                if ens.get("gene"):
                    sym_to_ensg[sym] = ens["gene"]
        except Exception as e:
            print(f"  {sym}: lookup fail {e}")
    print(f"sym→ENSG: {len(sym_to_ensg)} / {len(target_symbols)} mapped")
    print(sym_to_ensg)

    # Build matrix (gene-symbol × sample)
    rows = []
    for sym, ensg in sym_to_ensg.items():
        row = {"symbol": sym}
        for sid, counts in sample_counts.items():
            row[sid] = counts.get(ensg, 0)
        rows.append(row)
    df = pd.DataFrame(rows).set_index("symbol")
    print(f"matrix: {df.shape}")
    df.to_csv(RES / "spark_korean_GSE213647_module_counts.tsv", sep="\t")

    # log1p + per-gene z-score → module mean
    df_log = np.log1p(df.astype(float))
    z = df_log.subtract(df_log.mean(axis=1), axis=0).div(df_log.std(axis=1) + 1e-6, axis=0)
    mods = {}
    for name, syms in MODULES.items():
        present = [s for s in syms if s in z.index]
        if len(present) < 2: continue
        mods[name] = z.loc[present].mean(axis=0)

    rows = []
    for n1 in ["CAF", "M1_M2", "TLS", "CD8"]:
        if n1 not in mods or "RAI" not in mods: continue
        common = mods[n1].index.intersection(mods["RAI"].index)
        if len(common) < 30: continue
        r = float(spearmanr(mods[n1].loc[common], mods["RAI"].loc[common]).statistic)
        rows.append({"cohort": "Korean_GSE213647", "module": n1, "n": len(common), "r_vs_RAI": r})
        print(f"  Korean (n={len(common)}): {n1} × RAI Spearman = {r:+.3f}")

    pd.DataFrame(rows).to_csv(RES / "spark_korean_GSE213647_replicate.tsv", sep="\t", index=False)
    print(f"\nwrote {(RES / 'spark_korean_GSE213647_replicate.tsv').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
