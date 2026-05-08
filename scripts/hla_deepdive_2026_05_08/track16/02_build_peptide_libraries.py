#!/usr/bin/env python3
"""
Track 16 Step 2: Sliding-window peptide libraries.
HLA-I = 9-mers, HLA-II = 15-mers (canonical core +/- registers averaged later).
Boundary: thyroid autoimmunity peptidomics — not cancer neoantigen prediction.
"""
from __future__ import annotations
import csv
import sys
from pathlib import Path

OUTDIR = Path("/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track16_peptidomics")
FASTA_DIR = OUTDIR / "fasta"
PANEL = OUTDIR / "tables" / "T01_antigen_panel.tsv"

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def parse_fasta_seq(path: Path) -> str:
    text = path.read_text().strip().splitlines()
    return "".join(l for l in text if not l.startswith(">"))


def sliding(seq: str, k: int):
    for i in range(len(seq) - k + 1):
        pep = seq[i:i + k]
        if all(a in VALID_AA for a in pep):
            yield i + 1, pep   # 1-based start position


def main() -> int:
    panel_rows = []
    with open(PANEL) as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for row in rdr:
            panel_rows.append(row)

    # HLA-I 9-mer library
    pepI_path = OUTDIR / "tables" / "T02_peptide_library_9mer.tsv"
    pepII_path = OUTDIR / "tables" / "T03_peptide_library_15mer.tsv"
    nI = nII = 0
    with open(pepI_path, "w") as fI, open(pepII_path, "w") as fII:
        fI.write("gene\taccession\tstart\tpeptide\n")
        fII.write("gene\taccession\tstart\tpeptide\n")
        for r in panel_rows:
            seq = parse_fasta_seq(FASTA_DIR / Path(r["fasta_path"]).name)
            for s, p in sliding(seq, 9):
                fI.write(f"{r['gene']}\t{r['accession']}\t{s}\t{p}\n")
                nI += 1
            for s, p in sliding(seq, 15):
                fII.write(f"{r['gene']}\t{r['accession']}\t{s}\t{p}\n")
                nII += 1
    print(f"wrote {pepI_path} (n={nI})")
    print(f"wrote {pepII_path} (n={nII})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
