#!/usr/bin/env python3
"""Wave 3 — prep inputs for off-the-shelf neoantigen scoring.

Bundle: project/results/p_neo_bayesian_2026_05_09/bundle.tsv
Filter to ITSNdb (ext_itsndb_main, ext_itsndb_val).
VENUS rows are full proteins (35-618 AA) → not 8-11mer → skipped for all
class-I tools (documented).
"""
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BUNDLE = ROOT.parent / "bundle.tsv"

df = pd.read_csv(BUNDLE, sep="\t")
print(f"bundle shape: {df.shape}")

# ITSNdb only (8-11mer, real HLA classI)
ext = df[df["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])].copy()
ext["len"] = ext["peptide"].str.len()
print(f"ITSNdb rows: {len(ext)}, peptide-length range: {ext['len'].min()}-{ext['len'].max()}")
print(ext[["split", "in_master"]].value_counts().sort_index())

# Sanity: drop any peptides with non-standard AA or wrong length (PRIME/BigMHC/MHCflurry expect 8-14)
AA = set("ACDEFGHIKLMNPQRSTVWY")
ext["clean"] = ext["peptide"].apply(lambda s: set(s) <= AA)
ext["len_ok"] = ext["len"].between(8, 14)
ext_clean = ext[ext["clean"] & ext["len_ok"]].copy()
print(f"after AA/len filter: {len(ext_clean)} (dropped {len(ext) - len(ext_clean)})")

# Save canonical scoring set
ext_clean[["peptide", "HLA_norm", "label", "source", "split", "in_master", "len"]].to_csv(
    ROOT / "scoring_set_itsndb.tsv", sep="\t", index=False
)

# === MHCflurry input ===
# Class1PresentationPredictor.predict accepts peptides + alleles directly via Python API.
# We'll feed it from Python. Just make a tidy CSV for tracking.
ext_clean[["peptide", "HLA_norm"]].to_csv(
    ROOT / "mhcflurry_input.csv", index=False
)

# === BigMHC input ===
# CSV with mhc,pep header (matches example1.csv format).
bmhc = ext_clean[["HLA_norm", "peptide", "label"]].copy()
bmhc.columns = ["mhc", "pep", "tgt"]
bmhc.to_csv(ROOT / "bigmhc_input.csv", index=False)

# === PRIME input ===
# PRIME wants ONE peptide per line (text), alleles on -a flag.
# Strategy: group peptides by their HLA, write one peptide-list per allele,
# then run PRIME per-allele. PRIME allele convention: A0101 (no HLA-, no *, no :).
def prime_allele(hla):
    # "HLA-A*02:01" -> "A0201"
    s = hla.replace("HLA-", "").replace("*", "").replace(":", "")
    return s

ext_clean["prime_allele"] = ext_clean["HLA_norm"].apply(prime_allele)
prime_dir = ROOT / "prime_inputs"
prime_dir.mkdir(exist_ok=True)
groups = ext_clean.groupby("prime_allele")
allele_summary = []
for allele, grp in groups:
    peps = grp["peptide"].tolist()
    out = prime_dir / f"peps_{allele}.txt"
    with open(out, "w") as fh:
        for p in peps:
            fh.write(p + "\n")
    allele_summary.append({"allele": allele, "n": len(peps), "file": str(out.name)})
pd.DataFrame(allele_summary).to_csv(ROOT / "prime_allele_index.tsv", sep="\t", index=False)
print(f"PRIME: wrote {len(allele_summary)} per-allele peptide files")

print("\nDONE prep")
print(f"scoring set rows: {len(ext_clean)}")
print(f"unique HLAs: {ext_clean['HLA_norm'].nunique()}")
