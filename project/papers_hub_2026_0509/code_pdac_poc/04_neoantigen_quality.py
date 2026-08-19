"""
Balachandran 2017 (Nature) neoantigen quality framework — toy/PoC implementation
using the public PDAC mutation profile from cBioPortal.

For each missense mutation in TCGA-PAAD samples, we:
  1) generate a window-of-9 mutant peptide centered on the variant aa
  2) compute D = sequence dissimilarity to the corresponding self peptide
     (1 - identity over the 9-mer; we have no full proteome lookup, so we use
     the variant alone as a stand-in for the differential — the formal pipeline
     would slide windows; this is a clearly-marked toy)
  3) compute R = MAX cross-reactivity to a curated IEDB-pathogen anchor set
     using a substitution-matrix alignment proxy (BLOSUM62 sum)
  4) NeoQ = R × D, sample-level top-NeoQ ranking + class summaries

Outputs:
  results/neoantigen_quality_table.tsv  — per-mutation
  results/neoantigen_quality_summary.json — per-sample aggregates

This is explicitly toy-grade for the PoC dashboard. Production-grade requires
NetMHCpan binding + TCRdist3 cross-reactivity + full proteome BLAST.
"""

import json
import math
from pathlib import Path
import numpy as np
import pandas as pd

RAW = Path("/data/pdac_poc/raw")
OUT = Path("/data/pdac_poc/processed")
RES = Path("/data/pdac_poc/results")
RES.mkdir(parents=True, exist_ok=True)

# tiny BLOSUM62 lookup (subset)
AA = "ARNDCQEGHILKMFPSTWYV"
B62 = np.array([
    [4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0],
    [-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3],
    [-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3],
    [-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3],
    [0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
    [-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2],
    [-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2],
    [0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3],
    [-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3],
    [-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3],
    [-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1],
    [-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2],
    [-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1],
    [-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1],
    [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2],
    [1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2],
    [0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0],
    [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3],
    [-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1],
    [0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4],
])
AAIX = {a: i for i, a in enumerate(AA)}


# Curated IEDB-pathogen anchor 9-mers (representative, public): EBV BMLF1
# GLCTLVAML, EBV BRLF1 YVLDHLIVV, CMV pp65 NLVPMVATV, Influenza M1 GILGFVFTL,
# HIV Nef KRWIILGLNK, HCV NS3 KLVALGINAV, HPV16 E7 YMLDLQPETT (for HLA-B*07
# variants) etc. We use a curated 16-peptide pool.
PATHOGEN_PEPTIDES = [
    "GLCTLVAML",  # EBV BMLF1 (HLA-A*02:01)
    "YVLDHLIVV",  # EBV BRLF1
    "NLVPMVATV",  # CMV pp65
    "GILGFVFTL",  # Influenza M1
    "KRWIILGLNK", # HIV Nef
    "KLVALGINAV", # HCV NS3
    "YMLDLQPETT", # HPV16 E7 (10-mer, trimmed)
    "RPHERNGFTV", # HCV
    "FLAFLDKEY",  # IAV
    "TPRVTGGGAM", # M.tuberculosis
    "CTELKLSDY",  # adenovirus E1A
    "KAFSPEVIPMF",# Mart-1 modified
    "ELAGIGILTV", # Mart-1
    "SLYNTVATL",  # HIV Gag
    "ITDQVPFSV",  # HCMV pp65
    "RAKFKQLL",   # EBV BZLF1
]


def blosum_score(a, b):
    """Sum of BLOSUM62 over aligned residues (no gaps)."""
    if len(a) != len(b):
        L = min(len(a), len(b))
        a, b = a[:L], b[:L]
    s = 0
    for x, y in zip(a, b):
        if x in AAIX and y in AAIX:
            s += B62[AAIX[x], AAIX[y]]
    return s


def cross_reactivity_R(neo_pep, k=4.86, a=26.0):
    """Balachandran 2017 R = sum exp(k*(a-s_iedb)) sigmoid, here we use the
    max-similarity heuristic to keep PoC tractable."""
    if not neo_pep:
        return 0.0
    scores = []
    for pp in PATHOGEN_PEPTIDES:
        s = blosum_score(neo_pep, pp[:len(neo_pep)])
        scores.append(s)
    s_max = max(scores)
    # logistic on alignment quality, normalized
    return 1.0 / (1.0 + math.exp(-(s_max - 10.0) / 4.0))


def dissimilarity_D(mutant, wildtype):
    """1 - fraction identity. Larger ⇒ more divergent from self."""
    L = min(len(mutant), len(wildtype))
    if L == 0:
        return 0.0
    matches = sum(1 for a, b in zip(mutant[:L], wildtype[:L]) if a == b)
    return 1.0 - matches / L


def parse_protein_change(pc):
    """e.g., 'G12D' or 'p.G12D' -> ('G', 12, 'D'). Returns None for non-missense."""
    if not isinstance(pc, str):
        return None
    pc = pc.strip()
    if pc.startswith("p."):
        pc = pc[2:]
    if len(pc) < 3:
        return None
    # missense form only: <AA><digits><AA>
    if pc[0] not in AA or pc[-1] not in AA:
        return None
    middle = pc[1:-1]
    if not middle.isdigit():
        return None
    return pc[0], int(middle), pc[-1]


# Tiny pseudo-context map for KRAS, TP53 hotspots — surrounding residues
# (real pipeline pulls full RefSeq protein). PoC only.
HOTSPOT_CONTEXT = {
    ("KRAS", 12): "MTEYKLVVVGAGGVGKSALTIQLIQ",  # 12 ~ position 13
    ("KRAS", 13): "MTEYKLVVVGAGGVGKSALTIQLIQ",
    ("KRAS", 61): "AGGVGKSALTIQLIQNHFVDEYDPT",
    ("KRAS", 117): "DTAGQEEYSAMRDQYMRTGEGFLCV",
    ("KRAS", 146): "EDAFYTLVREIRQHKLRKLNPPDESG",
    ("TP53", 175): "VRVCACPGRDRRTEEENLRKKGEPHHELP",
    ("TP53", 245): "MGGMNRRPILTIITLEDSSGNLLGRNSF",
    ("TP53", 248): "RPILTIITLEDSSGNLLGRNSFEVRVC",
    ("TP53", 249): "PILTIITLEDSSGNLLGRNSFEVRVCA",
    ("TP53", 273): "GTRVRAMAIYKQSQHMTEVVRRCPHHE",
    ("TP53", 282): "AIYKQSQHMTEVVRRCPHHERCSDSDG",
    ("SMAD4", 361): "GSPNARRFSDQLNVTPLNVMASLPNYG",
    ("CDKN2A", 80): "RRAGSSMEPSADWLATAAARGRVEEVR",
}


def build_window(gene, pos, mut_aa, win=4):
    ctx = HOTSPOT_CONTEXT.get((gene, pos))
    if not ctx:
        return None, None
    # Build a 9-mer centered on the hotspot. We assume ctx is the local stretch
    # starting at pos-4 (PoC convention).
    if len(ctx) < 9:
        return None, None
    wt = ctx[:9]
    mut = wt[:4] + mut_aa + wt[5:]
    return wt, mut


def main():
    mut = pd.read_csv(RAW / "mutations.tsv", sep="\t")
    mut = mut[mut["mutationType"].fillna("").str.contains("Missense")]
    print(f"[04] missense rows: {len(mut)}")

    rows = []
    for _, r in mut.iterrows():
        parsed = parse_protein_change(r["proteinChange"])
        if not parsed:
            continue
        wt_aa, pos, mu_aa = parsed
        gene = r.get("hugo")
        wt9, mu9 = build_window(gene, pos, mu_aa)
        if not wt9:
            continue
        D = dissimilarity_D(mu9, wt9)
        R = cross_reactivity_R(mu9)
        NeoQ = R * D
        rows.append({
            "sampleId": r["sampleId"],
            "gene": gene,
            "proteinChange": r["proteinChange"],
            "wt_peptide": wt9, "mu_peptide": mu9,
            "D": D, "R": R, "NeoQ": NeoQ,
        })

    df = pd.DataFrame(rows)
    df.to_csv(RES / "neoantigen_quality_table.tsv", sep="\t", index=False)
    print(f"[04] scored peptides: {len(df)}")

    # per-sample aggregate
    df_sorted = df.sort_values(["sampleId", "NeoQ"], ascending=[True, False])
    top_per_sample = df_sorted.groupby("sampleId").head(1).set_index("sampleId")
    agg = (df.groupby("sampleId")
           .agg(n_missense_hot=("gene", "size"),
                top_NeoQ=("NeoQ", "max"),
                mean_NeoQ=("NeoQ", "mean"))
           .reset_index())
    agg["top_gene"] = agg["sampleId"].map(top_per_sample["gene"])
    agg["top_change"] = agg["sampleId"].map(top_per_sample["proteinChange"])
    agg.to_csv(RES / "neoantigen_quality_per_sample.tsv", sep="\t", index=False)

    summary = {
        "n_samples_scored": int(agg["sampleId"].nunique()),
        "n_peptides_scored": int(len(df)),
        "kras_g12_count": int(((df["gene"] == "KRAS") &
                                df["proteinChange"]
                                  .str.match(r"^G12[A-Z]$", na=False)).sum()),
        "kras_g12_share_of_kras": float(
            ((df["gene"] == "KRAS") &
             df["proteinChange"].str.match(r"^G12[A-Z]$", na=False)).sum() /
            max(1, (df["gene"] == "KRAS").sum())),
        "tp53_R175H_count": int(((df["gene"] == "TP53") &
                                  (df["proteinChange"] == "R175H")).sum()),
        "tp53_R248_count": int(((df["gene"] == "TP53") &
                                 df["proteinChange"]
                                   .str.match(r"^R248[A-Z]$", na=False)).sum()),
        "tp53_R273_count": int(((df["gene"] == "TP53") &
                                 df["proteinChange"]
                                   .str.match(r"^R273[A-Z]$", na=False)).sum()),
        "top_neo_genes": (df.groupby("gene")["NeoQ"]
                          .max().sort_values(ascending=False).head(10).to_dict()),
        "median_top_NeoQ": float(agg["top_NeoQ"].median()) if len(agg) else None,
    }
    with open(OUT / "neoantigen_quality_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[04] samples scored: {summary['n_samples_scored']}; "
          f"median top NeoQ: {summary['median_top_NeoQ']}")


if __name__ == "__main__":
    main()
