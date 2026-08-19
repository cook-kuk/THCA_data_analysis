"""
PDAC vaccine target prioritization.

Combines:
  - mutation prevalence (TCGA-PAAD frequency)
  - NeoQ score from script 04
  - tumor-vs-normal differential (proxy via mutation existence in tumor + GTEx
    normal pancreas reference for self-antigens; here we use a small curated
    self-antigen list)
  - HLA allele frequency overlap (Korean / pan-Asian / European)
  - public-vs-private flag (mKRAS = public ⇒ off-the-shelf candidate;
    private = personalized neoantigen)

Output:
  results/vaccine_targets.tsv    — ranked antigen panel
  processed/vaccine_priority_summary.json
"""

import json
from pathlib import Path
import pandas as pd
import numpy as np

RAW = Path("/data/pdac_poc/raw")
RES = Path("/data/pdac_poc/results")
OUT = Path("/data/pdac_poc/processed")

# Public neoantigen catalog - PDAC dominant + companion targets
# HLA frequencies from AFND (Korean / Pan-Asian / European pooled estimates)
PUBLIC_ANTIGENS = [
    {"target": "KRAS G12D", "class": "public_neoantigen", "vaccine_format": "mRNA / long-peptide / amphiphile",
     "hla_examples": ["A*11:01", "A*03:01", "B*07:02", "DRB1*04"],
     "afnd_korean_pct": 31.5,  # aggregated common alleles for these restrictions
     "afnd_panasian_pct": 33.0, "afnd_european_pct": 26.4,
     "evidence": "Pant 2024 ELI-002 amphiphile + Rojas 2023 BNT122 (subset)"},
    {"target": "KRAS G12V", "class": "public_neoantigen", "vaccine_format": "mRNA / peptide",
     "hla_examples": ["A*11:01", "B*07:02", "DRB1*0701"],
     "afnd_korean_pct": 28.1, "afnd_panasian_pct": 30.5, "afnd_european_pct": 23.0,
     "evidence": "Wang 2024 mKRAS peptide + Veatch JCI 2018 CD4 framework"},
    {"target": "KRAS G12R", "class": "public_neoantigen", "vaccine_format": "mRNA / peptide",
     "hla_examples": ["B*07:02", "DRB1*04:01"],
     "afnd_korean_pct": 8.9, "afnd_panasian_pct": 9.5, "afnd_european_pct": 16.6,
     "evidence": "PDAC-enriched; ELI-002 includes G12R cassette"},
    {"target": "KRAS G12C", "class": "public_neoantigen", "vaccine_format": "post-sotorasib vaccine",
     "hla_examples": ["A*02:01", "B*07:02"],
     "afnd_korean_pct": 22.0, "afnd_panasian_pct": 25.0, "afnd_european_pct": 38.0,
     "evidence": "rare in PDAC (~1-2%) but sotorasib synergy candidate"},
    {"target": "TP53 R175H", "class": "public_neoantigen", "vaccine_format": "mRNA / TCR-T",
     "hla_examples": ["A*02:01", "DRB1*15:01"],
     "afnd_korean_pct": 23.0, "afnd_panasian_pct": 25.0, "afnd_european_pct": 33.0,
     "evidence": "Lo Lab/NCI TCR clinical (Lowery JCI 2022)"},
    {"target": "TP53 R248Q", "class": "public_neoantigen", "vaccine_format": "TCR-T / vaccine cocktail",
     "hla_examples": ["A*02:01", "B*15:01"],
     "afnd_korean_pct": 18.5, "afnd_panasian_pct": 20.0, "afnd_european_pct": 27.0,
     "evidence": "Hot-spot mut, defined R-D class"},
    {"target": "Mesothelin", "class": "tumor_associated_self", "vaccine_format": "CAR-T / DC vaccine",
     "hla_examples": ["A*02:01", "DRB1*04:01"],
     "afnd_korean_pct": 22.0, "afnd_panasian_pct": 24.0, "afnd_european_pct": 33.0,
     "evidence": "CRS-207 listeria (Le 2015 JCO); benign pleura/peritoneum risk"},
    {"target": "WT1", "class": "tumor_associated_self", "vaccine_format": "DC peptide",
     "hla_examples": ["A*24:02", "A*02:01", "DRB1*04:05"],
     "afnd_korean_pct": 60.0, "afnd_panasian_pct": 55.0, "afnd_european_pct": 25.0,
     "evidence": "Higashihara 2014 PDAC trial; A*24:02-rich Asian populations"},
    {"target": "MUC1", "class": "tumor_associated_self", "vaccine_format": "MUC1-VNTR vaccine",
     "hla_examples": ["A*02:01", "DRB1*04"],
     "afnd_korean_pct": 23.0, "afnd_panasian_pct": 24.0, "afnd_european_pct": 32.0,
     "evidence": "Tecemotide / L-BLP25 historical; PDAC immunopeptidome support"},
    {"target": "Survivin (BIRC5)", "class": "tumor_associated_self", "vaccine_format": "peptide",
     "hla_examples": ["A*24:02", "A*02:01"],
     "afnd_korean_pct": 60.0, "afnd_panasian_pct": 55.0, "afnd_european_pct": 33.0,
     "evidence": "Validated CTL epitopes; A*24:02 Korean enrichment"},
    {"target": "Personalized neoantigen panel", "class": "private_neoantigen",
     "vaccine_format": "personalized mRNA (BNT122 class)",
     "hla_examples": ["per-patient"], "afnd_korean_pct": None,
     "afnd_panasian_pct": None, "afnd_european_pct": None,
     "evidence": "Rojas 2023 (PDAC) + Sethna 2025 long-FU"},
]


def main():
    # mutation prevalence from tcga-paad
    mut = pd.read_csv(RAW / "mutations.tsv", sep="\t")
    n_samples = pd.read_csv(RAW / "clinical.tsv", sep="\t")["sampleId"].nunique()

    def prev(gene, regex):
        sub = mut[(mut["hugo"] == gene) &
                  mut["proteinChange"].fillna("").str.match(regex, na=False)]
        return sub["sampleId"].nunique() / n_samples * 100

    prevalence = {
        "KRAS G12D": prev("KRAS", r"^G12D$"),
        "KRAS G12V": prev("KRAS", r"^G12V$"),
        "KRAS G12R": prev("KRAS", r"^G12R$"),
        "KRAS G12C": prev("KRAS", r"^G12C$"),
        "TP53 R175H": prev("TP53", r"^R175H$"),
        "TP53 R248Q": prev("TP53", r"^R248Q$"),
        "TP53 R273H": prev("TP53", r"^R273H$"),
    }

    # neoQ
    neo_summary = json.load(open(OUT / "neoantigen_quality_summary.json"))
    neoq_by_gene = neo_summary["top_neo_genes"]

    rows = []
    for ant in PUBLIC_ANTIGENS:
        target = ant["target"]
        gene = target.split()[0] if target.split() else ""
        prev_pct = prevalence.get(target, None)
        neoq = neoq_by_gene.get(gene)

        # priority score: 0.5 * prevalence/100 + 0.3 * korean_HLA/100 + 0.2 * neoQ
        prev_term = (prev_pct or 0) / 100.0
        kor_term = (ant.get("afnd_korean_pct") or 0) / 100.0
        neoq_term = (neoq or 0)
        priority = 0.5 * prev_term + 0.3 * kor_term + 0.2 * neoq_term

        rows.append({
            "target": target,
            "class": ant["class"],
            "vaccine_format": ant["vaccine_format"],
            "tcga_prevalence_pct": round(prev_pct, 2) if prev_pct is not None else None,
            "afnd_korean_pct": ant["afnd_korean_pct"],
            "afnd_panasian_pct": ant["afnd_panasian_pct"],
            "afnd_european_pct": ant["afnd_european_pct"],
            "neoQ_top": round(neoq, 4) if neoq is not None else None,
            "priority_score": round(priority, 4),
            "hla_examples": ",".join(ant["hla_examples"]),
            "evidence": ant["evidence"],
        })

    df = pd.DataFrame(rows).sort_values("priority_score", ascending=False)
    df.to_csv(RES / "vaccine_targets.tsv", sep="\t", index=False)

    # Korean coverage simulation: probability that a Korean PDAC patient carries
    # at least one targetable allele × mutation pair
    kor_covers = []
    public_targets = [r for r in rows if r["class"] == "public_neoantigen"
                      and r["tcga_prevalence_pct"]]
    for r in public_targets:
        p_carry_mut = (r["tcga_prevalence_pct"] or 0) / 100.0
        p_hla = (r["afnd_korean_pct"] or 0) / 100.0
        kor_covers.append((1 - p_carry_mut * p_hla))
    kor_combined = 1 - np.prod(kor_covers) if kor_covers else 0
    eu_covers = []
    for r in public_targets:
        p_carry_mut = (r["tcga_prevalence_pct"] or 0) / 100.0
        p_hla = (r["afnd_european_pct"] or 0) / 100.0
        eu_covers.append((1 - p_carry_mut * p_hla))
    eu_combined = 1 - np.prod(eu_covers) if eu_covers else 0

    summary = {
        "n_targets_ranked": int(len(df)),
        "tcga_prevalence_pct": prevalence,
        "korean_off_the_shelf_population_coverage_pct":
            round(kor_combined * 100, 2),
        "european_off_the_shelf_population_coverage_pct":
            round(eu_combined * 100, 2),
        "top5_priority": df.head(5)[["target", "priority_score",
                                      "tcga_prevalence_pct"]].to_dict(orient="records"),
    }
    with open(OUT / "vaccine_priority_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[06] {len(df)} targets ranked")
    for r in summary["top5_priority"]:
        print(f"  {r['target']:30s}  prev={r['tcga_prevalence_pct']}%  pri={r['priority_score']:.3f}")
    print(f"[06] Korean off-the-shelf pop coverage: {summary['korean_off_the_shelf_population_coverage_pct']}%")
    print(f"[06] European off-the-shelf pop coverage: {summary['european_off_the_shelf_population_coverage_pct']}%")


if __name__ == "__main__":
    main()
