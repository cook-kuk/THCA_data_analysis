#!/usr/bin/env python3
"""
Track 2 / Step 01 — Parse arcasHLA genotype JSONs for the K2 PRJEB11591 cohort
into a wide sample x locus table at 2-field (4-digit) resolution.

Boundary statement (CRITICAL):
  This script characterizes the HLA frequency distribution in K2 as a Korean
  RNA-seq imputation snapshot. It does NOT condition on tumor/normal/disease
  state. It does NOT test cancer or autoimmune association. See
  project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md.
"""
import json
import os
import re
from pathlib import Path

import pandas as pd

ARCAS_DIR = Path("/data/thca/_repo_offload/arcasHLA")
OUT_DIR = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "hla_deepdive_2026_05_08/track2_k2_baseline/tables"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

LOCI = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
NONCLASSICAL = ["E", "F", "G"]  # arcasHLA emits these only if reference panel includes them


def to_2field(allele: str) -> str:
    """Convert e.g. 'A*11:303' or 'A*32:01:01' to 'A*11:303' / 'A*32:01' (2 fields).

    Rules:
      - Strip trailing :NN groups beyond the first two colon-separated fields.
      - Preserve null-suffix letters (N, L, S, Q) when present at the *last kept* field.
      - 'A*11:303' -> 'A*11:303' (already 2 fields).
    """
    if not isinstance(allele, str) or "*" not in allele:
        return allele
    locus, rest = allele.split("*", 1)
    parts = rest.split(":")
    if len(parts) <= 2:
        return f"{locus}*{':'.join(parts)}"
    # Keep first two fields. Preserve any expression-suffix letter on the second field.
    f2 = parts[1]
    # If field 3+ has a suffix letter at the end and field 2 is purely numeric, that suffix
    # logically belongs to the truncated allele — keep its expression letter on field 2.
    return f"{locus}*{parts[0]}:{f2}"


def sample_id_from_filename(fname: str) -> str:
    # arcasHLA filenames look like ERR1518628_1.genotype.json
    m = re.match(r"(ERR\d+|SRR\d+)(?:_\d+)?\.genotype\.json$", fname)
    if m:
        return m.group(1)
    return fname.replace(".genotype.json", "")


def main() -> None:
    rows = []
    files = sorted(ARCAS_DIR.glob("*.genotype.json"))
    for fp in files:
        try:
            with open(fp) as fh:
                data = json.load(fh)
        except Exception as exc:
            print(f"[WARN] could not parse {fp.name}: {exc}")
            continue

        run_acc = sample_id_from_filename(fp.name)
        row = {"run_accession": run_acc, "source_file": fp.name}

        for locus in LOCI + NONCLASSICAL:
            calls = data.get(locus, [])
            if isinstance(calls, list) and len(calls) >= 1:
                a1 = to_2field(calls[0]) if calls[0] else None
                a2 = to_2field(calls[1]) if len(calls) > 1 and calls[1] else None
            else:
                a1 = a2 = None
            row[f"{locus}_a1"] = a1
            row[f"{locus}_a2"] = a2
        rows.append(row)

    df = pd.DataFrame(rows).sort_values("run_accession").reset_index(drop=True)

    # Merge in K2 metadata (sample_alias, etc.) from the run sheet for traceability —
    # but we DO NOT use any disease/histology label downstream (boundary).
    runs_path = Path(
        "/home/seungho/personal/THCA_data_analysis/project/results/v17_korean/"
        "K1A_prjeb11591_runs.tsv"
    )
    if runs_path.exists():
        runs = pd.read_csv(runs_path, sep="\t")[
            ["run_accession", "sample_alias", "sample_title", "library_strategy"]
        ]
        df = df.merge(runs, on="run_accession", how="left")

    out = OUT_DIR / "T01_K2_genotypes_2field.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"[OK] wrote {out} (n={len(df)} samples, {df.shape[1]} cols)")

    # Quick parse-quality summary
    summary = []
    for locus in LOCI + NONCLASSICAL:
        a1c = df[f"{locus}_a1"].notna().sum()
        a2c = df[f"{locus}_a2"].notna().sum()
        het = ((df[f"{locus}_a1"] != df[f"{locus}_a2"]) &
               df[f"{locus}_a1"].notna() & df[f"{locus}_a2"].notna()).sum()
        both = (df[f"{locus}_a1"].notna() & df[f"{locus}_a2"].notna()).sum()
        summary.append({
            "locus": locus,
            "n_samples_with_a1": a1c,
            "n_samples_with_both": both,
            "n_heterozygous": het,
            "het_rate_observed": het / both if both else None,
        })
    psum = pd.DataFrame(summary)
    psum.to_csv(OUT_DIR / "T02_K2_parse_quality.tsv", sep="\t", index=False)
    print(psum)


if __name__ == "__main__":
    main()
