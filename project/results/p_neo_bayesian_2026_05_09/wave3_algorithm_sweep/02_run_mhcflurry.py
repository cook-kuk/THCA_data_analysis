#!/usr/bin/env python3
"""Score ITSNdb scoring set with MHCflurry 2.0 Class1PresentationPredictor."""
import pandas as pd
from pathlib import Path
from mhcflurry import Class1PresentationPredictor

ROOT = Path(__file__).resolve().parent
df = pd.read_csv(ROOT / "scoring_set_itsndb.tsv", sep="\t")
print(f"scoring n={len(df)}")

predictor = Class1PresentationPredictor.load()
sup = set(predictor.supported_alleles)
print(f"mhcflurry supports {len(sup)} alleles; first 10: {list(sorted(sup))[:10]}")

# Filter to supported alleles
df["mhcflurry_supported"] = df["HLA_norm"].isin(sup)
n_drop = (~df["mhcflurry_supported"]).sum()
print(f"unsupported: {n_drop}")

df_ok = df[df["mhcflurry_supported"]].copy()
peps = df_ok["peptide"].tolist()
alleles = df_ok["HLA_norm"].tolist()

# Per-peptide-per-allele prediction.
# Use sample_names + dict alleles so each peptide is scored ONLY against its assigned allele.
sample_ids = [f"s{i}" for i in range(len(peps))]
alleles_dict = {sid: [allele] for sid, allele in zip(sample_ids, alleles)}

out = predictor.predict(
    peptides=peps,
    sample_names=sample_ids,
    alleles=alleles_dict,
    verbose=0,
)
print("output columns:", list(out.columns))
print(out.head(3))

# Re-align with df_ok order (peptide_num column gives original index)
df_ok = df_ok.reset_index(drop=True)
df_ok["sample_id"] = sample_ids
out_keyed = out.set_index(["sample_name", "peptide_num"])[
    ["presentation_score", "affinity", "processing_score"]
]
df_ok["mhcflurry_presentation"] = [
    out_keyed.loc[(sid, i), "presentation_score"]
    for i, sid in enumerate(sample_ids)
]
df_ok["mhcflurry_affinity"] = [
    out_keyed.loc[(sid, i), "affinity"]
    for i, sid in enumerate(sample_ids)
]
df_ok["mhcflurry_processing"] = [
    out_keyed.loc[(sid, i), "processing_score"]
    for i, sid in enumerate(sample_ids)
]

# Re-merge with full set including unsupported (they get NaN scores)
df_full = df.merge(
    df_ok[["peptide", "HLA_norm", "mhcflurry_presentation", "mhcflurry_affinity", "mhcflurry_processing"]],
    on=["peptide", "HLA_norm"], how="left"
)

df_full.to_csv(ROOT / "mhcflurry_scored.tsv", sep="\t", index=False)
print(f"saved {ROOT / 'mhcflurry_scored.tsv'}")
print(f"scored: {df_full['mhcflurry_presentation'].notna().sum()} / {len(df_full)}")
