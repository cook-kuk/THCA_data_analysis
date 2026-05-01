"""P2-B Strategy C — Yoo 2016 PLoS Genet S6 → K2 mutation calls + NBNR labels."""
from pathlib import Path
import re
import json
import pandas as pd

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/dark_matter_phase2")
SRC = OUT / "yoo2016_S6_clinical.xlsx"

df = pd.read_excel(SRC, sheet_name="S6 Table", header=1)
df.columns = [c.split("\n")[0].strip() for c in df.columns]
df = df[df["SampleID"].notna() & (df["SampleID"] != "SampleID")].copy()
print(f"Patients: {len(df)}")
print(f"Columns: {list(df.columns)[:8]}")

def parse_driver(s):
    if pd.isna(s) or str(s).strip() in ("", "."):
        return pd.Series({"genes": "", "has_braf_v600e": False, "has_ras": False,
                          "has_tert": False, "has_dicer1": False, "has_eif1ax": False,
                          "has_fusion": False, "alt_drivers": ""})
    s = str(s)
    has_braf = bool(re.search(r"BRAF.*V600E|BRAF\s*p\.?V600E", s, re.IGNORECASE))
    has_ras = bool(re.search(r"\b(NRAS|HRAS|KRAS)\b", s))
    has_tert = bool(re.search(r"TERT", s, re.IGNORECASE))
    has_dicer1 = bool(re.search(r"DICER1", s, re.IGNORECASE))
    has_eif1ax = bool(re.search(r"EIF1AX", s, re.IGNORECASE))
    has_fusion = bool(re.search(r"\b(RET|NTRK\d?|ALK|ROS1|ETV6)\b.*(fusion|/|::)|fusion|PAX8.PPARG", s, re.IGNORECASE))
    gene_matches = re.findall(r"\b(BRAF|NRAS|HRAS|KRAS|TERT|DICER1|EIF1AX|RET|NTRK\d?|ALK|ROS1|PAX8|PPARG|SOS1|SPOP|IDH1|PTEN|TP53|CHEK2|PPM1D)\b", s)
    genes = sorted(set(gene_matches))
    alt = [g for g in genes if g not in {"BRAF", "NRAS", "HRAS", "KRAS", "TERT"}]
    return pd.Series({"genes": ";".join(genes), "has_braf_v600e": has_braf, "has_ras": has_ras,
                      "has_tert": has_tert, "has_dicer1": has_dicer1, "has_eif1ax": has_eif1ax,
                      "has_fusion": has_fusion, "alt_drivers": ";".join(alt)})

parsed = df["Driver mutation"].apply(parse_driver)
df = pd.concat([df, parsed], axis=1)
df["is_dark_matter"] = (~df["has_braf_v600e"]) & (~df["has_ras"])
df["mol_subtype_label"] = df["Molecular subtype"].map({1: "NBNR", 2: "RAS-like", 3: "BRAF-like",
                                                       1.0: "NBNR", 2.0: "RAS-like", 3.0: "BRAF-like"}).fillna("Unknown")

n = len(df); n_dm = int(df["is_dark_matter"].sum())
print(f"\n=== K2 (Yoo 2016, n={n}) Mutation Landscape ===")
print(f"  BRAF V600E+: {int(df['has_braf_v600e'].sum())} ({df['has_braf_v600e'].sum()/n*100:.1f}%)")
print(f"  RAS+: {int(df['has_ras'].sum())} ({df['has_ras'].sum()/n*100:.1f}%)")
print(f"  TERT+: {int(df['has_tert'].sum())}")
print(f"  DICER1+: {int(df['has_dicer1'].sum())}")
print(f"  EIF1AX+: {int(df['has_eif1ax'].sum())}")
print(f"  Fusion+: {int(df['has_fusion'].sum())}")
print(f"  **Dark Matter (BRAF−/RAS−): {n_dm} ({n_dm/n*100:.1f}%)** | TCGA ref: 28.4%")

print("\nMolecular subtype × DM:")
print(pd.crosstab(df["mol_subtype_label"], df["is_dark_matter"]))
print("\nPathology distribution:")
print(df["Pathology"].value_counts())

dm = df[df["is_dark_matter"]]
print(f"\n=== K2 Dark Matter (n={len(dm)}) ===")
print(f"  DICER1+: {int(dm['has_dicer1'].sum())} ({dm['has_dicer1'].sum()/len(dm)*100:.1f}%)")
print(f"  EIF1AX+: {int(dm['has_eif1ax'].sum())} ({dm['has_eif1ax'].sum()/len(dm)*100:.1f}%)")
print(f"  Fusion+: {int(dm['has_fusion'].sum())}")
print(f"  TERT+: {int(dm['has_tert'].sum())}")
print(f"\nDM × Pathology:\n{pd.crosstab(df['is_dark_matter'], df['Pathology'])}")

df.to_csv(OUT / "k2_yoo2016_mutations_parsed.tsv", sep="\t", index=False)

summary = {
    "source": "Yoo SK et al. PLoS Genet 2016 PMID 27494611, Supplementary S6",
    "n_total": int(n),
    "mutations": {k: int(df[f"has_{k.lower()}"].sum() if k.lower() in ("braf_v600e","ras","tert","dicer1","eif1ax","fusion") else 0)
                  for k in ["BRAF_V600E", "RAS", "TERT", "DICER1", "EIF1AX", "Fusion"]},
    "dark_matter": {"n": int(n_dm), "pct": round(n_dm/n*100, 2), "tcga_reference_pct": 28.42},
    "molecular_subtype": {str(k): int(v) for k, v in df["mol_subtype_label"].value_counts().to_dict().items()},
    "dark_matter_alt_driver": {
        "DICER1": int(dm["has_dicer1"].sum()), "EIF1AX": int(dm["has_eif1ax"].sum()),
        "Fusion": int(dm["has_fusion"].sum()), "TERT": int(dm["has_tert"].sum()),
    },
}
(OUT / "k2_mutation_summary.json").write_text(json.dumps(summary, indent=2))
print(f"\nSaved: k2_yoo2016_mutations_parsed.tsv + k2_mutation_summary.json")
