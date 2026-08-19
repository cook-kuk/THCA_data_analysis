#!/usr/bin/env python3
"""Independent re-derivation of the two conflicting GSE151179 RAI claims.

Reads the raw series matrix + platform annotation directly (no reuse of any
prior repo output file) and reproduces:

  Claim 1 (manuscript / Fig 6D):  whole-series "collection before/after rai"
           contrast on thyroid_differentiation module, all tissue types pooled.
  Claim 2 (2026-08-06 audit):     "rai uptake at the metastatic site" contrast,
           tumour-only, patient-level, with tumour-purity adjustment.

Also runs the tissue-type decomposition that explains why Claim 1 and Claim 2
diverge.
"""
import gzip
import re
import numpy as np
import pandas as pd
from scipy import stats

RAW = "/data/rai_atlas/raw/GSE151179"
SM = f"{RAW}/GSE151179_series_matrix.txt.gz"
PROBE2GENE = f"{RAW}/GPL23159_probe2gene.tsv"

# ---------------------------------------------------------------------------
# 1. Parse sample characteristics from the series matrix header block
# ---------------------------------------------------------------------------
def parse_series_matrix_meta(path):
    meta_lines = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_geo_accession"):
                meta_lines["geo_accession"] = line.rstrip("\n").split("\t")[1:]
            if line.startswith("!Sample_characteristics_ch1"):
                meta_lines.setdefault("characteristics", []).append(
                    line.rstrip("\n").split("\t")[1:]
                )
            if line.startswith("!Sample_title"):
                meta_lines["title"] = line.rstrip("\n").split("\t")[1:]
    gsm = [x.strip('"') for x in meta_lines["geo_accession"]]
    title = [x.strip('"') for x in meta_lines.get("title", [""] * len(gsm))]
    n = len(gsm)
    df = pd.DataFrame({"gsm": gsm, "title": title})
    for row in meta_lines["characteristics"]:
        row = [x.strip('"') for x in row]
        # field label is before the first ": "
        labels = set()
        for cell in row:
            if ": " in cell:
                labels.add(cell.split(": ", 1)[0])
        assert len(labels) <= 1, f"mixed labels in one characteristics row: {labels}"
        if not labels:
            continue
        label = labels.pop()
        colname = label.strip().lower().replace(" ", "_").replace("/", "_")
        vals = []
        for cell in row:
            if ": " in cell:
                vals.append(cell.split(": ", 1)[1])
            else:
                vals.append(np.nan)
        assert len(vals) == n
        df[colname] = vals
    return df


meta = parse_series_matrix_meta(SM)
print("=== Sample characteristic fields found ===")
for c in meta.columns:
    if c in ("gsm", "title"):
        continue
    print(f"  {c}: {sorted(meta[c].dropna().unique().tolist())}")
print(f"\nTotal samples (GSMs): {len(meta)}")

import os
_OUTDIR = os.path.dirname(os.path.abspath(__file__))
meta.to_csv(os.path.join(_OUTDIR, "gse151179_full_metadata_2026_08_19.tsv"), sep="\t", index=False)

# ---------------------------------------------------------------------------
# 2. Parse expression matrix, restrict to the panel genes we need
# ---------------------------------------------------------------------------
probe2gene = pd.read_csv(PROBE2GENE, sep="\t")
probe2gene = probe2gene.dropna(subset=["gene_symbol"])
p2g = dict(zip(probe2gene["ID"], probe2gene["gene_symbol"]))

THYROID_DIFF_MODULE = ["SLC5A5", "TPO", "TSHR", "TG", "PAX8", "NKX2-1", "DIO1", "DIO2"]
DM8_PANEL = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "DIO1", "DIO2", "FOXE1"]
need_genes = set(THYROID_DIFF_MODULE) | set(DM8_PANEL)

wanted_probes = {pid for pid, g in p2g.items() if g in need_genes}
print(f"\nProbes mapping to needed genes: {len(wanted_probes)} (need {len(need_genes)} genes)")

rows = []
gsm_order = None
with gzip.open(SM, "rt") as fh:
    in_table = False
    for line in fh:
        if line.startswith("!series_matrix_table_begin"):
            in_table = True
            header = next(fh)
            gsm_order = [x.strip().strip('"') for x in header.rstrip("\n").split("\t")[1:]]
            continue
        if line.startswith("!series_matrix_table_end"):
            break
        if in_table:
            parts = line.rstrip("\n").split("\t")
            pid = parts[0].strip('"')
            if pid in wanted_probes:
                vals = [float(x) for x in parts[1:]]
                rows.append((pid, vals))

expr = pd.DataFrame({pid: vals for pid, vals in rows}, index=gsm_order).T
expr["gene_symbol"] = expr.index.map(p2g)
# collapse multiple probes per gene by mean (matches conventional practice)
gene_expr = expr.groupby("gene_symbol").mean(numeric_only=True).T  # samples x genes
gene_expr.index.name = "gsm"
print(f"\nGene-level matrix: {gene_expr.shape[0]} samples x {gene_expr.shape[1]} genes")
print("Genes recovered:", sorted(gene_expr.columns.tolist()))

# within-cohort z-score per gene (matches repo convention: "within-cohort z-mean")
gene_z = (gene_expr - gene_expr.mean()) / gene_expr.std(ddof=1)


def module_score(gene_z_df, genes):
    keep = [g for g in genes if g in gene_z_df.columns]
    missing = [g for g in genes if g not in gene_z_df.columns]
    if missing:
        print(f"  [module_score] missing genes on array: {missing}")
    return gene_z_df[keep].mean(axis=1)


meta = meta.set_index("gsm")
meta["thyroid_diff_module_z"] = module_score(gene_z, THYROID_DIFF_MODULE)
meta["dm8_panel_z"] = module_score(gene_z, DM8_PANEL)


def cohens_d(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    sa, sb = np.var(a, ddof=1), np.var(b, ddof=1)
    pooled = np.sqrt(((len(a) - 1) * sa + (len(b) - 1) * sb) / (len(a) + len(b) - 2))
    return (np.mean(a) - np.mean(b)) / pooled


def report_contrast(name, high, low, col, test="mwu"):
    a = meta.loc[high, col].dropna()
    b = meta.loc[low, col].dropna()
    d = cohens_d(a, b)
    if test == "mwu":
        _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    else:
        _, p = stats.ttest_ind(a, b, equal_var=False)
    print(f"{name}: n_high={len(a)} n_low={len(b)}  d={d:.3f}  P({test})={p:.4g}")
    return {"contrast": name, "n_high": len(a), "n_low": len(b), "d": d, "p": p, "test": test}


print("\n" + "=" * 80)
print("CLAIM 1 — manuscript Fig 6D wording: 'collection before/after rai', ALL samples, thyroid_diff module")
print("=" * 80)
before_mask = meta["collection_before_after_rai"].str.contains("before", case=False, na=False)
after_mask = meta["collection_before_after_rai"].str.contains("after", case=False, na=False)
print("before n =", before_mask.sum(), " after n =", after_mask.sum())
print("tissue_type crosstab (before/after x tissue):")
print(pd.crosstab(
    np.where(before_mask, "before", np.where(after_mask, "after", "?")),
    meta["tissue_type"],
))

c1_t = report_contrast(
    "C1a all-tissue after-vs-before, thyroid_diff_module, Welch t-test (matches repo script exactly)",
    meta.index[after_mask], meta.index[before_mask], "thyroid_diff_module_z", test="ttest",
)
c1_mwu = report_contrast(
    "C1b all-tissue after-vs-before, thyroid_diff_module, Mann-Whitney (matches manuscript's stated test)",
    meta.index[after_mask], meta.index[before_mask], "thyroid_diff_module_z", test="mwu",
)

print("\n--- decomposition: does the effect survive holding tissue type constant? ---")
tumor_mask = meta["tissue_type"].isin(["primary tumour", "primary tumor"]) | meta["tissue_type"].str.contains("primary", case=False, na=False)
lnmet_mask = meta["tissue_type"].str.contains("lymph node metastasis", case=False, na=False)
normal_mask = meta["tissue_type"].str.contains("non-neoplastic|normal", case=False, na=False)
print("primary tumour n =", tumor_mask.sum(), " LN-met (any) n =", lnmet_mask.sum(), " normal thyroid n =", normal_mask.sum())

# tumours-only (drop normal thyroid) before/after
tumors_only = meta.index[(before_mask | after_mask) & ~normal_mask]
c1_tumoronly = report_contrast(
    "C1c tumours-only (drop non-neoplastic thyroid) after-vs-before, thyroid_diff_module, MWU",
    meta.index[after_mask & meta.index.isin(tumors_only)],
    meta.index[before_mask & meta.index.isin(tumors_only)],
    "thyroid_diff_module_z", test="mwu",
)

# LN-met only (tissue held constant): pre-RAI synchronous LN met vs post-RAI LN met
presync_lnmet = meta.index[before_mask & lnmet_mask]
postlnmet = meta.index[after_mask & lnmet_mask]
c1_lnmetonly = report_contrast(
    "C1d LN-metastasis-only (tissue held constant) after-vs-before, thyroid_diff_module, MWU",
    postlnmet, presync_lnmet, "thyroid_diff_module_z", test="mwu",
)

# normal thyroid vs LN met (the confound)
c1_confound = report_contrast(
    "C1e CONFOUND: non-neoplastic thyroid vs lymph-node-metastasis (any), thyroid_diff_module, MWU",
    meta.index[lnmet_mask], meta.index[normal_mask], "thyroid_diff_module_z", test="mwu",
)

print("\n" + "=" * 80)
print("CLAIM 2 — audit wording: 'rai uptake at the metastatic site', tumour-only, patient-level")
print("=" * 80)
print("rai_uptake_at_the_metastatic_site value counts (all samples):")
print(meta["rai_uptake_at_the_metastatic_site"].value_counts(dropna=False))

uptake_col = "rai_uptake_at_the_metastatic_site"
has_uptake = meta[uptake_col].notna()
print("\nTissue type of samples carrying an uptake label:")
print(meta.loc[has_uptake, "tissue_type"].value_counts())
print("\nBefore/after status of samples carrying an uptake label:")
print(meta.loc[has_uptake, "collection_before_after_rai"].value_counts())

yes_mask = meta[uptake_col].str.lower().eq("yes")
no_mask = meta[uptake_col].str.lower().eq("no")
c2_sample = report_contrast(
    "C2a sample-level No-uptake vs Yes-uptake, dm8_panel_z, MWU",
    meta.index[no_mask], meta.index[yes_mask], "dm8_panel_z", test="mwu",
)
c2_sample_diff = report_contrast(
    "C2b sample-level No-uptake vs Yes-uptake, thyroid_diff_module_z, MWU",
    meta.index[no_mask], meta.index[yes_mask], "thyroid_diff_module_z", test="mwu",
)

# patient-level de-duplication (one row per patient_id, since some patients
# contribute >1 tumour specimen with the same label per audit 04b)
patient_col = "patient_id"
pat = meta.reset_index()
pat_level = (
    pat[pat[uptake_col].notna()]
    .groupby(patient_col)
    .agg({uptake_col: "first", "dm8_panel_z": "mean", "thyroid_diff_module_z": "mean",
          "tumor_purity_class_by_cibersort": "first"})
)
yes_p = pat_level[pat_level[uptake_col].str.lower() == "yes"]
no_p = pat_level[pat_level[uptake_col].str.lower() == "no"]
d_pat = cohens_d(no_p["dm8_panel_z"], yes_p["dm8_panel_z"])
_, p_pat = stats.mannwhitneyu(no_p["dm8_panel_z"].dropna(), yes_p["dm8_panel_z"].dropna(), alternative="two-sided")
print(f"\nC2c PATIENT-LEVEL (one row/patient) No-uptake vs Yes-uptake, dm8_panel_z:")
print(f"  n_no={len(no_p)} n_yes={len(yes_p)}  d={d_pat:.3f}  MWU P={p_pat:.4g}")

print("\ntumor purity class by uptake status:")
print(pd.crosstab(pat_level[uptake_col], pat_level["tumor_purity_class_by_cibersort"]))

# ---------------------------------------------------------------------------
# C2d: exact independent replication of audit script
# gse151179_uptake_at_met_site_2026_08_06.py -- PANEL_8, tumour-only
# (Primary tumor + all lymph-node-metastasis categories), patient-level.
# ---------------------------------------------------------------------------
print("\n--- C2d exact-methodology replication (tumour-only incl. LN-mets, patient-level) ---")
AUDIT_PANEL_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
audit_panel_z = module_score(gene_z, AUDIT_PANEL_8)
meta["audit_panel_z"] = audit_panel_z
is_tumor = meta["tissue_type"].str.contains("primary tumor|lymph node metastasis", case=False, na=False)
tum = meta[is_tumor].copy()
print(f"tumour samples (primary tumor + any LN-met category): n={len(tum)} "
      f"(expect 39 per audit 04b)")
pat = tum.reset_index().groupby(["patient_id", uptake_col], as_index=False)["audit_panel_z"].mean()
dup = pat.patient_id.duplicated(keep=False)
if dup.any():
    print("  discordant-label patients:", pat.loc[dup, "patient_id"].tolist())
else:
    print("  no patient carries a discordant uptake label (confirms patient-level field)")
yes_v = pat.loc[pat[uptake_col] == "Yes", "audit_panel_z"].values
no_v = pat.loc[pat[uptake_col] == "No", "audit_panel_z"].values
d_exact = cohens_d(yes_v, no_v)
_, p_exact = stats.mannwhitneyu(yes_v, no_v, alternative="two-sided")
print(f"  n_patients Yes={len(yes_v)} No={len(no_v)}  d(Yes-No)={d_exact:.3f}  MWU P={p_exact:.4g}")
print("  (compare to audit 04b: n=32 patients, d=0.373, P=0.534)")

# negative control: non-neoplastic thyroid samples carrying the same patient-level label
print("\n--- negative control: non-neoplastic thyroid samples by (patient-level) uptake label ---")
normal = meta[normal_mask].copy()
yes_n = normal.loc[normal[uptake_col] == "Yes", "audit_panel_z"].dropna().values
no_n = normal.loc[normal[uptake_col] == "No", "audit_panel_z"].dropna().values
if len(yes_n) >= 2 and len(no_n) >= 2:
    d_neg = cohens_d(yes_n, no_n)
    _, p_neg = stats.mannwhitneyu(yes_n, no_n, alternative="two-sided")
    print(f"  n_yes={len(yes_n)} n_no={len(no_n)}  d(Yes-No)={d_neg:.3f}  MWU P={p_neg:.4g}")
    print("  -> non-tumour tissue should show no biological association with met-site uptake;")
    print("     a nonzero effect here demonstrates the label is patient-level, not specimen-level")

print("\nDone. Full metadata written to gse151179_full_metadata.tsv")
