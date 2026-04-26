#!/usr/bin/env python3
"""v5.2 Task A2 — Table 1 schema fix.

Separately tracks n_rnaseq, n_maf, n_class_A, n_class_B, n_usable per
(cancer, cohort). Resolves the v5.1 inconsistency where LGG showed
n_samples=260 with n_class_A=414 (impossible).

TCGA: RNA-seq = tcga_samples.tsv (patient-level uniq);
      MAF      = tcga_maf_combined.tsv.gz patient_id uniq;
      labels   = tcga_driver_labels.tsv.
GEO:  RNA-seq = geo/{GSE}_expr.tsv columns;
      labels   = geo/{GSE}_pheno.tsv rows with non-null label.
THCA TCGA uses the shortcut (data_processed/bulk_rnaseq/...) — no MAF here.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/opt/thyroid-dash/project/notebooks_or_scripts")
from v5p1_common import COHORTS, CANCERS, DATA_RAW_V5, PROJECT

RES = PROJECT / "results" / "v5p2_fix"
RPT = PROJECT / "reports" / "v5p2"
RES.mkdir(parents=True, exist_ok=True)
RPT.mkdir(parents=True, exist_ok=True)


def tcga_counts(cancer: str) -> dict:
    base = DATA_RAW_V5 / cancer

    if cancer == "THCA":
        # Shortcut path — labels+RNA come from processed master, no MAF recount
        lab = pd.read_csv(base / "tcga_samples_labels.tsv", sep="\t")
        lab["patient_id"] = lab["sample_id"].astype(str).str[:12]
        n_rnaseq = lab["patient_id"].nunique()
        class_a = COHORTS[cancer]["class_a"]
        class_b = COHORTS[cancer]["class_b"]
        n_a = int((lab["label"] == class_a).sum())
        n_b = int((lab["label"] == class_b).sum())
        return dict(n_rnaseq=n_rnaseq, n_maf=None, n_class_A=n_a, n_class_B=n_b,
                    n_usable=n_rnaseq)

    samples = pd.read_csv(base / "tcga_samples.tsv", sep="\t")
    maf = pd.read_csv(base / "tcga_maf_combined.tsv.gz", sep="\t",
                      compression="gzip", low_memory=False)
    labels = pd.read_csv(base / "tcga_driver_labels.tsv", sep="\t")

    n_rnaseq = samples["patient_id"].nunique()
    n_maf = maf["patient_id"].nunique()

    class_a = COHORTS[cancer]["class_a"]
    class_b = COHORTS[cancer]["class_b"]
    n_a = int((labels["label"] == class_a).sum())
    n_b = int((labels["label"] == class_b).sum())

    rnaseq_pids = set(samples["patient_id"])
    labeled_pids = set(labels["patient_id"])
    n_usable = len(rnaseq_pids & labeled_pids)
    return dict(n_rnaseq=n_rnaseq, n_maf=n_maf, n_class_A=n_a, n_class_B=n_b,
                n_usable=n_usable)


def geo_counts(cancer: str, gse: str) -> dict:
    base = DATA_RAW_V5 / cancer / "geo"
    expr_path = base / f"{gse}_expr.tsv"
    pheno_path = base / f"{gse}_pheno.tsv"
    if not expr_path.exists() or not pheno_path.exists():
        return dict(n_rnaseq=0, n_maf=None, n_class_A=0, n_class_B=0, n_usable=0,
                    status="excluded")
    expr = pd.read_csv(expr_path, sep="\t", index_col=0, nrows=1)
    sids = expr.columns.tolist()
    pheno = pd.read_csv(pheno_path, sep="\t")
    n_rnaseq = len(sids)
    class_a = COHORTS[cancer]["class_a"]
    class_b = COHORTS[cancer]["class_b"]
    n_a = int((pheno["label"] == class_a).sum())
    n_b = int((pheno["label"] == class_b).sum())
    labeled = pheno[pheno["label"].notna()]
    n_usable = int(labeled["sample_id"].isin(sids).sum())
    return dict(n_rnaseq=n_rnaseq, n_maf=None, n_class_A=n_a, n_class_B=n_b,
                n_usable=n_usable, status="included")


def main():
    rows = []
    for cancer in CANCERS:
        t = tcga_counts(cancer)
        rows.append(dict(cancer=cancer, cohort_type="tcga",
                         cohort_id=COHORTS[cancer]["tcga"], **t))
        for gse in COHORTS[cancer]["geo"]:
            g = geo_counts(cancer, gse)
            g.pop("status", None)
            rows.append(dict(cancer=cancer, cohort_type="geo", cohort_id=gse, **g))

    df = pd.DataFrame(rows, columns=[
        "cancer", "cohort_type", "cohort_id",
        "n_rnaseq", "n_maf", "n_class_A", "n_class_B", "n_usable",
    ])
    out = RES / "cohort_availability_detailed.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"wrote {out} ({len(df)} rows)")

    # TOTALs + DIAL-usable rows per cancer (only cohorts actually used in harmonization)
    harm = pd.read_csv(PROJECT / "results" / "v5" / "v5p1_harmonization.tsv", sep="\t")

    # Redesigned Table 1 combining detail + harmonization TOTAL
    table_rows = []
    for cancer in CANCERS:
        sub = df[df["cancer"] == cancer]
        for _, r in sub.iterrows():
            nm = r["n_maf"]
            nm_s = "—" if pd.isna(nm) else str(int(nm))
            table_rows.append({
                "Cancer": cancer,
                "Cohort": r["cohort_id"],
                "Type": r["cohort_type"].upper(),
                "RNA-seq n": int(r["n_rnaseq"]),
                "MAF n": nm_s,
                "Class A": int(r["n_class_A"]),
                "Class B": int(r["n_class_B"]),
                "Usable (RNA ∩ label)": int(r["n_usable"]),
            })
        # harmonization row (total used in DIAL)
        hrow = harm[harm["cancer"] == cancer]
        if not hrow.empty:
            hr = hrow.iloc[0]
            table_rows.append({
                "Cancer": cancer,
                "Cohort": "DIAL cohort (post-harmonize)",
                "Type": "TOTAL",
                "RNA-seq n": int(hr["n_samples"]),
                "MAF n": "—",
                "Class A": int(hr["n_class_A"]),
                "Class B": int(hr["n_class_B"]),
                "Usable (RNA ∩ label)": int(hr["n_samples"]),
            })

    tbl = pd.DataFrame(table_rows)
    tbl_out = RES / "table1_redesigned.tsv"
    tbl.to_csv(tbl_out, sep="\t", index=False)
    print(f"wrote {tbl_out} ({len(tbl)} rows)")

    # LaTeX version
    tex_lines = [
        r"% Table 1 (v5.2) — cohort availability with separated counts.",
        r"% Columns: RNA-seq n, MAF n, Class A, Class B, Usable (RNA∩label).",
        r"\begin{tabular}{llrrrrr}",
        r"\hline",
        r"Cancer & Cohort & RNA-seq & MAF & Class A & Class B & Usable \\",
        r"\hline",
    ]
    for _, r in tbl.iterrows():
        cohort = str(r["Cohort"]).replace("_", r"\_")
        tex_lines.append(
            f'{r["Cancer"]} & {cohort} & {r["RNA-seq n"]} & {r["MAF n"]} '
            f'& {r["Class A"]} & {r["Class B"]} & {r["Usable (RNA ∩ label)"]} \\\\'
        )
    tex_lines += [
        r"\hline",
        r"\end{tabular}",
        "",
        r"\vspace{0.5em}",
        r"\noindent\textit{Note: n\_usable (intersection of RNA-seq and labelled patients) "
        r"is used for DIAL computation. RNA-seq and MAF counts are reported for cohort "
        r"transparency; they differ because GDC distributes RNA-seq and MAF under "
        r"separate cohort budgets. For TCGA-LGG the GEO cohorts were excluded due to "
        r"platform mapping failures, so the DIAL cohort reported is the TCGA-LGG "
        r"TSS split (Tissue Source Site) with the four largest TSS codes as batches.}",
    ]
    tex_out = RPT / "table1_redesigned.tex"
    tex_out.write_text("\n".join(tex_lines))
    print(f"wrote {tex_out}")

    # Show to console
    print("\n=== Table 1 redesigned ===")
    print(tbl.to_string(index=False))


if __name__ == "__main__":
    main()
