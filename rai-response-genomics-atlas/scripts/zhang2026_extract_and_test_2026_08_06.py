#!/usr/bin/env python3
"""Zhang 2026 (Cell Rep Med) — extract the per-patient advanced-DTC table and test drivers vs RAI.

Zhang T, Cui W, Tang H, ... Shi X. "Proteogenomic characterization delineates clinically
relevant subtypes of advanced differentiated thyroid cancer." Cell Reports Medicine
2026;7(3):102661, doi 10.1016/j.xcrm.2026.102661, PMID 41794039, PMCID PMC13006413.

The proteomic matrix (iProX IPX0011848000) and sequencing (GSA HRA011340) are both gated,
but the supplemental PDF is open and contains genuine per-patient tables:

  Table S1  113 patients: CC subtype, radioiodine therapy given, **RAI sensitivity
            (Avid / Refractory)**, **biochemical response to RAI (G1/G2/G3)**, histology,
            age, sex, ECOG, disease status.
  Table S2  per-patient somatic mutations (gene, cDNA, protein, VAF).

Retrieval that worked without contacting anyone:
    https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13006413/supplementaryFiles
returns a zip with mmc1.pdf; `pdftotext -layout` renders each patient as one line.

This script extracts both tables and runs the same driver-versus-refractoriness test that
was run on Siraj 2022, giving a second independent advanced-disease cohort.

Outputs:
  results/tables/zhang2026_patient_table_2026_08_06.tsv
  results/tables/zhang2026_mutations_2026_08_06.tsv
  results/tables/zhang2026_driver_vs_rai_2026_08_06.tsv
  results/figures/figure_zhang2026_driver_rai_2026_08_06.{png,pdf}
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
PDF = Path("/data/rai_atlas/external/supp_bulk/mmc1.pdf")
FIG = ROOT / "results" / "figures"
TAB = ROOT / "results" / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)
STAMP = "2026_08_06"

S1_COLS = ["patient", "age", "sex", "ecog", "immune_status", "disease_status", "cc_subtype",
           "proteomics", "ngs", "digital_pathology", "rai_therapy", "rai_sensitivity",
           "biochemical_response", "histology"]
S2_COLS = ["patient", "gene", "cdna", "protein", "vaf"]
RAS_GENES = {"NRAS", "HRAS", "KRAS"}


def pdf_text():
    if not PDF.exists():
        raise SystemExit(f"{PDF} not found — fetch the Europe PMC supplementary zip first")
    return subprocess.run(["pdftotext", "-layout", str(PDF), "-"],
                          capture_output=True, text=True, timeout=300).stdout


def parse_tables(txt):
    s1_start, s2_start = txt.index("Table S1."), txt.index("Table S2.")
    s3_start = txt.index("Table S3.")

    def rows(block):
        return [re.split(r"\s{2,}", l.strip())
                for l in block.split("\n") if re.match(r"\s*FUSCC-\d+", l)]

    s1 = [r for r in rows(txt[s1_start:s2_start]) if len(r) == len(S1_COLS)]
    s1 = pd.DataFrame(s1, columns=S1_COLS)

    s2_raw = rows(txt[s2_start:s3_start])
    s2 = pd.DataFrame([r for r in s2_raw if len(r) == len(S2_COLS)], columns=S2_COLS)
    return s1, s2


def driver_class(s2, patients):
    m = s2.copy()
    m["gene"] = m["gene"].str.upper().str.strip()
    m["protein"] = m["protein"].str.upper().str.strip()
    braf = set(m.loc[m.gene.eq("BRAF") & m.protein.str.contains("V600E"), "patient"])
    ras = set(m.loc[m.gene.isin(RAS_GENES)
                    & m.protein.str.contains(r"G12|G13|Q61", regex=True), "patient"])
    out = {}
    for p in patients:
        out[p] = "BRAF V600E" if p in braf else ("RAS hotspot" if p in ras
                                                 else "BRAF/RAS-negative")
    return pd.Series(out)


def main():
    txt = pdf_text()
    s1, s2 = parse_tables(txt)
    print(f"Table S1 parsed: {len(s1)} patients")
    print(f"Table S2 parsed: {len(s2)} mutation rows, "
          f"{s2['patient'].nunique()} patients with sequencing")

    for c in ["age"]:
        s1[c] = pd.to_numeric(s1[c], errors="coerce")
    s1["rai_sensitivity"] = s1["rai_sensitivity"].str.strip()
    s1["biochemical_response"] = s1["biochemical_response"].str.strip()
    print("\nRAI sensitivity:", s1["rai_sensitivity"].value_counts(dropna=False).to_dict())
    print("Biochemical response:", s1["biochemical_response"].value_counts(dropna=False).to_dict())
    print("CC subtype:", s1["cc_subtype"].value_counts().to_dict())
    print("Histology:", s1["histology"].value_counts().to_dict())

    s1["driver"] = s1["patient"].map(driver_class(s2, s1["patient"]))
    print("Driver class:", s1["driver"].value_counts().to_dict())

    s1.to_csv(TAB / f"zhang2026_patient_table_{STAMP}.tsv", sep="\t", index=False)
    s2.to_csv(TAB / f"zhang2026_mutations_{STAMP}.tsv", sep="\t", index=False)

    ev = s1[s1["rai_sensitivity"].isin(["Avid", "Refractory"])].copy()
    ev["refractory"] = (ev["rai_sensitivity"] == "Refractory").astype(int)
    print(f"\nevaluable for RAI sensitivity: {len(ev)} "
          f"({ev['rai_sensitivity'].value_counts().to_dict()})")

    rows = []

    ct = pd.crosstab(ev["driver"], ev["rai_sensitivity"])
    if ct.shape[0] > 1:
        chi = stats.chi2_contingency(ct.values)
        rows.append(dict(test="driver class x RAI refractoriness (chi-square)",
                         detail=str(ct.to_dict()), n=int(ct.values.sum()),
                         statistic=float(chi.statistic), p=float(chi.pvalue)))
        print("\ndriver x RAI:\n", ct.to_string())
        for d in ct.index:
            a = int(ct.loc[d].get("Refractory", 0)); b = int(ct.loc[d].get("Avid", 0))
            orr, p = stats.fisher_exact([[a, b],
                                         [int(ct["Refractory"].sum()) - a,
                                          int(ct["Avid"].sum()) - b]])
            rows.append(dict(test=f"{d} vs rest — refractory", detail=f"{a}/{a+b}",
                             n=int(ct.values.sum()), statistic=float(orr), p=float(p)))
            print(f"  {d:20s} refractory {a}/{a+b} ({100*a/max(a+b,1):.0f}%)  "
                  f"OR={orr:.2f} P={p:.4f}")

    ct2 = pd.crosstab(ev["cc_subtype"], ev["rai_sensitivity"])
    if ct2.shape[0] > 1:
        chi2 = stats.chi2_contingency(ct2.values)
        rows.append(dict(test="proteomic CC subtype x RAI refractoriness (chi-square)",
                         detail=str(ct2.to_dict()), n=int(ct2.values.sum()),
                         statistic=float(chi2.statistic), p=float(chi2.pvalue)))
        print("\nCC subtype x RAI:\n", ct2.to_string())
        print(f"  chi-square P = {chi2.pvalue:.3g}")

    ct3 = pd.crosstab(ev["histology"], ev["rai_sensitivity"])
    if ct3.shape[0] > 1:
        chi3 = stats.chi2_contingency(ct3.values)
        rows.append(dict(test="histology x RAI refractoriness (chi-square)",
                         detail=str(ct3.to_dict()), n=int(ct3.values.sum()),
                         statistic=float(chi3.statistic), p=float(chi3.pvalue)))
        print("\nhistology x RAI:\n", ct3.to_string())

    res = pd.DataFrame(rows)
    res.to_csv(TAB / f"zhang2026_driver_vs_rai_{STAMP}.tsv", sep="\t", index=False)

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0), facecolor="white", layout="constrained")
    cR = "#9c4742"

    def bars(ax, table, title, order=None):
        idx = order or list(table.index)
        idx = [i for i in idx if i in table.index]
        frac = [table.loc[i].get("Refractory", 0) / table.loc[i].sum() * 100 for i in idx]
        ns = [int(table.loc[i].sum()) for i in idx]
        ax.bar(np.arange(len(idx)), frac, color=cR, alpha=0.85, edgecolor="white")
        overall = table["Refractory"].sum() / table.values.sum() * 100
        ax.axhline(overall, color="#444", ls="--", lw=0.9)
        ax.set_xticks(np.arange(len(idx)))
        ax.set_xticklabels([str(i).replace(" ", "\n") for i in idx], fontsize=8.5)
        for k, (f_, n_) in enumerate(zip(frac, ns)):
            ax.text(k, f_ + 1.5, f"{f_:.0f}%\nn={n_}", ha="center", fontsize=8.5, color="#333")
        ax.set_ylabel("% radioiodine-refractory")
        ax.set_ylim(0, max(frac) * 1.28 + 6)
        ax.set_title(title, fontsize=10.5, loc="left", fontweight="bold")
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
        ax.text(0.99, 0.02, f"cohort {overall:.0f}%", transform=ax.transAxes, ha="right",
                fontsize=8.5, color="#444")

    def pv(name):
        r = res[res.test.str.startswith(name)]
        return float(r["p"].iloc[0]) if len(r) else float("nan")

    bars(axes[0], ct, f"a · By driver mutation class\nchi-square P = {pv('driver class'):.3g}",
         ["BRAF V600E", "RAS hotspot", "BRAF/RAS-negative"])
    bars(axes[1], ct2, "b · By proteomic consensus subtype\n"
                       f"chi-square P = {pv('proteomic CC subtype'):.3g}")
    bars(axes[2], ct3, f"c · By histology\nchi-square P = {pv('histology'):.3g}")

    fig.suptitle("Zhang 2026 advanced differentiated thyroid cancer (n = %d) — "
                 "radioiodine refractoriness by molecular and histological class" % len(ev),
                 fontsize=12.5, fontweight="bold")
    fig.savefig(FIG / f"figure_zhang2026_driver_rai_{STAMP}.png", dpi=180,
                bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"figure_zhang2026_driver_rai_{STAMP}.pdf",
                bbox_inches="tight", facecolor="white")
    print(f"\nwrote {FIG / f'figure_zhang2026_driver_rai_{STAMP}.png'}")


if __name__ == "__main__":
    main()
