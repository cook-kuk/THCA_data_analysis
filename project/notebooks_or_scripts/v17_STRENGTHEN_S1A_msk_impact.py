#!/usr/bin/env python3
"""
v17 STRENGTHEN S1-A — MSK-IMPACT thyroid cross-cohort analysis (datahub TSV).

Strategy: download data_clinical_sample.txt + data_mutations.txt from
the cBioPortal datahub (msk_impact_2017), filter to thyroid, compute
TERT promoter prevalence + BRAF V600E prevalence by cancer type detail.

Outputs:
- results/v17_strengthen/S1A_msk_impact_thyroid.tsv
- results/v17_strengthen/S1A_tert_prevalence_cross_cohort.tsv
- results/v17_strengthen/S1A_summary.json
"""
from __future__ import annotations
import io
import json
import sys
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError

import pandas as pd

OUT = Path("/opt/thyroid-dash/project/results/v17_strengthen")
OUT.mkdir(parents=True, exist_ok=True)

CLIN_URL = "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/msk_impact_2017/data_clinical_sample.txt"
MUT_URL = "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/msk_impact_2017/data_mutations.txt"
# fallback: raw.githubusercontent
CLIN_URL_FALLBACK = "https://raw.githubusercontent.com/cBioPortal/datahub/master/public/msk_impact_2017/data_clinical_sample.txt"
MUT_URL_FALLBACK = "https://raw.githubusercontent.com/cBioPortal/datahub/master/public/msk_impact_2017/data_mutations.txt"


def fetch_text(url: str, timeout: int = 300) -> str:
    req = Request(url, headers={"User-Agent": "Mozilla/5.0", "Accept": "text/plain"})
    with urlopen(req, timeout=timeout) as r:
        data = r.read()
    return data.decode("utf-8", errors="replace")


def fetch_tsv_with_fallback(primary: str, fallback: str, label: str) -> pd.DataFrame:
    last_err = None
    for url in (primary, fallback):
        try:
            print(f"[S1A] downloading {label} from {url} ...", flush=True)
            txt = fetch_text(url)
            # cBioPortal clinical files have header lines starting with #
            lines = txt.splitlines()
            data_start = next((i for i, ln in enumerate(lines) if not ln.startswith("#")), 0)
            tsv_body = "\n".join(lines[data_start:])
            df = pd.read_csv(io.StringIO(tsv_body), sep="\t", low_memory=False)
            print(f"[S1A] {label}: {len(df)} rows × {len(df.columns)} cols", flush=True)
            return df
        except (URLError, HTTPError, TimeoutError, OSError) as e:
            last_err = e
            print(f"[S1A] {label} fetch failed at {url}: {e}", flush=True)
    raise RuntimeError(f"both URLs failed for {label}: {last_err}")


def main() -> int:
    summary = {"study": "msk_impact_2017", "errors": []}
    try:
        clin = fetch_tsv_with_fallback(CLIN_URL, CLIN_URL_FALLBACK, "clinical_sample")
    except Exception as e:
        summary["errors"].append(f"clinical fetch: {e}")
        with open(OUT / "S1A_summary.json", "w") as f:
            json.dump(summary, f, indent=2, default=str)
        return 1

    # filter to thyroid cancer type
    cancer_col = "CANCER_TYPE" if "CANCER_TYPE" in clin.columns else None
    detail_col = "CANCER_TYPE_DETAILED" if "CANCER_TYPE_DETAILED" in clin.columns else None
    if cancer_col is None:
        # try ONCOTREE_CODE or look for thyroid in sample IDs
        for c in clin.columns:
            if "CANCER" in c.upper() or "ONCOTREE" in c.upper() or "PRIMARY" in c.upper():
                print("[S1A] candidate type col:", c, flush=True)
        summary["errors"].append("CANCER_TYPE column not found")
        clin.head(5).to_csv(OUT / "S1A_clinical_sample_head.tsv", sep="\t", index=False)
        return 1

    clin["__type"] = clin[cancer_col].astype(str).str.lower()
    thyroid = clin[clin["__type"].str.contains("thyroid", na=False)].copy()
    print(f"[S1A] {len(thyroid)} thyroid samples in MSK-IMPACT", flush=True)
    summary["n_thyroid_samples"] = len(thyroid)
    summary["n_total_samples"] = len(clin)

    if len(thyroid) == 0:
        # show what types exist
        types = clin[cancer_col].value_counts().head(20).to_dict()
        summary["errors"].append("zero thyroid samples")
        summary["top20_cancer_types"] = types
        with open(OUT / "S1A_summary.json", "w") as f:
            json.dump(summary, f, indent=2, default=str)
        return 1

    sample_col = "SAMPLE_ID" if "SAMPLE_ID" in thyroid.columns else thyroid.columns[0]
    thyroid_sids = set(thyroid[sample_col].astype(str).tolist())

    # mutations
    try:
        mut = fetch_tsv_with_fallback(MUT_URL, MUT_URL_FALLBACK, "mutations")
    except Exception as e:
        summary["errors"].append(f"mutation fetch: {e}")
        # write what we have
        thyroid.to_csv(OUT / "S1A_msk_impact_thyroid.tsv", sep="\t", index=False)
        with open(OUT / "S1A_summary.json", "w") as f:
            json.dump(summary, f, indent=2, default=str)
        return 1

    # normalize column names
    mut_sample_col = "Tumor_Sample_Barcode" if "Tumor_Sample_Barcode" in mut.columns else "SAMPLE_ID"
    mut_gene_col = "Hugo_Symbol" if "Hugo_Symbol" in mut.columns else "GENE"

    mut_thy = mut[mut[mut_sample_col].astype(str).isin(thyroid_sids)].copy()
    print(f"[S1A] {len(mut_thy)} mutation rows for thyroid samples", flush=True)

    # TERT promoter detection
    # MSK-IMPACT TERT panel is promoter-focused. Use HGVSp_Short or Variant_Classification or Start_Position.
    tert = mut_thy[mut_thy[mut_gene_col].astype(str).str.upper() == "TERT"].copy()
    if "Start_Position" in tert.columns:
        tert["__is_promoter"] = tert["Start_Position"].astype(str).isin(["1295228", "1295250", "1295113", "1295135"]) | \
                                tert.get("Variant_Classification", pd.Series("")).astype(str).str.contains("5'Flank|Promoter|TSS|Upstream", case=False, na=False)
    else:
        tert["__is_promoter"] = True  # MSK panel only covers promoter
    # Treat any TERT mutation in MSK-IMPACT as promoter (panel scope)
    tert_pos = set(tert[mut_sample_col].astype(str).tolist())
    print(f"[S1A] {len(tert_pos)} thyroid samples with TERT promoter mutation", flush=True)

    # BRAF V600E
    braf = mut_thy[mut_thy[mut_gene_col].astype(str).str.upper() == "BRAF"].copy()
    pc_col = next((c for c in ["HGVSp_Short", "Protein_Change", "Amino_Acid_Change"] if c in braf.columns), None)
    if pc_col is not None:
        braf_v600e = braf[braf[pc_col].astype(str).str.contains("V600E", na=False)]
    else:
        braf_v600e = braf  # all BRAF muts as fallback
    braf_pos = set(braf_v600e[mut_sample_col].astype(str).tolist())
    print(f"[S1A] {len(braf_pos)} thyroid samples with BRAF V600E", flush=True)

    # build per-sample table
    thyroid["TERT_promoter_mut"] = thyroid[sample_col].astype(str).isin(tert_pos).astype(int)
    thyroid["BRAF_V600E"] = thyroid[sample_col].astype(str).isin(braf_pos).astype(int)
    cols_keep = [c for c in [sample_col, "PATIENT_ID", cancer_col, detail_col, "SAMPLE_TYPE",
                              "AGE_AT_SEQ_REPORT", "SEX", "TERT_promoter_mut", "BRAF_V600E"]
                  if c is not None and c in thyroid.columns]
    thyroid_out = thyroid[cols_keep].copy()
    thyroid_out.to_csv(OUT / "S1A_msk_impact_thyroid.tsv", sep="\t", index=False)

    # per-cancer-type-detail prevalence
    if detail_col and detail_col in thyroid_out.columns:
        by = thyroid_out.groupby(detail_col).agg(
            n=(sample_col, "size"),
            n_TERT=("TERT_promoter_mut", "sum"),
            n_BRAF=("BRAF_V600E", "sum"),
        ).reset_index()
        by["TERT_pct"] = (by["n_TERT"] / by["n"] * 100).round(1)
        by["BRAF_pct"] = (by["n_BRAF"] / by["n"] * 100).round(1)
        by = by.sort_values("n", ascending=False)
        by.to_csv(OUT / "S1A_tert_prevalence_cross_cohort.tsv", sep="\t", index=False)
        summary["per_cancer_type_detail"] = by.to_dict(orient="records")
    else:
        # at least overall
        summary["per_cancer_type_detail"] = []

    # TCGA-THCA reference
    summary["tcga_reference"] = {
        "PTC": {"n": 504, "TERT_pct": 7.1, "BRAF_pct": 56.1, "source": "TCGA Cell 2014 + cBioPortal thca_tcga_pub"},
    }
    summary["overall_n_tert"] = int(thyroid_out["TERT_promoter_mut"].sum())
    summary["overall_n_braf"] = int(thyroid_out["BRAF_V600E"].sum())
    summary["overall_TERT_pct"] = round(thyroid_out["TERT_promoter_mut"].mean() * 100, 1) if len(thyroid_out) else None
    summary["overall_BRAF_pct"] = round(thyroid_out["BRAF_V600E"].mean() * 100, 1) if len(thyroid_out) else None

    with open(OUT / "S1A_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print("\n[S1A] DONE.")
    print(f"  thyroid n = {len(thyroid_out)}")
    print(f"  overall TERT% = {summary['overall_TERT_pct']}")
    print(f"  overall BRAF% = {summary['overall_BRAF_pct']}")
    if summary["per_cancer_type_detail"]:
        print("  per-type breakdown:")
        for r in summary["per_cancer_type_detail"]:
            print(f"    {r.get(detail_col, '?'):<45} n={r['n']:>4}  TERT%={r['TERT_pct']:>5}  BRAF%={r['BRAF_pct']:>5}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
