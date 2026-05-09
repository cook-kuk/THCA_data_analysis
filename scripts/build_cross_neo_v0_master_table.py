#!/usr/bin/env python3
"""Build CROSS-Neo v0 master table."""

from __future__ import annotations

import json

import pandas as pd

from cross_neo_v0_common import INPUT, OUT, ensure_dirs, hla_supertype, peptide_cluster


def main() -> None:
    ensure_dirs()
    bundle = pd.read_csv(INPUT / "bundle.tsv", sep="\t")
    bundle["HLA_norm"] = bundle["HLA_norm"].fillna("")
    # Keep peptide-HLA rows only; Venus rows are protein-level windows with no HLA.
    bundle = bundle[bundle["split"].isin(["train", "ext_itsndb_main", "ext_itsndb_val"])].copy()
    bundle = bundle[bundle["HLA_norm"].astype(str).str.len() > 0].reset_index(drop=True)

    strict = pd.read_csv(INPUT / "wave8b_strict_structure_2026_05_09/strict_no_tcr_self_exact_input.tsv", sep="\t")
    strict_keys = set(zip(strict["peptide"], strict["hla"]))

    wave8 = pd.read_csv(INPUT / "wave8/wave8_combined_features.tsv", sep="\t")
    wave8 = wave8.drop_duplicates(["peptide", "hla"])
    wave8_cols = [
        "peptide",
        "hla",
        "tcr_motif_score",
        "tcr_motif_count_log",
        "tcr_class_score",
        "self_exact_match",
        "self_hamming1_count_log",
        "self_hamming2_count_log",
        "self_blosum_max",
    ]
    m = bundle.merge(
        wave8[wave8_cols],
        left_on=["peptide", "HLA_norm"],
        right_on=["peptide", "hla"],
        how="left",
    )

    rows = []
    for i, r in m.iterrows():
        pep = str(r["peptide"])
        hla = str(r["HLA_norm"])
        study = str(r["source"])
        year = ""
        if study.startswith("ITSNdb"):
            year = "2024"
        row = {
            "sample_id": f"CNV0_{i:05d}",
            "peptide_mut": pep,
            "peptide_wt": "",
            "hla": hla,
            "hla_supertype": hla_supertype(hla),
            "source_protein": str(r.get("protein_id", "") or ""),
            "source_window_15aa": "",
            "source_window_30aa": "",
            "label": int(r["label"]),
            "study": study,
            "patient_id": "",
            "assay_type": "",
            "date_or_publication_year": year,
            "tcr_beta": "",
            "tcr_alpha": "",
            "strict_set_flag": (pep, hla) in strict_keys,
            "public_overlap_flags": "in_master" if bool(r.get("in_master", False)) else "no_inhouse_master_overlap",
            "split": str(r["split"]),
            "in_master": bool(r.get("in_master", False)),
            "near_peptide_cluster": peptide_cluster(pep),
            "tcr_motif_score": float(r.get("tcr_motif_score", 0.0) or 0.0),
            "self_exact_match": float(r.get("self_exact_match", 0.0) or 0.0),
            "self_hamming1_count_log": float(r.get("self_hamming1_count_log", 0.0) or 0.0),
            "self_hamming2_count_log": float(r.get("self_hamming2_count_log", 0.0) or 0.0),
            "self_blosum_max": float(r.get("self_blosum_max", 0.0) or 0.0),
        }
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "master_table.tsv", sep="\t", index=False)

    report = {
        "n_rows": int(len(out)),
        "n_strict": int(out["strict_set_flag"].sum()),
        "n_pos": int(out["label"].sum()),
        "split_counts": out["split"].value_counts().to_dict(),
        "study_counts": out["study"].value_counts().to_dict(),
        "strict_pos": int(out.loc[out["strict_set_flag"], "label"].sum()),
        "missing_columns_reason": {
            "peptide_wt": "not available in current bundle; counterfactual branch falls back to zero WT/delta flags",
            "source windows": "not available in current bundle for ITSNdb/train rows",
            "patient/tcr": "not available in current bundle",
        },
    }
    (OUT / "master_table_report.json").write_text(json.dumps(report, indent=2))
    lines = [
        "# CROSS-Neo v0 Master Table Report",
        "",
        f"Rows: {len(out)}",
        f"Strict rows: {int(out['strict_set_flag'].sum())}; strict positives: {int(out.loc[out['strict_set_flag'], 'label'].sum())}",
        "",
        "Missing WT/source-window/patient/TCR fields are preserved as empty columns because the current bundle does not carry them.",
        "The strict set is the no-overlap, no TCR/self exact, HLA-pseudo-available subset.",
    ]
    (OUT / "master_table_report.md").write_text("\n".join(lines) + "\n")
    print(f"[master] wrote {OUT / 'master_table.tsv'} rows={len(out)} strict={int(out['strict_set_flag'].sum())}")


if __name__ == "__main__":
    main()
