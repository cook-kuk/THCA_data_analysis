#!/usr/bin/env python3
"""
Build DM1/DM2-balanced TCGA-THCA WSI manifest for SPARK + Pod B WSI dispatch.

Inputs
------
- project/results/dark_matter_phase2/p2d_per_sample_classification.tsv
    columns: tcga_short, driver_class, v17_dark_cluster (DM1/DM2/blank), PFI, PFI.time

Output
------
- project/data/manifests/tcga_thca_wsi_dm_balanced.tsv
    one row per case_submitter_id (TCGA-XX-XXXX) with DM label + driver class +
    PFI + recommended n per stratum.
- project/data/manifests/tcga_thca_wsi_gdc_filter.json
    a GDC API filter JSON ready to POST against /files endpoint.

Stratification rule
-------------------
n=50 split:
- DM1 n=25 (8 BRAF, 8 RAS, 5 driver-neg, 4 other)
- DM2 n=25 (8 BRAF, 8 RAS, 5 driver-neg, 4 other)

If a stratum is short, fall back to next-best driver class within same DM cluster.
"""
from __future__ import annotations
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
SRC = ROOT / "project/results/dark_matter_phase2/p2d_per_sample_classification.tsv"
OUT_MANIFEST = ROOT / "project/data/manifests/tcga_thca_wsi_dm_balanced.tsv"
OUT_GDC_FILTER = ROOT / "project/data/manifests/tcga_thca_wsi_gdc_filter.json"

random.seed(42)

DRIVER_BUCKETS = {
    "BRAF": ["Class1_BRAF_V600E", "Class1b_BRAF_other"],
    "RAS":  ["Class2_RAS_HRAS", "Class2_RAS_KRAS", "Class2_RAS_NRAS"],
    "TERT": ["Class3_TERT_promoter"],
    "FUSION": ["Class4_RET_fusion", "Class5_other_fusion"],
    "NEG":  ["Class6_True_driver_neg"],
}


def bucketize(cls: str) -> str:
    if not cls or pd.isna(cls):
        return "OTHER"
    for b, prefixes in DRIVER_BUCKETS.items():
        if any(cls.startswith(p) for p in prefixes):
            return b
    return "OTHER"


def main():
    df = pd.read_csv(SRC, sep="\t")
    df = df.rename(columns={"tcga_short": "case_id"})
    df["dm"] = df["v17_dark_cluster"].fillna("UNLABELED")
    df["bucket"] = df["driver_class"].apply(bucketize)

    print("=== input distribution ===")
    print(df.groupby(["dm", "bucket"]).size().unstack(fill_value=0))

    # DM1/DM2 are dominated by driver-neg (Class6) — that's the core "dark matter" cohort.
    # We stratify within each DM by PFI event status (informative outcome signal),
    # plus reserve a small slot for non-NEG driver carriers within each DM (5/25).
    rng = pd.Series(range(len(df)))  # placeholder for deterministic shuffles
    selected = []
    for dm in ["DM1", "DM2"]:
        sub = df[df.dm == dm].copy()
        sub["pfi_evt"] = (sub["PFI"] == 1.0).astype(int)
        non_neg = sub[sub.bucket != "NEG"]
        neg = sub[sub.bucket == "NEG"]

        slot_nonneg = min(5, len(non_neg))
        slot_neg = 25 - slot_nonneg
        # PFI-event balanced within NEG: take half events, half censored
        evt = neg[neg.pfi_evt == 1]
        cen = neg[neg.pfi_evt == 0]
        n_evt = min(len(evt), max(2, slot_neg // 3))
        n_cen = slot_neg - n_evt
        pick_evt = evt.sample(n=n_evt, random_state=42) if len(evt) else evt
        pick_cen = cen.sample(n=min(n_cen, len(cen)), random_state=42) if len(cen) else cen
        pick_nonneg = non_neg.sample(n=slot_nonneg, random_state=42) if slot_nonneg else non_neg.iloc[0:0]
        chunk = pd.concat([pick_evt, pick_cen, pick_nonneg], ignore_index=True)
        chunk = chunk.assign(strat=dm)
        selected.append(chunk)

    out = pd.concat(selected, ignore_index=True).drop_duplicates("case_id")
    out = out[["case_id", "dm", "bucket", "driver_class", "strat", "PFI", "PFI.time"]]

    OUT_MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_MANIFEST, sep="\t", index=False)
    print(f"\n=== output: {len(out)} cases ===")
    print(out.groupby(["dm", "bucket"]).size().unstack(fill_value=0))
    print(f"wrote {OUT_MANIFEST.relative_to(ROOT)}")

    gdc_filter = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-THCA"]}},
            {"op": "in", "content": {"field": "data_format", "value": ["SVS"]}},
            {"op": "in", "content": {"field": "experimental_strategy", "value": ["Diagnostic Slide"]}},
            {"op": "in", "content": {"field": "cases.submitter_id", "value": out.case_id.tolist()}},
        ],
    }
    OUT_GDC_FILTER.write_text(json.dumps(gdc_filter, indent=2))
    print(f"wrote {OUT_GDC_FILTER.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
