#!/usr/bin/env python3
"""Create CROSS-Neo 2.0 locked split files."""

from __future__ import annotations

import json
import pandas as pd

from common import OUT, V0, ensure_dirs, write_json


def main() -> None:
    ensure_dirs()
    reg = pd.read_csv(OUT / "canonical_registry.tsv", sep="\t")
    split_dir = OUT / "splits"
    rows = []
    if (V0 / "folds.tsv").exists():
        folds = pd.read_csv(V0 / "folds.tsv", sep="\t")
        strict_ids = set(reg.loc[reg["strict_set_flag"].astype(bool), "row_id"])
        for (split, fold), g in folds.groupby(["split_name", "fold_id"]):
            test = set(g["sample_id"].astype(str))
            train = strict_ids - test
            for rid in train:
                rows.append({"split_name": split, "fold_id": fold, "row_id": rid, "role": "train"})
            for rid in test:
                rows.append({"split_name": split, "fold_id": fold, "row_id": rid, "role": "test"})
    sources = ["CEDAR", "NEPdb", "TESLA_mmc4", "TESLA_mmc7_validation"]
    pool = reg[reg["source_dataset"].isin(sources)]
    for src in sources:
        split = f"source_heldout_{src}"
        for rid in pool.loc[~pool["source_dataset"].eq(src), "row_id"]:
            rows.append({"split_name": split, "fold_id": src, "row_id": rid, "role": "train"})
        for rid in pool.loc[pool["source_dataset"].eq(src), "row_id"]:
            rows.append({"split_name": split, "fold_id": src, "row_id": rid, "role": "test"})
    # Low-prevalence stress: TESLA sources as test.
    tesla = set(pool.loc[pool["source_dataset"].astype(str).str.contains("TESLA"), "row_id"])
    for rid in set(pool["row_id"]) - tesla:
        rows.append({"split_name": "low_prevalence_stress_split", "fold_id": "tesla_all", "row_id": rid, "role": "train"})
    for rid in tesla:
        rows.append({"split_name": "low_prevalence_stress_split", "fold_id": "tesla_all", "row_id": rid, "role": "test"})
    df = pd.DataFrame(rows)
    for split, g in df.groupby("split_name"):
        g.to_csv(split_dir / f"{split}.tsv", sep="\t", index=False, na_rep="NA")
    manifest = {
        "seed": 20260509,
        "splits": {
            split: {
                "folds": int(g["fold_id"].nunique()),
                "train_rows": int((g["role"] == "train").sum()),
                "test_rows": int((g["role"] == "test").sum()),
            }
            for split, g in df.groupby("split_name")
        },
        "selection_rule": "outer test labels are not used for model/fusion selection",
    }
    write_json(split_dir / "split_manifest.json", manifest)
    lines = ["# CROSS-Neo 2.0 Split Report", "", pd.DataFrame(manifest["splits"]).T.reset_index(names="split_name").to_markdown(index=False)]
    (split_dir / "split_report.md").write_text("\n".join(lines) + "\n")
    print(f"[v2-splits] rows={len(df)} splits={df['split_name'].nunique()}")


if __name__ == "__main__":
    main()
