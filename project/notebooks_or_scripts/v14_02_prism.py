"""
v14 Task 4 — PRISM Repurposing drug sensitivity cross-reference.

Uses PRISM Repurposing 19Q4 primary screen (Broad/figshare 9393293):
  primary-screen-replicate-collapsed-logfold-change.csv  (41 MB, rows = cells by
     DepMap/Broad ID, cols = compound×dose replicate-collapsed log-fold-change)
  primary-screen-replicate-collapsed-treatment-info.csv  (compound metadata)

Lower logFC = more cell killing. Per-compound AUC across the (2-3) doses in the
primary screen gives a simple sensitivity proxy per cell line.

Output:
  $RES/ccle_drug_sensitivity.tsv   — long-form per compound × line
  $RES/ccle_drug_sensitivity_by_subtype.tsv — BRAF vs RAS differential
"""
from __future__ import annotations
import json
from pathlib import Path
import pandas as pd
import numpy as np

PROJECT = Path("/opt/thyroid-dash/project")
RES = PROJECT / "results" / "v14_ccle"
PRISM = RES / "prism_cache"

# v7 top candidates + v13 Tier-A anchors of interest
CANDIDATES = [
    "cannabidiol", "cannabinol", "quercetin", "luteolin", "resveratrol",
    "fisetin", "kaempferol", "baicalein",                       # CYP1B1 flavonoids
    "diazepam", "midazolam", "clonazepam", "etomidate", "propofol",
    "zolpidem", "flumazenil",                                   # GABRB2
    "simvastatin", "atorvastatin", "rosuvastatin", "ezetimibe", # LDLR
    "topotecan", "irinotecan",                                  # TACSTD2 payload
    "nafamostat", "camostat", "gabexate",                       # TMPRSS4
    "thalidomide", "lenalidomide", "pomalidomide",              # PLEKHA6 CRBN
    "indomethacin", "dronabinol",                               # v13 Tier-A
]


def main():
    # Treatment info
    tinfo = pd.read_csv(PRISM / "primary-screen-replicate-collapsed-treatment-info.csv")
    print(f"PRISM treatments: {len(tinfo)}")

    mask = tinfo["name"].str.lower().isin([c.lower() for c in CANDIDATES])
    hits = tinfo[mask].copy()
    print(f"Candidates matched: {hits['name'].nunique()} compounds, {len(hits)} (compound×dose) rows")
    print(hits[["name", "broad_id", "dose", "moa", "phase"]].head(30).to_string(index=False))

    # Cell-line info
    cinfo = pd.read_csv(PRISM / "primary-screen-cell-line-info.csv")
    thy = cinfo[cinfo["primary_tissue"] == "thyroid"].copy()
    print(f"\nPRISM thyroid cells: {len(thy)}")
    print(thy[["depmap_id", "ccle_name"]].to_string(index=False))
    thy_ids = thy["depmap_id"].tolist()

    # Stream-parse the logfold-change matrix; rows are cells (ACH-*), cols are column_name keys.
    lfc_path = PRISM / "primary-screen-replicate-collapsed-logfold-change.csv"
    hdr = pd.read_csv(lfc_path, nrows=0).columns.tolist()
    print(f"LFC matrix has {len(hdr)-1} treatment columns")

    # keep only candidate columns
    cand_cols = set(hits["column_name"].tolist())
    keep_cols = [hdr[0]] + [c for c in hdr[1:] if c in cand_cols]
    print(f"Columns to keep: {len(keep_cols)-1} / {len(cand_cols)} candidate (compound×dose)")

    lfc = pd.read_csv(lfc_path, usecols=keep_cols, low_memory=False)
    lfc = lfc.rename(columns={lfc.columns[0]: "depmap_id"})
    thy_lfc = lfc[lfc["depmap_id"].isin(thy_ids)].copy()
    print(f"Thyroid-line rows: {len(thy_lfc)}")

    # Merge with our v14 BRS predictions
    brs = pd.read_csv(RES / "ccle_brs52_validation.tsv", sep="\t")
    brs = brs.merge(
        cinfo[["depmap_id", "ccle_name"]],
        left_on="sample", right_on="ccle_name", how="left"
    )
    # Some samples have slight name drift: map via depmap_id from v8.1 metadata if present
    meta81 = pd.read_csv(RES / "ccle_thyroid_metadata.tsv", sep="\t")
    brs = brs.merge(meta81[["sample_id", "depmap_id"]], left_on="sample", right_on="sample_id",
                    how="left", suffixes=("", "_meta"))
    brs["depmap_id"] = brs["depmap_id_meta"].fillna(brs["depmap_id"])
    thy_lfc = thy_lfc.merge(brs[["depmap_id", "sample", "name", "v14_prediction", "mutation_label"]],
                             on="depmap_id", how="left")

    # Melt to long form
    long_rows = []
    for _, r in thy_lfc.iterrows():
        for col in keep_cols[1:]:
            val = r[col]
            if pd.isna(val):
                continue
            info = hits[hits["column_name"] == col]
            if info.empty:
                continue
            info = info.iloc[0]
            long_rows.append({
                "depmap_id": r["depmap_id"],
                "cell_line": r["name"] if pd.notna(r["name"]) else r["depmap_id"],
                "v14_prediction": r.get("v14_prediction"),
                "mutation_label": r.get("mutation_label"),
                "compound": info["name"],
                "broad_id": info["broad_id"],
                "dose_uM": float(info["dose"]),
                "moa": info["moa"],
                "phase": info["phase"],
                "logfold_change": float(val),
            })
    long_df = pd.DataFrame(long_rows)
    long_df.to_csv(RES / "ccle_drug_sensitivity.tsv", sep="\t", index=False)
    print(f"Wrote ccle_drug_sensitivity.tsv  ({len(long_df)} rows)")

    # Per-compound mean LFC per subtype
    if len(long_df):
        # collapse (compound, cell_line) across doses by median
        per_cell = long_df.groupby(["compound", "cell_line", "v14_prediction", "mutation_label"],
                                    dropna=False)["logfold_change"].median().reset_index()
        # subtype averages
        subtype = per_cell.groupby(["compound", "v14_prediction"])["logfold_change"].mean().unstack()
        if "RAS_like" in subtype.columns and "BRAF_like" in subtype.columns:
            subtype["braf_minus_ras"] = (subtype["BRAF_like"] - subtype["RAS_like"]).round(3)
            subtype = subtype.sort_values("braf_minus_ras")
        subtype.to_csv(RES / "ccle_drug_sensitivity_by_subtype.tsv", sep="\t")
        print("\nCompound × v14-subtype mean LFC (lower = more killing):")
        print(subtype.round(3).to_string())

        # Headline: which compounds show BRAF-selective sensitivity?
        if "braf_minus_ras" in subtype.columns:
            sel = subtype.dropna(subset=["braf_minus_ras"]).copy()
            top_braf_sel = sel[sel["braf_minus_ras"] < -0.5]
            top_ras_sel = sel[sel["braf_minus_ras"] > 0.5]
            print(f"\nCompounds with BRAF-selective killing (ΔLFC < -0.5): {len(top_braf_sel)}")
            print(top_braf_sel.round(3).to_string())
            print(f"\nCompounds with RAS-selective killing (ΔLFC > +0.5): {len(top_ras_sel)}")
            print(top_ras_sel.round(3).to_string())

    # CBD specifically (user called out as key check)
    cbd = long_df[long_df["compound"].str.lower() == "cannabidiol"] if len(long_df) else None
    if cbd is not None and len(cbd):
        print("\nCannabidiol per-line logFC (PRISM):")
        print(cbd[["cell_line", "dose_uM", "logfold_change", "v14_prediction", "mutation_label"]]
              .sort_values("logfold_change").to_string(index=False))
    else:
        print("\nCannabidiol: not present in PRISM 19Q4 primary screen (confirmed via zero-match above).")

    # update summary
    sumfile = RES / "v14_summary.json"
    summary = json.loads(sumfile.read_text()) if sumfile.exists() else {}
    summary["prism_compounds_matched"] = int(hits["name"].nunique())
    summary["prism_thyroid_lines"] = int(len(thy_lfc))
    summary["prism_rows"] = int(len(long_df))
    sumfile.write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
