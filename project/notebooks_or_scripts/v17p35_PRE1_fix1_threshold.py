"""
v17p35 PRE-1 — FIX1 threshold rerun.

Problem: friend's FIX1 used a calibrated LogReg classifier that returned
prob_dm1 ~ 1.0 for all 13 CCLE thyroid lines → no DM2-like comparison
group → drug actionability undefined.

Fix: use raw gene-set mean z-scores (DM1-up minus DM2-up) and median-split
the 13 thyroid lines. Validate with BRAF mutation status (sanity check),
then re-run PRISM Mann-Whitney.

Outputs:
  results/v17p35/tables/FIX1_celline_dm_scores_v2.tsv
  results/v17p35/tables/FIX1_top_drugs_dm1_selective_v2.tsv
  results/v17p35/tables/FIX1_top_drugs_dm2_selective_v2.tsv
  results/v17p35/tables/FIX1_threshold_method_audit.tsv
  results/v17p35/figs/FIX1_drug_volcano_v2.html
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ranksums

PROJECT = Path("/opt/thyroid-dash/project")
RES = PROJECT / "results" / "v17p35"
TAB = RES / "tables"
FIG = RES / "figs"
PRISM = PROJECT / "results" / "v14_ccle" / "prism_cache"

np.random.seed(42)


def load_dm_markers():
    df = pd.read_csv(PROJECT / "results/v17/tables/dark_matter_cluster_markers.tsv", sep="\t")
    # cluster 0 = DM1 (MAPK-up, dedifferentiated)
    # cluster 1 = DM2 (well-differentiated)
    dm1_up = df[(df["cluster"] == 0) & (df["log2fc"] > 0)]["gene"].tolist()
    dm2_up = df[(df["cluster"] == 1) & (df["log2fc"] > 0)]["gene"].tolist()
    return dm1_up, dm2_up


def main():
    expr = pd.read_csv(PROJECT / "results/v14_ccle/ccle_thyroid_expression.tsv",
                       sep="\t", index_col=0)
    print(f"CCLE expression: {expr.shape}")

    # DM markers
    dm1_up, dm2_up = load_dm_markers()
    print(f"DM1-up markers: {len(dm1_up)}  DM2-up markers: {len(dm2_up)}")

    # Within-cohort z-score
    Z = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1).replace(0, 1), axis=0)
    dm1_in = [g for g in dm1_up if g in Z.index]
    dm2_in = [g for g in dm2_up if g in Z.index]
    print(f"  DM1 markers present: {len(dm1_in)}  DM2 markers present: {len(dm2_in)}")

    dm1_score = Z.loc[dm1_in].mean(axis=0)
    dm2_score = Z.loc[dm2_in].mean(axis=0)
    delta = dm1_score - dm2_score   # >0 = DM1-leaning

    # Median split
    median = delta.median()
    dm_label = (delta > median).map({True: "DM1_like", False: "DM2_like"})

    # BRAF sanity check
    truth = pd.read_csv(PROJECT / "results/v14_ccle/ccle_thyroid_mutation_truth.tsv", sep="\t")
    truth_map = dict(zip(truth["sample_id"], truth["mutation_label"]))

    score_df = pd.DataFrame({
        "sample": delta.index,
        "dm1_score_z": dm1_score.round(3),
        "dm2_score_z": dm2_score.round(3),
        "delta_z": delta.round(3),
        "dm_like_v2": dm_label,
        "mutation_label": [truth_map.get(s, "other") for s in delta.index],
    })
    score_df = score_df.sort_values("delta_z", ascending=False).reset_index(drop=True)

    # Sanity: do BRAF mutants concentrate on DM1 side?
    braf_lines = score_df[score_df["mutation_label"] == "BRAF"]
    braf_dm1 = (braf_lines["dm_like_v2"] == "DM1_like").sum()
    print(f"\n=== BRAF sanity ===")
    print(f"BRAF-mutant lines: {len(braf_lines)}, in DM1 cluster: {braf_dm1}/{len(braf_lines)}")
    print(score_df.to_string(index=False))

    score_df.to_csv(TAB / "FIX1_celline_dm_scores_v2.tsv", sep="\t", index=False)

    # method audit
    audit = pd.DataFrame([{
        "method": "median_split_on_delta_z(dm1_z - dm2_z)",
        "n_dm1": int((dm_label == "DM1_like").sum()),
        "n_dm2": int((dm_label == "DM2_like").sum()),
        "median_delta_z": round(median, 3),
        "braf_in_dm1": int(braf_dm1),
        "braf_total": int(len(braf_lines)),
        "ras_in_dm1": int(((score_df["mutation_label"] == "RAS") & (score_df["dm_like_v2"] == "DM1_like")).sum()),
        "ras_total": int((score_df["mutation_label"] == "RAS").sum()),
    }])
    audit.to_csv(TAB / "FIX1_threshold_method_audit.tsv", sep="\t", index=False)
    print(f"\nAudit:\n{audit.to_string(index=False)}")

    # ---- PRISM Mann-Whitney ----
    print("\n=== PRISM Mann-Whitney (DM1 vs DM2) ===")
    cellinfo = pd.read_csv(PRISM / "primary-screen-cell-line-info.csv")
    tinfo = pd.read_csv(PRISM / "primary-screen-replicate-collapsed-treatment-info.csv")

    # Map sample → depmap_id
    sid_map = dict(zip(cellinfo["ccle_name"], cellinfo["depmap_id"]))
    score_df["depmap_id"] = score_df["sample"].map(sid_map)
    valid = score_df.dropna(subset=["depmap_id"]).copy()
    print(f"  Lines with depmap_id: {len(valid)}/{len(score_df)}")

    # Load LFC for these cells (filtered)
    lfc_path = PRISM / "primary-screen-replicate-collapsed-logfold-change.csv"
    hdr = pd.read_csv(lfc_path, nrows=0).columns.tolist()
    lfc = pd.read_csv(lfc_path, usecols=hdr, low_memory=False)
    lfc = lfc.rename(columns={lfc.columns[0]: "depmap_id"})
    lfc_thy = lfc[lfc["depmap_id"].isin(valid["depmap_id"])].set_index("depmap_id")
    print(f"  LFC rows for thyroid lines: {len(lfc_thy)}")

    # group cells by DM label (using delta_z median)
    valid_with_label = valid[valid["depmap_id"].isin(lfc_thy.index)].set_index("depmap_id")
    dm1_cells = valid_with_label[valid_with_label["dm_like_v2"] == "DM1_like"].index.tolist()
    dm2_cells = valid_with_label[valid_with_label["dm_like_v2"] == "DM2_like"].index.tolist()
    print(f"  DM1 cells in PRISM: {len(dm1_cells)}  DM2 cells in PRISM: {len(dm2_cells)}")

    if len(dm1_cells) < 2 or len(dm2_cells) < 2:
        print("  ! Not enough cells in one group for Mann-Whitney")
        return

    # Per-compound, collapse across doses to median LFC, then Mann-Whitney
    drug_results = []
    tinfo_map = tinfo.set_index("column_name")
    for col in lfc_thy.columns:
        if col not in tinfo_map.index:
            continue
        meta = tinfo_map.loc[col]
        if isinstance(meta, pd.DataFrame):
            meta = meta.iloc[0]
        dm1_vals = lfc_thy.loc[dm1_cells, col].dropna().values
        dm2_vals = lfc_thy.loc[dm2_cells, col].dropna().values
        if len(dm1_vals) < 2 or len(dm2_vals) < 2:
            continue
        try:
            stat, p = mannwhitneyu(dm1_vals, dm2_vals, alternative="two-sided")
        except Exception:
            continue
        drug_results.append({
            "column_name": col,
            "compound": meta["name"] if "name" in meta else "",
            "broad_id": meta.get("broad_id", ""),
            "dose_uM": meta.get("dose", np.nan),
            "moa": meta.get("moa", ""),
            "phase": meta.get("phase", ""),
            "dm1_n": len(dm1_vals),
            "dm2_n": len(dm2_vals),
            "dm1_mean_lfc": float(np.mean(dm1_vals)),
            "dm2_mean_lfc": float(np.mean(dm2_vals)),
            "delta_lfc": float(np.mean(dm1_vals) - np.mean(dm2_vals)),
            "pvalue": float(p),
        })

    drug_df = pd.DataFrame(drug_results)
    print(f"  Drugs tested: {len(drug_df)}")
    if len(drug_df) == 0:
        return
    # BH correction
    from scipy.stats import false_discovery_control
    try:
        drug_df["fdr"] = false_discovery_control(drug_df["pvalue"].values, method="bh")
    except Exception:
        # manual BH
        m = len(drug_df)
        rk = drug_df["pvalue"].rank(method="min")
        drug_df["fdr"] = (drug_df["pvalue"] * m / rk).clip(upper=1.0)

    # Collapse per compound (multiple doses → take min p)
    drug_df = drug_df.sort_values("pvalue")
    drug_df_compound = drug_df.groupby("compound", as_index=False).first().sort_values("pvalue")

    n_fdr_lt_01 = int((drug_df_compound["fdr"] < 0.1).sum())
    n_fdr_lt_05 = int((drug_df_compound["fdr"] < 0.05).sum())
    print(f"\n  Compounds: {len(drug_df_compound)}, FDR<0.1: {n_fdr_lt_01}, FDR<0.05: {n_fdr_lt_05}")

    # Top DM1-selective (lower mean_lfc in DM1 → more killing)
    dm1_sel = drug_df_compound[drug_df_compound["delta_lfc"] < 0].sort_values("pvalue").head(50)
    dm2_sel = drug_df_compound[drug_df_compound["delta_lfc"] > 0].sort_values("pvalue").head(50)

    dm1_sel.to_csv(TAB / "FIX1_top_drugs_dm1_selective_v2.tsv", sep="\t", index=False)
    dm2_sel.to_csv(TAB / "FIX1_top_drugs_dm2_selective_v2.tsv", sep="\t", index=False)

    # MOA enrichment (simple count)
    moa_dm1 = dm1_sel.head(20)["moa"].value_counts().head(15).to_frame("dm1_count")
    moa_dm2 = dm2_sel.head(20)["moa"].value_counts().head(15).to_frame("dm2_count")
    moa_combined = pd.concat([moa_dm1, moa_dm2], axis=1).fillna(0).astype(int).reset_index().rename(columns={"index": "moa"})
    moa_combined.to_csv(TAB / "FIX1_moa_enrichment_v2.tsv", sep="\t", index=False)

    # Sanity: check MEK / BRAF / Topo-I family
    print("\n=== Mechanism sanity ===")
    for keyword, label in [("MEK", "MEK inhibitor"), ("BRAF", "BRAF inhibitor"),
                            ("topoisomerase", "Topo-I"), ("HMGCR", "Statin")]:
        m = drug_df_compound[drug_df_compound["moa"].str.contains(keyword, case=False, na=False)]
        if len(m) > 0:
            print(f"  {label}: {len(m)} compounds, top-3:")
            print(m.head(3)[["compound", "delta_lfc", "pvalue", "moa"]].to_string(index=False))

    # Save summary v2
    summary = {
        "method": "median_split_delta_z",
        "n_cell_lines": int(len(valid_with_label)),
        "n_dm1_cells": int(len(dm1_cells)),
        "n_dm2_cells": int(len(dm2_cells)),
        "n_drugs_fdr_lt_0_1": n_fdr_lt_01,
        "n_drugs_fdr_lt_0_05": n_fdr_lt_05,
        "braf_in_dm1": int(braf_dm1),
        "braf_total": int(len(braf_lines)),
        "top_dm1_drug": dm1_sel.iloc[0]["compound"] if len(dm1_sel) else None,
        "top_dm1_moa": dm1_sel.iloc[0]["moa"] if len(dm1_sel) else None,
        "top_dm2_drug": dm2_sel.iloc[0]["compound"] if len(dm2_sel) else None,
    }
    (TAB / "FIX1_summary_v2.json").write_text(json.dumps(summary, indent=2))
    print(f"\n=== SUMMARY ===\n{json.dumps(summary, indent=2)}")


if __name__ == "__main__":
    main()
