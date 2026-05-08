#!/usr/bin/env python3
"""
COMMOT (Cao et al., Nat Methods 2023) cell-cell communication on GSE250521 Visium ST.

Upgrade from our hand-rolled lag correlation: COMMOT runs collective optimal
transport over a curated L-R database (CellChatDB) so we get formal sender/receiver
flux for each L-R pair at single-spot resolution, with spatial distance as cost.

For each slide: top-3 L-R pairs by total communication strength, per-stage rollup.

Outputs:
  spark_commot_per_slide_LR_strength.tsv     — slide × L-R pair × total strength
  spark_commot_per_stage_summary.tsv         — per-stage top L-R pairs ranking
  spark_commot_DM1_quartile_diff.tsv         — DM1 top-25% vs bot-25% L-R activity diff
"""
from __future__ import annotations
import sys, os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import commot as ct

ROOT = Path(__file__).resolve().parent.parent.parent.parent
GS = ROOT / "project/data/processed/GSE250521"
RES = ROOT / "project/results/03_pathology_poc"

# Subset DB to most relevant pathways for tumor-immune crosstalk
PATHWAY_FOCUS = [
    "TGFb", "PD-L1", "CXCL", "CCL", "IFN-II", "IL10", "TNF", "FGF", "WNT",
    "VEGF", "EGF", "PDGF", "GAS6", "LGALS9", "FN1", "COLLAGEN", "LAMININ",
]


def main():
    # CellChatDB — try human secreted signaling
    try:
        df_lig_rec = ct.pp.ligand_receptor_database(species="human", database="CellChat", signaling_type="Secreted Signaling")
    except Exception as e:
        print(f"DB fail: {e}; falling back to built-in")
        df_lig_rec = ct.pp.ligand_receptor_database(species="human")
    print(f"L-R database: {len(df_lig_rec)} entries, columns: {list(df_lig_rec.columns)}")
    if "pathway_name" in df_lig_rec.columns:
        df_lig_rec = df_lig_rec[df_lig_rec.pathway_name.isin(PATHWAY_FOCUS)]
        print(f"after pathway filter: {len(df_lig_rec)}")

    sids = sorted([d.name for d in GS.iterdir() if (d / f"{d.name}.scored.h5ad").exists()])
    rows_strength = []
    rows_dm1 = []

    for sid in sids:
        a = ad.read_h5ad(GS / sid / f"{sid}.scored.h5ad")
        stage = a.obs.stage.iloc[0] if "stage" in a.obs.columns else "?"
        if "spatial" not in a.obsm or a.n_obs < 50:
            continue
        a.obsm["spatial"] = np.asarray(a.obsm["spatial"], dtype=float)
        # Filter to genes present in L-R DB
        lr_genes = set()
        for c in df_lig_rec.columns[:2]:  # first two cols typically ligand/receptor
            for entry in df_lig_rec[c].astype(str):
                for g in entry.split("_"):
                    lr_genes.add(g)
        keep_var = [g for g in a.var.index if g in lr_genes]
        if len(keep_var) < 10:
            print(f"  {sid}: not enough L-R genes; skip")
            continue
        # Run COMMOT
        try:
            ct.tl.spatial_communication(
                a, database_name="CellChat", df_ligrec=df_lig_rec,
                dis_thr=500, heteromeric=True,
                pathway_sum=True,
            )
        except Exception as e:
            print(f"  {sid}: COMMOT fail ({type(e).__name__}: {str(e)[:100]}); skip")
            continue
        print(f"  {sid} ({stage}, {a.n_obs} spots): COMMOT OK")

        # Extract per-pathway sender / receiver totals from a.obsp / a.uns
        # Each pathway → 's-pathway_name' and 'r-pathway_name' obs columns
        s_cols = [c for c in a.obs.columns if c.startswith("s-CellChat-")]
        r_cols = [c for c in a.obs.columns if c.startswith("r-CellChat-")]
        v_dm1 = a.obs["DM1_like_score"].values
        q = pd.qcut(v_dm1, 4, labels=False, duplicates="drop")
        for c in s_cols:
            pname = c.replace("s-CellChat-", "")
            send = a.obs[c].values
            recv = a.obs[c.replace("s-", "r-")].values if c.replace("s-", "r-") in a.obs.columns else None
            tot = float(np.nansum(send))
            tot_recv = float(np.nansum(recv)) if recv is not None else 0.0
            rows_strength.append({
                "sample_id": sid, "stage": stage, "pathway": pname,
                "n_spots": int(a.n_obs),
                "total_send": tot, "total_recv": tot_recv,
                "mean_per_spot": float(np.nanmean(send)),
            })
            # DM1 quartile diff
            top = send[q == 3]
            bot = send[q == 0]
            if len(top) >= 5 and len(bot) >= 5:
                rows_dm1.append({
                    "sample_id": sid, "stage": stage, "pathway": pname,
                    "DM1_top_send_mean": float(np.nanmean(top)),
                    "DM1_bot_send_mean": float(np.nanmean(bot)),
                    "DM1_top_minus_bot": float(np.nanmean(top) - np.nanmean(bot)),
                })

    if rows_strength:
        df = pd.DataFrame(rows_strength)
        df.to_csv(RES / "spark_commot_per_slide_LR_strength.tsv", sep="\t", index=False)
        # Per-stage top pathways
        stg = df.groupby(["stage", "pathway"]).agg(
            n_slides=("sample_id", "count"),
            total_send_mean=("total_send", "mean"),
            mean_per_spot_mean=("mean_per_spot", "mean"),
        ).reset_index().sort_values(["stage", "total_send_mean"], ascending=[True, False])
        stg.to_csv(RES / "spark_commot_per_stage_summary.tsv", sep="\t", index=False)
        print("\n=== top 3 pathways per stage (COMMOT total_send) ===")
        for s in ["PT", "PTC", "LPTC", "ATC"]:
            sub = stg[stg.stage == s].head(5)
            if len(sub):
                print(f"\n[{s}]")
                print(sub[["pathway", "n_slides", "total_send_mean", "mean_per_spot_mean"]].to_string(index=False))

    if rows_dm1:
        df_dm1 = pd.DataFrame(rows_dm1)
        df_dm1.to_csv(RES / "spark_commot_DM1_quartile_diff.tsv", sep="\t", index=False)
        # per-stage mean
        stg_dm1 = df_dm1.groupby(["stage", "pathway"]).agg(
            n_slides=("sample_id", "count"),
            top_minus_bot=("DM1_top_minus_bot", "mean"),
        ).reset_index().sort_values(["stage", "top_minus_bot"], ascending=[True, False])
        print("\n=== top 5 pathways enriched in DM1-top vs DM1-bot per stage ===")
        for s in ["PT", "PTC", "LPTC", "ATC"]:
            sub = stg_dm1[stg_dm1.stage == s].head(5)
            if len(sub):
                print(f"\n[{s}]")
                print(sub.to_string(index=False))


if __name__ == "__main__":
    main()
