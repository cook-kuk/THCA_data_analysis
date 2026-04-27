#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
import scanpy as sc

from v17p35_common import FIG35, LOG35, P2TAB, SCRNA, TAB35, log_line, write_json

IMMUNE_MARKERS = {
    "T_cell": ["CD3D", "CD3E", "TRBC1", "IL7R"],
    "B_cell": ["MS4A1", "CD79A", "CD79B", "CD74"],
    "Macrophage": ["LST1", "TYROBP", "FCER1G", "CTSB"],
    "NK_cell": ["NKG7", "GNLY", "KLRD1", "TRAC"],
    "Epithelial": ["EPCAM", "KRT18", "KRT19", "KRT8"],
    "Stromal": ["COL1A1", "COL1A2", "DCN", "LUM"],
}


def main() -> None:
    log_line(LOG35, "v17p35 FIX5 start")
    scdf = pd.read_csv(P2TAB / "scrna_cell_dm_score.tsv", sep="\t")
    ad = sc.read_h5ad(SCRNA)
    ad.var["gene_symbol"] = ad.var["gene_symbol"].astype(str).str.upper()

    score_df = pd.DataFrame(index=ad.obs_names)
    for ctype, genes in IMMUNE_MARKERS.items():
        idx = [i for i, g in enumerate(ad.var["gene_symbol"]) if g in genes]
        vals = np.asarray(ad[:, idx].X.mean(axis=1)).ravel() if idx else np.zeros(ad.n_obs)
        score_df[ctype] = vals
    ann = pd.DataFrame(index=ad.obs_names)
    ann["cell_type"] = score_df.idxmax(axis=1)
    ann["patient_id_ann"] = ad.obs["patient_id"].astype(str).values
    ann["gsm"] = ad.obs["gsm"].astype(str).values
    ann["sample"] = ad.obs["sample"].astype(str).values
    ann["mutation_status"] = ad.obs["mutation_status"].astype(str).values

    cell = scdf.set_index("cell_id").join(ann, how="inner").reset_index().rename(columns={"index": "cell_id"})
    if "patient_id" not in cell.columns and "patient_id_ann" in cell.columns:
        cell["patient_id"] = cell["patient_id_ann"]
    per_pt = cell.groupby("patient_id").agg(
        n_cells=("cell_id", "size"),
        dm1_dominant_cells=("dominant_dm", lambda s: int((s == "DM1").sum())),
        dm2_dominant_cells=("dominant_dm", lambda s: int((s == "DM2").sum())),
        mutation_status=("mutation_status", lambda s: str(pd.Series(s).mode().iloc[0]) if len(pd.Series(s).mode()) else "NA"),
    ).reset_index()
    per_pt["dominance_ratio"] = per_pt[["dm1_dominant_cells", "dm2_dominant_cells"]].max(axis=1) / per_pt["n_cells"]
    per_pt["dominant_phenotype"] = np.where(per_pt["dm1_dominant_cells"] >= per_pt["dm2_dominant_cells"], "DM1", "DM2")
    per_pt.to_csv(TAB35 / "FIX5_per_patient_full.tsv", sep="\t", index=False)

    immune = cell[cell["cell_type"].isin(["T_cell", "B_cell", "Macrophage", "NK_cell"])].copy()
    breakdown = immune.groupby(["dominant_dm", "cell_type"]).size().reset_index(name="n")
    breakdown["fraction_within_dm"] = breakdown.groupby("dominant_dm")["n"].transform(lambda s: s / s.sum())
    breakdown.to_csv(TAB35 / "FIX5_immune_cell_type_breakdown.tsv", sep="\t", index=False)

    per_cell = cell[["cell_id", "patient_id", "dominant_dm", "DM1_score", "DM2_score", "cell_type", "mutation_status"]].copy()
    per_cell.to_csv(TAB35 / "FIX5_per_cell_dm_signature.tsv", sep="\t", index=False)

    px.bar(breakdown, x="cell_type", y="fraction_within_dm", color="dominant_dm", barmode="group", title="Immune cell composition by DM state").write_html(
        FIG35 / "FIX5_immune_composition_stacked.html", include_plotlyjs="cdn"
    )
    px.bar(per_pt, x="patient_id", y="dominance_ratio", color="dominant_phenotype", hover_data=["mutation_status", "n_cells"], title="Per-patient DM landscape").write_html(
        FIG35 / "FIX5_per_patient_dm_landscape.html", include_plotlyjs="cdn"
    )
    plot_df = per_cell.sample(min(30000, len(per_cell)), random_state=42)
    px.scatter(plot_df, x="DM1_score", y="DM2_score", color="cell_type", hover_data=["patient_id", "dominant_dm", "mutation_status"], title="Single-cell DM score UMAP proxy").write_html(
        FIG35 / "FIX5_single_cell_dm_score_umap.html", include_plotlyjs="cdn"
    )

    payload = {
        "patients": int(per_pt.shape[0]),
        "immune_rows": int(breakdown.shape[0]),
        "mutation_status_values": sorted(per_pt["mutation_status"].dropna().unique().tolist()),
    }
    write_json(TAB35 / "FIX5_summary.json", payload)
    log_line(LOG35, f"v17p35 FIX5 done {payload}")


if __name__ == "__main__":
    main()
