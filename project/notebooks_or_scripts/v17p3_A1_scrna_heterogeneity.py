#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
import scanpy as sc

from v17p3_common import FIG, LOG, P2TAB, SCRNA, TAB, log_line, write_json

CELLTYPE_MARKERS = {
    "epithelial": ["EPCAM", "KRT18", "KRT19"],
    "immune": ["PTPRC", "LST1", "TYROBP"],
    "stromal": ["COL1A1", "COL1A2", "DCN"],
}


def main() -> None:
    log_line(LOG, "v17p3 A1 start")
    scdf = pd.read_csv(P2TAB / "scrna_cell_dm_score.tsv", sep="\t")
    per_pt = scdf.groupby("patient_id").agg(
        n_cells=("cell_id", "size"),
        dm1_dominance=("dominant_dm", lambda s: float((s == "DM1").mean())),
        dm2_dominance=("dominant_dm", lambda s: float((s == "DM2").mean())),
    ).reset_index()
    per_pt["dominance_ratio"] = per_pt[["dm1_dominance", "dm2_dominance"]].max(axis=1)
    per_pt.to_csv(TAB / "A1_per_patient_dominance.tsv", sep="\t", index=False)

    ad = sc.read_h5ad(SCRNA)
    ad.var["gene_symbol"] = ad.var["gene_symbol"].astype(str).str.upper()
    ctype = pd.DataFrame(index=ad.obs_names)
    for name, genes in CELLTYPE_MARKERS.items():
        idx = [i for i, g in enumerate(ad.var["gene_symbol"]) if g in genes]
        vals = np.asarray(ad[:, idx].X.mean(axis=1)).ravel() if idx else np.zeros(ad.n_obs)
        ctype[name] = vals
    celltype = ctype.idxmax(axis=1)
    ann = pd.DataFrame({"patient_id_ann": ad.obs["patient_id"].astype(str).values, "cell_type": celltype.values}, index=ad.obs_names)
    cell = scdf.set_index("cell_id").join(ann, how="inner").reset_index()
    if "patient_id" not in cell.columns:
        cell["patient_id"] = cell["patient_id_ann"]
    ctab = cell.groupby(["cell_type", "dominant_dm"]).size().reset_index(name="n")
    ctab.to_csv(TAB / "A1_cell_type_per_cluster.tsv", sep="\t", index=False)
    dom = per_pt.copy()
    dom["dominant_label"] = np.where(dom["dm2_dominance"] > dom["dm1_dominance"], "DM2", "DM1")
    dom.assign(rank=np.arange(len(dom)) + 1).to_csv(TAB / "A1_dominance_vs_patient_features.tsv", sep="\t", index=False)

    px.bar(dom, x="patient_id", y="dominance_ratio", color="dominant_label", title="DM dominance per patient").write_html(
        FIG / "A1_dominance_per_patient_bar.html", include_plotlyjs="cdn"
    )
    px.scatter(dom, x="n_cells", y="dominance_ratio", color="dominant_label", hover_name="patient_id", title="Dominance vs patient cell count").write_html(
        FIG / "A1_dominance_clinical_correlation.html", include_plotlyjs="cdn"
    )
    px.bar(ctab, x="cell_type", y="n", color="dominant_dm", barmode="stack", title="Cell type by DM state").write_html(
        FIG / "A1_scrna_celltype_dm_overlap.html", include_plotlyjs="cdn"
    )
    plot_df = scdf.sample(min(20000, len(scdf)), random_state=42)
    px.scatter(plot_df, x="DM1_score", y="DM2_score", color="dominant_dm", hover_data=["patient_id"], title="Single-cell DM trajectory proxy").write_html(
        FIG / "A1_single_cell_trajectory.html", include_plotlyjs="cdn"
    )

    summary = {
        "patients": int(len(dom)),
        "mixed_patients": int((dom["dominance_ratio"] < 0.99).sum()),
        "dominance_ratio_range": [float(dom["dominance_ratio"].min()), float(dom["dominance_ratio"].max())],
    }
    write_json(TAB / "A1_summary.json", summary)
    log_line(LOG, f"v17p3 A1 done {summary}")


if __name__ == "__main__":
    main()
