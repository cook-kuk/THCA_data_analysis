#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px

from v17_common import align_dataset_expr_meta
from v17p3_common import FIG, LOG, TAB, THYROID_DIFF, TIERA67, V17TAB, dm_labels, fit_dm_classifier, load_dark_tiera, log_line, score_expr, write_json


def main() -> None:
    log_line(LOG, "v17p3 A2 start")
    tcga_dark, dark_meta = load_dark_tiera()
    tcga_dark.columns = tcga_dark.columns.astype(str).str.upper()
    marker_df = pd.read_csv(V17TAB / "dark_matter_cluster_markers.tsv", sep="\t")
    genes = sorted(set([g.upper() for g in TIERA67] + marker_df["gene"].astype(str).str.upper().tolist()))
    clf, use = fit_dm_classifier(tcga_dark, dark_meta, genes)

    expr_all, meta_all = align_dataset_expr_meta("TCGA-THCA", genes=use, tumor_only=True)
    expr_all.columns = expr_all.columns.astype(str).str.upper()
    meta_all = meta_all.copy()
    prob_dm2 = clf.predict_proba(expr_all[use])[:, 1]
    meta_all["prob_dm2"] = prob_dm2
    meta_all["prob_dm1"] = 1.0 - prob_dm2
    meta_all["dm_like"] = np.where(prob_dm2 >= 0.5, "DM2_like", "DM1_like")
    meta_all["rai_score_recalc"] = score_expr(expr_all, [g.upper() for g in THYROID_DIFF]).values
    meta_all.to_csv(TAB / "A2_dm_score_full_cohort.tsv", sep="\t", index=False)

    braf = meta_all[meta_all["driver_anchor"] == "BRAF"].copy()
    ras = meta_all[meta_all["driver_anchor"] == "RAS"].copy()
    braf_int = braf.groupby("dm_like").agg(n=("sample_id", "size"), age_mean=("age", "mean"), tds_mean=("tds_score", "mean"), rai_mean=("rai_score_recalc", "mean")).reset_index()
    ras_int = ras.groupby("dm_like").agg(n=("sample_id", "size"), age_mean=("age", "mean"), tds_mean=("tds_score", "mean"), rai_mean=("rai_score_recalc", "mean")).reset_index()
    pd.concat([braf_int.assign(driver="BRAF"), ras_int.assign(driver="RAS")], ignore_index=True).to_csv(TAB / "A2_braf_internal_stratification.tsv", sep="\t", index=False)

    meta_all["driver_dm_group"] = meta_all["driver_anchor"].astype(str) + "/" + meta_all["dm_like"].astype(str)
    keep = meta_all[meta_all["driver_anchor"].isin(["BRAF", "RAS"])].copy()
    grp = keep.groupby("driver_dm_group").agg(n=("sample_id", "size"), age_mean=("age", "mean"), tds_mean=("tds_score", "mean"), rai_mean=("rai_score_recalc", "mean")).reset_index()
    grp.to_csv(TAB / "A2_4group_clinical.tsv", sep="\t", index=False)

    px.violin(keep, x="driver_anchor", y="prob_dm2", color="dm_like", box=True, points="all", title="DM score by driver").write_html(
        FIG / "A2_dm_score_by_driver.html", include_plotlyjs="cdn"
    )
    px.scatter(keep, x="tds_score", y="prob_dm2", color="driver_anchor", symbol="dm_like", title="BRAF/RAS internal DM axis").write_html(
        FIG / "A2_braf_internal_umap.html", include_plotlyjs="cdn"
    )
    px.bar(grp, x="driver_dm_group", y="rai_mean", color="driver_dm_group", title="4-group outcomes proxy").write_html(
        FIG / "A2_4group_outcomes.html", include_plotlyjs="cdn"
    )

    corr = float(np.corrcoef((meta_all["driver_anchor"] == "BRAF").astype(int), meta_all["prob_dm2"])[0, 1])
    summary = {"driver_dm_correlation_braf_vs_probdm2": corr, "braf_dm2_like_n": int((braf["dm_like"] == "DM2_like").sum()), "ras_dm1_like_n": int((ras["dm_like"] == "DM1_like").sum())}
    write_json(TAB / "A2_summary.json", summary)
    log_line(LOG, f"v17p3 A2 done {summary}")


if __name__ == "__main__":
    main()

