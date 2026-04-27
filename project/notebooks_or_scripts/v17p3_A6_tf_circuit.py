#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
from sklearn.metrics import roc_auc_score

from v17_common import align_dataset_expr_meta
from v17p3_common import FIG, LOG, RAI_GENES, TAB, TF_TARGETS, THYROID_DIFF, load_dark, log_line, score_expr, write_json


def main() -> None:
    log_line(LOG, "v17p3 A6 start")
    expr, meta = load_dark()
    expr.columns = expr.columns.astype(str).str.upper()
    tf_rows = []
    tf_overlap = []
    for tf, genes in TF_TARGETS.items():
        keep = [g.upper() for g in genes if g.upper() in expr.columns]
        tf_overlap.append({"tf": tf, "n_targets_present": len(keep), "targets": ",".join(keep)})
        vals = expr[keep].mean(axis=1) if keep else pd.Series(np.nan, index=expr.index)
        tf_rows.append(pd.DataFrame({"sample_id": expr.index, "cluster": meta["v17_dark_cluster"].values, "tf": tf, "activity": vals.values}))
    tf_df = pd.concat(tf_rows, ignore_index=True)
    tf_df.to_csv(TAB / "A6_tf_activity_per_sample.tsv", sep="\t", index=False)
    pd.DataFrame(tf_overlap).to_csv(TAB / "A6_tf_target_overlap.tsv", sep="\t", index=False)
    rai = pd.DataFrame({"sample_id": expr.index, "cluster": meta["v17_dark_cluster"].values, "rai_score": score_expr(expr, [g.upper() for g in RAI_GENES]).values})
    rai.to_csv(TAB / "A6_rai_uptake_score.tsv", sep="\t", index=False)

    g760_expr, g760_meta = align_dataset_expr_meta("GSE76039", genes=[g.upper() for g in THYROID_DIFF], tumor_only=True)
    g760_expr.columns = g760_expr.columns.astype(str).str.upper()
    g760_meta = g760_meta.copy()
    g760_meta["rai_score"] = score_expr(g760_expr, [g.upper() for g in THYROID_DIFF]).values
    g760_meta["is_pdtc"] = (g760_meta["histology_subtype"] == "PDTC").astype(int)
    auc = float(roc_auc_score(g760_meta["is_pdtc"], -g760_meta["rai_score"])) if len(np.unique(g760_meta["is_pdtc"])) > 1 else np.nan
    g760_meta.to_csv(TAB / "A6_pdtc_rai_validation.tsv", sep="\t", index=False)

    px.bar(tf_df.groupby(["tf", "cluster"])["activity"].mean().reset_index(), x="tf", y="activity", color="cluster", barmode="group", title="TF circuit activity").write_html(
        FIG / "A6_tf_circuit_diagram.html", include_plotlyjs="cdn"
    )
    px.violin(rai, x="cluster", y="rai_score", color="cluster", box=True, points="all", title="RAI score distribution").write_html(
        FIG / "A6_rai_score_distribution.html", include_plotlyjs="cdn"
    )
    px.box(g760_meta, x="histology_subtype", y="rai_score", color="histology_subtype", points="all", title="PDTC/ATC RAI validation").write_html(
        FIG / "A6_pdtc_rai_validation.html", include_plotlyjs="cdn"
    )
    summary = {"pdtc_auc_using_low_rai": auc}
    write_json(TAB / "A6_summary.json", summary)
    log_line(LOG, f"v17p3 A6 done {summary}")


if __name__ == "__main__":
    main()

