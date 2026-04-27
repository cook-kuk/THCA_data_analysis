#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
from sklearn.cluster import KMeans

from v17p3_common import FIG, LOG, TAB, THYROID_DIFF, V17TAB, read_cross_cancer, log_line, write_json


def main() -> None:
    log_line(LOG, "v17p3 A4 start")
    markers = pd.read_csv(V17TAB / "dark_matter_cluster_markers.tsv", sep="\t")
    dm1 = markers[markers["cluster"] == 0]["gene"].astype(str).str.upper().head(20).tolist()
    dm2 = markers[markers["cluster"] == 1]["gene"].astype(str).str.upper().head(20).tolist()
    rows = []
    gene_rows = []
    for cancer in ["THCA", "LUAD", "COAD", "LGG", "SKCM"]:
        try:
            X, y, b, genes = read_cross_cancer(cancer)
        except Exception:
            continue
        genes_up = [g.upper() for g in genes]
        df = pd.DataFrame(X, columns=genes_up)
        keep1 = [g for g in dm1 if g in df.columns]
        keep2 = [g for g in dm2 if g in df.columns]
        if len(keep1) < 5 or len(keep2) < 5:
            continue
        df["DM1_score"] = df[keep1].mean(axis=1)
        df["DM2_score"] = df[keep2].mean(axis=1)
        km = KMeans(n_clusters=2, random_state=42, n_init=25).fit(df[["DM1_score", "DM2_score"]])
        grp = pd.DataFrame({"cluster": km.labels_, "DM1_score": df["DM1_score"], "DM2_score": df["DM2_score"]})
        rows.append({
            "cancer": cancer,
            "n": int(len(df)),
            "dm1_overlap": int(len(keep1)),
            "dm2_overlap": int(len(keep2)),
            "dm_score_corr": float(np.corrcoef(df["DM1_score"], df["DM2_score"])[0, 1]),
            "cluster_delta": float(grp.groupby("cluster")["DM2_score"].mean().diff().abs().max()),
        })
        for g in sorted(set(keep1 + keep2)):
            gene_rows.append({"cancer": cancer, "gene": g, "present": 1})
    out = pd.DataFrame(rows)
    out.to_csv(TAB / "A4_pancancer_dm_signature_transfer.tsv", sep="\t", index=False)
    cons = pd.DataFrame(gene_rows).groupby("gene").agg(cancer_count=("cancer", "nunique")).reset_index().sort_values("cancer_count", ascending=False)
    cons.to_csv(TAB / "A4_conserved_genes_across_cancers.tsv", sep="\t", index=False)
    if not out.empty:
        px.bar(out, x="cancer", y="cluster_delta", color="dm_score_corr", title="Pan-cancer DM applicability").write_html(FIG / "A4_pancancer_heatmap.html", include_plotlyjs="cdn")
        px.line_polar(out.assign(score=out["cluster_delta"]), r="score", theta="cancer", line_close=True, title="Conservation radar").write_html(FIG / "A4_conservation_radar.html", include_plotlyjs="cdn")
    summary = {"applicable_cancers": int(len(out)), "max_cluster_delta": float(out["cluster_delta"].max()) if not out.empty else None}
    write_json(TAB / "A4_summary.json", summary)
    log_line(LOG, f"v17p3 A4 done {summary}")


if __name__ == "__main__":
    main()

