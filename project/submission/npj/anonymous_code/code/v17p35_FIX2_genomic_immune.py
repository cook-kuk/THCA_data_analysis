#!/usr/bin/env python3
from __future__ import annotations

import tarfile

import numpy as np
import pandas as pd
import plotly.express as px
from scipy import stats

from v17_common import ROOT
from v17p2_common import tcga12
from v17p35_common import FIG35, LOG35, TAB35, load_dark, log_line, write_json


GISTIC = ROOT / "data_raw" / "v3_ext" / "tcga_scna" / "gistic2.tar.gz"
PREFIX = "gdac.broadinstitute.org_THCA-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/"


def read_member(member: str) -> pd.DataFrame:
    with tarfile.open(GISTIC, "r:gz") as tf:
        f = tf.extractfile(PREFIX + member)
        if f is None:
            raise FileNotFoundError(member)
        return pd.read_csv(f, sep="\t")


def main() -> None:
    log_line(LOG35, "v17p35 FIX2 start")
    _, meta = load_dark()
    meta = meta.copy()
    meta["tcga12"] = meta["sample_id"].astype(str).map(tcga12)

    gene_df = read_member("all_thresholded.by_genes.txt")
    gene_df = gene_df.rename(columns={"Gene Symbol": "gene_symbol"})
    sample_cols = [c for c in gene_df.columns if c.startswith("TCGA-")]
    gene_num = gene_df[sample_cols].apply(pd.to_numeric, errors="coerce")
    gene_summary = pd.DataFrame(
        {
            "sample_barcode": sample_cols,
            "cnv_nonzero_frac": (gene_num.ne(0).sum(axis=0) / gene_num.shape[0]).values,
            "amp_frac": (gene_num.gt(0).sum(axis=0) / gene_num.shape[0]).values,
            "del_frac": (gene_num.lt(0).sum(axis=0) / gene_num.shape[0]).values,
            "cnv_abs_mean": gene_num.abs().mean(axis=0).values,
        }
    )
    gene_summary["tcga12"] = gene_summary["sample_barcode"].map(tcga12)

    arm_df = read_member("broad_values_by_arm.txt")
    arm_sample_cols = [c for c in arm_df.columns if c.startswith("TCGA-")]
    arm_num = arm_df[arm_sample_cols].apply(pd.to_numeric, errors="coerce")
    arm_summary = pd.DataFrame(
        {
            "sample_barcode": arm_sample_cols,
            "arm_abs_mean": arm_num.abs().mean(axis=0).values,
            "arm_gain_frac": (arm_num.gt(0.2).sum(axis=0) / arm_num.shape[0]).values,
            "arm_loss_frac": (arm_num.lt(-0.2).sum(axis=0) / arm_num.shape[0]).values,
            "wgd_proxy": (arm_num.abs().gt(0.2).sum(axis=0) >= (arm_num.shape[0] * 0.5)).astype(int).values,
        }
    )
    arm_summary["tcga12"] = arm_summary["sample_barcode"].map(tcga12)

    merged = meta.merge(gene_summary.drop(columns=["sample_barcode"]), on="tcga12", how="left").merge(
        arm_summary.drop(columns=["sample_barcode"]), on="tcga12", how="left"
    )
    merged.to_csv(TAB35 / "FIX2_aneuploidy_per_sample.tsv", sep="\t", index=False)

    summ_rows = []
    for metric in ["cnv_nonzero_frac", "amp_frac", "del_frac", "cnv_abs_mean", "arm_abs_mean", "arm_gain_frac", "arm_loss_frac", "wgd_proxy"]:
        a = pd.to_numeric(merged.loc[merged["v17_dark_cluster"] == "DM1", metric], errors="coerce")
        b = pd.to_numeric(merged.loc[merged["v17_dark_cluster"] == "DM2", metric], errors="coerce")
        t = stats.ttest_ind(a, b, equal_var=False, nan_policy="omit")
        summ_rows.append(
            {
                "metric": metric,
                "dm1_mean": float(a.mean()),
                "dm2_mean": float(b.mean()),
                "dm2_minus_dm1": float(b.mean() - a.mean()),
                "pvalue": float(t.pvalue),
            }
        )
    pd.DataFrame(summ_rows).to_csv(TAB35 / "FIX2_genomic_summary_per_dm.tsv", sep="\t", index=False)

    pd.DataFrame([{"sample_id": np.nan, "immune_subtype": "MISSING", "source": "Thorsson subtype table not available locally in this sprint"}]).to_csv(
        TAB35 / "FIX2_thorsson_subtype.tsv", sep="\t", index=False
    )

    px.box(merged, x="v17_dark_cluster", y="arm_abs_mean", color="v17_dark_cluster", points="all", title="Aneuploidy proxy by DM cluster").write_html(
        FIG35 / "FIX2_aneuploidy_violin.html", include_plotlyjs="cdn"
    )
    px.box(merged, x="v17_dark_cluster", y="cnv_nonzero_frac", color="v17_dark_cluster", points="all", title="CNV burden proxy by DM cluster").write_html(
        FIG35 / "FIX2_cnv_genome_view.html", include_plotlyjs="cdn"
    )
    px.bar(pd.DataFrame([{"cluster": "DM1", "status": "MISSING"}, {"cluster": "DM2", "status": "MISSING"}]), x="cluster", y=[1, 1], color="status", title="Thorsson distribution placeholder").write_html(
        FIG35 / "FIX2_thorsson_distribution.html", include_plotlyjs="cdn"
    )

    payload = {
        "cnv_nonzero_frac_dm2_minus_dm1": float(pd.DataFrame(summ_rows).set_index("metric").loc["cnv_nonzero_frac", "dm2_minus_dm1"]),
        "thorsson_status": "MISSING",
    }
    write_json(TAB35 / "FIX2_summary.json", payload)
    log_line(LOG35, f"v17p35 FIX2 done {payload}")


if __name__ == "__main__":
    main()
