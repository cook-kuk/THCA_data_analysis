#!/usr/bin/env python3
from __future__ import annotations

import re

import numpy as np
import pandas as pd
import plotly.express as px

from v17p35_common import FIG35, LOG35, ROOT, TAB35, log_line, write_json


def normalize_pathway_name(x: str) -> str:
    s = str(x).strip().upper()
    s = re.sub(r"^HALLMARK[_ :]+", "", s)
    s = s.replace("HALLMARK_", "")
    s = re.sub(r"[^A-Z0-9]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s


def main() -> None:
    log_line(LOG35, "v17p35 FIX4 start")
    proper_path = ROOT / "results" / "v17p3" / "tables" / "F2_gsea_hallmark_proper.tsv"
    proxy_path = ROOT / "results" / "v17p2" / "tables" / "gsea_dm1_vs_dm2.tsv"

    proper = pd.read_csv(proper_path, sep="\t")
    proxy = pd.read_csv(proxy_path, sep="\t")

    proper_name_col = "pathway.1" if "pathway.1" in proper.columns else "pathway"
    proxy_name_col = "Term" if "Term" in proxy.columns else ("pathway" if "pathway" in proxy.columns else proxy.columns[0])
    proper_nes_col = "NES"
    proper_fdr_col = "FDR q-val"
    proxy_nes_col = "NES"
    proxy_fdr_col = "FDR q-val"

    proper_sub = proper[[proper_name_col, proper_nes_col, proper_fdr_col]].copy()
    proxy_sub = proxy[[proxy_name_col, proxy_nes_col, proxy_fdr_col]].copy()
    proper_sub.columns = ["pathway_raw", "NES_proper", "FDR_proper"]
    proxy_sub.columns = ["pathway_raw", "NES_proxy", "FDR_proxy"]

    proper_sub["pathway_norm"] = proper_sub["pathway_raw"].map(normalize_pathway_name)
    proxy_sub["pathway_norm"] = proxy_sub["pathway_raw"].map(normalize_pathway_name)

    proper_sub = proper_sub.drop_duplicates("pathway_norm", keep="first")
    proxy_sub = proxy_sub.drop_duplicates("pathway_norm", keep="first")

    proper_sub = proper_sub.rename(columns={"pathway_raw": "pathway_raw_proper"})
    proxy_sub = proxy_sub.rename(columns={"pathway_raw": "pathway_raw_proxy"})
    merged = proper_sub.merge(proxy_sub, on="pathway_norm", how="outer")
    merged["proper_detected_fdr_lt_0_05"] = merged["FDR_proper"].astype(float) < 0.05
    merged["proxy_detected_fdr_lt_0_05"] = merged["FDR_proxy"].astype(float) < 0.05
    merged["sign_agree"] = np.sign(merged["NES_proper"].fillna(0.0)) == np.sign(merged["NES_proxy"].fillna(0.0))
    merged["detected_both_fdr_lt_0_05"] = merged["proper_detected_fdr_lt_0_05"] & merged["proxy_detected_fdr_lt_0_05"]
    merged = merged.sort_values(["pathway_norm"]).reset_index(drop=True)
    merged.to_csv(TAB35 / "FIX4_proxy_vs_proper_full.tsv", sep="\t", index=False)

    both = merged.dropna(subset=["NES_proper", "NES_proxy"]).copy()
    summary = pd.DataFrame(
        [
            {"metric": "n_proper_rows", "value": int(proper_sub.shape[0])},
            {"metric": "n_proxy_rows", "value": int(proxy_sub.shape[0])},
            {"metric": "n_overlap_rows", "value": int(both.shape[0])},
            {"metric": "sign_agreement_rate", "value": float(both["sign_agree"].mean()) if not both.empty else np.nan},
            {
                "metric": "top10_overlap_fdr_lt_0_05",
                "value": int(
                    merged.sort_values("FDR_proper", na_position="last")
                    .head(10)["proxy_detected_fdr_lt_0_05"]
                    .fillna(False)
                    .sum()
                ),
            },
        ]
    )
    summary.to_csv(TAB35 / "FIX4_agreement_summary.tsv", sep="\t", index=False)

    if not both.empty:
        px.scatter(
            both,
            x="NES_proxy",
            y="NES_proper",
            hover_name="pathway_norm",
            color="sign_agree",
            title="Proxy vs proper Hallmark NES",
        ).write_html(FIG35 / "FIX4_proxy_proper_scatter.html", include_plotlyjs="cdn")

        top = merged.sort_values("FDR_proper", na_position="last").head(10).copy()
        top["pathway_display"] = top["pathway_norm"]
        long = top.melt(
            id_vars=["pathway_display"],
            value_vars=["NES_proper", "NES_proxy"],
            var_name="source",
            value_name="NES",
        )
        px.bar(
            long,
            x="pathway_display",
            y="NES",
            color="source",
            barmode="group",
            title="Top pathway proxy/proper concordance",
        ).write_html(FIG35 / "FIX4_top_pathway_concordance.html", include_plotlyjs="cdn")

    payload = {
        "n_overlap_rows": int(both.shape[0]),
        "sign_agreement_rate": float(both["sign_agree"].mean()) if not both.empty else None,
    }
    write_json(TAB35 / "FIX4_summary.json", payload)
    log_line(LOG35, f"v17p35 FIX4 done {payload}")


if __name__ == "__main__":
    main()
