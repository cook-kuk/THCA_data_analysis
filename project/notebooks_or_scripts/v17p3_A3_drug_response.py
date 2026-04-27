#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
import requests

from v17p3_common import FIG, LOG, TAB, TIERA67, log_line, write_json


def map_entrez_to_symbol(ids: list[str]) -> pd.DataFrame:
    rows = []
    chunk = 800
    for i in range(0, len(ids), chunk):
        sub = ids[i:i + chunk]
        r = requests.post("https://mygene.info/v3/query", data={"q": ",".join(sub), "scopes": "entrezgene", "fields": "symbol", "species": "human"}, timeout=120)
        r.raise_for_status()
        rows.extend(r.json())
    out = []
    for rec in rows:
        q = str(rec.get("query", ""))
        sym = rec.get("symbol")
        if q and sym:
            out.append({"entrezGeneId": q, "symbol": str(sym).upper()})
    return pd.DataFrame(out, columns=["entrezGeneId", "symbol"]).drop_duplicates()


def main() -> None:
    log_line(LOG, "v17p3 A3 start")
    expr = pd.read_csv("project/results/v14_ccle/ccle_thyroid_expression.tsv", sep="\t")
    expr["entrezGeneId"] = expr["entrezGeneId"].astype(str)
    looks_like_symbol = expr["entrezGeneId"].str.fullmatch(r"[A-Za-z0-9_.-]+").mean() > 0.95 and not expr["entrezGeneId"].str.fullmatch(r"\d+").mean() > 0.5
    if looks_like_symbol:
        ex = expr.rename(columns={"entrezGeneId": "symbol"}).copy()
        ex["symbol"] = ex["symbol"].str.upper()
        ex = ex.drop_duplicates("symbol").set_index("symbol")
    else:
        mapping = map_entrez_to_symbol(expr["entrezGeneId"].tolist())
        ex = expr.merge(mapping, on="entrezGeneId", how="left").dropna(subset=["symbol"]).drop_duplicates("symbol").set_index("symbol")
        ex = ex.drop(columns=["entrezGeneId"])
    val = pd.read_csv("project/results/v14_ccle/ccle_brs52_validation.tsv", sep="\t")
    samples = [c for c in ex.columns if c in set(val["sample"])]
    ex = ex[samples].T
    dm1_mark = ["DUSP5", "DUSP6", "DUSP4", "MET", "FOXP3"]
    dm2_mark = ["TPO", "DIO1", "DIO2", "SLC5A8", "TG", "TSHR"]
    score = pd.DataFrame({"sample": ex.index, "DM1_score": ex[[g for g in dm1_mark if g in ex.columns]].mean(axis=1), "DM2_score": ex[[g for g in dm2_mark if g in ex.columns]].mean(axis=1)})
    score["dm_like"] = np.where(score["DM2_score"] > score["DM1_score"], "DM2_like", "DM1_like")
    score = score.merge(val[["sample", "v14_prediction", "mutation_label"]], on="sample", how="left")
    score.to_csv(TAB / "A3_celline_dm_scores.tsv", sep="\t", index=False)

    drug = pd.read_csv("project/results/v14_ccle/ccle_drug_sensitivity.tsv", sep="\t")
    label_map = score.groupby("v14_prediction")[["DM1_score", "DM2_score"]].mean().reset_index()
    label_map["dm_like"] = np.where(label_map["DM2_score"] > label_map["DM1_score"], "DM2_like", "DM1_like")
    drug = drug.merge(label_map, on="v14_prediction", how="left")
    rows = []
    from scipy import stats
    for compound, sub in drug.groupby("compound"):
        a = sub.loc[sub["dm_like"] == "DM1_like", "logfold_change"].astype(float)
        b = sub.loc[sub["dm_like"] == "DM2_like", "logfold_change"].astype(float)
        if a.notna().sum() >= 2 and b.notna().sum() >= 2:
            p = float(stats.ttest_ind(a, b, equal_var=False, nan_policy="omit").pvalue)
            rows.append({"compound": compound, "dm1_mean_lfc": float(a.mean()), "dm2_mean_lfc": float(b.mean()), "delta_dm1_minus_dm2": float(a.mean() - b.mean()), "pvalue": p})
    out = pd.DataFrame(rows, columns=["compound", "dm1_mean_lfc", "dm2_mean_lfc", "delta_dm1_minus_dm2", "pvalue"])
    out["fdr"] = np.nan
    if not out.empty:
        from v17_common import bh_fdr
        out["fdr"] = bh_fdr(out["pvalue"].to_numpy(float))
    out.sort_values("delta_dm1_minus_dm2").head(20).to_csv(TAB / "A3_top_drugs_dm1_selective.tsv", sep="\t", index=False)
    out.sort_values("delta_dm1_minus_dm2", ascending=False).head(20).to_csv(TAB / "A3_top_drugs_dm2_selective.tsv", sep="\t", index=False)
    coherence = out[out["compound"].str.contains("topotecan|irinotecan|simvastatin|atorvastatin|trametinib|dabrafenib|selumetinib", case=False, na=False)].copy()
    coherence.to_csv(TAB / "A3_pathway_drug_coherence.tsv", sep="\t", index=False)
    pd.DataFrame().to_csv(TAB / "A3_tumor_predicted_response.tsv", sep="\t", index=False)

    if not out.empty:
        volc = out.copy()
        volc["mlog10fdr"] = -np.log10(np.clip(volc["fdr"].fillna(1.0), 1e-300, 1.0))
        px.scatter(volc, x="delta_dm1_minus_dm2", y="mlog10fdr", hover_name="compound", title="Drug selectivity volcano").write_html(
            FIG / "A3_drug_volcano_dm1_vs_dm2.html", include_plotlyjs="cdn"
        )
        top = pd.concat([out.nsmallest(10, "delta_dm1_minus_dm2"), out.nlargest(10, "delta_dm1_minus_dm2")], ignore_index=True)
        px.bar(top, x="compound", y="delta_dm1_minus_dm2", color="fdr", title="Top selective drugs").write_html(
            FIG / "A3_top_drugs_heatmap.html", include_plotlyjs="cdn"
        )
        if not coherence.empty:
            px.line_polar(coherence.assign(abs_delta=coherence["delta_dm1_minus_dm2"].abs()), r="abs_delta", theta="compound", line_close=True, title="Pathway-drug coherence").write_html(
                FIG / "A3_pathway_drug_radar.html", include_plotlyjs="cdn"
            )
    summary = {"n_drugs_fdr_lt_0_1": int((out["fdr"] < 0.1).sum()) if not out.empty else 0, "best_dm1_drug": None if out.empty else str(out.sort_values("delta_dm1_minus_dm2").iloc[0]["compound"])}
    write_json(TAB / "A3_summary.json", summary)
    log_line(LOG, f"v17p3 A3 done {summary}")


if __name__ == "__main__":
    main()
