#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from v17p35_common import FIG35, LOG35, ROOT, TAB35, dm_labels, load_dark_tiera, log_line, write_json


def main() -> None:
    log_line(LOG35, "v17p35 FIX1 start")
    expr, meta = load_dark_tiera()
    expr.columns = expr.columns.astype(str).str.upper()
    y = dm_labels(meta)
    clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=4000, random_state=42))
    clf.fit(expr, y)

    ccle = pd.read_csv(ROOT / "results" / "v14_ccle" / "ccle_thyroid_expression.tsv", sep="\t").rename(columns={"entrezGeneId": "symbol"})
    ccle["symbol"] = ccle["symbol"].astype(str).str.upper()
    ccle = ccle.drop_duplicates("symbol").set_index("symbol")
    X = ccle.reindex(expr.columns).fillna(0.0).T
    probs = clf.predict_proba(X)[:, 1]
    scores = pd.DataFrame({"sample": X.index, "prob_dm2": probs, "prob_dm1": 1.0 - probs})
    scores["dm_like"] = np.where(scores["prob_dm2"] >= 0.5, "DM2_like", "DM1_like")
    val = pd.read_csv(ROOT / "results" / "v14_ccle" / "ccle_brs52_validation.tsv", sep="\t")
    scores = scores.merge(val[["sample", "mutation_label", "v14_prediction", "name"]], on="sample", how="left")
    scores.to_csv(TAB35 / "FIX1_celline_dm_scores.tsv", sep="\t", index=False)

    drug = pd.read_csv(ROOT / "results" / "v14_ccle" / "ccle_drug_sensitivity.tsv", sep="\t")
    dm_map = scores.groupby("v14_prediction")["prob_dm2"].mean().to_dict()
    drug["prob_dm2"] = drug["v14_prediction"].map(dm_map)
    rows = []
    for compound, sub in drug.groupby("compound"):
        ok = sub["prob_dm2"].notna() & sub["logfold_change"].notna()
        if ok.sum() < 4:
            continue
        rho, p = stats.spearmanr(sub.loc[ok, "prob_dm2"], sub.loc[ok, "logfold_change"])
        rows.append(
            {
                "compound": compound,
                "n_points": int(ok.sum()),
                "spearman_rho": float(rho),
                "pvalue": float(p),
                "mean_logfold_change": float(sub.loc[ok, "logfold_change"].mean()),
                "moa": str(sub["moa"].dropna().iloc[0]) if sub["moa"].notna().any() else np.nan,
            }
        )
    out = pd.DataFrame(rows).sort_values("pvalue")
    if not out.empty:
        from v17_common import bh_fdr
        out["fdr"] = bh_fdr(out["pvalue"].to_numpy(float))
    else:
        out["fdr"] = []
    out.sort_values("spearman_rho").head(50).to_csv(TAB35 / "FIX1_top_drugs_dm1_selective.tsv", sep="\t", index=False)
    out.sort_values("spearman_rho", ascending=False).head(50).to_csv(TAB35 / "FIX1_top_drugs_dm2_selective.tsv", sep="\t", index=False)

    if not out.empty:
        moa = out.groupby("moa").agg(n=("compound", "size"), mean_rho=("spearman_rho", "mean"), min_fdr=("fdr", "min")).reset_index().sort_values(["min_fdr", "n"])
    else:
        moa = pd.DataFrame(columns=["moa", "n", "mean_rho", "min_fdr"])
    moa.to_csv(TAB35 / "FIX1_moa_enrichment.tsv", sep="\t", index=False)

    sm = pd.read_csv(ROOT / "results" / "v17p3" / "tables" / "A2_dm_score_full_cohort.tsv", sep="\t")
    tumor_pred = sm[["sample_id", "driver_anchor", "prob_dm2", "prob_dm1", "dm_like"]].copy()
    if not out.empty:
        top_dm1 = out.sort_values("spearman_rho").head(3)["compound"].tolist()
        top_dm2 = out.sort_values("spearman_rho", ascending=False).head(3)["compound"].tolist()
        tumor_pred["top3_predicted_drugs"] = np.where(tumor_pred["prob_dm2"] < 0.5, ",".join(top_dm1), ",".join(top_dm2))
        tumor_pred["response_score"] = np.where(tumor_pred["prob_dm2"] < 0.5, 1.0 - tumor_pred["prob_dm2"], tumor_pred["prob_dm2"])
    tumor_pred.to_csv(TAB35 / "FIX1_tumor_predicted_response.tsv", sep="\t", index=False)

    if not out.empty:
        volc = out.copy()
        volc["mlog10p"] = -np.log10(np.clip(volc["pvalue"], 1e-300, 1.0))
        px.scatter(volc, x="spearman_rho", y="mlog10p", hover_name="compound", color="fdr", title="Drug selectivity by DM2 probability correlation").write_html(
            FIG35 / "FIX1_drug_volcano.html", include_plotlyjs="cdn"
        )
        top = pd.concat([out.nsmallest(15, "spearman_rho"), out.nlargest(15, "spearman_rho")], ignore_index=True)
        px.bar(top, x="compound", y="spearman_rho", color="fdr", title="Top drugs by DM correlation").write_html(
            FIG35 / "FIX1_top_drugs_heatmap.html", include_plotlyjs="cdn"
        )
        if not moa.empty:
            px.line_polar(moa.head(12).assign(abs_rho=lambda d: d["mean_rho"].abs()), r="abs_rho", theta="moa", line_close=True, title="MOA correlation radar").write_html(
                FIG35 / "FIX1_moa_pathway_radar.html", include_plotlyjs="cdn"
            )
        px.histogram(tumor_pred, x="response_score", color="dm_like", nbins=40, title="Tumor predicted response score distribution").write_html(
            FIG35 / "FIX1_tumor_prediction_dist.html", include_plotlyjs="cdn"
        )

    payload = {
        "n_drugs_fdr_lt_0_1": int((out["fdr"] < 0.1).sum()) if not out.empty else 0,
        "n_cell_lines": int(scores.shape[0]),
        "n_dm2_like_cell_lines": int((scores["dm_like"] == "DM2_like").sum()),
    }
    write_json(TAB35 / "FIX1_summary.json", payload)
    log_line(LOG35, f"v17p35 FIX1 done {payload}")


if __name__ == "__main__":
    main()
