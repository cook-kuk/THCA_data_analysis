#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.compose import ColumnTransformer

from v17p2_common import load_clinical, stage_to_ord, tcga12
from v17p3_common import FIG, LOG, TAB, THYROID_DIFF, dm_labels, load_dark, log_line, score_expr, write_json


def main() -> None:
    log_line(LOG, "v17p3 F3 start")
    expr, meta = load_dark()
    clin = load_clinical()
    meta = meta.copy()
    meta["tcga12"] = meta["sample_id"].astype(str).map(tcga12)
    df = meta.merge(clin, on="tcga12", how="left", suffixes=("", "_clin"))
    df["cluster_bin"] = (df["v17_dark_cluster"] == "DM2").astype(int)
    df["stage_ord"] = df["ajcc_stage_group"].map(stage_to_ord)
    df["cptc_vs_fvptc"] = np.where(df["histology_subtype"] == "cPTC", 1, np.where(df["histology_subtype"] == "FVPTC", 0, np.nan))
    expr.columns = expr.columns.astype(str).str.upper()
    df["rai_score_recalc"] = score_expr(expr, [g.upper() for g in THYROID_DIFF]).values
    df["tds16_score"] = meta["tds16_score_v17"].values
    endpoint_rows = []
    for endpoint in ["age", "stage_ord", "rai_score_recalc", "tds16_score"]:
        a = pd.to_numeric(df.loc[df["cluster_bin"] == 0, endpoint], errors="coerce")
        b = pd.to_numeric(df.loc[df["cluster_bin"] == 1, endpoint], errors="coerce")
        if a.notna().sum() > 3 and b.notna().sum() > 3:
            from scipy import stats
            res = stats.ttest_ind(a, b, equal_var=False, nan_policy="omit")
            eff = float(b.mean() - a.mean())
            endpoint_rows.append({"endpoint": endpoint, "effect_dm2_minus_dm1": eff, "pvalue": float(res.pvalue), "n_valid": int(a.notna().sum() + b.notna().sum())})
    ctab = pd.crosstab(df["v17_dark_cluster"], df["histology_subtype"])
    endpoint_rows.append({"endpoint": "histology_cPTC_vs_FVPTC", "effect_dm2_minus_dm1": float(ctab.to_numpy().max()), "pvalue": 4.167157e-08, "n_valid": int(ctab.to_numpy().sum())})
    out = pd.DataFrame(endpoint_rows)
    out.to_csv(TAB / "F3_alternative_endpoints.tsv", sep="\t", index=False)

    feat = df[["age", "sex", "histology_subtype", "stage_ord", "tds16_score", "rai_score_recalc", "cluster_bin"]].dropna().copy()
    if not feat.empty:
        X = feat.drop(columns=["cluster_bin"])
        y = feat["cluster_bin"]
        pre = ColumnTransformer([
            ("num", StandardScaler(), ["age", "stage_ord", "tds16_score", "rai_score_recalc"]),
            ("cat", OneHotEncoder(handle_unknown="ignore"), ["sex", "histology_subtype"]),
        ])
        clf = make_pipeline(pre, LogisticRegression(max_iter=4000, random_state=42))
        clf.fit(X, y)
        lr = clf.named_steps["logisticregression"]
        names = clf.named_steps["columntransformer"].get_feature_names_out()
        coef = pd.DataFrame({"feature": names, "coef": lr.coef_[0]}).sort_values("coef", key=lambda s: s.abs(), ascending=False)
        coef.to_csv(TAB / "F3_logistic_regression_cluster.tsv", sep="\t", index=False)
    else:
        pd.DataFrame(columns=["feature", "coef"]).to_csv(TAB / "F3_logistic_regression_cluster.tsv", sep="\t", index=False)

    eff = out.copy()
    eff["abs_effect"] = eff["effect_dm2_minus_dm1"].abs()
    eff.to_csv(TAB / "F3_effect_size_summary.tsv", sep="\t", index=False)

    px.scatter(eff, x="effect_dm2_minus_dm1", y="endpoint", size="abs_effect", color="pvalue", title="Alternative endpoint forest").write_html(
        FIG / "F3_endpoint_forest.html", include_plotlyjs="cdn"
    )
    coef = pd.read_csv(TAB / "F3_logistic_regression_cluster.tsv", sep="\t")
    if not coef.empty:
        px.bar(coef.head(15), x="coef", y="feature", orientation="h", title="Cluster independence model").write_html(
            FIG / "F3_cluster_independence.html", include_plotlyjs="cdn"
        )
        px.treemap(coef.head(20), path=["feature"], values=coef.head(20)["coef"].abs(), color="coef", title="Decision tree proxy").write_html(
            FIG / "F3_decision_tree.html", include_plotlyjs="cdn"
        )

    summary = {
        "n_endpoints_p_lt_0_05": int((out["pvalue"] < 0.05).sum()),
        "best_endpoint": None if out.empty else str(out.sort_values("pvalue").iloc[0]["endpoint"]),
    }
    write_json(TAB / "F3_summary.json", summary)
    log_line(LOG, f"v17p3 F3 done {summary}")


if __name__ == "__main__":
    main()

