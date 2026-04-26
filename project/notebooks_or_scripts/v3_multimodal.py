"""v3 STEP 7 — Multimodal (best-effort).

Task A (fusion): no local data -> skipped stub
Task B (SCNA): no local data -> skipped stub
Task C (methylation, GSE97466): within-cohort analysis only — tumor vs normal
       because we have no mutation labels.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, "/opt/thyroid-dash/project/notebooks_or_scripts")
from v3_common import (FIGS, GSE97466_METH, RESULTS_TABLES, load_expr,
                        small_n_warning, write_skipped, write_tsv)


def task_fusion():
    out = RESULTS_TABLES / "v3_fusion_results.json"
    write_skipped(out, "No fusion callset locally available.",
                  {"note": "Would require e.g. STAR-Fusion/Arriba output."})
    return {"task": "A_fusion", "status": "skipped", "out": str(out)}


def task_scna():
    out = RESULTS_TABLES / "v3_scna_results.json"
    write_skipped(out, "No SCNA (copy-number) data locally available.",
                  {"note": "Would require GISTIC2 or ASCAT calls."})
    return {"task": "B_scna", "status": "skipped", "out": str(out)}


def task_methylation():
    # GSE97466 beta top5000 — genes unknown at mutation level, so
    # tumor-vs-normal classifier using provided sample_master.
    meth = pd.read_csv(GSE97466_METH, sep="\t")
    meth = meth.rename(columns={meth.columns[0]: "probe"}).set_index("probe")

    sm = pd.read_csv("/opt/thyroid-dash/project/metadata/sample_master.tsv",
                     sep="\t")
    sm = sm[sm.dataset == "GSE97466"].copy()
    sm = sm[sm["normal_vs_tumor"].isin(["tumor", "normal"])]
    samps = [s for s in sm.sample_id if s in meth.columns]
    sm = sm[sm.sample_id.isin(samps)]
    y = (sm.set_index("sample_id").loc[samps, "normal_vs_tumor"] ==
         "tumor").astype(int).to_numpy()
    X = meth[samps].T.to_numpy(float)

    rows = []
    if len(np.unique(y)) >= 2 and y.sum() >= 5 and (len(y) - y.sum()) >= 5:
        skf = StratifiedKFold(n_splits=min(5, int(y.sum()),
                                           int(len(y) - y.sum())),
                              shuffle=True, random_state=42)
        fold_aucs = []
        for tr, te in skf.split(X, y):
            sc = StandardScaler().fit(X[tr])
            clf = LogisticRegression(penalty="l2", max_iter=800,
                                     solver="liblinear",
                                     class_weight="balanced")
            clf.fit(sc.transform(X[tr]), y[tr])
            p = clf.predict_proba(sc.transform(X[te]))[:, 1]
            try:
                fold_aucs.append(roc_auc_score(y[te], p))
            except Exception:
                pass
        mean_auc = float(np.mean(fold_aucs)) if fold_aucs else np.nan
        rows.append(dict(task="methylation_within_gse97466",
                          label="tumor_vs_normal",
                          n_samples=int(len(y)),
                          n_positive=int(y.sum()),
                          mean_cv_auc=mean_auc,
                          note=("⚠ small n" if len(y) < 20 else "")
                               + "; LODO infeasible: no TCGA methylation"))
    else:
        rows.append(dict(task="methylation_within_gse97466",
                          label="tumor_vs_normal",
                          n_samples=int(len(y)),
                          n_positive=int(y.sum()),
                          mean_cv_auc=np.nan,
                          note="insufficient class balance"))

    df = pd.DataFrame(rows)
    out = RESULTS_TABLES / "v3_methylation_results.tsv"
    write_tsv(df, out)

    fig = go.Figure()
    fig.add_bar(x=df["task"], y=df["mean_cv_auc"],
                marker_color="#5eead4",
                text=[f"{v:.3f}" if not pd.isna(v) else "NA"
                      for v in df["mean_cv_auc"]],
                textposition="outside")
    fig.update_layout(
        title=f"GSE97466 methylation within-cohort CV AUC (n={int(df['n_samples'].iloc[0])})",
        yaxis_range=[0, 1.05], yaxis_title="5-fold CV AUC",
        template="plotly_dark", height=360)
    fig.write_html(FIGS / "v3_multimodal_methylation.html",
                    include_plotlyjs="cdn")
    return {"task": "C_methylation", "out": str(out),
             "note": "LODO impossible: no TCGA methylation locally"}


def run():
    return [task_fusion(), task_scna(), task_methylation()]


def main():
    res = run()
    print(json.dumps(res, indent=2, default=str))


if __name__ == "__main__":
    main()
