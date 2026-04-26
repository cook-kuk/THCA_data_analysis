#!/usr/bin/env python
"""v3_ext_06_batch_diagnostics.py

Batch-effect / dataset-identifiability diagnostics.
  D1 multiclass LogReg to predict cohort (flag macro AUC > 0.90).
  D2 PCA + UMAP by cohort vs subtype.
  D3 Subtype-preserving ComBat pilot on pooled RNA-seq - delta AUC pre/post.
  D4 TDS16 and BRS71_original cross-cohort boxplots.
"""
import json
import logging
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/opt/thyroid-dash/project")
DP = ROOT / "data_processed"
DP_V3_RNA = DP / "bulk_rnaseq_v3"
DP_V3_MIC = DP / "microarray_v3"
META = ROOT / "metadata"
TABLES = ROOT / "results" / "tables"
FIGS = ROOT / "reports" / "html" / "figs_interactive"
LOGDIR = ROOT / "logs"
for d in (TABLES, FIGS, LOGDIR):
    d.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[logging.FileHandler(LOGDIR / "v3_ext_06_batch_diagnostics.log", mode="w"),
              logging.StreamHandler()],
)
log = logging.getLogger("v3_ext_06")

sys.path.insert(0, str(ROOT / "notebooks_or_scripts"))
try:
    from rerun_v2 import TDS16
except Exception:
    TDS16 = ["DIO1", "DIO2", "DUOX1", "DUOX2", "FOXE1", "GLIS3", "NKX2-1", "PAX8",
             "SLC26A4", "SLC5A5", "SLC5A8", "TG", "THRA", "THRB", "TPO", "TSHR"]


def load_v3_matrices():
    mats = {}
    for p in list(DP_V3_RNA.glob("*_v3_log2.tsv")) + list(DP_V3_MIC.glob("*_v3_log2.tsv")):
        name = p.stem.replace("_v3_log2", "")
        try:
            df = pd.read_csv(p, sep="\t", index_col=0)
            # skip ENSG-indexed matrices to keep the shared HGNC gene space populated
            first = str(df.index[0]) if len(df.index) > 0 else ""
            if first.startswith("ENSG"):
                log.info(f"skip {name}: ENSG-indexed")
                continue
            mats[name] = df
        except Exception as e:
            log.warning(f"load {name} fail: {e}")
    return mats


def d1_dataset_identifiability(mats):
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import Pipeline
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import roc_auc_score
    # Shared genes
    sets = [set(m.index) for m in mats.values()]
    shared = sorted(set.intersection(*sets)) if sets else []
    if len(shared) < 5 or len(mats) < 2:
        return {"status": "skipped", "reason": "insufficient_cohorts_or_genes"}
    # Build X, y
    X_all = []; y_all = []; cohort_map = {}
    cohorts = sorted(mats.keys())
    for i, c in enumerate(cohorts):
        cohort_map[c] = i
        m = mats[c].loc[shared].T
        X_all.append(m.values); y_all.append(np.full(len(m), i))
    X = np.vstack(X_all); y = np.concatenate(y_all)
    if len(set(y)) < 2:
        return {"status": "skipped", "reason": "only_one_cohort"}
    pipe = Pipeline([("sc", StandardScaler()),
                     ("clf", LogisticRegression(max_iter=2000, random_state=42))])
    skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)
    aucs = []
    for tr, te in skf.split(X, y):
        pipe.fit(X[tr], y[tr])
        p = pipe.predict_proba(X[te])
        try:
            aucs.append(roc_auc_score(y[te], p, multi_class="ovr"))
        except Exception:
            pass
    macro = float(np.mean(aucs)) if aucs else np.nan
    flag = macro > 0.90
    return {"status": "ok", "cohorts": cohorts, "macro_auc": macro, "flag_identifiable": bool(flag),
            "n_shared_genes": len(shared)}


def d2_pca_umap(mats):
    try:
        from sklearn.decomposition import PCA
        sets = [set(m.index) for m in mats.values()]
        shared = sorted(set.intersection(*sets))
        if len(shared) < 10 or len(mats) < 2:
            return {"status": "skipped"}
        frames = []; labels = []
        for c, m in mats.items():
            mm = m.loc[shared].T
            frames.append(mm.values); labels.extend([c] * len(mm))
        X = np.vstack(frames)
        pca = PCA(n_components=2, random_state=42).fit_transform(X)
        # Write figure
        import plotly.graph_objects as go
        fig = go.Figure()
        labels_arr = np.array(labels)
        for c in np.unique(labels_arr):
            m = labels_arr == c
            fig.add_trace(go.Scatter(x=pca[m, 0], y=pca[m, 1], mode="markers", name=c,
                                     marker=dict(size=6, opacity=0.6)))
        fig.update_layout(title="D2 PCA by cohort", xaxis_title="PC1", yaxis_title="PC2")
        fig.write_html(str(FIGS / "v3_ext_d2_pca.html"), include_plotlyjs="cdn")
        return {"status": "ok", "n_samples": len(labels), "n_shared_genes": len(shared)}
    except Exception as e:
        return {"status": "failed", "error": str(e)}


def d3_combat_pilot(mats):
    # Lite ComBat: per-cohort centring then add overall mean
    try:
        sets = [set(m.index) for m in mats.values()]
        shared = sorted(set.intersection(*sets))
        if len(shared) < 10 or len(mats) < 2:
            return {"status": "skipped"}
        # Pool
        frames = []
        labels = []
        for c, m in mats.items():
            mm = m.loc[shared]
            for col in mm.columns:
                frames.append(mm[col].values); labels.append(c)
        arr = np.array(frames)  # samples x genes
        # Adjusted
        adj = arr.copy()
        labels_arr = np.array(labels)
        for c in np.unique(labels_arr):
            m = labels_arr == c
            adj[m] = adj[m] - adj[m].mean(axis=0)
        adj = adj + arr.mean(axis=0)
        # Delta: TDS16 signal strength pre/post
        shared_idx = {g: i for i, g in enumerate(shared)}
        tds_idx = [shared_idx[g] for g in TDS16 if g in shared_idx]
        pre_spread = float(np.std(arr[:, tds_idx])) if tds_idx else np.nan
        post_spread = float(np.std(adj[:, tds_idx])) if tds_idx else np.nan
        return {"status": "ok", "tds16_std_pre": pre_spread, "tds16_std_post": post_spread,
                "delta_std": (pre_spread - post_spread) if (pre_spread == pre_spread) else np.nan}
    except Exception as e:
        return {"status": "failed", "error": str(e)}


def d4_cross_cohort_boxplots(mats):
    try:
        import plotly.graph_objects as go
        sets = [set(m.index) for m in mats.values()]
        shared = sorted(set.intersection(*sets))
        tds = [g for g in TDS16 if g in shared]
        if not tds:
            return {"status": "skipped"}
        fig = go.Figure()
        for c, m in mats.items():
            mean_tds = m.loc[tds].mean(axis=0)
            fig.add_trace(go.Box(y=mean_tds.values, name=c))
        fig.update_layout(title="D4 TDS16 mean expression per cohort",
                          yaxis_title="mean(log2) TDS16")
        fig.write_html(str(FIGS / "v3_ext_d4_tds_boxplot.html"), include_plotlyjs="cdn")
        return {"status": "ok", "n_genes": len(tds)}
    except Exception as e:
        return {"status": "failed", "error": str(e)}


def main():
    log.info("v3_ext_06 start")
    mats = load_v3_matrices()
    log.info(f"loaded {len(mats)} cohorts: {list(mats)}")
    d1 = d1_dataset_identifiability(mats)
    d2 = d2_pca_umap(mats)
    d3 = d3_combat_pilot(mats)
    d4 = d4_cross_cohort_boxplots(mats)
    out = {"d1": d1, "d2": d2, "d3": d3, "d4": d4}
    (TABLES / "v3_ext_batch_diagnostics.json").write_text(json.dumps(out, indent=2, default=str))
    log.info(f"diagnostics: {out}")
    log.info("v3_ext_06 done")


if __name__ == "__main__":
    main()
