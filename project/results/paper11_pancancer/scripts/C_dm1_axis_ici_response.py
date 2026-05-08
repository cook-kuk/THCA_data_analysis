#!/usr/bin/env python3
"""Paper 11 — Phase C: DM1-axis inflammation signature × ICI response.

Re-uses existing Track B-lite per-cohort module scores (HLA-I/II, IFNG,
TLS, checkpoint, myeloid, thyroid_diff) and tests whether a
**lineage-agnostic** DM1 inflammation axis predicts ICI response.

DM1 inflammation composite (no thyroid effectors):
  DM1_inflam = mean(IFNG_T_cell_inflamed, myeloid_suppressive,
                     checkpoint_exhaustion, HLA_class_II)

Outputs:
  - phase_C_ICI/dm1_inflam_per_sample.tsv
  - phase_C_ICI/dm1_inflam_response_per_cohort.tsv (Cohen's d, AUC, p)
  - phase_C_ICI/forest_plot.png
  - phase_C_ICI/summary.json
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd

TRACK_B = Path("/data/thca/repo_results/paper3_ici_track_b_lite")
OUT = Path("/data/thca/repo_results/paper11_pancancer/phase_C_ICI")

INFLAM_MODULES = ["IFNG_T_cell_inflamed", "myeloid_suppressive",
                  "checkpoint_exhaustion", "HLA_class_II"]

# Cohorts where ICI response is meaningful (skip pure thyroid cohorts)
ICI_COHORTS = {
    "Hugo_GSE78220": "Hugo_GSE78220_module_scores.tsv",
    "Riaz_GSE91061_pre": "Riaz_GSE91061_pre_module_scores.tsv",
}


def cohens_d(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2: return np.nan
    s = np.sqrt((np.var(a, ddof=1) + np.var(b, ddof=1)) / 2 + 1e-12)
    return (np.mean(a) - np.mean(b)) / s


def auc_score(y, score):
    try:
        from sklearn.metrics import roc_auc_score
        return float(roc_auc_score(y, score))
    except Exception:
        return float("nan")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows_per_sample = []
    for cohort, fname in ICI_COHORTS.items():
        sp = TRACK_B / "scores_per_cohort" / fname
        if not sp.exists():
            print(f"[C] skip {cohort}: missing {sp}")
            continue
        df = pd.read_csv(sp, sep="\t")
        # composite
        df["DM1_inflam"] = df[INFLAM_MODULES].mean(axis=1)
        rows_per_sample.append(df.assign(cohort=cohort))
    if not rows_per_sample:
        print("[C] no ICI cohorts found — abort")
        return
    all_scores = pd.concat(rows_per_sample, ignore_index=True)
    all_scores.to_csv(OUT / "dm1_inflam_per_sample.tsv", sep="\t", index=False)

    # Response data: pull from existing dial_lite_hugo_riaz.tsv (response class
    # mapping is encoded by aggregate module x response there). We need
    # per-sample response — which is in riaz_pre_to_on_delta_with_response.tsv
    # for Riaz; Hugo lives in track_b_lite_phase2.py logic. Rebuild minimally:
    response_map = {}

    # Riaz pre: response is "PD" / "SD" / "PR" / "CR" — RECIST. Bin to responder.
    riaz_resp = TRACK_B / "riaz_pre_to_on_delta_with_response.tsv"
    if riaz_resp.exists():
        rd = pd.read_csv(riaz_resp, sep="\t")
        # patient col, resp_class col (1.0/0.0)
        for _, r in rd.iterrows():
            v = r.get("resp_class")
            if pd.isna(v): continue
            # Riaz patient column has format Pt101 — but score sample IDs are Pt101_Pre_*
            response_map[("Riaz_GSE91061_pre", str(r["patient"]))] = int(v)
    # Hugo: load from sample metadata if present
    pheno = Path("/home/seungho/personal/THCA_data_analysis/project/results/paper3_ici_data_registry/paper3_ici_public_sample_metadata.tsv")
    if pheno.exists():
        ph = pd.read_csv(pheno, sep="\t")
        for _, r in ph.iterrows():
            ds = str(r.get("dataset_id", ""))
            title = str(r.get("sample_title", ""))   # e.g. "Pt1"
            resp = str(r.get("response", "")).strip().upper()
            if ds == "GSE78220":
                # Sample id in score table is "Pt1", "Pt2" — matches sample_title
                if resp in {"COMPLETE RESPONSE", "PARTIAL RESPONSE"}:
                    response_map[("Hugo_GSE78220", title)] = 1
                elif resp in {"STABLE DISEASE", "PROGRESSIVE DISEASE"}:
                    response_map[("Hugo_GSE78220", title)] = 0

    # Per-cohort effect — for Riaz, the score sample id like "Pt101_Pre_AD..." needs
    # mapping to "Pt101" for the response map.
    def lookup(cohort, sid):
        if (cohort, sid) in response_map:
            return response_map[(cohort, sid)]
        if cohort.startswith("Riaz"):
            # strip _Pre_xxx suffix
            head = str(sid).split("_")[0]
            return response_map.get((cohort, head), np.nan)
        return np.nan

    rows = []
    for cohort, sub in all_scores.groupby("cohort"):
        sub = sub.copy()
        sub["resp"] = sub.apply(lambda r: lookup(cohort, r["sample_id"]), axis=1)
        sub_resp = sub.dropna(subset=["resp"])
        if len(sub_resp) < 5:
            rows.append({"cohort": cohort, "n_total": len(sub),
                         "n_with_response": len(sub_resp),
                         "n_resp": int(sub_resp["resp"].sum()) if len(sub_resp) else 0,
                         "warning": "no/insufficient response labels"})
            continue
        n_resp = int(sub_resp["resp"].sum())
        n_non = int((sub_resp["resp"] == 0).sum())
        d = cohens_d(sub_resp.loc[sub_resp.resp == 1, "DM1_inflam"],
                     sub_resp.loc[sub_resp.resp == 0, "DM1_inflam"])
        auc = auc_score(sub_resp["resp"].astype(int), sub_resp["DM1_inflam"])
        from scipy import stats as st
        try:
            t, p = st.ttest_ind(
                sub_resp.loc[sub_resp.resp == 1, "DM1_inflam"],
                sub_resp.loc[sub_resp.resp == 0, "DM1_inflam"],
                equal_var=False)
        except Exception:
            t, p = np.nan, np.nan
        rows.append({"cohort": cohort, "n_total": len(sub),
                     "n_with_response": len(sub_resp),
                     "n_resp": n_resp, "n_non": n_non,
                     "cohens_d": float(d) if not np.isnan(d) else np.nan,
                     "auc": auc,
                     "t": float(t) if not np.isnan(t) else np.nan,
                     "p": float(p) if not np.isnan(p) else np.nan})
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "dm1_inflam_response_per_cohort.tsv", sep="\t", index=False)
    print("[C] per-cohort effect:")
    print(res.to_string(index=False))

    # forest
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plot_df = res.dropna(subset=["cohens_d"])
        if len(plot_df) > 0:
            fig, ax = plt.subplots(figsize=(7, 3))
            y = np.arange(len(plot_df))
            ax.errorbar(plot_df["cohens_d"], y, fmt="o", capsize=3,
                        color="black")
            ax.axvline(0, color="grey", ls=":")
            ax.set_yticks(y)
            ax.set_yticklabels(
                [f"{r.cohort} (n_R={int(r.n_resp)}, n_NR={int(r.n_non)})"
                 for r in plot_df.itertuples()])
            ax.set_xlabel("Cohen's d (DM1 inflam, R vs NR)")
            ax.set_title("DM1 inflammation axis × ICI response")
            plt.tight_layout()
            plt.savefig(OUT / "forest_plot.png", dpi=140)
            plt.close()
    except Exception as e:
        print(f"[C] plot failed: {e}")

    summary = {
        "ici_cohorts_tested": list(ICI_COHORTS),
        "modules_in_composite": INFLAM_MODULES,
        "per_cohort_results": res.to_dict(orient="records"),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"[C] summary → {OUT / 'summary.json'}")


if __name__ == "__main__":
    main()
