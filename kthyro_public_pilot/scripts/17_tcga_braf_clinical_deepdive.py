#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import kruskal, mannwhitneyu, spearmanr

from pilot_utils import fdr_bh, pilot_root, setup_logging


AXES = [
    "rai_differentiation_score",
    "hla_i_apm_score",
    "immune_visibility_score",
    "cd8_exclusion_proxy",
    "myeloid_caf_barrier_score",
    "drug_delivery_failure_proxy",
    "aggressive_dedifferentiation_score",
    "proliferation_score",
]

AXIS_LABELS = {
    "rai_differentiation_score": "RAI differentiation",
    "hla_i_apm_score": "HLA-I/APM",
    "immune_visibility_score": "Immune visibility",
    "cd8_exclusion_proxy": "CD8 exclusion",
    "myeloid_caf_barrier_score": "Myeloid/CAF barrier",
    "drug_delivery_failure_proxy": "Drug-delivery failure proxy",
    "aggressive_dedifferentiation_score": "Aggressive dedifferentiation",
    "proliferation_score": "Proliferation",
}


def cohens_d(a: pd.Series, b: pd.Series) -> float:
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) < 2 or len(b) < 2:
        return np.nan
    pooled = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return float((a.mean() - b.mean()) / pooled) if np.isfinite(pooled) and pooled > 0 else np.nan


def stage_ord(x: object) -> float:
    s = str(x).upper()
    if not s or s == "NAN":
        return np.nan
    if "STAGE IV" in s:
        return 4
    if "STAGE III" in s:
        return 3
    if "STAGE II" in s:
        return 2
    if "STAGE I" in s:
        return 1
    m = re.search(r"\b([IVX]+)\b", s)
    if m:
        token = m.group(1)
        return {"I": 1, "II": 2, "III": 3, "IV": 4}.get(token, np.nan)
    return np.nan


def binary_event_col(df: pd.DataFrame, event_col: str) -> pd.Series:
    if event_col not in df.columns:
        return pd.Series(np.nan, index=df.index)
    s = df[event_col]
    if s.dtype.kind in "biufc":
        return pd.to_numeric(s, errors="coerce")
    return s.astype(str).str.upper().map({"1": 1, "TRUE": 1, "DECEASED": 1, "DEAD": 1, "0": 0, "FALSE": 0, "ALIVE": 0})


def braf_deepdive(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    braf = df[df["driver_anchor"].fillna("").eq("BRAF")].copy()
    label_counts = pd.crosstab(braf["primary_vulnerability_label"], columns="n").reset_index()
    label_counts["fraction_of_braf"] = label_counts["n"] / label_counts["n"].sum()
    label_counts["claim_boundary"] = "BRAF-only heterogeneity; mutation status is not dismissed, but it is insufficient to determine vulnerability state."

    rows = []
    for label, ldf in braf.groupby("primary_vulnerability_label", sort=False):
        row = {"primary_vulnerability_label": label, "n_braf_patients": len(ldf)}
        for axis in AXES:
            row[f"{axis}_mean"] = pd.to_numeric(ldf[axis], errors="coerce").mean()
            row[f"{axis}_median"] = pd.to_numeric(ldf[axis], errors="coerce").median()
        rows.append(row)
    score_summary = pd.DataFrame(rows)

    test_rows = []
    for axis in AXES:
        groups = [g[axis].dropna().to_numpy() for _, g in braf.groupby("primary_vulnerability_label") if len(g[axis].dropna()) >= 2]
        p = kruskal(*groups).pvalue if len(groups) >= 2 else np.nan
        test_rows.append({"axis": axis, "axis_label": AXIS_LABELS.get(axis, axis), "kruskal_p_across_braf_labels": p})
    tests = pd.DataFrame(test_rows)
    tests["fdr_q"] = fdr_bh(tests["kruskal_p_across_braf_labels"])
    score_summary = score_summary.merge(pd.DataFrame([{"dummy": 1}]), how="cross").drop(columns="dummy")
    return label_counts, tests


def clinical_associations(df: pd.DataFrame) -> pd.DataFrame:
    work = df.copy()
    work["stage_order"] = work["stage"].combine_first(work.get("ajcc_pathologic_tumor_stage", pd.Series(index=work.index))).map(stage_ord)
    work["stage_high_III_IV"] = work["stage_order"].ge(3)
    if "tumor_status" in work.columns:
        work["tumor_present_status"] = work["tumor_status"].astype(str).str.upper().str.contains("WITH|TUMOR", regex=True)
        work.loc[work["tumor_status"].isna(), "tumor_present_status"] = np.nan
    else:
        work["tumor_present_status"] = np.nan
    rows = []
    for axis in AXES:
        vals = pd.to_numeric(work[axis], errors="coerce")
        st = pd.DataFrame({"stage": work["stage_order"], "value": vals}).dropna()
        rho, sp_p = (np.nan, np.nan)
        if st["stage"].nunique() > 1 and st["value"].nunique() > 1:
            rho, sp_p = spearmanr(st["stage"], st["value"])
        high = vals[work["stage_high_III_IV"].eq(True)]
        low = vals[work["stage_order"].isin([1, 2])]
        mw_p = mannwhitneyu(high.dropna(), low.dropna(), alternative="two-sided").pvalue if len(high.dropna()) >= 2 and len(low.dropna()) >= 2 else np.nan
        rows.append(
            {
                "endpoint": "stage_III_IV_vs_I_II",
                "axis": axis,
                "axis_label": AXIS_LABELS.get(axis, axis),
                "n_high_or_event": int(high.dropna().shape[0]),
                "n_low_or_nonevent": int(low.dropna().shape[0]),
                "effect_cohens_d_high_minus_low": cohens_d(high, low),
                "spearman_stage_rho": rho,
                "spearman_p": sp_p,
                "mannwhitney_p": mw_p,
                "claim_boundary": "Exploratory TCGA clinical association; not a validated prognostic or predictive biomarker.",
            }
        )
        present = vals[work["tumor_present_status"].eq(True)]
        absent = vals[work["tumor_present_status"].eq(False)]
        mw_p2 = mannwhitneyu(present.dropna(), absent.dropna(), alternative="two-sided").pvalue if len(present.dropna()) >= 2 and len(absent.dropna()) >= 2 else np.nan
        rows.append(
            {
                "endpoint": "tumor_present_vs_tumor_free_status",
                "axis": axis,
                "axis_label": AXIS_LABELS.get(axis, axis),
                "n_high_or_event": int(present.dropna().shape[0]),
                "n_low_or_nonevent": int(absent.dropna().shape[0]),
                "effect_cohens_d_high_minus_low": cohens_d(present, absent),
                "spearman_stage_rho": np.nan,
                "spearman_p": np.nan,
                "mannwhitney_p": mw_p2,
                "claim_boundary": "Exploratory TCGA clinical association; tumor status fields are heterogeneous.",
            }
        )
    out = pd.DataFrame(rows)
    out["fdr_q_mannwhitney"] = fdr_bh(out["mannwhitney_p"])
    out["fdr_q_spearman"] = fdr_bh(out["spearman_p"])
    return out.sort_values("fdr_q_mannwhitney", na_position="last")


def survival_exploratory(df: pd.DataFrame) -> pd.DataFrame:
    try:
        from lifelines import CoxPHFitter
        from lifelines.statistics import logrank_test
    except Exception:
        return pd.DataFrame([{"status": "not_run", "reason": "lifelines not available"}])
    rows = []
    outcomes = [("PFI", "PFI.time", "PFI"), ("OS", "OS.time", "OS"), ("DSS", "DSS.time", "DSS")]
    for event_col, time_col, endpoint in outcomes:
        if event_col not in df.columns or time_col not in df.columns:
            continue
        event = pd.to_numeric(df[event_col], errors="coerce")
        time = pd.to_numeric(df[time_col], errors="coerce")
        for axis in AXES:
            vals = pd.to_numeric(df[axis], errors="coerce")
            frame = pd.DataFrame({"time": time, "event": event, "score": vals}).replace([np.inf, -np.inf], np.nan).dropna()
            frame = frame[(frame["time"] > 0) & frame["event"].isin([0, 1])]
            if len(frame) < 50 or frame["event"].sum() < 10 or frame["score"].nunique() < 3:
                continue
            frame["score_z"] = (frame["score"] - frame["score"].mean()) / frame["score"].std(ddof=0)
            try:
                cph = CoxPHFitter()
                cph.fit(frame[["time", "event", "score_z"]], duration_col="time", event_col="event")
                hr = float(np.exp(cph.params_["score_z"]))
                p = float(cph.summary.loc["score_z", "p"])
            except Exception:
                hr, p = np.nan, np.nan
            frame["high"] = frame["score"] >= frame["score"].median()
            hi = frame[frame["high"]]
            lo = frame[~frame["high"]]
            try:
                lr_p = logrank_test(hi["time"], lo["time"], hi["event"], lo["event"]).p_value
            except Exception:
                lr_p = np.nan
            rows.append(
                {
                    "endpoint": endpoint,
                    "axis": axis,
                    "axis_label": AXIS_LABELS.get(axis, axis),
                    "n": len(frame),
                    "n_events": int(frame["event"].sum()),
                    "cox_hr_per_sd": hr,
                    "cox_p": p,
                    "median_split_logrank_p": lr_p,
                    "claim_boundary": "Exploratory survival association only; not validated and not treatment-predictive.",
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty and "status" not in out.columns:
        out["cox_fdr_q"] = fdr_bh(out["cox_p"])
        out["logrank_fdr_q"] = fdr_bh(out["median_split_logrank_p"])
    return out


def plot_braf(df: pd.DataFrame, fig_dir: Path) -> None:
    braf = df[df["driver_anchor"].fillna("").eq("BRAF")].copy()
    means = braf.groupby("primary_vulnerability_label")[AXES].mean()
    means = means.rename(columns=AXIS_LABELS)
    fig, ax = plt.subplots(figsize=(12.5, 5.7))
    sns.heatmap(means, cmap="vlag", center=0, linewidths=0.4, linecolor="#222", ax=ax)
    ax.set_title("BRAF-mutant TCGA tumors split into distinct vulnerability states", fontsize=14, weight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(fig_dir / "tcga_braf_only_vulnerability_axis_heatmap.png", dpi=300, bbox_inches="tight")
    fig.savefig(fig_dir / "tcga_braf_only_vulnerability_axis_heatmap.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_clinical(assoc: pd.DataFrame, fig_dir: Path) -> None:
    sub = assoc[assoc["endpoint"].eq("stage_III_IV_vs_I_II")].copy()
    if sub.empty:
        return
    sub = sub.sort_values("effect_cohens_d_high_minus_low")
    fig, ax = plt.subplots(figsize=(9.2, 5.8))
    sns.barplot(data=sub, x="effect_cohens_d_high_minus_low", y="axis_label", color="#fb6a4a", ax=ax)
    ax.axvline(0, color="#333333", lw=1)
    ax.set_xlabel("Cohen's d: stage III/IV minus stage I/II")
    ax.set_ylabel("")
    ax.set_title("TCGA clinical-risk association of K-Thyro axes (exploratory)", fontsize=14, weight="bold")
    fig.tight_layout()
    fig.savefig(fig_dir / "tcga_clinical_stage_axis_associations.png", dpi=300, bbox_inches="tight")
    fig.savefig(fig_dir / "tcga_clinical_stage_axis_associations.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=pilot_root())
    args = parser.parse_args()
    logger = setup_logging("17_tcga_braf_clinical_deepdive")
    out = args.root / "results" / "extra_analyses"
    tables = out / "tables"
    figs = out / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.root / "results" / "tables" / "tcga_thca_patient_vulnerability_scores.tsv", sep="\t")
    label_counts, braf_tests = braf_deepdive(df)
    clinical = clinical_associations(df)
    surv = survival_exploratory(df)
    label_counts.to_csv(tables / "tcga_braf_only_label_distribution.tsv", sep="\t", index=False)
    braf_tests.to_csv(tables / "tcga_braf_only_axis_heterogeneity_tests.tsv", sep="\t", index=False)
    clinical.to_csv(tables / "tcga_clinical_axis_associations_expanded.tsv", sep="\t", index=False)
    surv.to_csv(tables / "tcga_survival_exploratory_axis_associations.tsv", sep="\t", index=False)
    plot_braf(df, figs)
    plot_clinical(clinical, figs)
    logger.info("Wrote BRAF/clinical deep dive outputs.")


if __name__ == "__main__":
    main()
