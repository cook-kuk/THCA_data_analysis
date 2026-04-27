#!/usr/bin/env python3
"""
v17 REAL FIX · R6 — DM1/DM2 + leak-free cluster survival analysis (OS).

Endpoints:
  - OS (overall survival) by DM1/DM2 (R1-A leak-free)
  - OS by combined DM × TERT axis
  - OS by stage-stratified DM
  - Univariate + multivariate Cox

Outputs:
  - results/v17_realfix/R6_survival_table.tsv
  - results/v17_realfix/R6_summary.json
  - results/v17_realfix/figures/FigR6_*.png
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test

OUT = Path("/opt/thyroid-dash/project/results/v17_realfix")
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

SAMPLE_MASTER = Path("/opt/thyroid-dash/project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv")
R1A_LABELS = OUT / "R1A_cluster_labels.tsv"

plt.rcParams.update({"font.family": "sans-serif", "font.size": 11, "figure.dpi": 130})


def load_data():
    sm = pd.read_csv(SAMPLE_MASTER, sep="\t", low_memory=False)
    print(f"[load] sample_master: {sm.shape}", flush=True)
    # filter: TCGA-THCA primary tumours with OS info
    sm = sm[sm["dataset"] == "TCGA-THCA"]
    sm = sm[sm["normal_vs_tumor"] == "tumor"]
    sm = sm.dropna(subset=["os_days", "os_event"])
    sm = sm[sm["os_days"] > 0]
    sm["os_days"] = pd.to_numeric(sm["os_days"], errors="coerce")
    sm["os_event"] = pd.to_numeric(sm["os_event"], errors="coerce")
    sm = sm.dropna(subset=["os_days", "os_event"])
    sm["os_event"] = sm["os_event"].astype(int)
    print(f"[load] after filter: n={len(sm)}, events={int(sm['os_event'].sum())}", flush=True)

    lbl = pd.read_csv(R1A_LABELS, sep="\t")
    lbl["dm_leak_free"] = lbl["cluster"].map(lambda x: "DM1" if x.startswith("DM1") else "DM2")
    print(f"[load] R1-A labels: n={len(lbl)} (DM1={(lbl['dm_leak_free']=='DM1').sum()}, DM2={(lbl['dm_leak_free']=='DM2').sum()})", flush=True)

    df = sm.merge(lbl[["sample_id", "dm_leak_free"]], on="sample_id", how="inner")
    print(f"[load] joined cohort: n={len(df)}, OS events={int(df['os_event'].sum())}", flush=True)
    return df


def km_2group(df, group_col, fig_path, title, group_labels):
    fig, ax = plt.subplots(figsize=(7, 5))
    kmfs = {}
    palette = {"DM1": "#dc2626", "DM2": "#1d4ed8"}
    for g in group_labels:
        mask = df[group_col] == g
        if mask.sum() < 2:
            continue
        kmf = KaplanMeierFitter()
        kmf.fit(df.loc[mask, "os_days"] / 365.25, df.loc[mask, "os_event"], label=f"{g} (n={int(mask.sum())}, ev={int(df.loc[mask,'os_event'].sum())})")
        kmf.plot_survival_function(ax=ax, ci_show=True, color=palette.get(g, "black"), lw=2.5)
        kmfs[g] = kmf
    if len(kmfs) >= 2:
        gA, gB = list(kmfs.keys())[:2]
        mA = df[group_col] == gA
        mB = df[group_col] == gB
        lr = logrank_test(df.loc[mA, "os_days"], df.loc[mB, "os_days"],
                          df.loc[mA, "os_event"], df.loc[mB, "os_event"])
        p_val = lr.p_value
    else:
        p_val = float("nan")
    ax.set_title(f"{title}\nlogrank p = {p_val:.3g}" if not np.isnan(p_val) else title)
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Overall survival probability")
    ax.set_ylim(0.5, 1.02)
    ax.legend(loc="lower left", fontsize=10)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight", dpi=180)
    fig.savefig(str(fig_path).replace(".png", ".pdf"), bbox_inches="tight")
    plt.close(fig)
    return p_val


def cox_uni_multi(df, group_col, group_labels):
    """Cox with DM1 = ref (0), DM2 = 1; univariate + multivariate (stage + age + sex)."""
    sub = df.copy()
    sub["dm_binary"] = (sub[group_col] == group_labels[1]).astype(int)
    sub["age"] = pd.to_numeric(sub["age_at_diagnosis"], errors="coerce")
    # stage encoding
    stage_map = {"Stage I": 1, "Stage II": 2, "Stage III": 3, "Stage IV": 4,
                 "Stage IVA": 4, "Stage IVB": 4, "Stage IVC": 4}
    sub["stage_int"] = sub["clinical_stage"].astype(str).map(lambda x: stage_map.get(x, np.nan))
    # sex encoding
    sub["sex_male"] = (sub.get("sex_clinical", "").astype(str).str.lower() == "male").astype(int)
    sub_uni = sub.dropna(subset=["dm_binary", "os_days", "os_event"]).copy()
    sub_mul = sub.dropna(subset=["dm_binary", "os_days", "os_event", "stage_int", "age", "sex_male"]).copy()

    out = {}
    if len(sub_uni) >= 30 and sub_uni["os_event"].sum() >= 3:
        try:
            cph_u = CoxPHFitter(penalizer=0.01)
            cph_u.fit(sub_uni[["dm_binary", "os_days", "os_event"]], duration_col="os_days", event_col="os_event")
            row = cph_u.summary.loc["dm_binary"]
            out["univariate"] = {
                "n": int(len(sub_uni)),
                "events": int(sub_uni["os_event"].sum()),
                "HR": float(row["exp(coef)"]),
                "HR_lo": float(row["exp(coef) lower 95%"]),
                "HR_hi": float(row["exp(coef) upper 95%"]),
                "p": float(row["p"]),
            }
        except Exception as e:
            out["univariate"] = {"error": str(e)}

    if len(sub_mul) >= 30 and sub_mul["os_event"].sum() >= 5:
        try:
            cph_m = CoxPHFitter(penalizer=0.01)
            cph_m.fit(sub_mul[["dm_binary", "stage_int", "age", "sex_male", "os_days", "os_event"]],
                      duration_col="os_days", event_col="os_event")
            row = cph_m.summary.loc["dm_binary"]
            out["multivariate"] = {
                "n": int(len(sub_mul)),
                "events": int(sub_mul["os_event"].sum()),
                "HR": float(row["exp(coef)"]),
                "HR_lo": float(row["exp(coef) lower 95%"]),
                "HR_hi": float(row["exp(coef) upper 95%"]),
                "p": float(row["p"]),
                "covariates": "stage_int + age + sex_male",
            }
        except Exception as e:
            out["multivariate"] = {"error": str(e)}

    return out


def km_4group_dm_tert(df, fig_path):
    """4-group: DM2/TERT-, DM2/TERT+, DM1/TERT-, DM1/TERT+."""
    sub = df.copy()
    # tert_status column = "wildtype"/"mutant" string, or numeric 0/1 in tert_promoter_integrated
    if "tert_status" in sub.columns:
        sub["tert_pos"] = (sub["tert_status"].astype(str).str.lower() == "mutant").astype(int)
    elif "tert_promoter_integrated" in sub.columns:
        sub["tert_pos"] = pd.to_numeric(sub["tert_promoter_integrated"], errors="coerce").fillna(0).astype(int)
    else:
        sub["tert_pos"] = 0
    sub["dm_tert"] = sub["dm_leak_free"].astype(str) + "/TERT" + sub["tert_pos"].map({0: "−", 1: "+"})
    fig, ax = plt.subplots(figsize=(8, 5))
    palette = {"DM1/TERT−": "#fbbf24", "DM1/TERT+": "#dc2626",
               "DM2/TERT−": "#1d4ed8", "DM2/TERT+": "#7c3aed"}
    groups = ["DM2/TERT−", "DM2/TERT+", "DM1/TERT−", "DM1/TERT+"]
    for g in groups:
        mask = sub["dm_tert"] == g
        if mask.sum() < 2:
            continue
        kmf = KaplanMeierFitter()
        kmf.fit(sub.loc[mask, "os_days"] / 365.25, sub.loc[mask, "os_event"],
                label=f"{g} (n={int(mask.sum())}, ev={int(sub.loc[mask,'os_event'].sum())})")
        kmf.plot_survival_function(ax=ax, ci_show=False, color=palette.get(g, "gray"), lw=2.5)

    # multivariate logrank
    sub_valid = sub.dropna(subset=["dm_tert", "os_days", "os_event"])
    sub_valid = sub_valid[sub_valid["dm_tert"].isin(groups)]
    if len(sub_valid["dm_tert"].unique()) >= 2:
        mlr = multivariate_logrank_test(sub_valid["os_days"], sub_valid["dm_tert"], sub_valid["os_event"])
        p_val = mlr.p_value
    else:
        p_val = float("nan")
    ax.set_title(f"Combined DM × TERT axis OS\nmultivariate logrank p = {p_val:.3g}")
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Overall survival probability")
    ax.set_ylim(0.5, 1.02)
    ax.legend(loc="lower left", fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight", dpi=180)
    fig.savefig(str(fig_path).replace(".png", ".pdf"), bbox_inches="tight")
    plt.close(fig)
    return p_val


def main():
    df = load_data()
    summary = {"n_total": len(df), "n_events": int(df["os_event"].sum())}

    # 1. DM1 vs DM2 OS (full cohort)
    print("\n=== R6.1 · DM1 vs DM2 OS (full cohort) ===", flush=True)
    p1 = km_2group(df, "dm_leak_free",
                    FIG / "FigR6_1_DM1_vs_DM2_OS.png",
                    "DM1 vs DM2 — OS (R1-A leak-free, n={})".format(len(df)),
                    ["DM1", "DM2"])
    cox1 = cox_uni_multi(df, "dm_leak_free", ["DM1", "DM2"])
    summary["R6_1_DM1_vs_DM2_full"] = {"logrank_p": float(p1), **cox1,
                                        "n_DM1": int((df["dm_leak_free"] == "DM1").sum()),
                                        "n_DM2": int((df["dm_leak_free"] == "DM2").sum())}
    print(f"  logrank p = {p1:.4g}", flush=True)
    print(f"  univariate: {cox1.get('univariate', {})}", flush=True)
    print(f"  multivariate: {cox1.get('multivariate', {})}", flush=True)

    # 2. DM1 vs DM2 OS — Stage III/IV only
    print("\n=== R6.2 · DM1 vs DM2 OS (Stage III/IV) ===", flush=True)
    df_advanced = df[df["clinical_stage"].astype(str).isin(["Stage III", "Stage IV", "Stage IVA", "Stage IVB", "Stage IVC"])]
    if len(df_advanced) >= 20:
        p2 = km_2group(df_advanced, "dm_leak_free",
                        FIG / "FigR6_2_DM1_vs_DM2_OS_StageIII_IV.png",
                        "DM1 vs DM2 — OS (Stage III/IV, n={})".format(len(df_advanced)),
                        ["DM1", "DM2"])
        cox2 = cox_uni_multi(df_advanced, "dm_leak_free", ["DM1", "DM2"])
        summary["R6_2_DM1_vs_DM2_advanced"] = {"logrank_p": float(p2), **cox2,
                                                "n_advanced": int(len(df_advanced))}
        print(f"  logrank p = {p2:.4g}", flush=True)
    else:
        summary["R6_2_DM1_vs_DM2_advanced"] = {"error": f"insufficient n={len(df_advanced)}"}

    # 3. DM1 vs DM2 OS — age >= 45
    print("\n=== R6.3 · DM1 vs DM2 OS (age ≥ 45) ===", flush=True)
    df["age_num"] = pd.to_numeric(df["age_at_diagnosis"], errors="coerce")
    df_old = df[df["age_num"] >= 45]
    if len(df_old) >= 30:
        p3 = km_2group(df_old, "dm_leak_free",
                        FIG / "FigR6_3_DM1_vs_DM2_OS_age45plus.png",
                        "DM1 vs DM2 — OS (age ≥ 45, n={})".format(len(df_old)),
                        ["DM1", "DM2"])
        cox3 = cox_uni_multi(df_old, "dm_leak_free", ["DM1", "DM2"])
        summary["R6_3_DM1_vs_DM2_age45plus"] = {"logrank_p": float(p3), **cox3,
                                                  "n_age45plus": int(len(df_old))}
        print(f"  logrank p = {p3:.4g}", flush=True)

    # 4. Combined DM × TERT 4-group
    print("\n=== R6.4 · DM × TERT 4-group OS ===", flush=True)
    p4 = km_4group_dm_tert(df, FIG / "FigR6_4_DM_TERT_4group_OS.png")
    summary["R6_4_DM_TERT_4group"] = {"multivariate_logrank_p": float(p4)}
    print(f"  multivariate logrank p = {p4:.4g}", flush=True)

    # write
    with open(OUT / "R6_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    # write table
    rows = []
    for key in ["R6_1_DM1_vs_DM2_full", "R6_2_DM1_vs_DM2_advanced", "R6_3_DM1_vs_DM2_age45plus"]:
        r = summary.get(key, {})
        rows.append({
            "analysis": key,
            "n": r.get("univariate", {}).get("n", "—"),
            "events": r.get("univariate", {}).get("events", "—"),
            "logrank_p": r.get("logrank_p", "—"),
            "univariate_HR": r.get("univariate", {}).get("HR", "—"),
            "univariate_p": r.get("univariate", {}).get("p", "—"),
            "multivariate_HR": r.get("multivariate", {}).get("HR", "—"),
            "multivariate_p": r.get("multivariate", {}).get("p", "—"),
        })
    pd.DataFrame(rows).to_csv(OUT / "R6_survival_table.tsv", sep="\t", index=False)

    print("\n=== R6 SUMMARY ===")
    print(json.dumps(summary, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
