"""v17 Q18 — combined BRAF x TERT survival meta-analysis.

Cohorts with patient-level OS:
  - TCGA-THCA            (n=502 in v3 cache, 16 events)
  - Landa 2016 / MSK     (n=43,  22 events)
  - Pozdeyev 2018        — no OS in supp / cohort_clinical => prevalence only
  - Lu 2023 (GSE193581)  — scRNA, no survival

Outputs:
  results/v17_survival_meta/q18_combined_survival_table.tsv
  results/v17_survival_meta/q18_summary.json
  results/v17_survival_meta/v17_q18_combined_KM.html
  results/v17_survival_meta/v17_q18_forest_HR.html
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from plotly.subplots import make_subplots
from scipy import stats

RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)

ROOT = Path("/opt/thyroid-dash/project")
SRC = ROOT / "results/v17_tert_recovery/v3/v3_external_combined_BRAF_TERT.tsv"
POZ = ROOT / "results/v17_pozdeyev/Pozdeyev2018_per_patient.tsv"
OUT = ROOT / "results/v17_survival_meta"
OUT.mkdir(parents=True, exist_ok=True)

BG = "#0b0e12"
GROUP_ORDER = ["BRAF-/TERT-", "BRAF+/TERT-", "BRAF-/TERT+", "BRAF+/TERT+"]
GROUP_COLORS = {
    "BRAF-/TERT-": "#7fdbff",  # cyan
    "BRAF+/TERT-": "#ffd166",  # amber
    "BRAF-/TERT+": "#a78bfa",  # violet
    "BRAF+/TERT+": "#ef476f",  # crimson
}


# --------------------------------------------------------------------------- #
# 1. load & prepare cohorts with OS
# --------------------------------------------------------------------------- #
df = pd.read_csv(SRC, sep="\t")
df = df.dropna(subset=["os_event", "os_days"]).copy()
df["os_event"] = df["os_event"].astype(int)
df["os_days"] = df["os_days"].astype(float)
df = df[df["os_days"] >= 0]

# BRAF V600E binary (V600E vs everything else)
df["braf_pos"] = (df["braf_class"] == "V600E").astype(int)
df["tert_pos"] = (df["tert"] == "mutated").astype(int)
df["group"] = np.where(
    df["braf_pos"] & df["tert_pos"], "BRAF+/TERT+",
    np.where(df["braf_pos"] & ~df["tert_pos"].astype(bool), "BRAF+/TERT-",
             np.where(~df["braf_pos"].astype(bool) & df["tert_pos"], "BRAF-/TERT+",
                      "BRAF-/TERT-")))

surv_table_rows = []
for cohort, sub in df.groupby("cohort"):
    for g in GROUP_ORDER:
        s = sub[sub["group"] == g]
        surv_table_rows.append(dict(
            cohort=cohort, group=g, n=len(s),
            events=int(s["os_event"].sum()),
            median_os_days=float(s["os_days"].median()) if len(s) else float("nan"),
        ))
# combined-cohort rows
for g in GROUP_ORDER:
    s = df[df["group"] == g]
    surv_table_rows.append(dict(
        cohort="COMBINED", group=g, n=len(s),
        events=int(s["os_event"].sum()),
        median_os_days=float(s["os_days"].median()) if len(s) else float("nan"),
    ))
surv_table = pd.DataFrame(surv_table_rows)
surv_table.to_csv(OUT / "q18_combined_survival_table.tsv", sep="\t", index=False)

# --------------------------------------------------------------------------- #
# 2. Cohort-stratified Cox (BRAF, TERT, BRAF*TERT) with cohort as strata
# --------------------------------------------------------------------------- #
def cox_fit(data: pd.DataFrame, formula_cols: list[str], strata: list[str] | None = None):
    cph = CoxPHFitter()
    keep = ["os_days", "os_event"] + formula_cols + (strata or [])
    cph.fit(data[keep].dropna(),
            duration_col="os_days", event_col="os_event",
            strata=strata)
    return cph


df["braf_x_tert"] = df["braf_pos"] * df["tert_pos"]

# combined Cox stratified on cohort
cph_combined = cox_fit(df, ["braf_pos", "tert_pos", "braf_x_tert"], strata=["cohort"])
hr_table_combined = cph_combined.summary[["coef", "exp(coef)",
                                          "exp(coef) lower 95%",
                                          "exp(coef) upper 95%", "p"]]
hr_table_combined.columns = ["coef", "HR", "HR_lo", "HR_hi", "p"]

# per-cohort Cox (no strata)
per_cohort_hr: dict[str, dict] = {}
for cohort, sub in df.groupby("cohort"):
    # need at least one event in BRAF+/TERT+ for the interaction; fall back to main effects
    try:
        cph = cox_fit(sub, ["braf_pos", "tert_pos"])
        for cov in ["braf_pos", "tert_pos"]:
            row = cph.summary.loc[cov]
            per_cohort_hr.setdefault(cov, {})[cohort] = dict(
                HR=float(row["exp(coef)"]),
                HR_lo=float(row["exp(coef) lower 95%"]),
                HR_hi=float(row["exp(coef) upper 95%"]),
                p=float(row["p"]),
                n=int(len(sub)), events=int(sub["os_event"].sum()),
            )
    except Exception as e:
        per_cohort_hr.setdefault("error", {})[cohort] = str(e)

# Cox for BRAF+/TERT+ vs the rest (most clinically meaningful flag)
df["braftert_pos"] = ((df["braf_pos"] == 1) & (df["tert_pos"] == 1)).astype(int)
cph_dual = cox_fit(df, ["braftert_pos"], strata=["cohort"])
dual_row = cph_dual.summary.loc["braftert_pos"]
dual_hr = dict(
    HR=float(dual_row["exp(coef)"]),
    HR_lo=float(dual_row["exp(coef) lower 95%"]),
    HR_hi=float(dual_row["exp(coef) upper 95%"]),
    p=float(dual_row["p"]),
    n=int(len(df)), events=int(df["os_event"].sum()),
)

# per-cohort dual flag HR (for forest)
per_cohort_dual: dict[str, dict] = {}
for cohort, sub in df.groupby("cohort"):
    sub2 = sub.copy()
    sub2["braftert_pos"] = ((sub2["braf_pos"] == 1) & (sub2["tert_pos"] == 1)).astype(int)
    if sub2["braftert_pos"].nunique() < 2 or sub2["os_event"].sum() < 2:
        per_cohort_dual[cohort] = dict(HR=float("nan"), HR_lo=float("nan"),
                                       HR_hi=float("nan"), p=float("nan"),
                                       n=int(len(sub2)), events=int(sub2["os_event"].sum()),
                                       note="insufficient events / no contrast")
        continue
    try:
        cph = cox_fit(sub2, ["braftert_pos"])
        row = cph.summary.loc["braftert_pos"]
        per_cohort_dual[cohort] = dict(
            HR=float(row["exp(coef)"]),
            HR_lo=float(row["exp(coef) lower 95%"]),
            HR_hi=float(row["exp(coef) upper 95%"]),
            p=float(row["p"]),
            n=int(len(sub2)), events=int(sub2["os_event"].sum()),
        )
    except Exception as e:
        per_cohort_dual[cohort] = dict(HR=float("nan"), HR_lo=float("nan"),
                                       HR_hi=float("nan"), p=float("nan"),
                                       n=int(len(sub2)), events=int(sub2["os_event"].sum()),
                                       note=str(e))

# 4-group multivariate logrank (combined cohort)
mlr = multivariate_logrank_test(df["os_days"], df["group"], df["os_event"])
logrank_p = float(mlr.p_value)
logrank_chi2 = float(mlr.test_statistic)

# Pozdeyev prevalence context
poz = pd.read_csv(POZ, sep="\t")
poz_overall_braf_tert = int(((poz["BRAF_V600E"] == 1) & (poz["TERT_promoter"] == 1)).sum())
poz_overall_n = int(len(poz))
poz_ptc = poz[poz["Tumor_Type"].str.contains("Papillary", na=False)]
poz_ptc_braf_tert = int(((poz_ptc["BRAF_V600E"] == 1) & (poz_ptc["TERT_promoter"] == 1)).sum())

# --------------------------------------------------------------------------- #
# 3. KM figure (4 groups, combined cohort)
# --------------------------------------------------------------------------- #
def km_curve(times, events):
    kmf = KaplanMeierFitter()
    kmf.fit(times, events)
    return kmf


fig_km = go.Figure()
for g in GROUP_ORDER:
    s = df[df["group"] == g]
    if len(s) == 0:
        continue
    kmf = km_curve(s["os_days"].values, s["os_event"].values)
    sf = kmf.survival_function_
    ci = kmf.confidence_interval_
    x = sf.index.values / 365.25  # years
    y = sf.iloc[:, 0].values
    lo = ci.iloc[:, 0].values
    hi = ci.iloc[:, 1].values
    color = GROUP_COLORS[g]
    label = f"{g}  (n={len(s)}, ev={int(s['os_event'].sum())})"
    fig_km.add_trace(go.Scatter(
        x=np.concatenate([x, x[::-1]]),
        y=np.concatenate([hi, lo[::-1]]),
        fill="toself", fillcolor=f"rgba({int(color[1:3], 16)},{int(color[3:5], 16)},{int(color[5:7], 16)},0.12)",
        line=dict(color="rgba(0,0,0,0)"), hoverinfo="skip",
        showlegend=False, name=g + " CI",
    ))
    fig_km.add_trace(go.Scatter(
        x=x, y=y, mode="lines",
        line=dict(color=color, width=2.4, shape="hv"),
        name=label,
    ))

ann = (f"4-group log-rank: chi2 = {logrank_chi2:.1f}, p = {logrank_p:.2e}<br>"
       f"BRAF+/TERT+ vs rest (cohort-stratified Cox): "
       f"HR = {dual_hr['HR']:.2f} ({dual_hr['HR_lo']:.2f}-{dual_hr['HR_hi']:.2f}), "
       f"p = {dual_hr['p']:.2e}<br>"
       f"Combined n = {len(df)}  events = {int(df['os_event'].sum())} "
       f"(TCGA-THCA + MSK-thyroid-2016)")

fig_km.update_layout(
    template="plotly_dark",
    paper_bgcolor=BG, plot_bgcolor=BG,
    title=dict(text="BRAF V600E x TERT promoter — overall survival (combined external cohorts)",
               x=0.02, xanchor="left", font=dict(size=16)),
    xaxis=dict(title="Years from diagnosis", gridcolor="#1c2230",
               zerolinecolor="#1c2230"),
    yaxis=dict(title="Overall survival probability", range=[0, 1.02],
               gridcolor="#1c2230", zerolinecolor="#1c2230"),
    legend=dict(bgcolor="rgba(0,0,0,0)", x=0.02, y=0.18, font=dict(size=11)),
    margin=dict(l=70, r=30, t=70, b=60),
    annotations=[dict(xref="paper", yref="paper", x=0.98, y=0.98,
                      xanchor="right", yanchor="top",
                      align="right", text=ann, showarrow=False,
                      bgcolor="rgba(11,14,18,0.6)",
                      bordercolor="#2a3142", borderwidth=1,
                      font=dict(size=11, color="#cbd5e1"))],
)
fig_km.write_html(OUT / "v17_q18_combined_KM.html", include_plotlyjs="cdn")

# --------------------------------------------------------------------------- #
# 4. Forest plot — BRAF+/TERT+ HR per cohort + meta
# --------------------------------------------------------------------------- #
forest_rows = []
for cohort in sorted(per_cohort_dual):
    rec = per_cohort_dual[cohort]
    forest_rows.append(dict(label=f"{cohort}  (n={rec['n']}, ev={rec['events']})",
                            HR=rec["HR"], lo=rec["HR_lo"], hi=rec["HR_hi"],
                            p=rec["p"], kind="cohort"))
forest_rows.append(dict(label=f"COMBINED (cohort-stratified)  (n={dual_hr['n']}, ev={dual_hr['events']})",
                        HR=dual_hr["HR"], lo=dual_hr["HR_lo"], hi=dual_hr["HR_hi"],
                        p=dual_hr["p"], kind="meta"))

# add Pozdeyev as prevalence-only context (no HR — annotation panel)
forest_df = pd.DataFrame(forest_rows)
labels = forest_df["label"].tolist()
y_pos = list(range(len(labels), 0, -1))

fig_forest = go.Figure()
for i, row in enumerate(forest_df.itertuples()):
    color = "#ef476f" if row.kind == "meta" else "#7fdbff"
    if not (math.isnan(row.HR) or math.isnan(row.lo) or math.isnan(row.hi)):
        fig_forest.add_trace(go.Scatter(
            x=[row.lo, row.hi], y=[y_pos[i], y_pos[i]],
            mode="lines", line=dict(color=color, width=2),
            showlegend=False, hoverinfo="skip",
        ))
        fig_forest.add_trace(go.Scatter(
            x=[row.HR], y=[y_pos[i]], mode="markers",
            marker=dict(color=color, size=14 if row.kind == "meta" else 11,
                        symbol="diamond" if row.kind == "meta" else "square",
                        line=dict(color="#0b0e12", width=1)),
            text=[f"HR={row.HR:.2f} ({row.lo:.2f}-{row.hi:.2f})  p={row.p:.2g}"],
            hovertemplate="%{text}<extra></extra>",
            showlegend=False,
        ))
    else:
        fig_forest.add_annotation(x=0, y=y_pos[i], xref="x", yref="y",
                                  text="(insufficient events to estimate)",
                                  showarrow=False, font=dict(color="#94a3b8", size=10),
                                  xanchor="left")

fig_forest.add_vline(x=1, line=dict(color="#94a3b8", dash="dash", width=1))

# right-side annotation: Pozdeyev prevalence + Lu 2023 note
context_text = (
    f"Cohorts without OS (context only):<br>"
    f"  Pozdeyev 2018 — BRAF+/TERT+ prevalence "
    f"{poz_overall_braf_tert}/{poz_overall_n} "
    f"({100*poz_overall_braf_tert/poz_overall_n:.1f}%); "
    f"PTC subset {poz_ptc_braf_tert}/{len(poz_ptc)} "
    f"({100*poz_ptc_braf_tert/len(poz_ptc):.1f}%)<br>"
    f"  Lu 2023 (GSE193581) — scRNA only, no survival data"
)

fig_forest.update_layout(
    template="plotly_dark",
    paper_bgcolor=BG, plot_bgcolor=BG,
    title=dict(text="Forest plot — BRAF+/TERT+ vs all others (overall survival HR)",
               x=0.02, xanchor="left", font=dict(size=16)),
    xaxis=dict(title="Hazard ratio (log scale)", type="log",
               gridcolor="#1c2230", zerolinecolor="#1c2230"),
    yaxis=dict(tickmode="array", tickvals=y_pos, ticktext=labels,
               gridcolor="#1c2230", zerolinecolor="#1c2230",
               range=[0.4, len(labels) + 0.6]),
    margin=dict(l=260, r=30, t=70, b=110),
    annotations=[dict(xref="paper", yref="paper", x=0.5, y=-0.22,
                      xanchor="center", yanchor="top",
                      text=context_text, showarrow=False,
                      bgcolor="rgba(11,14,18,0.6)",
                      bordercolor="#2a3142", borderwidth=1,
                      align="left",
                      font=dict(size=11, color="#cbd5e1"))],
    height=420,
)
fig_forest.write_html(OUT / "v17_q18_forest_HR.html", include_plotlyjs="cdn")

# --------------------------------------------------------------------------- #
# 5. JSON summary
# --------------------------------------------------------------------------- #
summary = dict(
    random_state=RANDOM_STATE,
    inputs=dict(
        v3_combined_tsv=str(SRC),
        pozdeyev_tsv=str(POZ),
    ),
    cohorts_with_OS=dict(
        TCGA_THCA=int((df["cohort"] == "TCGA-THCA").sum()),
        MSK_thyroid_2016=int((df["cohort"] == "MSK-thyroid-2016").sum()),
        combined_n=int(len(df)),
        combined_events=int(df["os_event"].sum()),
    ),
    cohorts_without_OS=dict(
        Pozdeyev_2018=dict(
            n=poz_overall_n,
            braf_v600e_n=int((poz["BRAF_V600E"] == 1).sum()),
            tert_promoter_n=int((poz["TERT_promoter"] == 1).sum()),
            braf_and_tert_n=poz_overall_braf_tert,
            ptc_n=int(len(poz_ptc)),
            ptc_braf_and_tert_n=poz_ptc_braf_tert,
            note="No patient-level OS in supplement-8 nor cohort_clinical.tsv; used as prevalence context.",
        ),
        Lu_2023_GSE193581=dict(note="scRNA-seq only; no clinical survival data."),
    ),
    four_group_logrank=dict(chi2=logrank_chi2, p=logrank_p),
    cox_combined_stratified_by_cohort=dict(
        covariates=hr_table_combined.to_dict(orient="index"),
    ),
    cox_per_cohort_main_effects=per_cohort_hr,
    cox_braf_tert_dual_vs_rest=dict(
        per_cohort=per_cohort_dual,
        meta_stratified=dual_hr,
    ),
    figures=dict(
        km=str(OUT / "v17_q18_combined_KM.html"),
        forest=str(OUT / "v17_q18_forest_HR.html"),
    ),
)
with open(OUT / "q18_summary.json", "w") as f:
    json.dump(summary, f, indent=2, default=lambda o: None if (isinstance(o, float) and math.isnan(o)) else o)

print(json.dumps(summary, indent=2, default=lambda o: None if (isinstance(o, float) and math.isnan(o)) else o))
