#!/usr/bin/env python3
"""
v4 Track B: strengthen the Bethesda-III/IV clinical decision simulation.

Extends the v3 Bethesda simulation to:
- 6 prevalence levels: 5%, 10%, 15%, 20%, 30%, 40%
- 20,000 simulated patients per prevalence
- Cost-utility in KRW (lobectomy 3.5M, FU 0.25M, panel 0.3M, missed-cancer penalty 30M)
- Decision-curve (net-benefit) analysis
- 3 clinical operating points: Se>=0.95, balanced, Sp>=0.95
- patient-impact JSON for Gemma hackathon

All cost assumptions explicitly labeled "Korean single-payer estimates, not generalizable".
Language kept as "decision-support prototype / retrospective computational triage".

Outputs:
  results/tables/v4_bethesda_sim.tsv
  results/tables/v4_bethesda_operating_points.tsv
  results/tables/v4_bethesda_cost_utility.tsv
  results/ml/v4_trackB_patient_impact.json
  reports/html/figs_interactive/v4_bethesda_fig{1..5}.{png,html,tsv}
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path("/opt/thyroid-dash/project")
OUT_ML = ROOT / "results" / "ml"
OUT_FIG = ROOT / "reports" / "html" / "figs_interactive"
OUT_TABLES = ROOT / "results" / "tables"
for d in (OUT_ML, OUT_FIG, OUT_TABLES):
    d.mkdir(parents=True, exist_ok=True)

PREVALENCES = [0.05, 0.10, 0.15, 0.20, 0.30, 0.40]
N_PER = 20_000
THRESHOLDS = np.arange(0.05, 0.91, 0.01)

# KRW cost model (explicit, labeled "Korean single-payer estimates, not generalizable")
COST = {
    "lobectomy_krw": 3_500_000,
    "followup_per_visit_krw": 250_000,
    "genomic_panel_krw": 300_000,
    "missed_cancer_penalty_krw": 30_000_000,
}


def log(msg: str) -> None:
    print(f"[trackB] {msg}", flush=True)


def simulate_one(prev: float, n: int, seed: int = 42) -> tuple[np.ndarray, np.ndarray]:
    """Simulate probability output of a moderately calibrated classifier.
    Malignant distribution Beta(5, 2); Benign distribution Beta(2, 5).
    """
    rng = np.random.default_rng(seed)
    n_pos = int(round(n * prev))
    n_neg = n - n_pos
    p_pos = rng.beta(5.0, 2.0, size=n_pos)
    p_neg = rng.beta(2.0, 5.0, size=n_neg)
    y = np.concatenate([np.ones(n_pos), np.zeros(n_neg)])
    p = np.concatenate([p_pos, p_neg])
    idx = rng.permutation(n)
    return y[idx], p[idx]


def metrics_at_threshold(y: np.ndarray, p: np.ndarray, thr: float, prev: float) -> dict:
    pred = (p >= thr).astype(int)
    tp = int(((pred == 1) & (y == 1)).sum())
    fp = int(((pred == 1) & (y == 0)).sum())
    tn = int(((pred == 0) & (y == 0)).sum())
    fn = int(((pred == 0) & (y == 1)).sum())
    sens = tp / max(tp + fn, 1)
    spec = tn / max(tn + fp, 1)
    ppv = tp / max(tp + fp, 1) if (tp + fp) else 0.0
    npv = tn / max(tn + fn, 1) if (tn + fn) else 0.0
    n = tp + fp + tn + fn
    # Net benefit formula (Vickers)
    pt = thr
    nb_model = (tp / n) - (fp / n) * (pt / (1 - pt)) if (1 - pt) > 0 else 0.0
    nb_all = prev - (1 - prev) * (pt / (1 - pt)) if (1 - pt) > 0 else 0.0
    nb_none = 0.0
    return {
        "threshold": float(thr), "prev": float(prev),
        "tp": tp, "fp": fp, "tn": tn, "fn": fn,
        "sens": float(sens), "spec": float(spec),
        "ppv": float(ppv), "npv": float(npv),
        "unnecessary_surgery_rate": float(fp / n),
        "missed_cancer_rate": float(fn / n),
        "nb_model": float(nb_model), "nb_all": float(nb_all), "nb_none": float(nb_none),
    }


def cost_per_patient(m: dict) -> float:
    n = m["tp"] + m["fp"] + m["tn"] + m["fn"]
    # Treat "positive" = lobectomy; "negative" = surveillance (6 follow-ups over 2y).
    # TP: lobectomy (correct)     = lobectomy + panel
    # FP: lobectomy (unnecessary) = lobectomy + panel
    # TN: follow-up (correct)     = 6 * FU + panel
    # FN: miss cancer             = 6*FU + panel + missed_cancer_penalty
    total = (
        (m["tp"] + m["fp"]) * (COST["lobectomy_krw"] + COST["genomic_panel_krw"]) +
        m["tn"] * (6 * COST["followup_per_visit_krw"] + COST["genomic_panel_krw"]) +
        m["fn"] * (6 * COST["followup_per_visit_krw"] + COST["genomic_panel_krw"] + COST["missed_cancer_penalty_krw"])
    )
    return float(total / max(n, 1))


def treat_all_cost(prev: float) -> float:
    # everyone gets lobectomy + panel
    return COST["lobectomy_krw"] + COST["genomic_panel_krw"]


def treat_none_cost(prev: float) -> float:
    # nobody operated: all get followups; malignant penalty for prev fraction
    return 6 * COST["followup_per_visit_krw"] + COST["genomic_panel_krw"] + prev * COST["missed_cancer_penalty_krw"]


def pick_operating_points(sim_df: pd.DataFrame) -> pd.DataFrame:
    """For each prevalence, pick (a) Se>=0.95, (b) balanced (Se+Sp max), (c) Sp>=0.95."""
    out = []
    for prev, grp in sim_df.groupby("prev"):
        grp = grp.sort_values("threshold")
        # Se>=0.95: largest threshold still with sens>=0.95
        se_high = grp[grp["sens"] >= 0.95]
        if len(se_high):
            r = se_high.iloc[-1].copy()
            r["operating_point"] = "Se>=0.95 (rule-out)"
            out.append(r)
        # Balanced: max sens+spec
        grp["sum_ss"] = grp["sens"] + grp["spec"]
        r = grp.sort_values("sum_ss", ascending=False).iloc[0].copy()
        r["operating_point"] = "balanced (max Se+Sp)"
        out.append(r)
        # Sp>=0.95: smallest threshold with spec>=0.95
        sp_high = grp[grp["spec"] >= 0.95]
        if len(sp_high):
            r = sp_high.iloc[0].copy()
            r["operating_point"] = "Sp>=0.95 (rule-in)"
            out.append(r)
    return pd.DataFrame(out).reset_index(drop=True)


def save_fig(path_base: Path, fig, tsv_df: pd.DataFrame | None = None) -> None:
    fig.savefig(path_base.with_suffix(".png"), dpi=110, bbox_inches="tight")
    plt.close(fig)
    if tsv_df is not None:
        tsv_df.to_csv(path_base.with_suffix(".tsv"), sep="\t", index=False)
    html = f"<html><body><img src='{path_base.name}.png' style='max-width:100%'></body></html>"
    path_base.with_suffix(".html").write_text(html)


def main() -> int:
    log(f"simulating {len(PREVALENCES)} prevalences x {N_PER} patients each...")
    sim_rows = []
    for prev in PREVALENCES:
        y, p = simulate_one(prev, N_PER, seed=42 + int(prev * 100))
        for thr in THRESHOLDS:
            row = metrics_at_threshold(y, p, float(thr), float(prev))
            row["cost_per_patient_krw"] = cost_per_patient(row)
            row["treat_all_cost_krw"] = treat_all_cost(float(prev))
            row["treat_none_cost_krw"] = treat_none_cost(float(prev))
            sim_rows.append(row)
    sim_df = pd.DataFrame(sim_rows)
    sim_df.to_csv(OUT_TABLES / "v4_bethesda_sim.tsv", sep="\t", index=False)
    log(f"sim rows: {len(sim_df)}")

    # Operating points
    op_df = pick_operating_points(sim_df)
    op_df.to_csv(OUT_TABLES / "v4_bethesda_operating_points.tsv", sep="\t", index=False)
    log(f"operating points: {len(op_df)}")

    # Cost-utility: at balanced point for each prevalence
    bal = op_df[op_df["operating_point"] == "balanced (max Se+Sp)"].copy()
    bal["savings_vs_treat_all_krw"] = bal["treat_all_cost_krw"] - bal["cost_per_patient_krw"]
    bal["savings_vs_treat_none_krw"] = bal["treat_none_cost_krw"] - bal["cost_per_patient_krw"]
    bal.to_csv(OUT_TABLES / "v4_bethesda_cost_utility.tsv", sep="\t", index=False)

    # Patient-impact JSON for Gemma
    impact = {
        "_note": (
            "Decision-support prototype / retrospective computational triage. "
            "NOT a diagnostic device. Costs are Korean single-payer estimates, "
            "NOT generalizable to other healthcare systems."
        ),
        "cost_assumptions_krw": COST,
        "prevalences_simulated": PREVALENCES,
        "n_patients_per_prev": N_PER,
        "operating_points": op_df.to_dict(orient="records"),
        "avoided_surgeries_per_1000_at_prev_0.20_balanced": None,
    }
    at20 = op_df[(op_df["prev"] == 0.20) & (op_df["operating_point"] == "balanced (max Se+Sp)")]
    if len(at20):
        r = at20.iloc[0]
        # at 20% prev and sens/spec; treat-all surgery rate = 1.0
        model_surgery_rate = (r["tp"] + r["fp"]) / (r["tp"] + r["fp"] + r["tn"] + r["fn"])
        avoided = (1.0 - model_surgery_rate) * 1000
        impact["avoided_surgeries_per_1000_at_prev_0.20_balanced"] = float(avoided)
    (OUT_ML / "v4_trackB_patient_impact.json").write_text(json.dumps(impact, indent=2))

    # --- Figures ---
    # Fig 1: ROC curves per prevalence (synthetic model → same ROC, so illustrate 1 curve +
    # show net benefit per prevalence instead).
    fig, ax = plt.subplots(figsize=(7, 5))
    for prev in PREVALENCES:
        g = sim_df[sim_df["prev"] == prev].sort_values("sens")
        ax.plot(1 - g["spec"], g["sens"], label=f"prev={prev:.0%}", lw=1.4)
    ax.plot([0, 1], [0, 1], "--", color="gray", lw=1)
    ax.set_xlabel("1 - specificity (FPR)")
    ax.set_ylabel("sensitivity (TPR)")
    ax.set_title("Track B: ROC curves across prevalences (synthetic)")
    ax.legend()
    save_fig(OUT_FIG / "v4_bethesda_fig1", fig, sim_df)

    # Fig 2: net-benefit vs threshold
    fig, ax = plt.subplots(figsize=(8, 5))
    for prev in [0.10, 0.20, 0.30]:
        g = sim_df[sim_df["prev"] == prev].sort_values("threshold")
        ax.plot(g["threshold"], g["nb_model"], label=f"model prev={prev:.0%}", lw=1.6)
        ax.plot(g["threshold"], g["nb_all"], ls="--", alpha=0.5, label=f"treat-all prev={prev:.0%}")
    ax.axhline(0, color="gray", lw=0.8)
    ax.set_xlabel("probability threshold")
    ax.set_ylabel("net benefit (Vickers)")
    ax.set_title("Track B: decision-curve analysis")
    ax.legend(fontsize=8)
    save_fig(OUT_FIG / "v4_bethesda_fig2", fig)

    # Fig 3: cost-per-patient vs prevalence
    fig, ax = plt.subplots(figsize=(7, 5))
    xs = bal["prev"].values
    ax.plot(xs, bal["cost_per_patient_krw"].values / 1e6, "-o", label="model (balanced)", lw=2)
    ax.plot(xs, bal["treat_all_cost_krw"].values / 1e6, "--s", label="treat-all (lobectomy for all)", lw=2)
    ax.plot(xs, bal["treat_none_cost_krw"].values / 1e6, "--^", label="treat-none (surveillance only)", lw=2)
    ax.set_xlabel("malignancy prevalence")
    ax.set_ylabel("expected cost per patient (M KRW)")
    ax.set_title("Track B: cost-utility (Korean single-payer est., not generalizable)")
    ax.legend()
    save_fig(OUT_FIG / "v4_bethesda_fig3", fig, bal)

    # Fig 4: operating points table as bar chart of sens/spec
    fig, ax = plt.subplots(figsize=(9, 5))
    op_df_local = op_df.copy()
    op_df_local["label"] = op_df_local["prev"].map(lambda x: f"{x:.0%}") + " | " + op_df_local["operating_point"]
    x = np.arange(len(op_df_local))
    w = 0.4
    ax.bar(x - w/2, op_df_local["sens"], w, label="sensitivity", color="#3498db")
    ax.bar(x + w/2, op_df_local["spec"], w, label="specificity", color="#e67e22")
    ax.set_xticks(x)
    ax.set_xticklabels(op_df_local["label"], rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("metric")
    ax.set_ylim(0, 1.05)
    ax.set_title("Track B: three operating points across prevalences")
    ax.legend()
    save_fig(OUT_FIG / "v4_bethesda_fig4", fig, op_df_local)

    # Fig 5: unnecessary surgery reduction vs missed cancer at balanced op
    fig, ax = plt.subplots(figsize=(7, 5))
    # at balanced op, x = prevalence, y1 = unnecessary surgery rate vs y2 = missed cancer rate
    ax.plot(bal["prev"], bal["unnecessary_surgery_rate"], "-o", color="#e74c3c", label="unnecessary surgery rate (FP/N)")
    ax.plot(bal["prev"], bal["missed_cancer_rate"], "-s", color="#8e44ad", label="missed cancer rate (FN/N)")
    ax.set_xlabel("prevalence")
    ax.set_ylabel("rate")
    ax.set_title("Track B: error trade-off at balanced operating point")
    ax.legend()
    save_fig(OUT_FIG / "v4_bethesda_fig5", fig, bal)

    log("done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
