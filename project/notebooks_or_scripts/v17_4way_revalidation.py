"""B-prompt: TERT × BRAF × RAS 4-way matrix re-validation with bootstrap CI.

Resolves the TERT-only paradox by exposing the (driver_anchor × TERT) 8-cell
breakdown that the prior 4-group flattening hid.

Outputs to project/results/audit_2026_04_29/4way/.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
from scipy import stats

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
SM_PATH = ROOT / "project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv"
OUT = ROOT / "project/results/audit_2026_04_29/4way"
OUT.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(42)


def load_master() -> pd.DataFrame:
    sm = pd.read_csv(SM_PATH, sep="\t", low_memory=False)
    sm = sm[sm["dataset"] == "TCGA-THCA"].copy()
    sm = sm[sm["normal_vs_tumor"] == "tumor"].copy()
    raw = sm["tert_promoter_integrated"].fillna("").astype(str).str.lower().str.strip()
    sm["tert_int"] = (raw == "mutated").astype(int)
    sm["os_event"] = pd.to_numeric(sm["os_event"], errors="coerce")
    sm["os_days"] = pd.to_numeric(sm["os_days"], errors="coerce")
    return sm


def make_cells(sm: pd.DataFrame) -> pd.DataFrame:
    """driver_anchor × tert_int 8-cell breakdown.

    driver_anchor in {BRAF, RAS, NTRK, fusion, unknown, ...} — mutually exclusive
    in the curated TCGA classification.
    """
    drv = sm["driver_anchor"].fillna("unknown").astype(str).str.upper()
    drv = drv.where(drv.isin(["BRAF", "RAS", "NTRK"]), other="OTHER")
    sm = sm.assign(driver_simple=drv)
    sm["cell"] = sm["driver_simple"] + ("_TERT+" if sm["tert_int"].iloc[0] == 1 else "")
    sm["cell"] = sm.apply(
        lambda r: f"{r['driver_simple']}_TERT+" if r["tert_int"] == 1 else f"{r['driver_simple']}_TERT-",
        axis=1,
    )
    return sm


def crosstab_summary(sm: pd.DataFrame) -> pd.DataFrame:
    rows = []
    surv = sm.dropna(subset=["os_event", "os_days"])
    for cell, grp in surv.groupby("cell"):
        rows.append(
            {
                "cell": cell,
                "n": len(grp),
                "events": int(grp["os_event"].sum()),
                "event_rate_pct": round(100 * grp["os_event"].mean(), 2),
                "median_followup_days": float(grp["os_days"].median()),
            }
        )
    df = pd.DataFrame(rows).sort_values("n", ascending=False)
    return df


def bootstrap_cox_hr(
    sm: pd.DataFrame,
    cell: str,
    reference_cell: str = "OTHER_TERT-",
    n_boot: int = 1000,
) -> dict:
    """Bootstrap Cox HR for cell vs reference. Handles 0-event cells."""
    df = sm.dropna(subset=["os_event", "os_days"]).copy()
    df = df[df["cell"].isin([cell, reference_cell])].copy()
    df["is_target"] = (df["cell"] == cell).astype(int)
    df = df[["os_days", "os_event", "is_target"]].rename(
        columns={"os_days": "T", "os_event": "E", "is_target": "G"}
    )

    n_target = int(df["G"].sum())
    n_ref = int(len(df) - n_target)
    e_target = int(df.loc[df["G"] == 1, "E"].sum())
    e_ref = int(df.loc[df["G"] == 0, "E"].sum())

    cph = CoxPHFitter(penalizer=0.01)  # ridge ~ Firth-like for small N
    try:
        cph.fit(df, duration_col="T", event_col="E", show_progress=False)
        hr = float(np.exp(cph.params_.loc["G"]))
    except Exception:
        hr = np.nan

    boot_hrs = []
    for _ in range(n_boot):
        idx = RNG.integers(0, len(df), size=len(df))
        bs = df.iloc[idx].copy()
        if bs["G"].nunique() < 2 or bs["E"].sum() == 0:
            continue
        try:
            c = CoxPHFitter(penalizer=0.01)
            c.fit(bs, duration_col="T", event_col="E", show_progress=False)
            boot_hrs.append(float(np.exp(c.params_.loc["G"])))
        except Exception:
            continue

    if len(boot_hrs) >= 100:
        ci_lo = float(np.percentile(boot_hrs, 2.5))
        ci_hi = float(np.percentile(boot_hrs, 97.5))
        median_hr = float(np.median(boot_hrs))
    else:
        ci_lo = ci_hi = median_hr = np.nan

    try:
        lr = logrank_test(
            df.loc[df["G"] == 1, "T"],
            df.loc[df["G"] == 0, "T"],
            df.loc[df["G"] == 1, "E"],
            df.loc[df["G"] == 0, "E"],
        )
        p = float(lr.p_value)
    except Exception:
        p = np.nan

    return {
        "cell": cell,
        "reference": reference_cell,
        "n_target": n_target,
        "n_ref": n_ref,
        "events_target": e_target,
        "events_ref": e_ref,
        "hr_point": round(hr, 3) if not np.isnan(hr) else None,
        "hr_boot_median": round(median_hr, 3) if not np.isnan(median_hr) else None,
        "hr_boot_ci_lo": round(ci_lo, 3) if not np.isnan(ci_lo) else None,
        "hr_boot_ci_hi": round(ci_hi, 3) if not np.isnan(ci_hi) else None,
        "n_boot_valid": len(boot_hrs),
        "logrank_p": round(p, 4) if not np.isnan(p) else None,
        "ci_crosses_1": (ci_lo < 1 < ci_hi) if not np.isnan(ci_lo) else None,
    }


def tert_subgroup_breakdown(sm: pd.DataFrame) -> pd.DataFrame:
    """How are TERT+ patients distributed across drivers? This is the smoking
    gun for the 'TERT-only triple-neg-ish is worst' paradox."""
    surv = sm.dropna(subset=["os_event", "os_days"])
    tert_plus = surv[surv["tert_int"] == 1]
    rows = []
    for drv in ["BRAF", "RAS", "NTRK", "OTHER"]:
        sub = tert_plus[tert_plus["driver_simple"] == drv]
        rows.append(
            {
                "tert_status": "TERT+",
                "driver": drv,
                "n": len(sub),
                "events": int(sub["os_event"].sum()),
                "frac_of_TERT+": round(100 * len(sub) / len(tert_plus), 1) if len(tert_plus) else 0,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    sm = load_master()
    sm = make_cells(sm)

    crosstab = crosstab_summary(sm)
    crosstab.to_csv(OUT / "8cell_crosstab.tsv", sep="\t", index=False)
    print("8-cell crosstab:\n", crosstab.to_string(index=False))

    tert_subgroup = tert_subgroup_breakdown(sm)
    tert_subgroup.to_csv(OUT / "tert_subgroup_breakdown.tsv", sep="\t", index=False)
    print("\nTERT+ subgroup breakdown:\n", tert_subgroup.to_string(index=False))

    cells = sorted([c for c in sm["cell"].dropna().unique() if c != "OTHER_TERT-"])
    forest_rows = []
    for cell in cells:
        try:
            forest_rows.append(bootstrap_cox_hr(sm, cell, "OTHER_TERT-", n_boot=1000))
        except Exception as e:
            print(f"  cell={cell} failed: {e}")

    forest_df = pd.DataFrame(forest_rows)
    forest_df.to_csv(OUT / "forest_HR_8cell.tsv", sep="\t", index=False)
    print("\nForest HR (vs OTHER_TERT-):\n", forest_df.to_string(index=False))

    surv = sm.dropna(subset=["os_event", "os_days"]).copy()
    try:
        ml = multivariate_logrank_test(surv["os_days"], surv["cell"], surv["os_event"])
        omnibus_p = float(ml.p_value)
    except Exception:
        omnibus_p = np.nan

    summary = {
        "n_total_tumor": int(len(sm)),
        "n_with_survival": int(len(surv)),
        "n_total_events": int(surv["os_event"].sum()),
        "omnibus_logrank_p_8cell": round(omnibus_p, 6) if not np.isnan(omnibus_p) else None,
        "key_finding": _interpret(crosstab, tert_subgroup, forest_df),
        "files": {
            "8cell_crosstab": str(OUT / "8cell_crosstab.tsv"),
            "tert_subgroup_breakdown": str(OUT / "tert_subgroup_breakdown.tsv"),
            "forest_HR_8cell": str(OUT / "forest_HR_8cell.tsv"),
        },
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("\nSummary:", json.dumps(summary, indent=2))


def _interpret(crosstab, tert_sub, forest):
    msgs = []
    tplus = crosstab[crosstab["cell"].str.endswith("_TERT+")]
    if not tplus.empty:
        max_cell = tplus.loc[tplus["event_rate_pct"].idxmax()]
        msgs.append(
            f"Highest event-rate cell: {max_cell['cell']} "
            f"(N={max_cell['n']}, events={max_cell['events']}, rate={max_cell['event_rate_pct']}%)"
        )
    braf_tert = tert_sub.loc[tert_sub["driver"] == "BRAF"]
    if not braf_tert.empty:
        msgs.append(
            f"BRAF+TERT+ accounts for {braf_tert.iloc[0]['frac_of_TERT+']}% "
            f"of all TERT+ patients (N={braf_tert.iloc[0]['n']})"
        )
    if not forest.empty:
        ci_cross = forest[forest["ci_crosses_1"] == True]
        msgs.append(
            f"Cells with bootstrap CI crossing 1 (small-N artifact warning): "
            f"{', '.join(ci_cross['cell'].tolist()) if len(ci_cross) else 'none'}"
        )
    return msgs


if __name__ == "__main__":
    main()
