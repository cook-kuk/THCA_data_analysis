#!/usr/bin/env python3
"""Paper 11 — Phase B: DM1 axis × pan-cancer survival Cox.

Reads:
  - pancan_dm1_scored.tsv      (sample-level DM1_like, TF_collapse, RAI_8, lineage)
  - /data/thca/repo_data/raw/TCGA_pancan/survival.tsv  (TCGA harmonized survival)

Per-lineage Cox: OS ~ DM1_like (+ age, + stage where parsable). Also DSS, PFI.
Outputs:
  - phase_B_survival/cox_per_lineage.tsv  (HR, CI, p, n, events)
  - phase_B_survival/forest_plot.png
  - phase_B_survival/km_top_hits.png       (KM for top 4 by |HR|, FDR<0.1)
  - phase_B_survival/summary.json

Treats clinical claims as exploratory.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
import pandas as pd

P11 = Path("/data/thca/repo_results/paper11_pancancer")
SCORES = P11 / "pancan_dm1_scored.tsv"
SURV = Path("/data/thca/repo_data/raw/TCGA_pancan/survival.tsv")
OUT = P11 / "phase_B_survival"


def parse_stage(s):
    if pd.isna(s) or not isinstance(s, str): return np.nan
    s = s.lower().strip()
    if "stage iv" in s: return 4
    if "stage iii" in s: return 3
    if "stage ii" in s: return 2
    if "stage i" in s: return 1
    return np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scores", default=SCORES, type=Path)
    ap.add_argument("--surv", default=SURV, type=Path)
    ap.add_argument("--out-dir", default=OUT, type=Path)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if not args.scores.exists():
        sys.exit(f"missing {args.scores}")
    if not args.surv.exists():
        sys.exit(f"missing {args.surv}")

    scores = pd.read_csv(args.scores, sep="\t")
    surv = pd.read_csv(args.surv, sep="\t", low_memory=False)
    print(f"[B] scores: {len(scores)}, surv: {len(surv)}")

    # sample id harmonization. scores.sample looks like TCGA-XX-XXXX-01
    # surv.sample same format
    df = scores.merge(surv, left_on="sample", right_on="sample", how="left")
    df["age"] = pd.to_numeric(df.get("age_at_initial_pathologic_diagnosis"), errors="coerce")
    df["stage_int"] = df.get("ajcc_pathologic_tumor_stage").apply(parse_stage)
    print(f"[B] joined: {len(df)} rows; n with OS.time: {df['OS.time'].notna().sum()}")

    try:
        from lifelines import CoxPHFitter
    except ImportError:
        sys.exit("[B] lifelines not installed; pip install lifelines")

    rows = []
    km_data = {}
    for lin, sub in df.groupby("lineage"):
        for outcome, t_col, e_col in [("OS", "OS.time", "OS"),
                                       ("DSS", "DSS.time", "DSS"),
                                       ("PFI", "PFI.time", "PFI")]:
            cox = sub[[outcome, t_col, "DM1_like", "age", "stage_int"]].copy()
            cox = cox.dropna(subset=[outcome, t_col, "DM1_like"])
            cox = cox[cox[t_col] > 0]
            if len(cox) < 20 or cox[outcome].sum() < 5:
                continue
            # Drop covariates that are entirely missing
            covars = ["DM1_like"]
            if cox["age"].notna().sum() / len(cox) > 0.5:
                covars.append("age")
            if cox["stage_int"].notna().sum() / len(cox) > 0.3:
                covars.append("stage_int")
            cox = cox[[outcome, t_col] + covars].dropna()
            if len(cox) < 20: continue
            cph = CoxPHFitter(penalizer=0.001)
            try:
                cph.fit(cox, duration_col=t_col, event_col=outcome)
                s = cph.summary.loc["DM1_like"]
                rows.append({
                    "lineage": lin, "outcome": outcome,
                    "n": len(cox), "events": int(cox[outcome].sum()),
                    "covariates": "+".join(covars),
                    "HR": float(s["exp(coef)"]),
                    "HR_lower95": float(s["exp(coef) lower 95%"]),
                    "HR_upper95": float(s["exp(coef) upper 95%"]),
                    "p": float(s["p"]),
                    "concordance": float(cph.concordance_index_),
                })
            except Exception as e:
                rows.append({"lineage": lin, "outcome": outcome, "n": len(cox),
                             "events": int(cox[outcome].sum()),
                             "covariates": "+".join(covars),
                             "HR": np.nan, "p": np.nan, "error": str(e)[:80]})
        # save subset for KM
        km_data[lin] = sub

    res = pd.DataFrame(rows)
    # FDR per outcome
    from scipy import stats as st
    out_rows = []
    for outc, g in res.groupby("outcome"):
        g = g.copy().sort_values("p")
        valid = g["p"].notna()
        if valid.sum() == 0:
            out_rows.append(g); continue
        # BH
        m = valid.sum()
        ranks = np.arange(1, m + 1)
        pvals = g.loc[valid, "p"].values
        bh = pvals * m / ranks
        bh = np.minimum.accumulate(bh[::-1])[::-1]
        g.loc[valid, "fdr"] = bh
        out_rows.append(g)
    res = pd.concat(out_rows, ignore_index=True).sort_values(["outcome", "fdr"])
    out_path = args.out_dir / "cox_per_lineage.tsv"
    res.to_csv(out_path, sep="\t", index=False)
    print(f"[B] cox table → {out_path}: {len(res)} rows ({(res['p'].notna()).sum()} fitted)")

    # Headline: count significant (FDR<0.1) per outcome
    summary = {}
    for outc, g in res.groupby("outcome"):
        sig = g[(g["fdr"].notna()) & (g["fdr"] < 0.1)].copy()
        summary[outc] = {
            "n_lineages_fit": int(g["p"].notna().sum()),
            "n_significant_fdr10": int(len(sig)),
            "top_hits": sig.head(10)[["lineage", "HR", "p", "fdr"]].to_dict(orient="records")
        }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str)[:1500])

    # ---- forest plot ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        os = res[res.outcome == "OS"].dropna(subset=["HR"]).copy()
        os = os.sort_values("HR", ascending=True)
        fig, ax = plt.subplots(figsize=(7, max(5, 0.22 * len(os))))
        y = np.arange(len(os))
        ax.errorbar(os["HR"], y,
                    xerr=[os["HR"] - os["HR_lower95"], os["HR_upper95"] - os["HR"]],
                    fmt="o", capsize=3, color="black", ecolor="grey")
        for i, r in enumerate(os.itertuples()):
            mk = "**" if (r.fdr < 0.05) else ("*" if (r.fdr < 0.1) else "")
            ax.text(max(r.HR_upper95, r.HR) * 1.05, i, mk, va="center")
        ax.axvline(1, color="r", ls=":")
        ax.set_yticks(y); ax.set_yticklabels(os["lineage"], fontsize=7)
        ax.set_xscale("log")
        ax.set_xlabel("HR (DM1_like, OS Cox; * FDR<0.1, ** FDR<0.05)")
        ax.set_title("Pan-cancer DM1 axis prognostic effect (OS, Cox)")
        plt.tight_layout()
        fp = args.out_dir / "forest_plot_OS.png"
        plt.savefig(fp, dpi=140)
        plt.close()
        print(f"[B] forest → {fp}")
    except Exception as e:
        print(f"[B] forest plot failed: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
