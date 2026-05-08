#!/usr/bin/env python3
"""Track 18 — Pan-cancer ICI response x HLA-I/II module x DM1 score.

Builds on Paper 11 Phase C v2 outputs (4 ICI cohorts, n=421 with binary response).
This script consumes the existing per-sample signature scores and re-computes:

  1. Per-cohort HLA-I + HLA-II module score table (deliverable 1)
  2. HLA-I module x ICI response per cohort: OR (per +1 SD), 95% CI, p, FDR (deliverable 2)
  3. HLA-II module x ICI response per cohort: same (deliverable 3, NEW)
  4. DM1 x HLA-I x ICI response interaction (logistic with DM1, HLA-I, DM1xHLA-I) (deliverable 4)
  5. DM1-high stratified ICI response by HLA-I tertile (deliverable 5)
  6. Survival in ICI cohorts (HLA-I, HLA-II x OS Cox HR) (deliverable 6)
  7. Composite DM1_inflam re-confirm + extension with HLA-II (deliverable 7)
  8. Single-cell ICI cohorts: documented absent in repo (deliverable 8)

HLA gene-expression module ONLY -- not allele genotype.
"""

from __future__ import annotations
from pathlib import Path
import json
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.exceptions import ConvergenceError, ConvergenceWarning

warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track18_ici_hla")
TBL = OUT / "tables"
FIG = OUT / "figures"
LOG = OUT / "logs"
TBL.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
LOG.mkdir(parents=True, exist_ok=True)

PHASEC = Path("/home/seungho/personal/THCA_data_analysis/project/results/paper11_pancancer/phase_C_ICI")
SCORES_PATH = PHASEC / "results" / "tables" / "signature_scores_per_sample.tsv"

CAPTION = "HLA gene-expression module - not allele genotype."


def main() -> None:
    log_lines: list[str] = []

    def log(msg: str) -> None:
        print(msg)
        log_lines.append(msg)

    log(f"[track18] loading per-sample scores from {SCORES_PATH}")
    df = pd.read_csv(SCORES_PATH, sep="\t")
    log(f"[track18] {len(df)} samples x {df.shape[1]} cols")
    log(f"[track18] cohorts: {sorted(df['cohort'].unique())}")

    # Restrict to pre-treatment with binary response (Phase C v2 inclusion)
    pre_mask = df["timepoint"].fillna("pre").str.lower().isin(["pre", "baseline", "screening", "unknown"])
    # GSE176307 / IMvigor210 have no timepoint splits in the data; keep all when no 'on' label
    df_pre = df[pre_mask | df["timepoint"].isna()].copy()
    df_resp = df_pre.dropna(subset=["response_binary_CRPR_vs_SD_PD"]).copy()
    df_resp["response"] = df_resp["response_binary_CRPR_vs_SD_PD"].astype(int)
    log(f"[track18] pre + response: n={len(df_resp)} ({df_resp['cohort'].value_counts().to_dict()})")

    # ------------------------------------------------------------------
    # Deliverable 1: per-cohort HLA-I + HLA-II module score TSV
    # ------------------------------------------------------------------
    keep_cols = ["sample_id", "cohort", "cancer_type", "therapy", "timepoint",
                 "HLA_class_I", "HLA_class_II",
                 "IFNG_T_cell_inflamed", "checkpoint_exhaustion",
                 "myeloid_suppressive", "cytolytic_GZMA_PRF1",
                 "DM1_inflam_composite", "lineage_portable_DM1",
                 "thyroid_differentiation", "TLS_CXCL13_like",
                 "response_binary_CRPR_vs_SD_PD",
                 "os_time", "os_event", "pfs_time", "pfs_event"]
    keep_cols = [c for c in keep_cols if c in df.columns]
    out1 = df[keep_cols].copy()
    out1.to_csv(TBL / "T1_per_cohort_hla_dm1_scores.tsv", sep="\t", index=False)
    log(f"[T1] {len(out1)} samples written -> T1_per_cohort_hla_dm1_scores.tsv")

    # ------------------------------------------------------------------
    # Helper: per-cohort logistic regression of response ~ score (z-scored within cohort already in pipeline)
    # ------------------------------------------------------------------
    def per_cohort_or(score_col: str) -> pd.DataFrame:
        rows = []
        for cohort, sub in df_resp.groupby("cohort"):
            x = sub[score_col].astype(float).values
            y = sub["response"].astype(int).values
            if np.isnan(x).any():
                m = ~np.isnan(x)
                x, y = x[m], y[m]
            # z within cohort
            if x.std() > 0:
                xz = (x - x.mean()) / x.std()
            else:
                rows.append({"cohort": cohort, "score": score_col, "n": len(x),
                             "n_R": int(y.sum()), "n_NR": int(len(y) - y.sum()),
                             "logOR": np.nan, "se": np.nan,
                             "OR": np.nan, "lo": np.nan, "hi": np.nan, "p": np.nan,
                             "cohens_d": np.nan, "auc": np.nan})
                continue
            try:
                X = sm.add_constant(xz)
                model = sm.Logit(y, X).fit(disp=0, maxiter=100)
                logOR = float(model.params[1])
                se = float(model.bse[1])
                p = float(model.pvalues[1])
                lo = np.exp(logOR - 1.96 * se)
                hi = np.exp(logOR + 1.96 * se)
            except Exception as e:
                logOR = se = p = np.nan
                lo = hi = np.nan
            # cohen's d (R vs NR)
            xR = xz[y == 1]; xNR = xz[y == 0]
            if len(xR) > 1 and len(xNR) > 1:
                pooled = np.sqrt(((len(xR)-1)*xR.var(ddof=1) + (len(xNR)-1)*xNR.var(ddof=1)) / (len(xR)+len(xNR)-2))
                d = (xR.mean() - xNR.mean()) / pooled if pooled > 0 else np.nan
                # Mann-Whitney U as AUC
                try:
                    auc = stats.mannwhitneyu(xR, xNR, alternative="two-sided").statistic / (len(xR)*len(xNR))
                except Exception:
                    auc = np.nan
            else:
                d = np.nan; auc = np.nan
            rows.append({"cohort": cohort, "score": score_col, "n": len(x),
                         "n_R": int(y.sum()), "n_NR": int(len(y) - y.sum()),
                         "logOR": logOR, "se": se,
                         "OR": float(np.exp(logOR)) if np.isfinite(logOR) else np.nan,
                         "lo": lo, "hi": hi, "p": p,
                         "cohens_d": d, "auc": auc})
        out = pd.DataFrame(rows)
        # FDR within score
        if out["p"].notna().any():
            out["fdr"] = np.nan
            mask = out["p"].notna()
            out.loc[mask, "fdr"] = multipletests(out.loc[mask, "p"], method="fdr_bh")[1]
        return out

    # ------------------------------------------------------------------
    # Deliverable 2: HLA-I module x ICI response per cohort + meta
    # Deliverable 3: HLA-II module x ICI response per cohort + meta (NEW)
    # ------------------------------------------------------------------
    def ivw_meta(df_or: pd.DataFrame) -> dict:
        m = df_or.dropna(subset=["logOR", "se"])
        if len(m) == 0:
            return {"k": 0, "logOR": np.nan, "se": np.nan, "OR": np.nan, "p": np.nan, "n_pos": 0}
        w = 1.0 / (m["se"] ** 2)
        logOR = (m["logOR"] * w).sum() / w.sum()
        se = np.sqrt(1.0 / w.sum())
        p = 2 * (1 - stats.norm.cdf(abs(logOR / se)))
        return {"k": int(len(m)), "logOR": float(logOR), "se": float(se),
                "OR": float(np.exp(logOR)), "lo": float(np.exp(logOR - 1.96*se)),
                "hi": float(np.exp(logOR + 1.96*se)), "p": float(p),
                "n_pos": int((m["logOR"] > 0).sum()), "n_total": int(len(m))}

    score_panel = ["HLA_class_I", "HLA_class_II", "IFNG_T_cell_inflamed",
                   "checkpoint_exhaustion", "cytolytic_GZMA_PRF1",
                   "myeloid_suppressive", "DM1_inflam_composite",
                   "lineage_portable_DM1", "TLS_CXCL13_like"]

    all_or_rows = []
    meta_rows = []
    for sc in score_panel:
        df_or = per_cohort_or(sc)
        all_or_rows.append(df_or)
        m = ivw_meta(df_or)
        m["score"] = sc
        meta_rows.append(m)
    or_table = pd.concat(all_or_rows, ignore_index=True)
    or_table.to_csv(TBL / "T2_per_cohort_OR_response.tsv", sep="\t", index=False)
    log(f"[T2] per-cohort OR table written: {len(or_table)} rows ({or_table['score'].nunique()} scores x {or_table['cohort'].nunique()} cohorts)")

    meta_df = pd.DataFrame(meta_rows)[["score", "k", "n_pos", "n_total", "logOR", "se", "OR", "lo", "hi", "p"]]
    meta_df.to_csv(TBL / "T3_meta_OR_response.tsv", sep="\t", index=False)
    log(f"[T3] meta OR (IVW): {meta_df.to_string(index=False)}")

    # Phase C v2 reconfirm summary
    hla1_meta = meta_df.query("score=='HLA_class_I'").iloc[0]
    hla2_meta = meta_df.query("score=='HLA_class_II'").iloc[0]
    ifng_meta = meta_df.query("score=='IFNG_T_cell_inflamed'").iloc[0]
    chk_meta  = meta_df.query("score=='checkpoint_exhaustion'").iloc[0]
    cyt_meta  = meta_df.query("score=='cytolytic_GZMA_PRF1'").iloc[0]
    log(f"[reconfirm] Phase C v2: IFNG OR={ifng_meta['OR']:.3f} p={ifng_meta['p']:.2e}; HLA-I OR={hla1_meta['OR']:.3f} p={hla1_meta['p']:.3f}; checkpoint OR={chk_meta['OR']:.3f} p={chk_meta['p']:.3f}; cytolytic OR={cyt_meta['OR']:.3f} p={cyt_meta['p']:.3f}")
    log(f"[NEW HLA-II] OR={hla2_meta['OR']:.3f} 95%CI [{hla2_meta['lo']:.3f}, {hla2_meta['hi']:.3f}] p={hla2_meta['p']:.3f} k={int(hla2_meta['k'])} n_pos={int(hla2_meta['n_pos'])}/{int(hla2_meta['n_total'])}")

    # ------------------------------------------------------------------
    # Deliverable 4: DM1 x HLA-I x response logistic regression interaction
    # Per cohort + pooled
    # ------------------------------------------------------------------
    def fit_interaction(sub: pd.DataFrame) -> dict:
        out = {"n": len(sub)}
        sub = sub.dropna(subset=["DM1_inflam_composite", "HLA_class_I", "response"]).copy()
        if len(sub) < 12 or sub["response"].nunique() < 2:
            return {**out, "error": "insufficient data or single-class outcome"}
        # z-score within sub
        for c in ["DM1_inflam_composite", "HLA_class_I"]:
            v = sub[c].values
            sub[c+"_z"] = (v - v.mean()) / (v.std() if v.std() > 0 else 1.0)
        sub["DM1xHLA1"] = sub["DM1_inflam_composite_z"] * sub["HLA_class_I_z"]
        X = sm.add_constant(sub[["DM1_inflam_composite_z", "HLA_class_I_z", "DM1xHLA1"]])
        y = sub["response"].astype(int).values
        try:
            model = sm.Logit(y, X).fit(disp=0, maxiter=200)
            return {**out,
                    "beta_DM1": float(model.params["DM1_inflam_composite_z"]),
                    "p_DM1":    float(model.pvalues["DM1_inflam_composite_z"]),
                    "beta_HLA1":float(model.params["HLA_class_I_z"]),
                    "p_HLA1":   float(model.pvalues["HLA_class_I_z"]),
                    "beta_int": float(model.params["DM1xHLA1"]),
                    "p_int":    float(model.pvalues["DM1xHLA1"]),
                    "OR_int":   float(np.exp(model.params["DM1xHLA1"])),
                    "se_int":   float(model.bse["DM1xHLA1"]),
                    "lo_int":   float(np.exp(model.params["DM1xHLA1"] - 1.96*model.bse["DM1xHLA1"])),
                    "hi_int":   float(np.exp(model.params["DM1xHLA1"] + 1.96*model.bse["DM1xHLA1"])),
                    "llf":      float(model.llf), "aic": float(model.aic)}
        except Exception as e:
            return {**out, "error": str(e)}

    inter_rows = []
    for cohort, sub in df_resp.groupby("cohort"):
        r = fit_interaction(sub)
        r["cohort"] = cohort
        inter_rows.append(r)
    inter_rows.append({**fit_interaction(df_resp), "cohort": "POOLED"})
    inter_df = pd.DataFrame(inter_rows)
    cols = ["cohort", "n", "beta_DM1", "p_DM1", "beta_HLA1", "p_HLA1",
            "beta_int", "se_int", "OR_int", "lo_int", "hi_int", "p_int", "llf", "aic", "error"]
    cols = [c for c in cols if c in inter_df.columns]
    inter_df = inter_df[cols]
    inter_df.to_csv(TBL / "T4_dm1_x_hla1_interaction.tsv", sep="\t", index=False)
    log(f"[T4] DM1 x HLA-I interaction (per-cohort + pooled):\n{inter_df.to_string(index=False)}")

    # ------------------------------------------------------------------
    # Deliverable 5: within-DM1-high samples, HLA-I tertile -> response
    # ------------------------------------------------------------------
    strat_rows = []
    for cohort, sub in df_resp.groupby("cohort"):
        sub = sub.dropna(subset=["DM1_inflam_composite", "HLA_class_I", "response"]).copy()
        if len(sub) < 12:
            continue
        # Top tertile of DM1
        dm1_q = sub["DM1_inflam_composite"].quantile([1/3, 2/3]).values
        sub["dm1_strata"] = pd.cut(sub["DM1_inflam_composite"], bins=[-np.inf, dm1_q[0], dm1_q[1], np.inf],
                                    labels=["low", "mid", "high"])
        for strat in ["high", "low"]:
            sg = sub[sub["dm1_strata"] == strat].copy()
            if len(sg) < 6 or sg["response"].nunique() < 2:
                strat_rows.append({"cohort": cohort, "dm1_strata": strat, "n": len(sg),
                                   "n_R": int(sg["response"].sum()),
                                   "auc_hla1": np.nan, "logOR_hla1": np.nan,
                                   "p_hla1": np.nan})
                continue
            xR = sg.loc[sg["response"]==1, "HLA_class_I"]
            xNR = sg.loc[sg["response"]==0, "HLA_class_I"]
            try:
                auc = stats.mannwhitneyu(xR, xNR, alternative="two-sided").statistic / (len(xR)*len(xNR))
            except Exception:
                auc = np.nan
            x = sg["HLA_class_I"].values
            xz = (x - x.mean()) / (x.std() if x.std() > 0 else 1.0)
            try:
                Xm = sm.add_constant(xz)
                m = sm.Logit(sg["response"].astype(int).values, Xm).fit(disp=0, maxiter=100)
                logOR = float(m.params[1]); p = float(m.pvalues[1])
            except Exception:
                logOR = np.nan; p = np.nan
            strat_rows.append({"cohort": cohort, "dm1_strata": strat, "n": len(sg),
                               "n_R": int(sg["response"].sum()),
                               "auc_hla1": auc, "logOR_hla1": logOR,
                               "OR_hla1": float(np.exp(logOR)) if np.isfinite(logOR) else np.nan,
                               "p_hla1": p})
    strat_df = pd.DataFrame(strat_rows)
    strat_df.to_csv(TBL / "T5_within_DM1_strata_HLA1_response.tsv", sep="\t", index=False)
    log(f"[T5] within-DM1-stratum HLA-I -> response (DM1-high vs DM1-low rows):\n{strat_df.to_string(index=False)}")

    # Pool DM1-high logORs
    dm1_high = strat_df[strat_df["dm1_strata"]=="high"].dropna(subset=["logOR_hla1"])
    if len(dm1_high) > 0:
        # need SE; refit pooled
        sub_all = df_resp.dropna(subset=["DM1_inflam_composite","HLA_class_I","response"]).copy()
        sub_all["dm1_high_global"] = (sub_all.groupby("cohort")["DM1_inflam_composite"]
                                      .transform(lambda x: x >= x.quantile(2/3)))
        sg = sub_all[sub_all["dm1_high_global"]].copy()
        if len(sg) > 10 and sg["response"].nunique() > 1:
            sg["hla1_z"] = sg.groupby("cohort")["HLA_class_I"].transform(lambda x: (x-x.mean())/x.std())
            sg["cohort_id"] = sg["cohort"].astype("category").cat.codes
            X = pd.get_dummies(sg["cohort"], drop_first=True).astype(float)
            X["hla1_z"] = sg["hla1_z"].values
            X = sm.add_constant(X)
            try:
                pooled_m = sm.Logit(sg["response"].astype(int).values, X.values.astype(float)).fit(disp=0, maxiter=200)
                idx = list(X.columns).index("hla1_z")
                hla1_logOR = float(pooled_m.params[idx])
                hla1_se = float(pooled_m.bse[idx])
                hla1_p = float(pooled_m.pvalues[idx])
                log(f"[T5 pooled DM1-high] HLA-I logOR_per_z={hla1_logOR:.3f} se={hla1_se:.3f} p={hla1_p:.4f} OR={np.exp(hla1_logOR):.3f}")
            except Exception as e:
                log(f"[T5 pooled DM1-high] fit failed: {e}")

    # ------------------------------------------------------------------
    # Deliverable 6: Survival -- HLA-I, HLA-II x OS Cox HR per cohort
    # ------------------------------------------------------------------
    surv_rows = []
    for cohort, sub in df_pre.groupby("cohort"):
        if "os_time" not in sub.columns or "os_event" not in sub.columns:
            continue
        sub = sub.dropna(subset=["os_time", "os_event"]).copy()
        if len(sub) < 12 or sub["os_event"].sum() < 5:
            continue
        for sc in ["HLA_class_I", "HLA_class_II", "DM1_inflam_composite",
                   "IFNG_T_cell_inflamed", "checkpoint_exhaustion",
                   "cytolytic_GZMA_PRF1", "lineage_portable_DM1"]:
            ss = sub.dropna(subset=[sc]).copy()
            if len(ss) < 12:
                continue
            v = ss[sc].values
            ss["score_z"] = (v - v.mean()) / (v.std() if v.std() > 0 else 1.0)
            cph = CoxPHFitter()
            try:
                cph.fit(ss[["os_time", "os_event", "score_z"]],
                        duration_col="os_time", event_col="os_event")
                hr = float(np.exp(cph.params_["score_z"]))
                lo = float(np.exp(cph.confidence_intervals_.iloc[0,0]))
                hi = float(np.exp(cph.confidence_intervals_.iloc[0,1]))
                p = float(cph.summary["p"].iloc[0])
            except (ConvergenceError, Exception):
                hr = lo = hi = p = np.nan
            surv_rows.append({"cohort": cohort, "score": sc, "n": len(ss),
                              "n_event": int(ss["os_event"].sum()),
                              "HR_per_z": hr, "HR_lo": lo, "HR_hi": hi, "p": p})
    surv_df = pd.DataFrame(surv_rows)
    surv_df.to_csv(TBL / "T6_survival_cox_per_cohort.tsv", sep="\t", index=False)
    log(f"[T6] survival cox: {len(surv_df)} rows ({surv_df.dropna(subset=['HR_per_z'])['cohort'].nunique()} cohorts with OS)")

    # ------------------------------------------------------------------
    # Deliverable 7: composite DM1_inflam re-confirmation + HLA-II augmented composite
    # ------------------------------------------------------------------
    df_resp["DM1_inflam_plus_HLA2"] = df_resp[["IFNG_T_cell_inflamed", "myeloid_suppressive",
                                                 "checkpoint_exhaustion", "HLA_class_II",
                                                 "HLA_class_I"]].mean(axis=1)
    df_resp["DM1_HLA2only"] = df_resp[["IFNG_T_cell_inflamed", "myeloid_suppressive",
                                         "checkpoint_exhaustion", "HLA_class_II"]].mean(axis=1)

    composite_rows = []
    for sc in ["DM1_inflam_composite", "DM1_HLA2only", "DM1_inflam_plus_HLA2"]:
        sign_pos = 0
        cohorts_used = 0
        per_cohort = []
        for cohort, sub in df_resp.groupby("cohort"):
            x = sub[sc].astype(float).values; y = sub["response"].astype(int).values
            m = ~np.isnan(x); x = x[m]; y = y[m]
            if x.std() == 0 or len(x) < 6:
                continue
            xz = (x - x.mean()) / x.std()
            try:
                model = sm.Logit(y, sm.add_constant(xz)).fit(disp=0, maxiter=100)
                logOR = float(model.params[1]); se = float(model.bse[1]); p = float(model.pvalues[1])
            except Exception:
                continue
            cohorts_used += 1
            if logOR > 0:
                sign_pos += 1
            per_cohort.append({"cohort": cohort, "score": sc, "logOR": logOR, "se": se, "p": p})
        meta = ivw_meta(pd.DataFrame(per_cohort)) if per_cohort else None
        if meta is not None:
            composite_rows.append({"score": sc, "k_cohorts": cohorts_used,
                                   "sign_positive": sign_pos,
                                   "pooled_OR": meta["OR"], "lo": meta["lo"], "hi": meta["hi"],
                                   "p": meta["p"]})
    comp_df = pd.DataFrame(composite_rows)
    comp_df.to_csv(TBL / "T7_composite_DM1_inflam_extended.tsv", sep="\t", index=False)
    log(f"[T7] composite extended:\n{comp_df.to_string(index=False)}")

    # ------------------------------------------------------------------
    # FIGURES
    # ------------------------------------------------------------------
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})

    # Figure 1: HLA-I module x response forest (per cohort + meta)
    def forest(score: str, fname: str, title: str) -> None:
        sub = or_table[or_table["score"] == score].copy().reset_index(drop=True)
        meta = meta_df[meta_df["score"] == score].iloc[0]
        fig, ax = plt.subplots(figsize=(7.5, 0.5*len(sub)+2.0))
        ys = np.arange(len(sub)+1)[::-1]
        for i, row in sub.iterrows():
            ax.errorbar(row["OR"], ys[i],
                        xerr=[[row["OR"]-row["lo"]], [row["hi"]-row["OR"]]],
                        fmt="o", color="#1f77b4", capsize=3, lw=1.4, ms=6)
            ax.text(0.02, ys[i], f"{row['cohort']}  (n={int(row['n'])}, R={int(row['n_R'])}/NR={int(row['n_NR'])})",
                    transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=8)
            ax.text(1.02, ys[i], f"OR={row['OR']:.2f} [{row['lo']:.2f},{row['hi']:.2f}]  p={row['p']:.2g}",
                    transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=8)
        # meta diamond
        ax.errorbar(meta["OR"], ys[-1],
                    xerr=[[meta["OR"]-meta["lo"]], [meta["hi"]-meta["OR"]]],
                    fmt="D", color="#d62728", capsize=3, lw=1.6, ms=8)
        ax.text(0.02, ys[-1], f"IVW META (k={int(meta['k'])})",
                transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=8, fontweight="bold")
        ax.text(1.02, ys[-1], f"OR={meta['OR']:.2f} [{meta['lo']:.2f},{meta['hi']:.2f}]  p={meta['p']:.2g}",
                transform=ax.get_yaxis_transform(), ha="left", va="center",
                fontsize=8, fontweight="bold")
        ax.axvline(1.0, color="grey", ls="--", lw=0.8)
        ax.set_xscale("log")
        ax.set_xlim(0.3, 3.0)
        ax.set_yticks([])
        ax.set_xlabel("OR per +1 SD (responder vs non-responder)")
        ax.set_title(f"{title}\n{CAPTION}", fontsize=10)
        plt.tight_layout()
        fig.savefig(FIG / fname, dpi=170, bbox_inches="tight")
        plt.close(fig)

    forest("HLA_class_I",  "F1_forest_HLA1_response.png",  "Track 18 - HLA-I module x ICI response (per cohort + IVW meta)")
    forest("HLA_class_II", "F2_forest_HLA2_response.png",  "Track 18 - HLA-II module x ICI response (per cohort + IVW meta)")
    forest("DM1_inflam_composite", "F3_forest_DM1_inflam_response.png", "Track 18 - DM1_inflam composite x ICI response")

    # Figure 4: meta OR heatmap (rows scores, cols cohorts) -- log10 OR
    pivot_logOR = or_table.pivot(index="score", columns="cohort", values="logOR")
    pivot_p = or_table.pivot(index="score", columns="cohort", values="p")
    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(pivot_logOR.values, cmap="RdBu_r", vmin=-1.0, vmax=1.0, aspect="auto")
    ax.set_xticks(range(pivot_logOR.shape[1])); ax.set_xticklabels(pivot_logOR.columns, rotation=30, ha="right")
    ax.set_yticks(range(pivot_logOR.shape[0])); ax.set_yticklabels(pivot_logOR.index)
    for i in range(pivot_logOR.shape[0]):
        for j in range(pivot_logOR.shape[1]):
            v = pivot_logOR.values[i,j]; p = pivot_p.values[i,j]
            star = "*" if (np.isfinite(p) and p < 0.05) else ""
            if np.isfinite(v):
                ax.text(j, i, f"{np.exp(v):.2f}{star}", ha="center", va="center", fontsize=8,
                        color="white" if abs(v) > 0.5 else "black")
    plt.colorbar(im, ax=ax, label="logOR per +1 SD (response)")
    ax.set_title(f"Track 18 - per-cohort OR heatmap (response, * p<0.05)\n{CAPTION}")
    plt.tight_layout()
    fig.savefig(FIG / "F4_OR_heatmap_per_cohort.png", dpi=170, bbox_inches="tight")
    plt.close(fig)

    # Figure 5: DM1 x HLA-I interaction scatter (pooled, color by response)
    sub_all = df_resp.dropna(subset=["DM1_inflam_composite","HLA_class_I","response"]).copy()
    sub_all["dm1_z"] = sub_all.groupby("cohort")["DM1_inflam_composite"].transform(lambda x: (x-x.mean())/x.std())
    sub_all["hla1_z"] = sub_all.groupby("cohort")["HLA_class_I"].transform(lambda x: (x-x.mean())/x.std())
    fig, axes = plt.subplots(1, len(sub_all["cohort"].unique()), figsize=(4*sub_all["cohort"].nunique(), 4),
                              sharex=True, sharey=True)
    if sub_all["cohort"].nunique() == 1:
        axes = [axes]
    for ax, (cohort, sg) in zip(axes, sub_all.groupby("cohort")):
        for resp_lab, color, marker in [(1, "#2ca02c", "o"), (0, "#d62728", "x")]:
            ss = sg[sg["response"] == resp_lab]
            ax.scatter(ss["dm1_z"], ss["hla1_z"], c=color, marker=marker, s=26,
                       alpha=0.75, label=f"{'R' if resp_lab==1 else 'NR'} (n={len(ss)})")
        ax.axhline(0, color="grey", ls="--", lw=0.6); ax.axvline(0, color="grey", ls="--", lw=0.6)
        ax.set_title(cohort, fontsize=9)
        ax.set_xlabel("DM1_inflam (z)"); ax.set_ylabel("HLA-I module (z)")
        ax.legend(loc="upper left", fontsize=7)
    fig.suptitle(f"Track 18 - DM1 x HLA-I per cohort, by response\n{CAPTION}", fontsize=10)
    plt.tight_layout()
    fig.savefig(FIG / "F5_DM1_x_HLA1_scatter_by_response.png", dpi=170, bbox_inches="tight")
    plt.close(fig)

    # Figure 6: within-DM1-high stratum, HLA-I distribution by response
    fig, axes = plt.subplots(1, sub_all["cohort"].nunique(), figsize=(4*sub_all["cohort"].nunique(), 4),
                              sharey=True)
    if sub_all["cohort"].nunique() == 1:
        axes = [axes]
    for ax, (cohort, sg) in zip(axes, sub_all.groupby("cohort")):
        thr = sg["DM1_inflam_composite"].quantile(2/3)
        sg_high = sg[sg["DM1_inflam_composite"] >= thr]
        data_R = sg_high.loc[sg_high["response"]==1, "HLA_class_I"]
        data_NR = sg_high.loc[sg_high["response"]==0, "HLA_class_I"]
        bp = ax.boxplot([data_NR.values, data_R.values], labels=["NR","R"], widths=0.5, patch_artist=True)
        for patch, c in zip(bp["boxes"], ["#d62728","#2ca02c"]):
            patch.set_facecolor(c); patch.set_alpha(0.5)
        ax.scatter(np.random.normal(1, 0.06, len(data_NR)), data_NR.values, color="#d62728", s=18, alpha=0.7)
        ax.scatter(np.random.normal(2, 0.06, len(data_R)),  data_R.values, color="#2ca02c", s=18, alpha=0.7)
        ax.set_title(f"{cohort}\nDM1-high (top tertile) n={len(sg_high)}", fontsize=9)
        ax.set_ylabel("HLA-I module score")
    fig.suptitle(f"Track 18 - within DM1-high: HLA-I module by response\n{CAPTION}", fontsize=10)
    plt.tight_layout()
    fig.savefig(FIG / "F6_within_DM1_high_HLA1_by_response.png", dpi=170, bbox_inches="tight")
    plt.close(fig)

    # Figure 7: survival forest -- HLA-I and HLA-II x OS, per cohort
    if len(surv_df.dropna(subset=["HR_per_z"])) > 0:
        s2 = surv_df[surv_df["score"].isin(["HLA_class_I","HLA_class_II","DM1_inflam_composite",
                                             "IFNG_T_cell_inflamed","checkpoint_exhaustion"])].dropna(subset=["HR_per_z"])
        fig, ax = plt.subplots(figsize=(7.5, 0.4*len(s2)+1.5))
        ys = np.arange(len(s2))[::-1]
        colors = {"HLA_class_I":"#1f77b4","HLA_class_II":"#9467bd",
                  "DM1_inflam_composite":"#ff7f0e","IFNG_T_cell_inflamed":"#2ca02c",
                  "checkpoint_exhaustion":"#8c564b"}
        for (i, (_, row)), y in zip(enumerate(s2.iterrows()), ys):
            row = row[1] if isinstance(row, tuple) else s2.iloc[i]
        for i in range(len(s2)):
            row = s2.iloc[i]
            ax.errorbar(row["HR_per_z"], ys[i],
                        xerr=[[row["HR_per_z"]-row["HR_lo"]], [row["HR_hi"]-row["HR_per_z"]]],
                        fmt="o", color=colors.get(row["score"], "black"), capsize=3, lw=1.3, ms=6)
            ax.text(0.02, ys[i], f"{row['cohort']}  {row['score']}",
                    transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=8)
            ax.text(1.02, ys[i], f"HR={row['HR_per_z']:.2f} [{row['HR_lo']:.2f},{row['HR_hi']:.2f}] p={row['p']:.2g}",
                    transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=8)
        ax.axvline(1.0, color="grey", ls="--", lw=0.8); ax.set_xscale("log"); ax.set_xlim(0.3, 3)
        ax.set_yticks([])
        ax.set_xlabel("HR per +1 SD (overall survival)")
        ax.set_title(f"Track 18 - survival in ICI cohorts: HLA-I, HLA-II, DM1 modules\n{CAPTION}", fontsize=10)
        plt.tight_layout()
        fig.savefig(FIG / "F7_survival_forest.png", dpi=170, bbox_inches="tight")
        plt.close(fig)

    # ------------------------------------------------------------------
    # Single-cell ICI cohort discovery (deliverable 8)
    # ------------------------------------------------------------------
    sc_rows = []
    candidates = ["sade-feldman", "GSE120575", "GSE123813", "GSE166181", "Yost", "Bassez"]
    for c in candidates:
        sc_rows.append({"candidate": c, "found_in_repo": False,
                        "note": "no scRNA ICI cohort in /data/thca/repo_results/paper11_pancancer or project/data/external"})
    pd.DataFrame(sc_rows).to_csv(TBL / "T8_singlecell_ici_audit.tsv", sep="\t", index=False)
    log("[T8] single-cell ICI: no scRNA ICI cohort located in repo (documented)")

    # ------------------------------------------------------------------
    # Save reconfirmation summary JSON
    # ------------------------------------------------------------------
    summary = {
        "phase_C_v2_reconfirm": {
            "IFNG":         {"OR": ifng_meta["OR"], "p": ifng_meta["p"], "k_pos": int(ifng_meta["n_pos"])},
            "HLA_I":        {"OR": hla1_meta["OR"], "p": hla1_meta["p"], "k_pos": int(hla1_meta["n_pos"])},
            "checkpoint":   {"OR": chk_meta["OR"],  "p": chk_meta["p"],  "k_pos": int(chk_meta["n_pos"])},
            "cytolytic":    {"OR": cyt_meta["OR"],  "p": cyt_meta["p"],  "k_pos": int(cyt_meta["n_pos"])},
        },
        "new_HLA_II": {
            "OR_meta": hla2_meta["OR"], "p_meta": hla2_meta["p"],
            "lo": hla2_meta["lo"], "hi": hla2_meta["hi"],
            "k_cohorts": int(hla2_meta["k"]),
            "n_pos_direction": int(hla2_meta["n_pos"])
        },
        "interaction_pooled":  inter_df[inter_df["cohort"]=="POOLED"].to_dict("records")[0],
        "n_cohorts_total": int(or_table["cohort"].nunique()),
        "n_samples_with_response": int(len(df_resp)),
        "boundary": CAPTION
    }
    (OUT / "track18_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    log(f"[summary] written -> track18_summary.json")

    (LOG / "track18_run.log").write_text("\n".join(log_lines))
    print(f"\n[track18] DONE. outputs at {OUT}")


if __name__ == "__main__":
    main()
