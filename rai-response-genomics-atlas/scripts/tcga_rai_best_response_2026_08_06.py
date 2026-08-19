#!/usr/bin/env python3
"""TCGA-THCA — 8-gene differentiation panel vs actual radioiodine treatment response.

Why this exists
---------------
The project has repeatedly stated that "no cohort in this study contains RAI outcome
data linked to per-sample DM1 calls". That is wrong. The GDC BCR Biotab file
`nationwidechildrens.org_clinical_radiation_thca.txt` records, per radiation course:

  radiation_adjuvant_units = mCi   (249 courses -> radioiodine, not external beam)
  radiation_total_dose             (cumulative activity in mCi)
  treatment_best_response          (Complete Response / Partial Response /
                                    Stable Disease / Radiographic Progressive Disease)

These fields are absent from every cBioPortal THCA study, which is why they were missed.
They give an RAI *response* label on the same patients for whom we already have the
8-gene panel score — i.e. the missing link between the molecular axis and RAI effect,
without invoking survival.

Endpoint choice
---------------
Primary endpoint is failure to achieve complete response to RAI (non-CR vs CR), NOT
overall survival: TCGA-THCA has 16 deaths / 500, which cannot support a survival claim
(and the ATA 2015 evidence base itself establishes an RAI survival benefit only for
T4 gross extrathyroidal extension and M1 disease, on observational data alone).

Outputs (new files only):
  results/tables/tcga_rai_response_{main,adjusted,dose,pfi,power}_2026_08_06.tsv
  results/figures/figure_tcga_rai_best_response_2026_08_06.{png,pdf}
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
FIG = ROOT / "results" / "figures"
TAB = ROOT / "results" / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

REPO = Path("/data/thca/repo_results")
PANEL = REPO / "r17_tcga_panel_d4p2_reconciliation" / "r17_per_sample_merged.tsv"
PURITY = REPO / "p2_braf_nature_sprint_2026_05_09" / "r6_purity_audit" / "r6_purity_per_sample.tsv"
CBIO = (REPO / "p2_braf_nature_sprint_2026_05_09" / "h24_survival_sensitivity" / "cache"
        / "cbio_thca_tcga_patient_clinical.tsv")
SURV = Path("/data/thca/repo_data/raw/TCGA_pancan/survival.tsv")
RADIO = DATA / "raw" / "TCGA_THCA_biotab" / "nationwidechildrens.org_clinical_radiation_thca.txt"

SEED = 20260806
STAMP = "2026_08_06"
RESP_ORDER = {"Complete Response": 0, "Partial Response": 1, "Stable Disease": 2,
              "Radiographic Progressive Disease": 3}


def cohens_d(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return (a.mean() - b.mean()) / sp if sp > 0 else np.nan


def boot_ci_d(a, b, n=5000, seed=SEED):
    rng = np.random.default_rng(seed)
    a, b = np.asarray(a, float), np.asarray(b, float)
    v = [cohens_d(rng.choice(a, len(a), True), rng.choice(b, len(b), True)) for _ in range(n)]
    v = np.asarray([x for x in v if np.isfinite(x)])
    return (np.percentile(v, 2.5), np.percentile(v, 97.5)) if len(v) else (np.nan, np.nan)


def auc_mw(pos, neg):
    if len(pos) < 1 or len(neg) < 1:
        return np.nan
    u = stats.mannwhitneyu(pos, neg, alternative="two-sided").statistic
    return u / (len(pos) * len(neg))


def build():
    panel = pd.read_csv(PANEL, sep="\t")
    panel["patient"] = panel["sampleId"].str.slice(0, 12)
    panel = panel.drop_duplicates("patient")

    rad = pd.read_csv(RADIO, sep="\t", skiprows=[1, 2], low_memory=False)
    rad["dose"] = pd.to_numeric(rad["radiation_total_dose"], errors="coerce")
    # radioiodine courses are the ones dosed in millicuries; Gy/cGy = external beam
    rai = rad[rad["radiation_adjuvant_units"].astype(str).str.lower().eq("mci")].copy()
    ebrt_pts = set(rad.loc[rad["radiation_adjuvant_units"].astype(str).str.lower()
                           .isin(["gy", "cgy"]), "bcr_patient_barcode"])

    rai["resp_rank"] = rai["treatment_best_response"].map(RESP_ORDER)
    agg = rai.groupby("bcr_patient_barcode").agg(
        rai_courses=("bcr_radiation_barcode", "count"),
        rai_dose_mci=("dose", "sum"),
        best_rank=("resp_rank", "min"),
        worst_rank=("resp_rank", "max"),
    ).reset_index().rename(columns={"bcr_patient_barcode": "patient"})
    inv = {v: k for k, v in RESP_ORDER.items()}
    agg["best_response"] = agg["best_rank"].map(inv)
    agg["worst_response"] = agg["worst_rank"].map(inv)
    agg["had_ebrt"] = agg["patient"].isin(ebrt_pts)

    df = panel.merge(agg, on="patient", how="inner")

    # This audit table stores composition as leukocyte fraction, not a "purity" column;
    # leuko_frac is the inverse-purity proxy and is the covariate that matters here,
    # because the eight panel genes are thyrocyte-specific.
    pur = pd.read_csv(PURITY, sep="\t")
    pcol = next((c for c in pur.columns if "purity" in c.lower()), None) or "leuko_frac"
    pur["patient"] = pur["sample_id"].astype(str).str.slice(0, 12)
    pur[pcol] = pd.to_numeric(pur[pcol], errors="coerce")
    pur = pur.dropna(subset=[pcol]).drop_duplicates("patient")
    df = df.merge(pur[["patient", pcol]].rename(columns={pcol: "purity"}), on="patient", how="left")
    print(f"composition covariate = {pcol}; non-missing {df['purity'].notna().sum()}/{len(df)}")

    cb = pd.read_csv(CBIO, sep="\t", low_memory=False)
    keep = ["patientId", "AGE", "AJCC_PATHOLOGIC_TUMOR_STAGE", "EXTRATHYROIDAL_EXTENSION", "SEX"]
    cb = cb[[c for c in keep if c in cb.columns]].rename(columns={"patientId": "patient"})
    df = df.merge(cb, on="patient", how="left")

    sv = pd.read_csv(SURV, sep="\t", low_memory=False)
    sv = sv[sv["cancer type abbreviation"] == "THCA"]
    sv["patient"] = sv["_PATIENT"]
    df = df.merge(sv[["patient", "PFI", "PFI.time", "OS", "OS.time"]].drop_duplicates("patient"),
                  on="patient", how="left")
    return df


def main():
    df = build()
    df["evaluable"] = df["best_response"].notna()
    ev = df[df["evaluable"]].copy()
    ev["nonCR"] = (ev["best_response"] != "Complete Response").astype(int)

    print(f"RAI-treated patients with panel score: {len(df)}")
    print(f"  with evaluable best response       : {len(ev)}")
    print("  response distribution:", ev["best_response"].value_counts().to_dict())
    print("  cumulative mCi: median %.1f (IQR %.1f-%.1f)" % (
        df["rai_dose_mci"].median(), df["rai_dose_mci"].quantile(.25), df["rai_dose_mci"].quantile(.75)))

    cr = ev.loc[ev.nonCR == 0, "RAI_8"].dropna().values
    ncr = ev.loc[ev.nonCR == 1, "RAI_8"].dropna().values
    d = cohens_d(ncr, cr)
    lo, hi = boot_ci_d(ncr, cr)
    p = stats.mannwhitneyu(ncr, cr, alternative="two-sided").pvalue

    main_rows = [dict(
        analysis="PRIMARY: panel z, non-CR vs CR (RAI-treated)", n_nonCR=len(ncr), n_CR=len(cr),
        median_nonCR=float(np.median(ncr)), median_CR=float(np.median(cr)),
        cohens_d=float(d), d_ci_low=float(lo), d_ci_high=float(hi),
        auc=float(auc_mw(ncr, cr)), p=float(p))]

    # sensitivity: drop anyone who also received external-beam radiation
    pure = ev[~ev.had_ebrt]
    a = pure.loc[pure.nonCR == 1, "RAI_8"].dropna().values
    b = pure.loc[pure.nonCR == 0, "RAI_8"].dropna().values
    if len(a) >= 2 and len(b) >= 2:
        lo2, hi2 = boot_ci_d(a, b)
        main_rows.append(dict(analysis="sensitivity: RAI only, no external beam",
                              n_nonCR=len(a), n_CR=len(b),
                              median_nonCR=float(np.median(a)), median_CR=float(np.median(b)),
                              cohens_d=float(cohens_d(a, b)), d_ci_low=float(lo2), d_ci_high=float(hi2),
                              auc=float(auc_mw(a, b)),
                              p=float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)))

    # DM call contingency
    ct = pd.crosstab(ev["DM_call"], ev["nonCR"])
    if ct.shape == (2, 2):
        orr, pf = stats.fisher_exact(ct.values)
        main_rows.append(dict(analysis=f"DM1 vs DM2 x non-CR (Fisher) {ct.values.tolist()}",
                              n_nonCR=int(ct[1].sum()), n_CR=int(ct[0].sum()),
                              median_nonCR=np.nan, median_CR=np.nan, cohens_d=float(orr),
                              d_ci_low=np.nan, d_ci_high=np.nan, auc=np.nan, p=float(pf)))
    main_df = pd.DataFrame(main_rows)
    main_df.to_csv(TAB / f"tcga_rai_response_main_{STAMP}.tsv", sep="\t", index=False)
    print("\n", main_df.to_string(index=False))

    # adjusted logistic regression
    try:
        import statsmodels.formula.api as smf
        m = ev.copy()
        m["stage_high"] = m["AJCC_PATHOLOGIC_TUMOR_STAGE"].astype(str).str.contains(
            "III|IV", regex=True).astype(int)
        m["braf"] = m["has_braf_v600e"].astype(str).str.lower().eq("true").astype(int)
        cols = ["RAI_8", "nonCR", "stage_high", "braf", "AGE", "purity", "rai_dose_mci"]
        m = m[cols].apply(pd.to_numeric, errors="coerce").dropna()
        f = "nonCR ~ RAI_8 + stage_high + braf + AGE + purity + rai_dose_mci"
        fit = smf.logit(f, data=m).fit(disp=0)
        adj = pd.DataFrame({"term": fit.params.index, "beta": fit.params.values,
                            "or": np.exp(fit.params.values),
                            "ci_low": np.exp(fit.conf_int()[0].values),
                            "ci_high": np.exp(fit.conf_int()[1].values),
                            "p": fit.pvalues.values})
        adj["n_used"] = len(m)
        adj.to_csv(TAB / f"tcga_rai_response_adjusted_{STAMP}.tsv", sep="\t", index=False)
        print("\nAdjusted logistic (non-CR):\n", adj.to_string(index=False))
    except Exception as exc:  # noqa: BLE001
        print(f"adjusted model skipped: {exc}")
        adj = pd.DataFrame()

    # dose association
    dd = df[["RAI_8", "rai_dose_mci"]].dropna()
    rho, prho = stats.spearmanr(dd["RAI_8"], dd["rai_dose_mci"])
    dose = pd.DataFrame([dict(analysis="panel z vs cumulative RAI activity (mCi)", n=len(dd),
                              spearman_rho=float(rho), p=float(prho),
                              note="dose is prescribed by stage/risk: confounded by indication")])
    dose.to_csv(TAB / f"tcga_rai_response_dose_{STAMP}.tsv", sep="\t", index=False)
    print("\n", dose.to_string(index=False))

    # PFI within RAI-treated patients
    pfi_rows = []
    try:
        from lifelines import CoxPHFitter
        c = df[["RAI_8", "PFI", "PFI.time"]].dropna()
        c = c[c["PFI.time"] > 0]
        cph = CoxPHFitter().fit(c, duration_col="PFI.time", event_col="PFI")
        s = cph.summary.loc["RAI_8"]
        pfi_rows.append(dict(model="PFI ~ panel z (RAI-treated only)", n=len(c),
                             events=int(c["PFI"].sum()), hr=float(s["exp(coef)"]),
                             ci_low=float(s["exp(coef) lower 95%"]),
                             ci_high=float(s["exp(coef) upper 95%"]), p=float(s["p"])))
        o = df[["RAI_8", "OS", "OS.time"]].dropna()
        o = o[o["OS.time"] > 0]
        pfi_rows.append(dict(model="OS events available (not modelled)", n=len(o),
                             events=int(o["OS"].sum()), hr=np.nan, ci_low=np.nan,
                             ci_high=np.nan, p=np.nan))
    except Exception as exc:  # noqa: BLE001
        print(f"Cox skipped: {exc}")
    pfi = pd.DataFrame(pfi_rows)
    if len(pfi):
        pfi.to_csv(TAB / f"tcga_rai_response_pfi_{STAMP}.tsv", sep="\t", index=False)
        print("\n", pfi.to_string(index=False))

    # power
    from scipy.stats import nct, t as tdist

    def power(n1, n2, dd_, alpha=0.05):
        dfree = n1 + n2 - 2
        ncp = dd_ * np.sqrt(n1 * n2 / (n1 + n2))
        crit = tdist.ppf(1 - alpha / 2, dfree)
        return 1 - nct.cdf(crit, dfree, ncp) + nct.cdf(-crit, dfree, ncp)

    det = float(next(x for x in np.arange(0.05, 3, 0.005) if power(len(ncr), len(cr), x) >= 0.80))
    pw = pd.DataFrame([
        dict(quantity="n non-CR / n CR", value=f"{len(ncr)} / {len(cr)}"),
        dict(quantity="observed d", value=round(float(d), 3)),
        dict(quantity="power at observed d", value=round(float(power(len(ncr), len(cr), abs(d))), 3)),
        dict(quantity="d detectable at 80% power", value=round(det, 3)),
        dict(quantity="TCGA-THCA OS events (why survival is not the endpoint)",
             value=int(df["OS"].dropna().sum()) if df["OS"].notna().any() else "NA"),
    ])
    pw.to_csv(TAB / f"tcga_rai_response_power_{STAMP}.tsv", sep="\t", index=False)
    print("\n", pw.to_string(index=False))

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.2), facecolor="white", layout="constrained",
                             gridspec_kw={"width_ratios": [1.0, 1.15, 1.0]})
    cCR, cN = "#2f6f4f", "#9c4742"

    ax = axes[0]
    groups = [cr, ncr]
    bp = ax.boxplot(groups, tick_labels=[f"Complete\nresponse", "Non-CR\n(PR / SD / PD)"],
                    showfliers=False, patch_artist=True, widths=0.55)
    for patch, c in zip(bp["boxes"], [cCR, cN]):
        patch.set_facecolor(c); patch.set_alpha(0.28); patch.set_edgecolor(c)
    rng = np.random.default_rng(SEED)
    for i, (vals, c) in enumerate(zip(groups, [cCR, cN]), start=1):
        ax.scatter(rng.normal(i, 0.06, len(vals)), vals, s=16, color=c, alpha=0.65,
                   edgecolor="white", linewidth=0.3, zorder=3)
        ax.text(i, 0.02, f"n = {len(vals)}", ha="center", va="bottom", fontsize=9,
                color="#444", transform=ax.get_xaxis_transform(),
                bbox=dict(facecolor="white", edgecolor="none", pad=1.0, alpha=0.85))
    ax.axhline(0, color="#666", lw=0.6, ls="--")
    ax.set_ylabel("8-gene panel z (higher = differentiation preserved)")
    ax.set_title(f"a · Response to radioiodine, TCGA-THCA\nd = {d:+.2f} [{lo:+.2f}, {hi:+.2f}] · P = {p:.3g}",
                 fontsize=10.5, loc="left", fontweight="bold", color="#2a2a2a")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)

    ax = axes[1]
    if len(adj):
        a2 = adj[adj.term != "Intercept"].copy()
        y = np.arange(len(a2))
        ax.errorbar(a2["or"], y, xerr=[a2["or"] - a2.ci_low, a2.ci_high - a2["or"]],
                    fmt="o", color="#37618e", ecolor="#8fa8c4", capsize=3, ms=6)
        ax.scatter(a2.loc[a2.term == "RAI_8", "or"], y[a2.term.values == "RAI_8"],
                   s=80, facecolor="#b4472f", edgecolor="white", zorder=4)
        ax.axvline(1, color="#666", lw=0.7, ls="--")
        ax.set_xscale("log")
        ax.set_yticks(y)
        ax.set_yticklabels([{"RAI_8": "8-gene panel z", "stage_high": "Stage III/IV",
                             "braf": "BRAF V600E", "AGE": "Age",
                             "purity": "Tumour purity",
                             "rai_dose_mci": "Cumulative RAI (mCi)"}.get(t, t)
                            for t in a2.term], fontsize=9)
        ax.set_xlabel("Odds ratio for failing to achieve CR (95% CI)")
        ax.set_title(f"b · Adjusted logistic model (n = {int(a2.n_used.iloc[0])})",
                     fontsize=10.5, loc="left", fontweight="bold", color="#2a2a2a")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)

    ax = axes[2]
    order = ["Complete Response", "Partial Response", "Stable Disease",
             "Radiographic Progressive Disease"]
    present = [o for o in order if (ev.best_response == o).any()]
    vals = [ev.loc[ev.best_response == o, "RAI_8"].dropna().values for o in present]
    bp2 = ax.boxplot(vals, tick_labels=[o.replace("Radiographic ", "").replace(" ", "\n")
                                        for o in present],
                     showfliers=False, patch_artist=True, widths=0.6)
    shades = ["#2f6f4f", "#7d9a55", "#cd8b3a", "#9c4742"]
    for patch, c in zip(bp2["boxes"], shades):
        patch.set_facecolor(c); patch.set_alpha(0.30); patch.set_edgecolor(c)
    for i, v in enumerate(vals, start=1):
        ax.text(i, 0.02, f"n={len(v)}", ha="center", va="bottom", fontsize=8.5,
                color="#444", transform=ax.get_xaxis_transform(),
                bbox=dict(facecolor="white", edgecolor="none", pad=1.0, alpha=0.85))
    ax.axhline(0, color="#666", lw=0.6, ls="--")
    ax.set_ylabel("8-gene panel z")
    rk = ev[["RAI_8", "best_rank"]].dropna()
    rho2, p2 = stats.spearmanr(rk["RAI_8"], rk["best_rank"])
    ax.set_title(f"c · Ordinal response gradient\nSpearman rho = {rho2:+.3f} · P = {p2:.3g}",
                 fontsize=10.5, loc="left", fontweight="bold", color="#2a2a2a")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)

    fig.suptitle("TCGA-THCA — 8-gene differentiation panel vs documented radioiodine treatment "
                 "response (GDC BCR Biotab, mCi-dosed courses)",
                 fontsize=12.5, fontweight="bold")
    fig.savefig(FIG / f"figure_tcga_rai_best_response_{STAMP}.png", dpi=180,
                bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"figure_tcga_rai_best_response_{STAMP}.pdf",
                bbox_inches="tight", facecolor="white")
    print(f"\nwrote {FIG / f'figure_tcga_rai_best_response_{STAMP}.png'}")
    print(f"ordinal Spearman rho = {rho2:+.4f}, P = {p2:.4g}, n = {len(rk)}")


if __name__ == "__main__":
    main()
