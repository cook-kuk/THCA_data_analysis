#!/usr/bin/env python3
"""TCGA-THCA — 8-gene panel vs persistent structural disease after radioiodine.

Companion to `tcga_rai_best_response_2026_08_06.py`. That script used the per-course
`treatment_best_response` field (RECIST-flavoured) and returned a clean null.

This script tests the endpoint that the ATA framework actually cares about and that
reviewers of a radioiodine paper will demand: **structural disease persisting after
treatment**, restricted to patients who received radioiodine. TCGA carries two proxies:

  TUMOR_STATUS                            WITH TUMOR vs TUMOR FREE at last follow-up
  CLINICAL_STATUS_WITHIN_3_MTHS_SURGERY   persistent locoregional / distant vs no evidence

Neither is an adjudicated ATA response-to-therapy category, and that is stated as a
limitation rather than papered over — no public dataset carries true ATA categories.

Outputs:
  results/tables/tcga_rai_structural_{main,adjusted}_2026_08_06.tsv
  results/figures/figure_tcga_rai_structural_2026_08_06.{png,pdf}
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "results" / "figures"
TAB = ROOT / "results" / "tables"
STAMP = "2026_08_06"
SEED = 20260806

_spec = importlib.util.spec_from_file_location(
    "rai_base", Path(__file__).resolve().parent / f"tcga_rai_best_response_{STAMP}.py")
_base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_base)

CBIO = _base.CBIO


def main():
    df = _base.build()

    cb = pd.read_csv(CBIO, sep="\t", low_memory=False)
    cols = ["patientId", "TUMOR_STATUS", "CLINICAL_STATUS_WITHIN_3_MTHS_SURGERY",
            "NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT"]
    cb = cb[[c for c in cols if c in cb.columns]].rename(columns={"patientId": "patient"})
    df = df.merge(cb, on="patient", how="left")

    rows = []

    def test(name, mask_pos, mask_neg, label_pos, label_neg):
        pos = df.loc[mask_pos, "RAI_8"].dropna().values
        neg = df.loc[mask_neg, "RAI_8"].dropna().values
        if len(pos) < 2 or len(neg) < 2:
            rows.append(dict(endpoint=name, group_pos=label_pos, group_neg=label_neg,
                             n_pos=len(pos), n_neg=len(neg), cohens_d=np.nan,
                             d_ci_low=np.nan, d_ci_high=np.nan, auc=np.nan, p=np.nan,
                             note="insufficient n"))
            return None, None
        lo, hi = _base.boot_ci_d(pos, neg)
        rows.append(dict(endpoint=name, group_pos=label_pos, group_neg=label_neg,
                         n_pos=len(pos), n_neg=len(neg),
                         cohens_d=float(_base.cohens_d(pos, neg)),
                         d_ci_low=float(lo), d_ci_high=float(hi),
                         auc=float(_base.auc_mw(pos, neg)),
                         p=float(stats.mannwhitneyu(pos, neg, alternative="two-sided").pvalue),
                         note="RAI-treated patients only"))
        return pos, neg

    # headline endpoint first: the only proxy that is unambiguously post-treatment
    nte = df["NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT"].astype(str).str.upper()
    event_v, noevent_v = test("new tumour event after initial treatment",
                              nte.eq("YES"), nte.eq("NO"), "YES", "NO")

    cs = df["CLINICAL_STATUS_WITHIN_3_MTHS_SURGERY"].astype(str).str.lower()
    test("persistent disease within 3 months of surgery",
         cs.str.contains("persistent"), cs.str.contains("no imaging evidence|no evidence"),
         "persistent locoregional/distant", "no imaging evidence")

    ts = df["TUMOR_STATUS"].astype(str).str.upper()
    test("structural disease at last follow-up",
         ts.eq("WITH TUMOR"), ts.eq("TUMOR FREE"), "WITH TUMOR", "TUMOR FREE")

    res = pd.DataFrame(rows)
    res.to_csv(TAB / f"tcga_rai_structural_main_{STAMP}.tsv", sep="\t", index=False)
    print(res.to_string(index=False))

    # Benjamini-Hochberg across the three endpoints tested here
    ok = res["p"].notna()
    ranked = res.loc[ok, "p"].rank(method="first")
    res.loc[ok, "q_bh"] = (res.loc[ok, "p"] * ok.sum() / ranked).clip(upper=1.0)
    res.to_csv(TAB / f"tcga_rai_structural_main_{STAMP}.tsv", sep="\t", index=False)

    # Adjusted models for EVERY endpoint. The GSE151179 analysis showed a nominal
    # association that vanished once thyrocyte content was held constant, so an
    # unadjusted effect here means nothing until the same covariates are applied.
    adj_all = []
    try:
        import statsmodels.formula.api as smf
        m = df.copy()
        m["stage_high"] = m["AJCC_PATHOLOGIC_TUMOR_STAGE"].astype(str).str.contains(
            "III|IV", regex=True).astype(int)
        m["braf"] = m["has_braf_v600e"].astype(str).str.lower().eq("true").astype(int)
        m["y_structural"] = np.where(ts.eq("WITH TUMOR"), 1,
                                     np.where(ts.eq("TUMOR FREE"), 0, np.nan))
        m["y_persistent"] = np.where(cs.str.contains("persistent"), 1,
                                     np.where(cs.str.contains("no imaging evidence|no evidence"),
                                              0, np.nan))
        m["y_newevent"] = np.where(nte.eq("YES"), 1, np.where(nte.eq("NO"), 0, np.nan))

        for y, nice in [("y_structural", "structural disease at last follow-up"),
                        ("y_persistent", "persistent disease within 3 months of surgery"),
                        ("y_newevent", "new tumour event after initial treatment")]:
            use = m[[y, "RAI_8", "stage_high", "braf", "AGE", "purity", "rai_dose_mci"]].apply(
                pd.to_numeric, errors="coerce").dropna()
            n_ev = int(use[y].sum())
            # keep events-per-variable defensible: drop covariates when events are scarce
            terms = ["RAI_8", "stage_high", "braf", "AGE", "purity", "rai_dose_mci"]
            while len(terms) > 1 and n_ev / len(terms) < 5:
                terms.pop()
            try:
                fit = smf.logit(f"{y} ~ " + " + ".join(terms), data=use).fit(disp=0)
            except Exception as exc:  # noqa: BLE001
                print(f"  {nice}: model failed ({exc})")
                continue
            a = pd.DataFrame({"endpoint": nice, "term": fit.params.index,
                              "or": np.exp(fit.params.values),
                              "ci_low": np.exp(fit.conf_int()[0].values),
                              "ci_high": np.exp(fit.conf_int()[1].values),
                              "p": fit.pvalues.values, "n_used": len(use), "n_events": n_ev,
                              "terms_kept": ",".join(terms)})
            adj_all.append(a)
            print(f"\nAdjusted logistic — {nice} (n={len(use)}, events={n_ev}, "
                  f"covariates kept for EPV>=5: {terms})")
            print(a[["term", "or", "ci_low", "ci_high", "p"]].to_string(index=False))
        if adj_all:
            adj = pd.concat(adj_all, ignore_index=True)
            adj.to_csv(TAB / f"tcga_rai_structural_adjusted_{STAMP}.tsv", sep="\t", index=False)
    except Exception as exc:  # noqa: BLE001
        print(f"adjusted models skipped: {exc}")

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), facecolor="white", layout="constrained",
                             gridspec_kw={"width_ratios": [1.0, 1.2]})
    cFree, cWith = "#2f6f4f", "#9c4742"

    ax = axes[0]
    if event_v is not None:
        groups = [noevent_v, event_v]
        bp = ax.boxplot(groups, tick_labels=["No new tumour event", "New tumour event"],
                        showfliers=False, patch_artist=True, widths=0.55)
        for patch, c in zip(bp["boxes"], [cFree, cWith]):
            patch.set_facecolor(c); patch.set_alpha(0.28); patch.set_edgecolor(c)
        rng = np.random.default_rng(SEED)
        for i, (vals, c) in enumerate(zip(groups, [cFree, cWith]), start=1):
            ax.scatter(rng.normal(i, 0.06, len(vals)), vals, s=16, color=c, alpha=0.65,
                       edgecolor="white", linewidth=0.3, zorder=3)
            ax.text(i, 0.02, f"n = {len(vals)}", ha="center", va="bottom", fontsize=9,
                    color="#444", transform=ax.get_xaxis_transform(),
                    bbox=dict(facecolor="white", edgecolor="none", pad=1.0, alpha=0.85))
        r0 = res.iloc[0]
        ax.axhline(0, color="#666", lw=0.6, ls="--")
        ax.set_ylabel("8-gene panel z (higher = differentiation preserved)")
        ax.set_title("a · Structural disease at last follow-up, RAI-treated\n"
                     f"d = {r0.cohens_d:+.2f} [{r0.d_ci_low:+.2f}, {r0.d_ci_high:+.2f}] · P = {r0.p:.3g}",
                     fontsize=10.5, loc="left", fontweight="bold", color="#2a2a2a")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)

    ax = axes[1]
    sub = res.dropna(subset=["cohens_d"])
    y = np.arange(len(sub))
    ax.errorbar(sub.cohens_d, y, xerr=[sub.cohens_d - sub.d_ci_low, sub.d_ci_high - sub.cohens_d],
                fmt="o", color="#37618e", ecolor="#8fa8c4", capsize=3, ms=7)
    ax.axvline(0, color="#666", lw=0.7, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([e.replace(" after initial treatment", "\nafter initial treatment")
                        .replace(" within 3 months of surgery", "\nwithin 3 months of surgery")
                        .replace(" at last follow-up", "\nat last follow-up")
                        for e in sub.endpoint], fontsize=8.5)
    for yi, (a_, b_) in enumerate(zip(sub.n_pos, sub.n_neg)):
        ax.text(0.99, yi + 0.28, f"{a_} vs {b_}", transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=7.5, color="#777")
    ax.set_xlabel("Cohen's d (disease present − disease absent), 95% bootstrap CI")
    ax.set_title("b · All structural-disease proxies in TCGA-THCA",
                 fontsize=10.5, loc="left", fontweight="bold", color="#2a2a2a")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)

    fig.suptitle("TCGA-THCA — 8-gene panel vs persistent structural disease after radioiodine",
                 fontsize=12.5, fontweight="bold")
    fig.savefig(FIG / f"figure_tcga_rai_structural_{STAMP}.png", dpi=180,
                bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"figure_tcga_rai_structural_{STAMP}.pdf",
                bbox_inches="tight", facecolor="white")
    print(f"\nwrote {FIG / f'figure_tcga_rai_structural_{STAMP}.png'}")


if __name__ == "__main__":
    main()
