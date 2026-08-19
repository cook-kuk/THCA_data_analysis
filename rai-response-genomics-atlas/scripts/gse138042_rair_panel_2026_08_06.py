#!/usr/bin/env python3
"""GSE138042 — 8-gene differentiation panel in radioiodine-refractory thyroid cancer.

Source: Sisdelli L / Cerutti JM group, *Heliyon* 2021 (PMID 33748479), GEO GSE138042.
The GEO series carries a processed mRNA count matrix (36,596 genes x 95 libraries) in which
13 columns are named `RAIR-*` — radioiodine-refractory tumours — and the remaining 82 are
the unselected surgical thyroid cohort. All eight panel genes are present, which makes this
the only public dataset found in this sweep that has BOTH a full transcriptome AND a
radioiodine-refractory label, and needs no data-access application.

Caveats built into the analysis rather than glossed:
  * the comparator group is an unselected surgical cohort (papillary, follicular, follicular
    adenoma, medullary, poorly differentiated), not a matched radioiodine-avid group, so this
    is a refractory-vs-unselected contrast;
  * the RAI label lives in the sample names and the paper, not in structured GEO metadata;
  * there is no tumour-purity estimate, so immune and stromal transcripts are used as
    composition proxies — the covariate that decided the GSE151179 and TCGA analyses.

Outputs:
  results/tables/gse138042_rair_panel_{main,per_gene}_2026_08_06.tsv
  results/figures/figure_gse138042_rair_panel_2026_08_06.{png,pdf}
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
MATRIX = Path("/data/rai_atlas/raw/GSE138042/GSE138042_mRNA_seq_thyroid.csv.gz")
# Per-patient clinical annotation, Supplementary File 1 of the source paper. Retrieved via
# the Europe PMC supplementaryFiles endpoint, which bypasses the PMC JavaScript gate:
#   https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7970325/supplementaryFiles
CLIN = Path("/data/rai_atlas/raw/GSE138042/supp/mmc1.xlsx")
FIG = ROOT / "results" / "figures"
TAB = ROOT / "results" / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

PANEL_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
IMMUNE = ["PTPRC", "CD3E", "CD68", "CD19"]
STROMA = ["COL1A1", "COL1A2", "ACTA2", "FN1"]
STAMP = "2026_08_06"
SEED = 20260806


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
    u = stats.mannwhitneyu(pos, neg, alternative="two-sided").statistic
    return u / (len(pos) * len(neg))


def score(lz, genes):
    g = [x for x in genes if x in lz.index]
    if not g:
        return None
    z = lz.loc[g]
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1, ddof=1) + 1e-9, axis=0)
    return z.mean(axis=0)


def main():
    counts = pd.read_csv(MATRIX, index_col=0)
    counts = counts.loc[counts.sum(axis=1) > 0]
    cpm = counts / counts.sum(axis=0) * 1e6
    lz = np.log2(cpm + 1)
    print(f"matrix {counts.shape[0]} genes x {counts.shape[1]} libraries")

    rair = [c for c in lz.columns if c.upper().startswith("RAIR")]
    other = [c for c in lz.columns if c not in rair]
    print(f"RAI-refractory {len(rair)} / unselected surgical cohort {len(other)}")

    panel = score(lz, PANEL_8)
    immune = score(lz, IMMUNE)
    stroma = score(lz, STROMA)

    df = pd.DataFrame({"panel_z": panel, "immune": immune, "stroma": stroma})
    df["rair"] = df.index.isin(rair).astype(int)

    # attach per-patient clinical annotation
    if CLIN.exists():
        cl = pd.read_excel(CLIN)
        cl.columns = [str(c).replace("\n", " ").strip() for c in cl.columns]
        tcol = next(c for c in cl.columns if "Iodine" in c)
        cl["key"] = cl["mRNAseq file name"].astype(str).str.replace(".fastq.gz", "", regex=False)
        cl = cl.set_index("key")
        df["rai_status"] = cl[tcol].reindex(df.index)
        df["cancer_type"] = cl["Cancer type"].reindex(df.index)
        df["braf"] = cl["BRAF mutation V600E"].reindex(df.index)
        df["age"] = pd.to_numeric(cl["age"].reindex(df.index), errors="coerce")
        df["sex"] = cl["sex"].reindex(df.index)
        matched = df["rai_status"].notna().sum()
        print(f"clinical annotation joined: {df['cancer_type'].notna().sum()}/{len(df)} samples; "
              f"radioiodine status recorded for {matched}")
        print("  ", df["rai_status"].value_counts(dropna=False).to_dict())
    df.to_csv(TAB / f"gse138042_rair_per_sample_{STAMP}.tsv", sep="\t")

    a = df.loc[df.rair == 1, "panel_z"].values
    b = df.loc[df.rair == 0, "panel_z"].values
    d = cohens_d(a, b)
    lo, hi = boot_ci_d(a, b)
    p = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue

    rows = [dict(analysis="context: panel z, RAI-refractory vs whole unselected cohort",
                 n_rair=len(a), n_other=len(b), median_rair=float(np.median(a)),
                 median_other=float(np.median(b)), cohens_d=float(d),
                 d_ci_low=float(lo), d_ci_high=float(hi),
                 auc=float(auc_mw(a, b)), p=float(p))]
    print(f"\ncontext  d = {d:+.2f} [{lo:+.2f}, {hi:+.2f}]  AUC = {auc_mw(a, b):.3f}  P = {p:.3g}")

    # ---- PRIMARY: the contrast the source paper actually defines ----
    # radioresistance vs radiosensitivity, i.e. carcinomas that all received radioiodine.
    # This removes the benign follicular adenomas and the never-treated tumours from the
    # comparator, which is the confound in the contrast above.
    clean = None
    if "rai_status" in df.columns and df["rai_status"].notna().any():
        st = df["rai_status"].astype(str).str.lower()
        res = df.loc[st.str.contains("resist"), "panel_z"].dropna().values
        sen = df.loc[st.str.contains("sensitiv"), "panel_z"].dropna().values
        if len(res) >= 2 and len(sen) >= 2:
            dc = cohens_d(res, sen)
            lc, hc = boot_ci_d(res, sen)
            pc = stats.mannwhitneyu(res, sen, alternative="two-sided").pvalue
            clean = (res, sen, dc, lc, hc, pc)
            rows.insert(0, dict(
                analysis="PRIMARY: radioresistant vs radiosensitive (both radioiodine-treated)",
                n_rair=len(res), n_other=len(sen), median_rair=float(np.median(res)),
                median_other=float(np.median(sen)), cohens_d=float(dc),
                d_ci_low=float(lc), d_ci_high=float(hc),
                auc=float(auc_mw(res, sen)), p=float(pc)))
            print(f"PRIMARY  radioresistant (n={len(res)}) vs radiosensitive (n={len(sen)}): "
                  f"d = {dc:+.2f} [{lc:+.2f}, {hc:+.2f}]  AUC = {auc_mw(res, sen):.3f}  P = {pc:.3g}")

            # power for this contrast
            from scipy.stats import nct, t as tdist

            def power(n1, n2, dd_, alpha=0.05):
                dfree = n1 + n2 - 2
                ncp = dd_ * np.sqrt(n1 * n2 / (n1 + n2))
                crit = tdist.ppf(1 - alpha / 2, dfree)
                return 1 - nct.cdf(crit, dfree, ncp) + nct.cdf(-crit, dfree, ncp)

            det = float(next(x for x in np.arange(0.05, 4, 0.005)
                             if power(len(res), len(sen), x) >= 0.80))
            print(f"         power at observed d = {power(len(res), len(sen), abs(dc)):.2f}; "
                  f"d detectable at 80% power = {det:.2f}")

    # malignant-only comparator (drop benign follicular adenoma and medullary carcinoma)
    if "cancer_type" in df.columns and df["cancer_type"].notna().any():
        ct = df["cancer_type"].astype(str).str.lower()
        keep = ~ct.str.contains("adenoma|medullar")
        sub = df[keep]
        x = sub.loc[sub.rair == 1, "panel_z"].values
        y = sub.loc[sub.rair == 0, "panel_z"].values
        if len(x) >= 2 and len(y) >= 2:
            l2, h2 = boot_ci_d(x, y)
            rows.append(dict(
                analysis="sensitivity: malignant comparator only (no adenoma/medullary)",
                n_rair=len(x), n_other=len(y), median_rair=float(np.median(x)),
                median_other=float(np.median(y)), cohens_d=float(cohens_d(x, y)),
                d_ci_low=float(l2), d_ci_high=float(h2), auc=float(auc_mw(x, y)),
                p=float(stats.mannwhitneyu(x, y, alternative="two-sided").pvalue)))
            print(f"  malignant-only comparator: d = {cohens_d(x, y):+.2f}, "
                  f"P = {stats.mannwhitneyu(x, y, alternative='two-sided').pvalue:.3g} "
                  f"(n = {len(x)} vs {len(y)})")

    # composition controls: are immune/stromal content different, and does adjustment matter?
    for col, name in [("immune", "immune content"), ("stroma", "stromal content")]:
        x = df.loc[df.rair == 1, col].values
        y = df.loc[df.rair == 0, col].values
        pp = stats.mannwhitneyu(x, y, alternative="two-sided").pvalue
        rows.append(dict(analysis=f"composition check: {name}, RAI-refractory vs rest",
                         n_rair=len(x), n_other=len(y), median_rair=float(np.median(x)),
                         median_other=float(np.median(y)), cohens_d=float(cohens_d(x, y)),
                         d_ci_low=np.nan, d_ci_high=np.nan, auc=np.nan, p=float(pp)))
        print(f"  {name}: d = {cohens_d(x, y):+.2f}, P = {pp:.3g}")

    adj = pd.DataFrame()
    try:
        import statsmodels.formula.api as smf
        fit = smf.logit("rair ~ panel_z + immune + stroma", data=df).fit(disp=0)
        adj = pd.DataFrame({"term": fit.params.index, "or": np.exp(fit.params.values),
                            "ci_low": np.exp(fit.conf_int()[0].values),
                            "ci_high": np.exp(fit.conf_int()[1].values),
                            "p": fit.pvalues.values, "n": len(df)})
        adj.to_csv(TAB / f"gse138042_rair_adjusted_{STAMP}.tsv", sep="\t", index=False)
        print("\nAdjusted logistic (odds of being RAI-refractory):")
        print(adj.to_string(index=False))
        r = adj[adj.term == "panel_z"].iloc[0]
        rows.append(dict(analysis="ADJUSTED logistic: panel z (immune + stroma held constant)",
                         n_rair=len(a), n_other=len(b), median_rair=np.nan, median_other=np.nan,
                         cohens_d=float(r["or"]), d_ci_low=float(r.ci_low),
                         d_ci_high=float(r.ci_high), auc=np.nan, p=float(r["p"])))
    except Exception as exc:  # noqa: BLE001
        print(f"adjusted model skipped: {exc}")

    res = pd.DataFrame(rows)
    res.to_csv(TAB / f"gse138042_rair_panel_main_{STAMP}.tsv", sep="\t", index=False)

    # per-gene, on the primary contrast when it exists
    if clean is not None:
        st = df["rai_status"].astype(str).str.lower()
        g_pos = df.index[st.str.contains("resist")].tolist()
        g_neg = df.index[st.str.contains("sensitiv")].tolist()
        gene_contrast = "radioresistant vs radiosensitive"
    else:
        g_pos, g_neg = rair, other
        gene_contrast = "RAI-refractory vs rest"
    pg = []
    for g in [x for x in PANEL_8 if x in lz.index]:
        v = lz.loc[g]
        x, y = v[g_pos].values, v[g_neg].values
        pg.append(dict(gene=g, contrast=gene_contrast, cohens_d=cohens_d(x, y), auc=auc_mw(x, y),
                       median_rair=float(np.median(x)), median_other=float(np.median(y)),
                       p=float(stats.mannwhitneyu(x, y, alternative="two-sided").pvalue)))
    pgd = pd.DataFrame(pg).sort_values("cohens_d")
    ranked = pgd["p"].rank(method="first")
    pgd["q_bh"] = (pgd["p"] * len(pgd) / ranked).clip(upper=1.0)
    pgd.to_csv(TAB / f"gse138042_rair_per_gene_{STAMP}.tsv", sep="\t", index=False)
    print("\nPer-gene:\n", pgd.to_string(index=False))

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.1), facecolor="white", layout="constrained",
                             gridspec_kw={"width_ratios": [0.95, 1.0, 1.15]})
    cR, cO = "#9c4742", "#2f6f4f"

    ax = axes[0]
    if clean is not None:
        res_v, sen_v, dc, lc, hc, pc = clean
        left, right = sen_v, res_v
        labels = ["Radiosensitive", "Radioresistant"]
        title = (f"a · Radioiodine-treated carcinomas only\n"
                 f"d = {dc:+.2f} [{lc:+.2f}, {hc:+.2f}] · P = {pc:.3g}")
    else:
        left, right = b, a
        labels = ["Unselected\nsurgical cohort", "RAI-refractory"]
        title = f"a · GSE138042\nd = {d:+.2f} [{lo:+.2f}, {hi:+.2f}] · P = {p:.3g}"
    bp = ax.boxplot([left, right], tick_labels=labels,
                    showfliers=False, patch_artist=True, widths=0.55)
    for patch, c in zip(bp["boxes"], [cO, cR]):
        patch.set_facecolor(c); patch.set_alpha(0.28); patch.set_edgecolor(c)
    rng = np.random.default_rng(SEED)
    for i, (v, c) in enumerate(zip([b, a], [cO, cR]), start=1):
        ax.scatter(rng.normal(i, 0.07, len(v)), v, s=26, color=c, alpha=0.75,
                   edgecolor="white", linewidth=0.4, zorder=3)
        ax.text(i, 0.02, f"n = {len(v)}", ha="center", va="bottom", fontsize=9, color="#444",
                transform=ax.get_xaxis_transform(),
                bbox=dict(facecolor="white", edgecolor="none", pad=1.0, alpha=0.85))
    ax.axhline(0, color="#666", lw=0.6, ls="--")
    ax.set_ylabel("8-gene panel z (higher = differentiation preserved)")
    ax.set_title(f"a · GSE138042\nd = {d:+.2f} [{lo:+.2f}, {hi:+.2f}] · P = {p:.3g}",
                 fontsize=10.5, loc="left", fontweight="bold")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)

    ax = axes[1]
    if len(adj):
        a2 = adj[adj.term != "Intercept"]
        y = np.arange(len(a2))
        ax.errorbar(a2["or"], y, xerr=[a2["or"] - a2.ci_low, a2.ci_high - a2["or"]],
                    fmt="o", color="#37618e", ecolor="#8fa8c4", capsize=3, ms=6)
        ax.scatter(a2.loc[a2.term == "panel_z", "or"], y[a2.term.values == "panel_z"],
                   s=80, facecolor="#b4472f", edgecolor="white", zorder=4)
        ax.axvline(1, color="#666", lw=0.7, ls="--")
        ax.set_xscale("log")
        ax.set_yticks(y)
        ax.set_yticklabels({"panel_z": "8-gene panel z", "immune": "Immune content",
                            "stroma": "Stromal content"}.get(t, t) for t in a2.term)
        ax.set_xlabel("Odds ratio for RAI-refractory (95% CI)")
        ax.set_title("b · Composition-adjusted model", fontsize=10.5, loc="left",
                     fontweight="bold")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)

    ax = axes[2]
    ax.barh(np.arange(len(pgd)), pgd.cohens_d,
            color=["#9c4742" if v < 0 else "#2f6f4f" for v in pgd.cohens_d],
            alpha=0.85, edgecolor="white")
    ax.set_yticks(np.arange(len(pgd)))
    ax.set_yticklabels(pgd.gene, fontsize=9)
    ax.axvline(0, color="#666", lw=0.7, ls="--")
    ax.set_xlabel("Cohen's d (RAI-refractory − rest)")
    ax.set_title("c · Per-gene effect", fontsize=10.5, loc="left", fontweight="bold")
    for i, q_ in enumerate(pgd.q_bh):
        ax.text(0.04, i, f"q={q_:.2g}", va="center", ha="left", fontsize=7.5, color="#333",
                transform=ax.get_yaxis_transform())
    ax.set_xlim(min(pgd.cohens_d) * 1.12, max(0.25, max(pgd.cohens_d) * 1.2))
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)

    fig.suptitle("GSE138042 — 8-gene differentiation panel in radioiodine-refractory "
                 "thyroid carcinoma", fontsize=12.5, fontweight="bold")
    fig.savefig(FIG / f"figure_gse138042_rair_panel_{STAMP}.png", dpi=180,
                bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"figure_gse138042_rair_panel_{STAMP}.pdf",
                bbox_inches="tight", facecolor="white")
    print(f"\nwrote {FIG / f'figure_gse138042_rair_panel_{STAMP}.png'}")


if __name__ == "__main__":
    main()
