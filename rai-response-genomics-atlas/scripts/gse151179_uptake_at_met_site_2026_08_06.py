#!/usr/bin/env python3
"""GSE151179 — 8-gene panel vs LESION-LEVEL RAI uptake at the metastatic site.

Prior atlas work tested the patient-level `patient rai responce` field
(Refractory 46 / Avid 6), which is severely unbalanced and returned a null.
This script instead uses the GEO characteristic

    `rai uptake at the metastatic site: Yes / No`   (Yes 25 / No 27)

which is the closest available direct measurement of whether radioiodine
actually reached and was retained by disease — i.e. RAI effect, not RAI receipt.

Design
------
Biological unit      patient (33 patients, 52 samples, repeated specimens present)
Primary test         patient-level mean panel z in tumour samples, Yes vs No
Secondary            sample-level test (non-independent; reported as such)
Adjustment           CIBERSORT tumour-purity class + lesion driver class
Negative control     non-neoplastic thyroid samples (should show no association)
Sensitivity          primary tumours only; pre-RAI collection only; leave-one-gene-out

Outputs (new files only, nothing overwritten):
  results/tables/gse151179_uptake_at_met_site_2026_08_06.tsv
  results/tables/gse151179_uptake_per_gene_2026_08_06.tsv
  results/figures/figure_gse151179_uptake_at_met_site_2026_08_06.{png,pdf}
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

PANEL_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
SEED = 20260806
STAMP = "2026_08_06"


def cohens_d(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return np.nan
    sp = np.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    return (a.mean() - b.mean()) / sp if sp > 0 else np.nan


def auc_mw(pos, neg):
    """AUC via Mann-Whitney U (probability a 'pos' sample scores higher)."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if len(pos) < 1 or len(neg) < 1:
        return np.nan
    u = stats.mannwhitneyu(pos, neg, alternative="two-sided").statistic
    return u / (len(pos) * len(neg))


def boot_ci_d(a, b, n=5000, seed=SEED):
    rng = np.random.default_rng(seed)
    a, b = np.asarray(a, float), np.asarray(b, float)
    vals = [cohens_d(rng.choice(a, len(a), replace=True),
                     rng.choice(b, len(b), replace=True)) for _ in range(n)]
    vals = np.asarray([v for v in vals if np.isfinite(v)])
    return (np.percentile(vals, 2.5), np.percentile(vals, 97.5)) if len(vals) else (np.nan, np.nan)


def raw_characteristics():
    """Parse !Sample_characteristics_ch1 from the series matrix into a GSM-indexed frame.

    The processed label table carries an all-NaN `tumor_purity`; the CIBERSORT purity
    class only exists in the raw GEO characteristics, and it is the key confounder here
    (panel genes are thyrocyte-specific, so thyrocyte content drives the score).
    """
    import gzip
    path = DATA / "raw" / "GSE151179" / "GSE151179_series_matrix.txt.gz"
    fields, gsms = {}, []
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_geo_accession"):
                gsms = [v.strip().strip('"') for v in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = [v.strip().strip('"') for v in line.rstrip("\n").split("\t")[1:]]
                key = vals[0].split(":")[0].strip()
                fields[key] = [v.split(":", 1)[1].strip() if ":" in v else np.nan for v in vals]
            elif line.startswith("!series_matrix_table_begin"):
                break
    return pd.DataFrame(fields, index=gsms)


def load_panel_matrix():
    """Return (panel gene x sample z-matrix, metadata) aligned on GSM ids."""
    expr = pd.read_csv(DATA / "interim" / "GSE151179_expression.tsv", sep="\t", index_col=0)
    p2g = pd.read_csv(DATA / "raw" / "GSE151179" / "GPL23159_probe2gene.tsv", sep="\t")
    gcol = next(c for c in p2g.columns if "gene" in c.lower() or "symbol" in c.lower())
    pcol = next(c for c in p2g.columns if "probe" in c.lower() or c.lower().startswith("id"))
    p2g = p2g[[pcol, gcol]].dropna()
    p2g.columns = ["probe_id", "symbol"]
    p2g["probe_id"] = p2g["probe_id"].astype(str)
    expr.index = expr.index.astype(str)

    keep = p2g[p2g["symbol"].isin(PANEL_8)]
    rows = []
    for g, sub in keep.groupby("symbol"):
        probes = [p for p in sub["probe_id"] if p in expr.index]
        if not probes:
            continue
        # collapse multi-probe genes by the probe with highest mean signal
        block = expr.loc[probes]
        best = block.mean(axis=1).idxmax()
        rows.append(pd.Series(block.loc[best], name=g))
    gm = pd.DataFrame(rows)
    # per-gene z across all 52 arrays (single platform, RMA-normalised)
    gz = gm.sub(gm.mean(axis=1), axis=0).div(gm.std(axis=1, ddof=1), axis=0)

    meta = pd.read_csv(DATA / "processed" / "GSE151179_label_joined.tsv", sep="\t")
    idcol = next(c for c in meta.columns if meta[c].astype(str).str.startswith("GSM").any())
    meta = meta.set_index(idcol)
    common = [s for s in gz.columns if s in meta.index]
    return gz[common], meta.loc[common]


def summarise(name, unit, pos, neg, direction_note=""):
    """pos = uptake Yes group values, neg = uptake No group values."""
    if len(pos) < 2 or len(neg) < 2:
        return dict(analysis=name, unit=unit, n_uptake_yes=len(pos), n_uptake_no=len(neg),
                    median_yes=np.nan, median_no=np.nan, cohens_d=np.nan,
                    d_ci_low=np.nan, d_ci_high=np.nan, auc=np.nan, p_mannwhitney=np.nan,
                    note="insufficient n")
    p = stats.mannwhitneyu(pos, neg, alternative="two-sided").pvalue
    d = cohens_d(pos, neg)
    lo, hi = boot_ci_d(pos, neg)
    return dict(analysis=name, unit=unit, n_uptake_yes=len(pos), n_uptake_no=len(neg),
                median_yes=float(np.median(pos)), median_no=float(np.median(neg)),
                cohens_d=float(d), d_ci_low=float(lo), d_ci_high=float(hi),
                auc=float(auc_mw(pos, neg)), p_mannwhitney=float(p), note=direction_note)


def main():
    gz, meta = load_panel_matrix()
    found = list(gz.index)
    print(f"panel genes mapped: {len(found)}/8 -> {found}")

    panel_z = gz.mean(axis=0)
    df = meta.copy()
    df["panel_z"] = panel_z.reindex(df.index)
    df["uptake"] = df["rai_uptake_met"].str.lower().map({"yes": "Yes", "no": "No"})
    df["is_tumor"] = df["Sample_source_name_ch1"].str.contains("papillary", case=False)
    df["patient"] = df["patient_id"] if "patient_id" in df.columns else df.filter(
        like="patient").iloc[:, 0]

    raw = raw_characteristics()
    pcol = next(c for c in raw.columns if "purity" in c.lower())
    df["purity_class"] = raw[pcol].reindex(df.index)
    print("purity classes:", df["purity_class"].value_counts(dropna=False).to_dict())

    tum = df[df["is_tumor"]].copy()
    norm = df[~df["is_tumor"]].copy()
    print(f"tumour samples {len(tum)} (Yes {(tum.uptake=='Yes').sum()} / "
          f"No {(tum.uptake=='No').sum()}) across {tum.patient.nunique()} patients")

    rows = []

    # --- PRIMARY: patient-level (one value per patient, tumour samples only) ---
    pat = tum.groupby(["patient", "uptake"], as_index=False)["panel_z"].mean()
    dup = pat.patient.duplicated(keep=False)
    if dup.any():
        print("WARNING: patients with discordant uptake labels:", pat[dup].patient.tolist())
    rows.append(summarise("PRIMARY patient-level (tumour samples)", "patient",
                          pat.loc[pat.uptake == "Yes", "panel_z"].values,
                          pat.loc[pat.uptake == "No", "panel_z"].values,
                          "one mean panel z per patient"))

    # --- Secondary: sample-level (non-independent) ---
    rows.append(summarise("sample-level (tumour samples, non-independent)", "sample",
                          tum.loc[tum.uptake == "Yes", "panel_z"].values,
                          tum.loc[tum.uptake == "No", "panel_z"].values,
                          "repeated specimens per patient not accounted"))

    # --- Sensitivity: primary tumours only ---
    pt = tum[tum["sample_type"].str.contains("primary", case=False)]
    rows.append(summarise("sensitivity: primary tumours only", "sample",
                          pt.loc[pt.uptake == "Yes", "panel_z"].values,
                          pt.loc[pt.uptake == "No", "panel_z"].values,
                          "one specimen per patient by construction"))

    # --- Sensitivity: specimens collected BEFORE RAI ---
    pre = tum[tum["collection_timing"].str.lower() == "before"]
    rows.append(summarise("sensitivity: pre-RAI specimens only", "sample",
                          pre.loc[pre.uptake == "Yes", "panel_z"].values,
                          pre.loc[pre.uptake == "No", "panel_z"].values,
                          "removes post-treatment sampling effect"))

    # --- Sensitivity: high-purity tumours only (thyrocyte-content confounder) ---
    hp = tum[tum["purity_class"].astype(str).str.contains("high", case=False, na=False)]
    rows.append(summarise("sensitivity: high-purity tumours only", "sample",
                          hp.loc[hp.uptake == "Yes", "panel_z"].values,
                          hp.loc[hp.uptake == "No", "panel_z"].values,
                          "controls thyrocyte-content confounding"))

    # --- Negative control: non-neoplastic thyroid ---
    rows.append(summarise("NEGATIVE CONTROL: non-neoplastic thyroid", "sample",
                          norm.loc[norm.uptake == "Yes", "panel_z"].values,
                          norm.loc[norm.uptake == "No", "panel_z"].values,
                          "expect no association"))

    res = pd.DataFrame(rows)

    # --- Adjusted model: panel_z ~ uptake + purity class + driver class ---
    try:
        import statsmodels.formula.api as smf
        m = pd.DataFrame({
            "panel_z": tum["panel_z"].astype(float),
            "purity": tum["purity_class"].astype(str),
            # collapse driver levels with <3 samples so the design stays full rank
            "driver": tum["lesion_class"].astype(str).where(
                tum["lesion_class"].map(tum["lesion_class"].value_counts()) >= 3, "other"),
            "uptake_bin": (tum["uptake"] == "Yes").astype(int),
        }).dropna()
        fit = smf.ols("panel_z ~ uptake_bin + C(purity) + C(driver)", data=m).fit()
        adj = dict(analysis="ADJUSTED OLS (purity + driver), sample-level", unit="sample",
                   n_uptake_yes=int(m.uptake_bin.sum()), n_uptake_no=int((1 - m.uptake_bin).sum()),
                   median_yes=np.nan, median_no=np.nan,
                   cohens_d=float(fit.params["uptake_bin"]),
                   d_ci_low=float(fit.conf_int().loc["uptake_bin", 0]),
                   d_ci_high=float(fit.conf_int().loc["uptake_bin", 1]),
                   auc=np.nan, p_mannwhitney=float(fit.pvalues["uptake_bin"]),
                   note="coefficient is beta in panel-z units, not Cohen's d")
        res = pd.concat([res, pd.DataFrame([adj])], ignore_index=True)
        print(fit.summary().tables[1])
    except Exception as exc:  # noqa: BLE001
        print(f"adjusted model skipped: {exc}")

    # --- Power: what this cohort could and could not have detected ---
    from scipy.stats import nct, t as tdist

    def power(n1, n2, d, alpha=0.05):
        dfree = n1 + n2 - 2
        ncp = d * np.sqrt(n1 * n2 / (n1 + n2))
        crit = tdist.ppf(1 - alpha / 2, dfree)
        return 1 - nct.cdf(crit, dfree, ncp) + nct.cdf(-crit, dfree, ncp)

    n_yes, n_no = int(res.iloc[0]["n_uptake_yes"]), int(res.iloc[0]["n_uptake_no"])
    obs_d = float(res.iloc[0]["cohens_d"])
    z_a, z_b = 1.959963985, 0.8416212336
    n_needed = int(np.ceil(2 * (z_a + z_b) ** 2 / obs_d ** 2)) if obs_d else -1
    pw = pd.DataFrame([
        dict(quantity="observed patient-level d", value=round(obs_d, 3)),
        dict(quantity="power to detect observed d", value=round(power(n_yes, n_no, obs_d), 3)),
        dict(quantity="d detectable at 80% power", value=round(
            float(next(x for x in np.arange(0.05, 3, 0.005) if power(n_yes, n_no, x) >= 0.80)), 3)),
        dict(quantity="patients per group needed for observed d at 80% power", value=n_needed),
        dict(quantity="total patients needed", value=2 * n_needed),
    ])
    pw.to_csv(TAB / f"gse151179_uptake_power_{STAMP}.tsv", sep="\t", index=False)
    print("\nPower analysis:\n", pw.to_string(index=False))

    res.to_csv(TAB / f"gse151179_uptake_at_met_site_{STAMP}.tsv", sep="\t", index=False)
    print("\n", res.to_string(index=False))

    # --- Per-gene effects (patient-level) ---
    per_gene = []
    for g in found:
        s = gz.loc[g]
        t = tum.copy()
        t["v"] = s.reindex(t.index)
        pg = t.groupby(["patient", "uptake"], as_index=False)["v"].mean()
        yes = pg.loc[pg.uptake == "Yes", "v"].values
        no = pg.loc[pg.uptake == "No", "v"].values
        per_gene.append(dict(gene=g, n_yes=len(yes), n_no=len(no),
                             cohens_d=cohens_d(yes, no), auc=auc_mw(yes, no),
                             p=stats.mannwhitneyu(yes, no, alternative="two-sided").pvalue))
    pg_df = pd.DataFrame(per_gene).sort_values("cohens_d", ascending=False)
    # Benjamini-Hochberg across the 8 panel genes
    ranked = pg_df["p"].rank(method="first")
    pg_df["q_bh"] = (pg_df["p"] * len(pg_df) / ranked).clip(upper=1.0)
    pg_df.to_csv(TAB / f"gse151179_uptake_per_gene_{STAMP}.tsv", sep="\t", index=False)
    print("\nPer-gene (patient-level):\n", pg_df.to_string(index=False))

    # --- Leave-one-gene-out on the primary test ---
    logo = []
    for drop in found:
        keep = [g for g in found if g != drop]
        s = gz.loc[keep].mean(axis=0)
        t = tum.copy()
        t["v"] = s.reindex(t.index)
        pgg = t.groupby(["patient", "uptake"], as_index=False)["v"].mean()
        yes = pgg.loc[pgg.uptake == "Yes", "v"].values
        no = pgg.loc[pgg.uptake == "No", "v"].values
        logo.append(dict(dropped=drop, cohens_d=cohens_d(yes, no), auc=auc_mw(yes, no),
                         p=stats.mannwhitneyu(yes, no, alternative="two-sided").pvalue))
    logo_df = pd.DataFrame(logo)
    print("\nLeave-one-gene-out (patient-level primary test):\n", logo_df.to_string(index=False))
    logo_df.to_csv(TAB / f"gse151179_uptake_logo_{STAMP}.tsv", sep="\t", index=False)

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.4), facecolor="white", layout="constrained",
                             gridspec_kw={"width_ratios": [0.9, 1.45, 1.0]})
    col = {"Yes": "#2f6f4f", "No": "#9c4742"}

    ax = axes[0]
    groups = [pat.loc[pat.uptake == "No", "panel_z"].values,
              pat.loc[pat.uptake == "Yes", "panel_z"].values]
    bp = ax.boxplot(groups, tick_labels=["No uptake", "RAI uptake"], showfliers=False,
                    patch_artist=True, widths=0.55)
    for patch, c in zip(bp["boxes"], [col["No"], col["Yes"]]):
        patch.set_facecolor(c); patch.set_alpha(0.30); patch.set_edgecolor(c)
    rng = np.random.default_rng(SEED)
    for i, vals in enumerate(groups, start=1):
        ax.scatter(rng.normal(i, 0.055, len(vals)), vals, s=34,
                   color=col["No"] if i == 1 else col["Yes"], alpha=0.9,
                   edgecolor="white", linewidth=0.6, zorder=3)
        ax.text(i, 0.02, f"n = {len(vals)}", ha="center", va="bottom", fontsize=9,
                color="#444", transform=ax.get_xaxis_transform(),
                bbox=dict(facecolor="white", edgecolor="none", pad=1.0, alpha=0.85))
    prim = res.iloc[0]
    ax.axhline(0, color="#666", lw=0.6, ls="--")
    ax.set_ylabel("8-gene panel z (patient mean)")
    ax.set_title(f"a · Patient-level, tumour specimens\nd = {prim.cohens_d:+.2f} "
                 f"[{prim.d_ci_low:+.2f}, {prim.d_ci_high:+.2f}] · P = {prim.p_mannwhitney:.3g}",
                 fontsize=10.5, loc="left", fontweight="bold", color="#2a2a2a")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    ax = axes[1]
    label_map = {
        "sample-level (tumour samples, non-independent)": "All tumours (sample-level)",
        "sensitivity: primary tumours only": "Primary tumours only",
        "sensitivity: pre-RAI specimens only": "Pre-RAI specimens only",
        "sensitivity: high-purity tumours only": "High-purity tumours only",
        "NEGATIVE CONTROL: non-neoplastic thyroid": "NEGATIVE CONTROL\nnon-neoplastic thyroid",
        "ADJUSTED OLS (purity + driver), sample-level": "ADJUSTED for purity + driver\n(β, panel-z units)",
    }
    sub = res[res.analysis.isin(label_map)]
    y = np.arange(len(sub))
    is_adj = sub.analysis.str.startswith("ADJUSTED").values
    ax.errorbar(sub.cohens_d, y,
                xerr=[sub.cohens_d - sub.d_ci_low, sub.d_ci_high - sub.cohens_d],
                fmt="o", color="#37618e", ecolor="#8fa8c4", capsize=3, ms=6, lw=1.4)
    ax.scatter(sub.cohens_d[is_adj], y[is_adj], s=70, facecolor="#b4472f",
               edgecolor="white", zorder=4)
    ax.axvline(0, color="#666", lw=0.7, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([label_map[s] for s in sub.analysis], fontsize=8.5)
    for yi, (n1, n2) in enumerate(zip(sub.n_uptake_yes, sub.n_uptake_no)):
        ax.text(0.99, yi + 0.30, f"{n1} vs {n2}", transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=7.5, color="#777")
    ax.set_xlabel("Cohen's d (uptake Yes − No), 95% bootstrap CI")
    ax.set_title("b · Sensitivity and negative control", fontsize=10.5, loc="left",
                 fontweight="bold", color="#2a2a2a")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    ax = axes[2]
    o = pg_df.sort_values("cohens_d")
    ax.barh(np.arange(len(o)), o.cohens_d,
            color=["#2f6f4f" if v > 0 else "#9c4742" for v in o.cohens_d],
            alpha=0.85, edgecolor="white")
    ax.set_yticks(np.arange(len(o)))
    ax.set_yticklabels(o.gene, fontsize=9)
    ax.axvline(0, color="#666", lw=0.7, ls="--")
    ax.set_xlabel("Cohen's d (patient-level)")
    ax.set_title("c · Per-gene effect, uptake Yes vs No", fontsize=10.5, loc="left",
                 fontweight="bold", color="#2a2a2a")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for i, (d_, q_) in enumerate(zip(o.cohens_d, o.q_bh)):
        ax.text(d_ + (0.03 if d_ >= 0 else -0.03), i, f"q={q_:.3f}",
                va="center", ha="left" if d_ >= 0 else "right", fontsize=7.5, color="#333")

    fig.suptitle("GSE151179 — 8-gene differentiation panel vs RAI uptake at the metastatic site "
                 "(lesion-level RAI effect label)", fontsize=12.5, fontweight="bold")
    fig.savefig(FIG / f"figure_gse151179_uptake_at_met_site_{STAMP}.png", dpi=180,
                bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"figure_gse151179_uptake_at_met_site_{STAMP}.pdf",
                bbox_inches="tight", facecolor="white")
    print(f"\nwrote {FIG / f'figure_gse151179_uptake_at_met_site_{STAMP}.png'}")


if __name__ == "__main__":
    main()
