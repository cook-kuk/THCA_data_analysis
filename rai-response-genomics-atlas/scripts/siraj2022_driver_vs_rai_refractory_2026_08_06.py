#!/usr/bin/env python3
"""Siraj 2022 (KFSHRC, n=158 PTC) — driver class vs radioiodine refractoriness.

Source: Siraj AK et al., "APOBEC SBS13 Mutational Signature - A Novel Predictor of
Radioactive Iodine Refractory Papillary Thyroid Carcinoma", Cancers 2022;14(6):1584,
doi 10.3390/cancers14061584. Supplementary tables are open (CC BY) and carry per-patient
rows, which is why this cohort is usable today while the matched WES (EGAS00001001788)
still needs an ICGC DACO application.

  Table S1  158 patients: RAI Classification (Refractory 66 / Avid 92), cumulative RAI
            activity, thyroglobulin after surgery and 6 months after RAI, PFS + censor.
  Table S3  4,788 somatic mutations tagged with Sample and RAI class.

What this can and cannot test
-----------------------------
The cohort has no transcriptome, so the 8-gene differentiation panel cannot be scored.
What it CAN test is the genomic premise of the manuscript: that the BRAF/RAS-negative
compartment is where radioiodine decisions are mechanistically untethered. If that framing
is right, driver-negative tumours should not simply be the low-risk remainder — they should
carry their own refractoriness burden.

Outputs:
  results/tables/siraj2022_driver_vs_rai_2026_08_06.tsv
  results/tables/siraj2022_patient_level_2026_08_06.tsv
  results/figures/figure_siraj2022_driver_rai_2026_08_06.{png,pdf}
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
XLSX = Path("/tmp/claude-1000/-home-seungho-personal-THCA-data-analysis/"
            "b9611370-1293-405b-bb35-984bcf91262c/scratchpad/rai_supp/siraj_extracted/"
            "SupplementaryTables.xlsx")
LOCAL = Path("/data/rai_atlas/external/siraj2022")
FIG = ROOT / "results" / "figures"
TAB = ROOT / "results" / "tables"
for d in (FIG, TAB, LOCAL):
    d.mkdir(parents=True, exist_ok=True)

STAMP = "2026_08_06"
SEED = 20260806
RAS_GENES = {"NRAS", "HRAS", "KRAS"}
FUSION_GENES = {"RET", "NTRK1", "NTRK3", "ALK", "BRAF"}   # BRAF here only as a fusion partner


def load():
    src = XLSX if XLSX.exists() else LOCAL / "SupplementaryTables.xlsx"
    if not src.exists():
        raise SystemExit(f"supplementary workbook not found at {src}")
    # keep a permanent copy next to the other external data
    if src != LOCAL / "SupplementaryTables.xlsx":
        (LOCAL / "SupplementaryTables.xlsx").write_bytes(src.read_bytes())

    s1 = pd.read_excel(src, "Supplementary Table S1", skiprows=2)
    s1.columns = [str(c).strip() for c in s1.columns]

    raw = pd.read_excel(src, "Supplementary Table S3", header=None)
    hdr = next(i for i in range(60)
               if any(str(x).strip().lower() == "sample" for x in raw.iloc[i].tolist()))
    s3 = pd.read_excel(src, "Supplementary Table S3", skiprows=hdr)
    s3.columns = [str(c).strip() for c in s3.columns]
    return s1, s3


def driver_class(s3):
    """Assign one driver class per patient from the somatic mutation table."""
    m = s3.copy()
    m["gene"] = m["Gene"].astype(str).str.upper()
    m["aa"] = m["AAChange.ensGene"].astype(str).str.upper()
    m["exonic"] = m["ExonicFunc.ensGene"].astype(str).str.lower()
    coding = m[~m["exonic"].isin(["nan", "synonymous snv"])]

    braf = set(coding.loc[coding.gene.eq("BRAF") & coding.aa.str.contains("V600E"), "Sample"])
    ras = set(coding.loc[coding.gene.isin(RAS_GENES)
                         & coding.aa.str.contains(r"G12|G13|Q61", regex=True), "Sample"])
    out = {}
    for s in m["Sample"].dropna().unique():
        if s in braf:
            out[s] = "BRAF V600E"
        elif s in ras:
            out[s] = "RAS hotspot"
        else:
            out[s] = "BRAF/RAS-negative"
    return pd.Series(out, name="driver")


def main():
    s1, s3 = load()
    s1 = s1.rename(columns={"SAMPLE": "Sample", "RAI Classification": "rai"})
    s1["rai"] = s1["rai"].astype(str).str.strip()
    s1 = s1[s1["rai"].isin(["Refractory", "Avid"])].copy()
    print(f"S1 patients: {len(s1)}  {s1['rai'].value_counts().to_dict()}")

    s1["driver"] = s1["Sample"].map(driver_class(s3))
    n_seq = s1["driver"].notna().sum()
    print(f"patients with somatic mutation data: {n_seq}/{len(s1)}")
    print("driver classes:", s1["driver"].value_counts(dropna=False).to_dict())

    dose_col = next(c for c in s1.columns if "CUMULATIVE RAI DOSE" in c.upper())
    pfs_col = next(c for c in s1.columns if "Progression-free" in c)
    cen_col = next(c for c in s1.columns if "Censor" in c)
    tg6_col = next((c for c in s1.columns if "6 months RAI" in c and "DATE" not in c.upper()), None)

    s1["dose"] = pd.to_numeric(s1[dose_col], errors="coerce")
    s1["pfs"] = pd.to_numeric(s1[pfs_col], errors="coerce")
    # `PFS Censor` is 1 = censored, 0 = progressed (verified: Avid 91 censored / 1 event,
    # Refractory 13 censored / 53 events). The event indicator is therefore its complement.
    s1["event"] = 1 - pd.to_numeric(s1[cen_col], errors="coerce")
    s1["tg6"] = pd.to_numeric(s1[tg6_col], errors="coerce") if tg6_col else np.nan
    s1["refractory"] = (s1["rai"] == "Refractory").astype(int)

    rows = []

    # --- driver class x refractoriness ---
    sub = s1.dropna(subset=["driver"])
    ct = pd.crosstab(sub["driver"], sub["rai"])
    print("\ndriver x RAI:\n", ct.to_string())
    chi = stats.chi2_contingency(ct.values)
    rows.append(dict(test="driver class x RAI refractoriness (chi-square)",
                     detail=str(ct.to_dict()), n=int(ct.values.sum()),
                     statistic=float(chi.statistic), p=float(chi.pvalue)))

    for d in ct.index:
        a = int(ct.loc[d, "Refractory"]) if "Refractory" in ct.columns else 0
        b = int(ct.loc[d, "Avid"]) if "Avid" in ct.columns else 0
        oth_r = int(ct["Refractory"].sum()) - a
        oth_a = int(ct["Avid"].sum()) - b
        orr, p = stats.fisher_exact([[a, b], [oth_r, oth_a]])
        rows.append(dict(test=f"{d} vs rest — refractory",
                         detail=f"{a}/{a+b} refractory in group; {oth_r}/{oth_r+oth_a} in rest",
                         n=int(ct.values.sum()), statistic=float(orr), p=float(p)))
        print(f"  {d:20s} refractory {a}/{a+b} ({100*a/max(a+b,1):.0f}%)  "
              f"OR={orr:.2f} P={p:.4f}")

    # --- cumulative RAI activity, thyroglobulin ---
    for col, name in [("dose", "cumulative RAI activity (mCi)"),
                      ("tg6", "thyroglobulin 6 months after RAI")]:
        a = s1.loc[s1.refractory == 1, col].dropna()
        b = s1.loc[s1.refractory == 0, col].dropna()
        if len(a) > 2 and len(b) > 2:
            u = stats.mannwhitneyu(a, b, alternative="two-sided")
            rows.append(dict(test=f"{name}: refractory vs avid",
                             detail=f"median {a.median():.1f} vs {b.median():.1f}",
                             n=len(a) + len(b), statistic=float(u.statistic), p=float(u.pvalue)))
            print(f"\n{name}: refractory median {a.median():.1f} (n={len(a)}) vs "
                  f"avid {b.median():.1f} (n={len(b)}), P={u.pvalue:.3g}")

    # --- PFS ---
    try:
        from lifelines import CoxPHFitter
        c = s1[["pfs", "event", "refractory"]].dropna()
        c = c[c["pfs"] > 0]
        cph = CoxPHFitter().fit(c, duration_col="pfs", event_col="event")
        r = cph.summary.loc["refractory"]
        rows.append(dict(test="PFS ~ RAI refractoriness (Cox) — PARTLY CIRCULAR",
                         detail=f"HR {r['exp(coef)']:.2f} "
                                f"({r['exp(coef) lower 95%']:.2f}-{r['exp(coef) upper 95%']:.2f}); "
                                "the refractoriness definition itself includes structural "
                                "progression, so this is not an independent outcome test",
                         n=len(c), statistic=float(r["exp(coef)"]), p=float(r["p"])))
        print(f"\nPFS Cox: refractory HR = {r['exp(coef)']:.2f} "
              f"({r['exp(coef) lower 95%']:.2f}-{r['exp(coef) upper 95%']:.2f}), "
              f"P = {r['p']:.3g}, events = {int(c['event'].sum())}/{len(c)}")

        cd = s1.dropna(subset=["driver"])[["pfs", "event", "driver"]].dropna()
        cd = cd[cd["pfs"] > 0]
        dm = pd.get_dummies(cd["driver"], drop_first=False).astype(int)
        if "BRAF/RAS-negative" in dm.columns:
            cd2 = pd.concat([cd[["pfs", "event"]], dm[["BRAF/RAS-negative"]]], axis=1)
            cd2.columns = ["pfs", "event", "driver_negative"]
            cph2 = CoxPHFitter().fit(cd2, duration_col="pfs", event_col="event")
            r2 = cph2.summary.loc["driver_negative"]
            rows.append(dict(test="PFS ~ BRAF/RAS-negative (Cox)",
                             detail=f"HR {r2['exp(coef)']:.2f} "
                                    f"({r2['exp(coef) lower 95%']:.2f}-"
                                    f"{r2['exp(coef) upper 95%']:.2f})",
                             n=len(cd2), statistic=float(r2["exp(coef)"]), p=float(r2["p"])))
            print(f"PFS Cox: BRAF/RAS-negative HR = {r2['exp(coef)']:.2f}, P = {r2['p']:.3g}")
    except Exception as exc:  # noqa: BLE001
        print(f"Cox skipped: {exc}")

    res = pd.DataFrame(rows)
    res.to_csv(TAB / f"siraj2022_driver_vs_rai_{STAMP}.tsv", sep="\t", index=False)
    s1[["Sample", "rai", "driver", "dose", "tg6", "pfs", "event"]].to_csv(
        TAB / f"siraj2022_patient_level_{STAMP}.tsv", sep="\t", index=False)
    print("\n", res.to_string(index=False))

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0), facecolor="white", layout="constrained")
    cR, cA = "#9c4742", "#2f6f4f"

    ax = axes[0]
    order = [d for d in ["BRAF V600E", "RAS hotspot", "BRAF/RAS-negative"] if d in ct.index]
    frac = [ct.loc[d, "Refractory"] / ct.loc[d].sum() * 100 for d in order]
    ns = [int(ct.loc[d].sum()) for d in order]
    ax.bar(np.arange(len(order)), frac, color=cR, alpha=0.85, edgecolor="white")
    overall = ct["Refractory"].sum() / ct.values.sum() * 100
    ax.axhline(overall, color="#444", ls="--", lw=0.9)
    ax.text(len(order) - 0.4, overall + 1.5, f"cohort {overall:.0f}%", fontsize=8.5, color="#444",
            ha="right")
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels([o.replace(" ", "\n") for o in order], fontsize=9)
    for i, (f_, n_) in enumerate(zip(frac, ns)):
        ax.text(i, f_ + 1.2, f"{f_:.0f}%\nn={n_}", ha="center", fontsize=8.5, color="#333")
    ax.set_ylabel("% radioiodine-refractory")
    ax.set_title("a · Refractoriness by driver class\nSiraj 2022, KFSHRC (n = %d)" % int(ct.values.sum()),
                 fontsize=10.5, loc="left", fontweight="bold")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)

    ax = axes[1]
    a = s1.loc[s1.refractory == 1, "dose"].dropna()
    b = s1.loc[s1.refractory == 0, "dose"].dropna()
    bp = ax.boxplot([b, a], tick_labels=["Avid", "Refractory"], showfliers=False,
                    patch_artist=True, widths=0.55)
    for patch, c in zip(bp["boxes"], [cA, cR]):
        patch.set_facecolor(c); patch.set_alpha(0.28); patch.set_edgecolor(c)
    rng = np.random.default_rng(SEED)
    for i, (v, c) in enumerate(zip([b, a], [cA, cR]), start=1):
        ax.scatter(rng.normal(i, 0.06, len(v)), v, s=18, color=c, alpha=0.7,
                   edgecolor="white", linewidth=0.3, zorder=3)
        ax.text(i, 0.02, f"n = {len(v)}", ha="center", va="bottom", fontsize=9, color="#444",
                transform=ax.get_xaxis_transform(),
                bbox=dict(facecolor="white", edgecolor="none", pad=1.0, alpha=0.85))
    ax.set_ylabel("Cumulative RAI activity (mCi)")
    ax.set_title("b · Treatment intensity", fontsize=10.5, loc="left", fontweight="bold")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)

    ax = axes[2]
    try:
        from lifelines import KaplanMeierFitter
        km = KaplanMeierFitter()
        for lab, c in [("Avid", cA), ("Refractory", cR)]:
            m = s1[(s1["rai"] == lab)][["pfs", "event"]].dropna()
            m = m[m["pfs"] > 0]
            km.fit(m["pfs"], m["event"], label=f"{lab} (n={len(m)})")
            km.plot_survival_function(ax=ax, color=c, ci_show=True, linewidth=1.8)
        ax.set_xlabel("Months since surgery")
        ax.set_ylabel("Progression-free")
        ax.set_ylim(0, 1.02)
        ax.legend(fontsize=8.5, frameon=False)
        ax.set_title("c · Progression-free survival\nNOT an independent test — structural "
                     "progression is part of\nthe refractoriness definition",
                     fontsize=10.5, loc="left", fontweight="bold")
    except Exception:
        ax.axis("off")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)

    fig.suptitle("Siraj 2022 KFSHRC cohort — driver class, treatment intensity and outcome "
                 "by radioiodine refractoriness", fontsize=12.5, fontweight="bold")
    fig.savefig(FIG / f"figure_siraj2022_driver_rai_{STAMP}.png", dpi=180,
                bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"figure_siraj2022_driver_rai_{STAMP}.pdf",
                bbox_inches="tight", facecolor="white")
    print(f"\nwrote {FIG / f'figure_siraj2022_driver_rai_{STAMP}.png'}")


if __name__ == "__main__":
    main()
