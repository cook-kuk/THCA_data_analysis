"""
H13 — causal mediation: MAPK output -> 8-gene methylation -> HT immune axis.

Three sets of variables on TCGA-THCA primary tumors:
  X = MAPK-output panel mean z-score (DUSP4/5/6, SPRY2/4, ETV4/5, PHLDA1, CCND1)
  M = 8-gene HM450 mean beta (DIO1, FOXE1, NKX2-1, PAX8, SLC5A5, TG, TPO, TSHR)
  Y = HT-13 panel mean z-score (HLA-DR/DP/DQ, CD79A/B, MS4A1, AICDA, CXCL13, CCR6, IFNG)

Strata:
  - BRAF_like + cPTC
  - RAS_like + FVPTC
  - All TCGA-THCA primary tumors with all three layers

Per stratum we run forward (X -> M -> Y) and reverse (Y -> M -> X via Y -> X -> M
ordering) mediation with 5000-bootstrap proportion-mediated CI.

Outputs:
  - h13_mediation_results.tsv
  - h13_dm_polarity_audit.tsv
  - h13_polarity_histograms.png
  - H13_REPORT.md (separately)
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/"
           "p2_braf_nature_sprint_2026_05_09/h13_mediation")
OUT.mkdir(parents=True, exist_ok=True)

MASTER = "/data/thca/repo_results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv"
EXPR = "/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_zscore.tsv"
METH = "/data/thca/repo_results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv"

MAPK_GENES = ["DUSP4", "DUSP5", "DUSP6", "SPRY2", "SPRY4",
              "ETV4", "ETV5", "PHLDA1", "CCND1"]
HT_GENES = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
            "HLA-DQA1", "HLA-DQB1", "CD79A", "CD79B", "MS4A1",
            "AICDA", "CXCL13", "CCR6", "IFNG"]
EIGHT_G = ["DIO1", "FOXE1", "NKX2-1", "PAX8",
           "SLC5A5", "TG", "TPO", "TSHR"]

RNG = np.random.default_rng(20260509)
NBOOT = 5000


# --------------------------------------------------------------------- loaders
def load_master() -> pd.DataFrame:
    df = pd.read_csv(MASTER, sep="\t", low_memory=False)
    # primary tumor only
    df = df[df["normal_vs_tumor"] == "tumor"].copy()
    # 12-char short id for methylation join
    df["meth_id"] = df["sample_short"].str.slice(0, 12)
    return df


def panel_mean_z(expr_path: str, genes: list[str]) -> pd.Series:
    """Read selected gene rows, exact-match symbol, then average per sample."""
    keep = set(genes)
    rows = []
    with open(expr_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            sym = line.split("\t", 1)[0]
            if sym in keep:
                vals = line.rstrip("\n").split("\t")
                rows.append(vals)
    df = pd.DataFrame(rows, columns=header).set_index("gene_symbol")
    df = df.apply(pd.to_numeric, errors="coerce")
    found = sorted(set(df.index) & keep)
    missing = sorted(keep - set(df.index))
    print(f"  panel: found {len(found)}/{len(keep)} (missing: {missing})")
    # gene-level z-score is already done -> just average
    return df.loc[found].mean(axis=0).rename("panel_mean")


def load_methylation() -> pd.DataFrame:
    return pd.read_csv(METH, sep="\t")


# ------------------------------------------------------------------ mediation
def ols_fit(y: np.ndarray, X: np.ndarray) -> sm.regression.linear_model.RegressionResultsWrapper:
    Xc = sm.add_constant(X, has_constant="add")
    return sm.OLS(y, Xc).fit()


def mediation_one(x: np.ndarray, m: np.ndarray, y: np.ndarray) -> dict:
    """Baron & Kenny mediation. Returns total/direct/indirect with bootstrap CI."""
    # standardise so coefficients are comparable across strata
    def z(v):
        v = np.asarray(v, dtype=float)
        return (v - np.nanmean(v)) / np.nanstd(v, ddof=1)

    x, m, y = z(x), z(m), z(y)

    # total: Y ~ X
    f_tot = ols_fit(y, x)
    c = f_tot.params[1]

    # M ~ X
    f_a = ols_fit(m, x)
    a = f_a.params[1]

    # Y ~ X + M
    f_b = ols_fit(y, np.column_stack([x, m]))
    c_prime = f_b.params[1]
    b = f_b.params[2]

    indirect = a * b
    prop_med = indirect / c if c != 0 else np.nan

    # bootstrap
    n = len(x)
    a_boot = np.empty(NBOOT)
    b_boot = np.empty(NBOOT)
    cprime_boot = np.empty(NBOOT)
    c_boot = np.empty(NBOOT)
    for i in range(NBOOT):
        idx = RNG.integers(0, n, n)
        xi, mi, yi = x[idx], m[idx], y[idx]
        try:
            a_boot[i] = ols_fit(mi, xi).params[1]
            f_y_full = ols_fit(yi, np.column_stack([xi, mi]))
            cprime_boot[i] = f_y_full.params[1]
            b_boot[i] = f_y_full.params[2]
            c_boot[i] = ols_fit(yi, xi).params[1]
        except Exception:
            a_boot[i] = b_boot[i] = cprime_boot[i] = c_boot[i] = np.nan
    ind_boot = a_boot * b_boot
    prop_boot = np.where(c_boot != 0, ind_boot / c_boot, np.nan)

    def ci(arr):
        arr = arr[~np.isnan(arr)]
        if len(arr) == 0:
            return (np.nan, np.nan)
        return tuple(np.quantile(arr, [0.025, 0.975]))

    # AIC/BIC for the "Y ~ X + M" full model
    aic = f_b.aic
    bic = f_b.bic

    return dict(
        n=n,
        a=a, b=b,
        total_c=c,
        direct_c_prime=c_prime,
        indirect_ab=indirect,
        prop_mediated=prop_med,
        a_p=f_a.pvalues[1],
        b_p=f_b.pvalues[2],
        total_p=f_tot.pvalues[1],
        direct_p=f_b.pvalues[1],
        ci_indirect_lo=ci(ind_boot)[0],
        ci_indirect_hi=ci(ind_boot)[1],
        ci_prop_lo=ci(prop_boot)[0],
        ci_prop_hi=ci(prop_boot)[1],
        aic=aic, bic=bic,
    )


# ---------------------------------------------------------------- driver code
def main():
    print("[1] master ...")
    master = load_master()
    print(f"  master tumor rows: {len(master)}")

    print("[2] expression panels ...")
    print(" MAPK")
    mapk = panel_mean_z(EXPR, MAPK_GENES)
    print(" HT")
    ht = panel_mean_z(EXPR, HT_GENES)

    print("[3] methylation ...")
    meth = load_methylation()
    print(f"  meth rows: {len(meth)}")

    # --- assemble per-sample table
    df = master[[
        "sample_id", "sample_short", "meth_id",
        "histology_subtype", "driver_anchor", "molecular_subtype",
    ]].copy()

    df["MAPK"] = df["sample_id"].map(mapk).astype(float)
    df["HT"] = df["sample_id"].map(ht).astype(float)
    df = df.merge(meth[["sample_short", "mean_8g_beta"]],
                  left_on="meth_id", right_on="sample_short",
                  how="left", suffixes=("", "_meth"))

    full = df.dropna(subset=["MAPK", "HT", "mean_8g_beta"]).copy()
    print(f"  samples with all three layers: {len(full)}")

    # ------------------------------------------------------------- strata
    strata = {
        "BRAF_like_cPTC": full[(full["molecular_subtype"] == "BRAF_like") &
                               (full["histology_subtype"] == "cPTC")],
        "RAS_like_FVPTC": full[(full["molecular_subtype"] == "RAS_like") &
                               (full["histology_subtype"] == "FVPTC")],
        "All_TCGA_THCA": full,
    }

    # ------------------------------------------------------------- mediation
    print("[4] mediation ...")
    rows = []
    for name, sub in strata.items():
        if len(sub) < 25:
            print(f"  skip {name}: n={len(sub)}")
            continue
        x = sub["MAPK"].to_numpy()
        m = sub["mean_8g_beta"].to_numpy()
        y = sub["HT"].to_numpy()

        fwd = mediation_one(x, m, y)
        fwd["stratum"] = name
        fwd["direction"] = "forward (MAPK->8gMeth->HT)"
        rows.append(fwd)

        # Reverse: HT -> 8gMeth -> MAPK  (HT is X, MAPK is Y)
        rev = mediation_one(y, m, x)
        rev["stratum"] = name
        rev["direction"] = "reverse (HT->8gMeth->MAPK)"
        rows.append(rev)

        # Alternate ordering check: MAPK -> HT -> 8gMeth
        # (i.e. is methylation downstream of immunity?)
        alt = mediation_one(x, y, m)
        alt["stratum"] = name
        alt["direction"] = "alt (MAPK->HT->8gMeth)"
        rows.append(alt)

        print(f"  {name} n={fwd['n']}  forward indirect={fwd['indirect_ab']:+.3f} "
              f"prop_med={fwd['prop_mediated']:+.2f} "
              f"CI=({fwd['ci_prop_lo']:+.2f},{fwd['ci_prop_hi']:+.2f})")

    res = pd.DataFrame(rows)
    cols_first = ["stratum", "direction", "n",
                  "total_c", "total_p",
                  "direct_c_prime", "direct_p",
                  "a", "a_p", "b", "b_p",
                  "indirect_ab", "ci_indirect_lo", "ci_indirect_hi",
                  "prop_mediated", "ci_prop_lo", "ci_prop_hi",
                  "aic", "bic"]
    res = res[cols_first]
    res.to_csv(OUT / "h13_mediation_results.tsv", sep="\t", index=False, float_format="%.4f")
    print(f"  wrote {OUT/'h13_mediation_results.tsv'}")

    # ------------------------------------------------------------- DM polarity
    print("[5] DM polarity audit ...")
    # use master-level dm column (paper-1 standard) joined to methylation
    dm_master = master[["meth_id", "dm",
                        "histology_subtype", "molecular_subtype"]].copy()
    dm_master = dm_master.merge(meth[["sample_short", "mean_8g_beta"] + EIGHT_G],
                                left_on="meth_id", right_on="sample_short",
                                how="inner")
    dm_master["dm"] = dm_master["dm"].astype(str)
    dm_master = dm_master[dm_master["dm"].isin(["dm1", "dm2", "DM1", "DM2", "1", "2"])]
    dm_master["dm"] = dm_master["dm"].str.lower().replace(
        {"1": "dm1", "2": "dm2"})

    audit_rows = []
    for name, sel in [
        ("BRAF_like_cPTC", (dm_master["molecular_subtype"] == "BRAF_like") &
                            (dm_master["histology_subtype"] == "cPTC")),
        ("RAS_like_FVPTC", (dm_master["molecular_subtype"] == "RAS_like") &
                            (dm_master["histology_subtype"] == "FVPTC")),
        ("All_TCGA_THCA", pd.Series([True] * len(dm_master),
                                     index=dm_master.index)),
    ]:
        sub = dm_master[sel]
        for dm_label in ["dm1", "dm2"]:
            grp = sub[sub["dm"] == dm_label]
            audit_rows.append(dict(
                stratum=name,
                dm=dm_label,
                n=len(grp),
                mean_8g_beta=grp["mean_8g_beta"].mean(),
                std_8g_beta=grp["mean_8g_beta"].std(),
                **{f"{g}_mean": grp[g].mean() for g in EIGHT_G}
            ))
        # contrast
        a = sub[sub["dm"] == "dm1"]["mean_8g_beta"].dropna()
        b = sub[sub["dm"] == "dm2"]["mean_8g_beta"].dropna()
        if len(a) > 1 and len(b) > 1:
            d = (a.mean() - b.mean()) / np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
            audit_rows.append(dict(
                stratum=name, dm="dm1_minus_dm2_cohend",
                n=len(a) + len(b),
                mean_8g_beta=a.mean() - b.mean(),
                std_8g_beta=d,
            ))

    aud = pd.DataFrame(audit_rows)
    aud.to_csv(OUT / "h13_dm_polarity_audit.tsv", sep="\t", index=False,
               float_format="%.4f")
    print(f"  wrote {OUT/'h13_dm_polarity_audit.tsv'}")

    # quick ascii print
    print("  -- per-stratum mean 8g-beta --")
    for name in ["BRAF_like_cPTC", "RAS_like_FVPTC", "All_TCGA_THCA"]:
        for dm in ["dm1", "dm2"]:
            r = aud[(aud["stratum"] == name) & (aud["dm"] == dm)]
            if not r.empty:
                print(f"    {name:18s} {dm}: n={int(r['n'].iat[0]):3d}  "
                      f"beta={r['mean_8g_beta'].iat[0]:.3f}")

    # ------------------------------------------------------------- histograms
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=False)
    for ax, name in zip(axes, ["BRAF_like_cPTC", "RAS_like_FVPTC",
                               "All_TCGA_THCA"]):
        sel = {
            "BRAF_like_cPTC": (dm_master["molecular_subtype"] == "BRAF_like") &
                              (dm_master["histology_subtype"] == "cPTC"),
            "RAS_like_FVPTC": (dm_master["molecular_subtype"] == "RAS_like") &
                              (dm_master["histology_subtype"] == "FVPTC"),
            "All_TCGA_THCA":  pd.Series([True] * len(dm_master),
                                        index=dm_master.index),
        }[name]
        sub = dm_master[sel]
        for dm_label, color in [("dm1", "tab:red"), ("dm2", "tab:blue")]:
            v = sub[sub["dm"] == dm_label]["mean_8g_beta"].dropna()
            ax.hist(v, bins=20, alpha=0.55, color=color,
                    label=f"{dm_label} (n={len(v)}, mean={v.mean():.2f})")
        ax.set_title(name)
        ax.set_xlabel("mean 8-gene HM450 beta")
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "h13_polarity_histograms.png", dpi=150)
    plt.close(fig)
    print(f"  wrote {OUT/'h13_polarity_histograms.png'}")

    # save short JSON summary for the report
    summary = {
        "n_full": int(len(full)),
        "strata_counts": {k: int(len(v)) for k, v in strata.items()},
        "key": [
            res.set_index(["stratum", "direction"]).loc[
                (s, d), ["n", "indirect_ab", "prop_mediated",
                        "ci_prop_lo", "ci_prop_hi", "aic"]].to_dict()
            for s in ["BRAF_like_cPTC", "RAS_like_FVPTC", "All_TCGA_THCA"]
            for d in ["forward (MAPK->8gMeth->HT)",
                      "reverse (HT->8gMeth->MAPK)",
                      "alt (MAPK->HT->8gMeth)"]
        ],
    }
    (OUT / "h13_summary.json").write_text(json.dumps(summary, indent=2,
                                                      default=str))
    print("[done]")


if __name__ == "__main__":
    main()
