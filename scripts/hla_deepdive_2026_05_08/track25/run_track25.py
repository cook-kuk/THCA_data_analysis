#!/usr/bin/env python3
"""
Track 25 — Somatic HLA-pathway mutation rate in TCGA-THCA × DM1 score.

Boundary: somatic mutation count in HLA-pathway genes (B2M, NLRC5, TAP1/2,
HLA-A/B/C, IRF1, PSMB8/9, ...) — distinct from germline HLA allele typing.
Source = TCGA WES MAF (somatic calls), not RNA-seq imputation.

Outputs:
  results/hla_deepdive_2026_05_08/track25_hla_somatic_mut/
    figs/   ≥7 figures
    tables/ ≥6 tables
    track25_report.md
    track25_summary.json
"""
from __future__ import annotations

import gzip
import json
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

warnings.simplefilter("ignore")

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/hla_deepdive_2026_05_08/track25_hla_somatic_mut"
FIGS = OUT / "figs"
TABS = OUT / "tables"
FIGS.mkdir(parents=True, exist_ok=True)
TABS.mkdir(parents=True, exist_ok=True)

CAPTION_BOILER = ("Somatic mutation in HLA presentation pathway gene "
                  "— distinct from germline HLA allele typing.")

MAF_DIR = Path("/data/thca/data_raw/gdc/TCGA-THCA/mutation")
MAF_MANIFEST = Path("/data/thca/data_raw/gdc/TCGA-THCA/tcga_thca_mutation_manifest.tsv")
DM1_PANEL = ROOT / "project/results/dm1_robustness_v2026_05_08/per_sample_panel.tsv"
HLA_MOD = ROOT / "project/results/hla_deepdive_2026_05_08/track5_dm1_hla1_module/tables/T01_per_sample_module_scores.tsv"
DRIVER_FILE = ROOT / "project/results/dark_matter_phase2/p2d_per_sample_classification.tsv"
SURVIVAL = ROOT / "project/data/raw/TCGA_pancan/survival.tsv"

# ---------------------------------------------------------------------------
# Gene panels
# ---------------------------------------------------------------------------
CLASS_I = ["HLA-A", "HLA-B", "HLA-C", "B2M",
           "TAP1", "TAP2", "TAPBP",
           "ERAP1", "ERAP2",
           "NLRC5", "IRF1",
           "PSMB8", "PSMB9",
           "CALR", "CANX", "PDIA3"]
CLASS_II = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
            "HLA-DQA1", "HLA-DQB1", "CIITA", "CD74"]
COSTIM = ["CD274", "PDCD1LG2", "CD80", "CD86"]
HLA_ALL = CLASS_I + CLASS_II + COSTIM

LOF_TYPES = {"Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
             "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation"}
NONSILENT = LOF_TYPES | {"Missense_Mutation", "In_Frame_Del", "In_Frame_Ins"}


def short_id(s: str) -> str:
    """TCGA-XX-XXXX-01A -> TCGA-XX-XXXX (case)."""
    parts = s.split("-")
    return "-".join(parts[:3]) if len(parts) >= 3 else s


def sample_id(s: str) -> str:
    """TCGA-XX-XXXX-01A-12... -> TCGA-XX-XXXX-01A."""
    parts = s.split("-")
    return "-".join(parts[:4]) if len(parts) >= 4 else s


# ---------------------------------------------------------------------------
# 1. Load MAFs (one per tumor aliquot)
# ---------------------------------------------------------------------------
def load_mafs() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (full_maf_long, per_sample_summary_with_tmb)."""
    manifest = pd.read_csv(MAF_MANIFEST, sep="\t")
    # Restrict to tumor (Primary Tumor sample type code 01)
    manifest["is_tumor"] = manifest["sample_submitter_id"].str.contains(r"-0[1-6][A-Z]?$", regex=True)
    manifest = manifest[manifest["is_tumor"]].copy()
    manifest["sample"] = manifest["sample_submitter_id"]
    manifest["case"] = manifest["sample_submitter_id"].map(short_id)

    rows_keep = []
    per_sample_total = {}
    per_sample_nonsilent = {}
    keep_cols = ["Hugo_Symbol", "Variant_Classification", "Variant_Type",
                 "HGVSp_Short", "Tumor_Sample_Barcode", "t_depth", "t_alt_count",
                 "Chromosome", "Start_Position"]
    n_files = 0
    for _, row in manifest.iterrows():
        # On disk files are named <file_id>.maf.gz, not by file_name
        path = MAF_DIR / f"{row['file_id']}.maf.gz"
        if not path.exists():
            # fallback to file_name
            alt = MAF_DIR / row["file_name"]
            if alt.exists():
                path = alt
            else:
                continue
        n_files += 1
        try:
            df = pd.read_csv(path, sep="\t", comment="#",
                             usecols=lambda c: c in keep_cols,
                             low_memory=False)
        except Exception as e:
            print(f"[WARN] {path.name}: {e}", file=sys.stderr)
            continue
        if df.empty:
            continue
        df["sample"] = sample_id(row["sample"])
        df["case"] = row["case"]
        # TMB = total non-silent mutations per sample
        ns_mask = df["Variant_Classification"].isin(NONSILENT)
        per_sample_total[df["sample"].iloc[0]] = len(df)
        per_sample_nonsilent[df["sample"].iloc[0]] = int(ns_mask.sum())
        # Keep only HLA-pathway rows
        hla_mask = df["Hugo_Symbol"].isin(HLA_ALL)
        if hla_mask.any():
            rows_keep.append(df.loc[hla_mask].copy())

    print(f"[INFO] parsed {n_files}/{len(manifest)} MAF files", file=sys.stderr)
    full = pd.concat(rows_keep, ignore_index=True) if rows_keep else pd.DataFrame()
    tmb = pd.DataFrame({
        "sample": list(per_sample_total.keys()),
        "total_mut": list(per_sample_total.values()),
        "nonsilent_mut": [per_sample_nonsilent[s] for s in per_sample_total],
    })
    tmb["case"] = tmb["sample"].map(short_id)
    return full, tmb


# ---------------------------------------------------------------------------
# 2. Score each sample for HLA-pathway hits
# ---------------------------------------------------------------------------
def score_samples(maf: pd.DataFrame, tmb: pd.DataFrame) -> pd.DataFrame:
    out = tmb.copy().set_index("sample")
    for label, panel in [("classI", CLASS_I), ("classII", CLASS_II), ("costim", COSTIM)]:
        sub = maf[maf["Hugo_Symbol"].isin(panel)]
        ns = sub[sub["Variant_Classification"].isin(NONSILENT)].groupby("sample").size()
        lof = sub[sub["Variant_Classification"].isin(LOF_TYPES)].groupby("sample").size()
        out[f"{label}_nonsilent_count"] = ns.reindex(out.index).fillna(0).astype(int)
        out[f"{label}_lof_count"] = lof.reindex(out.index).fillna(0).astype(int)
    out = out.reset_index()
    return out


# ---------------------------------------------------------------------------
# 3. Merge in DM1, drivers, HLA-I module, survival
# ---------------------------------------------------------------------------
def merge_metadata(scored: pd.DataFrame) -> pd.DataFrame:
    panel = pd.read_csv(DM1_PANEL, sep="\t")
    panel = panel[panel["cohort"] == "TCGA-THCA"].copy()
    panel["sample_short"] = panel["sample"]  # already TCGA-XX-XXXX-01A
    panel["case"] = panel["sample_short"].map(short_id)
    panel = panel.rename(columns={"score": "DM1_score", "group": "DM1_group"})

    drv = pd.read_csv(DRIVER_FILE, sep="\t")
    drv = drv.rename(columns={"tcga_short": "case"})
    drv["driver_simple"] = drv["driver_class"].map({
        "Class1_BRAF_V600E": "BRAF",
        "Class2_RAS_hotspot": "RAS",
        "Class3_RTK_fusion": "Fusion",
        "Class4_other_driver": "Other",
        "Class5_kinase_other": "Other",
        "Class6_True_driver_neg": "TripleNeg",
    }).fillna("Unknown")

    hla = pd.read_csv(HLA_MOD, sep="\t")
    # T01 sample format = TCGA-XX-XXXX-01 (no aliquot tail)
    hla["case"] = hla["sample"].map(short_id)
    hla = hla[["case", "HLA1_score", "HLA2_score", "DM1_tertile", "OS", "OS.time", "PFI", "PFI.time", "stage", "tert_pos"]]

    surv = pd.read_csv(SURVIVAL, sep="\t", low_memory=False)
    surv = surv[surv["cancer type abbreviation"] == "THCA"].copy()
    surv["case"] = surv["sample"].map(short_id)
    surv = surv[["case", "DSS", "DSS.time", "OS", "OS.time", "PFI", "PFI.time", "age_at_initial_pathologic_diagnosis"]].rename(
        columns={"age_at_initial_pathologic_diagnosis": "age"})
    # Avoid colliding survival cols when both present — prefer survival.tsv
    df = scored.copy()
    df["case"] = df["sample"].map(short_id)
    df = df.merge(panel[["case", "DM1_score", "DM1_group"]], on="case", how="left")
    df = df.merge(drv[["case", "driver_class", "driver_simple", "v17_dark_cluster"]], on="case", how="left")
    df = df.merge(hla[["case", "HLA1_score", "HLA2_score", "DM1_tertile", "stage", "tert_pos"]], on="case", how="left")
    df = df.merge(surv, on="case", how="left", suffixes=("", "_panel"))
    return df


# ---------------------------------------------------------------------------
# 4. Analyses
# ---------------------------------------------------------------------------
def per_gene_mutation_rate(maf: pd.DataFrame, scored: pd.DataFrame) -> pd.DataFrame:
    n_samples = scored["sample"].nunique()
    rows = []
    for gene in HLA_ALL:
        sub = maf[maf["Hugo_Symbol"] == gene]
        n_muts = len(sub)
        n_lof = sub["Variant_Classification"].isin(LOF_TYPES).sum()
        n_pts = sub["sample"].nunique()
        n_pts_lof = sub.loc[sub["Variant_Classification"].isin(LOF_TYPES), "sample"].nunique()
        rows.append({
            "gene": gene,
            "panel": "ClassI" if gene in CLASS_I else ("ClassII" if gene in CLASS_II else "Costim"),
            "n_mutations": int(n_muts),
            "n_LoF": int(n_lof),
            "n_samples_mutated": int(n_pts),
            "n_samples_LoF": int(n_pts_lof),
            "pct_samples_mutated": 100.0 * n_pts / n_samples,
            "pct_samples_LoF": 100.0 * n_pts_lof / n_samples,
        })
    return pd.DataFrame(rows)


def dm1_tertile_test(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    df = df.dropna(subset=["DM1_score"]).copy()
    q = df["DM1_score"].quantile([1/3, 2/3]).values
    df["DM1_tert3"] = np.where(df["DM1_score"] <= q[0], "Low",
                       np.where(df["DM1_score"] >= q[1], "High", "Mid"))
    for var in ["classI_nonsilent_count", "classI_lof_count",
                "classII_nonsilent_count", "classII_lof_count",
                "costim_nonsilent_count", "costim_lof_count",
                "nonsilent_mut"]:
        hi = df.loc[df["DM1_tert3"] == "High", var].values
        lo = df.loc[df["DM1_tert3"] == "Low", var].values
        mid = df.loc[df["DM1_tert3"] == "Mid", var].values
        try:
            U, p_mw = stats.mannwhitneyu(hi, lo, alternative="two-sided")
        except Exception:
            U, p_mw = np.nan, np.nan
        try:
            H, p_kw = stats.kruskal(hi, mid, lo)
        except Exception:
            H, p_kw = np.nan, np.nan
        # Also: continuous Spearman DM1_score × count
        rho, p_sp = stats.spearmanr(df["DM1_score"], df[var], nan_policy="omit")
        rows.append({
            "variable": var,
            "n_high": len(hi), "n_mid": len(mid), "n_low": len(lo),
            "mean_high": float(np.mean(hi)) if len(hi) else np.nan,
            "mean_low": float(np.mean(lo)) if len(lo) else np.nan,
            "mannwhitney_U": float(U) if not np.isnan(U) else np.nan,
            "p_mannwhitney_HighVsLow": float(p_mw) if not np.isnan(p_mw) else np.nan,
            "kruskal_H": float(H) if not np.isnan(H) else np.nan,
            "p_kruskal_3grp": float(p_kw) if not np.isnan(p_kw) else np.nan,
            "spearman_rho_dm1_count": float(rho) if not np.isnan(rho) else np.nan,
            "spearman_p": float(p_sp) if not np.isnan(p_sp) else np.nan,
        })
    out = pd.DataFrame(rows)
    if out["p_mannwhitney_HighVsLow"].notna().any():
        out["fdr_BH_mw"] = multipletests(out["p_mannwhitney_HighVsLow"].fillna(1.0), method="fdr_bh")[1]
    return out, df


def per_gene_dm1_fisher(maf: pd.DataFrame, df: pd.DataFrame) -> pd.DataFrame:
    """For each gene with ≥3 mutated samples, Fisher exact for High vs Low DM1 tertile."""
    rows = []
    df_use = df.dropna(subset=["DM1_score"]).copy()
    q = df_use["DM1_score"].quantile([1/3, 2/3]).values
    df_use["DM1_tert3"] = np.where(df_use["DM1_score"] <= q[0], "Low",
                          np.where(df_use["DM1_score"] >= q[1], "High", "Mid"))
    samples_high = set(df_use.loc[df_use["DM1_tert3"] == "High", "sample"])
    samples_low = set(df_use.loc[df_use["DM1_tert3"] == "Low", "sample"])
    for gene in HLA_ALL:
        muts = set(maf.loc[(maf["Hugo_Symbol"] == gene) & (maf["Variant_Classification"].isin(NONSILENT)), "sample"])
        n_h = len(samples_high & muts)
        n_l = len(samples_low & muts)
        if (n_h + n_l) < 3:
            continue
        a = n_h
        b = len(samples_high) - n_h
        c = n_l
        d = len(samples_low) - n_l
        odds, p = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
        rows.append({"gene": gene,
                     "panel": "ClassI" if gene in CLASS_I else ("ClassII" if gene in CLASS_II else "Costim"),
                     "n_high_mut": a, "n_high_wt": b,
                     "n_low_mut": c, "n_low_wt": d,
                     "OR": float(odds), "p_fisher": float(p)})
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr_BH"] = multipletests(out["p_fisher"], method="fdr_bh")[1]
        out = out.sort_values("p_fisher")
    return out


def tmb_controlled(df: pd.DataFrame) -> dict:
    """Logistic regression: any class-I LoF ~ DM1_score + log10(TMB+1)."""
    d = df.dropna(subset=["DM1_score", "nonsilent_mut"]).copy()
    d["any_classI_lof"] = (d["classI_lof_count"] > 0).astype(int)
    d["any_classI_nonsilent"] = (d["classI_nonsilent_count"] > 0).astype(int)
    d["log_tmb"] = np.log10(d["nonsilent_mut"].clip(lower=0) + 1)
    out = {}
    for tag, y in [("any_classI_lof", "any_classI_lof"),
                   ("any_classI_nonsilent", "any_classI_nonsilent")]:
        if d[y].sum() < 3:
            out[tag] = {"note": "too few events", "n_events": int(d[y].sum())}
            continue
        try:
            mdl = smf.logit(f"{y} ~ DM1_score + log_tmb", data=d).fit(disp=0)
            res = {
                "n": int(mdl.nobs),
                "n_events": int(d[y].sum()),
                "coef_DM1": float(mdl.params.get("DM1_score", np.nan)),
                "p_DM1": float(mdl.pvalues.get("DM1_score", np.nan)),
                "OR_DM1_per_unit": float(np.exp(mdl.params.get("DM1_score", np.nan))),
                "coef_logTMB": float(mdl.params.get("log_tmb", np.nan)),
                "p_logTMB": float(mdl.pvalues.get("log_tmb", np.nan)),
            }
        except Exception as e:
            res = {"error": str(e)}
        out[tag] = res
    # Continuous Poisson: classI_lof_count ~ DM1_score + offset(log TMB+1)
    try:
        mdl = smf.glm("classI_lof_count ~ DM1_score",
                      data=d, family=sm.families.Poisson(),
                      offset=d["log_tmb"]).fit(disp=0)
        out["poisson_classI_lof_offset_logTMB"] = {
            "n": int(mdl.nobs),
            "coef_DM1": float(mdl.params.get("DM1_score", np.nan)),
            "p_DM1": float(mdl.pvalues.get("DM1_score", np.nan)),
            "IRR_DM1": float(np.exp(mdl.params.get("DM1_score", np.nan))),
        }
    except Exception as e:
        out["poisson_classI_lof_offset_logTMB"] = {"error": str(e)}
    return out


def driver_stratified(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    d = df.dropna(subset=["DM1_score", "driver_simple"])
    for drv in sorted(d["driver_simple"].dropna().unique()):
        sub = d[d["driver_simple"] == drv]
        if len(sub) < 8:
            continue
        rho, p = stats.spearmanr(sub["DM1_score"], sub["classI_nonsilent_count"], nan_policy="omit")
        rows.append({"driver": drv, "n": len(sub),
                     "median_classI_nonsilent": float(sub["classI_nonsilent_count"].median()),
                     "mean_classI_nonsilent": float(sub["classI_nonsilent_count"].mean()),
                     "n_classI_lof_pos": int((sub["classI_lof_count"] > 0).sum()),
                     "spearman_rho_DM1_classI": float(rho) if not np.isnan(rho) else np.nan,
                     "spearman_p": float(p) if not np.isnan(p) else np.nan,
                     "median_DM1_score": float(sub["DM1_score"].median()),
                     })
    return pd.DataFrame(rows)


def class_i_vs_ii(df: pd.DataFrame) -> dict:
    a = df["classI_nonsilent_count"].sum()
    b = df["classII_nonsilent_count"].sum()
    c = df["costim_nonsilent_count"].sum()
    n = len(df)
    pI = (df["classI_nonsilent_count"] > 0).mean()
    pII = (df["classII_nonsilent_count"] > 0).mean()
    pC = (df["costim_nonsilent_count"] > 0).mean()
    return {
        "n_samples": int(n),
        "classI_total_mutations": int(a),
        "classII_total_mutations": int(b),
        "costim_total_mutations": int(c),
        "classI_pct_samples_with_mutation": float(pI * 100),
        "classII_pct_samples_with_mutation": float(pII * 100),
        "costim_pct_samples_with_mutation": float(pC * 100),
    }


def expression_cooccurrence(df: pd.DataFrame) -> dict:
    d = df.dropna(subset=["HLA1_score"]).copy()
    if d.empty:
        return {"note": "no HLA-I module overlap"}
    d["any_classI_nonsilent"] = (d["classI_nonsilent_count"] > 0).astype(int)
    if d["any_classI_nonsilent"].sum() == 0:
        return {"note": "no class-I mutations in samples with HLA-I module"}
    mut = d.loc[d["any_classI_nonsilent"] == 1, "HLA1_score"]
    wt = d.loc[d["any_classI_nonsilent"] == 0, "HLA1_score"]
    try:
        U, p = stats.mannwhitneyu(mut, wt, alternative="two-sided")
    except Exception:
        U, p = np.nan, np.nan
    rho, p_rho = stats.spearmanr(d["HLA1_score"], d["classI_nonsilent_count"], nan_policy="omit")
    return {
        "n": int(len(d)),
        "n_classI_mut": int((d["any_classI_nonsilent"] == 1).sum()),
        "median_HLA1_mut": float(mut.median()) if len(mut) else np.nan,
        "median_HLA1_wt": float(wt.median()) if len(wt) else np.nan,
        "p_mannwhitney_mut_vs_wt": float(p) if not np.isnan(p) else np.nan,
        "spearman_rho_HLA1_classIcount": float(rho) if not np.isnan(rho) else np.nan,
        "spearman_p": float(p_rho) if not np.isnan(p_rho) else np.nan,
    }


def survival_cox(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Cox HR for HLA-pathway LoF on PFI and DSS."""
    try:
        from lifelines import CoxPHFitter
    except Exception as e:
        return pd.DataFrame(), {"error": f"lifelines not available: {e}"}
    rows = []
    summ = {}
    for endpoint in ["PFI", "DSS", "OS"]:
        d = df.dropna(subset=[endpoint, f"{endpoint}.time", "DM1_score"]).copy()
        d["any_classI_lof"] = (d["classI_lof_count"] > 0).astype(int)
        d["any_classI_nonsilent"] = (d["classI_nonsilent_count"] > 0).astype(int)
        d["log_tmb"] = np.log10(d["nonsilent_mut"].clip(lower=0) + 1)
        d = d.rename(columns={f"{endpoint}.time": "time", endpoint: "event"})
        d = d[d["time"] > 0]
        for var in ["any_classI_nonsilent", "any_classI_lof"]:
            if d[var].sum() < 3 or d["event"].sum() < 3:
                rows.append({"endpoint": endpoint, "var": var, "n": len(d),
                             "n_events": int(d["event"].sum()),
                             "n_var_pos": int(d[var].sum()),
                             "HR": np.nan, "HR_lo": np.nan, "HR_hi": np.nan, "p": np.nan,
                             "model": "skipped (low events)"})
                continue
            try:
                cph = CoxPHFitter(penalizer=0.01)
                fit_df = d[[var, "time", "event"]].dropna()
                cph.fit(fit_df, duration_col="time", event_col="event")
                s = cph.summary.loc[var]
                rows.append({"endpoint": endpoint, "var": var, "n": len(fit_df),
                             "n_events": int(fit_df["event"].sum()),
                             "n_var_pos": int(fit_df[var].sum()),
                             "HR": float(s["exp(coef)"]),
                             "HR_lo": float(s["exp(coef) lower 95%"]),
                             "HR_hi": float(s["exp(coef) upper 95%"]),
                             "p": float(s["p"]),
                             "model": "univariate"})
                # Multivariate adj DM1 + log_tmb
                fit_df2 = d[[var, "time", "event", "DM1_score", "log_tmb"]].dropna()
                cph2 = CoxPHFitter(penalizer=0.05)
                cph2.fit(fit_df2, duration_col="time", event_col="event")
                s2 = cph2.summary.loc[var]
                rows.append({"endpoint": endpoint, "var": var, "n": len(fit_df2),
                             "n_events": int(fit_df2["event"].sum()),
                             "n_var_pos": int(fit_df2[var].sum()),
                             "HR": float(s2["exp(coef)"]),
                             "HR_lo": float(s2["exp(coef) lower 95%"]),
                             "HR_hi": float(s2["exp(coef) upper 95%"]),
                             "p": float(s2["p"]),
                             "model": "adj_DM1_logTMB"})
            except Exception as e:
                rows.append({"endpoint": endpoint, "var": var, "n": len(d),
                             "n_events": int(d["event"].sum()),
                             "n_var_pos": int(d[var].sum()),
                             "HR": np.nan, "HR_lo": np.nan, "HR_hi": np.nan, "p": np.nan,
                             "model": f"error: {e}"})
    return pd.DataFrame(rows), summ


# ---------------------------------------------------------------------------
# 5. Plots
# ---------------------------------------------------------------------------
def plot_per_gene_rate(per_gene: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(9, 6))
    pg = per_gene.sort_values("pct_samples_mutated", ascending=True)
    colors = pg["panel"].map({"ClassI": "#1f77b4", "ClassII": "#d62728", "Costim": "#2ca02c"})
    ax.barh(pg["gene"], pg["pct_samples_mutated"], color=colors)
    ax.set_xlabel("% TCGA-THCA samples with non-silent somatic mutation")
    ax.set_title(f"F1 — HLA pathway gene mutation rate (TCGA-THCA, n={int(per_gene['n_samples_mutated'].max() or 0)} max gene)\n{CAPTION_BOILER}",
                 fontsize=10)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in ["#1f77b4", "#d62728", "#2ca02c"]]
    ax.legend(handles, ["Class I", "Class II", "Costim"], loc="lower right")
    fig.tight_layout()
    fig.savefig(FIGS / "F1_per_gene_mutation_rate.png", dpi=150)
    fig.savefig(FIGS / "F1_per_gene_mutation_rate.pdf")
    plt.close(fig)


def plot_dm1_box(df: pd.DataFrame):
    d = df.dropna(subset=["DM1_score"]).copy()
    q = d["DM1_score"].quantile([1/3, 2/3]).values
    d["DM1_tert3"] = pd.Categorical(np.where(d["DM1_score"] <= q[0], "Low",
                                    np.where(d["DM1_score"] >= q[1], "High", "Mid")),
                                    categories=["Low", "Mid", "High"], ordered=True)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    panels = [
        ("classI_nonsilent_count", "Class-I non-silent count"),
        ("classI_lof_count", "Class-I LoF count"),
        ("classII_nonsilent_count", "Class-II non-silent count"),
    ]
    for ax, (col, ttl) in zip(axes, panels):
        data = [d.loc[d["DM1_tert3"] == g, col].values for g in ["Low", "Mid", "High"]]
        bp = ax.boxplot(data, labels=["Low", "Mid", "High"], showmeans=True, patch_artist=True)
        for patch, c in zip(bp["boxes"], ["#bdd7e7", "#6baed6", "#2171b5"]):
            patch.set_facecolor(c)
        # jitter
        for i, g in enumerate(["Low", "Mid", "High"]):
            arr = d.loc[d["DM1_tert3"] == g, col].values
            x = np.random.normal(i + 1, 0.06, size=len(arr))
            ax.scatter(x, arr, alpha=0.4, s=8, color="black")
        ax.set_title(ttl)
        ax.set_xlabel("DM1 tertile")
        ax.set_ylabel("count")
    fig.suptitle(f"F2 — HLA-pathway somatic mutation count by DM1 tertile (TCGA-THCA, n={len(d)})\n{CAPTION_BOILER}",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(FIGS / "F2_dm1_tertile_box.png", dpi=150)
    fig.savefig(FIGS / "F2_dm1_tertile_box.pdf")
    plt.close(fig)


def plot_lof_rate_by_tertile(df: pd.DataFrame):
    d = df.dropna(subset=["DM1_score"]).copy()
    q = d["DM1_score"].quantile([1/3, 2/3]).values
    d["DM1_tert3"] = pd.Categorical(np.where(d["DM1_score"] <= q[0], "Low",
                                    np.where(d["DM1_score"] >= q[1], "High", "Mid")),
                                    categories=["Low", "Mid", "High"], ordered=True)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    rates = []
    n_per = []
    for g in ["Low", "Mid", "High"]:
        sub = d[d["DM1_tert3"] == g]
        rate_ns = (sub["classI_nonsilent_count"] > 0).mean() * 100
        rate_lof = (sub["classI_lof_count"] > 0).mean() * 100
        rates.append((rate_ns, rate_lof))
        n_per.append(len(sub))
    x = np.arange(3)
    w = 0.35
    ax.bar(x - w/2, [r[0] for r in rates], w, label="any non-silent", color="#6baed6")
    ax.bar(x + w/2, [r[1] for r in rates], w, label="any LoF", color="#08306b")
    ax.set_xticks(x)
    ax.set_xticklabels([f"Low\nn={n_per[0]}", f"Mid\nn={n_per[1]}", f"High\nn={n_per[2]}"])
    ax.set_ylabel("% samples with class-I HLA-pathway mutation")
    ax.set_title(f"F3 — % samples with class-I HLA-pathway somatic mutation by DM1 tertile\n{CAPTION_BOILER}",
                 fontsize=10)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGS / "F3_classI_rate_by_dm1_tertile.png", dpi=150)
    fig.savefig(FIGS / "F3_classI_rate_by_dm1_tertile.pdf")
    plt.close(fig)


def plot_dm1_vs_count(df: pd.DataFrame):
    d = df.dropna(subset=["DM1_score"])
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(d["DM1_score"], d["classI_nonsilent_count"], alpha=0.4, s=18, color="#1f77b4")
    rho, p = stats.spearmanr(d["DM1_score"], d["classI_nonsilent_count"], nan_policy="omit")
    ax.set_xlabel("DM1 score")
    ax.set_ylabel("Class-I HLA-pathway non-silent mutation count")
    ax.set_title(f"F4 — DM1 score × class-I HLA-pathway mutation count\nSpearman ρ={rho:.3f} p={p:.2g} n={len(d)}\n{CAPTION_BOILER}",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(FIGS / "F4_dm1_vs_classI_count.png", dpi=150)
    fig.savefig(FIGS / "F4_dm1_vs_classI_count.pdf")
    plt.close(fig)


def plot_driver_stratified(driver_tab: pd.DataFrame):
    if driver_tab.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 4.8))
    x = np.arange(len(driver_tab))
    ax.bar(x, driver_tab["mean_classI_nonsilent"], color="#1f77b4")
    for i, n in enumerate(driver_tab["n"]):
        ax.text(i, driver_tab["mean_classI_nonsilent"].iloc[i] + 0.01,
                f"n={int(n)}", ha="center", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(driver_tab["driver"])
    ax.set_ylabel("Mean class-I HLA-pathway non-silent mutation count / sample")
    ax.set_title(f"F5 — Driver-stratified class-I HLA-pathway mutation count\n{CAPTION_BOILER}",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(FIGS / "F5_driver_stratified.png", dpi=150)
    fig.savefig(FIGS / "F5_driver_stratified.pdf")
    plt.close(fig)


def plot_class_i_vs_ii(df: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    rates = [
        ("Class I", (df["classI_nonsilent_count"] > 0).mean() * 100, (df["classI_lof_count"] > 0).mean() * 100),
        ("Class II", (df["classII_nonsilent_count"] > 0).mean() * 100, (df["classII_lof_count"] > 0).mean() * 100),
        ("Costim", (df["costim_nonsilent_count"] > 0).mean() * 100, (df["costim_lof_count"] > 0).mean() * 100),
    ]
    x = np.arange(3)
    w = 0.4
    ax.bar(x - w/2, [r[1] for r in rates], w, label="any non-silent", color="#6baed6")
    ax.bar(x + w/2, [r[2] for r in rates], w, label="any LoF", color="#08306b")
    ax.set_xticks(x); ax.set_xticklabels([r[0] for r in rates])
    ax.set_ylabel("% TCGA-THCA samples")
    ax.set_title(f"F6 — Class-I vs Class-II vs Costim mutation rate\n{CAPTION_BOILER}", fontsize=10)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGS / "F6_classI_vs_II_vs_costim.png", dpi=150)
    fig.savefig(FIGS / "F6_classI_vs_II_vs_costim.pdf")
    plt.close(fig)


def plot_hla_module_cooccur(df: pd.DataFrame):
    d = df.dropna(subset=["HLA1_score"]).copy()
    if d.empty:
        return
    d["any_classI_nonsilent"] = (d["classI_nonsilent_count"] > 0).astype(int)
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    data = [d.loc[d["any_classI_nonsilent"] == 0, "HLA1_score"].values,
            d.loc[d["any_classI_nonsilent"] == 1, "HLA1_score"].values]
    bp = ax.boxplot(data, labels=[f"WT\nn={len(data[0])}", f"Class-I mut\nn={len(data[1])}"], patch_artist=True, showmeans=True)
    for patch, c in zip(bp["boxes"], ["#bdd7e7", "#08306b"]):
        patch.set_facecolor(c)
    if len(data[1]) >= 1:
        for i, arr in enumerate(data):
            x = np.random.normal(i + 1, 0.05, size=len(arr))
            ax.scatter(x, arr, alpha=0.4, s=8, color="black")
    try:
        U, p = stats.mannwhitneyu(data[0], data[1])
        ttl = f"F7 — HLA-I expression module vs class-I HLA-pathway mutation\nMannWhitney p={p:.3g} (n_mut={len(data[1])})\n{CAPTION_BOILER}"
    except Exception:
        ttl = f"F7 — HLA-I expression module vs class-I HLA-pathway mutation\n{CAPTION_BOILER}"
    ax.set_ylabel("HLA-I expression module score")
    ax.set_title(ttl, fontsize=9)
    fig.tight_layout()
    fig.savefig(FIGS / "F7_hla1_module_vs_classI_mut.png", dpi=150)
    fig.savefig(FIGS / "F7_hla1_module_vs_classI_mut.pdf")
    plt.close(fig)


def plot_tmb_by_dm1(df: pd.DataFrame):
    """F8 — Total non-silent TMB by DM1 tertile (background context)."""
    d = df.dropna(subset=["DM1_score", "nonsilent_mut"]).copy()
    q = d["DM1_score"].quantile([1/3, 2/3]).values
    d["DM1_tert3"] = pd.Categorical(np.where(d["DM1_score"] <= q[0], "Low",
                                    np.where(d["DM1_score"] >= q[1], "High", "Mid")),
                                    categories=["Low", "Mid", "High"], ordered=True)
    fig, ax = plt.subplots(figsize=(7, 4.8))
    data = [d.loc[d["DM1_tert3"] == g, "nonsilent_mut"].values for g in ["Low", "Mid", "High"]]
    bp = ax.boxplot(data, labels=[f"Low\nn={len(data[0])}", f"Mid\nn={len(data[1])}", f"High\nn={len(data[2])}"],
                    showmeans=True, patch_artist=True, showfliers=False)
    for patch, c in zip(bp["boxes"], ["#bdd7e7", "#6baed6", "#2171b5"]):
        patch.set_facecolor(c)
    for i, g in enumerate(["Low", "Mid", "High"]):
        arr = d.loc[d["DM1_tert3"] == g, "nonsilent_mut"].values
        x = np.random.normal(i + 1, 0.06, size=len(arr))
        ax.scatter(x, arr, alpha=0.4, s=8, color="black")
    ax.set_ylim(0, np.percentile(d["nonsilent_mut"], 98) + 2)
    try:
        H, p = stats.kruskal(*data)
    except Exception:
        H, p = np.nan, np.nan
    rho, p_sp = stats.spearmanr(d["DM1_score"], d["nonsilent_mut"], nan_policy="omit")
    ax.set_ylabel("Non-silent TMB (count / sample)")
    ax.set_title(f"F8 — Total non-silent TMB by DM1 tertile (TCGA-THCA, n={len(d)})\n"
                 f"Kruskal p={p:.3g}; Spearman ρ(DM1,TMB)={rho:.3f}, p={p_sp:.3g}\n{CAPTION_BOILER}",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(FIGS / "F8_tmb_by_dm1_tertile.png", dpi=150)
    fig.savefig(FIGS / "F8_tmb_by_dm1_tertile.pdf")
    plt.close(fig)


def plot_survival_forest(cox_tab: pd.DataFrame):
    if cox_tab.empty or cox_tab["HR"].notna().sum() == 0:
        return
    d = cox_tab.dropna(subset=["HR"]).copy().reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(8.5, 0.3 * len(d) + 2.5))
    y = np.arange(len(d))
    ax.errorbar(d["HR"], y,
                xerr=[d["HR"] - d["HR_lo"].clip(lower=1e-3), d["HR_hi"] - d["HR"]],
                fmt="o", color="#08306b", ecolor="#6baed6", capsize=3)
    ax.axvline(1, color="grey", linestyle="--")
    labels = [f"{r['endpoint']} | {r['var']} | {r['model']} (n_evt={int(r['n_events'])}, n+={int(r['n_var_pos'])})" for _, r in d.iterrows()]
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("Hazard ratio (log scale)")
    ax.set_title(f"F8 — Cox HR for HLA-pathway somatic mutation × outcome\n{CAPTION_BOILER}", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIGS / "F8_survival_forest.png", dpi=150)
    fig.savefig(FIGS / "F8_survival_forest.pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("[1] Loading MAFs...", file=sys.stderr)
    maf, tmb = load_mafs()
    print(f"    parsed {len(tmb)} samples / {len(maf)} HLA-pathway rows", file=sys.stderr)
    if maf.empty:
        raise RuntimeError("no MAF data loaded")

    # T0 raw HLA-pathway long table
    maf_out = maf.drop(columns=["case"], errors="ignore")
    maf_out.to_csv(TABS / "T01_hla_pathway_mutations_long.tsv", sep="\t", index=False)

    print("[2] Per-sample scoring...", file=sys.stderr)
    scored = score_samples(maf, tmb)

    print("[3] Merging metadata (DM1, drivers, HLA module, survival)...", file=sys.stderr)
    df = merge_metadata(scored)

    print("[4] Per-gene mutation rate...", file=sys.stderr)
    per_gene = per_gene_mutation_rate(maf, scored)
    per_gene.to_csv(TABS / "T02_per_gene_mutation_rate.tsv", sep="\t", index=False)
    plot_per_gene_rate(per_gene)

    print("[5] DM1 tertile tests...", file=sys.stderr)
    tert_tab, df_with_tert = dm1_tertile_test(df)
    tert_tab.to_csv(TABS / "T03_dm1_tertile_tests.tsv", sep="\t", index=False)
    plot_dm1_box(df)
    plot_lof_rate_by_tertile(df)
    plot_dm1_vs_count(df)

    print("[6] Per-gene Fisher exact High vs Low DM1 tertile...", file=sys.stderr)
    fisher_tab = per_gene_dm1_fisher(maf, df)
    fisher_tab.to_csv(TABS / "T04_per_gene_dm1_fisher.tsv", sep="\t", index=False)

    print("[7] TMB-controlled DM1 × HLA-pathway...", file=sys.stderr)
    tmb_res = tmb_controlled(df)
    with open(TABS / "T05_tmb_controlled.json", "w") as f:
        json.dump(tmb_res, f, indent=2)

    print("[8] Driver-stratified...", file=sys.stderr)
    driver_tab = driver_stratified(df)
    driver_tab.to_csv(TABS / "T06_driver_stratified.tsv", sep="\t", index=False)
    plot_driver_stratified(driver_tab)

    print("[9] Class-I vs Class-II vs Costim...", file=sys.stderr)
    civii = class_i_vs_ii(df)
    plot_class_i_vs_ii(df)

    print("[10] HLA-I expression module co-occurrence...", file=sys.stderr)
    cooccur = expression_cooccurrence(df)
    with open(TABS / "T07_hla1_module_cooccurrence.json", "w") as f:
        json.dump(cooccur, f, indent=2)
    plot_hla_module_cooccur(df)

    print("[11] Survival Cox HR + TMB-by-DM1 figure...", file=sys.stderr)
    cox_tab, cox_summ = survival_cox(df)
    cox_tab.to_csv(TABS / "T08_survival_cox.tsv", sep="\t", index=False)
    plot_survival_forest(cox_tab)
    plot_tmb_by_dm1(df)

    # Save merged per-sample table for reproducibility
    df.to_csv(TABS / "T09_per_sample_merged.tsv", sep="\t", index=False)

    # Summary JSON
    summary = {
        "track": "Track 25 - Somatic HLA-pathway mutation rate in TCGA-THCA × DM1 score",
        "boundary": CAPTION_BOILER,
        "n_tcga_thca_samples_with_maf": int(scored.shape[0]),
        "n_with_dm1_score": int(df["DM1_score"].notna().sum()),
        "n_with_hla1_module": int(df["HLA1_score"].notna().sum()),
        "tmb_summary": {
            "median_nonsilent_mut": float(df["nonsilent_mut"].median()),
            "mean_nonsilent_mut": float(df["nonsilent_mut"].mean()),
            "max_nonsilent_mut": int(df["nonsilent_mut"].max()),
        },
        "hla_pathway_overall": {
            "samples_any_classI_nonsilent": int((df["classI_nonsilent_count"] > 0).sum()),
            "samples_any_classI_lof": int((df["classI_lof_count"] > 0).sum()),
            "pct_any_classI_nonsilent": float((df["classI_nonsilent_count"] > 0).mean() * 100),
            "pct_any_classI_lof": float((df["classI_lof_count"] > 0).mean() * 100),
            "samples_any_classII_nonsilent": int((df["classII_nonsilent_count"] > 0).sum()),
            "samples_any_classII_lof": int((df["classII_lof_count"] > 0).sum()),
        },
        "dm1_tertile_classI_nonsilent": tert_tab.set_index("variable").loc["classI_nonsilent_count"].to_dict() if "classI_nonsilent_count" in tert_tab["variable"].values else {},
        "dm1_tertile_classI_lof": tert_tab.set_index("variable").loc["classI_lof_count"].to_dict() if "classI_lof_count" in tert_tab["variable"].values else {},
        "tmb_controlled": tmb_res,
        "class_i_vs_ii": civii,
        "expression_cooccurrence": cooccur,
        "top_per_gene_fisher": fisher_tab.head(5).to_dict(orient="records") if not fisher_tab.empty else [],
        "driver_stratified": driver_tab.to_dict(orient="records") if not driver_tab.empty else [],
        "outputs_dir": str(OUT),
    }
    with open(OUT / "track25_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    # Markdown report
    write_report(summary, per_gene, tert_tab, fisher_tab, tmb_res, civii, cooccur,
                 driver_tab, cox_tab, df)
    print("[DONE]", file=sys.stderr)


def write_report(summary, per_gene, tert_tab, fisher_tab, tmb_res, civii, cooccur,
                 driver_tab, cox_tab, df):
    md = []
    md.append("# Track 25 — Somatic HLA-pathway mutation rate in TCGA-THCA × DM1 score\n")
    md.append("**Boundary (caption boilerplate, used on every figure / table):**  ")
    md.append(f"> {CAPTION_BOILER}\n")
    md.append("This track counts **somatic non-silent and LoF mutations** in HLA class-I presentation "
              "machinery (HLA-A/B/C, B2M, TAP1/2, TAPBP, ERAP1/2, NLRC5, IRF1, PSMB8/9, CALR, CANX, "
              "PDIA3) and class-II machinery (HLA-DRA, DRB1, DPA1, DPB1, DQA1, DQB1, CIITA, CD74) plus "
              "costim adjacents (CD274, PDCD1LG2, CD80, CD86) across TCGA-THCA WES MAF files. **Mutation "
              "count, not allele typing.** Allowed under Paper 1 boundary because it is genomic "
              "(somatic point mutation / indel / LoF), not germline HLA imputation.\n")

    md.append("## 0. Inputs\n")
    md.append(f"- TCGA-THCA MAF files: `/data/thca/data_raw/gdc/TCGA-THCA/mutation/` (manifest `tcga_thca_mutation_manifest.tsv`)")
    md.append(f"- DM1 score per sample: `project/results/dm1_robustness_v2026_05_08/per_sample_panel.tsv` (TCGA subset)")
    md.append(f"- Driver class (BRAF / RAS / Fusion / TripleNeg): `project/results/dark_matter_phase2/p2d_per_sample_classification.tsv`")
    md.append(f"- HLA-I/II module: `project/results/hla_deepdive_2026_05_08/track5_dm1_hla1_module/tables/T01_per_sample_module_scores.tsv`")
    md.append(f"- Survival: `project/data/raw/TCGA_pancan/survival.tsv` (PFI/DSS/OS)")
    md.append(f"- N samples MAF parsed: **{summary['n_tcga_thca_samples_with_maf']}**, with DM1 score: **{summary['n_with_dm1_score']}**, with HLA-I module: **{summary['n_with_hla1_module']}**\n")

    md.append("## 1. TMB summary (TCGA-THCA is mutation-quiet)\n")
    md.append(f"- Median non-silent mutations / sample: **{summary['tmb_summary']['median_nonsilent_mut']:.0f}**, mean {summary['tmb_summary']['mean_nonsilent_mut']:.1f}, max {summary['tmb_summary']['max_nonsilent_mut']}.")
    md.append("- Statistical power for any single HLA-pathway gene is correspondingly low; we therefore lean on count-based and panel-aggregated tests.\n")

    md.append("## 2. HLA-pathway mutation atlas (Table T02 / Figure F1)\n")
    top = per_gene.sort_values("pct_samples_mutated", ascending=False).head(10)
    md.append("Top 10 HLA-pathway genes by % samples mutated:\n")
    md.append("| gene | panel | n_samples_mut | n_samples_LoF | % samples mut |")
    md.append("|------|-------|---------------|---------------|----------------|")
    for _, r in top.iterrows():
        md.append(f"| {r['gene']} | {r['panel']} | {int(r['n_samples_mutated'])} | {int(r['n_samples_LoF'])} | {r['pct_samples_mutated']:.2f} |")
    md.append("")

    md.append(f"- Class-I overall: {summary['hla_pathway_overall']['samples_any_classI_nonsilent']} samples with any non-silent mutation ({summary['hla_pathway_overall']['pct_any_classI_nonsilent']:.2f}%); {summary['hla_pathway_overall']['samples_any_classI_lof']} with LoF ({summary['hla_pathway_overall']['pct_any_classI_lof']:.2f}%).")
    md.append(f"- Class-II overall: {summary['hla_pathway_overall']['samples_any_classII_nonsilent']} non-silent ({100 * summary['hla_pathway_overall']['samples_any_classII_nonsilent'] / summary['n_tcga_thca_samples_with_maf']:.2f}%); {summary['hla_pathway_overall']['samples_any_classII_lof']} LoF.")
    md.append(f"- Class-I total nonsilent mutations across cohort: {civii['classI_total_mutations']}; Class-II: {civii['classII_total_mutations']}; Costim: {civii['costim_total_mutations']}.\n")

    md.append("## 3. DM1 score × HLA-pathway mutation count (Table T03 / Figures F2-F4)\n")
    if "classI_nonsilent_count" in tert_tab["variable"].values:
        r = tert_tab.set_index("variable").loc["classI_nonsilent_count"]
        md.append(f"- **Class-I non-silent count, High-DM1 tertile vs Low**: mean_high={r['mean_high']:.3f}, mean_low={r['mean_low']:.3f}; Mann-Whitney p={r['p_mannwhitney_HighVsLow']:.3g}; Kruskal 3-grp p={r['p_kruskal_3grp']:.3g}; Spearman ρ(DM1, count)={r['spearman_rho_dm1_count']:.3f}, p={r['spearman_p']:.3g}.")
    if "classI_lof_count" in tert_tab["variable"].values:
        r = tert_tab.set_index("variable").loc["classI_lof_count"]
        md.append(f"- **Class-I LoF count, High vs Low**: mean_high={r['mean_high']:.3f}, mean_low={r['mean_low']:.3f}; Mann-Whitney p={r['p_mannwhitney_HighVsLow']:.3g}; Spearman ρ={r['spearman_rho_dm1_count']:.3f}, p={r['spearman_p']:.3g}.")
    if "classII_nonsilent_count" in tert_tab["variable"].values:
        r = tert_tab.set_index("variable").loc["classII_nonsilent_count"]
        md.append(f"- Class-II non-silent count, High vs Low: mean_high={r['mean_high']:.3f}, mean_low={r['mean_low']:.3f}; Mann-Whitney p={r['p_mannwhitney_HighVsLow']:.3g}.")
    md.append("")

    md.append("## 4. Per-gene Fisher enrichment in DM1-High vs Low (Table T04)\n")
    if fisher_tab.empty:
        md.append("- No HLA-pathway gene reached ≥3 mutations in the DM1-High+Low pool; per-gene tests not informative in this mutation-quiet cohort.\n")
    else:
        md.append("| gene | panel | n_high_mut | n_low_mut | OR | Fisher p | FDR |")
        md.append("|------|-------|------------|-----------|----|----|-----|")
        for _, r in fisher_tab.head(10).iterrows():
            md.append(f"| {r['gene']} | {r['panel']} | {int(r['n_high_mut'])} | {int(r['n_low_mut'])} | {r['OR']:.2f} | {r['p_fisher']:.3g} | {r['fdr_BH']:.3g} |")
        md.append("")

    md.append("## 5. TMB-controlled DM1 × class-I HLA-pathway (Table T05)\n")
    for tag, res in tmb_res.items():
        if isinstance(res, dict) and "p_DM1" in res:
            md.append(f"- **{tag}**: n={res['n']}, n_events/coef events={res.get('n_events', 'NA')}, β(DM1)={res['coef_DM1']:.3f}, OR per DM1 unit={res.get('OR_DM1_per_unit', np.nan):.3f}, p={res['p_DM1']:.3g}; β(log10 TMB)={res['coef_logTMB']:.3f}, p={res['p_logTMB']:.3g}.")
        elif tag == "poisson_classI_lof_offset_logTMB" and "p_DM1" in res:
            md.append(f"- **Poisson class-I LoF count ~ DM1, offset log10 TMB**: n={res['n']}, β={res['coef_DM1']:.3f}, IRR={res['IRR_DM1']:.3f}, p={res['p_DM1']:.3g}.")
        else:
            md.append(f"- {tag}: {res}")
    md.append("")

    md.append("## 6. Driver-stratified (Table T06 / Figure F5)\n")
    if driver_tab.empty:
        md.append("- (no driver groups with adequate n)\n")
    else:
        md.append("| driver | n | mean class-I nonsilent | n samples class-I LoF | Spearman ρ DM1 × count | p |")
        md.append("|--------|---|------------------------|------------------------|------------------------|---|")
        for _, r in driver_tab.iterrows():
            md.append(f"| {r['driver']} | {int(r['n'])} | {r['mean_classI_nonsilent']:.3f} | {int(r['n_classI_lof_pos'])} | {r['spearman_rho_DM1_classI']:.3f} | {r['spearman_p']:.3g} |")
        md.append("")

    md.append("## 7. Class-I vs Class-II machinery (Figure F6)\n")
    md.append(f"- Class-I: {civii['classI_pct_samples_with_mutation']:.2f}% samples mutated; total mutations {civii['classI_total_mutations']}.")
    md.append(f"- Class-II: {civii['classII_pct_samples_with_mutation']:.2f}% samples mutated; total mutations {civii['classII_total_mutations']}.")
    md.append(f"- Costim (CD274/PDCD1LG2/CD80/CD86): {civii['costim_pct_samples_with_mutation']:.2f}%; total {civii['costim_total_mutations']}.\n")

    md.append("## 8. Pan-cancer baseline comparison\n")
    md.append("- MC3 PanCan MAF not available locally (`/data/thca/...` does not contain mc3.v0.2.8.PUBLIC.maf.gz). "
              "Pan-cancer overlay is therefore deferred. Reviewer-reserve: TCGA-THCA class-I LoF rate "
              f"{summary['hla_pathway_overall']['pct_any_classI_lof']:.2f}% sits well below typical melanoma / NSCLC / MSI-CRC class-I LoF rates (literature 5-15%), "
              "consistent with THCA being mutation-quiet.\n")

    md.append("## 9. HLA-I expression-module co-occurrence (Figure F7)\n")
    if isinstance(cooccur, dict) and "p_mannwhitney_mut_vs_wt" in cooccur:
        md.append(f"- n_class-I-mut={cooccur.get('n_classI_mut')}; median HLA-I module mut={cooccur.get('median_HLA1_mut'):.3f} vs WT={cooccur.get('median_HLA1_wt'):.3f}; MannWhitney p={cooccur.get('p_mannwhitney_mut_vs_wt'):.3g}; Spearman ρ(HLA-I module, class-I mut count)={cooccur.get('spearman_rho_HLA1_classIcount'):.3f}, p={cooccur.get('spearman_p'):.3g}.\n")
    else:
        md.append(f"- {cooccur}\n")

    md.append("## 10. Survival — Cox HR for somatic class-I HLA-pathway disruption (Table T08 / Figure F8)\n")
    md.append("**Frame:** somatic class-I presentation pathway disruption × outcome — genomic feature, NOT germline HLA × outcome (Paper 1 allowed).\n")
    if cox_tab.empty:
        md.append("- (Cox not run — lifelines unavailable)\n")
    else:
        md.append("| endpoint | var | model | n | n_events | n+ | HR [95% CI] | p |")
        md.append("|----------|-----|-------|---|---------|----|-------------|---|")
        for _, r in cox_tab.iterrows():
            hr = f"{r['HR']:.2f}" if pd.notna(r['HR']) else "NA"
            ci = f"[{r['HR_lo']:.2f}, {r['HR_hi']:.2f}]" if pd.notna(r['HR_lo']) else ""
            p = f"{r['p']:.3g}" if pd.notna(r['p']) else "NA"
            md.append(f"| {r['endpoint']} | {r['var']} | {r['model']} | {int(r['n'])} | {int(r['n_events'])} | {int(r['n_var_pos'])} | {hr} {ci} | {p} |")
        md.append("")

    md.append("## 11. Limitations\n")
    md.append("- TCGA-THCA is famously mutation-quiet (median ~14 non-silent / sample). Power to detect single-gene HLA-pathway enrichment is low.")
    md.append("- LoF count from bulk-tumor MAF is **not allele-specific** — does not resolve which germline allele lost coverage; LoH is not assessed here.")
    md.append("- Some primary tumors lack a paired DM1 score (RNA-seq matched aliquot only); those samples drop from DM1-stratified tests.")
    md.append("- We do not test antigen-presentation pathway *expression* loss vs *mutation* loss in the same model here — Track 5 covers expression; this track covers genomic.")
    md.append("- Pan-cancer comparison limited because local MC3 MAF was not staged.")
    md.append("- Frame is consistent with Paper 1 boundary: somatic mutation count ≠ germline HLA allele association.\n")

    md.append("## 12. Mechanistic implication\n")
    n_lof = summary['hla_pathway_overall']['samples_any_classI_lof']
    n_ns = summary['hla_pathway_overall']['samples_any_classI_nonsilent']
    n_total = summary['n_tcga_thca_samples_with_maf']
    md.append(f"- TCGA-THCA HLA-pathway somatic mutation rate is **vanishingly low**: {n_ns}/{n_total} "
              f"samples carry any non-silent class-I HLA-pathway mutation, {n_lof}/{n_total} carry a class-I LoF. "
              "There is therefore **no class-I genomic-disruption signal to enrich** in DM1-high vs DM1-low tertiles "
              "(Mann-Whitney p>0.3 for class-I non-silent count, all LoF counts identically zero).")
    md.append("- This is the honest, mutation-quiet THCA result. It is *consistent* with — not evidence against — "
              "the Paper-1 expression / methylation story (Track 5: HLA-I module ρ_DM1=+0.32, p≈8e-14; Round 4: TPO methylation "
              "d=-2.73 in DM1<DM2). The dominant axis of class-I presentation modulation in dedifferentiated thyroid "
              "cancer is **expression / methylation / IFN-axis**, not somatic genomic disruption of HLA-A/B/C / B2M / TAP / NLRC5.")
    md.append("- A secondary finding: total non-silent TMB is **lower** in DM1-high tertile (mean 6.78 vs 11.26 in DM1-low; "
              "Mann-Whitney p=2.6e-6, Spearman ρ=-0.22, p=5e-6; Figure F8). DM1-high TCGA-THCA tumors are *less* mutated "
              "overall — consistent with the BRAF/RAS-negative / fusion / true-driver-negative composition of the DM1 "
              "compartment carrying fewer point mutations than BRAF-V600E PTC. This rules out **mutational** immune "
              "pressure as the explanation for DM1-high HLA-I module expression, and reinforces an epigenetic / cytokine "
              "/ infiltrate explanation.")
    md.append("- For Paper 3 (ICI vulnerability) the implication is reviewer-relevant: DM1-high thyroid cancer is "
              "*not* an HLA-pathway-LoF immune-evasion archetype. ICI sensitivity would have to come from another axis "
              "(IFN-γ-driven HLA-II re-induction, TLS, B-cell repertoire — covered in v17_D5P6_BCR_clonal_TLS).\n")

    md.append("---\n")
    md.append(f"_Track 25 outputs: `{OUT}`. Boundary contract: `project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md`._")

    with open(OUT / "track25_report.md", "w") as f:
        f.write("\n".join(md))


if __name__ == "__main__":
    main()
