#!/usr/bin/env python3
"""Track 6 — Pan-cancer HLA module x DM1 score deep dive.

Boundary: HLA-I/II as transcriptomic gene-expression module signatures only.
NO allele genotype from cancer data. Caption every figure accordingly.

Pipeline:
  1. Stream HLA-I + HLA-II module genes from TCGA pancan expression
  2. Build per-sample HLA-I and HLA-II module scores (within-lineage z-mean)
  3. DM1 x HLA module Spearman per lineage (forest, FDR/Bonferroni)
  4. Sign-coherence binomial test across 33 lineages
  5. Lineage-stratified scatter grid
  6. HLA module x survival (Cox per lineage, OS + DSS)
  7. Three-way Cox: DM1 + HLA module simultaneously
  8. DepMap connection: MYC/NAMPT-high vs HLA-I expression in cell lines
  9. Hallmark cross-talk: DM1 vs HLA-I vs Allograft Rejection triangle
 10. DM1-high x HLA-I-low quadrant per lineage table
 11. Methylation x HLA-I (documented limitation - HM450 fetch blocked per memory)

Outputs:
  /home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track6_pancan_hla_dm1/
"""
from __future__ import annotations

import gzip
import json
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
DATA_PANCAN = ROOT / "project/data/raw/TCGA_pancan"
PANCAN_EXPR = DATA_PANCAN / "pancan_geneExp.gz"
SURV = DATA_PANCAN / "survival.tsv"
P11 = ROOT / "project/results/paper11_pancancer"
DM1_FILE = P11 / "phase_E_lineage_specific/lineage_portable_dm1_per_sample.tsv"
HALLMARK_LONG = P11 / "phase_G_hallmark/hallmark_dm1_corr_long.tsv"
DEPMAP_RAW = Path("/data/thca/repo_results/p3_p9_full_execution/paper9/raw")
DEPMAP_EXPR_FILE = DEPMAP_RAW / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
CCLE_STATE = Path("/data/thca/repo_results/paper9_sl_first_pass/ccle_dm1_state.tsv")

OUT = ROOT / "project/results/hla_deepdive_2026_05_08/track6_pancan_hla_dm1"
FIG_DIR = OUT / "figures"
TBL_DIR = OUT / "tables"
GRID_DIR = OUT / "scatter_grid"
for d in (OUT, FIG_DIR, TBL_DIR, GRID_DIR):
    d.mkdir(parents=True, exist_ok=True)

CAPTION = "HLA gene-expression module - not allele genotype."

HLA1_CORE = [
    "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP", "NLRC5",
    "IRF1", "PSMB8", "PSMB9", "ERAP1", "ERAP2", "HLA-E", "HLA-F", "HLA-G",
    "CALR", "CANX", "PDIA3",
]
HLA2_CORE = [
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
    "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "CIITA", "CD74",
]

MIN_N = 30  # lineage filter


def stream_pancan_expr(path: Path, target_genes: list[str]) -> pd.DataFrame:
    """Stream pancan expression matrix, return DataFrame samples x genes."""
    print(f"[stream] {path} for {len(target_genes)} genes ...", flush=True)
    target = set(target_genes)
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        sample_ids = header[1:]
        out = {}
        n = 0
        for line in f:
            n += 1
            if n % 5000 == 0:
                print(f"  ... scanned {n} rows, matched {len(out)}", flush=True)
            tab = line.find("\t")
            if tab < 0:
                continue
            gene = line[:tab].strip()
            if gene in target:
                vals = line[tab + 1 :].rstrip("\n").split("\t")
                out[gene] = pd.to_numeric(pd.Series(vals), errors="coerce").values
                if len(out) == len(target):
                    break
    print(f"[stream] matched {len(out)}/{len(target)} genes")
    df = pd.DataFrame(out, index=sample_ids)
    df.index.name = "sample"
    return df


def within_lineage_zscore(df: pd.DataFrame, lineage_col: str, gene_cols: list[str]) -> pd.DataFrame:
    """Z-score each gene within each lineage."""
    out = df.copy()
    for lin, idx in df.groupby(lineage_col).groups.items():
        sub = df.loc[idx, gene_cols]
        mu = sub.mean(axis=0)
        sd = sub.std(axis=0).replace(0, np.nan)
        out.loc[idx, gene_cols] = (sub - mu) / sd
    return out


# ------------------------------------------------------------
# Step 1+2: load DM1, expression, build HLA-I and HLA-II module scores
# ------------------------------------------------------------
def step_load_and_score():
    print("\n=== Step 1+2: Load + HLA module scoring ===")
    dm1 = pd.read_csv(DM1_FILE, sep="\t")
    print(f"DM1 portable samples: {len(dm1)}, lineages: {dm1['lineage'].nunique()}")

    target = HLA1_CORE + HLA2_CORE
    expr = stream_pancan_expr(PANCAN_EXPR, target)
    expr_long = expr.reset_index()

    df = dm1.merge(expr_long, on="sample", how="inner")
    print(f"Merged: {len(df)} samples")

    # Filter lineages with N >= 30
    counts = df["lineage"].value_counts()
    keep = counts[counts >= MIN_N].index.tolist()
    df = df[df["lineage"].isin(keep)].reset_index(drop=True)
    print(f"After N>={MIN_N} filter: {len(df)} samples, {len(keep)} lineages")

    hla1_present = [g for g in HLA1_CORE if g in df.columns]
    hla2_present = [g for g in HLA2_CORE if g in df.columns]
    print(f"HLA-I genes present: {len(hla1_present)}/{len(HLA1_CORE)}: {hla1_present}")
    print(f"HLA-II genes present: {len(hla2_present)}/{len(HLA2_CORE)}: {hla2_present}")

    # Within-lineage z-score then mean -> module score
    df = within_lineage_zscore(df, "lineage", hla1_present + hla2_present)
    df["HLA1_module"] = df[hla1_present].mean(axis=1)
    df["HLA2_module"] = df[hla2_present].mean(axis=1)

    # Save per-sample module scores
    out_cols = ["sample", "lineage", "DM1_portable", "HLA1_module", "HLA2_module"]
    df[out_cols].to_csv(TBL_DIR / "T01_pancan_hla_module_per_sample.tsv", sep="\t", index=False)
    df.to_csv(TBL_DIR / "T01b_pancan_full_module_genes_zscored.tsv", sep="\t", index=False)
    print(f"Saved T01: {len(df)} samples")
    return df, hla1_present, hla2_present


# ------------------------------------------------------------
# Step 3: DM1 x HLA module per-lineage Spearman forest
# ------------------------------------------------------------
def step_per_lineage_spearman(df: pd.DataFrame):
    print("\n=== Step 3: Per-lineage Spearman ===")

    def _do(module_col: str, label: str):
        rows = []
        for lin, sub in df.groupby("lineage"):
            x = sub["DM1_portable"].values
            y = sub[module_col].values
            mask = np.isfinite(x) & np.isfinite(y)
            if mask.sum() < 10:
                continue
            rho, p = stats.spearmanr(x[mask], y[mask])
            rows.append({"lineage": lin, "n": int(mask.sum()), "rho": rho, "p": p})
        out = pd.DataFrame(rows)
        out["fdr_bh"] = multipletests(out["p"].fillna(1), method="fdr_bh")[1]
        out["bonferroni"] = np.minimum(out["p"] * len(out), 1.0)
        out = out.sort_values("rho")
        out.to_csv(TBL_DIR / f"T02_dm1_{label}_per_lineage.tsv", sep="\t", index=False)
        print(f"  {label}: {len(out)} lineages; median rho = {out['rho'].median():.3f}")
        return out

    h1 = _do("HLA1_module", "hla1")
    h2 = _do("HLA2_module", "hla2")
    return h1, h2


def plot_forest(df: pd.DataFrame, out_stem: Path, title: str):
    df = df.sort_values("rho").reset_index(drop=True)
    n = len(df)
    fig, ax = plt.subplots(figsize=(8, max(5, n * 0.28)))
    y = np.arange(n)
    colors = ["#cc4c33" if (r > 0 and p < 0.05) else
              "#3577b3" if (r < 0 and p < 0.05) else "#888"
              for r, p in zip(df["rho"], df["p"])]
    ax.barh(y, df["rho"], color=colors, alpha=0.85)
    sig_marker = ["***" if f < 0.001 else "**" if f < 0.01 else "*" if f < 0.05 else ""
                  for f in df["fdr_bh"]]
    for i, (rho, m) in enumerate(zip(df["rho"], sig_marker)):
        if m:
            ax.text(rho + (0.01 if rho >= 0 else -0.01), i, m,
                    ha="left" if rho >= 0 else "right", va="center", fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{lin}  (n={n})" for lin, n in zip(df["lineage"], df["n"])],
                       fontsize=7)
    ax.axvline(0, color="black", lw=0.5)
    ax.set_xlabel("Spearman rho (DM1 vs HLA module)")
    ax.set_title(title + f"\n{CAPTION}", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_stem.with_suffix(".png"), dpi=180)
    fig.savefig(out_stem.with_suffix(".pdf"))
    plt.close(fig)


# ------------------------------------------------------------
# Step 4: sign coherence
# ------------------------------------------------------------
def step_sign_coherence(h1: pd.DataFrame, h2: pd.DataFrame):
    print("\n=== Step 4: Sign coherence ===")
    rows = []
    for label, df in [("HLA-I", h1), ("HLA-II", h2)]:
        n = len(df)
        npos = int((df["rho"] > 0).sum())
        nneg = int((df["rho"] < 0).sum())
        majority = max(npos, nneg)
        # Two-sided binomial
        p_two = stats.binomtest(majority, n, p=0.5).pvalue
        rows.append({"module": label, "n_lineages": n, "n_pos": npos, "n_neg": nneg,
                     "fraction_majority": majority / n,
                     "majority_sign": "+" if npos >= nneg else "-",
                     "binom_p_two_sided": p_two})
        print(f"  {label}: n={n}, +={npos}, -={nneg}, majority={majority}, p={p_two:.2e}")
    out = pd.DataFrame(rows)
    out.to_csv(TBL_DIR / "T03_sign_coherence.tsv", sep="\t", index=False)
    return out


# ------------------------------------------------------------
# Step 5: scatter grid
# ------------------------------------------------------------
def step_scatter_grid(df: pd.DataFrame, h1: pd.DataFrame, module_col: str, label: str):
    print(f"\n=== Step 5: Scatter grid {label} ===")
    lineages = sorted(df["lineage"].unique())
    n = len(lineages)
    cols = 6
    rows = (n + cols - 1) // cols
    fig = plt.figure(figsize=(cols * 2.6, rows * 2.5))
    gs = GridSpec(rows, cols, figure=fig, hspace=0.55, wspace=0.45)
    for i, lin in enumerate(lineages):
        ax = fig.add_subplot(gs[i // cols, i % cols])
        sub = df[df["lineage"] == lin]
        x = sub["DM1_portable"].values
        y = sub[module_col].values
        mask = np.isfinite(x) & np.isfinite(y)
        ax.scatter(x[mask], y[mask], s=4, alpha=0.4, c="#445", edgecolors="none")
        # Linear fit
        if mask.sum() >= 5:
            try:
                m, b = np.polyfit(x[mask], y[mask], 1)
                xs = np.linspace(x[mask].min(), x[mask].max(), 30)
                ax.plot(xs, m * xs + b, color="#cc4c33", lw=1)
            except Exception:
                pass
        rho_row = h1[h1["lineage"] == lin]
        if len(rho_row):
            r = rho_row["rho"].iloc[0]
            p = rho_row["p"].iloc[0]
            ax.set_title(f"{lin[:28]}\nrho={r:+.2f} p={p:.1e}", fontsize=7)
        else:
            ax.set_title(lin[:30], fontsize=7)
        ax.tick_params(labelsize=6)
        if i // cols == rows - 1:
            ax.set_xlabel("DM1", fontsize=7)
        if i % cols == 0:
            ax.set_ylabel(label, fontsize=7)
    fig.suptitle(f"DM1 vs {label} per lineage  -  {CAPTION}", fontsize=11, y=0.998)
    fig.savefig(GRID_DIR / f"F_grid_{label.replace('-', '')}.png", dpi=160, bbox_inches="tight")
    fig.savefig(GRID_DIR / f"F_grid_{label.replace('-', '')}.pdf", bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------
# Step 6: HLA module x survival (Cox per lineage)
# ------------------------------------------------------------
def parse_stage(s):
    if pd.isna(s) or not isinstance(s, str):
        return np.nan
    s = s.lower().strip()
    if "stage iv" in s:
        return 4
    if "stage iii" in s:
        return 3
    if "stage ii" in s:
        return 2
    if "stage i" in s:
        return 1
    return np.nan


def step_survival(df: pd.DataFrame):
    print("\n=== Step 6: HLA module x survival Cox ===")
    try:
        from lifelines import CoxPHFitter
    except ImportError:
        print("[step6] lifelines not installed; skipping")
        return None
    surv = pd.read_csv(SURV, sep="\t", low_memory=False)
    j = df.merge(surv, on="sample", how="left")
    j["age"] = pd.to_numeric(j.get("age_at_initial_pathologic_diagnosis"), errors="coerce")
    j["stage_int"] = j.get("ajcc_pathologic_tumor_stage").apply(parse_stage)

    rows = []
    for lin, sub in j.groupby("lineage"):
        for module_col, mod_label in [("HLA1_module", "HLA-I"), ("HLA2_module", "HLA-II")]:
            for outcome, t_col, e_col in [("OS", "OS.time", "OS"), ("DSS", "DSS.time", "DSS")]:
                cox = sub[[outcome, t_col, module_col, "age", "stage_int"]].copy()
                cox = cox.dropna(subset=[outcome, t_col, module_col])
                cox = cox[cox[t_col] > 0]
                if len(cox) < 25 or cox[outcome].sum() < 5:
                    continue
                covars = [module_col]
                if cox["age"].notna().sum() / len(cox) > 0.5:
                    covars.append("age")
                if cox["stage_int"].notna().sum() / len(cox) > 0.3:
                    covars.append("stage_int")
                cox = cox[[outcome, t_col] + covars].dropna()
                if len(cox) < 25:
                    continue
                cph = CoxPHFitter(penalizer=0.001)
                try:
                    cph.fit(cox, duration_col=t_col, event_col=outcome)
                    s = cph.summary.loc[module_col]
                    rows.append({
                        "lineage": lin, "module": mod_label, "outcome": outcome,
                        "n": len(cox), "events": int(cox[outcome].sum()),
                        "covariates": "+".join(covars),
                        "HR": float(s["exp(coef)"]),
                        "HR_lower95": float(s["exp(coef) lower 95%"]),
                        "HR_upper95": float(s["exp(coef) upper 95%"]),
                        "p": float(s["p"]),
                        "concordance": float(cph.concordance_index_),
                    })
                except Exception as e:
                    print(f"   skip {lin}/{mod_label}/{outcome}: {e}")
    out = pd.DataFrame(rows)
    if len(out):
        # FDR within (module x outcome) groups
        out["fdr_bh"] = np.nan
        for (mod, oc), g in out.groupby(["module", "outcome"]):
            out.loc[g.index, "fdr_bh"] = multipletests(g["p"].fillna(1), method="fdr_bh")[1]
    out.to_csv(TBL_DIR / "T04_hla_survival_cox.tsv", sep="\t", index=False)
    print(f"  Cox rows: {len(out)}")
    if len(out):
        sig = out[(out["fdr_bh"] < 0.1)]
        print(f"  FDR<0.1 hits: {len(sig)}")
    return out


def plot_survival_forest(cox_df: pd.DataFrame, module: str, outcome: str):
    if cox_df is None or not len(cox_df):
        return
    sub = cox_df[(cox_df["module"] == module) & (cox_df["outcome"] == outcome)].copy()
    if not len(sub):
        return
    sub = sub.sort_values("HR")
    n = len(sub)
    fig, ax = plt.subplots(figsize=(8, max(5, n * 0.28)))
    y = np.arange(n)
    log_hr = np.log(sub["HR"].values)
    log_lo = np.log(sub["HR_lower95"].clip(lower=1e-6).values)
    log_hi = np.log(sub["HR_upper95"].clip(upper=1e6).values)
    err = np.array([log_hr - log_lo, log_hi - log_hr])
    colors = ["#cc4c33" if (h > 1 and p < 0.05) else
              "#3577b3" if (h < 1 and p < 0.05) else "#888"
              for h, p in zip(sub["HR"], sub["p"])]
    ax.errorbar(log_hr, y, xerr=err, fmt="o", color="black", ecolor="#888",
                markerfacecolor="white", markersize=4, lw=0.8, capsize=2)
    for i, (lh, c) in enumerate(zip(log_hr, colors)):
        ax.plot(lh, i, "o", color=c, markersize=5, markeredgecolor="black", markeredgewidth=0.4)
    ax.axvline(0, color="black", lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{lin}  (n={n2}, e={e})" for lin, n2, e in
                        zip(sub["lineage"], sub["n"], sub["events"])], fontsize=7)
    ax.set_xlabel("log HR (per 1 unit module score)")
    ax.set_title(f"{module} module x {outcome} - Cox forest\n{CAPTION}", fontsize=10)
    fig.tight_layout()
    stem = FIG_DIR / f"F06_survival_forest_{module.replace('-','')}_{outcome}"
    fig.savefig(stem.with_suffix(".png"), dpi=180)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


# ------------------------------------------------------------
# Step 7: Three-way Cox: DM1 + HLA module simultaneously
# ------------------------------------------------------------
def step_three_way(df: pd.DataFrame):
    print("\n=== Step 7: Three-way DM1 + HLA Cox ===")
    try:
        from lifelines import CoxPHFitter
    except ImportError:
        return None
    surv = pd.read_csv(SURV, sep="\t", low_memory=False)
    j = df.merge(surv, on="sample", how="left")
    j["age"] = pd.to_numeric(j.get("age_at_initial_pathologic_diagnosis"), errors="coerce")
    j["stage_int"] = j.get("ajcc_pathologic_tumor_stage").apply(parse_stage)

    rows = []
    for lin, sub in j.groupby("lineage"):
        for module_col, mod_label in [("HLA1_module", "HLA-I"), ("HLA2_module", "HLA-II")]:
            for outcome, t_col in [("OS", "OS.time"), ("DSS", "DSS.time")]:
                cox = sub[[outcome, t_col, "DM1_portable", module_col, "age", "stage_int"]].copy()
                cox = cox.dropna(subset=[outcome, t_col, "DM1_portable", module_col])
                cox = cox[cox[t_col] > 0]
                if len(cox) < 30 or cox[outcome].sum() < 5:
                    continue
                covars = ["DM1_portable", module_col]
                if cox["age"].notna().sum() / len(cox) > 0.5:
                    covars.append("age")
                if cox["stage_int"].notna().sum() / len(cox) > 0.3:
                    covars.append("stage_int")
                cox = cox[[outcome, t_col] + covars].dropna()
                if len(cox) < 30:
                    continue
                cph = CoxPHFitter(penalizer=0.001)
                # Also fit DM1-only baseline
                base = cox[[outcome, t_col, "DM1_portable"] + [c for c in covars if c not in {"DM1_portable", module_col}]].copy()
                try:
                    cph.fit(cox, duration_col=t_col, event_col=outcome)
                    s_dm1 = cph.summary.loc["DM1_portable"]
                    s_mod = cph.summary.loc[module_col]
                    cph_b = CoxPHFitter(penalizer=0.001)
                    cph_b.fit(base, duration_col=t_col, event_col=outcome)
                    s_dm1_alone = cph_b.summary.loc["DM1_portable"]
                    rows.append({
                        "lineage": lin, "module": mod_label, "outcome": outcome,
                        "n": len(cox), "events": int(cox[outcome].sum()),
                        "DM1_alone_HR": float(s_dm1_alone["exp(coef)"]),
                        "DM1_alone_p": float(s_dm1_alone["p"]),
                        "DM1_adj_HR": float(s_dm1["exp(coef)"]),
                        "DM1_adj_p": float(s_dm1["p"]),
                        "module_adj_HR": float(s_mod["exp(coef)"]),
                        "module_adj_p": float(s_mod["p"]),
                        "DM1_HR_attenuation": float(s_dm1["exp(coef)"]) / float(s_dm1_alone["exp(coef)"]),
                    })
                except Exception:
                    continue
    out = pd.DataFrame(rows)
    out.to_csv(TBL_DIR / "T05_three_way_dm1_hla_cox.tsv", sep="\t", index=False)
    print(f"  Three-way Cox rows: {len(out)}")
    return out


def plot_three_way_attenuation(threeway: pd.DataFrame):
    if threeway is None or not len(threeway):
        return
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    for ax, mod in zip(axes, ["HLA-I", "HLA-II"]):
        sub = threeway[(threeway["module"] == mod) & (threeway["outcome"] == "OS")].copy()
        if not len(sub):
            ax.set_title(f"{mod} (no data)")
            continue
        sub = sub.sort_values("DM1_HR_attenuation")
        y = np.arange(len(sub))
        ax.barh(y, sub["DM1_HR_attenuation"] - 1, color="#3577b3", alpha=0.7)
        ax.axvline(0, color="black", lw=0.5)
        ax.set_yticks(y)
        ax.set_yticklabels([s[:28] for s in sub["lineage"]], fontsize=7)
        ax.set_xlabel("DM1 HR attenuation when adding HLA module\n(>0 = HLA strengthens DM1; <0 = HLA absorbs DM1)")
        ax.set_title(f"{mod} OS")
    fig.suptitle(f"DM1 prognostic effect change after adding HLA module - {CAPTION}", fontsize=10)
    fig.tight_layout()
    stem = FIG_DIR / "F07_dm1_hr_attenuation"
    fig.savefig(stem.with_suffix(".png"), dpi=180)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


# ------------------------------------------------------------
# Step 8: DepMap connection
# ------------------------------------------------------------
def step_depmap_hla(hla1_present: list[str], hla2_present: list[str]):
    print("\n=== Step 8: DepMap HLA-I expression vs MYC/NAMPT ===")
    if not DEPMAP_EXPR_FILE.exists():
        print(f"  DepMap expr missing: {DEPMAP_EXPR_FILE}")
        return None
    if not CCLE_STATE.exists():
        print(f"  CCLE DM1 state missing: {CCLE_STATE}")
        return None

    state = pd.read_csv(CCLE_STATE, sep="\t")
    state = state[["line", "DM1_like_score", "DM1_high"]].copy()
    state.columns = ["ModelID_alias", "DM1_like_score", "DM1_high"]

    # Read DepMap expression header
    print("  Reading DepMap expression header...")
    expr_header = pd.read_csv(DEPMAP_EXPR_FILE, nrows=0).columns.tolist()
    import re
    sym2col = {}
    target_genes = hla1_present + hla2_present + ["MYC", "NAMPT"]
    for c in expr_header:
        m = re.match(r"^([A-Z0-9\-]+)\s*\(\d+\)$", c)
        if m and m.group(1) in target_genes:
            sym2col[m.group(1)] = c
    first_col = expr_header[0]
    cols = [first_col] + list(sym2col.values())
    expr = pd.read_csv(DEPMAP_EXPR_FILE, usecols=cols)
    expr = expr.rename(columns={first_col: "ModelID"})
    expr = expr.rename(columns={v: k for k, v in sym2col.items()})
    print(f"  DepMap expr loaded: {expr.shape}, matched {len(sym2col)} target genes")

    # Match CCLE state line to DepMap ModelID via Model.csv (uses StrippedCellLineName)
    model = pd.read_csv(DEPMAP_RAW / "Model.csv", usecols=["ModelID", "StrippedCellLineName"])
    state["StrippedCellLineName"] = state["ModelID_alias"].str.split("_").str[0].str.upper()
    state2 = state.merge(model, on="StrippedCellLineName", how="inner")
    df = expr.merge(state2[["ModelID", "DM1_like_score", "DM1_high"]], on="ModelID", how="inner")
    print(f"  Cell lines with both DM1 score + expression: {len(df)}")

    hla1_in_dm = [g for g in hla1_present if g in df.columns]
    hla2_in_dm = [g for g in hla2_present if g in df.columns]
    if not len(hla1_in_dm):
        print("  No HLA-I genes in DepMap; bailing")
        return None
    # Z-score globally and average
    df["HLA1_dm_module"] = df[hla1_in_dm].apply(stats.zscore, nan_policy="omit").mean(axis=1)
    df["HLA2_dm_module"] = df[hla2_in_dm].apply(stats.zscore, nan_policy="omit").mean(axis=1) if hla2_in_dm else np.nan

    # Test: MYC-high lines have low HLA-I?
    rows = []
    for gene in ["MYC", "NAMPT"]:
        if gene not in df.columns:
            continue
        hi = df[df[gene] >= df[gene].median()]
        lo = df[df[gene] < df[gene].median()]
        for mod_label, mod_col in [("HLA-I", "HLA1_dm_module"), ("HLA-II", "HLA2_dm_module")]:
            if mod_col not in df.columns:
                continue
            x = hi[mod_col].dropna().values
            y = lo[mod_col].dropna().values
            if len(x) < 5 or len(y) < 5:
                continue
            t, p = stats.ttest_ind(x, y, equal_var=False)
            d = (np.mean(x) - np.mean(y)) / np.sqrt((np.var(x, ddof=1) + np.var(y, ddof=1)) / 2)
            rho_g, p_g = stats.spearmanr(df[gene], df[mod_col], nan_policy="omit")
            rows.append({"gene": gene, "module": mod_label,
                         "n_hi": len(x), "n_lo": len(y),
                         "mean_hi_module": float(np.mean(x)),
                         "mean_lo_module": float(np.mean(y)),
                         "cohens_d": float(d), "ttest_p": float(p),
                         "spearman_rho": float(rho_g), "spearman_p": float(p_g)})
    # Also DM1 vs HLA module
    rho_dm1_h1, p_dm1_h1 = stats.spearmanr(df["DM1_like_score"], df["HLA1_dm_module"], nan_policy="omit")
    rho_dm1_h2, p_dm1_h2 = stats.spearmanr(df["DM1_like_score"], df["HLA2_dm_module"], nan_policy="omit") \
        if "HLA2_dm_module" in df.columns else (np.nan, np.nan)
    rows.append({"gene": "DM1_score", "module": "HLA-I",
                 "n_hi": np.nan, "n_lo": np.nan,
                 "mean_hi_module": np.nan, "mean_lo_module": np.nan,
                 "cohens_d": np.nan, "ttest_p": np.nan,
                 "spearman_rho": float(rho_dm1_h1), "spearman_p": float(p_dm1_h1)})
    rows.append({"gene": "DM1_score", "module": "HLA-II",
                 "n_hi": np.nan, "n_lo": np.nan,
                 "mean_hi_module": np.nan, "mean_lo_module": np.nan,
                 "cohens_d": np.nan, "ttest_p": np.nan,
                 "spearman_rho": float(rho_dm1_h2), "spearman_p": float(p_dm1_h2)})
    out = pd.DataFrame(rows)
    out.to_csv(TBL_DIR / "T06_depmap_myc_nampt_vs_hla.tsv", sep="\t", index=False)

    # Plot scatter: MYC vs HLA-I in DepMap, colored by DM1
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    for ax, gene in zip(axes, ["MYC", "NAMPT", "DM1_like_score"]):
        if gene not in df.columns:
            continue
        ax.scatter(df[gene], df["HLA1_dm_module"], s=10, c=df["DM1_like_score"],
                   cmap="coolwarm", alpha=0.6, edgecolors="none")
        ax.set_xlabel(gene)
        ax.set_ylabel("HLA-I module (z-mean)")
        rho_v, p_v = stats.spearmanr(df[gene], df["HLA1_dm_module"], nan_policy="omit")
        ax.set_title(f"{gene} vs HLA-I  rho={rho_v:+.2f} p={p_v:.1e}", fontsize=9)
    fig.suptitle(f"DepMap cell lines: DM1-axis genes vs HLA-I module - {CAPTION}", fontsize=10)
    fig.tight_layout()
    stem = FIG_DIR / "F08_depmap_myc_nampt_dm1_vs_hla1"
    fig.savefig(stem.with_suffix(".png"), dpi=180)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)
    return out


# ------------------------------------------------------------
# Step 9: Hallmark cross-talk
# ------------------------------------------------------------
def step_hallmark_triangle(per_lineage_h1: pd.DataFrame, per_lineage_h2: pd.DataFrame):
    print("\n=== Step 9: Hallmark triangle ===")
    if not HALLMARK_LONG.exists():
        print("  Hallmark file missing")
        return None
    hall = pd.read_csv(HALLMARK_LONG, sep="\t")
    print(f"  Hallmark long table: {hall.shape}; cols: {hall.columns.tolist()}")
    # Filter to allograft rejection + a few canonical immune sets
    targets = ["Allograft Rejection", "Inflammatory Response",
               "Interferon Gamma Response", "Interferon Alpha Response",
               "TNF-alpha Signaling via NF-kB", "Complement",
               "TGF-beta Signaling", "Apoptosis"]
    hall["match"] = hall["hallmark"].astype(str).apply(
        lambda x: next((t for t in targets if t.lower() in x.lower()), None))
    hall_imm = hall[hall["match"].notna()].copy()

    # Triangle: per lineage, DM1-HLA1 rho vs DM1-Allograft rho
    # Pivot per lineage / hallmark
    if "spearman_r" in hall.columns:
        rcol = "spearman_r"
    elif "r" in hall.columns:
        rcol = "r"
    else:
        rcol = [c for c in hall.columns if "spearman" in c.lower() or c == "r"][0]

    pivot = hall_imm.pivot_table(index="lineage", columns="match", values=rcol, aggfunc="mean")
    # Merge with HLA-I per-lineage
    h1m = per_lineage_h1[["lineage", "rho"]].rename(columns={"rho": "rho_DM1_HLA1"})
    h2m = per_lineage_h2[["lineage", "rho"]].rename(columns={"rho": "rho_DM1_HLA2"})
    triangle = pivot.merge(h1m, on="lineage").merge(h2m, on="lineage")
    triangle.to_csv(TBL_DIR / "T07_hallmark_triangle.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, mod_col, mod_label in zip(axes, ["rho_DM1_HLA1", "rho_DM1_HLA2"], ["HLA-I", "HLA-II"]):
        if "Allograft Rejection" in triangle.columns:
            x = triangle[mod_col].values
            y = triangle["Allograft Rejection"].values
            mask = np.isfinite(x) & np.isfinite(y)
            ax.scatter(x[mask], y[mask], s=30, alpha=0.7, c="#445")
            for _, row in triangle.iterrows():
                if np.isfinite(row[mod_col]) and np.isfinite(row["Allograft Rejection"]):
                    ax.annotate(row["lineage"][:14], (row[mod_col], row["Allograft Rejection"]),
                                fontsize=5, alpha=0.6)
            if mask.sum() >= 5:
                rho, p = stats.spearmanr(x[mask], y[mask])
                ax.set_title(f"DM1-{mod_label} rho vs DM1-AllograftRej rho\nlineage Spearman={rho:+.2f} p={p:.2e}",
                             fontsize=9)
            ax.set_xlabel(f"DM1-{mod_label} rho per lineage")
            ax.set_ylabel("DM1-Allograft Rejection rho per lineage")
            ax.axvline(0, color="black", lw=0.4)
            ax.axhline(0, color="black", lw=0.4)
    fig.suptitle(f"Triangle: DM1 - HLA module - Allograft Rejection (Hallmark) - {CAPTION}", fontsize=10)
    fig.tight_layout()
    stem = FIG_DIR / "F09_triangle_dm1_hla_allograft"
    fig.savefig(stem.with_suffix(".png"), dpi=180)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)
    return triangle


# ------------------------------------------------------------
# Step 10: DM1-high + HLA-I-low quadrant
# ------------------------------------------------------------
def step_quadrant(df: pd.DataFrame):
    print("\n=== Step 10: DM1-high x HLA-I-low quadrant ===")
    rows = []
    rows_per_lin = []
    for lin, sub in df.groupby("lineage"):
        n = len(sub)
        # Per-lineage thresholds
        dm1_med = sub["DM1_portable"].median()
        h1_med = sub["HLA1_module"].median()
        h2_med = sub["HLA2_module"].median()
        q_h1cold = ((sub["DM1_portable"] > dm1_med) & (sub["HLA1_module"] < h1_med)).sum()
        q_h1hot = ((sub["DM1_portable"] > dm1_med) & (sub["HLA1_module"] >= h1_med)).sum()
        q_dm1lo_h1cold = ((sub["DM1_portable"] <= dm1_med) & (sub["HLA1_module"] < h1_med)).sum()
        q_dm1lo_h1hot = ((sub["DM1_portable"] <= dm1_med) & (sub["HLA1_module"] >= h1_med)).sum()
        q_h2cold = ((sub["DM1_portable"] > dm1_med) & (sub["HLA2_module"] < h2_med)).sum()
        rows_per_lin.append({"lineage": lin, "n": n,
                              "DM1hi_HLA1lo": int(q_h1cold),
                              "DM1hi_HLA1hi": int(q_h1hot),
                              "DM1lo_HLA1lo": int(q_dm1lo_h1cold),
                              "DM1lo_HLA1hi": int(q_dm1lo_h1hot),
                              "DM1hi_HLA2lo": int(q_h2cold),
                              "frac_DM1hi_HLA1lo": q_h1cold / n})
    quad = pd.DataFrame(rows_per_lin).sort_values("frac_DM1hi_HLA1lo", ascending=False)
    quad.to_csv(TBL_DIR / "T08_quadrant_per_lineage.tsv", sep="\t", index=False)

    # Heatmap of fractions
    fig, ax = plt.subplots(figsize=(6, max(5, len(quad) * 0.28)))
    ax.barh(np.arange(len(quad)), quad["frac_DM1hi_HLA1lo"], color="#3577b3", alpha=0.85)
    ax.set_yticks(np.arange(len(quad)))
    ax.set_yticklabels([f"{l}  (n={n})" for l, n in zip(quad["lineage"], quad["n"])], fontsize=7)
    ax.axvline(0.25, color="#cc4c33", lw=1, ls="--")
    ax.set_xlabel("Fraction in DM1-high & HLA-I-low quadrant")
    ax.set_title(f"Immune-cold DM1-high quadrant per lineage - {CAPTION}", fontsize=10)
    fig.tight_layout()
    stem = FIG_DIR / "F10_quadrant_dm1hi_hla1lo"
    fig.savefig(stem.with_suffix(".png"), dpi=180)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)
    return quad


# ------------------------------------------------------------
# Step 11: HLA-I module vs methylation proxy (epi_silencing_index)
# ------------------------------------------------------------
def step_methylation_proxy(df: pd.DataFrame):
    print("\n=== Step 11: HLA-I module vs epi_silencing_index (RNA proxy) ===")
    epi_file = P11 / "phase_A_epigenetic" / "epi_index_per_sample.tsv"
    if not epi_file.exists():
        print("  Epigenetic file missing; documenting limitation")
        return None
    epi = pd.read_csv(epi_file, sep="\t")
    j = df.merge(epi[["sample", "epi_silencing_index", "epi_writers", "epi_PRC2",
                       "epi_HDACs", "epi_demethylase"]], on="sample", how="inner")
    rows = []
    for lin, sub in j.groupby("lineage"):
        for ax_col in ["epi_silencing_index", "epi_writers", "epi_PRC2", "epi_HDACs"]:
            x = sub[ax_col].values
            y = sub["HLA1_module"].values
            mask = np.isfinite(x) & np.isfinite(y)
            if mask.sum() < 10:
                continue
            r, p = stats.spearmanr(x[mask], y[mask])
            rows.append({"lineage": lin, "axis": ax_col, "n": int(mask.sum()),
                         "rho": r, "p": p})
    out = pd.DataFrame(rows)
    out["fdr_bh"] = np.nan
    for ax, g in out.groupby("axis"):
        out.loc[g.index, "fdr_bh"] = multipletests(g["p"].fillna(1), method="fdr_bh")[1]
    out.to_csv(TBL_DIR / "T09_epi_silencing_vs_hla1.tsv", sep="\t", index=False)
    print(f"  Rows: {len(out)}")

    # Scatter: epi_silencing vs HLA1 module across all samples (pooled, color by lineage)
    fig, ax = plt.subplots(figsize=(8, 6))
    rho_pool, p_pool = stats.spearmanr(j["epi_silencing_index"], j["HLA1_module"], nan_policy="omit")
    ax.scatter(j["epi_silencing_index"], j["HLA1_module"], s=4, c=j["DM1_portable"],
               cmap="coolwarm", alpha=0.4, edgecolors="none")
    ax.set_xlabel("epi_silencing_index (RNA proxy of writer/PRC2/HDAC over demethylase)")
    ax.set_ylabel("HLA-I module (within-lineage z-mean)")
    ax.set_title(f"Pan-cancer HLA-I module vs epigenetic silencing proxy\npooled rho={rho_pool:+.3f} p={p_pool:.1e} - {CAPTION}", fontsize=10)
    fig.tight_layout()
    stem = FIG_DIR / "F11_epi_silencing_vs_hla1"
    fig.savefig(stem.with_suffix(".png"), dpi=180)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)
    return out


# ------------------------------------------------------------
# Master HLA1 vs HLA2 comparison (extra figure)
# ------------------------------------------------------------
def plot_module_comparison(h1: pd.DataFrame, h2: pd.DataFrame):
    merged = h1[["lineage", "rho"]].rename(columns={"rho": "rho_HLA1"}).merge(
        h2[["lineage", "rho"]].rename(columns={"rho": "rho_HLA2"}), on="lineage")
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(merged["rho_HLA1"], merged["rho_HLA2"], s=40, c="#445")
    for _, r in merged.iterrows():
        ax.annotate(r["lineage"][:18], (r["rho_HLA1"], r["rho_HLA2"]), fontsize=6, alpha=0.7)
    if len(merged) >= 5:
        rho, p = stats.spearmanr(merged["rho_HLA1"], merged["rho_HLA2"])
        ax.set_title(f"DM1 rho: HLA-I vs HLA-II per lineage\nlineage Spearman={rho:+.3f} p={p:.1e} - {CAPTION}",
                     fontsize=10)
    ax.axhline(0, color="black", lw=0.5)
    ax.axvline(0, color="black", lw=0.5)
    ax.set_xlabel("rho (DM1 vs HLA-I)")
    ax.set_ylabel("rho (DM1 vs HLA-II)")
    fig.tight_layout()
    stem = FIG_DIR / "F12_hla1_vs_hla2_per_lineage"
    fig.savefig(stem.with_suffix(".png"), dpi=180)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------
def main():
    print("=" * 60)
    print("Track 6 - Pan-cancer HLA module x DM1")
    print("=" * 60)

    df, hla1_present, hla2_present = step_load_and_score()

    h1, h2 = step_per_lineage_spearman(df)
    plot_forest(h1, FIG_DIR / "F02_dm1_hla1_forest", "DM1 vs HLA-I module per lineage")
    plot_forest(h2, FIG_DIR / "F03_dm1_hla2_forest", "DM1 vs HLA-II module per lineage")

    sign_df = step_sign_coherence(h1, h2)

    step_scatter_grid(df, h1, "HLA1_module", "HLA-I")
    step_scatter_grid(df, h2, "HLA2_module", "HLA-II")

    cox = step_survival(df)
    if cox is not None and len(cox):
        for mod in ["HLA-I", "HLA-II"]:
            for outc in ["OS", "DSS"]:
                plot_survival_forest(cox, mod, outc)

    threeway = step_three_way(df)
    plot_three_way_attenuation(threeway)

    depmap = step_depmap_hla(hla1_present, hla2_present)

    triangle = step_hallmark_triangle(h1, h2)

    quad = step_quadrant(df)

    epi_out = step_methylation_proxy(df)

    plot_module_comparison(h1, h2)

    # Master summary
    summary = {
        "n_samples_total": int(len(df)),
        "n_lineages_total": int(df["lineage"].nunique()),
        "hla1_genes_present": hla1_present,
        "hla2_genes_present": hla2_present,
        "hla1": {
            "median_rho": float(h1["rho"].median()),
            "n_lineages_pos": int((h1["rho"] > 0).sum()),
            "n_lineages_neg": int((h1["rho"] < 0).sum()),
            "n_FDR10_pos": int(((h1["rho"] > 0) & (h1["fdr_bh"] < 0.1)).sum()),
            "n_FDR10_neg": int(((h1["rho"] < 0) & (h1["fdr_bh"] < 0.1)).sum()),
        },
        "hla2": {
            "median_rho": float(h2["rho"].median()),
            "n_lineages_pos": int((h2["rho"] > 0).sum()),
            "n_lineages_neg": int((h2["rho"] < 0).sum()),
            "n_FDR10_pos": int(((h2["rho"] > 0) & (h2["fdr_bh"] < 0.1)).sum()),
            "n_FDR10_neg": int(((h2["rho"] < 0) & (h2["fdr_bh"] < 0.1)).sum()),
        },
        "sign_coherence": sign_df.to_dict(orient="records"),
        "top5_hla1_pos": h1.sort_values("rho", ascending=False).head(5).to_dict(orient="records"),
        "top5_hla1_neg": h1.sort_values("rho", ascending=True).head(5).to_dict(orient="records"),
        "top5_hla2_pos": h2.sort_values("rho", ascending=False).head(5).to_dict(orient="records"),
        "survival_FDR10": int((cox["fdr_bh"] < 0.1).sum()) if cox is not None and len(cox) else 0,
        "depmap_summary": depmap.to_dict(orient="records") if depmap is not None else None,
        "top_immune_cold_dm1high_lineages": quad.head(8).to_dict(orient="records"),
        "boundary": CAPTION,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"\n[done] summary -> {OUT / 'summary.json'}")
    return summary


if __name__ == "__main__":
    main()
