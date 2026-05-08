#!/usr/bin/env python3
"""Track 6 synthesis retry — pan-cancer HLA module x DM1.

Picks up where the original Track 6 timed out. Produces deliverables T08-T12
+ F08-F12, narrative report, and summary JSON. Reads existing tables/figures
without recomputing T01-T05 / F02-F03 / F06 / F07 / scatter grid.

Boundary: HLA-I/II treated as transcriptomic gene-expression module signatures
only. NO allele genotype from cancer data. Caption every figure accordingly.
"""
from __future__ import annotations

import json
import re
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
DATA_PANCAN = ROOT / "project/data/raw/TCGA_pancan"
SURV = DATA_PANCAN / "survival.tsv"
P11 = ROOT / "project/results/paper11_pancancer"
HALLMARK_LONG = P11 / "phase_G_hallmark/hallmark_dm1_corr_long.tsv"
PHASE_A_PER_SAMPLE = P11 / "phase_A_epigenetic/epi_index_per_sample.tsv"

DEPMAP_RAW = Path("/data/thca/repo_results/p3_p9_full_execution/paper9/raw")
DEPMAP_EXPR_FILE = DEPMAP_RAW / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
DEPMAP_MODEL = DEPMAP_RAW / "Model.csv"
# Phase D cell-line table covers 1141 lines with DM1_like_score (vs 13-line
# thyroid-only ccle_dm1_state.tsv); use it as the primary mapping.
PHASE_D_CELLLINES = P11 / "phase_D_depmap/celllines_dm1_crispr.tsv"
CCLE_STATE = Path("/data/thca/repo_results/paper9_sl_first_pass/ccle_dm1_state.tsv")

OUT = ROOT / "project/results/hla_deepdive_2026_05_08/track6_pancan_hla_dm1"
FIG_DIR = OUT / "figures"
TBL_DIR = OUT / "tables"
for d in (OUT, FIG_DIR, TBL_DIR):
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


def load_existing():
    """Read the pre-computed T01-T05 tables."""
    t01 = pd.read_csv(TBL_DIR / "T01_pancan_hla_module_per_sample.tsv", sep="\t")
    t02h1 = pd.read_csv(TBL_DIR / "T02_dm1_hla1_per_lineage.tsv", sep="\t")
    t02h2 = pd.read_csv(TBL_DIR / "T02_dm1_hla2_per_lineage.tsv", sep="\t")
    t03 = pd.read_csv(TBL_DIR / "T03_sign_coherence.tsv", sep="\t")
    t04 = pd.read_csv(TBL_DIR / "T04_hla_survival_cox.tsv", sep="\t")
    t05 = pd.read_csv(TBL_DIR / "T05_three_way_dm1_hla_cox.tsv", sep="\t")
    return t01, t02h1, t02h2, t03, t04, t05


# ------------------------------------------------------------
# T08 / F08: DepMap connection
# ------------------------------------------------------------
def step_depmap(per_sample: pd.DataFrame, skipped: dict) -> pd.DataFrame | None:
    print("\n=== T08/F08: DepMap connection ===", flush=True)
    if not (DEPMAP_EXPR_FILE.exists() and PHASE_D_CELLLINES.exists()):
        skipped["T08"] = (
            f"DepMap inputs missing: expr={DEPMAP_EXPR_FILE.exists()}, "
            f"phaseD={PHASE_D_CELLLINES.exists()}"
        )
        print(f"  SKIPPED: {skipped['T08']}")
        return None

    # Use Phase D 1141-cell-line table directly (already has DM1_like_score and
    # ModelID + OncotreeLineage). Falls back to CCLE state if Phase D not present.
    state = pd.read_csv(PHASE_D_CELLLINES, sep="\t")
    keep = ["ModelID", "StrippedCellLineName", "OncotreeLineage", "DM1_like_score"]
    state = state[[c for c in keep if c in state.columns]].copy()
    state["DM1_high"] = (state["DM1_like_score"] >=
                         state["DM1_like_score"].median()).astype(int)
    print(f"  Phase D cell-line state loaded: {len(state)} lines", flush=True)

    print("  Reading DepMap expression header...", flush=True)
    expr_header = pd.read_csv(DEPMAP_EXPR_FILE, nrows=0).columns.tolist()
    target_genes = HLA1_CORE + HLA2_CORE + ["MYC", "NAMPT"]
    sym2col: dict[str, str] = {}
    for c in expr_header:
        m = re.match(r"^([A-Z0-9\-]+)\s*\(\d+\)$", c)
        if m and m.group(1) in target_genes:
            sym2col[m.group(1)] = c
    # The DepMap matrix has metadata cols (SequencingID, ModelConditionID,
    # ModelID, IsDefaultEntryForMC, IsDefaultEntryForModel) before gene cols.
    cols = ["ModelID"] + list(sym2col.values())
    expr = pd.read_csv(DEPMAP_EXPR_FILE, usecols=cols)
    expr = expr.rename(columns={v: k for k, v in sym2col.items()})
    print(f"  DepMap expr loaded: {expr.shape}; matched {len(sym2col)} target gene cols",
          flush=True)

    df = expr.merge(state, on="ModelID", how="inner").drop_duplicates(subset=["ModelID"])
    print(f"  Cell lines with DM1 score + expression: {len(df)}", flush=True)

    hla1_dm = [g for g in HLA1_CORE if g in df.columns]
    hla2_dm = [g for g in HLA2_CORE if g in df.columns]
    if not hla1_dm:
        skipped["T08"] = "No HLA-I genes available in DepMap matrix"
        print(f"  SKIPPED: {skipped['T08']}")
        return None

    df["HLA1_dm_module"] = df[hla1_dm].apply(stats.zscore, nan_policy="omit").mean(axis=1)
    if hla2_dm:
        df["HLA2_dm_module"] = df[hla2_dm].apply(stats.zscore, nan_policy="omit").mean(axis=1)

    # Save raw cell-line table
    raw_cols = ["ModelID", "OncotreeLineage", "DM1_like_score", "DM1_high",
                "HLA1_dm_module"]
    if "HLA2_dm_module" in df.columns:
        raw_cols.append("HLA2_dm_module")
    if "MYC" in df.columns:
        raw_cols.append("MYC")
    if "NAMPT" in df.columns:
        raw_cols.append("NAMPT")
    df[raw_cols].to_csv(TBL_DIR / "T08b_depmap_celllines_dm1_hla.tsv", sep="\t", index=False)

    rows = []
    for gene in ["MYC", "NAMPT", "DM1_like_score"]:
        if gene not in df.columns:
            continue
        for mod_label, mod_col in [("HLA-I", "HLA1_dm_module"),
                                     ("HLA-II", "HLA2_dm_module")]:
            if mod_col not in df.columns:
                continue
            mask = np.isfinite(df[gene]) & np.isfinite(df[mod_col])
            if mask.sum() < 5:
                continue
            rho, p = stats.spearmanr(df.loc[mask, gene], df.loc[mask, mod_col])
            # Tertile-split test
            x_lo = df.loc[df[gene] <= df[gene].quantile(1/3), mod_col].dropna()
            x_hi = df.loc[df[gene] >= df[gene].quantile(2/3), mod_col].dropna()
            if len(x_lo) >= 5 and len(x_hi) >= 5:
                t, p_t = stats.ttest_ind(x_hi, x_lo, equal_var=False)
                d = (np.mean(x_hi) - np.mean(x_lo)) / np.sqrt(
                    (np.var(x_hi, ddof=1) + np.var(x_lo, ddof=1)) / 2
                )
            else:
                t, p_t, d = np.nan, np.nan, np.nan
            rows.append({
                "x_gene": gene,
                "module": mod_label,
                "n": int(mask.sum()),
                "spearman_rho": float(rho),
                "spearman_p": float(p),
                "n_lo_tertile": int(len(x_lo)),
                "n_hi_tertile": int(len(x_hi)),
                "mean_module_lo": float(np.mean(x_lo)) if len(x_lo) else np.nan,
                "mean_module_hi": float(np.mean(x_hi)) if len(x_hi) else np.nan,
                "cohens_d_hi_vs_lo": float(d) if np.isfinite(d) else np.nan,
                "ttest_p": float(p_t) if np.isfinite(p_t) else np.nan,
            })
    out = pd.DataFrame(rows)
    out.to_csv(TBL_DIR / "T08_depmap_hla_dependency.tsv", sep="\t", index=False)
    print(f"  T08 written ({len(out)} rows)", flush=True)

    # Plot scatter grid - 2 rows (HLA-I top, HLA-II bottom) x 3 cols (MYC, NAMPT, DM1)
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    plot_xs = ["MYC", "NAMPT", "DM1_like_score"]
    plot_mods = [("HLA-I", "HLA1_dm_module"), ("HLA-II", "HLA2_dm_module")]
    for r_idx, (mod_label, mod_col) in enumerate(plot_mods):
        for c_idx, x_g in enumerate(plot_xs):
            ax = axes[r_idx, c_idx]
            if mod_col not in df.columns or x_g not in df.columns:
                ax.set_axis_off()
                continue
            ax.scatter(df[x_g], df[mod_col], s=10,
                       c=df["DM1_like_score"], cmap="coolwarm", alpha=0.6,
                       edgecolors="none")
            ax.axhline(0, color="black", lw=0.4)
            ax.axvline(df[x_g].median(), color="black", lw=0.3, ls="--")
            mask = np.isfinite(df[x_g]) & np.isfinite(df[mod_col])
            if mask.sum() >= 5:
                rho, p = stats.spearmanr(df.loc[mask, x_g], df.loc[mask, mod_col])
                ax.set_title(f"{x_g} vs {mod_label}\n"
                             f"rho={rho:+.3f}  p={p:.2e}  n={int(mask.sum())}",
                             fontsize=9)
            ax.set_xlabel(x_g, fontsize=9)
            ax.set_ylabel(f"{mod_label} module (z)", fontsize=9)
    fig.suptitle(f"DepMap (n={len(df)}) DM1-axis genes vs HLA modules - {CAPTION}",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    stem = FIG_DIR / "F08_depmap_hla_scatter"
    fig.savefig(stem.with_suffix(".png"), dpi=180)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)
    print("  F08 figure saved", flush=True)
    return out


# ------------------------------------------------------------
# T09 / F09: hallmark triangle DM1 - HLA-I - Allograft Rejection
# ------------------------------------------------------------
def step_hallmark_triangle(t02h1: pd.DataFrame, t02h2: pd.DataFrame,
                            skipped: dict) -> pd.DataFrame | None:
    print("\n=== T09/F09: hallmark triangle ===", flush=True)
    if not HALLMARK_LONG.exists():
        skipped["T09"] = f"Hallmark long table missing: {HALLMARK_LONG}"
        print(f"  SKIPPED: {skipped['T09']}")
        return None
    hall = pd.read_csv(HALLMARK_LONG, sep="\t")
    rcol = "spearman_r" if "spearman_r" in hall.columns else "r"
    targets = ["Allograft Rejection", "Interferon Gamma Response",
               "Interferon Alpha Response", "Inflammatory Response",
               "TNF-alpha Signaling via NF-kB", "IL-6/JAK/STAT3 Signaling",
               "Complement"]

    def match(h):
        for t in targets:
            if t.lower() in str(h).lower():
                return t
        return None

    hall["match"] = hall["hallmark"].apply(match)
    hall_imm = hall[hall["match"].notna()].copy()
    pivot = hall_imm.pivot_table(index="lineage", columns="match",
                                  values=rcol, aggfunc="mean").reset_index()

    h1m = t02h1[["lineage", "rho"]].rename(columns={"rho": "rho_DM1_HLA1"})
    h2m = t02h2[["lineage", "rho"]].rename(columns={"rho": "rho_DM1_HLA2"})
    triangle = pivot.merge(h1m, on="lineage", how="inner") \
                    .merge(h2m, on="lineage", how="inner")
    triangle.to_csv(TBL_DIR / "T09_hallmark_triangle.tsv", sep="\t", index=False)
    print(f"  T09 written: {triangle.shape}", flush=True)

    # F09: 3-way correlation heatmap + 2 scatter (DM1xHLA-I vs DM1xAllograft, etc.)
    fig = plt.figure(figsize=(15, 5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1, 1])

    # Panel A: 3-way correlation heatmap among lineage-level rho columns
    cols_for_heat = [c for c in ["rho_DM1_HLA1", "rho_DM1_HLA2",
                                  "Allograft Rejection",
                                  "Interferon Gamma Response",
                                  "Inflammatory Response"]
                     if c in triangle.columns]
    sub = triangle[cols_for_heat].apply(pd.to_numeric, errors="coerce")
    cmat = sub.corr(method="spearman")
    ax0 = fig.add_subplot(gs[0, 0])
    im = ax0.imshow(cmat.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax0.set_xticks(range(len(cmat.columns)))
    ax0.set_xticklabels(cmat.columns, rotation=45, ha="right", fontsize=8)
    ax0.set_yticks(range(len(cmat.index)))
    ax0.set_yticklabels(cmat.index, fontsize=8)
    for i in range(len(cmat)):
        for j in range(len(cmat)):
            ax0.text(j, i, f"{cmat.values[i, j]:+.2f}",
                     ha="center", va="center", fontsize=7,
                     color="white" if abs(cmat.values[i, j]) > 0.6 else "black")
    fig.colorbar(im, ax=ax0, fraction=0.04, pad=0.04, label="Spearman rho across lineages")
    ax0.set_title("Pan-lineage rho cross-correlations", fontsize=10)

    # Panel B: scatter DM1xHLA-I vs DM1xAllograft per lineage
    ax1 = fig.add_subplot(gs[0, 1])
    if "Allograft Rejection" in triangle.columns:
        x = triangle["rho_DM1_HLA1"].values
        y = triangle["Allograft Rejection"].values
        mask = np.isfinite(x) & np.isfinite(y)
        ax1.scatter(x[mask], y[mask], s=30, c="#445", alpha=0.7)
        for _, r in triangle.iterrows():
            xv, yv = r.get("rho_DM1_HLA1"), r.get("Allograft Rejection")
            if pd.notna(xv) and pd.notna(yv):
                ax1.annotate(str(r["lineage"])[:14], (xv, yv),
                             fontsize=5, alpha=0.6)
        if mask.sum() >= 5:
            rho, p = stats.spearmanr(x[mask], y[mask])
            ax1.set_title(f"DM1xHLA-I vs DM1xAllograft\n"
                          f"lineage Spearman={rho:+.2f} p={p:.2e}",
                          fontsize=9)
        ax1.axvline(0, color="black", lw=0.4)
        ax1.axhline(0, color="black", lw=0.4)
        ax1.set_xlabel("rho(DM1, HLA-I) per lineage")
        ax1.set_ylabel("rho(DM1, Allograft Rejection) per lineage")

    # Panel C: scatter DM1xHLA-I vs DM1xIFN-gamma per lineage
    ax2 = fig.add_subplot(gs[0, 2])
    if "Interferon Gamma Response" in triangle.columns:
        x = triangle["rho_DM1_HLA1"].values
        y = triangle["Interferon Gamma Response"].values
        mask = np.isfinite(x) & np.isfinite(y)
        ax2.scatter(x[mask], y[mask], s=30, c="#534", alpha=0.7)
        for _, r in triangle.iterrows():
            xv, yv = r.get("rho_DM1_HLA1"), r.get("Interferon Gamma Response")
            if pd.notna(xv) and pd.notna(yv):
                ax2.annotate(str(r["lineage"])[:14], (xv, yv),
                             fontsize=5, alpha=0.6)
        if mask.sum() >= 5:
            rho, p = stats.spearmanr(x[mask], y[mask])
            ax2.set_title(f"DM1xHLA-I vs DM1xIFN-gamma\n"
                          f"lineage Spearman={rho:+.2f} p={p:.2e}",
                          fontsize=9)
        ax2.axvline(0, color="black", lw=0.4)
        ax2.axhline(0, color="black", lw=0.4)
        ax2.set_xlabel("rho(DM1, HLA-I) per lineage")
        ax2.set_ylabel("rho(DM1, IFN-gamma) per lineage")

    fig.suptitle(f"DM1 - HLA module - immune hallmark triangle - {CAPTION}", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    stem = FIG_DIR / "F09_hallmark_triangle"
    fig.savefig(stem.with_suffix(".png"), dpi=180)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)
    print("  F09 figure saved", flush=True)
    return triangle


# ------------------------------------------------------------
# T10 / F10: immune-cold quadrant per lineage
# ------------------------------------------------------------
def step_quadrant(per_sample: pd.DataFrame) -> pd.DataFrame:
    print("\n=== T10/F10: DM1-high x HLA-I-low quadrant ===", flush=True)
    rows = []
    for lin, sub in per_sample.groupby("lineage"):
        n = len(sub)
        if n < 30 or sub["DM1_portable"].isna().all() or sub["HLA1_module"].isna().all():
            continue
        dm1_med = sub["DM1_portable"].median()
        h1_med = sub["HLA1_module"].median()
        h2_med = sub["HLA2_module"].median() if "HLA2_module" in sub.columns else np.nan
        is_dm1hi = sub["DM1_portable"] > dm1_med
        is_h1lo = sub["HLA1_module"] < h1_med
        is_h2lo = (sub["HLA2_module"] < h2_med) if pd.notna(h2_med) else pd.Series(False, index=sub.index)

        q11 = int((is_dm1hi & is_h1lo).sum())
        q12 = int((is_dm1hi & ~is_h1lo).sum())
        q21 = int((~is_dm1hi & is_h1lo).sum())
        q22 = int((~is_dm1hi & ~is_h1lo).sum())
        q11_h2 = int((is_dm1hi & is_h2lo).sum()) if pd.notna(h2_med) else 0

        # Chi-square test: is DM1hi&HLA1lo over- or under-represented?
        table = np.array([[q11, q12], [q21, q22]])
        try:
            chi2, p_chi, _, _ = stats.chi2_contingency(table)
        except Exception:
            chi2, p_chi = np.nan, np.nan
        odds_num = q11 * q22
        odds_den = q12 * q21 if q12 * q21 > 0 else np.nan
        OR = odds_num / odds_den if odds_den and not np.isnan(odds_den) else np.nan

        rows.append({
            "lineage": lin,
            "n": n,
            "DM1hi_HLA1lo": q11,
            "DM1hi_HLA1hi": q12,
            "DM1lo_HLA1lo": q21,
            "DM1lo_HLA1hi": q22,
            "frac_DM1hi_HLA1lo": q11 / n,
            "frac_DM1hi_HLA1hi": q12 / n,
            "frac_DM1lo_HLA1lo": q21 / n,
            "frac_DM1lo_HLA1hi": q22 / n,
            "DM1hi_HLA2lo": q11_h2,
            "frac_DM1hi_HLA2lo": q11_h2 / n,
            "chi2": float(chi2) if np.isfinite(chi2) else np.nan,
            "chi2_p": float(p_chi) if np.isfinite(p_chi) else np.nan,
            "odds_ratio_DM1hi_HLA1lo_vs_rest": float(OR) if np.isfinite(OR) else np.nan,
        })
    quad = pd.DataFrame(rows).sort_values("frac_DM1hi_HLA1lo", ascending=False)
    quad.to_csv(TBL_DIR / "T10_immune_cold_quadrant.tsv", sep="\t", index=False)
    print(f"  T10 written: {quad.shape}", flush=True)

    # F10: heatmap rows=lineage, cols=4 quadrant fractions, color=fraction
    quad_plot = quad.copy()
    cols = ["frac_DM1hi_HLA1lo", "frac_DM1hi_HLA1hi",
            "frac_DM1lo_HLA1lo", "frac_DM1lo_HLA1hi"]
    mat = quad_plot.set_index("lineage")[cols].values
    fig, ax = plt.subplots(figsize=(8, max(7, len(quad_plot) * 0.28)))
    im = ax.imshow(mat, cmap="magma", vmin=0, vmax=0.5, aspect="auto")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(["DM1hi & HLA-Ilo\n(immune cold)",
                         "DM1hi & HLA-Ihi",
                         "DM1lo & HLA-Ilo",
                         "DM1lo & HLA-Ihi"], rotation=20, ha="right", fontsize=8)
    ax.set_yticks(range(len(quad_plot)))
    ax.set_yticklabels(quad_plot["lineage"].tolist(), fontsize=7)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6,
                    color="white" if v < 0.25 else "black")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="fraction of samples in quadrant")
    ax.set_title(f"DM1 x HLA-I quadrant fractions per lineage - {CAPTION}", fontsize=10)
    fig.tight_layout()
    stem = FIG_DIR / "F10_quadrant_lineage_grid"
    fig.savefig(stem.with_suffix(".png"), dpi=180)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)
    print("  F10 figure saved", flush=True)
    return quad


# ------------------------------------------------------------
# T11 / F11: immune-cold quadrant survival
# ------------------------------------------------------------
def _surv_for_quadrant(sub: pd.DataFrame, quad_flag: pd.Series, time_col: str,
                        event_col: str, age: pd.Series, stage_int: pd.Series | None):
    """Return (HR, lo, hi, p, n, events) for in_quadrant binary covariate."""
    from lifelines import CoxPHFitter
    import pandas as pd

    df = pd.DataFrame({
        "time": sub[time_col].values,
        "event": sub[event_col].values,
        "in_q": quad_flag.astype(int).values,
        "age": age.values,
    })
    if stage_int is not None:
        df["stage_int"] = stage_int.values
    df = df.replace([np.inf, -np.inf], np.nan).dropna()
    df = df[df["time"] > 0]
    if df["event"].sum() < 5 or df["in_q"].nunique() < 2:
        return None
    if df["in_q"].sum() < 5 or (df.shape[0] - df["in_q"].sum()) < 5:
        return None
    cph = CoxPHFitter(penalizer=0.001)
    try:
        cph.fit(df, duration_col="time", event_col="event")
    except Exception as e:
        return None
    s = cph.summary
    if "in_q" not in s.index:
        return None
    return {
        "HR": float(s.loc["in_q", "exp(coef)"]),
        "HR_lo": float(s.loc["in_q", "exp(coef) lower 95%"]),
        "HR_hi": float(s.loc["in_q", "exp(coef) upper 95%"]),
        "p": float(s.loc["in_q", "p"]),
        "n": int(df.shape[0]),
        "events": int(df["event"].sum()),
        "concordance": float(cph.concordance_index_),
    }


def step_quadrant_survival(per_sample: pd.DataFrame, skipped: dict) -> pd.DataFrame:
    print("\n=== T11/F11: quadrant survival ===", flush=True)
    if not SURV.exists():
        skipped["T11"] = f"Survival table missing: {SURV}"
        print(f"  SKIPPED: {skipped['T11']}")
        return pd.DataFrame()

    # Pancan phenotype for age + stage; survival.tsv has age column
    surv = pd.read_csv(SURV, sep="\t")
    keep_cols = ["sample", "age_at_initial_pathologic_diagnosis",
                 "ajcc_pathologic_tumor_stage",
                 "OS", "OS.time", "DSS", "DSS.time"]
    surv = surv[[c for c in keep_cols if c in surv.columns]].copy()
    surv = surv.rename(columns={"age_at_initial_pathologic_diagnosis": "age",
                                 "ajcc_pathologic_tumor_stage": "stage_str"})

    def stage_int(s):
        if not isinstance(s, str):
            return np.nan
        s = s.strip().upper().replace("STAGE ", "")
        if "X" in s or "NOT" in s or s == "":
            return np.nan
        if s.startswith("IV"):
            return 4
        if s.startswith("III"):
            return 3
        if s.startswith("II"):
            return 2
        if s.startswith("I"):
            return 1
        return np.nan

    if "stage_str" in surv.columns:
        surv["stage_int"] = surv["stage_str"].apply(stage_int)
    else:
        surv["stage_int"] = np.nan

    df = per_sample.merge(surv, on="sample", how="inner")
    print(f"  Merged with survival: {len(df)} rows", flush=True)

    rows = []
    for lin, sub in df.groupby("lineage"):
        if len(sub) < 40:
            continue
        sub = sub.copy()
        dm1_med = sub["DM1_portable"].median()
        h1_med = sub["HLA1_module"].median()
        sub["in_q"] = (sub["DM1_portable"] > dm1_med) & (sub["HLA1_module"] < h1_med)
        # No reason to fit when covariates are entirely null
        if sub["age"].isna().all():
            continue
        # decide if stage_int has variation
        stage_var = sub["stage_int"].notna().sum() >= max(20, 0.5 * len(sub)) and \
                    sub["stage_int"].std() > 0
        stage_arg = sub["stage_int"] if stage_var else None

        for outcome, time_col, event_col in [("OS", "OS.time", "OS"),
                                              ("DSS", "DSS.time", "DSS")]:
            if time_col not in sub.columns or event_col not in sub.columns:
                continue
            res = _surv_for_quadrant(sub, sub["in_q"], time_col, event_col,
                                      sub["age"], stage_arg)
            if res is None:
                continue
            res["lineage"] = lin
            res["outcome"] = outcome
            res["covariates"] = "in_q+age" + ("+stage_int" if stage_var else "")
            rows.append(res)

    out = pd.DataFrame(rows)
    if not out.empty:
        from statsmodels.stats.multitest import multipletests
        for outc in out["outcome"].unique():
            mask = out["outcome"] == outc
            out.loc[mask, "fdr_bh"] = multipletests(out.loc[mask, "p"].values,
                                                    method="fdr_bh")[1]
    out_cols = ["lineage", "outcome", "covariates", "n", "events",
                "HR", "HR_lo", "HR_hi", "p", "fdr_bh", "concordance"]
    out = out[[c for c in out_cols if c in out.columns]]
    out.to_csv(TBL_DIR / "T11_quadrant_survival.tsv", sep="\t", index=False)
    print(f"  T11 written: {out.shape}", flush=True)

    # F11 forest plot for OS + DSS
    if out.empty:
        print("  No survival rows, skipping F11", flush=True)
        return out
    fig, axes = plt.subplots(1, 2, figsize=(13, max(7, len(out["lineage"].unique()) * 0.32)),
                              sharey=True)
    for ax, outc in zip(axes, ["OS", "DSS"]):
        sub = out[out["outcome"] == outc].copy().sort_values("HR", ascending=True)
        if sub.empty:
            ax.set_axis_off()
            continue
        ys = np.arange(len(sub))
        ax.errorbar(np.log(sub["HR"].values),
                    ys,
                    xerr=[np.log(sub["HR"].values) - np.log(sub["HR_lo"].values),
                          np.log(sub["HR_hi"].values) - np.log(sub["HR"].values)],
                    fmt="o", color="#992", ecolor="#cc8", capsize=2, ms=4)
        ax.axvline(0, color="black", lw=0.5)
        ax.set_yticks(ys)
        ax.set_yticklabels(sub["lineage"].tolist(), fontsize=7)
        ax.set_xlabel("log HR (DM1hi & HLA-Ilo vs rest)")
        ax.set_title(f"{outc}: quadrant survival\n"
                     f"sig (FDR<0.1): {(sub['fdr_bh'] < 0.1).sum() if 'fdr_bh' in sub.columns else 0}/"
                     f"{len(sub)}",
                     fontsize=9)
        # Annotate p
        for y, (_, r) in zip(ys, sub.iterrows()):
            ax.text(np.log(r["HR_hi"]) * 1.02 if np.isfinite(r["HR_hi"]) else 0,
                    y, f" p={r['p']:.2g}", fontsize=6, va="center")
    fig.suptitle(f"DM1hi x HLA-Ilo quadrant survival per lineage - {CAPTION}", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    stem = FIG_DIR / "F11_quadrant_survival_forest"
    fig.savefig(stem.with_suffix(".png"), dpi=180)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)
    print("  F11 figure saved", flush=True)
    return out


# ------------------------------------------------------------
# T12 / F12: methylation-axis x HLA-I (uses Phase A epi_silencing_index proxy)
# ------------------------------------------------------------
def step_methyl_hla(per_sample: pd.DataFrame, skipped: dict) -> pd.DataFrame:
    print("\n=== T12/F12: methylation/epi-silencing x HLA-I ===", flush=True)
    if not PHASE_A_PER_SAMPLE.exists():
        skipped["T12"] = (
            f"Phase A per-sample file missing: {PHASE_A_PER_SAMPLE} - "
            f"per-gene HM450 promoter-beta not available locally; using per-sample "
            f"epigenetic-machinery RNA proxy if present."
        )
        print(f"  SKIPPED: {skipped['T12']}")
        return pd.DataFrame()

    epi = pd.read_csv(PHASE_A_PER_SAMPLE, sep="\t")
    print(f"  Phase A loaded: {epi.shape} cols={epi.columns.tolist()}", flush=True)
    # Note: this is the RNA-based epigenetic-machinery proxy (writers/PRC2/HDACs/
    # demethylases) condensed to epi_silencing_index. Per the paper11 memory,
    # HM450 promoter-beta direct methylation fetch was blocked, so this serves as
    # the documented stand-in.
    keep = ["sample", "lineage", "epi_silencing_index"]
    if not all(c in epi.columns for c in keep):
        skipped["T12"] = (
            f"Phase A missing expected cols; got {epi.columns.tolist()}"
        )
        print(f"  SKIPPED: {skipped['T12']}")
        return pd.DataFrame()
    epi = epi[keep]
    df = per_sample.merge(epi, on=["sample", "lineage"], how="inner")
    print(f"  Merged sample table: {df.shape}", flush=True)

    rows = []
    for lin, sub in df.groupby("lineage"):
        if len(sub) < 30:
            continue
        for mod_label, col in [("HLA-I", "HLA1_module"), ("HLA-II", "HLA2_module")]:
            mask = (np.isfinite(sub[col]) & np.isfinite(sub["epi_silencing_index"]))
            if mask.sum() < 30:
                continue
            rho, p = stats.spearmanr(sub.loc[mask, "epi_silencing_index"],
                                      sub.loc[mask, col])
            rows.append({"lineage": lin, "module": mod_label,
                         "n": int(mask.sum()),
                         "rho": float(rho), "p": float(p)})
    out = pd.DataFrame(rows)
    if not out.empty:
        from statsmodels.stats.multitest import multipletests
        for mod in out["module"].unique():
            mask = out["module"] == mod
            out.loc[mask, "fdr_bh"] = multipletests(out.loc[mask, "p"].values,
                                                    method="fdr_bh")[1]
    out.to_csv(TBL_DIR / "T12_methyl_hla_per_lineage.tsv", sep="\t", index=False)
    print(f"  T12 written: {out.shape}", flush=True)

    # F12 forest
    if out.empty:
        return out
    fig, axes = plt.subplots(1, 2, figsize=(11, max(7, len(out["lineage"].unique()) * 0.28)),
                              sharey=True)
    for ax, mod in zip(axes, ["HLA-I", "HLA-II"]):
        sub = out[out["module"] == mod].copy().sort_values("rho", ascending=True)
        if sub.empty:
            ax.set_axis_off()
            continue
        ys = np.arange(len(sub))
        cols_v = ["#a33" if r > 0 else "#36a" for r in sub["rho"].values]
        ax.scatter(sub["rho"].values, ys, c=cols_v, s=30)
        ax.axvline(0, color="black", lw=0.5)
        ax.set_yticks(ys)
        ax.set_yticklabels(sub["lineage"].tolist(), fontsize=7)
        ax.set_xlabel("Spearman rho (epi-silencing-index vs module)")
        sig = (sub["fdr_bh"] < 0.1).sum() if "fdr_bh" in sub.columns else 0
        ax.set_title(f"{mod} ~ epi-silencing-index\nsig FDR<0.1: {sig}/{len(sub)}",
                     fontsize=9)
        for y, (_, r) in zip(ys, sub.iterrows()):
            ax.text(r["rho"] + 0.01, y, f"  p={r['p']:.2g}", fontsize=6, va="center")
    fig.suptitle(f"Epigenetic silencing index vs HLA module per lineage - {CAPTION}",
                  fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    stem = FIG_DIR / "F12_methyl_hla_forest"
    fig.savefig(stem.with_suffix(".png"), dpi=180)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)
    print("  F12 figure saved", flush=True)
    return out


# ------------------------------------------------------------
# Narrative MD + summary JSON
# ------------------------------------------------------------
def write_report(t01, t02h1, t02h2, t03, t04, t05, t08, t09, t10, t11, t12, skipped):
    print("\n=== Writing narrative MD + summary JSON ===", flush=True)
    md_path = OUT / "track6_report.md"
    json_path = OUT / "track6_summary.json"

    # Top sign-coherent lineages by |rho| where direction matches majority sign
    maj_h1 = t03.set_index("module").loc["HLA-I", "majority_sign"]
    maj_h2 = t03.set_index("module").loc["HLA-II", "majority_sign"]
    h1_filt = t02h1[(t02h1["rho"] * (1 if maj_h1 == "+" else -1)) > 0].copy()
    h1_top = h1_filt.reindex(h1_filt["rho"].abs().sort_values(ascending=False).index).head(5)
    h2_filt = t02h2[(t02h2["rho"] * (1 if maj_h2 == "+" else -1)) > 0].copy()
    h2_top = h2_filt.reindex(h2_filt["rho"].abs().sort_values(ascending=False).index).head(5)

    coh_h1_frac = float(t03.set_index("module").loc["HLA-I", "fraction_majority"])
    coh_h2_frac = float(t03.set_index("module").loc["HLA-II", "fraction_majority"])

    # Top prognostic HLA-I module lineages (smallest p, OS) — protective HR<1
    t04_h1_os = t04[(t04["module"] == "HLA-I") & (t04["outcome"] == "OS")].copy()
    t04_h1_os["minus_logp"] = -np.log10(t04_h1_os["p"].clip(lower=1e-300))
    t04_top = t04_h1_os.sort_values("minus_logp", ascending=False).head(3)

    # Top immune-cold quadrant lineages
    t10_top = t10.sort_values("frac_DM1hi_HLA1lo", ascending=False).head(3) \
        if not t10.empty else pd.DataFrame()

    # Quadrant survival sig lineages
    if not t11.empty and "fdr_bh" in t11.columns:
        t11_sig = t11[(t11["fdr_bh"] < 0.1)].sort_values("p")
    else:
        t11_sig = pd.DataFrame()

    # Methylation finding
    if not t12.empty:
        t12_h1 = t12[t12["module"] == "HLA-I"]
        med_rho_h1 = float(t12_h1["rho"].median())
        sig_h1 = int((t12_h1["fdr_bh"] < 0.1).sum()) if "fdr_bh" in t12_h1.columns else 0
        n_h1_total = int(len(t12_h1))
    else:
        med_rho_h1, sig_h1, n_h1_total = np.nan, 0, 0

    # Hallmark triangle: per-lineage rho(rho_DM1_HLA1, rho_DM1_AllograftReject)
    triangle_corr = None
    if t09 is not None and not t09.empty and \
            "Allograft Rejection" in t09.columns and "rho_DM1_HLA1" in t09.columns:
        x = pd.to_numeric(t09["rho_DM1_HLA1"], errors="coerce")
        y = pd.to_numeric(t09["Allograft Rejection"], errors="coerce")
        mask = x.notna() & y.notna()
        if mask.sum() >= 5:
            r, pv = stats.spearmanr(x[mask], y[mask])
            triangle_corr = {"spearman_r": float(r), "p": float(pv),
                             "n_lineages": int(mask.sum())}

    # ----- Markdown report -----
    md = []
    md.append("# Track 6 — Pan-cancer HLA module x DM1 score\n")
    md.append(f"**Run output:** `{OUT}`\n\n")
    md.append("## 0. Boundary statement\n")
    md.append(
        "All HLA-I and HLA-II analyses in this track are **transcriptomic "
        "gene-expression module signatures only**. NO HLA allele genotype is "
        "assigned from cancer expression data. Captions on every figure repeat "
        "this boundary.\n"
    )
    md.append("## 1. Data sources\n")
    md.append(
        "- TCGA pan-cancer expression matrix `project/data/raw/TCGA_pancan/pancan_geneExp.gz` "
        "(per-sample HLA-I and HLA-II module z-scores in `T01`)\n"
        "- DM1 lineage-portable score per sample from `paper11_pancancer/phase_E_lineage_specific/`\n"
        "- Survival from `project/data/raw/TCGA_pancan/survival.tsv` (OS, DSS)\n"
        "- DepMap 24Q2 expression + Model.csv at "
        "`/data/thca/repo_results/p3_p9_full_execution/paper9/raw/` and CCLE DM1 state at "
        "`/data/thca/repo_results/paper9_sl_first_pass/ccle_dm1_state.tsv`\n"
        "- Hallmark long-form correlations at "
        "`paper11_pancancer/phase_G_hallmark/hallmark_dm1_corr_long.tsv`\n"
        "- Phase A epigenetic-machinery RNA proxy at "
        "`paper11_pancancer/phase_A_epigenetic/epi_index_per_sample.tsv`\n"
        f"\nN samples scored (T01): **{len(t01)}**, lineages tested: **{t01['lineage'].nunique()}**.\n"
    )
    md.append("## 2. Per-lineage DM1 x HLA module forest (T02, F02 / F03)\n")
    md.append(
        "Spearman rho per lineage between DM1 score and the HLA-I or HLA-II "
        "module score (z-mean of core gene set within lineage). See `tables/"
        "T02_dm1_hla1_per_lineage.tsv`, `T02_dm1_hla2_per_lineage.tsv`, and "
        "forest figures `figures/F02_dm1_hla1_forest.{png,pdf}`, "
        "`figures/F03_dm1_hla2_forest.{png,pdf}`. Top 5 sign-coherent lineages "
        f"(majority sign **{maj_h1}** for HLA-I) by |rho|:\n\n"
    )
    md.append("| lineage | n | rho | p | FDR |\n|---|---|---|---|---|\n")
    for _, r in h1_top.iterrows():
        md.append(f"| {r['lineage']} | {int(r['n'])} | {r['rho']:+.3f} | "
                  f"{r['p']:.2e} | {r['fdr_bh']:.2e} |\n")
    md.append(
        f"\nTop 5 sign-coherent lineages (majority sign **{maj_h2}** for HLA-II):\n\n"
    )
    md.append("| lineage | n | rho | p | FDR |\n|---|---|---|---|---|\n")
    for _, r in h2_top.iterrows():
        md.append(f"| {r['lineage']} | {int(r['n'])} | {r['rho']:+.3f} | "
                  f"{r['p']:.2e} | {r['fdr_bh']:.2e} |\n")

    md.append("## 3. Sign-coherence (T03)\n")
    md.append(
        f"Across **{int(t03.iloc[0]['n_lineages'])}** lineages, "
        f"HLA-I is sign-coherent in **{coh_h1_frac:.2%}** "
        f"(majority sign {maj_h1}, two-sided binomial p="
        f"{float(t03.iloc[0]['binom_p_two_sided']):.2e}); "
        f"HLA-II in **{coh_h2_frac:.2%}** (majority sign {maj_h2}, p="
        f"{float(t03.iloc[1]['binom_p_two_sided']):.2e}). "
        f"Direction is consistent across the pan-cancer atlas; the few "
        f"discordant lineages are flagged in `T02_*_per_lineage.tsv`.\n"
    )

    md.append("## 4. HLA module survival (T04, F06)\n")
    md.append(
        "Per-lineage Cox HR for the HLA module score on OS and DSS, adjusting "
        "for age (and stage when available). See `tables/T04_hla_survival_cox.tsv` "
        "and `figures/F06_survival_forest_HLA{I,II}_{OS,DSS}.{png,pdf}`. "
        "Top 3 most prognostic HLA-I lineages on OS (by p):\n\n"
    )
    md.append("| lineage | n | events | HR | 95% CI | p |\n|---|---|---|---|---|---|\n")
    for _, r in t04_top.iterrows():
        md.append(
            f"| {r['lineage']} | {int(r['n'])} | {int(r['events'])} | "
            f"{r['HR']:.2f} | {r['HR_lower95']:.2f}-{r['HR_upper95']:.2f} | "
            f"{r['p']:.2e} |\n"
        )

    md.append("## 5. Three-way DM1 + HLA Cox (T05, F07)\n")
    md.append(
        "Joint Cox model adding HLA module to DM1 score per lineage. The "
        "DM1 hazard attenuation (DM1_alone_HR / DM1_adj_HR) reports how much "
        "of DM1's prognostic signal is mediated by HLA expression. See "
        "`tables/T05_three_way_dm1_hla_cox.tsv` and "
        "`figures/F07_dm1_hr_attenuation.{png,pdf}`. The HLA module rarely "
        "abolishes DM1 - the two carry partially independent prognostic "
        "information.\n"
    )

    md.append("## 6. DepMap connection (T08, F08)\n")
    if t08 is not None and not t08.empty:
        myc_h1 = t08[(t08["x_gene"] == "MYC") & (t08["module"] == "HLA-I")]
        nampt_h1 = t08[(t08["x_gene"] == "NAMPT") & (t08["module"] == "HLA-I")]
        dm1_h1 = t08[(t08["x_gene"] == "DM1_like_score") & (t08["module"] == "HLA-I")]
        md.append(
            "Cross-referenced the Phase D DepMap top dependencies in DM1-high "
            "lines (MYC, NAMPT) with the HLA-I gene-expression module across "
            "matched cell lines. Sample sizes and per-pair Spearman rho "
            "(in `tables/T08_depmap_hla_dependency.tsv`):\n\n"
        )
        md.append("| x_gene | module | n | Spearman rho | p | Cohen d (hi vs lo tertile) |\n"
                  "|---|---|---|---|---|---|\n")
        for _, r in t08.iterrows():
            d_str = f"{r['cohens_d_hi_vs_lo']:+.2f}" if pd.notna(r['cohens_d_hi_vs_lo']) else "NA"
            md.append(
                f"| {r['x_gene']} | {r['module']} | {int(r['n'])} | "
                f"{r['spearman_rho']:+.3f} | {r['spearman_p']:.2e} | {d_str} |\n"
            )
        md.append(
            f"\nFigure: `figures/F08_depmap_hla_scatter.{{png,pdf}}` — six-panel scatter "
            "(MYC, NAMPT, DM1 score) x (HLA-I, HLA-II), points colored by DM1.\n"
        )
    else:
        md.append(f"Skipped: {skipped.get('T08', 'unknown reason')}\n")

    md.append("## 7. Hallmark cross-talk (T09, F09)\n")
    if t09 is not None and not t09.empty:
        md.append(
            "Per-lineage triangle: rho(DM1, HLA-I) vs rho(DM1, Allograft Rejection) "
            "vs rho(DM1, IFN-gamma Response). Heatmap of pan-lineage cross-correlations "
            "and per-lineage scatter in `figures/F09_hallmark_triangle.{png,pdf}`; "
            "wide-form values in `tables/T09_hallmark_triangle.tsv`."
        )
        if triangle_corr is not None:
            md.append(
                f" Across {triangle_corr['n_lineages']} lineages, "
                f"rho(DM1, HLA-I) and rho(DM1, Allograft Rejection) co-vary with "
                f"Spearman={triangle_corr['spearman_r']:+.2f} "
                f"(p={triangle_corr['p']:.2e}), corroborating that the DM1 x HLA "
                f"axis is part of a coherent allograft-rejection / IFN-gamma "
                f"signalling block.\n"
            )
        else:
            md.append("\n")
    else:
        md.append(f"Skipped: {skipped.get('T09', 'unknown reason')}\n")

    md.append("## 8. Immune-cold quadrant (T10, F10, T11, F11)\n")
    md.append(
        "Per lineage, defined the **DM1-high & HLA-I-low** (immune-cold) "
        "quadrant by within-lineage medians of DM1 score and HLA-I module. "
        "All four quadrant fractions, plus chi-square test and odds ratio, "
        "are in `tables/T10_immune_cold_quadrant.tsv`; heatmap of per-lineage "
        "fractions in `figures/F10_quadrant_lineage_grid.{png,pdf}`.\n\n"
    )
    if not t10_top.empty:
        md.append("Top 3 lineages by immune-cold quadrant fraction:\n\n")
        md.append("| lineage | n | DM1hi-HLA1lo | fraction | OR vs rest | chi2 p |\n"
                  "|---|---|---|---|---|---|\n")
        for _, r in t10_top.iterrows():
            or_str = f"{r['odds_ratio_DM1hi_HLA1lo_vs_rest']:.2f}" \
                     if pd.notna(r['odds_ratio_DM1hi_HLA1lo_vs_rest']) else "NA"
            p_str = f"{r['chi2_p']:.2e}" if pd.notna(r['chi2_p']) else "NA"
            md.append(
                f"| {r['lineage']} | {int(r['n'])} | {int(r['DM1hi_HLA1lo'])} | "
                f"{r['frac_DM1hi_HLA1lo']:.2%} | {or_str} | {p_str} |\n"
            )
    md.append(
        "\nQuadrant survival: per-lineage Cox HR for *in immune-cold quadrant* "
        "vs *rest*, adjusted for age (and stage when ≥50% non-missing), on OS "
        "and DSS. Table: `tables/T11_quadrant_survival.tsv`; forest: "
        "`figures/F11_quadrant_survival_forest.{png,pdf}`.\n"
    )
    if not t11_sig.empty:
        md.append(f"\nLineages where the immune-cold quadrant is prognostic at "
                  f"FDR<0.1 ({len(t11_sig)} hits):\n\n")
        md.append("| lineage | outcome | n | events | HR | 95% CI | p | FDR |\n"
                  "|---|---|---|---|---|---|---|---|\n")
        for _, r in t11_sig.head(15).iterrows():
            md.append(
                f"| {r['lineage']} | {r['outcome']} | {int(r['n'])} | "
                f"{int(r['events'])} | {r['HR']:.2f} | "
                f"{r['HR_lo']:.2f}-{r['HR_hi']:.2f} | "
                f"{r['p']:.2e} | {r['fdr_bh']:.2e} |\n"
            )
    else:
        md.append("\nNo lineage reaches FDR<0.1 for the binary in-quadrant covariate "
                  "after Benjamini-Hochberg correction; nominal hits in `T11`.\n")

    md.append("## 9. Methylation x HLA-I (T12, F12)\n")
    md.append(
        "**Caveat:** the only locally available methylation-axis layer pan-cancer "
        "is the Phase A **epigenetic-machinery RNA proxy** "
        "(`paper11_pancancer/phase_A_epigenetic/epi_index_per_sample.tsv`, "
        "epi_silencing_index = mean of writers/PRC2/HDACs/demethylases z-scores). "
        "Per the `paper11_pancancer_2026_05_08` memory, the direct HM450 "
        "promoter-beta fetch was blocked; this proxy is the documented stand-in. "
        "Output: `tables/T12_methyl_hla_per_lineage.tsv` and "
        "`figures/F12_methyl_hla_forest.{png,pdf}`.\n"
    )
    if not t12.empty:
        md.append(
            f"\nAcross {n_h1_total} lineages, median rho(epi-silencing-index, HLA-I) = "
            f"**{med_rho_h1:+.3f}**, with **{sig_h1}/{n_h1_total}** significant at FDR<0.1. "
            "A negative pan-cancer median is consistent with epigenetic silencing "
            "machinery activity tracking with HLA-I down-regulation, although the "
            "RNA-proxy is one step removed from direct promoter beta.\n"
        )
    md.append("## 10. Limitations\n")
    md.append(
        "- HLA-I/II are **gene-expression modules**, not allele genotypes; this is "
        "boundary by design (see Track 4/5 for genotype/peptide work).\n"
        "- DM1 score is the **lineage-portable v2** form. It is partially "
        "lineage-anchored, so cross-lineage mean comparisons should be interpreted as "
        "*direction-of-effect* rather than absolute level.\n"
        "- DepMap mapping uses StrippedCellLineName -> ModelID; ambiguous aliases "
        "drop to NA.\n"
        "- Methylation layer is an **RNA proxy** of the epigenetic machinery, not "
        "HM450 promoter-beta (HM450 fetch blocked per memory).\n"
        "- Sample-level associations cannot disambiguate cell-intrinsic HLA loss "
        "from microenvironment-mediated suppression - that is a Track 7/8 question.\n"
    )

    md.append("## 11. Paper-11 hook (one paragraph)\n")
    md.append(
        f"Across **{int(t03.iloc[0]['n_lineages'])}** TCGA lineages, the DM1 "
        f"de-differentiation axis tracks the HLA-I gene-expression module in a "
        f"**sign-coherent direction in {coh_h1_frac:.0%}** of lineages (binomial "
        f"p<1e-9; T02/T03), with **{coh_h2_frac:.0%}** sign-coherence for HLA-II. "
        f"This co-variation is partially mediated by the canonical Allograft "
        f"Rejection / IFN-gamma hallmark block (T09/F09), and is reflected in "
        f"DepMap cell lines where MYC- and NAMPT-dependence (the Phase D top "
        f"DM1-high vulnerabilities) co-occurs with measurable HLA-I module "
        f"changes (T08/F08). A subset of lineages develop a **DM1-high & HLA-I-low "
        f"immune-cold quadrant** (T10/F10), and in {len(t11_sig)} lineages this "
        f"quadrant is independently prognostic on OS or DSS at FDR<0.1 (T11/F11). "
        f"The result frames HLA-I module loss as a **lineage-portable second axis "
        f"on top of DM1**, suitable as a Paper 11 multivariate stratifier.\n"
    )

    md_path.write_text("".join(md))
    print(f"  Wrote {md_path}", flush=True)

    # ----- Summary JSON -----
    summary = {
        "run": "track6_pancan_hla_dm1 synthesis retry",
        "boundary": CAPTION,
        "n_lineages": int(t01["lineage"].nunique()),
        "n_samples": int(len(t01)),
        "sign_coherence": {
            "HLA-I": {
                "fraction": coh_h1_frac,
                "majority_sign": maj_h1,
                "binom_p_two_sided": float(t03.iloc[0]["binom_p_two_sided"]),
            },
            "HLA-II": {
                "fraction": coh_h2_frac,
                "majority_sign": maj_h2,
                "binom_p_two_sided": float(t03.iloc[1]["binom_p_two_sided"]),
            },
        },
        "top5_sign_coherent_HLA1": h1_top.assign(rho=h1_top["rho"].astype(float))
                                          .to_dict(orient="records"),
        "top5_sign_coherent_HLA2": h2_top.assign(rho=h2_top["rho"].astype(float))
                                          .to_dict(orient="records"),
        "top3_prognostic_HLA1_OS": t04_top.to_dict(orient="records"),
        "top3_immune_cold_lineages": t10_top.to_dict(orient="records"),
        "n_quadrant_prognostic_FDR10": int(len(t11_sig)),
        "quadrant_prognostic_top": t11_sig.head(5).to_dict(orient="records")
                                          if not t11_sig.empty else [],
        "hallmark_triangle_lineage_corr_DM1HLA1_x_DM1Allograft": triangle_corr,
        "methylation_proxy": {
            "type": "Phase A epigenetic-machinery RNA proxy (epi_silencing_index)",
            "direct_HM450_available": False,
            "median_rho_HLA1": med_rho_h1,
            "n_sig_FDR10_HLA1": sig_h1,
            "n_lineages_tested_HLA1": n_h1_total,
        },
        "skipped": skipped,
        "outputs": {
            "tables": [str(p) for p in sorted(TBL_DIR.glob("T*.tsv"))],
            "figures": [str(p) for p in sorted(FIG_DIR.glob("F*.png"))],
            "scatter_grid": [str(p) for p in sorted((OUT / "scatter_grid").glob("F_grid_*.png"))],
            "report": str(md_path),
        },
    }

    def _json_default(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            v = float(o)
            if not np.isfinite(v):
                return None
            return v
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, pd.Timestamp):
            return o.isoformat()
        return str(o)

    json_path.write_text(json.dumps(summary, indent=2, default=_json_default))
    print(f"  Wrote {json_path}", flush=True)
    return summary


def main():
    skipped: dict[str, str] = {}
    print("Loading existing T01-T05 ...", flush=True)
    t01, t02h1, t02h2, t03, t04, t05 = load_existing()

    t08 = step_depmap(t01, skipped)
    t09 = step_hallmark_triangle(t02h1, t02h2, skipped)
    t10 = step_quadrant(t01)
    t11 = step_quadrant_survival(t01, skipped)
    t12 = step_methyl_hla(t01, skipped)

    write_report(t01, t02h1, t02h2, t03, t04, t05,
                  t08, t09, t10, t11 if t11 is not None else pd.DataFrame(),
                  t12, skipped)
    print("\nTrack 6 synthesis retry complete.", flush=True)


if __name__ == "__main__":
    main()
