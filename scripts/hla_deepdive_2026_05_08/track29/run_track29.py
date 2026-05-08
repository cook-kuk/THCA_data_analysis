"""Track 29 — TCGA-THCA HLA module x neoantigen burden x DM1 (Paper 1, expression-only).

Boundary: HLA gene-expression module ONLY, no germline allele genotype, no carrier
frequency. Uses Thorsson 2018 published per-sample TMB and predicted-neoantigen
counts (computed by Thorsson with reference HLA-A/B/C calls; we re-use them as
a published feature, NOT re-call HLA from cancer BAMs). See
/home/seungho/personal/THCA_data_analysis/project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md.
"""
from __future__ import annotations
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm
from lifelines import CoxPHFitter

warnings.filterwarnings("ignore")

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
THORSSON = ROOT / "project/data/external/thorsson_2018"
TRACK5 = ROOT / "project/results/hla_deepdive_2026_05_08/track5_dm1_hla1_module/tables"
TRACK25 = ROOT / "project/results/hla_deepdive_2026_05_08/track25_hla_somatic_mut/tables"
TRACK6 = ROOT / "project/results/hla_deepdive_2026_05_08/track6_pancan_hla_dm1"
PANCAN = ROOT / "project/results/paper11_pancancer"
OUT = ROOT / "project/results/hla_deepdive_2026_05_08/track29_neoantigen_hla"
FIG = OUT / "figs"
TAB = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

CAPTION = (
    "HLA gene-expression module only (Paper 1 residualization control). "
    "No germline allele genotype, no carrier-frequency claim. "
    "Per HLA_CANCER_SEPARATION_RULES.md."
)


def short_id(s: str, n: int = 12) -> str:
    return str(s)[:n].upper()


def sample_id(s: str, n: int = 15) -> str:
    return str(s)[:n].upper()


def spearman(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 5:
        return np.nan, np.nan, int(m.sum())
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


def cohen_d(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    s = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    if s == 0:
        return np.nan
    return float((a.mean() - b.mean()) / s)


# -----------------------------------------------------------------------------
# 1. Load per-sample THCA panel (Track 5 modules + Track 25 TMB + driver/survival)
# -----------------------------------------------------------------------------
print("[1] Loading THCA per-sample modules...")
t5 = pd.read_csv(TRACK5 / "T01_per_sample_module_scores.tsv", sep="\t")
t5["sample15"] = t5["sample"].astype(str).str.upper().str.slice(0, 15)
t5["patient"] = t5["sample"].astype(str).str.slice(0, 12).str.upper()
n_t5_raw = len(t5)
t5 = t5.drop_duplicates(subset=["sample15"], keep="first").reset_index(drop=True)
print(f"  Track 5 n_raw={n_t5_raw} -> dedup={len(t5)}; cols={list(t5.columns)[:8]}...")

# Track 26: IFN-g signatures
t26 = pd.read_csv(
    ROOT / "project/results/hla_deepdive_2026_05_08/track26_ifng_hla_dm1/tables/T10_per_sample_combined.tsv",
    sep="\t",
)
keep26 = ["sample", "IFNG_hallmark", "IFNG_ayers", "HLA1_recompute", "HLA2_recompute"]
t26 = t26[[c for c in keep26 if c in t26.columns]].copy()
t26["sample15"] = t26["sample"].astype(str).str.upper().str.slice(0, 15)
t26 = t26.drop_duplicates(subset=["sample15"], keep="first").reset_index(drop=True)

# Track 25: TMB (re-computed from MAF) and HLA-pathway mut counts
t25 = pd.read_csv(TRACK25 / "T09_per_sample_merged.tsv", sep="\t")
t25_keep = ["sample", "total_mut", "nonsilent_mut", "classI_nonsilent_count",
            "classI_lof_count", "age"]
t25 = t25[[c for c in t25_keep if c in t25.columns]].copy()
t25["sample15"] = t25["sample"].astype(str).str.upper().str.slice(0, 15)
t25 = t25.drop_duplicates(subset=["sample15"], keep="first").reset_index(drop=True)

# A5 immune (cytolytic_score, HLA-A/B/C, B2M)
a5 = pd.read_csv(ROOT / "project/results/v17p3/tables/A5_immune_subtype_distribution.tsv", sep="\t")
# also cytolytic + HLA expr is in A5_tmb_msi_per_cluster
a5_full = pd.read_csv(ROOT / "project/results/v17p3/tables/A5_tmb_msi_per_cluster.tsv", sep="\t")
a5_keep = ["sample_id", "cytolytic_score", "CD274", "PDCD1", "CTLA4", "IDO1",
           "HLA-A", "HLA-B", "HLA-C", "B2M", "total_mutations"]
a5_full = a5_full[[c for c in a5_keep if c in a5_full.columns]].copy()
a5_full["sample15"] = a5_full["sample_id"].astype(str).str.upper().str.slice(0, 15)
a5_full = a5_full.drop(columns=["sample_id"]).rename(
    columns={"HLA-A": "HLA_A_tpm", "HLA-B": "HLA_B_tpm", "HLA-C": "HLA_C_tpm",
             "B2M": "B2M_tpm", "cytolytic_score": "cyt_a5", "total_mutations": "tmb_a5"}
)
a5_full = a5_full.drop_duplicates(subset=["sample15"], keep="first").reset_index(drop=True)

# -----------------------------------------------------------------------------
# 2. Load Thorsson 2018: TMB, SNV neoantigens, indel neoantigens, leukocyte fraction
# -----------------------------------------------------------------------------
print("[2] Loading Thorsson 2018 panimmune tables...")

mut = pd.read_csv(THORSSON / "mutation_load.tsv", sep="\t")
mut = mut[mut["Cohort"] == "THCA"].copy()
mut["sample15"] = mut["Tumor_Sample_ID"].astype(str).str.upper().str.slice(0, 15)
mut = mut.rename(columns={"Non-silent per Mb": "thorsson_nonsil_perMb",
                          "Silent per Mb": "thorsson_silent_perMb"})
print(f"  Thorsson THCA mutation rows: {len(mut)}")

snv = pd.read_csv(THORSSON / "snv_neoantigens.tsv", sep="\t")
snv["sample15"] = snv["barcode"].astype(str).str.upper().str.slice(0, 15)
snv = snv.rename(
    columns={
        "numberOfNonSynonymousSNP": "n_nonsyn_snv",
        "numberOfImmunogenicMutation": "n_immunogenic_mut",
        "numberOfPeptideTested": "n_peptide_tested",
        "numberOfBindingPMHC": "n_binding_pMHC",
        "numberOfBindingExpressedPMHC": "n_binding_expr_pMHC",
    }
)
indel = pd.read_csv(THORSSON / "indel_neoantigens.tsv", sep="\t")
indel = indel[indel["cancer_type"] == "THCA"].copy()
# Indel file uses 12-char patient IDs, not full sample barcodes — join on patient
indel["patient"] = indel["sample"].astype(str).str.upper().str.slice(0, 12)
indel = indel.rename(columns={"indel_num": "n_indel",
                              "immunogenic_indel_num": "n_immunogenic_indel",
                              "neoantigen_num": "n_neoantigen_indel"})
indel = indel.drop_duplicates(subset=["patient"], keep="first").reset_index(drop=True)

leuk = pd.read_csv(THORSSON / "leukocyte_fraction.tsv", sep="\t",
                   header=None, names=["cohort", "barcode", "leukocyte_frac"])
leuk = leuk[leuk["cohort"] == "THCA"].copy()
leuk["sample15"] = leuk["barcode"].astype(str).str.upper().str.slice(0, 15)

print(f"  THCA: mut={len(mut)}, snv={snv[snv['sample15'].str.contains('TCGA')].shape[0]} (filter), "
      f"indel={len(indel)}, leuk={len(leuk)}")

# Filter snv neoantigens to THCA via merge with mut
snv_thca = snv.merge(mut[["sample15"]], on="sample15", how="inner")
print(f"  Thorsson SNV-neoantigen THCA rows after THCA filter: {len(snv_thca)}")

# -----------------------------------------------------------------------------
# 3. Build merged THCA per-sample frame
# -----------------------------------------------------------------------------
print("[3] Building merged frame...")
df = t5.copy()
for src in [t26.drop(columns=["sample"]), t25.drop(columns=["sample"]),
            a5_full,
            mut[["sample15", "thorsson_nonsil_perMb", "thorsson_silent_perMb"]],
            snv_thca[["sample15", "n_nonsyn_snv", "n_immunogenic_mut",
                      "n_peptide_tested", "n_binding_pMHC", "n_binding_expr_pMHC"]],
            leuk[["sample15", "leukocyte_frac"]]]:
    df = df.merge(src, on="sample15", how="left")
# indel join is patient-level
df = df.merge(indel[["patient", "n_indel", "n_immunogenic_indel",
                     "n_neoantigen_indel"]], on="patient", how="left")

# log neoantigen burden
for c in ["n_immunogenic_mut", "n_binding_pMHC", "n_binding_expr_pMHC",
          "n_nonsyn_snv", "n_neoantigen_indel"]:
    df[c + "_log1p"] = np.log1p(df[c].astype(float))

# Total neoantigen = SNV-binding + indel-binding (simple sum)
df["total_neoantigen_pMHC"] = df["n_binding_pMHC"].fillna(0) + df["n_neoantigen_indel"].fillna(0)
df["total_neoantigen_log1p"] = np.log1p(df["total_neoantigen_pMHC"])

# Crude proxy fallback (will report alongside Thorsson actual)
df["neoantigen_proxy_TMB05"] = df["nonsilent_mut"].astype(float) * 0.05

# Cytolytic from A5 (already GZMA+PRF1 geom-mean per repo convention)
# also recompute from HLA-A/B/C if needed; A5 cyt_a5 is fine.
df["cytolytic"] = df["cyt_a5"]
print(f"  Merged frame: n={len(df)}, with TMB={df['nonsilent_mut'].notna().sum()}, "
      f"with Thorsson SNV neoag={df['n_binding_pMHC'].notna().sum()}, "
      f"with indel neoag={df['n_neoantigen_indel'].notna().sum()}, "
      f"with cytolytic={df['cytolytic'].notna().sum()}, "
      f"with leukocyte={df['leukocyte_frac'].notna().sum()}")

df.to_csv(TAB / "T01_per_sample_merged.tsv", sep="\t", index=False)


def fmt_p(p):
    if p is None or not np.isfinite(p):
        return "NA"
    return f"{p:.2e}" if p < 0.01 else f"{p:.3f}"


# -----------------------------------------------------------------------------
# Deliverable 1: TMB distribution + DM1 tertile
# -----------------------------------------------------------------------------
print("[4] TMB x DM1...")
res = {}
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
ax = axes[0]
tmb = df["nonsilent_mut"].dropna()
ax.hist(tmb, bins=40, color="#5b8def")
ax.set_xlabel("Non-silent mutations (Track 25 MAF)")
ax.set_ylabel("# samples")
ax.set_title(f"THCA TMB distribution (median={tmb.median():.1f}, max={tmb.max():.0f})")
ax = axes[1]
groups = ["Low", "Mid", "High"]
data = [df.loc[df["DM1_tertile"] == g, "nonsilent_mut"].dropna() for g in groups]
ax.boxplot(data, labels=groups)
ax.set_xlabel("DM1 tertile")
ax.set_ylabel("Non-silent mutations")
r, p, n = spearman(df["DM1_use"].values, df["nonsilent_mut"].values)
ax.set_title(f"DM1 vs TMB  rho={r:+.3f}  p={fmt_p(p)}  n={n}")
fig.suptitle("Track 29 Fig 1 — TMB x DM1 in TCGA-THCA  |  " + CAPTION, fontsize=8)
fig.tight_layout()
fig.savefig(FIG / "F01_tmb_distribution_dm1.png", dpi=170, bbox_inches="tight")
plt.close(fig)

res["tmb_dm1"] = {"spearman_rho": r, "p": p, "n": n,
                  "median_low": float(np.nanmedian(data[0])),
                  "median_mid": float(np.nanmedian(data[1])),
                  "median_high": float(np.nanmedian(data[2])),
                  "kruskal_p": float(stats.kruskal(*data).pvalue)}

# -----------------------------------------------------------------------------
# Deliverable 2-3: Neoantigen burden x DM1
# -----------------------------------------------------------------------------
print("[5] Neoantigen x DM1...")
fig, axes = plt.subplots(2, 2, figsize=(11, 8))
metrics = [
    ("n_binding_pMHC_log1p", "log1p(SNV pMHC binding count)"),
    ("n_binding_expr_pMHC_log1p", "log1p(SNV expressed pMHC)"),
    ("n_immunogenic_mut_log1p", "log1p(SNV immunogenic mut)"),
    ("total_neoantigen_log1p", "log1p(SNV+indel neoantigens)"),
]
neoag_rows = []
for ax, (col, lab) in zip(axes.flatten(), metrics):
    valid = df[["DM1_use", col]].dropna()
    if len(valid) < 10:
        ax.set_title(f"{lab}: n<10")
        continue
    ax.scatter(valid["DM1_use"], valid[col], alpha=0.5, s=14, color="#7e3ff2")
    rho, p, n = spearman(valid["DM1_use"].values, valid[col].values)
    pe, pep = stats.pearsonr(valid["DM1_use"], valid[col])
    ax.set_xlabel("DM1 score")
    ax.set_ylabel(lab)
    ax.set_title(f"rho={rho:+.3f} p={fmt_p(p)} n={n} | r={pe:+.3f}")
    neoag_rows.append({"metric": col, "n": n, "spearman_rho": rho, "spearman_p": p,
                       "pearson_r": float(pe), "pearson_p": float(pep)})
fig.suptitle("Track 29 Fig 2 — DM1 x predicted neoantigen burden (Thorsson)  |  " + CAPTION,
             fontsize=8)
fig.tight_layout()
fig.savefig(FIG / "F02_dm1_x_neoantigen_burden.png", dpi=170, bbox_inches="tight")
plt.close(fig)

pd.DataFrame(neoag_rows).to_csv(TAB / "T02_dm1_x_neoantigen.tsv", sep="\t", index=False)
res["neoantigen_dm1"] = neoag_rows

# Tertile contrast for primary metric
prim = "n_binding_pMHC_log1p"
data = [df.loc[df["DM1_tertile"] == g, prim].dropna() for g in groups]
res["neoantigen_dm1_tertile"] = {
    "metric": prim,
    "n_low": int(len(data[0])), "n_mid": int(len(data[1])), "n_high": int(len(data[2])),
    "median_low": float(np.nanmedian(data[0])) if len(data[0]) else np.nan,
    "median_mid": float(np.nanmedian(data[1])) if len(data[1]) else np.nan,
    "median_high": float(np.nanmedian(data[2])) if len(data[2]) else np.nan,
    "cohens_d_high_vs_low": cohen_d(data[2], data[0]),
    "kruskal_p": float(stats.kruskal(*data).pvalue) if all(len(d) > 1 for d in data) else np.nan,
}

# -----------------------------------------------------------------------------
# Deliverable 4: Cytolytic x DM1, Cytolytic x neoantigen, Cytolytic x TMB
# -----------------------------------------------------------------------------
print("[6] Cytolytic chain...")
fig, axes = plt.subplots(1, 3, figsize=(14, 4))
chain = [
    ("DM1_use", "cytolytic", "DM1 vs Cytolytic"),
    ("nonsilent_mut", "cytolytic", "TMB vs Cytolytic"),
    ("n_binding_pMHC_log1p", "cytolytic", "log Neoantigen vs Cytolytic"),
]
chain_rows = []
for ax, (x, y, lab) in zip(axes, chain):
    sub = df[[x, y]].dropna()
    if len(sub) < 5:
        continue
    ax.scatter(sub[x], sub[y], alpha=0.5, s=14, color="#16a085")
    rho, p, n = spearman(sub[x].values, sub[y].values)
    ax.set_xlabel(x); ax.set_ylabel(y)
    ax.set_title(f"{lab}\nrho={rho:+.3f} p={fmt_p(p)} n={n}")
    chain_rows.append({"x": x, "y": y, "n": n, "rho": rho, "p": p})
fig.suptitle("Track 29 Fig 3 — Classical chain pieces  |  " + CAPTION, fontsize=8)
fig.tight_layout()
fig.savefig(FIG / "F03_cytolytic_chain.png", dpi=170, bbox_inches="tight")
plt.close(fig)

pd.DataFrame(chain_rows).to_csv(TAB / "T03_chain_pieces.tsv", sep="\t", index=False)
res["chain_pieces"] = chain_rows

# Also DM1 -> HLA-I module
rho, p, n = spearman(df["DM1_use"].values, df["HLA1_score"].values)
res["dm1_x_hla1"] = {"rho": rho, "p": p, "n": n}
rho2, p2, n2 = spearman(df["DM1_use"].values, df["HLA2_score"].values)
res["dm1_x_hla2"] = {"rho": rho2, "p": p2, "n": n2}

# -----------------------------------------------------------------------------
# Deliverable 5: Mediation chain DM1 -> TMB -> neoantigen -> cytolytic -> HLA-I
# -----------------------------------------------------------------------------
print("[7] Sequential mediation...")

# Use OLS partial regressions; report coefficient of DM1 on each next step
# controlled by previous. DM1 -> TMB (no controls), DM1 -> neo | TMB,
# DM1 -> cyt | TMB, neo, DM1 -> HLA1 | TMB, neo, cyt.
def beta_partial(y, x_list, df_):
    sub = df_[[*x_list, y]].dropna()
    if len(sub) < 25:
        return {"n": int(len(sub)), "beta": np.nan, "p": np.nan, "r2": np.nan,
                "coefs": {}}
    Y = sub[y].astype(float).values
    X = sub[x_list].astype(float).values
    Xc = sm.add_constant(X)
    m = sm.OLS(Y, Xc).fit()
    coefs = dict(zip(["const"] + x_list, m.params.tolist()))
    pvals = dict(zip(["const"] + x_list, m.pvalues.tolist()))
    return {"n": int(len(sub)), "beta_DM1": float(coefs[x_list[0]]),
            "p_DM1": float(pvals[x_list[0]]),
            "r2": float(m.rsquared), "coefs": coefs, "pvals": pvals}


# Standardize chain variables for comparable betas
zdf = df.copy()
for c in ["DM1_use", "nonsilent_mut", "n_binding_pMHC_log1p", "cytolytic",
          "HLA1_score", "HLA2_score", "leukocyte_frac"]:
    if c in zdf.columns:
        v = zdf[c].astype(float)
        zdf[c + "_z"] = (v - v.mean()) / v.std(ddof=0)

mediation = {}
mediation["step1_DM1_to_TMB"] = beta_partial("nonsilent_mut_z", ["DM1_use_z"], zdf)
mediation["step2_DM1_to_neo_given_TMB"] = beta_partial(
    "n_binding_pMHC_log1p_z", ["DM1_use_z", "nonsilent_mut_z"], zdf)
mediation["step3_DM1_to_cyt_given_TMB_neo"] = beta_partial(
    "cytolytic_z", ["DM1_use_z", "nonsilent_mut_z", "n_binding_pMHC_log1p_z"], zdf)
mediation["step4_DM1_to_HLA1_given_TMB_neo_cyt"] = beta_partial(
    "HLA1_score_z", ["DM1_use_z", "nonsilent_mut_z", "n_binding_pMHC_log1p_z", "cytolytic_z"], zdf)
mediation["step4_DM1_to_HLA2_given_TMB_neo_cyt"] = beta_partial(
    "HLA2_score_z", ["DM1_use_z", "nonsilent_mut_z", "n_binding_pMHC_log1p_z", "cytolytic_z"], zdf)
# Direct: DM1 -> HLA1 unconditional
mediation["direct_DM1_to_HLA1"] = beta_partial("HLA1_score_z", ["DM1_use_z"], zdf)
mediation["direct_DM1_to_HLA2"] = beta_partial("HLA2_score_z", ["DM1_use_z"], zdf)

med_rows = []
for k, v in mediation.items():
    med_rows.append({"step": k, "n": v["n"], "beta_DM1": v.get("beta_DM1", np.nan),
                     "p_DM1": v.get("p_DM1", np.nan), "r2": v.get("r2", np.nan)})
pd.DataFrame(med_rows).to_csv(TAB / "T04_mediation_chain.tsv", sep="\t", index=False)
res["mediation"] = med_rows

# Visualize beta cascade
fig, ax = plt.subplots(figsize=(8, 4))
labels = [r["step"].replace("_", "\n") for r in med_rows]
betas = [r["beta_DM1"] for r in med_rows]
colors = ["#1f77b4" if (b is not None and np.isfinite(b) and b > 0) else "#d62728" for b in betas]
ax.bar(range(len(labels)), betas, color=colors)
for i, r in enumerate(med_rows):
    if np.isfinite(r["p_DM1"]):
        ax.text(i, r["beta_DM1"], f"p={fmt_p(r['p_DM1'])}", ha="center",
                va="bottom" if r["beta_DM1"] > 0 else "top", fontsize=8)
ax.axhline(0, color="black", lw=0.5)
ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, fontsize=7, rotation=20, ha="right")
ax.set_ylabel(r"$\beta_{DM1}$ (z-scored)")
ax.set_title("Track 29 Fig 4 — DM1 effect through classical antigen-presentation chain\n"
             "Direct DM1->HLA-I/II independence test  |  " + CAPTION, fontsize=8)
fig.tight_layout()
fig.savefig(FIG / "F04_mediation_beta_cascade.png", dpi=170, bbox_inches="tight")
plt.close(fig)

# -----------------------------------------------------------------------------
# Deliverable 6: Driver-stratified
# -----------------------------------------------------------------------------
print("[8] Driver-stratified...")
drv_rows = []
for drv in ["BRAF", "RAS", "TripleNeg", "Fusion"]:
    sub = df[df["driver_simple"] == drv].copy()
    if len(sub) < 20:
        continue
    for ymet, lab in [("nonsilent_mut", "TMB"),
                      ("n_binding_pMHC_log1p", "neoantigen"),
                      ("cytolytic", "cytolytic"),
                      ("HLA1_score", "HLA1"),
                      ("HLA2_score", "HLA2")]:
        rho, p, n = spearman(sub["DM1_use"].values, sub[ymet].values)
        drv_rows.append({"driver": drv, "y": lab, "n": n, "rho": rho, "p": p,
                         "median_y": float(np.nanmedian(sub[ymet]))})
pd.DataFrame(drv_rows).to_csv(TAB / "T05_driver_stratified.tsv", sep="\t", index=False)
res["driver_stratified"] = drv_rows

# Heatmap-style figure
piv = pd.DataFrame(drv_rows).pivot(index="driver", columns="y", values="rho")
fig, ax = plt.subplots(figsize=(7, 3.5))
im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels(piv.columns)
ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index)
for i in range(piv.shape[0]):
    for j in range(piv.shape[1]):
        v = piv.values[i, j]
        if np.isfinite(v):
            ax.text(j, i, f"{v:+.2f}", ha="center", va="center",
                    color="white" if abs(v) > 0.35 else "black", fontsize=9)
plt.colorbar(im, ax=ax, label="Spearman rho (vs DM1)")
ax.set_title("Track 29 Fig 5 — Driver-stratified DM1 correlations  |  " + CAPTION, fontsize=8)
fig.tight_layout()
fig.savefig(FIG / "F05_driver_stratified_heatmap.png", dpi=170, bbox_inches="tight")
plt.close(fig)

# -----------------------------------------------------------------------------
# Deliverable 7: Inflamed/cold quadrant
# -----------------------------------------------------------------------------
print("[9] Inflamed/cold quadrant...")
quad = df.dropna(subset=["HLA1_score", "n_binding_pMHC_log1p", "cytolytic", "DM1_use"]).copy()
quad["hla_hi"] = quad["HLA1_score"] > quad["HLA1_score"].median()
quad["neo_hi"] = quad["n_binding_pMHC_log1p"] > quad["n_binding_pMHC_log1p"].median()


def label_quad(row):
    if row["hla_hi"] and row["neo_hi"]:
        return "Hot_HLAhi_NEOhi"
    if row["hla_hi"] and not row["neo_hi"]:
        return "Mismatch_HLAhi_NEOlo"
    if not row["hla_hi"] and row["neo_hi"]:
        return "Mismatch_HLAlo_NEOhi"
    return "Cold_HLAlo_NEOlo"


quad["quadrant"] = quad.apply(label_quad, axis=1)

quad_summary = (
    quad.groupby("quadrant")
    .agg(n=("DM1_use", "size"),
         median_DM1=("DM1_use", "median"),
         median_TMB=("nonsilent_mut", "median"),
         median_neo=("n_binding_pMHC", "median"),
         median_cyt=("cytolytic", "median"),
         median_HLA1=("HLA1_score", "median"),
         median_HLA2=("HLA2_score", "median"),
         pct_BRAF=("driver_simple", lambda s: 100 * (s == "BRAF").mean()),
         pct_TripleNeg=("driver_simple", lambda s: 100 * (s == "TripleNeg").mean()))
    .round(3)
)
quad_summary.to_csv(TAB / "T06_quadrant_summary.tsv", sep="\t")
res["quadrant_summary"] = quad_summary.reset_index().to_dict(orient="records")

fig, ax = plt.subplots(figsize=(7, 6))
sc = ax.scatter(quad["n_binding_pMHC_log1p"], quad["HLA1_score"],
                c=quad["DM1_use"], cmap="viridis", s=22, alpha=0.75, edgecolor="k", lw=0.2)
ax.axhline(quad["HLA1_score"].median(), ls="--", color="grey")
ax.axvline(quad["n_binding_pMHC_log1p"].median(), ls="--", color="grey")
ax.set_xlabel("log1p(SNV pMHC binding count)")
ax.set_ylabel("HLA-I module score")
plt.colorbar(sc, ax=ax, label="DM1 score")
ax.set_title("Track 29 Fig 6 — Antigen presentation quadrant (color=DM1)  |  " + CAPTION,
             fontsize=8)
fig.tight_layout()
fig.savefig(FIG / "F06_quadrant_neo_x_hla1_dm1.png", dpi=170, bbox_inches="tight")
plt.close(fig)

# -----------------------------------------------------------------------------
# Deliverable 8: Pan-cancer sanity (top 5 DM1 x HLA-I sign-coherent lineages)
# -----------------------------------------------------------------------------
print("[10] Pan-cancer sanity...")
pancan_neoag = None
try:
    pcd = pd.read_csv(PANCAN / "pancan_dm1_scored.tsv", sep="\t")
    pcd["sample15"] = pcd["sample"].astype(str).str.upper().str.slice(0, 15)
    # join Thorsson SNV neoag pancan
    pcs = snv.copy()
    pcs["sample15"] = pcs["barcode"].astype(str).str.upper().str.slice(0, 15)
    merged = pcd.merge(pcs[["sample15", "n_binding_pMHC", "n_immunogenic_mut"]],
                       on="sample15", how="inner")
    merged["log_pMHC"] = np.log1p(merged["n_binding_pMHC"])
    pan_rows = []
    for lin, sub in merged.groupby("lineage"):
        if len(sub) < 30:
            continue
        rho, p, n = spearman(sub["DM1_like"].values, sub["log_pMHC"].values)
        rho2, p2, _ = spearman(sub["DM1_like"].values, sub["n_immunogenic_mut"].values)
        pan_rows.append({"lineage": lin, "n": n,
                         "rho_DM1_log_pMHC": rho, "p_DM1_log_pMHC": p,
                         "rho_DM1_immunogenic": rho2, "p_DM1_immunogenic": p2})
    pan_df = pd.DataFrame(pan_rows).sort_values("rho_DM1_log_pMHC")
    pan_df.to_csv(TAB / "T07_pancan_dm1_neoantigen.tsv", sep="\t", index=False)
    pancan_neoag = pan_df

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.barh(pan_df["lineage"], pan_df["rho_DM1_log_pMHC"],
            color=["#1f77b4" if r > 0 else "#d62728" for r in pan_df["rho_DM1_log_pMHC"]])
    ax.axvline(0, color="black", lw=0.5)
    ax.set_xlabel("Spearman rho: DM1 vs log1p(SNV pMHC)")
    ax.set_title("Track 29 Fig 7 — Pan-cancer DM1 x neoantigen burden\n"
                 + "(THCA highlighted; expect coherent negative if THCA chain holds)\n  |  "
                 + CAPTION, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "F07_pancan_dm1_x_neoantigen.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    res["pancan_dm1_neoantigen_top"] = pan_df.head(10).to_dict(orient="records")
    res["pancan_dm1_neoantigen_bottom"] = pan_df.tail(10).to_dict(orient="records")
except Exception as e:
    print(f"  pancan sanity failed: {e}")
    res["pancan_dm1_neoantigen_error"] = str(e)

# -----------------------------------------------------------------------------
# Deliverable 9: HLA-I / neoantigen ratio — antigen-presentation efficiency
# -----------------------------------------------------------------------------
print("[11] Antigen-presentation efficiency...")
eff = df.dropna(subset=["HLA1_score", "n_binding_pMHC"]).copy()
# ratio = HLA1 (z-scored) / log1p(neoantigen + 1)
eff["pres_efficiency"] = eff["HLA1_score"] / (np.log1p(eff["n_binding_pMHC"]) + 1.0)
rho_eff, p_eff, n_eff = spearman(eff["DM1_use"].values, eff["pres_efficiency"].values)

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
ax = axes[0]
ax.scatter(eff["DM1_use"], eff["pres_efficiency"], s=14, alpha=0.5, color="#e67e22")
ax.set_xlabel("DM1 score"); ax.set_ylabel("HLA1 / log1p(neoantigen+1)")
ax.set_title(f"Antigen-presentation efficiency vs DM1\nrho={rho_eff:+.3f} p={fmt_p(p_eff)} n={n_eff}")
ax = axes[1]
data = [eff.loc[eff["DM1_tertile"] == g, "pres_efficiency"].dropna() for g in groups]
ax.boxplot(data, labels=groups)
ax.set_xlabel("DM1 tertile"); ax.set_ylabel("Presentation efficiency")
kp = stats.kruskal(*data).pvalue if all(len(d) > 1 for d in data) else np.nan
d_hl = cohen_d(data[2], data[0])
ax.set_title(f"By tertile  Kruskal p={fmt_p(kp)}  d(High-Low)={d_hl:+.3f}")
fig.suptitle("Track 29 Fig 8 — HLA-I / neoantigen ratio (presentation efficiency)  |  " + CAPTION,
             fontsize=8)
fig.tight_layout()
fig.savefig(FIG / "F08_presentation_efficiency.png", dpi=170, bbox_inches="tight")
plt.close(fig)

eff_tab = pd.DataFrame({
    "tertile": groups,
    "n": [len(d) for d in data],
    "median_efficiency": [float(np.nanmedian(d)) if len(d) else np.nan for d in data],
    "iqr_25": [float(np.nanpercentile(d, 25)) if len(d) else np.nan for d in data],
    "iqr_75": [float(np.nanpercentile(d, 75)) if len(d) else np.nan for d in data],
})
eff_tab.to_csv(TAB / "T08_presentation_efficiency.tsv", sep="\t", index=False)
res["presentation_efficiency"] = {
    "rho_DM1": rho_eff, "p": p_eff, "n": n_eff,
    "kruskal_p": float(kp) if np.isfinite(kp) else None,
    "cohens_d_high_vs_low": d_hl,
    "median_low": eff_tab.loc[0, "median_efficiency"],
    "median_high": eff_tab.loc[2, "median_efficiency"],
}

# -----------------------------------------------------------------------------
# Deliverable 10: Survival — DM1 x neoantigen interaction
# -----------------------------------------------------------------------------
print("[12] Survival...")
surv = df.dropna(subset=["OS", "OS.time", "DM1_use", "n_binding_pMHC"]).copy()
surv = surv[surv["OS.time"] > 0].copy()
surv["log_neo"] = np.log1p(surv["n_binding_pMHC"])
surv["DM1_z"] = (surv["DM1_use"] - surv["DM1_use"].mean()) / surv["DM1_use"].std(ddof=0)
surv["log_neo_z"] = (surv["log_neo"] - surv["log_neo"].mean()) / surv["log_neo"].std(ddof=0)
surv["DM1_x_neo"] = surv["DM1_z"] * surv["log_neo_z"]
surv["age_z"] = ((surv["age"].fillna(surv["age"].median()) - surv["age"].mean()) /
                 surv["age"].std(ddof=0))

surv_rows = []
for ev_col, t_col, label in [("OS", "OS.time", "OS"), ("PFI", "PFI.time", "PFI")]:
    sub = surv.dropna(subset=[ev_col, t_col]).copy()
    sub = sub[sub[t_col] > 0]
    if sub[ev_col].sum() < 5:
        surv_rows.append({"endpoint": label, "n": int(len(sub)), "events": int(sub[ev_col].sum()),
                          "note": "too few events"})
        continue
    cph = CoxPHFitter(penalizer=0.01)
    cox_df = sub[[t_col, ev_col, "DM1_z", "log_neo_z", "DM1_x_neo", "age_z"]].rename(
        columns={t_col: "T", ev_col: "E"}).dropna()
    try:
        cph.fit(cox_df, duration_col="T", event_col="E")
        out = cph.summary
        for var in ["DM1_z", "log_neo_z", "DM1_x_neo", "age_z"]:
            if var in out.index:
                surv_rows.append({
                    "endpoint": label,
                    "n": int(len(cox_df)),
                    "events": int(cox_df["E"].sum()),
                    "var": var,
                    "HR": float(out.loc[var, "exp(coef)"]),
                    "HR_ci_low": float(out.loc[var, "exp(coef) lower 95%"]),
                    "HR_ci_high": float(out.loc[var, "exp(coef) upper 95%"]),
                    "p": float(out.loc[var, "p"]),
                })
    except Exception as e:
        surv_rows.append({"endpoint": label, "n": int(len(cox_df)), "note": f"cox failed: {e}"})

pd.DataFrame(surv_rows).to_csv(TAB / "T09_survival_interaction.tsv", sep="\t", index=False)
res["survival"] = surv_rows

# Survival figure: KM by DM1 x neoantigen quadrant
try:
    from lifelines import KaplanMeierFitter
    kmf = KaplanMeierFitter()
    s2 = surv.copy()
    s2["DM1_hi"] = s2["DM1_use"] > s2["DM1_use"].median()
    s2["NEO_hi"] = s2["log_neo"] > s2["log_neo"].median()
    s2["grp"] = s2.apply(lambda r: f"{'DM1hi' if r['DM1_hi'] else 'DM1lo'}_"
                                   f"{'NEOhi' if r['NEO_hi'] else 'NEOlo'}", axis=1)
    fig, ax = plt.subplots(figsize=(7, 5))
    for g, sub in s2.groupby("grp"):
        if len(sub) < 5:
            continue
        kmf.fit(sub["OS.time"] / 30.4, sub["OS"], label=f"{g} (n={len(sub)})")
        kmf.plot_survival_function(ax=ax, ci_show=False)
    ax.set_xlabel("Months from diagnosis"); ax.set_ylabel("OS prob")
    ax.set_title("Track 29 Fig 9 — KM by DM1 x neoantigen burden  |  " + CAPTION, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "F09_KM_dm1_x_neoantigen.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
except Exception as e:
    print(f"  KM failed: {e}")

# -----------------------------------------------------------------------------
# Deliverable 11: Crude proxy comparison (TMB*0.05 vs Thorsson actual)
# -----------------------------------------------------------------------------
proxy = df.dropna(subset=["nonsilent_mut", "n_binding_pMHC"]).copy()
rho_proxy, p_proxy, n_proxy = spearman(
    proxy["neoantigen_proxy_TMB05"].values, proxy["n_binding_pMHC"].values)
res["proxy_vs_thorsson"] = {"rho": rho_proxy, "p": p_proxy, "n": n_proxy,
                            "comment": "Crude TMB*0.05 vs Thorsson published pMHC count."}

# -----------------------------------------------------------------------------
# Final summary
# -----------------------------------------------------------------------------
res["track"] = "Track 29 - DM1 × neoantigen × HLA-I in TCGA-THCA"
res["boundary"] = CAPTION
res["n_thca"] = int(len(df))
res["n_with_thorsson_neo"] = int(df["n_binding_pMHC"].notna().sum())
res["n_with_tmb"] = int(df["nonsilent_mut"].notna().sum())
res["n_with_cytolytic"] = int(df["cytolytic"].notna().sum())
res["outputs_dir"] = str(OUT)


def _convert(o):
    import numpy as _np
    if isinstance(o, (_np.integer,)):
        return int(o)
    if isinstance(o, (_np.floating,)):
        return float(o)
    if isinstance(o, _np.ndarray):
        return o.tolist()
    if isinstance(o, dict):
        return {k: _convert(v) for k, v in o.items()}
    if isinstance(o, list):
        return [_convert(v) for v in o]
    return o


(OUT / "track29_summary.json").write_text(json.dumps(_convert(res), indent=2, default=str))

print("\n=== Track 29 done ===")
print(f"  out: {OUT}")
print(f"  n_thca: {res['n_thca']}")
print(f"  n_thorsson_neo: {res['n_with_thorsson_neo']}")
print(f"  DM1 x TMB rho: {res['tmb_dm1']['spearman_rho']:.3f}")
if neoag_rows:
    print(f"  DM1 x neoantigen (pMHC log) rho: {neoag_rows[0]['spearman_rho']:.3f}")
print(f"  DM1 x cytolytic rho: {chain_rows[0]['rho']:.3f}")
print(f"  DM1 x HLA1 (direct) rho: {res['dm1_x_hla1']['rho']:.3f}")
print(f"  DM1 -> HLA1 | TMB,neo,cyt residual beta: "
      f"{mediation['step4_DM1_to_HLA1_given_TMB_neo_cyt']['beta_DM1']:.3f} "
      f"p={mediation['step4_DM1_to_HLA1_given_TMB_neo_cyt']['p_DM1']:.2e}")
print(f"  Presentation efficiency High vs Low Cohen's d: "
      f"{res['presentation_efficiency']['cohens_d_high_vs_low']:.3f}")
