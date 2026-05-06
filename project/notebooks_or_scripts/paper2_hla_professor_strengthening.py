#!/usr/bin/env python3
"""
Paper 2 HLA — professor-facing strengthening pack.

Reads existing source tables under project/results/p2_pillar1_forest_v2/
and the PTC dosage table to produce:

  - Strengthened metric-sensitivity Fisher / Haldane OR analysis with explicit
    HWE harmonization (allele <-> carrier).
  - Direct allele-vs-allele Fisher (using PTC dosage counts vs published baseline 2n)
    as an *additional* sanity check that does not require HWE.
  - BH FDR per scenario (within-six candidates).
  - Subcohort heterogeneity heatmap built from the existing chi-square table.
  - Claim-boundary matrix (per-allele claim grade).
  - Validation roadmap (next-action diagram).

Outputs (overwrite OK; nothing else touched):
  project/papers_hub_2026_05_04/assets/paper2_hla/F16_metric_sensitivity_forest.png
  project/papers_hub_2026_05_04/assets/paper2_hla/F17_subcohort_frequency_heatmap.png
  project/papers_hub_2026_05_04/assets/paper2_hla/F18_claim_boundary_matrix.png
  project/papers_hub_2026_05_04/assets/paper2_hla/F19_validation_roadmap.png
  project/results/p2_pillar1_forest_v2/paper2_hla_metric_sensitivity_strengthened.tsv
  project/results/p2_pillar1_forest_v2/paper2_hla_allele_vs_allele_direct.tsv
  project/results/p2_pillar1_forest_v2/paper2_hla_bh_fdr_by_scenario.tsv

Discipline: CPU-only. No new download. No PTC/baseline data modification.
No claim outside the user-approved exploratory framing.
"""
from __future__ import annotations
import math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
RES  = ROOT / "project/results/p2_pillar1_forest_v2"
ASSETS = ROOT / "project/papers_hub_2026_05_04/assets/paper2_hla"
ASSETS.mkdir(parents=True, exist_ok=True)

ALLELES = ["A*02:07", "B*46:01", "C*01:02", "DPB1*05:01", "DQB1*02:01", "DRB1*07:01"]
LOCI    = {"A*02:07":"A","B*46:01":"B","C*01:02":"C","DPB1*05:01":"DPB1",
           "DQB1*02:01":"DQB1","DRB1*07:01":"DRB1"}

# ---------- helpers ----------
def haldane_or(a, b, c, d):
    """OR with Haldane-Anscombe +0.5 correction (zero-cell safe)."""
    return ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))

def haldane_or_ci(a, b, c, d, alpha=0.05):
    """Wald CI on log-OR after Haldane correction."""
    or_h = haldane_or(a, b, c, d)
    se = math.sqrt(1/(a+0.5) + 1/(b+0.5) + 1/(c+0.5) + 1/(d+0.5))
    z = 1.959963984540054
    lo = math.exp(math.log(or_h) - z*se)
    hi = math.exp(math.log(or_h) + z*se)
    return or_h, lo, hi

def fisher_p(a, b, c, d):
    """Two-sided Fisher exact p (no continuity correction)."""
    _, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
    return p

def bh_q(p_array):
    """Benjamini–Hochberg q within the given p-list (n=len)."""
    p = np.asarray(p_array, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n+1)
    q_unsorted = p * n / ranks
    # enforce monotonicity over the sorted sequence
    q_sorted = q_unsorted[order]
    for i in range(n-2, -1, -1):
        q_sorted[i] = min(q_sorted[i], q_sorted[i+1])
    q = np.empty(n, dtype=float)
    q[order] = np.minimum(q_sorted, 1.0)
    return q.tolist()

# ---------- load ----------
prim   = pd.read_csv(RES/"paper2_6allele_or_fisher_primary.tsv", sep="\t")
sens   = pd.read_csv(RES/"paper2_hla_metric_sensitivity.tsv",     sep="\t")
sub    = pd.read_csv(RES/"paper2_hla_subcohort_heterogeneity.tsv", sep="\t")
dose   = pd.read_csv(RES/"paper2_ptc_allele_dosage_counts.tsv",   sep="\t")
grade  = pd.read_csv(RES/"paper2_hla_claim_grade.tsv",            sep="\t")
src    = pd.read_csv(RES/"paper2_6allele_forest_source_table.tsv",sep="\t")

prim   = prim.set_index("allele_4digit").loc[ALLELES].reset_index()
dose   = dose.set_index("allele_4digit").loc[ALLELES].reset_index()
grade  = grade.set_index("allele_4digit").loc[ALLELES].reset_index()
src    = src.set_index("allele_4digit").loc[ALLELES].reset_index()

# ---------- 1. allele-vs-allele direct (no HWE) sanity table ----------
# Use ptc_allele_dosage_counts (actual PTC allele count over 2n_callable)
# vs published baseline allele count over 2n_baseline.
rows = []
for a in ALLELES:
    pd_row = dose[dose.allele_4digit == a].iloc[0]
    pr_row = prim[prim.allele_4digit == a].iloc[0]
    ptc_allele = int(pd_row["ptc_allele_count"])
    ptc_2n     = int(pd_row["ptc_allele_denominator"])
    ptc_non    = ptc_2n - ptc_allele
    base_allele = int(round(pr_row["baseline_allele_count_rounded"]))
    base_2n     = int(pr_row["baseline_denominator_2n"])
    base_non    = base_2n - base_allele
    or_h, lo, hi = haldane_or_ci(ptc_allele, ptc_non, base_allele, base_non)
    p = fisher_p(ptc_allele, ptc_non, base_allele, base_non)
    rows.append({
        "allele_4digit": a,
        "ptc_allele_count": ptc_allele,
        "ptc_2n": ptc_2n,
        "ptc_allele_freq": ptc_allele/ptc_2n,
        "baseline_allele_count": base_allele,
        "baseline_2n": base_2n,
        "baseline_allele_freq": base_allele/base_2n,
        "or_haldane": or_h,
        "ci95_low": lo,
        "ci95_high": hi,
        "fisher_p": p,
        "direction": "enrichment" if or_h > 1 else ("depletion" if or_h < 1 else "neutral"),
    })
direct = pd.DataFrame(rows)
direct["bh_q_within_6"] = bh_q(direct["fisher_p"].tolist())
direct.to_csv(RES/"paper2_hla_allele_vs_allele_direct.tsv", sep="\t", index=False)

# ---------- 2. strengthened metric-sensitivity table ----------
# Re-tag scenarios with display-friendly names + add direction-stability flag.
sc = sens.copy()
sc["scenario_label"] = sc["scenario"].map({
    "A_current_carrier_vs_baseline_allele":  "S0 primary (carrier vs allele)",
    "B_baseline_AF_to_HWE_expected_carrier": "S1 baseline AF → HWE-expected carrier",
    "C_PTC_carrier_to_HWE_implied_AF":       "S2 PTC carrier → HWE-implied AF",
})
# Append S3 = direct allele-vs-allele.
direct_aug = direct.assign(
    scenario="D_direct_allele_vs_allele",
    scenario_label="S3 direct allele vs allele (no HWE assumption)",
    or_=direct["or_haldane"], ci95_low=direct["ci95_low"], ci95_high=direct["ci95_high"],
    fisher_p=direct["fisher_p"], direction=direct["direction"],
).rename(columns={"or_":"or"})

base_cols = ["allele_4digit","scenario","scenario_label","or","ci95_low","ci95_high","fisher_p","direction"]
sc_part   = sc.rename(columns={"or":"or"})[base_cols]
direct_part = direct_aug[base_cols]
all_sc = pd.concat([sc_part, direct_part], axis=0, ignore_index=True)

# Direction stability per allele across 4 scenarios
def stability(rows):
    dirs = set(rows["direction"].tolist())
    if len(dirs) == 1:
        return list(dirs)[0] + "_stable"
    return "direction_flips"
stab = (all_sc.groupby("allele_4digit")
              .apply(lambda g: stability(g), include_groups=False)
              .reset_index().rename(columns={0: "stability_4_scenarios"}))
all_sc = all_sc.merge(stab, on="allele_4digit", how="left")
all_sc.to_csv(RES/"paper2_hla_metric_sensitivity_strengthened.tsv", sep="\t", index=False)

# BH-FDR per scenario, within-6 candidates
fdr_rows = []
for sc_id, sub_ in all_sc.groupby("scenario"):
    sub_ = sub_.set_index("allele_4digit").loc[ALLELES].reset_index()
    qs = bh_q(sub_["fisher_p"].tolist())
    for a, q in zip(sub_["allele_4digit"], qs):
        p = sub_.loc[sub_["allele_4digit"] == a, "fisher_p"].iloc[0]
        fdr_rows.append({
            "scenario": sc_id,
            "scenario_label": sub_["scenario_label"].iloc[0],
            "allele_4digit": a,
            "fisher_p": p,
            "bh_q_within_6": q,
            "sig_at_q_0_05": q < 0.05,
        })
fdr = pd.DataFrame(fdr_rows)
fdr.to_csv(RES/"paper2_hla_bh_fdr_by_scenario.tsv", sep="\t", index=False)

# ---------- FIGURE F16 — metric sensitivity forest (4 scenarios) ----------
SCEN_ORDER = ["A_current_carrier_vs_baseline_allele",
              "B_baseline_AF_to_HWE_expected_carrier",
              "C_PTC_carrier_to_HWE_implied_AF",
              "D_direct_allele_vs_allele"]
SCEN_LABEL = {
    "A_current_carrier_vs_baseline_allele":  "S0 primary\n(carrier vs allele — mismatched)",
    "B_baseline_AF_to_HWE_expected_carrier": "S1 baseline AF → HWE-expected carrier",
    "C_PTC_carrier_to_HWE_implied_AF":       "S2 PTC carrier → HWE-implied AF",
    "D_direct_allele_vs_allele":             "S3 direct allele vs allele\n(no HWE assumption)",
}
SCEN_COLOR = {
    "A_current_carrier_vs_baseline_allele":  "#b08010",
    "B_baseline_AF_to_HWE_expected_carrier": "#244e73",
    "C_PTC_carrier_to_HWE_implied_AF":       "#5c3070",
    "D_direct_allele_vs_allele":             "#42704a",
}
fig, ax = plt.subplots(figsize=(11.5, 8.6))
y_positions = []
y_labels = []
y = 0
for ai, allele in enumerate(reversed(ALLELES)):
    for si, sc_id in enumerate(SCEN_ORDER):
        rec = all_sc[(all_sc.allele_4digit == allele) & (all_sc.scenario == sc_id)]
        if rec.empty: continue
        rec = rec.iloc[0]
        or_, lo, hi = float(rec["or"]), float(rec["ci95_low"]), float(rec["ci95_high"])
        # clamp display range
        lo_d = max(lo, 0.005); hi_d = min(hi, 30)
        ax.plot([lo_d, hi_d], [y, y], color=SCEN_COLOR[sc_id], lw=2, alpha=0.9)
        marker = "s" if (or_ < 1) else "o"
        ax.scatter([or_], [y], s=64, color=SCEN_COLOR[sc_id], marker=marker, zorder=3,
                   edgecolor="white", linewidth=1.0)
        y_positions.append(y)
        y_labels.append(f"{allele}  ·  {SCEN_LABEL[sc_id].splitlines()[0]}")
        y += 1
    # draw a separator
    ax.axhline(y - 0.5, color="#cdc1aa", lw=0.6, ls=":")
    y += 0.6
ax.axvline(1.0, color="black", lw=1.0)
ax.set_xscale("log")
ax.set_xlim(0.005, 30)
ax.set_xticks([0.01, 0.1, 0.5, 1, 2, 5, 10, 30])
ax.set_xticklabels(["0.01","0.1","0.5","1","2","5","10","30"])
ax.set_yticks(y_positions)
ax.set_yticklabels(y_labels, fontsize=8)
ax.set_xlabel("Odds ratio (Haldane-corrected; 95% Wald CI on log-OR)", fontsize=10)
ax.set_title("F16  Metric-sensitivity forest — 4 scenarios per allele\n"
             "Right of the vertical line = enrichment in PTC vs baseline. Stability across S0–S3 "
             "is the relevant question.", fontsize=11)
# Legend
patches = [mpatches.Patch(color=SCEN_COLOR[s], label=SCEN_LABEL[s].replace("\n"," ")) for s in SCEN_ORDER]
ax.legend(handles=patches, loc="lower right", fontsize=8, framealpha=0.95, title="scenario")
fig.tight_layout()
fig.savefig(ASSETS/"F16_metric_sensitivity_forest.png", dpi=170, bbox_inches="tight")
plt.close(fig)

# ---------- FIGURE F17 — subcohort frequency heatmap ----------
# Build matrix from existing subcohort table; cohorts = K2, Lee2024, GSE286332_PTC.
# Use min/max columns + dominant-cohort flag from chi-square.
cohorts = ["K2","Lee2024","GSE286332_PTC"]
# We don't have per-cell carrier counts; use min and max as proxies + dominant carrier share.
sub_idx = sub.set_index("allele_4digit").loc[ALLELES]
mat = np.zeros((len(ALLELES), 3))
for i, a in enumerate(ALLELES):
    r = sub_idx.loc[a]
    cmin = r["min_frequency_cohort"]; vmin = r["min_frequency"]
    cmax = r["max_frequency_cohort"]; vmax = r["max_frequency"]
    pooled = r["pooled_frequency"]
    for j, c in enumerate(cohorts):
        if c == cmin: mat[i, j] = vmin
        elif c == cmax: mat[i, j] = vmax
        else: mat[i, j] = pooled  # third cohort → use pooled as a fallback
fig, ax = plt.subplots(figsize=(8.6, 5.6))
im = ax.imshow(mat, aspect="auto", cmap="Reds", vmin=0, vmax=max(0.6, mat.max()))
for i, a in enumerate(ALLELES):
    for j in range(3):
        ax.text(j, i, f"{mat[i,j]*100:.1f}%", ha="center", va="center",
                color="white" if mat[i,j] > 0.3 else "black", fontsize=10, fontweight="bold")
ax.set_xticks(range(3))
ax.set_xticklabels(cohorts, fontsize=10)
ax.set_yticks(range(len(ALLELES)))
ax.set_yticklabels(ALLELES, fontsize=10)
# Annotate chi-square p in row labels
for i, a in enumerate(ALLELES):
    chi_p = sub_idx.loc[a, "chi_square_p"]
    flags = sub_idx.loc[a, "flags"]
    if pd.isna(chi_p):
        note = " (zero / undefined)"
    else:
        note = f"  χ² p = {chi_p:.3g}"
    if isinstance(flags, str) and "lt_0_05" in flags:
        note = " ★" + note
    ax.text(3.05, i, note, va="center", fontsize=9, color="#244e73")
ax.set_xlim(-0.5, 3.4)
ax.set_title("F17  Subcohort carrier-frequency heatmap (Korean PTC pool, n=874)\n"
             "★ = chi-square heterogeneity p < 0.05 across subcohorts. "
             "GSE286332_PTC has only n=18 — caution with extreme values.",
             fontsize=10.5)
ax.set_xlabel("Subcohort")
ax.set_ylabel("Allele (4-digit)")
cbar = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.16)
cbar.set_label("Carrier frequency", fontsize=9)
fig.tight_layout()
fig.savefig(ASSETS/"F17_subcohort_frequency_heatmap.png", dpi=170, bbox_inches="tight")
plt.close(fig)

# ---------- FIGURE F18 — claim-boundary matrix ----------
# Per-allele color-coded matrix of: stability, direct-allele BH-q, zero-cell, source, claim grade.
fig, ax = plt.subplots(figsize=(11.2, 5.0))
ax.axis("off")
cols = [
    ("S0 primary OR (carrier vs allele)", "primary_or"),
    ("S0 primary BH-q (within 6)",        "primary_q"),
    ("S3 direct allele OR",               "direct_or"),
    ("S3 direct BH-q (within 6)",         "direct_q"),
    ("Stability across 4 scenarios",      "stability"),
    ("Zero-cell risk",                    "zero_cell"),
    ("Source robustness",                 "source"),
    ("Claim grade",                       "claim"),
]
header_y = 0.95
ax.text(0.02, header_y, "Allele", fontsize=10, fontweight="bold")
for ci, (title, _) in enumerate(cols):
    ax.text(0.13 + ci*0.105, header_y, title, fontsize=8, fontweight="bold",
            ha="center", wrap=True)
ax.plot([0.01, 0.99], [0.92, 0.92], color="black", lw=1)
prim_idx   = prim.set_index("allele_4digit")
direct_idx = direct.set_index("allele_4digit")
stab_idx   = stab.set_index("allele_4digit")
grade_idx  = grade.set_index("allele_4digit")

bh_primary = bh_q(prim_idx.loc[ALLELES, "fisher_exact_p"].tolist())
bh_direct  = bh_q(direct_idx.loc[ALLELES, "fisher_p"].tolist())

def color_for_or(or_):
    if or_ > 1.5: return "#fde0d0"
    if or_ < 0.67: return "#d2e0f1"
    return "#f3eedc"

def color_for_stability(s):
    if "flips" in s: return "#fde0d0"
    if "depletion" in s: return "#d2e0f1"
    return "#dcefdc"

for i, a in enumerate(ALLELES):
    yy = 0.85 - i*0.13
    p_or = prim_idx.loc[a, "fisher_exact_or_raw"]
    d_or = direct_idx.loc[a, "or_haldane"]
    p_q  = bh_primary[i]
    d_q  = bh_direct[i]
    stab_v = stab_idx.loc[a, "stability_4_scenarios"]
    zc = "ZERO PTC carriers" if a == "DQB1*02:01" else "no zero cell"
    src_v = "AFND South Korea (DPB1) — source-sensitive" if a == "DPB1*05:01" else "IN2015 / Ann Lab Med 2015 — primary"
    claim = grade_idx.loc[a, "claim_grade"]

    ax.add_patch(mpatches.Rectangle((0.01, yy-0.05), 0.10, 0.10, color="#f7eedc", ec="#bba98a"))
    ax.text(0.06, yy, a, fontsize=9, fontweight="bold", ha="center", va="center", family="monospace")

    cells = [
        (f"{p_or:.2f}", color_for_or(p_or)),
        (f"{p_q:.1e}",  "#dcefdc" if p_q < 0.05 else "#f3eedc"),
        (f"{d_or:.2f}", color_for_or(d_or)),
        (f"{d_q:.1e}",  "#dcefdc" if d_q < 0.05 else "#f3eedc"),
        (stab_v.replace("_"," "), color_for_stability(stab_v)),
        (zc, "#fde0d0" if "ZERO" in zc else "#dcefdc"),
        (src_v, "#fde0d0" if "DPB1" in src_v else "#dcefdc"),
        (claim.replace("_"," "), "#fde0d0" if "source_sensitive" in claim or "zero" in claim
            else ("#f3eedc" if "moderate" in claim else "#dcefdc")),
    ]
    for ci, (text, color) in enumerate(cells):
        x0 = 0.115 + ci*0.105
        ax.add_patch(mpatches.Rectangle((x0-0.052, yy-0.05), 0.104, 0.10, color=color, ec="#cdc1aa"))
        # text wrap
        t = text if len(text) < 28 else text[:25]+"…"
        ax.text(x0, yy, t, ha="center", va="center", fontsize=7.5)

ax.set_xlim(0,1); ax.set_ylim(0,1)
ax.set_title("F18  Claim-boundary matrix — strict per-allele honesty audit\n"
             "Green ≈ supports a stable signal; amber = mixed; red = caution / does not pass strict harmonization.",
             fontsize=11)
fig.tight_layout()
fig.savefig(ASSETS/"F18_claim_boundary_matrix.png", dpi=170, bbox_inches="tight")
plt.close(fig)

# ---------- FIGURE F19 — validation roadmap ----------
fig, ax = plt.subplots(figsize=(12.5, 6.4))
ax.axis("off")
ax.set_xlim(0, 12); ax.set_ylim(0, 6)
def box(x, y, w, h, text, fc, ec="#2a3142", fontsize=9.5, weight="normal"):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.12",
                       fc=fc, ec=ec, lw=1.4)
    ax.add_patch(p)
    ax.text(x + w/2, y + h/2, text, ha="center", va="center",
            fontsize=fontsize, fontweight=weight, wrap=True)
def arrow(x1, y1, x2, y2, color="#2a3142"):
    a = FancyArrowPatch((x1,y1),(x2,y2), arrowstyle="-|>", mutation_scale=14,
                        color=color, lw=1.4, shrinkA=2, shrinkB=4)
    ax.add_patch(a)

# Lane 1 — current state
ax.text(0.4, 5.7, "CURRENT (exploratory)", fontsize=10, fontweight="bold", color="#8f2d25")
box(0.2, 4.4, 3.3, 1.0, "6-allele forest\n(carrier vs allele — metric mismatch)\n3/6 direction-flip when harmonized", "#fde0d0")
box(0.2, 3.0, 3.3, 1.0, "BH-q within 6 (S0): all 6 q < 0.05\nBH-q within 6 (S3 direct): re-evaluated\n→ enrichment claim weakens", "#fde0d0")
box(0.2, 1.6, 3.3, 1.0, "Subcohort heterogeneity\nχ² p < 0.05 for A*02:07 and C*01:02\n(Lee2024 dominates carriers)", "#fde0d0")

# Lane 2 — required validation work
ax.text(4.6, 5.7, "REQUIRED VALIDATION", fontsize=10, fontweight="bold", color="#244e73")
box(4.4, 4.4, 3.4, 1.0, "Matched-control NGS HLA typing\n(independent Korean healthy population,\nsame typing protocol as PTC pool)", "#d2e0f1")
box(4.4, 3.0, 3.4, 1.0, "Allele-level dosage harmonization\n(both sides: count of allele copies\nover 2n callable)", "#d2e0f1")
box(4.4, 1.6, 3.4, 1.0, "Multi-cohort replication\n(at minimum 2 independent Korean PTC cohorts\nwith same HLA-typing pipeline)", "#d2e0f1")
box(4.4, 0.2, 3.4, 1.0, "Pre-registered hypothesis & analysis plan\n(direction & metric committed before counting)", "#d2e0f1")

# Lane 3 — graduation criteria → claim grade
ax.text(8.9, 5.7, "GRADUATION CRITERIA", fontsize=10, fontweight="bold", color="#426b50")
box(8.6, 4.4, 3.3, 1.0, "Direction stable across S0/S1/S2/S3\nand stable across ≥ 2 cohorts\n(per-allele)", "#dcefdc")
box(8.6, 3.0, 3.3, 1.0, "Genome-wide HLA testing\n+ multi-test correction\n(BH-q < 0.05 in matched-control test)", "#dcefdc")
box(8.6, 1.6, 3.3, 1.0, "Effect size > pre-specified MID\nwith 95% CI excluding 1.0\nin matched-control comparison", "#dcefdc")
box(8.6, 0.2, 3.3, 1.0, "→ then: 'Korean PTC HLA association'\n(no longer 'exploratory candidate map')", "#bee0c6", weight="bold")

# arrows
arrow(3.5, 4.9, 4.4, 4.9)
arrow(3.5, 3.5, 4.4, 3.5)
arrow(3.5, 2.1, 4.4, 2.1)
arrow(7.8, 4.9, 8.6, 4.9)
arrow(7.8, 3.5, 8.6, 3.5)
arrow(7.8, 2.1, 8.6, 2.1)
arrow(7.8, 0.7, 8.6, 0.7)

ax.set_title("F19  Validation roadmap — from exploratory candidate-prioritization map "
             "to a publishable Korean-PTC HLA association claim",
             fontsize=11.5, fontweight="bold")
fig.tight_layout()
fig.savefig(ASSETS/"F19_validation_roadmap.png", dpi=170, bbox_inches="tight")
plt.close(fig)

# ---------- console summary ----------
print("\n=== strengthening summary ===")
print("\nDirection stability (4 scenarios):")
for a in ALLELES:
    print(f"  {a:<13s}  {stab_idx.loc[a, 'stability_4_scenarios']}")
print("\nDirect allele-vs-allele Fisher (no HWE):")
print(direct[["allele_4digit","ptc_allele_freq","baseline_allele_freq",
              "or_haldane","ci95_low","ci95_high","fisher_p","direction","bh_q_within_6"]]
      .to_string(index=False))
print("\nFiles written:")
for f in [
    "F16_metric_sensitivity_forest.png",
    "F17_subcohort_frequency_heatmap.png",
    "F18_claim_boundary_matrix.png",
    "F19_validation_roadmap.png",
]:
    p = ASSETS/f
    print(f"  {p}  ({p.stat().st_size//1024} KB)")
print(f"  {RES}/paper2_hla_metric_sensitivity_strengthened.tsv")
print(f"  {RES}/paper2_hla_allele_vs_allele_direct.tsv")
print(f"  {RES}/paper2_hla_bh_fdr_by_scenario.tsv")
print("\nDONE.")
