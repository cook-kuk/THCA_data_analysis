#!/usr/bin/env python3
"""
Paper 2 HLA — DQB1*02:01 zero-cell QC audit.

Goal: determine whether DQB1*02:01 = 0 carriers in the Korean PTC pool (n=874)
is biologically plausible or a technical artifact.

Inputs (read-only):
  project/results/d4p1_panasian_meta/korean_PTC_pool_n908.tsv          (per-sample 4-digit calls)
  project/results/v17_korean/arcasHLA/K2_arcasHLA_genotypes.tsv         (K2 raw arcasHLA)
  project/results/v17_korean/arcasHLA_GSE286332/                       (per-sample arcasHLA outputs)
  project/results/v17_korean/arcasHLA/                                  (K2 *.genotype.json if any)
  project/results/d4p1_panasian_meta/GSE286332_arcasHLA_genotypes.tsv

Outputs:
  project/results/p2_pillar1_forest_v2/paper2_dqb1_zero_cell_qc_audit.tsv
  project/reports/2026_05_06_paper2_dqb1_zero_cell_qc_report.md
  project/papers_hub_2026_05_04/assets/paper2_hla/F20_DQB1_zero_cell_qc_matrix.png
  project/papers_hub_2026_05_04/assets/paper2_hla/F21_DQB1_locus_callability_barplot.png

Discipline: CPU-only, no new download, no claim outside exploratory framing.
No connection to AITD/Paper 9. No edit to Paper 1.
"""
from __future__ import annotations
import json
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
RES  = ROOT / "project/results/p2_pillar1_forest_v2"
ASSETS = ROOT / "project/papers_hub_2026_05_04/assets/paper2_hla"
ASSETS.mkdir(parents=True, exist_ok=True)

POOL_TSV = ROOT / "project/results/d4p1_panasian_meta/korean_PTC_pool_n908.tsv"
K2_RAW   = ROOT / "project/results/v17_korean/arcasHLA/K2_arcasHLA_genotypes.tsv"
GSE286_RAW = ROOT / "project/results/d4p1_panasian_meta/GSE286332_arcasHLA_genotypes.tsv"
ARCAS_K2  = ROOT / "project/results/v17_korean/arcasHLA"
ARCAS_286 = ROOT / "project/results/v17_korean/arcasHLA_GSE286332"

LOCI_NAMES = ["A","B","C","DRB1","DQB1","DPB1"]

# ---------- 1. load per-sample pool ----------
print(f"[load] {POOL_TSV}")
pool = pd.read_csv(POOL_TSV, sep="\t")
print(f"  rows={len(pool)}, cohorts={dict(Counter(pool['cohort']))}")

# treat empty/whitespace as NA
for c in pool.columns:
    if c == "run" or c == "cohort": continue
    pool[c] = pool[c].replace({"": pd.NA, "nan": pd.NA, "NaN": pd.NA, "None": pd.NA})

# ---------- 2. per-locus callability ----------
def callable_n(df, locus):
    a1 = df[f"{locus}_a1_4d"].notna()
    a2 = df[f"{locus}_a2_4d"].notna()
    return int((a1 & a2).sum()), int(a1.sum() | a2.sum()), int((~(a1 | a2)).sum())

per_locus_rows = []
for locus in LOCI_NAMES:
    both, any_, miss = callable_n(pool, locus)
    n_total = len(pool)
    per_locus_rows.append({
        "locus": locus,
        "n_total_samples": n_total,
        "callable_both_alleles": both,
        "callable_at_least_one": any_,
        "missing_both_alleles": miss,
        "callable_pct_both": both / n_total * 100,
    })
    print(f"  {locus:>5}  callable_both={both:>4}  any={any_:>4}  miss={miss:>4}  pct={both/n_total*100:5.1f}%")
locus_df = pd.DataFrame(per_locus_rows)
locus_df.to_csv(RES/"_paper2_dqb1_audit_locus_callability.tsv", sep="\t", index=False)

# ---------- 3. DQB1 detail: 4-digit allele distribution + 2-digit family + subcohort ----------
def all_alleles_long(df, locus):
    alleles = pd.concat([df[[f"{locus}_a1_4d","cohort","run"]].rename(columns={f"{locus}_a1_4d":"allele"}),
                         df[[f"{locus}_a2_4d","cohort","run"]].rename(columns={f"{locus}_a2_4d":"allele"})],
                         axis=0, ignore_index=True)
    alleles = alleles[alleles["allele"].notna()]
    return alleles

dqb1_long = all_alleles_long(pool, "DQB1")
print(f"\n[DQB1] total non-null alleles = {len(dqb1_long)}; expected (callable both × 2) = {(pool['DQB1_a1_4d'].notna() & pool['DQB1_a2_4d'].notna()).sum()*2}")

# 4-digit unique
dqb1_4d_counts = dqb1_long["allele"].value_counts()
print(f"[DQB1] unique 4-digit alleles = {len(dqb1_4d_counts)}")

# 2-digit family extraction
dqb1_long["allele_2d"] = dqb1_long["allele"].apply(lambda s: s.split(":")[0] if isinstance(s, str) else None)
dqb1_2d_counts = dqb1_long["allele_2d"].value_counts()
print(f"[DQB1] 2-digit family distribution:")
for fam, n in dqb1_2d_counts.head(15).items():
    print(f"   {fam:<14s} {n}")

# DQB1*02 family detail (split into *02:01 vs *02:02 vs other 02:xx)
dqb1_02_mask = dqb1_long["allele_2d"] == "DQB1*02"
dqb1_02 = dqb1_long.loc[dqb1_02_mask, "allele"].value_counts()
print(f"\n[DQB1*02 family detail]")
for a, n in dqb1_02.items():
    print(f"   {a:<14s} {n}")

# Per-subcohort DQB1*02:01 carrier counts
def carrier_count(df, locus, allele_4d):
    a1 = df[f"{locus}_a1_4d"] == allele_4d
    a2 = df[f"{locus}_a2_4d"] == allele_4d
    return int((a1 | a2).sum())

cohorts = sorted(pool["cohort"].dropna().unique().tolist())
subcohort_rows = []
for c in cohorts:
    sub = pool[pool["cohort"] == c]
    n = len(sub)
    callable_both = int((sub["DQB1_a1_4d"].notna() & sub["DQB1_a2_4d"].notna()).sum())
    cc_0201 = carrier_count(sub, "DQB1", "DQB1*02:01")
    cc_0202 = carrier_count(sub, "DQB1", "DQB1*02:02")
    cc_02_any = int((sub["DQB1_a1_4d"].astype(str).str.startswith("DQB1*02") |
                     sub["DQB1_a2_4d"].astype(str).str.startswith("DQB1*02")).sum())
    cc_0501 = carrier_count(sub, "DQB1", "DQB1*05:01")
    cc_0301 = carrier_count(sub, "DQB1", "DQB1*03:01")
    subcohort_rows.append({
        "cohort": c,
        "n_total": n,
        "DQB1_callable_both": callable_both,
        "DQB1_callable_pct": callable_both/n*100 if n else 0,
        "carriers_DQB1*02:01": cc_0201,
        "carriers_DQB1*02:02": cc_0202,
        "carriers_DQB1*02_any": cc_02_any,
        "carriers_DQB1*05:01": cc_0501,
        "carriers_DQB1*03:01": cc_0301,
    })
sub_df = pd.DataFrame(subcohort_rows)
sub_df.to_csv(RES/"_paper2_dqb1_audit_subcohort.tsv", sep="\t", index=False)
print("\n[DQB1 subcohort detail]")
print(sub_df.to_string(index=False))

# ---------- 4. allele diversity per locus ----------
diversity_rows = []
for locus in LOCI_NAMES:
    long = all_alleles_long(pool, locus)
    n_unique = long["allele"].nunique()
    top1 = long["allele"].value_counts().head(1)
    top1_share = (top1.iloc[0] / len(long)) if len(long) else 0
    long["allele_2d"] = long["allele"].apply(lambda s: s.split(":")[0] if isinstance(s, str) else None)
    n_unique_2d = long["allele_2d"].nunique()
    diversity_rows.append({
        "locus": locus,
        "n_observations": len(long),
        "n_unique_4digit": n_unique,
        "n_unique_2digit_family": n_unique_2d,
        "top1_allele_share_pct": top1_share*100,
    })
div_df = pd.DataFrame(diversity_rows)
div_df.to_csv(RES/"_paper2_dqb1_audit_diversity.tsv", sep="\t", index=False)
print("\n[allele diversity per locus]")
print(div_df.to_string(index=False))

# ---------- 5. DQB1 detailed flagging ----------
# Detect "rare/impossible" allele patterns: 4-digit with very high digit code (e.g. *03:98 - flagged)
def numeric_subtype(allele):
    # parse "DQB1*03:98" -> 98
    try:
        sub = allele.split(":")[1].split(":")[0]
        return int(sub)
    except Exception:
        return None

dqb1_long["subtype_num"] = dqb1_long["allele"].apply(numeric_subtype)
high_subtype = dqb1_long[dqb1_long["subtype_num"].fillna(0) >= 50]
print(f"\n[DQB1] subtype number >= 50 (rare/impossible flag): {len(high_subtype)} observations across {high_subtype['allele'].nunique()} unique alleles")
print(high_subtype["allele"].value_counts().head(10))

# ---------- 6. arcasHLA per-sample QC inspection (K2 SRR/ERR + GSE286332 SRR) ----------
def scan_arcas_dir(dirpath):
    geno_jsons = list(dirpath.glob("**/*.genotype.json"))
    genes_jsons = list(dirpath.glob("**/*.genes.json"))
    sample_summary = []
    for gj in geno_jsons:
        sample = gj.name.replace(".genotype.json","")
        try:
            with open(gj) as fh: data = json.load(fh)
        except Exception:
            sample_summary.append({"sample":sample, "loci_present":[], "DQB1_in_genotype":False, "issue":"unreadable"})
            continue
        loci = list(data.keys())
        # arcasHLA stores list of length 2 per gene (or empty / ["",""])
        DQB1 = data.get("DQB1", None)
        sample_summary.append({
            "sample": sample,
            "loci_present": ",".join(loci),
            "DQB1_in_genotype": DQB1 is not None and len(DQB1) == 2 and all(DQB1),
            "DQB1_call": str(DQB1) if DQB1 else "",
        })
    return geno_jsons, genes_jsons, sample_summary

arcas_286_g, arcas_286_genes, arcas_286_summary = scan_arcas_dir(ARCAS_286)
print(f"\n[arcasHLA GSE286332] genotype.json files = {len(arcas_286_g)}, genes.json files = {len(arcas_286_genes)}")
arcas_286_df = pd.DataFrame(arcas_286_summary) if arcas_286_summary else pd.DataFrame()
if len(arcas_286_df):
    pos = int(arcas_286_df["DQB1_in_genotype"].sum())
    print(f"  DQB1 in genotype.json: {pos}/{len(arcas_286_df)}")

# Look into GSE286332 raw genotype TSV to see DQB1 distribution among the 18 samples
print(f"\n[GSE286332 raw genotypes — DQB1 distribution]")
gse286_raw = pd.read_csv(GSE286_RAW, sep="\t")
print(f"  rows={len(gse286_raw)}, columns include DQB1_a1_4d? {'DQB1_a1_4d' in gse286_raw.columns}")
gse286_dqb1_a1 = gse286_raw["DQB1_a1_4d"].fillna("(missing)").value_counts()
gse286_dqb1_a2 = gse286_raw["DQB1_a2_4d"].fillna("(missing)").value_counts()
print(f"  a1 distribution:\n{gse286_dqb1_a1.to_string()}")
print(f"  a2 distribution:\n{gse286_dqb1_a2.to_string()}")

# K2 raw genotypes for DQB1
print(f"\n[K2 raw genotypes — DQB1 distribution]")
k2_raw = pd.read_csv(K2_RAW, sep="\t")
print(f"  rows={len(k2_raw)}")
def k2_dqb1_status(row):
    a1 = str(row.get("DQB1_a1_4digit","")).strip()
    a2 = str(row.get("DQB1_a2_4digit","")).strip()
    if not a1 and not a2: return "both_missing"
    if not a1 or not a2: return "one_missing"
    return "both_called"
k2_raw["DQB1_status"] = k2_raw.apply(k2_dqb1_status, axis=1)
print("  K2 DQB1 call status:")
print(k2_raw["DQB1_status"].value_counts().to_string())

# how does K2 missing DQB1 relate to other loci on same sample?
k2_partial = k2_raw[k2_raw["DQB1_status"] != "both_called"]
print(f"\n  K2 samples with DQB1 missing/partial: {len(k2_partial)}")
if len(k2_partial):
    other_loci_status = []
    for _, r in k2_partial.iterrows():
        ok = []
        for L in ["A","B","C","DRB1","DPB1"]:
            v1 = str(r.get(f"{L}_a1_4digit","")).strip()
            v2 = str(r.get(f"{L}_a2_4digit","")).strip()
            if v1 and v2: ok.append(L)
        other_loci_status.append({"run": r["run"], "DQB1_status": r["DQB1_status"],
                                   "other_loci_called": ",".join(ok),
                                   "other_loci_called_n": len(ok)})
    pd.DataFrame(other_loci_status).to_csv(RES/"_paper2_dqb1_audit_K2_missing_samples.tsv", sep="\t", index=False)
    print(f"  → saved to _paper2_dqb1_audit_K2_missing_samples.tsv")

# ---------- 7. consolidated audit TSV ----------
audit_rows = []
for locus in LOCI_NAMES:
    rec = {"section":"locus_callability", "metric":locus,
           "value":locus_df[locus_df.locus == locus]["callable_pct_both"].iloc[0],
           "note":f"callable both alleles = {locus_df[locus_df.locus == locus]['callable_both_alleles'].iloc[0]} / {len(pool)}"}
    audit_rows.append(rec)

# DQB1 specific findings
audit_rows.append({"section":"DQB1_specific","metric":"unique_4digit_alleles","value":len(dqb1_4d_counts),
                   "note":"distinct 4-digit DQB1 alleles observed in callable samples"})
audit_rows.append({"section":"DQB1_specific","metric":"unique_2digit_families","value":dqb1_long["allele_2d"].nunique(),
                   "note":"distinct 2-digit DQB1 families observed"})
audit_rows.append({"section":"DQB1_specific","metric":"DQB1*02_family_total_observations",
                   "value":int(dqb1_02.sum()),
                   "note":f"sum across {len(dqb1_02)} 4-digit subtypes"})
for a, n in dqb1_02.items():
    audit_rows.append({"section":"DQB1_specific","metric":f"observations_{a}","value":n,
                       "note":"observed allele count in pool"})
audit_rows.append({"section":"DQB1_specific","metric":"observations_DQB1*02:01",
                   "value":int(dqb1_02.get("DQB1*02:01", 0)),
                   "note":"target allele observed count"})
audit_rows.append({"section":"DQB1_specific","metric":"high_subtype_number_count",
                   "value":int(len(high_subtype)),
                   "note":"observations with subtype number ≥ 50 (potential rare/typing-suspect)"})

# subcohort DQB1*02:01
for _, r in sub_df.iterrows():
    audit_rows.append({"section":"DQB1_subcohort","metric":f"{r['cohort']}_DQB1*02:01_carriers",
                       "value":int(r["carriers_DQB1*02:01"]),
                       "note":f"n_callable_both={int(r['DQB1_callable_both'])} of n_total={int(r['n_total'])} ({r['DQB1_callable_pct']:.1f}%)"})
    audit_rows.append({"section":"DQB1_subcohort","metric":f"{r['cohort']}_DQB1*02_any_carriers",
                       "value":int(r["carriers_DQB1*02_any"]),
                       "note":"any DQB1*02:xx carrier (regardless of subtype)"})

# diversity comparison
for _, r in div_df.iterrows():
    audit_rows.append({"section":"locus_diversity","metric":f"{r['locus']}_unique_4digit",
                       "value":int(r["n_unique_4digit"]),
                       "note":f"top-1 share = {r['top1_allele_share_pct']:.1f}%"})

audit_df = pd.DataFrame(audit_rows)
audit_df.to_csv(RES/"paper2_dqb1_zero_cell_qc_audit.tsv", sep="\t", index=False)
print(f"\n[save] {RES/'paper2_dqb1_zero_cell_qc_audit.tsv'}")

# ---------- FIGURE F20 — DQB1 zero-cell QC matrix (multipanel) ----------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6),
                         gridspec_kw=dict(width_ratios=[1.1, 1.4, 1.5]))
# panel A: DQB1*02 family breakdown
axA = axes[0]
fams = dqb1_02.head(8).reset_index()
fams.columns = ["allele","count"]
colors = ["#8f2d25" if a == "DQB1*02:01" else ("#b58534" if a.startswith("DQB1*02") else "#244e73")
          for a in fams["allele"]]
axA.barh(range(len(fams))[::-1], fams["count"].values, color=colors)
axA.set_yticks(range(len(fams))[::-1])
axA.set_yticklabels(fams["allele"], fontsize=10, family="monospace")
axA.set_xlabel("Observations (allele copies)")
axA.set_title(f"A. DQB1*02 family detail\n(target DQB1*02:01 highlighted; n={int(dqb1_02.sum())} total *02 obs)",
              fontsize=10.5)
for i, (a, c) in enumerate(zip(fams["allele"], fams["count"])):
    axA.text(c + 1, len(fams)-1-i, str(int(c)), va="center", fontsize=9)

# panel B: subcohort DQB1*02:01 vs *02_any vs *05:01 vs *03:01
axB = axes[1]
cohorts_order = sub_df["cohort"].tolist()
x = np.arange(len(cohorts_order))
w = 0.2
axB.bar(x - 1.5*w, sub_df["carriers_DQB1*02:01"].values, w, label="DQB1*02:01 (target)", color="#8f2d25")
axB.bar(x - 0.5*w, sub_df["carriers_DQB1*02_any"].values, w, label="DQB1*02:xx (any)", color="#b58534")
axB.bar(x + 0.5*w, sub_df["carriers_DQB1*05:01"].values, w, label="DQB1*05:01", color="#244e73")
axB.bar(x + 1.5*w, sub_df["carriers_DQB1*03:01"].values, w, label="DQB1*03:01", color="#426b50")
axB.set_xticks(x)
labs = [f"{c}\nn_callable={int(r)}\n/{int(t)} ({p:.0f}%)" for c, r, t, p in
        zip(sub_df["cohort"], sub_df["DQB1_callable_both"],
            sub_df["n_total"], sub_df["DQB1_callable_pct"])]
axB.set_xticklabels(labs, fontsize=8.5)
axB.set_ylabel("Carriers")
axB.set_title("B. DQB1 carrier counts per subcohort\n(any *02:xx vs target *02:01 vs reference alleles)",
              fontsize=10.5)
axB.legend(fontsize=8, loc="upper right", framealpha=0.95)

# panel C: per-locus callability bar (callable_pct_both)
axC = axes[2]
order = ["A","B","C","DRB1","DPB1","DQB1"]
ld = locus_df.set_index("locus").loc[order]
cols = ["#244e73"]*5 + ["#8f2d25"]
axC.bar(range(len(order)), ld["callable_pct_both"], color=cols)
axC.set_xticks(range(len(order))); axC.set_xticklabels(order, fontsize=10)
for i, (loc, pct, n) in enumerate(zip(order, ld["callable_pct_both"], ld["callable_both_alleles"])):
    axC.text(i, pct + 1, f"{pct:.1f}%\n({int(n)}/{len(pool)})", ha="center", fontsize=8.5)
axC.set_ylim(0, 110)
axC.set_ylabel("Callable both alleles (%)")
axC.set_title("C. Per-locus callability — DQB1 lowest by far\n(72.2% vs A/B/C ≥ 99.9%)",
              fontsize=10.5)
axC.axhline(100, color="black", lw=0.6, ls=":")

fig.suptitle("F20  DQB1*02:01 zero-cell QC matrix — Korean PTC pool n=874",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.95])
fig.savefig(ASSETS/"F20_DQB1_zero_cell_qc_matrix.png", dpi=170, bbox_inches="tight")
plt.close(fig)
print(f"[save] {ASSETS/'F20_DQB1_zero_cell_qc_matrix.png'}")

# ---------- FIGURE F21 — locus callability bar (with diversity overlay) ----------
fig, ax1 = plt.subplots(figsize=(10.5, 5.0))
order = ["A","B","C","DRB1","DPB1","DQB1"]
ld = locus_df.set_index("locus").loc[order]
dd = div_df.set_index("locus").loc[order]
cols = ["#244e73"]*5 + ["#8f2d25"]
xs = np.arange(len(order))
b1 = ax1.bar(xs, ld["callable_pct_both"], color=cols, label="callable both (%)")
ax1.set_xticks(xs); ax1.set_xticklabels(order, fontsize=11)
ax1.set_ylabel("Callable both alleles (%)", color="#244e73", fontsize=11)
ax1.set_ylim(0, 115)
for i, (pct, n) in enumerate(zip(ld["callable_pct_both"], ld["callable_both_alleles"])):
    ax1.text(i, pct + 2, f"{pct:.1f}%\nn={int(n)}", ha="center", fontsize=9, color="#1a1a1a")

ax2 = ax1.twinx()
ax2.plot(xs, dd["n_unique_4digit"], "o-", color="#b58534", lw=2, label="n unique 4-digit alleles")
ax2.set_ylabel("Unique 4-digit alleles observed", color="#b58534", fontsize=11)
for i, n in enumerate(dd["n_unique_4digit"]):
    ax2.text(i, n + 1.5, str(int(n)), ha="center", color="#7a5018", fontsize=9, fontweight="bold")

# annotate concern below the DQB1 bar
ax1.annotate(
    "DQB1: lowest call rate +\nmiddle-tier diversity →\nzero-cell signal needs technical audit",
    xy=(5, ld['callable_pct_both'].iloc[5]), xytext=(3.0, 35),
    fontsize=9, color="#8f2d25",
    arrowprops=dict(arrowstyle="->", color="#8f2d25", lw=1.2),
    bbox=dict(boxstyle="round,pad=0.4", fc="#fff0ed", ec="#8f2d25"))

ax1.set_title("F21  Per-locus callability + allele diversity — Korean PTC pool (n=874)\n"
              "DQB1 calls miss ~28% of samples; A/B/C ≥ 99.9%; DPB1 ~89.1%; DRB1 ~98.5%.",
              fontsize=11.5, fontweight="bold")
fig.tight_layout()
fig.savefig(ASSETS/"F21_DQB1_locus_callability_barplot.png", dpi=170, bbox_inches="tight")
plt.close(fig)
print(f"[save] {ASSETS/'F21_DQB1_locus_callability_barplot.png'}")

# ---------- FINAL VERDICT ----------
n_callable = int(locus_df.loc[locus_df.locus=="DQB1","callable_both_alleles"].iloc[0])
n_total    = int(locus_df.loc[locus_df.locus=="DQB1","n_total_samples"].iloc[0])
n_unique_dqb1 = int(div_df.loc[div_df.locus=="DQB1","n_unique_4digit"].iloc[0])
n_dqb1_02_obs = int(dqb1_02.sum())
n_target_obs = int(dqb1_02.get("DQB1*02:01", 0))
n_target_carriers = sum(sub_df["carriers_DQB1*02:01"].values)

print("\n=== DQB1*02:01 ZERO-CELL QC VERDICT (intermediate) ===")
print(f"DQB1 callable: {n_callable}/{n_total} ({n_callable/n_total*100:.1f}%) — vs A/B/C ≥ 99.9%")
print(f"DQB1 unique 4-digit alleles observed: {n_unique_dqb1}")
print(f"DQB1*02 family observations: {n_dqb1_02_obs}")
print(f"DQB1*02:01 observations: {n_target_obs}")
print(f"DQB1*02:01 carriers: {n_target_carriers}")
print(f"DQB1*02:xx (any) carriers per subcohort:")
for _, r in sub_df.iterrows():
    print(f"   {r['cohort']:>14s}  any-*02 = {int(r['carriers_DQB1*02_any']):>3}  target-*02:01 = {int(r['carriers_DQB1*02:01']):>3}  callable={int(r['DQB1_callable_both'])}/{int(r['n_total'])}")
print("DONE.")
