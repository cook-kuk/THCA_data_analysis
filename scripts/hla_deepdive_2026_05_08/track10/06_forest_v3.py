#!/usr/bin/env python3
"""
Track 10 — Step 6: Build forest_panasian_GD_v3 by merging Track 1 v2 source
rows with Track 10 Korean extractions.

Rules:
  - Drop the v2 Park 2005 DPB1*05:01 row (citation error — Park 2005 reports
    DR/DQ only, never DPB1).
  - Replace the v2 Shin 2019 rows for B*46:01, C*01:02, DPB1*05:01 with the
    new accurate full-CI-derived rows from S1 Table.
  - Add new Shin 2019 rows for DRB1*08:02, DRB1*15:01, DRB1*16:02, A*02:07,
    C*03:02, DQB1*03:02 (the focus alleles previously D-grade).
  - Add Park 2005 rows for DRB1*08:03, DRB1*16:02, DRB1*03:01 (overall stratum;
    we exclude the male-only stratum from main forest).
  - Add Jang 2011 rows for DRB1*03:01, DRB1*08:02, DRB1*14:03, DRB1*07:01,
    DRB1*13:02.

Then re-run DerSimonian-Laird random-effects meta on the focus alleles plus
the supplementary informative ones, save v3 forest plot in {png,pdf,svg}.
"""
from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
T1   = ROOT / "project/results/hla_deepdive_2026_05_08/track1_paper4_gd_panasian/tables"
OUT  = ROOT / "project/results/hla_deepdive_2026_05_08/track10_korean_lit"
TBL  = OUT / "tables"
PLT  = OUT / "plots"
PER  = PLT / "per_allele"
for d in (TBL, PLT, PER):
    d.mkdir(parents=True, exist_ok=True)

# ===== Load v2 anchor and master extraction =====
v2 = pd.read_csv(T1 / "T01_panasian_GD_source_rows_v2.tsv", sep="\t")
master = pd.read_csv(TBL / "T03_korean_lit_extraction_master.tsv", sep="\t")

# ===== Helpers =====
def ci_to_se(lo, hi):
    if not (np.isfinite(lo) and np.isfinite(hi)) or lo <= 0 or hi <= 0:
        return np.nan
    return (math.log(hi) - math.log(lo)) / (2 * 1.959964)

def or_p_to_se(or_val, p_val):
    if not (np.isfinite(or_val) and np.isfinite(p_val)) or p_val <= 0 or p_val >= 1:
        return np.nan
    z = stats.norm.isf(p_val / 2.0)
    if z == 0:
        return np.nan
    return abs(math.log(or_val)) / z

def ds_l(log_or, se):
    log_or = np.asarray(log_or, float); se = np.asarray(se, float)
    mask = np.isfinite(log_or) & np.isfinite(se) & (se > 0)
    log_or = log_or[mask]; se = se[mask]
    k = len(log_or)
    if k == 0:
        return None
    w_fe = 1.0 / (se ** 2)
    fe = (w_fe * log_or).sum() / w_fe.sum()
    Q = (w_fe * (log_or - fe) ** 2).sum()
    df = k - 1
    if df <= 0:
        tau2 = 0.0
    else:
        c = w_fe.sum() - (w_fe ** 2).sum() / w_fe.sum()
        tau2 = max(0.0, (Q - df) / c)
    w_re = 1.0 / (se ** 2 + tau2)
    pooled = (w_re * log_or).sum() / w_re.sum()
    pooled_se = math.sqrt(1.0 / w_re.sum())
    ci_lo = math.exp(pooled - 1.959964 * pooled_se)
    ci_hi = math.exp(pooled + 1.959964 * pooled_se)
    z = pooled / pooled_se
    p_random = 2 * (1 - stats.norm.cdf(abs(z)))
    I2 = max(0.0, (Q - df) / Q * 100) if Q > 0 else 0.0
    if df >= 1 and tau2 > 0:
        t_crit = stats.t.isf(0.025, df=df)
        pi_lo = math.exp(pooled - t_crit * math.sqrt(pooled_se ** 2 + tau2))
        pi_hi = math.exp(pooled + t_crit * math.sqrt(pooled_se ** 2 + tau2))
    else:
        pi_lo = pi_hi = math.nan
    return dict(k=k, pooled_or=math.exp(pooled), ci_lo=ci_lo, ci_hi=ci_hi,
                p_random=p_random, tau2=tau2, Q=Q, df=df, I2_pct=I2,
                pi_lo=pi_lo, pi_hi=pi_hi)

# ===== Build v3 source rows =====

# 1) Drop the citation-error Park 2005 DPB1*05:01 v2 row.
keep_v2 = ~((v2["source"] == "Park 2005") & (v2["allele_4d"] == "DPB1*05:01"))
v3 = v2[keep_v2].copy()

# 2) Drop the prior Shin 2019 rows that we'll replace with full-CI derived.
shin_to_replace = ["B*46:01", "C*01:02", "DPB1*05:01"]
v3 = v3[~((v3["source"] == "Shin 2019") & (v3["allele_4d"].isin(shin_to_replace)))].copy()

# 3) Build new Korean source rows from master (Controls_vs_GD only, focus + replaced).
NEW_ALLELES = {
    "Shin_2019": ["DPB1*05:01", "B*46:01", "DRB1*08:02", "DRB1*15:01", "DRB1*16:02",
                  "A*02:07", "C*03:02", "DQB1*03:02", "C*01:02"],
    "Park_2005": ["DRB1*08:03", "DRB1*16:02", "DRB1*03:01"],
    "Jang_2011": ["DRB1*03:01", "DRB1*08:02", "DRB1*14:03", "DRB1*07:01", "DRB1*13:02"],
}

new_rows = []
for paper_id, alleles in NEW_ALLELES.items():
    for a in alleles:
        if paper_id == "Shin_2019":
            sub = master[(master["paper_id"] == "Shin_2019") &
                          (master["comparison"] == "Controls_vs_GD") &
                          (master["allele"] == a)]
        else:
            sub = master[(master["paper_id"] == paper_id) &
                          (master["comparison"] == "Controls_vs_GD") &
                          (master["allele"] == a)]
        if len(sub) == 0:
            continue
        r = sub.iloc[0]
        OR = r["OR"]; lo = r["OR_lo"]; hi = r["OR_hi"]
        se = ci_to_se(lo, hi) if (np.isfinite(lo) and np.isfinite(hi)) else r.get("se_log_or", np.nan)
        if not np.isfinite(se) or se <= 0:
            continue
        log_or = math.log(OR)
        country = "Korea"
        ancestry = "Korean children" if paper_id == "Shin_2019" else "Korean adult"
        source_label = paper_id.replace("_", " ")
        nc_a = r.get("n_case", np.nan); nc_b = r.get("n_control", np.nan)
        new_rows.append(dict(
            source=source_label, country=country, ancestry=ancestry, allele_4d=a,
            or_value=OR, ci_lo_derived=math.exp(log_or - 1.959964*se),
            ci_hi_derived=math.exp(log_or + 1.959964*se),
            p_value=(float(r["p"]) if pd.notna(r["p"]) and isinstance(r["p"], (int, float)) and np.isfinite(r["p"]) else np.nan),
            n_case=nc_a, n_control=nc_b,
            log_or=log_or, se=se, se_source="reported_CI",
            weight_fe=np.nan, weight_fe_norm=np.nan,
            caveat=f"Track-10 extraction: {r['cohort']}; {r['typing']}; freq={r['frequency_type']}",
        ))
new_df = pd.DataFrame(new_rows)
v3 = pd.concat([v3, new_df], ignore_index=True, sort=False)

# Recompute weight columns
v3["weight_fe"] = 1.0 / v3["se"] ** 2
v3["weight_fe_norm"] = v3.groupby("allele_4d")["weight_fe"].transform(lambda x: x / x.sum())

v3 = v3.sort_values(["allele_4d", "source"])
v3.to_csv(TBL / "T04_panasian_GD_source_rows_v3.tsv", sep="\t", index=False)
print(f"saved -> {TBL/'T04_panasian_GD_source_rows_v3.tsv'}  ({len(v3)} rows)")
print(f"\nSource rows by allele (v3):")
print(v3.groupby("allele_4d")["source"].apply(lambda x: ";".join(x)).to_string())

# ===== Run meta v3 =====
FOCUS = ["DPB1*05:01", "B*46:01", "DRB1*08:02", "DRB1*15:01", "DRB1*16:02",
         "A*02:07", "C*03:02", "DQB1*03:02",
         "C*01:02", "DQB1*02:01", "DRB1*07:01",
         "DRB1*08:03", "DRB1*03:01", "DRB1*14:03", "DRB1*13:02"]

meta_rows = []
for a in FOCUS:
    sub = v3[v3["allele_4d"] == a]
    if len(sub) == 0:
        meta_rows.append(dict(allele=a, k=0, sources="(no data)", pooled_or=np.nan,
                              ci_lo=np.nan, ci_hi=np.nan, p_random=np.nan,
                              tau2=np.nan, Q=np.nan, df=np.nan, I2_pct=np.nan,
                              pi_lo=np.nan, pi_hi=np.nan, status="NO_DATA"))
        continue
    res = ds_l(sub["log_or"].values, sub["se"].values)
    meta_rows.append(dict(
        allele=a, k=res["k"], sources="; ".join(sub["source"].astype(str).tolist()),
        pooled_or=res["pooled_or"], ci_lo=res["ci_lo"], ci_hi=res["ci_hi"],
        p_random=res["p_random"], tau2=res["tau2"], Q=res["Q"], df=res["df"],
        I2_pct=res["I2_pct"], pi_lo=res["pi_lo"], pi_hi=res["pi_hi"],
        status=("single_source" if res["k"] == 1 else "random_effects_DL"),
    ))
meta_df = pd.DataFrame(meta_rows).sort_values("p_random")
meta_df.to_csv(TBL / "T05_panasian_GD_DL_random_effects_v3.tsv", sep="\t", index=False)
print(f"\nsaved -> {TBL/'T05_panasian_GD_DL_random_effects_v3.tsv'}")
print(meta_df[["allele","k","pooled_or","ci_lo","ci_hi","I2_pct","p_random"]].to_string(index=False))

# ===== v2 vs v3 diff =====
v2_meta = pd.read_csv(T1 / "T02_panasian_GD_DL_random_effects_v2.tsv", sep="\t")
diff_rows = []
for a in FOCUS:
    v2r = v2_meta[v2_meta["allele"] == a]
    v3r = meta_df[meta_df["allele"] == a]
    if len(v2r) == 0 or len(v3r) == 0:
        continue
    v2r = v2r.iloc[0]; v3r = v3r.iloc[0]
    diff_rows.append(dict(
        allele=a,
        k_v2=v2r["k"], k_v3=v3r["k"], delta_k=v3r["k"] - v2r["k"],
        OR_v2=v2r["pooled_or"], OR_v3=v3r["pooled_or"],
        delta_log_OR=math.log(v3r["pooled_or"] / v2r["pooled_or"])
            if (np.isfinite(v3r["pooled_or"]) and np.isfinite(v2r["pooled_or"]) and v3r["pooled_or"] > 0 and v2r["pooled_or"] > 0) else np.nan,
        ci_lo_v2=v2r["ci_lo"], ci_hi_v2=v2r["ci_hi"],
        ci_lo_v3=v3r["ci_lo"], ci_hi_v3=v3r["ci_hi"],
        I2_v2=v2r["I2_pct"], I2_v3=v3r["I2_pct"],
        delta_I2=v3r["I2_pct"] - v2r["I2_pct"],
        sources_v3=v3r["sources"],
    ))
diff_df = pd.DataFrame(diff_rows)
diff_df.to_csv(TBL / "T06_v2_vs_v3_meta_diff.tsv", sep="\t", index=False)
print(f"\nsaved -> {TBL/'T06_v2_vs_v3_meta_diff.tsv'}")
print("\nv2->v3 OR shift:")
print(diff_df[["allele","k_v2","k_v3","OR_v2","OR_v3","I2_v2","I2_v3"]].to_string(index=False))

# ===== Forest plot v3 =====
def forest_master(meta_df, out_stem):
    md = meta_df[(meta_df["k"] >= 1) & np.isfinite(meta_df["pooled_or"])].copy().reset_index(drop=True)
    md["log_or"] = np.log(md["pooled_or"])
    md["log_lo"] = np.log(md["ci_lo"])
    md["log_hi"] = np.log(md["ci_hi"])
    fig, ax = plt.subplots(figsize=(13, max(6, 0.55 * len(md) + 2)), dpi=150)
    yticks = np.arange(len(md))
    risk = md["pooled_or"] > 1
    ax.errorbar(md["log_or"], yticks,
                xerr=[md["log_or"] - md["log_lo"], md["log_hi"] - md["log_or"]],
                fmt="none", color="black", capsize=3, lw=1.2)
    for i, (lor, c) in enumerate(zip(md["log_or"], np.where(risk, "#cc3333", "#3366cc"))):
        ax.plot(lor, i, "s", color=c, markersize=11, markeredgecolor="black")
    for i, row in md.iterrows():
        if np.isfinite(row.get("pi_lo", np.nan)) and np.isfinite(row.get("pi_hi", np.nan)):
            ax.hlines(i, math.log(row["pi_lo"]), math.log(row["pi_hi"]),
                      color="grey", lw=0.8, ls=":")
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(yticks); ax.set_yticklabels(md["allele"], fontsize=11)
    ax.invert_yaxis()
    xt = [math.log(x) for x in [0.2, 0.5, 1, 2, 3, 5, 10, 25]]
    ax.set_xticks(xt); ax.set_xticklabels(["0.2","0.5","1","2","3","5","10","25"])
    ax.set_xlabel("Pooled OR (random-effects, DerSimonian-Laird)", fontsize=12)
    txt_x = math.log(35)
    ax.set_xlim(math.log(0.1), math.log(50))
    ax.set_title(
        "Pan-Asian Graves' disease — HLA susceptibility random-effects meta-analysis (v3)\n"
        "(autoimmune-only; Track 10 Korean literature integration)\n",
        fontsize=12,
    )
    for i, row in md.iterrows():
        s = f"OR {row['pooled_or']:.2f}  [{row['ci_lo']:.2f}, {row['ci_hi']:.2f}]   k={row['k']}   I²={row['I2_pct']:.0f}%"
        ax.text(txt_x, i, s, fontsize=8.5, va="center", family="monospace")
    plt.tight_layout()
    fig.savefig(f"{out_stem}.png", dpi=150, bbox_inches="tight")
    fig.savefig(f"{out_stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{out_stem}.svg", bbox_inches="tight")
    plt.close(fig)

forest_master(meta_df, PLT / "forest_panasian_GD_v3")
print(f"\nsaved -> {PLT/'forest_panasian_GD_v3.{png,pdf,svg}'}")

# ===== Per-allele forests for the 8 focus alleles =====
def forest_per_allele(allele, src_df, meta_row, out_path):
    sub = src_df[src_df["allele_4d"] == allele].sort_values(["country", "source"]).reset_index(drop=True)
    if len(sub) == 0:
        return
    fig, ax = plt.subplots(figsize=(11, max(3.5, 0.6 * len(sub) + 2.5)), dpi=150)
    yticks = np.arange(len(sub))
    log_or = sub["log_or"].astype(float).values
    se = sub["se"].astype(float).values
    log_lo = log_or - 1.959964 * se
    log_hi = log_or + 1.959964 * se
    weights = 1.0 / (se ** 2); weights = weights / weights.sum()
    sizes = 80 + 700 * weights
    colors = ["#cc3333" if (lo > 0) else ("#3366cc" if hi < 0 else "#888888")
              for lo, hi in zip(log_lo, log_hi)]
    ax.errorbar(log_or, yticks, xerr=[log_or - log_lo, log_hi - log_or],
                fmt="none", color="black", capsize=3, lw=1.0)
    for i, (lor, c, sz) in enumerate(zip(log_or, colors, sizes)):
        ax.scatter([lor], [i], s=sz, color=c, edgecolor="black", zorder=5)
    if meta_row is not None and np.isfinite(meta_row["pooled_or"]):
        diamond_y = len(sub) + 0.5
        cx = math.log(meta_row["pooled_or"])
        clo = math.log(meta_row["ci_lo"]); chi = math.log(meta_row["ci_hi"])
        diamond = plt.Polygon([(clo, diamond_y), (cx, diamond_y - 0.4),
                               (chi, diamond_y), (cx, diamond_y + 0.4)],
                              facecolor="black", edgecolor="black", zorder=10)
        ax.add_patch(diamond)
        ax.text(chi + 0.05, diamond_y,
                f"RE pooled OR={meta_row['pooled_or']:.2f}  [{meta_row['ci_lo']:.2f},{meta_row['ci_hi']:.2f}]  "
                f"I²={meta_row['I2_pct']:.0f}%  τ²={meta_row['tau2']:.3f}  k={meta_row['k']}",
                fontsize=9, va="center", family="monospace")
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(list(yticks) + [len(sub) + 0.5])
    ylabels = []
    for _, r in sub.iterrows():
        n_str = (f"n_case={int(r['n_case'])}/n_ctrl={int(r['n_control'])}"
                 if np.isfinite(r["n_case"]) and np.isfinite(r["n_control"]) else "n=NA")
        ylabels.append(f"{r['source']} ({r['country']}, {n_str})")
    ylabels.append("Pan-Asian RE pooled (v3)")
    ax.set_yticklabels(ylabels, fontsize=9)
    ax.invert_yaxis()
    xt = [math.log(x) for x in [0.2, 0.5, 1, 2, 3, 5, 10]]
    ax.set_xticks(xt); ax.set_xticklabels(["0.2","0.5","1","2","3","5","10"])
    ax.set_xlim(math.log(0.1), math.log(15))
    ax.set_xlabel("OR (95% CI), Graves' disease vs ancestry-matched control")
    ax.set_title(f"Pan-Asian GD HLA — {allele} (v3)\n(autoimmune Graves' disease; no cancer claim)", fontsize=12)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

cnt = 0
for a in FOCUS:
    sub = v3[v3["allele_4d"] == a]
    if len(sub) == 0:
        continue
    mr = meta_df[meta_df["allele"] == a]
    mr = mr.iloc[0] if len(mr) else None
    safe = a.replace("*", "_").replace(":", "_")
    forest_per_allele(a, v3, mr, PER / f"forest_{safe}_v3.png")
    cnt += 1
print(f"saved {cnt} per-allele v3 forest pages -> {PER}/")

# ===== Korean-only sub-meta v2 =====
korean = v3[v3["country"].astype(str).str.contains("Korea", case=False)]
korean_meta_rows = []
for a in FOCUS:
    sub = korean[korean["allele_4d"] == a]
    if len(sub) == 0:
        continue
    res = ds_l(sub["log_or"].values, sub["se"].values)
    korean_meta_rows.append(dict(
        allele=a, k=res["k"], sources="; ".join(sub["source"].astype(str).tolist()),
        pooled_or=res["pooled_or"], ci_lo=res["ci_lo"], ci_hi=res["ci_hi"],
        p_random=res["p_random"], I2_pct=res["I2_pct"],
        tau2=res["tau2"], Q=res["Q"], df=res["df"],
    ))
kor_meta = pd.DataFrame(korean_meta_rows).sort_values("p_random")
kor_meta.to_csv(TBL / "T07_korean_only_sub_meta_v2.tsv", sep="\t", index=False)
print(f"\nsaved -> {TBL/'T07_korean_only_sub_meta_v2.tsv'}")
print(kor_meta.to_string(index=False))

# Korean-only forest plot
forest_master(kor_meta, PLT / "forest_korean_only_v2")
print(f"saved -> {PLT/'forest_korean_only_v2.{png,pdf,svg}'}")
