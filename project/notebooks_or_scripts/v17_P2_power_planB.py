#!/usr/bin/env python3
"""P2 — GSE286332 power analysis + 시나리오 2 Plan B decision-tree map.

Power: minimum detectable effect size at n=18 (9+9), n=36 (18+18), n=54, n=108.
Plan B: 4 scenarios mapped against Bundang outreach response.
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats

PROJ = Path("/opt/thyroid-dash/project")
RES = PROJ / "results/p2_power_planB"
RES.mkdir(parents=True, exist_ok=True)

# ============================================================
# 1. Power calc — Welch's t / Mann-Whitney equivalent
# ============================================================
print("=== 1. Power table ===")

def two_sided_power_d(d, n_per_group, alpha=0.05):
    """Approx power for Welch's t with effect size d, balanced n."""
    from scipy.stats import nct, t
    df = 2 * (n_per_group - 1)
    nc = d * np.sqrt(n_per_group / 2.0)
    crit = t.ppf(1 - alpha/2, df)
    pw = 1 - nct.cdf(crit, df, nc) + nct.cdf(-crit, df, nc)
    return pw

def min_detectable_d(n_per_group, alpha=0.05, target_pw=0.8):
    """Bisection: smallest d achieving target_pw."""
    lo, hi = 0.05, 5.0
    for _ in range(40):
        mid = (lo + hi) / 2
        if two_sided_power_d(mid, n_per_group, alpha) < target_pw:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2

power_rows = []
for n in [9, 18, 25, 50, 75, 100]:
    for d in [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]:
        pw = two_sided_power_d(d, n, alpha=0.05)
        power_rows.append(dict(n_per_group=n, cohen_d=d, power=round(pw, 3)))
power_df = pd.DataFrame(power_rows)
power_table = power_df.pivot(index="n_per_group", columns="cohen_d", values="power")
print("Power table (rows=n/group, cols=Cohen d, alpha=0.05):")
print(power_table.to_string())
power_table.to_csv(RES / "power_table.tsv", sep="\t")

mdd = []
for n in [9, 12, 18, 25, 50, 75, 100]:
    for pw in [0.8, 0.9]:
        d = min_detectable_d(n, alpha=0.05, target_pw=pw)
        mdd.append(dict(n_per_group=n, target_power=pw, min_detectable_d=round(d, 3)))
mdd_df = pd.DataFrame(mdd)
print("\nMinimum detectable Cohen d:")
print(mdd_df.to_string(index=False))
mdd_df.to_csv(RES / "min_detectable_d.tsv", sep="\t", index=False)

# ============================================================
# 2. Observed effect sizes from P3 (GSE286332)
# ============================================================
print("\n=== 2. Observed effects from P3 (GSE286332 9+9) ===")
p3 = json.loads((PROJ / "results/p3_gse286332/P3_summary.json").read_text())
obs = {
    "8-gene RAI score (panel z-mean)": abs(p3["panel_8gene_compare"]["cohen_d"]),
    "HLA-I module": abs(p3["HLA_I_cohen_d"]),
    "HLA-II module": abs(p3["HLA_II_cohen_d"]),
    "P(DM1) shift": None,  # not direct d; reported MW p
}
print("Observed |Cohen d| at n=9+9:")
for k, v in obs.items():
    if v is not None:
        pw = two_sided_power_d(v, 9, alpha=0.05)
        print(f"  {k}: |d|={v:.2f}, achieved power={pw:.3f}")

# Power surplus for replication at n=50+50
print("\nIf Bundang adds n=50 PTC + 50 PTC+HT (matched), achievable power for these effects:")
for k, v in obs.items():
    if v is not None:
        for n_repl in [25, 50, 100]:
            pw = two_sided_power_d(v / 2, n_repl, alpha=0.05)  # half effect (conservative)
            print(f"  {k} | half-effect d={v/2:.2f} | n={n_repl}: power={pw:.3f}")

# ============================================================
# 3. Plan B decision tree (Bundang outreach)
# ============================================================
print("\n=== 3. Plan B map (Bundang outreach scenarios) ===")
plan_b = {
    "scenario_A_strong_GO": {
        "trigger": "Bundang has Graves' n>50 AND PTC+HT n>30 with RNA-seq/SNP",
        "timeline_months": "6-9",
        "venue": ["Cell Rep Med (IF 14)", "JCI Insight (IF 8)", "Nat Commun reach"],
        "deliverables": [
            "GSE286332 discovery + Bundang validation cohort",
            "Pan-Asian HLA meta (Han + K2 + Lee + Bundang)",
            "PTC+HT molecular axis paper",
        ],
        "must_have": ["Bundang RNA-seq access agreement", "Bundang HLA imputation (cookHLA on SNP) or arcasHLA on RNA-seq"],
        "risks": ["IRB delay", "metadata richness (TRAb/anti-TPO/anti-Tg titers needed)"],
    },
    "scenario_B_moderate_GO": {
        "trigger": "Bundang has PTC cohort (any n) but NO Graves'/HT",
        "timeline_months": "5-7",
        "venue": ["JCI Insight (IF 8)", "Genome Medicine (IF 11)"],
        "deliverables": [
            "K2+Lee+TCGA Hashimoto-like meta",
            "GSE286332 as orthogonal molecular validation",
            "8-gene panel external validation on Bundang PTC",
        ],
        "must_have": ["Bundang PTC RNA-seq or NanoString 8-gene"],
        "drop": ["Graves' arm — defer to Phase 2 paper"],
    },
    "scenario_C_method_paper": {
        "trigger": "Bundang non-responsive 6+ weeks OR all access denied",
        "timeline_months": "3-4",
        "venue": ["Brief Bioinform (IF 7)", "Bioinformatics (IF 4-7)"],
        "deliverables": [
            "DIAL audit framework + agent forensic case study (8-gene Task B)",
            "arcasHLA + cookHLA cross-platform pipeline",
            "GSE286332 as included case",
        ],
        "drop": ["Bundang validation"],
    },
    "scenario_D_minimum_viable": {
        "trigger": "Manuscript freeze required by date (paper conflict, ARIA pull-back, etc.)",
        "timeline_months": "2-3",
        "venue": ["Sci Rep (IF 4)", "Endocrine-Related Cancer (IF 5)"],
        "deliverables": [
            "Minimal: 8-gene panel + DM1/DM2 + K2 validation only",
            "GSE286332 as small Suppl figure",
        ],
        "drop": ["HLA arm", "Hashimoto axis"],
    },
}

(RES / "plan_B_map.json").write_text(json.dumps(plan_b, indent=2))
for k, v in plan_b.items():
    print(f"\n--- {k} ---")
    print(f"  Trigger: {v['trigger']}")
    print(f"  Timeline: {v['timeline_months']} months")
    print(f"  Venue: {', '.join(v['venue'])}")

# ============================================================
# 4. Bundang outreach query template (paper-ready)
# ============================================================
print("\n=== 4. Bundang outreach query template ===")
template = """[분당 outreach query — Day 2 PM 추가 안건]

분당서울대 cohort 와 관련해 다음 4가지 항목을 부탁드립니다:

1. PTC + Hashimoto's thyroiditis 동반 환자
   - 가능 sample 수
   - modality (RNA-seq / NanoString / WES / SNP genotype)
   - clinical metadata: anti-TPO, anti-Tg titer, TSH, T4, BRAFV600E status

2. Graves' disease (양성 자가면역 hyperthyroidism)
   - 가능 sample 수 (어떤 modality라도)
   - HLA typing 또는 SNP imputation 가능 여부
   - clinical metadata: TRAb titer, TSH, FT3/FT4, 치료력 (methimazole, RAI ablation, surgery)

3. 협력 가능한 modality
   - 만약 RNA-seq가 어렵다면 NanoString 8-gene panel만 측정해도 매우 유효합니다 (validation cohort)
   - 또는 SNP genotype만 있어도 cookHLA로 4-digit HLA imputation 가능

4. Timeline
   - 6-9 month 내 협력 가능 여부 (paper 제출 기한)
   - IRB 절차에 필요한 기간

--- 우리 측 제공 ---
- GSE286332 (Korean PTC+HT, n=18) discovery 완료: 10,380 DEGs, IFN-γ + HLA-II + dedifferentiation 강력 신호
- K2 (PRJEB11591, n=260) + Lee 2024 (GSE213647, n=630) arcasHLA HLA imputation 완료
- 8-gene panel 외부 validation 가능 pipeline
- cookHLA + arcasHLA 둘 다 본인 first-author 도구

--- 기대 outcome ---
- 분당 시료 + 우리 데이터 = Pan-Korean autoimmune-PTC-axis paper 가능
- Cell Rep Med / JCI Insight 수준 reach 가능 (Phase 2)
"""
(RES / "bundang_outreach_query.txt").write_text(template)
print(template)

# ============================================================
# 5. Summary JSON
# ============================================================
summary = {
    "power_at_n9_d_minus1.6_alpha05": round(two_sided_power_d(1.6, 9), 3),
    "power_at_n9_d_3.6_alpha05": round(two_sided_power_d(3.6, 9), 3),
    "power_at_n50_d_0.5_alpha05": round(two_sided_power_d(0.5, 50), 3),
    "power_at_n50_d_0.8_alpha05": round(two_sided_power_d(0.8, 50), 3),
    "min_detectable_d_n9_pw80": round(min_detectable_d(9, target_pw=0.8), 3),
    "min_detectable_d_n18_pw80": round(min_detectable_d(18, target_pw=0.8), 3),
    "min_detectable_d_n50_pw80": round(min_detectable_d(50, target_pw=0.8), 3),
    "observed_effects_GSE286332": {
        "panel_8gene": p3["panel_8gene_compare"]["cohen_d"],
        "hla_I": p3["HLA_I_cohen_d"],
        "hla_II": p3["HLA_II_cohen_d"],
    },
    "plan_B": plan_b,
}
(RES / "P2_summary.json").write_text(json.dumps(summary, indent=2, default=str))
print(f"\n✓ All outputs saved to {RES}")
