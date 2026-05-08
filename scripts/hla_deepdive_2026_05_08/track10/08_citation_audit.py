#!/usr/bin/env python3
"""
Track 10 — Step 8: Citation correction audit.

Cross-checks the prior Track 1 v2 source rows against the resolved Track 10
PMIDs/DOIs/abstracts and flags any inconsistency.
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT  = ROOT / "project/results/hla_deepdive_2026_05_08/track10_korean_lit"
TBL  = OUT / "tables"
T1   = ROOT / "project/results/hla_deepdive_2026_05_08/track1_paper4_gd_panasian/tables"

audits = []

# Audit 1: Track 1 v2 had Park 2005 DPB1*05:01 row.
audits.append(dict(
    audit_id="A1_park2005_DPB1",
    track="Track 1 v2",
    finding="Park 2005 PMID 15993720 (Park MH et al, Hum Immunol 66:741-7) reports HLA-DRB1 and HLA-DQB1 ONLY. The abstract verbatim: 'analyzed HLA-DRB1 and -DQB1 associations with GD in 198 Korean patients compared with 200 healthy controls'. No DPB1 typing was performed.",
    track1_v2_row="Park 2005 row for DPB1*05:01 with OR=2.05, p=0.003, n_case=88, n_ctrl=104",
    severity="HIGH — citation error",
    action="DROPPED in Track 10 v3 forest. Likely the source of those numbers was a different paper (perhaps Inoue 1992 / Onuma 1994 / a Korean DPB1 paper that was not Park 2005). Future re-curation needed.",
    status="resolved",
))

# Audit 2: Track 1 v2 Park 2005 cohort size mismatch
audits.append(dict(
    audit_id="A2_park2005_cohort",
    track="Track 1 v2",
    finding="Track 1 v2 row had n_case=88, n_control=104 for Park 2005. Park 2005 abstract states n_case=198, n_control=200.",
    track1_v2_row="Park 2005 (n_case=88, n_control=104)",
    severity="MEDIUM — cohort mis-sized",
    action="Track 10 v3 uses correct n_case=198, n_control=200 for Park 2005 DRB1 rows.",
    status="resolved",
))

# Audit 3: Shin 2019 DPB1*05:01 Pc-derived SE in Track 1 v2
audits.append(dict(
    audit_id="A3_shin2019_se_source",
    track="Track 1 v2",
    finding="Shin 2019 rows in v2 derived SE from corrected p (Pc) instead of from the reported 95% CI. This biases SE estimates because Pc reflects multiple-comparison correction, not the underlying test precision.",
    track1_v2_row="Shin 2019 rows for B*46:01, C*01:02, DPB1*05:01 (se_source='corrected_p_approx')",
    severity="MEDIUM — SE biased",
    action="Track 10 v3 replaces all Shin 2019 rows with full S1-Table-derived 95% CIs (raw p values, accurate Wald SE).",
    status="resolved",
))

# Audit 4: Manifest hint for Baek 2021 incorrect
audits.append(dict(
    audit_id="A4_baek2021_locus",
    track="data_registry manifest",
    finding="The data registry hint described Baek 2021 as 'Korean class-II NGS reference' implying DPB1/DQB1/DPA1/DQA1 coverage. The actual paper (HLA 97:112-126, n=26,202) is HLA-A, -B, -DRB1 only.",
    track1_v2_row="N/A (manifest, not yet ingested)",
    severity="LOW — manifest annotation only",
    action="Manifest text should be updated to: 'Korean class-I + DRB1 NGS reference (HLA-A/B/DRB1, n=26,202)'.",
    status="open — flagged for manifest update",
))

# Audit 5: Cho 2011 abstract reports cw* and PCR-SSP (low-resolution); subgroup ns differ
audits.append(dict(
    audit_id="A5_cho2011_resolution",
    track="Track 10 (new)",
    finding="Cho 2011 PMID 21952423 uses PCR-SSP at 2-digit allele family resolution (HLA-A*02 etc), not 4-digit. Cohort: n_AITD=73 (n_HD=32 + n_GD=41) vs n_CTRL=159. The abstract reports DIRECTIONAL associations only without ORs/CIs, so we recorded 18 directional rows (no numeric meta inputs).",
    track1_v2_row="N/A (Cho 2011 was missing from Track 1 v2 entirely)",
    severity="INFO — resolution mismatch",
    action="Cho 2011 cannot be added to the 4-digit pan-Asian meta as-is. Manual table extraction from the paywalled PDF is required for inclusion. Documented for future re-curation.",
    status="open — extract manually if PDF acquired",
))

# Audit 6: Jang 2011 6-digit -> 4-digit collapse
audits.append(dict(
    audit_id="A6_jang2011_resolution",
    track="Track 10 (new)",
    finding="Jang 2011 PMID 21062236 reports DRB1 at 6-digit (e.g. DRB1*030101). Collapsing to 4-digit (DRB1*03:01) loses sub-allele information but enables alignment with Track 1 v2 4-digit anchors. Frequencies are per-chromosome (allele frequencies), not carrier frequencies.",
    track1_v2_row="N/A (Jang 2011 was missing from Track 1 v2 entirely)",
    severity="INFO — collapse documented",
    action="Track 10 v3 carries an allele_6d column and a 4-digit-collapsed allele column. Forest uses 4-digit. Future high-res forest can use the 6-digit version.",
    status="resolved",
))

df = pd.DataFrame(audits)
df.to_csv(TBL / "T12_citation_audit.tsv", sep="\t", index=False)
print(f"saved -> {TBL/'T12_citation_audit.tsv'}")
print()
print(df[["audit_id","severity","status"]].to_string(index=False))
