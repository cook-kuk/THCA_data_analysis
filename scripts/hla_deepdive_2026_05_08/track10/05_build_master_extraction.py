#!/usr/bin/env python3
"""
Track 10 — Step 5: Unified korean_lit_extraction_master.tsv

Merge Shin 2019 (full S1 Table), Park 2005 (3 abstract-derived rows), Jang 2011
(5 abstract-derived rows). Cho 2011 directional rows (no OR) are kept as a
separate qualitative table.
"""
from __future__ import annotations
import math
from pathlib import Path
import pandas as pd
import numpy as np
from scipy import stats

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track10_korean_lit")
TBL = OUT / "tables"
WS  = OUT / "extraction_worksheets"

shin = pd.read_csv(WS / "shin2019_long.tsv", sep="\t")
park = pd.read_csv(WS / "park2005_long.tsv", sep="\t")
jang = pd.read_csv(WS / "jang2011_long.tsv", sep="\t")

# Canonical schema
COLS = [
    "paper_id", "pmid", "doi", "ancestry", "cohort", "typing", "comparison",
    "allele", "frequency_type",
    "n_case", "n_control", "n_case_carrier", "n_control_carrier",
    "n_case_chrom", "n_ctrl_chrom", "twoN_case", "twoN_ctrl",
    "OR", "OR_lo", "OR_hi", "se_log_or", "p", "Pc",
    "note",
]

# --- Shin 2019: keep all 4-digit rows for AITD, GD, HD comparisons ---
shin_out = shin.rename(columns={
    "OR": "OR", "OR_lo": "OR_lo", "OR_hi": "OR_hi", "p": "p", "Pc": "Pc",
}).copy()
# compute SE(log OR) from CI
def ci_to_se(lo, hi):
    if not (np.isfinite(lo) and np.isfinite(hi)) or lo <= 0 or hi <= 0:
        return np.nan
    return (math.log(hi) - math.log(lo)) / (2 * 1.959964)
shin_out["se_log_or"] = [ci_to_se(lo, hi) for lo, hi in zip(shin_out["OR_lo"], shin_out["OR_hi"])]
shin_out["n_case_chrom"] = np.nan
shin_out["n_ctrl_chrom"] = np.nan
shin_out["twoN_case"] = np.nan
shin_out["twoN_ctrl"] = np.nan
shin_out["note"] = "Shin 2019 S1 Table; carrier-frequency"

# --- Park 2005 ---
park_out = park.rename(columns={
    "OR_recomputed": "OR", "OR_lo_recomputed": "OR_lo", "OR_hi_recomputed": "OR_hi",
    "se_log_or_recomputed": "se_log_or", "p_published": "p", "Pc_published": "Pc",
}).copy()
park_out["n_case_carrier"] = np.nan
park_out["n_control_carrier"] = np.nan

# --- Jang 2011 ---
jang_out = jang.rename(columns={
    "OR_recomputed": "OR", "OR_lo_recomputed": "OR_lo", "OR_hi_recomputed": "OR_hi",
    "se_log_or_recomputed": "se_log_or",
}).copy()
jang_out["Pc"] = np.nan
jang_out["n_case_carrier"] = np.nan
jang_out["n_control_carrier"] = np.nan

master = pd.concat([shin_out, park_out, jang_out], ignore_index=True, sort=False)

# Reorder
for c in COLS:
    if c not in master.columns:
        master[c] = np.nan
master = master[COLS + [c for c in master.columns if c not in COLS]]
master.to_csv(TBL / "T03_korean_lit_extraction_master.tsv", sep="\t", index=False)
print(f"saved -> {TBL/'T03_korean_lit_extraction_master.tsv'}  ({len(master)} rows)")

# Print breakdown
print("\nBreakdown by paper × comparison:")
print(master.groupby(["paper_id", "comparison"]).size().reset_index(name="n_rows"))
