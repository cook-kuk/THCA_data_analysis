#!/usr/bin/env python3
"""
Track 10 — Step 3: Extract Shin 2019 (PMID 31091281) S1 Table.

S1 Table reports carrier-counts (n) and carrier-frequencies (%) for HLA-A/B/C/
DRB1/DQB1/DPB1 4-digit alleles in:
    Controls (n=142),  AITD (n=116),  GD (n=71),  HD (n=45)
plus 4 OR / 95% CI / p / Pc columns:
    Controls vs AITD,  Controls vs GD,  Controls vs HD,  GD vs HD.

This script flattens that into:
    - one row per (allele, comparison)
    - canonical columns matching Track-1 source schema for downstream meta-analysis.
"""
from __future__ import annotations
import re, zipfile, json
from pathlib import Path
import pandas as pd

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track10_korean_lit")
TBL = OUT / "tables"
WS  = OUT / "extraction_worksheets"
WS.mkdir(parents=True, exist_ok=True)

DOCX = OUT / "source_pdfs" / "Shin_2019_S1_table.docx"

with zipfile.ZipFile(DOCX) as z:
    xml = z.read("word/document.xml").decode("utf-8", errors="ignore")
text = re.sub(r"<[^>]+>", " ", xml)
text = re.sub(r"\s+", " ", text)

# Find the start of the table data (after the headers); skip everything before "Locus Alleles n"
m = re.search(r"Locus Alleles n = 142.*?Pc n = 71.*?Pc n = 45.*?Pc p value Pc", text)
if m is None:
    # Fallback: simpler anchor
    m = re.search(r"Locus Alleles n = 142", text)
    if m is None:
        raise SystemExit("S1 Table header not found")
start = m.end()
# The S1 Table block ends at the start of the next heading; safest cutoff: when "Pc" occurs followed by a non-numeric run with too few patterns
chunk = text[start:start + 200_000]

# Each allele row pattern (most rows). Loci collected: A, B, C, DRB1, DQB1, DPB1
# Row layout (when no Pc reported on a comparison):
#   <Locus> <4digit> <n_ctrl> (<pct>) <n_aitd> (<pct>) <OR_aitd> (<lo>-<hi>) <p_aitd> [<Pc_aitd>] <n_gd> (<pct>) <OR_gd> (<lo>-<hi>) <p_gd> [<Pc_gd>] <n_hd> (<pct>) <OR_hd> (<lo>-<hi>) <p_hd> [<Pc_hd>] <p_gd_vs_hd> [<Pc_gd_vs_hd>]
# Pc only when reported (some rows include Pc, others use "-" or omit). The table prints "-" for unestimable.
# We'll parse line-by-line by anchoring on a locus token.

LOCI = r"(?:DRB1|DQB1|DPB1|A|B|C)"
ALLELE_TOK = r"\d{2}:\d{2}(?:[NLQS])?"

# Numbers can be ints, decimals, or token "-"
NUM = r"(?:\d+\.\d+|\d+|-)"
PCT = r"\d+\.\d+|\d+\.0|\d+|0"
ORTOK = r"(?:\d+\.\d+|\d+|-)"
CITOK = r"(?:[\d\.\-]+|-)"

# Generic token for OR/CI sometimes "(0.05-3.07)"
CI = r"\([\d\.\-]+\s*-\s*[\d\.\-]+\)"

# A more lenient token-by-token approach: split chunk by the locus tokens
# But that breaks because "C" overlaps with "Cw" etc. Only A/B/C/DRB1/DQB1/DPB1 appear.
# Safer: walk allele patterns: locus space allele space rest until next locus
# We'll detect rows by regex anchored on "<LOCUS> <DIGITS>:<DIGITS>".

ROW_RE = re.compile(
    rf"({LOCI})\s+({ALLELE_TOK})\s+(\d+)\s*\(\s*([\d\.]+)\s*\)\s+(\d+)\s*\(\s*([\d\.]+)\s*\)\s+"
    # AITD comparison
    rf"({NUM})\s+\(([\d\.\-]+)\s*-\s*([\d\.\-]+)\)\s+([\d\.\-]+)(?:\s+([\d\.\-]+))?\s+"
    # GD arm: count (pct) OR (lo-hi) p [Pc]
    rf"(\d+)\s*\(\s*([\d\.]+)\s*\)\s+({NUM})\s+\(([\d\.\-]+)\s*-\s*([\d\.\-]+)\)\s+([\d\.\-]+)(?:\s+([\d\.\-]+))?\s+"
    # HD arm: count (pct) OR (lo-hi) p [Pc]
    rf"(\d+)\s*\(\s*([\d\.]+)\s*\)\s+({NUM})\s+\(([\d\.\-]+)\s*-\s*([\d\.\-]+)\)\s+([\d\.\-]+)(?:\s+([\d\.\-]+))?\s+"
    # GD vs HD: p [Pc]
    rf"([\d\.\-]+)(?:\s+([\d\.\-]+))?"
)

rows = []
for m in ROW_RE.finditer(chunk):
    g = m.groups()
    rec = dict(
        locus=g[0],
        allele=f"{g[0]}*{g[1]}",
        n_ctrl_carrier=int(g[2]), pct_ctrl=float(g[3]),
        n_aitd_carrier=int(g[4]), pct_aitd=float(g[5]),
        OR_aitd=g[6], OR_aitd_lo=g[7], OR_aitd_hi=g[8],
        p_aitd=g[9], Pc_aitd=g[10],
        n_gd_carrier=int(g[11]), pct_gd=float(g[12]),
        OR_gd=g[13], OR_gd_lo=g[14], OR_gd_hi=g[15],
        p_gd=g[16], Pc_gd=g[17],
        n_hd_carrier=int(g[18]), pct_hd=float(g[19]),
        OR_hd=g[20], OR_hd_lo=g[21], OR_hd_hi=g[22],
        p_hd=g[23], Pc_hd=g[24],
        p_gd_vs_hd=g[25], Pc_gd_vs_hd=g[26],
    )
    rows.append(rec)

print(f"parsed {len(rows)} S1 Table rows from Shin 2019")

if len(rows) < 30:
    # Fallback: print last unmatched portion
    print("\nWARNING: low row count. Tail of chunk:")
    print(chunk[-2000:])

df = pd.DataFrame(rows)
# Coerce numeric-as-string fields, treating "-" as NaN
for c in [k for k in df.columns if k.startswith(("OR_", "p_", "Pc_"))]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

df["n_ctrl"] = 142
df["n_aitd"] = 116
df["n_gd"]   = 71
df["n_hd"]   = 45
df["paper_id"] = "Shin_2019"
df["pmid"] = "31091281"
df["doi"]  = "10.1371/journal.pone.0216941"
df["typing_method"] = "Sanger SBT (HLA-A/B/C/DRB1/DQB1/DPB1 high-resolution)"
df["frequency_type"] = "carrier_frequency"  # carriers per subject (0/1/2 collapsed to ≥1)
df["table_source"] = "S1 Table (DOCX); columns Controls/AITD/GD/HD with 4 ORs"

df.to_csv(WS / "shin2019_S1_table_full.tsv", sep="\t", index=False)
print(f"saved -> {WS/'shin2019_S1_table_full.tsv'}  ({len(df)} rows)")

# Build a 'long' worksheet: one row per (allele, comparison)
def emit(allele, comp, n_case, n_ctrl, n_case_carrier, n_ctrl_carrier, OR, lo, hi, p, Pc):
    return dict(
        paper_id="Shin_2019", pmid="31091281", doi="10.1371/journal.pone.0216941",
        cohort="Korean pediatric (mean age ~12 y; SNUH/Catholic Univ)",
        typing="Sanger SBT high-resolution", ancestry="Korean",
        comparison=comp,
        allele=allele,
        n_case=n_case, n_control=n_ctrl,
        n_case_carrier=n_case_carrier, n_control_carrier=n_ctrl_carrier,
        OR=OR, OR_lo=lo, OR_hi=hi, p=p, Pc=Pc,
        frequency_type="carrier_frequency",
    )

long_rows = []
for _, r in df.iterrows():
    a = r["allele"]
    long_rows.append(emit(a, "Controls_vs_AITD", 116, 142, r["n_aitd_carrier"], r["n_ctrl_carrier"],
                          r["OR_aitd"], r["OR_aitd_lo"], r["OR_aitd_hi"], r["p_aitd"], r["Pc_aitd"]))
    long_rows.append(emit(a, "Controls_vs_GD",  71, 142, r["n_gd_carrier"], r["n_ctrl_carrier"],
                          r["OR_gd"], r["OR_gd_lo"], r["OR_gd_hi"], r["p_gd"], r["Pc_gd"]))
    long_rows.append(emit(a, "Controls_vs_HD",  45, 142, r["n_hd_carrier"], r["n_ctrl_carrier"],
                          r["OR_hd"], r["OR_hd_lo"], r["OR_hd_hi"], r["p_hd"], r["Pc_hd"]))
ldf = pd.DataFrame(long_rows)
ldf.to_csv(WS / "shin2019_long.tsv", sep="\t", index=False)
print(f"saved -> {WS/'shin2019_long.tsv'}  ({len(ldf)} rows)")

# Print headline numbers for the focus alleles
focus = ["DPB1*05:01", "B*46:01", "DRB1*08:02", "DRB1*15:01", "DRB1*16:02",
         "A*02:07", "C*03:02", "DQB1*03:02"]
print("\nFocus-allele carrier rows (Controls vs GD):")
for a in focus:
    sub = ldf[(ldf["allele"] == a) & (ldf["comparison"] == "Controls_vs_GD")]
    if len(sub):
        s = sub.iloc[0]
        print(f"  {a:14s}  GD {s['n_case_carrier']}/{s['n_case']}  CTRL {s['n_control_carrier']}/{s['n_control']}  OR={s['OR']}  ({s['OR_lo']}-{s['OR_hi']}) p={s['p']}  Pc={s['Pc']}")
    else:
        print(f"  {a:14s}  (not in S1 Table)")
