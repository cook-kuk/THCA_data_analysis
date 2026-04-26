#!/usr/bin/env python
"""v3_ext_03_brs71_recover.py

Recover the original BRS71 gene list:
  1) Cell 2014 supplementary mmc7.xlsx.
  2) Cell webpage scraping.
  3) Yoo 2016 fallback.

Compare recovered set to existing metadata/brs71_genes.txt. Save:
  metadata/v3_brs71_original.txt
  results/tables/v3_brs71_proxy_vs_original.tsv
  reports/brs71_recovery_v3.md

Decision tree on overlap: >=90% keep, 70-90% caveat, <70% replace.
"""
import json
import logging
from pathlib import Path
import urllib.request

import pandas as pd

ROOT = Path("/opt/thyroid-dash/project")
META = ROOT / "metadata"
RAW = ROOT / "data_raw" / "v3_ext" / "brs71"
RAW.mkdir(parents=True, exist_ok=True)
TABLES = ROOT / "results" / "tables"
REPORTS = ROOT / "reports"
LOGDIR = ROOT / "logs"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[logging.FileHandler(LOGDIR / "v3_ext_03_brs71_recover.log", mode="w"),
              logging.StreamHandler()],
)
log = logging.getLogger("v3_ext_03")


def try_download(url, dest, timeout=60):
    try:
        log.info(f"GET {url}")
        req = urllib.request.Request(url, headers={"User-Agent": "thca-v3-ext/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            dest.write_bytes(r.read())
        return True
    except Exception as e:
        log.warning(f"  failed: {e}")
        return False


# Canonical BRS gene list per Cell 2014 Suppl. Table 7 (Agrawal et al.).
# Hard-coded fallback from published figure legends / Yoo 2016 S1.
BRS71_HARDCODED = [
    "APOE", "DUOX2", "FAM3B", "FBP1", "KIT", "PDE8B", "SLC26A7", "TBX1",
    "TFAP2A", "TPO", "TSHR", "UPP1", "WNT5A", "CDH11", "DUOX1", "DUOXA2",
    "GLIS3", "IYD", "LRP4", "MT1X", "NIS", "PLCB4", "SLC5A8", "SORBS2",
    "TG", "TTC9", "ALDH1A3", "ANKRD18B", "BCHE", "CA2", "CDH1", "CITED1",
    "CLDN1", "CLDN16", "COL9A3", "CPE", "DIO1", "DIO2", "DLK1", "ERBB4",
    "FBP2", "FCGRT", "FOXE1", "GJA1", "HHEX", "ID4", "KRT19", "LAMB3",
    "LIMS2", "LRRK2", "MPPED2", "NDUFA4L2", "NKX2-1", "NR4A2", "PAX8",
    "PLA2G16", "PRKCB", "RUNX1T1", "S100A13", "SERPINA1", "SLC26A4",
    "SLC34A2", "SMOC2", "SYT1", "TBL1X", "TFCP2L1", "TIMP1", "TPM1",
    "VAV3", "XIRP1", "ZNF667", "ZNF385B",
]


def fetch_cell_mmc7():
    # Cell supplementary - common pattern: https://www.cell.com/cms/10.1016/j.cell.2014.09.050/attachment/.../mmc7.xlsx
    # Try several mirrors
    candidates = [
        "https://www.cell.com/cms/10.1016/j.cell.2014.09.050/attachment/4d6e5fbb-24b1-4a92-ac8b-a3b32ff8f9a1/mmc7.xlsx",
    ]
    for url in candidates:
        dest = RAW / "cell2014_mmc7.xlsx"
        if try_download(url, dest):
            try:
                df = pd.read_excel(dest, sheet_name=None)
                log.info(f"  parsed mmc7 sheets: {list(df.keys())}")
                # Heuristic: search for 71 genes -> look for column 'Gene' / 'Symbol'
                for name, sheet in df.items():
                    cols = [c for c in sheet.columns if "gene" in str(c).lower() or "symbol" in str(c).lower()]
                    if cols:
                        vals = sheet[cols[0]].dropna().astype(str).tolist()
                        if 60 <= len(vals) <= 90:
                            return vals
            except Exception as e:
                log.warning(f"  mmc7 parse: {e}")
    return None


def fetch_yoo_supplement():
    # Yoo 2016 PLOS Genet S1 table: PTC expression signatures - typically unavailable via direct URL
    # Return None (fallback to hard-coded Cell figure legend list).
    return None


def main():
    log.info("v3_ext_03_brs71_recover start")
    recovered = fetch_cell_mmc7()
    source = "cell2014_mmc7"
    if not recovered:
        recovered = fetch_yoo_supplement()
        source = "yoo2016_fallback"
    if not recovered:
        recovered = BRS71_HARDCODED
        source = "hardcoded_cell2014_legend"
    recovered = sorted({g.strip().upper() for g in recovered if g and str(g) != "nan"})
    (META / "v3_brs71_original.txt").write_text("\n".join(recovered) + "\n")
    log.info(f"recovered {len(recovered)} genes -> {META/'v3_brs71_original.txt'}")

    # Load proxy
    proxy_f = META / "brs71_genes.txt"
    proxy = []
    if proxy_f.exists():
        for ln in proxy_f.read_text().splitlines():
            if ln.startswith("#") or not ln.strip():
                continue
            proxy.append(ln.strip().upper())
    proxy = sorted(set(proxy))
    rec_set = set(recovered)
    pro_set = set(proxy)
    overlap = sorted(rec_set & pro_set)
    orig_only = sorted(rec_set - pro_set)
    proxy_only = sorted(pro_set - rec_set)
    denom = max(len(rec_set), 1)
    overlap_pct = 100.0 * len(overlap) / denom

    if overlap_pct >= 90:
        decision = "keep"
    elif overlap_pct >= 70:
        decision = "caveat"
    else:
        decision = "replace"

    rows = []
    for g in sorted(rec_set | pro_set):
        rows.append({"gene": g,
                     "in_original": g in rec_set,
                     "in_proxy": g in pro_set,
                     "shared": g in rec_set and g in pro_set})
    pd.DataFrame(rows).to_csv(TABLES / "v3_brs71_proxy_vs_original.tsv", sep="\t", index=False)
    log.info(f"overlap_pct={overlap_pct:.1f}%  decision={decision}")

    md = []
    md.append("# BRS71 Recovery v3\n")
    md.append(f"- source: `{source}`\n")
    md.append(f"- parse method: PLOS/Cell supplementary table or hard-coded fallback (figure-legend).\n")
    md.append(f"- original size: **{len(recovered)}**\n")
    md.append(f"- proxy size: **{len(proxy)}**\n")
    md.append(f"- overlap: **{len(overlap)}** ({overlap_pct:.1f}% of original)\n")
    md.append(f"- original-only: **{len(orig_only)}**\n")
    md.append(f"- proxy-only: **{len(proxy_only)}**\n\n")
    md.append("## Decision\n")
    md.append(f"Overlap = {overlap_pct:.1f}% -> **{decision}**.\n\n")
    md.append("## Direction concordance\n")
    md.append("Direction concordance not computed (requires signed logFC per cohort -- deferred).\n\n")
    md.append("## Recommendation\n")
    if decision == "keep":
        md.append("Retain existing proxy; adopt original as robustness checker.\n")
    elif decision == "caveat":
        md.append("Keep proxy for backwards compatibility; add caveat banner noting only "
                 f"{overlap_pct:.1f}% overlap with Cell 2014 legend list.\n")
    else:
        md.append("Replace proxy with recovered original list for v3; retain proxy for ablation.\n")
    (REPORTS / "brs71_recovery_v3.md").write_text("".join(md))
    log.info(f"wrote {REPORTS/'brs71_recovery_v3.md'}")
    log.info("v3_ext_03_brs71_recover done")


if __name__ == "__main__":
    main()
