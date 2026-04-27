#!/usr/bin/env python
"""v17 ULTIMATE U3A — Mu 2024 + RAI ground-truth landscape search.

Light wrapper over WebFetch/WebSearch results that were collected in the
parent agent turn (online lookups already executed). This script bakes
the harvested findings into three on-disk artifacts:

    project/results/v17_ultimate/U3A_mu2024_data_status.md
    project/results/v17_ultimate/U3A_search_log.md
    project/results/v17_ultimate/U3A_rai_landscape.tsv

NO paywall bypass. NO Sci-Hub. Honest accounting only.

Author: Seungho Cook
Date:   2026-04-27
"""
from __future__ import annotations

import csv
import datetime as _dt
from pathlib import Path

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project" / "results" / "v17_ultimate"
OUT.mkdir(parents=True, exist_ok=True)

TODAY = _dt.date.today().isoformat()

# ---------------------------------------------------------------------------
# 1. Findings harvested via WebFetch / WebSearch (already executed)
# ---------------------------------------------------------------------------
MU2024 = {
    "title": (
        "Characterizing Genetic Alterations Related to Radioiodine "
        "Avidity in Metastatic Thyroid Cancer"
    ),
    "authors": "Mu Z, Zhang X, Sun D, Sun Y, Shi C, et al.",
    "journal": "J Clin Endocrinol Metab",
    "year": 2024,
    "volume_pages": "109(5):1231-1240",
    "pmid": "38060243",
    "pmcid": "PMC11031230",
    "doi": "10.1210/clinem/dgad697",
    "n_patients": 214,
    "n_samples": 281,
    "rai_labels": "I-RAIA (n=134) vs I-RAIR (n=80); subgroups C-RAIA / P-RAIR / G-RAIR",
    "data_accession": "HRA004166 (NGDC GSA-Human)",
    "access": "CONTROLLED — DAC application required",
    "dac_contact": "Lin Yansong, linys@pumch.cn (Peking Union Medical College Hospital)",
    "supp_accessible": False,
    "reason": (
        "Data Availability: 'Restrictions apply ... to preserve patient "
        "confidentiality'. NGS deposited under HRA004166 at NGDC requires "
        "Data Access Committee approval; no openly downloadable variant table "
        "or supplementary mutation matrix on Oxford Academic / PMC."
    ),
}

# Cohort, n, RAI_label_type, access, score (1-5; 5=best for our pipeline)
RAI_COHORTS = [
    {
        "cohort": "GSE151181 (SuperSeries: GSE151179 mRNA + GSE151180 miRNA)",
        "n": 99,
        "rai_label_type": "RAI-avid vs RAI-refractory PTC (paired primary + LNM)",
        "platform": "Affymetrix Clariom S + Agilent miRNA V21",
        "pmid": "33198784",
        "access": "OPEN (Public; CEL + series matrix downloadable)",
        "score": 5,
        "notes": "Already targeted by U1B. Best public RAI ground-truth.",
    },
    {
        "cohort": "GSE151179 (mRNA subseries)",
        "n": 52,
        "rai_label_type": "17 primary PTC + 5 synchronous LN + 13 normal + 17 post-RAI refractory LNM",
        "platform": "Affymetrix Clariom S (GPL23159)",
        "pmid": "33198784",
        "access": "OPEN",
        "score": 5,
        "notes": "Direct mRNA arm — pairs with our v17 expression workflow.",
    },
    {
        "cohort": "GSE151180 (miRNA subseries)",
        "n": 47,
        "rai_label_type": "RAI-avid vs RAI-refractory (miRNA companion)",
        "platform": "Agilent-070156 miRNA V21 (GPL21575)",
        "pmid": "33198784",
        "access": "OPEN",
        "score": 4,
        "notes": "Useful for miRNA cross-modality validation.",
    },
    {
        "cohort": "GSE299988",
        "n": 14,
        "rai_label_type": "PTC w/ LN metastasis, RAI-refractory",
        "platform": "RNA-seq",
        "pmid": "n/a",
        "access": "OPEN",
        "score": 2,
        "notes": "Small n; refractory-only (no avid comparator).",
    },
    {
        "cohort": "GSE171473",
        "n": 9,
        "rai_label_type": "RAI-refractory; YK-4-279 / TERT context",
        "platform": "RNA-seq",
        "pmid": "n/a",
        "access": "OPEN",
        "score": 2,
        "notes": "TERT mechanistic; too small as standalone test set.",
    },
    {
        "cohort": "GSE162525 / 162523 / 162522 / 147475-78 (SWI/SNF series)",
        "n": 156,
        "rai_label_type": "Redifferentiation-resistance (RAI-related but not direct RAIA/RAIR)",
        "platform": "RNA-seq + ChIP-seq + ATAC-seq",
        "pmid": "n/a",
        "access": "OPEN",
        "score": 3,
        "notes": "Mechanistic multi-omic; useful as orthogonal validation.",
    },
    {
        "cohort": "Mu 2024 / HRA004166",
        "n": 214,
        "rai_label_type": "I-RAIA vs I-RAIR (gold-standard clinical avidity)",
        "platform": "Targeted NGS panel (BRAF/RAS/TERT/TP53 + fusions)",
        "pmid": "38060243",
        "access": "CONTROLLED (NGDC DAC; not paywalled — credentialed)",
        "score": 5,
        "notes": "Best clinical labels but requires DAC approval; consider formal request.",
    },
    {
        "cohort": "GSE190966",
        "n": 12,
        "rai_label_type": "NOT THYROID — Bos indicus muscle methylation",
        "platform": "RRBS",
        "pmid": "n/a",
        "access": "OPEN but irrelevant",
        "score": 0,
        "notes": "Accession exists but is bovine; exclude.",
    },
]

URLS_TRIED = [
    ("PubMed Mu 2024", "WebSearch query 'Mu Z 2024 thyroid radioiodine avidity JCEM 109 1231'", "FOUND PMID 38060243"),
    ("Oxford Academic JCEM", "https://academic.oup.com/jcem/article/109/5/1231/7460630", "Found article landing page (paywalled body)"),
    ("PMC full text", "https://pmc.ncbi.nlm.nih.gov/articles/PMC11031230/", "Open-access PMC version available; Data Availability extracted"),
    ("PubMed", "https://pubmed.ncbi.nlm.nih.gov/38060243/", "Confirmed PMID + abstract"),
    ("Ovid mirror", "https://www.ovid.com/journals/jceme/fulltext/10.1210/clinem/dgad697...", "Paywalled"),
    ("NGDC GSA-Human HRA004166", "https://ngdc.cncb.ac.cn/gsa-human/browse/HRA004166", "CONTROLLED ACCESS; DAC contact: linys@pumch.cn"),
    ("GEO GSE151181", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151181", "OPEN; n=99 SuperSeries"),
    ("GEO GSE151179", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151179", "OPEN; n=52 mRNA"),
    ("GEO GSE151180", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151180", "OPEN; n=47 miRNA"),
    ("GEO GSE190966", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190966", "OPEN but bovine — exclude"),
    ("GEO GDS landscape", "https://www.ncbi.nlm.nih.gov/gds/?term=thyroid+AND+radioiodine+AND+expression", "15 thyroid+RAI series enumerated"),
    ("Sci-Hub / preprint", "(skipped per ethics)", "NOT ATTEMPTED"),
]

# ---------------------------------------------------------------------------
# 2. Emit U3A_mu2024_data_status.md
# ---------------------------------------------------------------------------
status_md = OUT / "U3A_mu2024_data_status.md"
status_md.write_text(
    f"""# U3A — Mu Z et al. JCEM 2024 data accessibility status

_Date: {TODAY}_  •  _Author: Seungho Cook_  •  _v17 ULTIMATE sprint_

## TL;DR

**Mu 2024 supplementary / raw data accessible? NO (controlled).**

The Oxford Academic / PMC version of the paper is open-access for *reading*,
but the underlying NGS data is deposited under **NGDC accession
`HRA004166`** with **controlled access** — a Data Access Committee (DAC)
application is required. The published supplement does not contain a
patient-level mutation × RAI-avidity table that would allow us to bypass
the DAC.

## Citation

- {MU2024['authors']}. *{MU2024['title']}*.
  **{MU2024['journal']}** {MU2024['year']};{MU2024['volume_pages']}.
- PMID: {MU2024['pmid']}  •  PMCID: {MU2024['pmcid']}  •  DOI: {MU2024['doi']}

## Cohort

- N patients: **{MU2024['n_patients']}** (DM-DTC)
- N samples : **{MU2024['n_samples']}**
- RAI labels: {MU2024['rai_labels']}

## Access verdict

- Data accession: **{MU2024['data_accession']}**
- Access tier  : **{MU2024['access']}**
- DAC contact  : {MU2024['dac_contact']}
- Reason       : {MU2024['reason']}

## What we did NOT do (and why)

- We did **not** attempt Sci-Hub or other paywall-bypass mirrors.
  The paper itself is already open-access on PMC; the controlled item is
  patient-level NGS, where bypass would be both unethical and illegal.
- We did **not** scrape NGDC. A DAC request is the correct path.

## Recommended next step

1. Lean on the open RAI-avid/refractory cohorts in `U3A_rai_landscape.tsv`
   (GSE151181 SuperSeries, n=99) as the primary external validation set —
   already in flight in U1B.
2. **Optionally** submit a DAC request to Lin Yansong (linys@pumch.cn) for
   HRA004166 if reviewers ask for clinical-grade RAI labels beyond GEO.
3. Treat Mu 2024 as **citation + concordance reference**, not as a
   downloadable test set, in the manuscript.
""",
    encoding="utf-8",
)

# ---------------------------------------------------------------------------
# 3. Emit U3A_search_log.md
# ---------------------------------------------------------------------------
log_md = OUT / "U3A_search_log.md"
lines = [
    "# U3A — Search log",
    "",
    f"_Date: {TODAY}_",
    "",
    "| # | Source | URL / query | Outcome |",
    "|---|--------|-------------|---------|",
]
for i, (src, url, outcome) in enumerate(URLS_TRIED, 1):
    lines.append(f"| {i} | {src} | `{url}` | {outcome} |")
lines += [
    "",
    "## Notes",
    "",
    "- All lookups executed via the agent's WebFetch / WebSearch tools on "
    f"{TODAY}.",
    "- Sci-Hub / preprint paywall-bypass routes were intentionally skipped.",
    "- HRA004166 (Mu 2024 NGS) is **controlled access** — DAC application "
    "required; no scraping.",
    "- GEO landscape query returned 15 thyroid+RAI series; the GSE151179/180/181 "
    "trio is the strongest open RAI-avid vs RAI-refractory ground truth "
    "(n=99, PMID 33198784).",
]
log_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

# ---------------------------------------------------------------------------
# 4. Emit U3A_rai_landscape.tsv
# ---------------------------------------------------------------------------
tsv_path = OUT / "U3A_rai_landscape.tsv"
fields = ["cohort", "n", "rai_label_type", "platform", "pmid", "access", "score", "notes"]
with tsv_path.open("w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for row in RAI_COHORTS:
        w.writerow({k: row.get(k, "") for k in fields})

# ---------------------------------------------------------------------------
# 5. Console summary
# ---------------------------------------------------------------------------
print("=" * 72)
print("U3A — Mu 2024 + RAI ground-truth landscape")
print("=" * 72)
print(f"Output dir : {OUT}")
print(f"Status MD  : {status_md.name}")
print(f"Search log : {log_md.name}")
print(f"Landscape  : {tsv_path.name}  ({len(RAI_COHORTS)} cohorts)")
print()
print(f"Mu 2024 supp accessible? NO  (controlled via NGDC HRA004166)")
print(f"Best OPEN alternative   : GSE151181 SuperSeries (n=99, PMID 33198784)")
print(f"Total open thyroid+RAI series catalogued: {sum(1 for r in RAI_COHORTS if 'OPEN' in r['access'])}")
print("Done.")
