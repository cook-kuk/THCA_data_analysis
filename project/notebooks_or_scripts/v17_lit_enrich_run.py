#!/usr/bin/env python3
"""
v17 lit_enrich runner — multi-source literature enrichment for the THCA papers.

Hits 12 free academic-data APIs (OpenAlex, CrossRef, PubMed, Europe PMC,
Semantic Scholar, arXiv, bioRxiv, ICite, AFND, ClinicalTrials.gov, Unpaywall,
CORE) and writes 7 reports + JSON data into manuscript_p2_brief/lit_enrich_<date>/.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "project" / "notebooks_or_scripts"))

from helpers.lit_enrich import tasks


def main() -> int:
    bib = ROOT / "project" / "manuscript_v8" / "03_intro_references.bib"
    out_dir = ROOT / "project" / "manuscript_p2_brief" / f"lit_enrich_{time.strftime('%Y_%m_%d')}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[lit_enrich] bib       = {bib}")
    print(f"[lit_enrich] out_dir   = {out_dir}\n")

    results: dict = {}

    # ── Task 1: verify_references ────────────────────────────────────────────
    print("[1/6] verify_references …")
    t0 = time.time()
    r1 = tasks.verify_references(bib, out_dir)
    print(f"      → {r1['incomplete']} incomplete / {r1['entries']} entries "
          f"({time.time()-t0:.1f}s)\n")
    results["verify_references"] = {k: v for k, v in r1.items() if k != "results"}

    # ── Task 2: competitive_landscape ────────────────────────────────────────
    claims = [
        {"topic": "HLA-II PTC autoimmunity (DPB1*05:01 etc.)",
         "query": "HLA DPB1 papillary thyroid Hashimoto autoimmune", "year_from": 2018},
        {"topic": "Hashimoto-like signature & PTC outcomes",
         "query": "Hashimoto thyroiditis papillary thyroid cancer transcriptome", "year_from": 2020},
        {"topic": "BCR repertoire / TLS in thyroid cancer",
         "query": "tertiary lymphoid structures B cell repertoire thyroid cancer", "year_from": 2020},
        {"topic": "BRAF/RAS-negative ('dark matter') PTC subtypes",
         "query": "BRAF RAS negative papillary thyroid driver-negative", "year_from": 2018},
        {"topic": "Single-cell PTC progression",
         "query": "single-cell papillary thyroid carcinoma progression", "year_from": 2021},
        {"topic": "Pan-Asian HLA fine-mapping (Graves' / autoimmune thyroid)",
         "query": "HLA fine-mapping Graves autoimmune thyroid Asian", "year_from": 2016},
        {"topic": "TIERA-like Korean PTC molecular cohort",
         "query": "Korean papillary thyroid cancer transcriptome molecular subtype", "year_from": 2018},
        {"topic": "8-gene / driver-excluded PTC stratification",
         "query": "papillary thyroid prognostic gene signature stratification", "year_from": 2018},
    ]
    print("[2/6] competitive_landscape …")
    t0 = time.time()
    r2 = tasks.competitive_landscape(claims, out_dir)
    print(f"      → {r2['claims']} claims mined ({time.time()-t0:.1f}s)\n")
    results["competitive_landscape"] = {"claims": r2["claims"]}

    # ── Task 3: hla_frequencies ──────────────────────────────────────────────
    alleles = ["DPB1*05:01", "DRB1*04:05", "DRB1*15:01", "DQB1*06:02", "B*46:01"]
    populations = ["Korea", "China", "Japan", "Taiwan"]
    print("[3/6] hla_frequencies …")
    t0 = time.time()
    r3 = tasks.hla_frequencies(alleles, populations, out_dir)
    print(f"      → {r3['alleles']} alleles × {r3['populations']} pops "
          f"({time.time()-t0:.1f}s)\n")
    results["hla_frequencies"] = {"alleles": r3["alleles"], "populations": r3["populations"]}

    # ── Task 4: misattribution_check ─────────────────────────────────────────
    print("[4/6] misattribution_check …")
    t0 = time.time()
    r4 = tasks.misattribution_check(bib, out_dir)
    print(f"      → {r4['issues']} mismatches in {r4['checked']} DOI'd entries "
          f"({time.time()-t0:.1f}s)\n")
    results["misattribution_check"] = r4

    # ── Task 5: clinical_landscape ───────────────────────────────────────────
    trial_queries = [
        "papillary thyroid cancer immunotherapy",
        "RET fusion thyroid",
        "BRAF V600E thyroid",
        "anaplastic thyroid carcinoma",
        "thyroid cancer Hashimoto",
    ]
    print("[5/6] clinical_landscape …")
    t0 = time.time()
    r5 = tasks.clinical_landscape(trial_queries, out_dir)
    print(f"      → {r5['queries']} trial queries ({time.time()-t0:.1f}s)\n")
    results["clinical_landscape"] = r5

    # ── Task 6: full_text_discovery ──────────────────────────────────────────
    print("[6/6] full_text_discovery …")
    t0 = time.time()
    r6 = tasks.full_text_discovery(bib, out_dir)
    print(f"      → {r6['oa_count']}/{r6['entries']} open-access "
          f"({time.time()-t0:.1f}s)\n")
    results["full_text_discovery"] = r6

    # ── Summary ──────────────────────────────────────────────────────────────
    tasks.write_summary(out_dir, results)
    print(f"\n[lit_enrich] DONE → {out_dir}")
    print("              Reports: 01_…06_…07_summary.md")
    print("              Data:    data/*.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
