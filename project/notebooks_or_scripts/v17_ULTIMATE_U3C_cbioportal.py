#!/usr/bin/env python
"""
v17 ULTIMATE U3C: Cross-cohort thyroid mutation prevalence from cBioPortal API.

Pulls TERT promoter (C228T/C250T), BRAF V600E, RAS hotspots (HRAS/KRAS/NRAS
codon 12/13/61), and TP53 mutations across all available thyroid studies on
cBioPortal, and computes prevalence per histology subtype.

Outputs:
  - U3C_cbioportal_studies_list.tsv
  - U3C_per_study_mutation_counts.tsv
  - U3C_supplementary_table_S20.tsv
  - U3C_summary.json

Author: Seungho Cook
Date: 2026-04-27
"""

from __future__ import annotations

import json
import os
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import requests

API_ROOT = "https://www.cbioportal.org/api"
OUT_DIR = Path("/home/seungho/personal/THCA_data_analysis/project/results/v17_ultimate")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Genes of interest. cBioPortal uses Entrez IDs for fetch endpoints.
# TERT=7015, BRAF=673, HRAS=3265, KRAS=3845, NRAS=4893, TP53=7157
GENES: Dict[str, int] = {
    "TERT": 7015,
    "BRAF": 673,
    "HRAS": 3265,
    "KRAS": 3845,
    "NRAS": 4893,
    "TP53": 7157,
}

# RAS hotspot codons
RAS_HOTSPOT_CODONS = {12, 13, 61}

# TERT promoter genomic coordinates (GRCh37/hg19) -- chr5
# C228T = chr5:1295228 G>A   (-124 from ATG)
# C250T = chr5:1295250 G>A   (-146 from ATG)
# (Reference convention varies; we use position +/- 5 bp window for tolerance.)
TERT_C228_POS = 1295228
TERT_C250_POS = 1295250
TERT_PROMOTER_WINDOW = (1295113, 1295300)

REQ_SLEEP = 0.5
MAX_RETRIES = 5


def request_with_retry(method: str, url: str, **kwargs) -> Optional[requests.Response]:
    """HTTP call with backoff on 429/5xx; returns Response or None on terminal failure."""
    backoff = 1.0
    for attempt in range(MAX_RETRIES):
        try:
            resp = requests.request(method, url, timeout=60, **kwargs)
        except requests.RequestException as exc:
            print(f"  ! request error ({exc}); retry {attempt+1}/{MAX_RETRIES}")
            time.sleep(backoff)
            backoff *= 2
            continue

        if resp.status_code == 200:
            time.sleep(REQ_SLEEP)
            return resp
        if resp.status_code in (429, 500, 502, 503, 504):
            ra = resp.headers.get("Retry-After")
            wait = float(ra) if ra and ra.isdigit() else backoff
            print(f"  ! {resp.status_code} on {url}; waiting {wait:.1f}s "
                  f"(attempt {attempt+1}/{MAX_RETRIES})")
            time.sleep(wait)
            backoff *= 2
            continue
        # Non-retryable
        print(f"  ! {resp.status_code} on {url}: {resp.text[:200]}")
        return resp
    return None


def list_thyroid_studies() -> List[Dict[str, Any]]:
    """Find all studies on cBioPortal whose name/id mentions thyroid/THCA."""
    print("[1/5] Fetching study list...")
    resp = request_with_retry("GET", f"{API_ROOT}/studies",
                              params={"projection": "SUMMARY", "pageSize": 10000})
    if resp is None or resp.status_code != 200:
        return []
    studies = resp.json()
    keywords = ("thyroid", "thca", "thyr")
    out = []
    for st in studies:
        blob = " ".join(str(st.get(k, "")) for k in ("studyId", "name",
                                                     "shortName", "description",
                                                     "cancerTypeId")).lower()
        if any(k in blob for k in keywords):
            out.append(st)
    # Always include curated thyroid studies even if name filter missed them.
    curated_ids = {"thca_tcga", "thca_tcga_pub", "thca_tcga_pan_can_atlas_2018",
                   "thca_msk_2016", "thca_msk_chord_2024"}
    have = {s["studyId"] for s in out}
    for sid in curated_ids - have:
        # Verify it exists.
        r = request_with_retry("GET", f"{API_ROOT}/studies/{sid}")
        if r is not None and r.status_code == 200:
            out.append(r.json())
    print(f"   -> {len(out)} thyroid-related studies found")
    return out


def get_mutation_profile_id(study_id: str) -> Optional[str]:
    """Find the mutation molecular profile for a study.

    The /molecular-profiles list endpoint does not honour the studyId param on
    every deployment, so we hit the per-study endpoint and filter ourselves.
    """
    resp = request_with_retry(
        "GET", f"{API_ROOT}/studies/{study_id}/molecular-profiles",
        params={"projection": "SUMMARY"})
    if resp is None or resp.status_code != 200:
        # Fall back to global list with explicit study filter
        resp = request_with_retry("GET", f"{API_ROOT}/molecular-profiles",
                                   params={"projection": "SUMMARY",
                                           "pageSize": 100000})
        if resp is None or resp.status_code != 200:
            return None
        candidates = [p for p in resp.json() if p.get("studyId") == study_id]
    else:
        candidates = resp.json()

    # Prefer the canonical "<study>_mutations" profile when present.
    canonical = f"{study_id}_mutations"
    for prof in candidates:
        if (prof.get("molecularAlterationType") == "MUTATION_EXTENDED"
                and prof.get("molecularProfileId") == canonical
                and prof.get("studyId") == study_id):
            return canonical
    for prof in candidates:
        if (prof.get("molecularAlterationType") == "MUTATION_EXTENDED"
                and prof.get("studyId") == study_id):
            return prof.get("molecularProfileId")
    return None


def get_sample_list_id(study_id: str, profile_id: str) -> Optional[str]:
    """Get the sequenced-samples sample-list id for the given mutation profile."""
    # Standard convention: <study>_sequenced
    cand = f"{study_id}_sequenced"
    r = request_with_retry("GET", f"{API_ROOT}/sample-lists/{cand}")
    if r is not None and r.status_code == 200:
        return cand
    # Otherwise list and look for "sequenced"
    r = request_with_retry("GET", f"{API_ROOT}/studies/{study_id}/sample-lists")
    if r is None or r.status_code != 200:
        return None
    for sl in r.json():
        cat = (sl.get("category") or "").lower()
        if "sequenc" in cat or "all_cases_with_mutation" in (sl.get("sampleListId") or ""):
            return sl["sampleListId"]
    # Fallback: all samples
    return f"{study_id}_all"


def fetch_clinical(study_id: str) -> Dict[str, Dict[str, str]]:
    """Return {sampleId: {attr: value}} for SAMPLE-level clinical data."""
    resp = request_with_retry(
        "GET",
        f"{API_ROOT}/studies/{study_id}/clinical-data",
        params={"clinicalDataType": "SAMPLE", "projection": "SUMMARY",
                "pageSize": 100000},
    )
    if resp is None or resp.status_code != 200:
        return {}
    out: Dict[str, Dict[str, str]] = defaultdict(dict)
    for row in resp.json():
        sid = row.get("sampleId")
        if not sid:
            continue
        out[sid][row["clinicalAttributeId"]] = row.get("value", "")
    return out


def histology_for(sample: str, clin: Dict[str, Dict[str, str]]) -> str:
    """Best-effort histology / cancer-type-detailed lookup."""
    c = clin.get(sample, {})
    for key in ("CANCER_TYPE_DETAILED", "HISTOLOGICAL_DIAGNOSIS",
                "TUMOR_TYPE", "CANCER_TYPE", "SUBTYPE",
                "ONCOTREE_CODE", "TUMOR_TISSUE_SITE"):
        v = c.get(key)
        if v:
            return v
    return "NA"


def is_thyroid_sample(sample: str, clin: Dict[str, Dict[str, str]]) -> bool:
    """For pan-cancer studies, check if a sample is thyroid origin."""
    c = clin.get(sample, {})
    blob = " ".join(str(v) for v in c.values()).lower()
    return ("thyroid" in blob or "thca" in blob
            or c.get("ONCOTREE_CODE", "").upper() in {"THPA", "THFO", "THAP",
                                                       "THHC", "THME", "THYC",
                                                       "THCA"})


def fetch_mutations(profile_id: str, sample_list_id: str,
                    entrez_ids: Iterable[int]) -> List[Dict[str, Any]]:
    """POST mutations/fetch for the given gene set, scoped to sample list."""
    body = {
        "entrezGeneIds": list(entrez_ids),
        "sampleListId": sample_list_id,
    }
    url = (f"{API_ROOT}/molecular-profiles/{profile_id}/mutations/fetch"
           f"?projection=DETAILED")
    resp = request_with_retry("POST", url, json=body,
                              headers={"Content-Type": "application/json"})
    if resp is None or resp.status_code != 200:
        return []
    return resp.json()


# ----- mutation classification helpers -----

def classify_tert(mut: Dict[str, Any]) -> Optional[str]:
    """Return 'C228T', 'C250T', 'TERT_promoter_other', or None for non-promoter."""
    pos = mut.get("startPosition") or mut.get("startPos")
    chrom = str(mut.get("chr") or "").lstrip("chr")
    pchg = (mut.get("proteinChange") or "").strip()

    # Promoter mutations are intronic/UTR-ish; protein change often empty.
    if pos and chrom in {"5", "5.0"}:
        try:
            ipos = int(pos)
        except (TypeError, ValueError):
            ipos = -1
        if abs(ipos - TERT_C228_POS) <= 2:
            return "C228T"
        if abs(ipos - TERT_C250_POS) <= 2:
            return "C250T"
        if TERT_PROMOTER_WINDOW[0] <= ipos <= TERT_PROMOTER_WINDOW[1]:
            return "TERT_promoter_other"

    # Some studies report HGVS like "C228T" / "C250T" in proteinChange / mutationType
    aa = pchg.upper()
    if "C228T" in aa:
        return "C228T"
    if "C250T" in aa:
        return "C250T"

    mtype = (mut.get("mutationType") or "").lower()
    if "5'flank" in mtype or "promoter" in mtype or "5'utr" in mtype:
        return "TERT_promoter_other"
    return None


def classify_braf(mut: Dict[str, Any]) -> Optional[str]:
    pchg = (mut.get("proteinChange") or "").upper()
    if "V600E" in pchg:
        return "V600E"
    if pchg.startswith("V600"):
        return f"V600_{pchg[4:] or 'other'}"
    return "BRAF_other"


def classify_ras(mut: Dict[str, Any], gene: str) -> Optional[str]:
    pchg = (mut.get("proteinChange") or "").upper().lstrip("P.")
    if not pchg:
        return None
    # Parse codon: e.g. Q61R -> 61
    digits = ""
    for ch in pchg[1:]:
        if ch.isdigit():
            digits += ch
        else:
            break
    try:
        codon = int(digits) if digits else -1
    except ValueError:
        codon = -1
    if codon in RAS_HOTSPOT_CODONS:
        return f"{gene}_codon{codon}"
    return f"{gene}_other"


def classify_tp53(mut: Dict[str, Any]) -> str:
    return "TP53_any"


# -------------- main -----------------

def process_study(study: Dict[str, Any]) -> Dict[str, Any]:
    """Return per-study mutation summary dict."""
    sid = study["studyId"]
    print(f"\n[Study] {sid} -- {study.get('name','')}")

    profile_id = get_mutation_profile_id(sid)
    if profile_id is None:
        print(f"   ! no MUTATION_EXTENDED profile for {sid}; skipping mutations")
        return {"studyId": sid, "name": study.get("name", ""),
                "n_samples": study.get("allSampleCount", 0),
                "has_mutations": False, "mutations": [], "samples": [],
                "clinical": {}}

    sample_list_id = get_sample_list_id(sid, profile_id) or f"{sid}_all"
    print(f"   profile={profile_id}  sample_list={sample_list_id}")

    clin = fetch_clinical(sid)
    print(f"   clinical rows for {len(clin)} samples")

    muts = fetch_mutations(profile_id, sample_list_id, GENES.values())
    print(f"   pulled {len(muts)} mutation records across {len(GENES)} genes")

    # If pan-cancer-style (msk_impact_2017), restrict to thyroid samples
    is_pan = "msk_impact_2017" in sid or "msk_chord" in sid
    if is_pan:
        thy = {s for s in clin if is_thyroid_sample(s, clin)}
        muts = [m for m in muts if m.get("sampleId") in thy]
        # samples tested = thyroid samples in sequenced list (use clin keys ∩ thy)
        sample_universe = thy
        print(f"   pan-cancer study; restricting to {len(thy)} thyroid samples")
    else:
        # Need the list of sequenced samples
        r = request_with_retry("GET",
                               f"{API_ROOT}/sample-lists/{sample_list_id}/sample-ids")
        if r is not None and r.status_code == 200:
            sample_universe = set(r.json())
        else:
            sample_universe = set(clin.keys())

    return {
        "studyId": sid,
        "name": study.get("name", ""),
        "n_samples": len(sample_universe),
        "has_mutations": True,
        "mutations": muts,
        "samples": sorted(sample_universe),
        "clinical": clin,
    }


def aggregate(study_results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Compute prevalence tables and summary."""
    # rows: dict keyed by (study, gene, hotspot, histology) -> {n_mut, n_test}
    rows: Dict[tuple, Dict[str, Any]] = defaultdict(lambda: {"mutated": set(),
                                                              "tested": set()})

    # Build: for every (study, histology), all tested samples
    for sres in study_results:
        if not sres.get("has_mutations"):
            continue
        sid = sres["studyId"]
        clin = sres["clinical"]
        for samp in sres["samples"]:
            hist = histology_for(samp, clin)
            # Initialise tested set for every gene/hotspot we report
            for gene in GENES:
                if gene == "TERT":
                    hotspots = ["C228T", "C250T", "TERT_promoter_other", "TERT_any"]
                elif gene == "BRAF":
                    hotspots = ["V600E", "BRAF_other", "BRAF_any"]
                elif gene in {"HRAS", "KRAS", "NRAS"}:
                    hotspots = [f"{gene}_codon12", f"{gene}_codon13",
                                f"{gene}_codon61", f"{gene}_other", f"{gene}_any"]
                elif gene == "TP53":
                    hotspots = ["TP53_any"]
                else:
                    hotspots = [f"{gene}_any"]
                for hs in hotspots:
                    rows[(sid, gene, hs, hist)]["tested"].add(samp)

        # Mark mutated samples
        for mut in sres["mutations"]:
            samp = mut.get("sampleId")
            entrez = mut.get("entrezGeneId")
            gene = next((g for g, e in GENES.items() if e == entrez), None)
            if gene is None:
                continue
            hist = histology_for(samp, clin)

            if gene == "TERT":
                cl = classify_tert(mut)
                if cl is None:
                    continue
                rows[(sid, "TERT", cl, hist)]["mutated"].add(samp)
                rows[(sid, "TERT", "TERT_any", hist)]["mutated"].add(samp)
            elif gene == "BRAF":
                cl = classify_braf(mut) or "BRAF_other"
                rows[(sid, "BRAF", cl, hist)]["mutated"].add(samp)
                rows[(sid, "BRAF", "BRAF_any", hist)]["mutated"].add(samp)
            elif gene in {"HRAS", "KRAS", "NRAS"}:
                cl = classify_ras(mut, gene) or f"{gene}_other"
                rows[(sid, gene, cl, hist)]["mutated"].add(samp)
                rows[(sid, gene, f"{gene}_any", hist)]["mutated"].add(samp)
            elif gene == "TP53":
                rows[(sid, "TP53", "TP53_any", hist)]["mutated"].add(samp)

    # Materialize TSV rows
    per_study: List[Dict[str, Any]] = []
    for (sid, gene, hs, hist), d in rows.items():
        n_test = len(d["tested"])
        n_mut = len(d["mutated"])
        if n_test == 0:
            continue
        per_study.append({
            "study": sid,
            "gene": gene,
            "hotspot": hs,
            "histology": hist,
            "n_mutated": n_mut,
            "n_tested": n_test,
            "prevalence": n_mut / n_test if n_test else float("nan"),
        })
    return {"per_study_rows": per_study}


def write_tsv(path: Path, rows: List[Dict[str, Any]],
              columns: List[str]) -> None:
    with path.open("w") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in columns) + "\n")


def build_supplementary_table(per_study_rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Pivot to histology x gene x hotspot rows, with study columns + meta-pool."""
    studies = sorted({r["study"] for r in per_study_rows})
    # Group key
    grp: Dict[tuple, Dict[str, Dict[str, int]]] = defaultdict(lambda: defaultdict(dict))
    for r in per_study_rows:
        key = (r["histology"], r["gene"], r["hotspot"])
        grp[key][r["study"]] = {"n_mut": r["n_mutated"], "n_test": r["n_tested"]}

    out_rows: List[Dict[str, Any]] = []
    for (hist, gene, hs), per in sorted(grp.items()):
        row: Dict[str, Any] = {"histology": hist, "gene": gene, "hotspot": hs}
        tot_mut = 0
        tot_test = 0
        for st in studies:
            v = per.get(st)
            if v:
                row[f"{st}__n_mut"] = v["n_mut"]
                row[f"{st}__n_test"] = v["n_test"]
                row[f"{st}__prev"] = (v["n_mut"] / v["n_test"]
                                      if v["n_test"] else "")
                tot_mut += v["n_mut"]
                tot_test += v["n_test"]
            else:
                row[f"{st}__n_mut"] = ""
                row[f"{st}__n_test"] = ""
                row[f"{st}__prev"] = ""
        row["pooled_n_mut"] = tot_mut
        row["pooled_n_test"] = tot_test
        row["pooled_prevalence"] = tot_mut / tot_test if tot_test else ""
        out_rows.append(row)

    cols = ["histology", "gene", "hotspot"]
    for st in studies:
        cols += [f"{st}__n_mut", f"{st}__n_test", f"{st}__prev"]
    cols += ["pooled_n_mut", "pooled_n_test", "pooled_prevalence"]
    return out_rows, cols


def fallback_summary(reason: str) -> Dict[str, Any]:
    return {
        "status": "FAILED",
        "error": reason,
        "fallback": {
            "tcga_thca_pub_TERT_promoter_n_mutated": 36,
            "source": "v17 TERT recovery v2 (cBioPortal thca_tcga_pub) -- "
                      "see project memory.",
        },
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }


def main() -> int:
    print("=== v17 ULTIMATE U3C: cBioPortal cross-cohort thyroid mutations ===")
    studies = list_thyroid_studies()
    if not studies:
        # Retry once
        print("First studies fetch returned empty; retrying once...")
        time.sleep(2.0)
        studies = list_thyroid_studies()
    if not studies:
        summary = fallback_summary("Could not list studies from cBioPortal API")
        (OUT_DIR / "U3C_summary.json").write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))
        return 1

    # Persist study list
    study_rows = [{
        "study_id": s["studyId"],
        "name": s.get("name", ""),
        "n_samples": s.get("allSampleCount", ""),
        "has_mutations": "",  # filled later
    } for s in studies]

    print(f"\n[2/5] Querying {len(studies)} studies for mutations...\n")
    results = []
    for st in studies:
        try:
            res = process_study(st)
        except Exception as exc:  # be tolerant; keep going
            print(f"   !! exception in {st['studyId']}: {exc}")
            res = {"studyId": st["studyId"], "name": st.get("name", ""),
                   "n_samples": st.get("allSampleCount", 0),
                   "has_mutations": False, "mutations": [], "samples": [],
                   "clinical": {}}
        results.append(res)
        for srow in study_rows:
            if srow["study_id"] == res["studyId"]:
                srow["has_mutations"] = res["has_mutations"]
                if res.get("samples"):
                    srow["n_samples"] = len(res["samples"])

    write_tsv(OUT_DIR / "U3C_cbioportal_studies_list.tsv", study_rows,
              ["study_id", "name", "n_samples", "has_mutations"])

    print("\n[3/5] Aggregating prevalence...")
    agg = aggregate(results)
    per_study_rows = agg["per_study_rows"]
    write_tsv(OUT_DIR / "U3C_per_study_mutation_counts.tsv", per_study_rows,
              ["study", "gene", "hotspot", "histology",
               "n_mutated", "n_tested", "prevalence"])

    print("[4/5] Building supplementary table S20...")
    sup_rows, sup_cols = build_supplementary_table(per_study_rows)
    write_tsv(OUT_DIR / "U3C_supplementary_table_S20.tsv", sup_rows, sup_cols)

    # ---- summary stats ----
    def prev_range(gene: str, hotspot_filter: Optional[Iterable[str]] = None):
        vals = []
        for r in per_study_rows:
            if r["gene"] != gene:
                continue
            if hotspot_filter and r["hotspot"] not in hotspot_filter:
                continue
            if r["n_tested"] >= 10:  # ignore tiny strata
                vals.append(r["prevalence"])
        if not vals:
            return None
        return {"min": min(vals), "max": max(vals), "n_strata": len(vals)}

    total_samples = sum(len(r["samples"]) for r in results if r.get("has_mutations"))
    summary = {
        "status": "OK",
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "n_studies": len(results),
        "n_studies_with_mutations": sum(1 for r in results if r.get("has_mutations")),
        "total_samples_queried": total_samples,
        "studies": [{"id": r["studyId"], "name": r["name"],
                     "n_samples": len(r["samples"]),
                     "n_mutations_pulled": len(r["mutations"])} for r in results],
        "TERT_any_prevalence_range": prev_range("TERT", {"TERT_any"}),
        "TERT_C228T_prevalence_range": prev_range("TERT", {"C228T"}),
        "TERT_C250T_prevalence_range": prev_range("TERT", {"C250T"}),
        "BRAF_V600E_prevalence_range": prev_range("BRAF", {"V600E"}),
        "RAS_hotspot_prevalence_range_HRAS": prev_range(
            "HRAS", {"HRAS_codon12", "HRAS_codon13", "HRAS_codon61"}),
        "RAS_hotspot_prevalence_range_KRAS": prev_range(
            "KRAS", {"KRAS_codon12", "KRAS_codon13", "KRAS_codon61"}),
        "RAS_hotspot_prevalence_range_NRAS": prev_range(
            "NRAS", {"NRAS_codon12", "NRAS_codon13", "NRAS_codon61"}),
        "TP53_prevalence_range": prev_range("TP53", {"TP53_any"}),
        "outputs": {
            "studies_list": str(OUT_DIR / "U3C_cbioportal_studies_list.tsv"),
            "per_study_counts": str(OUT_DIR / "U3C_per_study_mutation_counts.tsv"),
            "supplementary_S20": str(OUT_DIR / "U3C_supplementary_table_S20.tsv"),
        },
    }
    (OUT_DIR / "U3C_summary.json").write_text(json.dumps(summary, indent=2,
                                                          default=str))

    print("\n[5/5] Done.")
    print(json.dumps({k: v for k, v in summary.items()
                      if k not in {"studies"}}, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
