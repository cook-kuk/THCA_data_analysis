"""Deep analysis of lit_enrich raw JSON — dedupe, rank, cross-reference, and emit
hard tables for manuscript_v8 / paper2_brief paste.

Inputs:  project/manuscript_p2_brief/lit_enrich_2026_05_02/data/*.json
Outputs: project/manuscript_p2_brief/lit_enrich_2026_05_02/deep_analysis/*

Run:    project/.venv/bin/python project/notebooks_or_scripts/v17_lit_enrich_deep.py
"""

from __future__ import annotations

import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
DATA = ROOT / "project/manuscript_p2_brief/lit_enrich_2026_05_02/data"
OUT = ROOT / "project/manuscript_p2_brief/lit_enrich_2026_05_02/deep_analysis"
OUT.mkdir(parents=True, exist_ok=True)


# ----------------------------------------------------------------------------- #
# Load
# ----------------------------------------------------------------------------- #
def load_json(name: str):
    return json.loads((DATA / f"{name}.json").read_text())


verified = load_json("verified_refs")
related = load_json("related_works")
afnd = load_json("afnd_alleles")
trials_raw = load_json("trials")
unpaywall = load_json("unpaywall_links")
misattr = load_json("misattribution")

# Quick shape report
shapes = {
    "verified_refs (entries)": len(verified) if isinstance(verified, list) else len(verified),
    "related_works (claims)": len(related) if isinstance(related, list) else len(related),
    "afnd_alleles (queries)": len(afnd) if isinstance(afnd, list) else len(afnd),
    "trials_raw (queries)": len(trials_raw) if isinstance(trials_raw, list) else len(trials_raw),
    "unpaywall (entries)": len(unpaywall) if isinstance(unpaywall, list) else len(unpaywall),
    "misattribution": len(misattr) if isinstance(misattr, list) else len(misattr),
}
print("=== shapes ===")
for k, v in shapes.items():
    print(f"  {k}: {v}")


# ----------------------------------------------------------------------------- #
# 1. verified_refs — triangulation status + completeness map
# ----------------------------------------------------------------------------- #
def normalize_doi(s: str | None) -> str | None:
    if not s:
        return None
    s = s.strip().lower()
    s = re.sub(r"^https?://(dx\.)?doi\.org/", "", s)
    return s or None


def reflect_verified(verified):
    """verified is a list of {key, incomplete, crossref_match, crossref_title, [pubmed_*], ...}"""
    rows = []
    for entry in verified:
        if not isinstance(entry, dict):
            continue
        key = entry.get("key")
        # Detect which source confirmed
        crossref_ok = entry.get("crossref_match")
        pubmed_ok = bool(entry.get("pubmed_match") or entry.get("pubmed_pmid"))
        openalex_ok = bool(entry.get("openalex_id") or entry.get("openalex_match"))
        epmc_ok = bool(entry.get("europepmc_id") or entry.get("epmc_match"))
        n_confirm = sum([bool(crossref_ok), pubmed_ok, openalex_ok, epmc_ok])
        rows.append(
            {
                "key": key,
                "incomplete": entry.get("incomplete"),
                "crossref": bool(crossref_ok),
                "pubmed": pubmed_ok,
                "openalex": openalex_ok,
                "europepmc": epmc_ok,
                "n_confirm": n_confirm,
                "crossref_title": (entry.get("crossref_title") or "")[:100],
                "pubmed_pmid": entry.get("pubmed_pmid"),
                "openalex_id": entry.get("openalex_id"),
                "doi": normalize_doi(entry.get("doi") or entry.get("crossref_doi")),
            }
        )
    df = pd.DataFrame(rows).sort_values(["incomplete", "n_confirm"], ascending=[False, False])
    return df


df_verify = reflect_verified(verified)
df_verify.to_csv(OUT / "01_verified_refs_triangulation.csv", index=False)
print(f"\n[01] verified_refs triangulation -> {len(df_verify)} entries")
print(df_verify.head(20).to_string(index=False))


# ----------------------------------------------------------------------------- #
# 2. related_works — cross-source dedupe + relevance rank per claim
# ----------------------------------------------------------------------------- #
def flatten_related(related: dict) -> pd.DataFrame:
    """Each claim has multiple sources (openalex/sem_scholar/europe_pmc/biorxiv/arxiv).
    Flatten to long format with one row per (claim, source, paper).
    """
    rows = []
    for claim, src_dict in related.items():
        if not isinstance(src_dict, dict):
            continue
        for src, hits in src_dict.items():
            if not isinstance(hits, list):
                continue
            for h in hits:
                if not isinstance(h, dict):
                    continue
                rows.append(
                    {
                        "claim": claim,
                        "source": src,
                        "doi": normalize_doi(h.get("doi")),
                        "title": (h.get("title") or "").strip(),
                        "year": h.get("year"),
                        "venue": h.get("venue") or h.get("journal"),
                        "first_author": h.get("first_author") or h.get("first") or h.get("author"),
                        "cited_by": h.get("cited_by") or h.get("cited") or h.get("citationCount"),
                        "influential": h.get("influential")
                        or h.get("influentialCitationCount"),
                    }
                )
    return pd.DataFrame(rows)


df_rel = flatten_related(related)
print(f"\n[02] related_works flat: {len(df_rel)} rows × source/claim")
print(df_rel.groupby(["claim", "source"]).size().unstack(fill_value=0))


def dedupe_rank(df: pd.DataFrame) -> pd.DataFrame:
    """For each claim, group identical papers across sources; rank by (cited_by desc, influential desc, year desc)."""
    df = df.copy()
    # dedup key: prefer DOI, fall back to (lowercased title prefix, year)
    df["dedup_key"] = df.apply(
        lambda r: r["doi"]
        if r["doi"]
        else (((r["title"] or "")[:60]).lower() + f"|{r['year'] or ''}"),
        axis=1,
    )
    df["cited_by"] = pd.to_numeric(df["cited_by"], errors="coerce").fillna(0)
    df["influential"] = pd.to_numeric(df["influential"], errors="coerce").fillna(0)
    df["year"] = pd.to_numeric(df["year"], errors="coerce")

    grouped = (
        df.groupby(["claim", "dedup_key"])
        .agg(
            sources=("source", lambda s: ",".join(sorted(set(s)))),
            n_sources=("source", "nunique"),
            doi=("doi", "first"),
            title=("title", lambda s: max(s, key=len) if len(s) else ""),
            year=("year", "max"),
            venue=("venue", lambda s: next((x for x in s if x), None)),
            first_author=("first_author", lambda s: next((x for x in s if x), None)),
            cited_by=("cited_by", "max"),
            influential=("influential", "max"),
        )
        .reset_index()
    )
    grouped["composite"] = (
        grouped["n_sources"] * 100  # cross-source corroboration
        + grouped["cited_by"].clip(upper=500)  # cap to avoid mega-citation skew
        + grouped["influential"] * 5
        + (grouped["year"].fillna(2000) - 2000)
    )
    return grouped.sort_values(["claim", "composite"], ascending=[True, False])


df_dedup = dedupe_rank(df_rel)
df_dedup.to_csv(OUT / "02_related_works_dedup_ranked.csv", index=False)

# Per-claim top 5 dedup table
top5 = (
    df_dedup.groupby("claim").head(5).reset_index(drop=True)
    [["claim", "n_sources", "year", "first_author", "title", "cited_by", "influential", "doi"]]
)
top5.to_csv(OUT / "02_related_works_top5_per_claim.csv", index=False)
print(f"\n[02] dedupe -> {len(df_dedup)} unique papers across {df_dedup['claim'].nunique()} claims")
print(f"     top-5 per claim emitted ({len(top5)} rows)")


# ----------------------------------------------------------------------------- #
# 3. trials.json — structured extraction
# ----------------------------------------------------------------------------- #
def flatten_trials(trials_raw) -> pd.DataFrame:
    rows = []
    container = trials_raw if isinstance(trials_raw, dict) else {"_": trials_raw}
    for query, payload in container.items():
        # payload may be dict with 'studies' (CT.gov v2 format)
        studies = []
        if isinstance(payload, dict):
            studies = payload.get("studies") or payload.get("results") or []
        elif isinstance(payload, list):
            studies = payload
        for st in studies:
            if not isinstance(st, dict):
                continue
            proto = st.get("protocolSection", {}) if "protocolSection" in st else st
            ident = proto.get("identificationModule", {})
            status = proto.get("statusModule", {})
            design = proto.get("designModule", {})
            sponsor = proto.get("sponsorCollaboratorsModule", {})
            arms = proto.get("armsInterventionsModule", {})
            cond = proto.get("conditionsModule", {})
            outcomes = proto.get("outcomesModule", {})
            elig = proto.get("eligibilityModule", {})

            interventions = arms.get("interventions", [])
            int_types = [i.get("type") for i in interventions if isinstance(i, dict)]
            int_names = [i.get("name") for i in interventions if isinstance(i, dict)]

            phases = design.get("phases", [])
            enroll = design.get("enrollmentInfo", {}).get("count")

            rows.append(
                {
                    "query": query,
                    "nct_id": ident.get("nctId"),
                    "title": ident.get("briefTitle") or ident.get("officialTitle"),
                    "status": status.get("overallStatus"),
                    "start_date": (status.get("startDateStruct") or {}).get("date"),
                    "phase": ",".join(phases) if phases else None,
                    "enrollment": enroll,
                    "study_type": design.get("studyType"),
                    "intervention_types": ",".join(sorted(set(t for t in int_types if t))),
                    "intervention_names": "; ".join(int_names) if int_names else None,
                    "conditions": ",".join(cond.get("conditions", [])),
                    "lead_sponsor": (sponsor.get("leadSponsor") or {}).get("name"),
                    "sponsor_class": (sponsor.get("leadSponsor") or {}).get("class"),
                    "primary_outcome": "; ".join(
                        o.get("measure") for o in outcomes.get("primaryOutcomes", []) if isinstance(o, dict) and o.get("measure")
                    ),
                    "min_age": elig.get("minimumAge"),
                    "sex": elig.get("sex"),
                }
            )
    return pd.DataFrame(rows)


df_tr = flatten_trials(trials_raw)
df_tr.to_csv(OUT / "03_trials_structured.csv", index=False)
print(f"\n[03] trials structured: {len(df_tr)} trials across {df_tr['query'].nunique()} queries")
print("     intervention_types distribution:")
print(df_tr["intervention_types"].value_counts().head(10).to_string())
print("     phase distribution:")
print(df_tr["phase"].value_counts().head(10).to_string())

# Targeted slices: BRAF V600E + immunotherapy (relevant to dark-matter Discussion)
sliced = df_tr[
    df_tr["intervention_names"].fillna("").str.contains(
        r"selpercatinib|pralsetinib|larotrectinib|dabrafenib|trametinib|encorafenib|vemurafenib|cobimetinib|pembrolizumab|nivolumab|atezolizumab|cemiplimab|lenvatinib|cabozantinib",
        case=False,
        regex=True,
    )
]
sliced.to_csv(OUT / "03_trials_targeted_drugs.csv", index=False)
print(f"     targeted-drug trials (RET/BRAF/IO/TKI): {len(sliced)}")


# ----------------------------------------------------------------------------- #
# 4. AFND HLA — full 5×4 matrix + forest-ready table
# ----------------------------------------------------------------------------- #
def flatten_afnd(afnd) -> pd.DataFrame:
    rows = []
    container = afnd if isinstance(afnd, dict) else {}
    for allele, pop_dict in container.items():
        if not isinstance(pop_dict, dict):
            continue
        for pop, hits in pop_dict.items():
            if isinstance(hits, list):
                for h in hits:
                    if not isinstance(h, dict):
                        continue
                    rows.append(
                        {
                            "allele": allele,
                            "population": pop,
                            "pop_name": h.get("population") or h.get("pop_name"),
                            "n": h.get("sample_size") or h.get("n"),
                            "frequency": h.get("allele_freq") or h.get("frequency") or h.get("af"),
                            "country": h.get("country") or pop,
                        }
                    )
            elif isinstance(hits, dict):
                # alt schema
                for k, v in hits.items():
                    if isinstance(v, dict):
                        rows.append(
                            {
                                "allele": allele,
                                "population": pop,
                                "pop_name": k,
                                "n": v.get("n"),
                                "frequency": v.get("frequency"),
                                "country": v.get("country"),
                            }
                        )
    return pd.DataFrame(rows)


df_hla = flatten_afnd(afnd)
df_hla["frequency"] = pd.to_numeric(df_hla["frequency"], errors="coerce")
df_hla["n"] = pd.to_numeric(df_hla["n"], errors="coerce")
df_hla.to_csv(OUT / "04_hla_long.csv", index=False)
print(f"\n[04] HLA AFND long: {len(df_hla)} (allele × pop × study)")

if len(df_hla):
    # Sample-weighted mean frequency per (allele, population)
    weighted = (
        df_hla.dropna(subset=["frequency", "n"])
        .groupby(["allele", "population"])
        .apply(lambda g: pd.Series({
            "weighted_freq": np.average(g["frequency"], weights=g["n"]) if g["n"].sum() else g["frequency"].mean(),
            "total_n": int(g["n"].sum()),
            "n_studies": len(g),
        }), include_groups=False)
        .reset_index()
    )
    weighted.to_csv(OUT / "04_hla_pop_weighted.csv", index=False)
    print("     weighted (allele × pop):")
    print(weighted.to_string(index=False))


# ----------------------------------------------------------------------------- #
# 5. Unpaywall — fetch actual abstracts/metadata for OA links
# ----------------------------------------------------------------------------- #
def reflect_unpaywall(unpaywall):
    """unpaywall is a list of {key, doi, oa_status, is_oa, pdf_url, landing_url, license, europepmc_url, pmcid}"""
    rows = []
    for entry in unpaywall:
        if not isinstance(entry, dict):
            continue
        rows.append(
            {
                "key": entry.get("key"),
                "doi": normalize_doi(entry.get("doi")),
                "is_oa": entry.get("is_oa"),
                "oa_status": entry.get("oa_status"),
                "pdf_url": entry.get("pdf_url"),
                "landing_url": entry.get("landing_url"),
                "license": entry.get("license"),
                "pmcid": entry.get("pmcid"),
                "europepmc_url": entry.get("europepmc_url"),
            }
        )
    return pd.DataFrame(rows)


df_up = reflect_unpaywall(unpaywall)
df_up.to_csv(OUT / "05_unpaywall_summary.csv", index=False)
print(f"\n[05] unpaywall summary: {len(df_up)} entries, {df_up['is_oa'].sum() if 'is_oa' in df_up else 0} OA")


# ----------------------------------------------------------------------------- #
# 6. Cross-reference: which dedup top-papers also appear in our verified bib?
# ----------------------------------------------------------------------------- #
verified_dois = set(df_verify["doi"].dropna()) | set(df_up["doi"].dropna())
df_dedup["in_our_bib"] = df_dedup["doi"].isin(verified_dois)
overlap = df_dedup[df_dedup["in_our_bib"]][["claim", "doi", "title", "first_author", "year"]]
overlap.to_csv(OUT / "06_overlap_with_bib.csv", index=False)
print(f"\n[06] dedup ∩ verified_bib: {len(overlap)} papers (both in our bib AND mined)")

# Conversely: top-ranked papers per claim NOT in bib (these are gap candidates)
gap = df_dedup[~df_dedup["in_our_bib"]].groupby("claim").head(3)
gap[["claim", "n_sources", "year", "first_author", "title", "cited_by", "doi"]].to_csv(
    OUT / "06_gap_candidates_top3_per_claim.csv", index=False
)
print(f"     gap candidates (top-3 per claim, not in bib): {len(gap)}")

print("\n=== DONE ===")
print(f"Outputs in: {OUT}")
for p in sorted(OUT.glob("*")):
    print(f"  {p.name}  ({p.stat().st_size/1024:.1f} KB)")
