"""
Six lit_enrich tasks. Each writes a markdown report + a JSON data file.
All paths are relative to an output directory passed by the caller.
"""
from __future__ import annotations

import json
import re
import time
from pathlib import Path
from typing import Any

import bibtexparser

from . import sources as S


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _write(out_dir: Path, name: str, content: str) -> Path:
    p = out_dir / name
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(content, encoding="utf-8")
    return p


def _dump_json(out_dir: Path, name: str, obj: Any) -> Path:
    p = out_dir / "data" / name
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")
    return p


def _norm_authors(record: dict) -> str:
    """Extract first-author surname from heterogeneous record shapes."""
    if not record:
        return ""
    # OpenAlex
    a = record.get("authorships") or []
    if a:
        name = (a[0].get("author") or {}).get("display_name", "")
        if name:
            return name.split()[-1]
    # CrossRef
    a = record.get("author") or []
    if a:
        return a[0].get("family") or a[0].get("name") or ""
    # PubMed esummary
    a = record.get("authors") or []
    if a and isinstance(a, list):
        first = a[0]
        if isinstance(first, dict):
            return first.get("name", "") or first.get("display_name", "")
        return str(first).split()[-1]
    # Semantic Scholar
    a = record.get("authors") or []
    if a and isinstance(a[0], dict):
        return (a[0].get("name") or "").split()[-1]
    return ""


def _entry_query(entry: dict) -> str:
    """Build a free-text query from a bib entry, mining hints from notes
    when the title field is a placeholder.
    """
    parts: list[str] = []
    title = entry.get("title", "").strip("{}").strip()
    if title:
        if "TBD" in title or "verify" in title.lower() or "본인" in title:
            # Skeleton title — extract just the topical prefix before the
            # first delimiter (— or full citation).
            head = re.split(r"[—\-]|full citation|본인", title, maxsplit=1)[0].strip()
            head = head.strip("{}").strip(":,. ")
            if head and len(head) > 4:
                parts.append(head)
        else:
            parts.append(title)

    author = entry.get("author", "")
    if author:
        first = author.split(" and ")[0].split(",")[0].strip()
        if first and "others" not in first.lower() and "{" not in first:
            parts.append(first)

    year = entry.get("year")
    if year:
        parts.append(str(year))

    note = entry.get("note", "")
    # 1) GEO / PMC / PRJ / SRA accessions are gold for finding the right paper
    for m in re.finditer(r"(GSE\d+|PRJ[A-Z]+\d+|SRP\d+|PMC\d+)", note):
        parts.append(m.group(1))
    # 2) Author-Year hints in the note ("Chu 2018", "Davies & Welch 2014")
    for m in re.finditer(r"([A-Z][a-z]+(?:\s*&\s*[A-Z][a-z]+)?\s+\d{4})", note):
        token = m.group(1)
        # Avoid duplicating the bib's own author/year
        if not (first and token.startswith(first)):
            parts.append(token)
    # 3) Journal hints ("JAMA Oncol", "Endocr Connect", "Cancer Cell")
    for m in re.finditer(r"\*([A-Z][\w\s]{2,30})\*", note):
        parts.append(m.group(1).strip())
    # 4) Topic keywords from "context: …"
    cm = re.search(r"context:\s*([A-Za-z0-9 ,]+?)(?:[?.\n]|$)", note)
    if cm:
        parts.append(cm.group(1).strip())

    return " ".join(parts)


def _load_bib(bib_path: Path) -> list[dict]:
    if not bib_path.exists():
        return []
    parser = bibtexparser.bparser.BibTexParser(common_strings=True)
    parser.ignore_nonstandard_types = False
    db = bibtexparser.load(bib_path.open(encoding="utf-8"), parser=parser)
    return db.entries


# ─────────────────────────────────────────────────────────────────────────────
# Task 1 — verify_references
# ─────────────────────────────────────────────────────────────────────────────

def verify_references(bib_path: Path, out_dir: Path) -> dict:
    entries = _load_bib(bib_path)
    oa = S.OpenAlex()
    pm = S.PubMed()
    cr = S.CrossRef()

    results: list[dict] = []
    md_lines = ["# 01 Reference verification\n",
                f"_Source bib: `{bib_path}` ({len(entries)} entries)_\n",
                "Triangulation across OpenAlex + PubMed + CrossRef. ",
                "Voice-protected: this is a report only — paste verified DOI/title/journal "
                "into the bib by hand.\n"]

    for e in entries:
        key = e.get("ID", "?")
        has_doi = bool(e.get("doi"))
        is_incomplete = (
            not has_doi
            or "TBD" in (e.get("journal") or "")
            or "verify" in (e.get("title") or "").lower()
            or "others" in (e.get("author") or "").lower()
        )
        row: dict[str, Any] = {"key": key, "incomplete": is_incomplete}

        if has_doi and not is_incomplete:
            # quick DOI sanity check via CrossRef
            cr_doc = cr.lookup(e["doi"])
            row["crossref_match"] = bool(cr_doc)
            row["crossref_title"] = cr_doc.get("title", [""])[0] if cr_doc else ""
            results.append(row)
            continue

        # incomplete → search
        q = _entry_query(e)
        row["query"] = q
        row["candidates"] = []

        oa_hits = oa.search(q, limit=3)
        time.sleep(0.2)
        pm_hits = pm.search_with_summary(q, limit=3)
        time.sleep(0.4)
        cr_hits = cr.search(q, limit=3)

        for h in oa_hits:
            row["candidates"].append({
                "source": "OpenAlex",
                "title": h.get("title", ""),
                "doi": (h.get("doi") or "").replace("https://doi.org/", ""),
                "journal": ((h.get("primary_location") or {}).get("source") or {}).get("display_name", ""),
                "year": h.get("publication_year"),
                "first_author": _norm_authors(h),
                "openalex_id": h.get("id", ""),
            })
        for h in pm_hits:
            row["candidates"].append({
                "source": "PubMed",
                "title": h.get("title", ""),
                "doi": next((x.get("value") for x in (h.get("articleids") or [])
                             if x.get("idtype") == "doi"), ""),
                "journal": h.get("source", ""),
                "year": (h.get("pubdate", "") or "")[:4],
                "first_author": _norm_authors(h),
                "pmid": h.get("uid", ""),
            })
        for h in cr_hits:
            row["candidates"].append({
                "source": "CrossRef",
                "title": (h.get("title") or [""])[0],
                "doi": h.get("DOI", ""),
                "journal": (h.get("container-title") or [""])[0],
                "year": ((h.get("issued") or {}).get("date-parts") or [[None]])[0][0],
                "first_author": _norm_authors(h),
            })
        results.append(row)

    # Markdown rendering
    incomplete = [r for r in results if r.get("incomplete")]
    md_lines.append(f"\n## Summary\n\n- Entries: {len(entries)}")
    md_lines.append(f"- Incomplete: **{len(incomplete)}**")
    md_lines.append(f"- Verified-with-DOI passing CrossRef sanity check: "
                    f"{sum(1 for r in results if r.get('crossref_match'))}\n")

    if incomplete:
        md_lines.append("\n## Incomplete entries — verified suggestions\n")
        for r in incomplete:
            md_lines.append(f"\n### `{r['key']}`\n")
            md_lines.append(f"_Query: `{r.get('query', '')}`_\n")
            cands = r.get("candidates", [])
            if not cands:
                md_lines.append("⚠️ No candidates found across OpenAlex / PubMed / CrossRef.\n")
                continue
            md_lines.append("| Src | First author | Year | Journal | Title | DOI |")
            md_lines.append("|---|---|---|---|---|---|")
            for c in cands[:9]:
                title = (c.get("title") or "")[:80].replace("|", "/")
                md_lines.append(
                    f"| {c['source']} | {c.get('first_author','')} | {c.get('year','')} "
                    f"| {(c.get('journal') or '')[:40]} | {title} | "
                    f"`{c.get('doi','')}` |"
                )
    _write(out_dir, "01_reference_verification.md", "\n".join(md_lines))
    _dump_json(out_dir, "verified_refs.json", results)
    return {"entries": len(entries), "incomplete": len(incomplete), "results": results}


# ─────────────────────────────────────────────────────────────────────────────
# Task 2 — competitive_landscape
# ─────────────────────────────────────────────────────────────────────────────

def competitive_landscape(claims: list[dict], out_dir: Path) -> dict:
    """
    claims: list of {topic, query, year_from, year_to}
    For each claim, hit OpenAlex + Semantic Scholar + Europe PMC + bioRxiv,
    rank by influential citations + recency.
    """
    oa = S.OpenAlex()
    ss = S.SemanticScholar()
    epmc = S.EuropePMC()
    bio = S.BioRxiv()
    arx = S.Arxiv()

    md_lines = ["# 02 Competitive landscape\n",
                "Mining OpenAlex + Semantic Scholar + Europe PMC + bioRxiv + arXiv "
                "for missing or competing citations on each main claim.\n"]
    out: dict[str, Any] = {}

    for claim in claims:
        topic = claim["topic"]
        query = claim["query"]
        yf = claim.get("year_from", 2018)
        md_lines.append(f"\n## {topic}\n")
        md_lines.append(f"_Query: `{query}` (since {yf})_\n")

        oa_hits = oa.search(query, limit=15, year_from=yf)
        time.sleep(0.3)
        ss_hits = ss.search(query, limit=15)
        time.sleep(1.0)  # SS is slow without key
        epmc_hits = epmc.search(query, limit=10, full_text_only=True)
        time.sleep(0.3)
        bio_hits = bio.search(query, limit=10)
        time.sleep(0.3)
        arx_hits = arx.search(query, limit=8)
        time.sleep(0.3)

        # Influential SS results first (sort by influentialCitationCount)
        ss_ranked = sorted(ss_hits, key=lambda x: x.get("influentialCitationCount", 0) or 0,
                           reverse=True)[:8]

        out[topic] = {
            "openalex": [{"title": h.get("title"), "year": h.get("publication_year"),
                          "doi": (h.get("doi") or "").replace("https://doi.org/", ""),
                          "cited_by": h.get("cited_by_count"),
                          "first_author": _norm_authors(h)} for h in oa_hits],
            "semantic_scholar_influential": [{
                "title": h.get("title"), "year": h.get("year"),
                "venue": h.get("venue"),
                "influential_cites": h.get("influentialCitationCount", 0),
                "total_cites": h.get("citationCount", 0),
                "doi": (h.get("externalIds") or {}).get("DOI", ""),
                "first_author": _norm_authors(h),
            } for h in ss_ranked],
            "europepmc_fulltext": [{"title": h.get("title"), "year": h.get("pubYear"),
                                    "journal": h.get("journalTitle"), "doi": h.get("doi"),
                                    "pmcid": h.get("pmcid")} for h in epmc_hits],
            "biorxiv_preprints": [{"title": h.get("title"), "year": h.get("pubYear"),
                                   "doi": h.get("doi")} for h in bio_hits],
            "arxiv": [{"title": h.get("title"), "published": h.get("published"),
                       "id": h.get("id"), "first_author": (h.get("authors") or [""])[0]}
                      for h in arx_hits],
        }

        # Render
        md_lines.append("\n**Top 8 (Semantic Scholar influential):**\n")
        if ss_ranked:
            md_lines.append("| Influential / Total | Year | Venue | First | Title | DOI |")
            md_lines.append("|---|---|---|---|---|---|")
            for h in ss_ranked:
                md_lines.append(
                    f"| {h.get('influentialCitationCount', 0)} / "
                    f"{h.get('citationCount', 0)} | {h.get('year','')} | "
                    f"{(h.get('venue') or '')[:30]} | "
                    f"{_norm_authors(h)} | "
                    f"{(h.get('title') or '')[:80].replace('|','/')} | "
                    f"`{(h.get('externalIds') or {}).get('DOI','')}` |"
                )
        else:
            md_lines.append("_(no Semantic Scholar hits — possibly rate-limited)_\n")

        md_lines.append("\n**OpenAlex top 10 (recency × cited):**\n")
        oa_ranked = sorted(oa_hits, key=lambda x: (x.get("publication_year") or 0,
                                                    x.get("cited_by_count") or 0),
                           reverse=True)[:10]
        if oa_ranked:
            md_lines.append("| Cited | Year | First | Title | DOI |")
            md_lines.append("|---|---|---|---|---|")
            for h in oa_ranked:
                md_lines.append(
                    f"| {h.get('cited_by_count', 0)} | "
                    f"{h.get('publication_year','')} | "
                    f"{_norm_authors(h)} | "
                    f"{(h.get('title') or '')[:80].replace('|','/')} | "
                    f"`{(h.get('doi') or '').replace('https://doi.org/','')}` |"
                )

        md_lines.append("\n**bioRxiv / medRxiv preprints (≥2024):**\n")
        bio_2024 = [h for h in bio_hits if (h.get("pubYear") or "") >= "2024"]
        if bio_2024:
            md_lines.append("| Year | Title | DOI |")
            md_lines.append("|---|---|---|")
            for h in bio_2024[:8]:
                md_lines.append(
                    f"| {h.get('pubYear','')} | "
                    f"{(h.get('title') or '')[:90].replace('|','/')} | "
                    f"`{h.get('doi','')}` |"
                )
        else:
            md_lines.append("_(no recent preprints found)_\n")

    _write(out_dir, "02_competitive_landscape.md", "\n".join(md_lines))
    _dump_json(out_dir, "related_works.json", out)
    return {"claims": len(claims), "results": out}


# ─────────────────────────────────────────────────────────────────────────────
# Task 3 — hla_frequencies (AFND)
# ─────────────────────────────────────────────────────────────────────────────

def hla_frequencies(alleles: list[str], populations: list[str], out_dir: Path) -> dict:
    afnd = S.AFND()
    md_lines = ["# 03 HLA allele frequencies (AFND)\n",
                "Direct lookup of risk alleles in Korean / Chinese / Japanese populations "
                "for Paper 2 Pillar 1 forest cross-check.\n"]
    out: dict[str, Any] = {}

    for allele in alleles:
        md_lines.append(f"\n## {allele}\n")
        out[allele] = {}
        for pop in populations:
            rows = afnd.allele_frequency(allele, population_filter=pop, limit=20)
            time.sleep(1.0)  # be polite to AFND scrape
            out[allele][pop] = rows
            md_lines.append(f"\n### {pop} (n={len(rows)} populations matched)\n")
            if rows:
                md_lines.append("| Population | Allele freq | Sample size |")
                md_lines.append("|---|---|---|")
                for r in rows[:10]:
                    md_lines.append(f"| {r['population']} | {r['allele_freq']} | {r['sample_size']} |")
            else:
                md_lines.append("_(no rows — AFND may return empty for indirect filters; verify manually)_\n")

    md_lines.append("\n---\n\n_AFND scrape is heuristic. Cross-check the populations of "
                    "interest at http://www.allelefrequencies.net before citing._\n")
    _write(out_dir, "03_hla_frequencies.md", "\n".join(md_lines))
    _dump_json(out_dir, "afnd_alleles.json", out)
    return {"alleles": len(alleles), "populations": len(populations), "results": out}


# ─────────────────────────────────────────────────────────────────────────────
# Task 4 — misattribution_check
# ─────────────────────────────────────────────────────────────────────────────

def misattribution_check(bib_path: Path, out_dir: Path) -> dict:
    """
    For every entry with a DOI: pull CrossRef + OpenAlex + PubMed metadata
    and flag mismatches in (year, first-author surname, journal).
    """
    entries = _load_bib(bib_path)
    cr = S.CrossRef()
    oa = S.OpenAlex()

    md_lines = ["# 04 Misattribution check\n",
                "DOI fingerprint cross-check (year, first-author surname, journal).\n"
                "Surfaces problems like the Krishnamoorthy / Landa misattribution found earlier.\n"]
    issues: list[dict] = []

    for e in entries:
        doi = e.get("doi")
        if not doi:
            continue
        cr_doc = cr.lookup(doi)
        time.sleep(0.2)
        oa_doc = oa.get_work(doi=doi)
        time.sleep(0.2)

        bib_year = str(e.get("year", ""))
        bib_journal = (e.get("journal") or "").lower()
        bib_first_author = ""
        if e.get("author"):
            bib_first_author = e["author"].split(" and ")[0].split(",")[0].strip()

        cr_year = str(((cr_doc.get("issued") or {}).get("date-parts") or [[""]])[0][0]) if cr_doc else ""
        cr_journal = (cr_doc.get("container-title") or [""])[0].lower() if cr_doc else ""
        cr_first = (cr_doc.get("author") or [{}])[0].get("family", "") if cr_doc else ""

        oa_year = str(oa_doc.get("publication_year", "")) if oa_doc else ""
        oa_journal = ((oa_doc.get("primary_location") or {}).get("source") or {}).get("display_name", "").lower()
        oa_first = _norm_authors(oa_doc) if oa_doc else ""

        flags = []
        if cr_year and bib_year and cr_year != bib_year:
            flags.append(f"year mismatch: bib={bib_year}, CrossRef={cr_year}")
        if cr_first and bib_first_author and cr_first.lower() not in bib_first_author.lower() \
                and bib_first_author.lower() not in cr_first.lower():
            flags.append(f"first author mismatch: bib={bib_first_author}, CrossRef={cr_first}")
        if cr_journal and bib_journal and cr_journal not in bib_journal and bib_journal not in cr_journal:
            # avoid false positives for abbreviations
            if not (cr_journal.startswith(bib_journal[:8]) or bib_journal.startswith(cr_journal[:8])):
                flags.append(f"journal mismatch: bib='{bib_journal}', CrossRef='{cr_journal}'")
        if oa_doc and oa_year and bib_year and oa_year != bib_year:
            flags.append(f"OpenAlex year={oa_year} vs bib={bib_year}")

        if flags:
            issues.append({
                "key": e.get("ID"), "doi": doi, "flags": flags,
                "bib": {"year": bib_year, "first_author": bib_first_author, "journal": bib_journal},
                "crossref": {"year": cr_year, "first_author": cr_first, "journal": cr_journal},
                "openalex": {"year": oa_year, "first_author": oa_first, "journal": oa_journal},
            })

    if issues:
        md_lines.append(f"\n## ⚠️ {len(issues)} entries with mismatches\n")
        for it in issues:
            md_lines.append(f"\n### `{it['key']}`  (DOI: `{it['doi']}`)\n")
            for f in it["flags"]:
                md_lines.append(f"- {f}")
            md_lines.append("")
            md_lines.append(f"  - bib       : {it['bib']}")
            md_lines.append(f"  - CrossRef  : {it['crossref']}")
            md_lines.append(f"  - OpenAlex  : {it['openalex']}")
    else:
        md_lines.append("\n✅ No DOI fingerprint mismatches found across CrossRef + OpenAlex.\n")

    _write(out_dir, "04_misattribution_check.md", "\n".join(md_lines))
    _dump_json(out_dir, "misattribution.json", issues)
    return {"checked": sum(1 for e in entries if e.get("doi")), "issues": len(issues)}


# ─────────────────────────────────────────────────────────────────────────────
# Task 5 — clinical_landscape
# ─────────────────────────────────────────────────────────────────────────────

def clinical_landscape(queries: list[str], out_dir: Path) -> dict:
    ct = S.ClinicalTrials()
    md_lines = ["# 05 Clinical trial landscape (ClinicalTrials.gov)\n",
                "Active and recruiting trials relevant to the manuscript's claims — "
                "useful for Discussion future-direction citations.\n"]
    out: dict[str, list] = {}

    for q in queries:
        rows = ct.search(q, limit=30)
        time.sleep(0.5)
        out[q] = rows
        md_lines.append(f"\n## `{q}` ({len(rows)} active/recruiting)\n")
        if not rows:
            md_lines.append("_(no matches)_\n")
            continue
        md_lines.append("| NCT ID | Phase | Status | Title | Sponsor |")
        md_lines.append("|---|---|---|---|---|")
        for s in rows[:20]:
            ps = s.get("protocolSection", {}) or {}
            ident = ps.get("identificationModule", {}) or {}
            status = ps.get("statusModule", {}) or {}
            design = ps.get("designModule", {}) or {}
            sponsor = (ps.get("sponsorCollaboratorsModule", {}) or {}).get("leadSponsor", {}) or {}
            md_lines.append(
                f"| {ident.get('nctId', '')} | "
                f"{','.join(design.get('phases', []) or [])} | "
                f"{status.get('overallStatus', '')} | "
                f"{(ident.get('briefTitle', '') or '')[:80].replace('|', '/')} | "
                f"{(sponsor.get('name', '') or '')[:30]} |"
            )
    _write(out_dir, "05_clinical_landscape.md", "\n".join(md_lines))
    _dump_json(out_dir, "trials.json", out)
    return {"queries": len(queries), "results_per_query": {q: len(v) for q, v in out.items()}}


# ─────────────────────────────────────────────────────────────────────────────
# Task 6 — full_text_discovery
# ─────────────────────────────────────────────────────────────────────────────

def full_text_discovery(bib_path: Path, out_dir: Path) -> dict:
    entries = _load_bib(bib_path)
    up = S.Unpaywall()
    epmc = S.EuropePMC()
    md_lines = ["# 06 Full-text discovery (Unpaywall + Europe PMC)\n",
                "Open-access PDF / full-text URLs for every DOI in the bib.\n"]
    rows: list[dict] = []

    for e in entries:
        doi = e.get("doi")
        if not doi:
            rows.append({"key": e.get("ID"), "doi": None, "oa_status": "no_doi"})
            continue
        u = up.lookup(doi)
        time.sleep(0.3)
        oa_loc = (u or {}).get("best_oa_location") or {}
        pmcid = (e.get("pmcid") or "")
        if not pmcid:
            # try Europe PMC search by DOI
            hits = epmc.search(f'DOI:"{doi}"', limit=1)
            if hits:
                pmcid = hits[0].get("pmcid", "")
            time.sleep(0.3)

        rows.append({
            "key": e.get("ID"),
            "doi": doi,
            "oa_status": (u or {}).get("oa_status", "unknown"),
            "is_oa": (u or {}).get("is_oa", False),
            "pdf_url": oa_loc.get("url_for_pdf") or "",
            "landing_url": oa_loc.get("url_for_landing_page") or "",
            "license": oa_loc.get("license") or "",
            "europepmc_url": epmc.get_full_text_url(pmcid) if pmcid else "",
            "pmcid": pmcid,
        })

    md_lines.append(f"\n## {len(rows)} entries · "
                    f"{sum(1 for r in rows if r.get('is_oa'))} open access\n")
    md_lines.append("\n| Key | OA | DOI | PDF | Europe PMC |")
    md_lines.append("|---|---|---|---|---|")
    for r in rows:
        oa_mark = "✅" if r.get("is_oa") else ("—" if r.get("doi") else "?")
        pdf = r.get("pdf_url", "")
        epmc_url = r.get("europepmc_url") or ""
        md_lines.append(
            f"| `{r['key']}` | {oa_mark} ({r.get('oa_status','')}) | "
            f"`{r.get('doi','')}` | "
            f"{('[PDF]('+pdf+')') if pdf else '—'} | "
            f"{('[PMC]('+epmc_url+')') if epmc_url else '—'} |"
        )
    _write(out_dir, "06_full_text_discovery.md", "\n".join(md_lines))
    _dump_json(out_dir, "unpaywall_links.json", rows)
    return {"entries": len(rows), "oa_count": sum(1 for r in rows if r.get('is_oa'))}


# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────

def write_summary(out_dir: Path, results: dict) -> Path:
    lines = [
        "# 07 lit_enrich summary\n",
        f"_Generated: {time.strftime('%Y-%m-%d %H:%M')}_\n",
        f"_Output dir: `{out_dir}`_\n\n",
        "Multi-source literature enrichment across 12 free academic-data APIs.\n",
        "## Tasks\n",
    ]
    for task, info in results.items():
        lines.append(f"- **{task}**: {info}")
    lines.append("\n## Sources used\n")
    lines.append("| # | Source | Purpose | Auth |")
    lines.append("|---|---|---|---|")
    sources = [
        ("OpenAlex", "Broad academic + abstracts + citations graph", "polite-pool email"),
        ("CrossRef", "DOI authority + metadata + references", "polite-pool email"),
        ("PubMed E-utilities", "Biomedical authority", "polite-pool email"),
        ("Europe PMC", "Full-text search for OA papers", "none"),
        ("Semantic Scholar", "Influential citations ranking", "none (rate limited)"),
        ("arXiv", "Preprint search", "none"),
        ("bioRxiv (via E-PMC)", "2024-2026 biomedical preprints", "none"),
        ("ICite (NIH)", "RCR citation impact", "none"),
        ("AFND", "HLA allele frequency (HTML scrape)", "none"),
        ("ClinicalTrials.gov v2", "Trial landscape", "none"),
        ("Unpaywall", "OA discovery", "polite-pool email"),
        ("CORE", "OA full-text (200M)", "free API key (skipped if absent)"),
    ]
    for i, (n, p, a) in enumerate(sources, 1):
        lines.append(f"| {i} | {n} | {p} | {a} |")

    lines.append("\n## How to re-run\n")
    lines.append("```bash")
    lines.append("project/.venv/bin/python project/notebooks_or_scripts/v17_lit_enrich_run.py")
    lines.append("```")
    lines.append("\nOutputs land back in this directory. Edit `claims`, `alleles`, "
                 "`trial_queries` in the runner to broaden scope.\n")

    p = out_dir / "07_summary.md"
    p.write_text("\n".join(lines), encoding="utf-8")
    return p
