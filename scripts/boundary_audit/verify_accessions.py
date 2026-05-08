#!/usr/bin/env python3
"""Verify accessions in master_dataset_catalog.csv against public APIs.

Hits NCBI eutils (GEO, BioProject), ENA browser API, and EBI biostudies
(ArrayExpress). Writes a verification report to
data_registry/reports/accession_verification.md and a per-row JSON cache so
re-runs are incremental.

Usage:
    python3 scripts/boundary_audit/verify_accessions.py
    python3 scripts/boundary_audit/verify_accessions.py --refresh   # ignore cache
    python3 scripts/boundary_audit/verify_accessions.py --limit 50  # cap rows
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Tuple

REPO = Path(__file__).resolve().parents[2]
MASTER = REPO / "data_registry" / "manifests" / "master_dataset_catalog.csv"
REPORT = REPO / "data_registry" / "reports" / "accession_verification.md"
CACHE = REPO / "data_registry" / "reports" / "_accession_verification_cache.json"

USER_AGENT = "THCA-data-registry/1.0 (research; contact: kukshomr@gmail.com)"
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENA_BASE = "https://www.ebi.ac.uk/ena/browser/api/xml"
BIOSTUDIES_BASE = "https://www.ebi.ac.uk/biostudies/api/v1/studies"
GDC_BASE = "https://api.gdc.cancer.gov/projects"
CBIO_BASE = "https://www.cbioportal.org/api/studies"
GWAS_BASE = "https://www.ebi.ac.uk/gwas/rest/api"
ICGC_BASE = "https://dcc.icgc.org/api/v1/projects"

# ~2.5 RPS NCBI without API key; ENA/biostudies are looser
NCBI_SLEEP = 0.4
OTHER_SLEEP = 0.15

GEO_RE = re.compile(r"^GSE\d+$")
ENA_RE = re.compile(r"^PRJEB\d+$")
SRA_RE = re.compile(r"^PRJNA\d+$")
AE_RE = re.compile(r"^E-MTAB-\d+$")
HRA_RE = re.compile(r"^HRA\d+$")
TCGA_RE = re.compile(r"^TCGA-[A-Z]{2,5}$")
ICGC_PROJECT_RE = re.compile(r"^(THCA|BRCA|LIRI|PACA|RECA|MELA|EOPC|GACA)-[A-Z]{2,3}$")
NCT_RE = re.compile(r"^NCT\d+$")

# Inherently not-API-verifiable sources (no per-entry endpoint)
MANUAL_SOURCES = {
    "AFND", "internal", "BMDW", "NMDP", "DepMap", "DICE", "Tools",
    "DSigDB", "MSigDB", "KEGG", "Reactome", "Tabula_Sapiens",
    "HPA", "DiseaseMeth", "MethBank", "Wanderer", "PanglaoDB",
    "CellTypist", "BLUEPRINT", "Roadmap_Epigenomics", "ENCODE",
    "TIMER2", "xCell", "CIBERSORTx", "Open_Targets", "DisGeNET",
    "HPO", "ClinGen", "SysteMHC", "TANTIGEN", "TSNAdb", "NEPdb",
    "dbPepNeo", "TESLA", "VDJdb", "McPAS-TCR", "IEDB", "MHCflurry",
    "NetMHCpan", "NetMHCIIpan", "HLA_dictionary", "HLA_Ligand_Atlas",
    "HLA-net", "HLA_Ligand_Atlas", "PRISM_Repurposing", "PRISM",
    "HuBMAP", "4DN", "Human_Cell_Atlas", "CellxGene", "Broad_SCP",
    "Pisco_Lab", "10x_Genomics", "Stereo_seq", "TICA", "PanImmune",
    "PanCanAtlas", "MET500", "HCMI", "Project_DRIVE", "GDSC",
    "GDSC2", "CCLE", "CCLE_Proteomics", "DepMap_Proteomics",
    "GDC_PDC", "ChEMBL", "DrugBank", "Open_PHACTS", "IMPC",
    "TumorFusionGeneDB", "FusionHub", "ChimerDB", "PCAWG",
    "ICGC_PCAWG", "Hartwig_Medical_Foundation", "TCIA",
    "Connectivity_Map", "CPTAC", "GTEx", "GTEx_Proteogenomic",
    "ProteomicsDB", "ARCHS4", "recount3", "PanNuke",
    "FinnGen", "UK_Biobank", "BioBank_Japan", "Pan_UKBB",
    "BBJ_PheWeb", "IEU_OpenGWAS", "Synapse", "ORIEN_TCC",
    "AACR_GENIE", "AACR_GENIE_BPC", "Korean_Network_PCT",
    "Korean_Tissue_Biobank", "Korean_Pathology_AI", "KOBIC",
    "KISTI", "KMDP", "KCDC", "KSHI", "APBMDR", "APHIA", "IHW",
    "KoreanGenomeProject", "KoGES", "KoGES_KARE", "KoGES_KOEX",
    "KoGES_HEXA", "KoGES_KBN", "Fukushima", "Chernobyl_Tissue_Bank",
    "Korean_Cancer_Genome_Atlas", "Korean_Hereditary_Cancer_Reg",
    "Korean_NM_Society", "Korean_Thyroid_Assoc", "KCGA",
    "BBMRI", "All_of_Us", "AMR", "Genes.ai", "TOPMed", "KOREA_Project",
    "AOTRC", "French_Cancer_Network", "Global_Thyroid_Network",
    "PheWAS_Catalog", "Mayo_Clinic", "ClinicalTrials_gov", "KNHANES",
    "ClinVar", "COSMIC", "COSMIC_Fusion", "dbSNP", "gnomAD",
    "Korean_donors_research",
    "ICGC_ARGO",
}

THYROID_TERMS = (
    "thyroid", "thyroïd", "thca", "ptc", "papillary",
    "follicular", "anaplastic", "atc", "pdtc", "medullary",
    "graves", "hashimoto", "tirads", "goiter", "thyroiditis",
)


def http_get(url: str, timeout: int = 15) -> Optional[str]:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json, application/xml, text/xml"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.read().decode("utf-8", errors="ignore")
    except urllib.error.HTTPError as e:
        return f"__HTTPERR__{e.code}"
    except Exception as e:
        return f"__ERR__{type(e).__name__}: {e}"


@dataclass
class Result:
    accession: str
    kind: str
    exists: Optional[bool] = None  # True/False/None(error)
    title: str = ""
    n_samples: str = ""
    platform: str = ""
    thyroid_match: bool = False
    note: str = ""


def is_thyroid(title: str) -> bool:
    t = title.lower()
    return any(term in t for term in THYROID_TERMS)


def verify_geo(acc: str) -> Result:
    r = Result(accession=acc, kind="GEO")
    q = urllib.parse.quote(f"{acc}[Accession]")
    body = http_get(f"{NCBI_BASE}/esearch.fcgi?db=gds&term={q}&retmode=json")
    if not body or body.startswith("__"):
        r.exists = None
        r.note = body or "empty response"
        return r
    try:
        data = json.loads(body)
        ids = data.get("esearchresult", {}).get("idlist", [])
        if not ids:
            r.exists = False
            r.note = "no UID in esearch"
            return r
        # Find the canonical Series record (UID starts with 200 = GSE)
        series_id = next((x for x in ids if x.startswith("200")), ids[0])
        time.sleep(NCBI_SLEEP)
        body2 = http_get(f"{NCBI_BASE}/esummary.fcgi?db=gds&id={series_id}&retmode=json")
        if not body2 or body2.startswith("__"):
            r.exists = True
            r.note = f"esummary failed: {body2}"
            return r
        data2 = json.loads(body2)
        rec = data2.get("result", {}).get(series_id, {})
        r.exists = True
        r.title = (rec.get("title") or "")[:160]
        r.n_samples = str(rec.get("n_samples", ""))
        r.platform = "/".join((rec.get("gpl") or "").split(";")[:3])
        r.thyroid_match = is_thyroid(r.title) or is_thyroid(rec.get("summary", ""))
    except Exception as e:
        r.exists = None
        r.note = f"parse error: {e}"
    return r


def verify_bioproject(acc: str) -> Result:
    r = Result(accession=acc, kind="SRA/BioProject")
    q = urllib.parse.quote(acc)
    body = http_get(f"{NCBI_BASE}/esearch.fcgi?db=bioproject&term={q}&retmode=json")
    if not body or body.startswith("__"):
        r.exists = None
        r.note = body or "empty"
        return r
    try:
        data = json.loads(body)
        ids = data.get("esearchresult", {}).get("idlist", [])
        if not ids:
            r.exists = False
            return r
        time.sleep(NCBI_SLEEP)
        body2 = http_get(f"{NCBI_BASE}/esummary.fcgi?db=bioproject&id={ids[0]}&retmode=json")
        if not body2 or body2.startswith("__"):
            r.exists = True
            return r
        rec = json.loads(body2).get("result", {}).get(ids[0], {})
        r.exists = True
        r.title = (rec.get("project_title") or rec.get("project_name") or "")[:160]
        r.thyroid_match = is_thyroid(r.title) or is_thyroid(rec.get("project_description", ""))
    except Exception as e:
        r.exists = None
        r.note = f"parse error: {e}"
    return r


def verify_ena(acc: str) -> Result:
    r = Result(accession=acc, kind="ENA")
    body = http_get(f"{ENA_BASE}/{acc}")
    if not body or body.startswith("__"):
        r.exists = None
        r.note = body or "empty"
        return r
    if "<PROJECT" in body or "<STUDY" in body:
        r.exists = True
        m = re.search(r"<TITLE>(.*?)</TITLE>", body, re.DOTALL)
        if m:
            r.title = m.group(1).strip()[:160]
        r.thyroid_match = is_thyroid(r.title) or is_thyroid(body[:2000])
    else:
        r.exists = False
        r.note = "no PROJECT/STUDY tag"
    return r


def verify_arrayexpress(acc: str) -> Result:
    r = Result(accession=acc, kind="ArrayExpress")
    body = http_get(f"{BIOSTUDIES_BASE}/{acc}")
    if not body or body.startswith("__"):
        if body and "404" in body:
            r.exists = False
            r.note = "404"
        else:
            r.exists = None
            r.note = body or "empty"
        return r
    try:
        data = json.loads(body)
        if data.get("accno") == acc:
            r.exists = True
            for attr in data.get("attributes", []) or []:
                if attr.get("name", "").lower() == "title":
                    r.title = (attr.get("value") or "")[:160]
                    break
            r.thyroid_match = is_thyroid(r.title) or is_thyroid(body[:2000])
        else:
            r.exists = False
    except Exception as e:
        r.exists = None
        r.note = f"parse error: {e}"
    return r


def verify_gdc_project(acc: str) -> Result:
    r = Result(accession=acc, kind="GDC")
    body = http_get(f"{GDC_BASE}/{acc}")
    if not body or body.startswith("__"):
        r.exists = None
        r.note = body or "empty"
        return r
    try:
        data = json.loads(body)
        d = data.get("data") or {}
        if d.get("project_id") == acc:
            r.exists = True
            r.title = (d.get("name") or "")[:160]
            sites = d.get("primary_site") or []
            if isinstance(sites, list):
                site_str = ", ".join(sites)[:80]
            else:
                site_str = str(sites)[:80]
            r.note = f"primary_site={site_str}"
            r.thyroid_match = is_thyroid(r.title) or is_thyroid(site_str)
        else:
            r.exists = False
    except Exception as e:
        r.exists = None
        r.note = f"parse error: {e}"
    return r


def verify_cbio_study(acc: str) -> Result:
    r = Result(accession=acc, kind="cBioPortal")
    body = http_get(f"{CBIO_BASE}/{acc}")
    if not body or body.startswith("__"):
        if body and "404" in body:
            r.exists = False
            r.note = "404"
        else:
            r.exists = None
            r.note = body or "empty"
        return r
    try:
        data = json.loads(body)
        if data.get("studyId") == acc:
            r.exists = True
            r.title = (data.get("name") or "")[:160]
            r.thyroid_match = is_thyroid(r.title) or is_thyroid(data.get("description", ""))
        else:
            r.exists = False
    except Exception as e:
        r.exists = None
        r.note = f"parse error: {e}"
    return r


def verify_icgc_project(acc: str) -> Result:
    r = Result(accession=acc, kind="ICGC")
    body = http_get(f"{ICGC_BASE}/{acc}")
    if not body or body.startswith("__"):
        if body and "404" in body:
            r.exists = False
            r.note = "404"
        else:
            r.exists = None
            r.note = body or "empty"
        return r
    try:
        data = json.loads(body)
        if data.get("id") == acc or data.get("projectId") == acc:
            r.exists = True
            r.title = (data.get("name") or data.get("primarySite") or "")[:160]
            r.thyroid_match = is_thyroid(r.title) or is_thyroid(str(data.get("primarySite", "")))
        else:
            r.exists = False
    except Exception as e:
        r.exists = None
        r.note = f"parse error: {e}"
    return r


def verify_pubmed_publication(title: str, year: str = "") -> Result:
    """Fuzzy-search PubMed for a publication title; verify the row points at a real paper."""
    r = Result(accession=title[:60], kind="PubMed")
    if not title:
        r.exists = False
        r.note = "no title"
        return r
    # PubMed esearch by title (use only first 80 chars to avoid query bloat)
    q = urllib.parse.quote(f"{title[:80]}[Title]" + (f" AND {year}[PDAT]" if year else ""))
    body = http_get(f"{NCBI_BASE}/esearch.fcgi?db=pubmed&term={q}&retmode=json&retmax=3")
    if not body or body.startswith("__"):
        r.exists = None
        r.note = body or "empty"
        return r
    try:
        ids = json.loads(body).get("esearchresult", {}).get("idlist", [])
        if not ids:
            r.exists = False
            r.note = "no PMID for title"
            return r
        time.sleep(NCBI_SLEEP)
        body2 = http_get(f"{NCBI_BASE}/esummary.fcgi?db=pubmed&id={ids[0]}&retmode=json")
        if body2 and not body2.startswith("__"):
            rec = json.loads(body2).get("result", {}).get(ids[0], {})
            r.title = (rec.get("title") or "")[:160]
        r.exists = True
        r.thyroid_match = is_thyroid(title) or is_thyroid(r.title)
    except Exception as e:
        r.exists = None
        r.note = f"parse error: {e}"
    return r


def dispatch(acc: str, *, source: str = "", title: str = "", year: str = "") -> Optional[Result]:
    if GEO_RE.match(acc):
        return verify_geo(acc)
    if ENA_RE.match(acc):
        return verify_ena(acc)
    if SRA_RE.match(acc):
        return verify_bioproject(acc)
    if AE_RE.match(acc):
        return verify_arrayexpress(acc)
    if TCGA_RE.match(acc):
        return verify_gdc_project(acc)
    if ICGC_PROJECT_RE.match(acc):
        return verify_icgc_project(acc)
    if source == "cBioPortal":
        return verify_cbio_study(acc)
    if source == "GDC" and acc.lower().startswith(("thca_", "thyroid_")):
        # cBioPortal-style study IDs that the catalog labeled as GDC
        return verify_cbio_study(acc)
    if source == "Publication" and title:
        return verify_pubmed_publication(title, year)
    return None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--refresh", action="store_true", help="Ignore on-disk cache")
    p.add_argument("--limit", type=int, default=None, help="Stop after N verifiable rows")
    p.add_argument("--only", choices=["GEO", "ENA", "SRA", "AE"], help="Restrict to one source")
    args = p.parse_args()

    cache = {}
    if CACHE.exists() and not args.refresh:
        try:
            cache = json.loads(CACHE.read_text())
        except Exception:
            cache = {}

    rows = list(csv.DictReader(MASTER.open()))
    verifiable = []
    for row in rows:
        acc = (row.get("accession_or_id") or "").strip()
        if not acc:
            continue
        # Skip DUPLICATE markers
        if "DUPLICATE" in row.get("title", ""):
            continue
        source = (row.get("source_database") or "").strip()
        title = (row.get("title") or "").strip()
        year = (row.get("year") or "").strip()
        # Determine if this row has a verifier path
        has_pattern = (
            GEO_RE.match(acc)
            or ENA_RE.match(acc)
            or SRA_RE.match(acc)
            or AE_RE.match(acc)
            or TCGA_RE.match(acc)
            or ICGC_PROJECT_RE.match(acc)
            or source in ("cBioPortal",)
            or (source == "Publication" and title)
        )
        if not has_pattern:
            continue
        if args.only:
            kind_map = {"GEO": GEO_RE, "ENA": ENA_RE, "SRA": SRA_RE, "AE": AE_RE}
            if args.only in kind_map:
                if not kind_map[args.only].match(acc):
                    continue
            elif args.only == "GDC" and not TCGA_RE.match(acc):
                continue
            elif args.only == "ICGC" and not ICGC_PROJECT_RE.match(acc):
                continue
            elif args.only == "PUB" and source != "Publication":
                continue
            elif args.only == "CBIO" and source != "cBioPortal":
                continue
        verifiable.append((row["dataset_id"], row["zone"], acc, source, title, year))

    if args.limit:
        verifiable = verifiable[: args.limit]

    print(f"verifying {len(verifiable)} accessions ...", file=sys.stderr)
    results = []
    n_done = 0
    for ds_id, zone, acc, source, title, year in verifiable:
        # Cache key includes source/title for publication entries (acc alone isn't unique)
        if source == "Publication":
            cache_key = f"PUB::{acc}::{title[:60]}"
        else:
            cache_key = f"{acc}"
        if cache_key in cache and not args.refresh:
            r = Result(**cache[cache_key])
        else:
            r = dispatch(acc, source=source, title=title, year=year)
            if r is None:
                continue
            cache[cache_key] = r.__dict__
            time.sleep(NCBI_SLEEP if r.kind in ("GEO", "SRA/BioProject", "PubMed") else OTHER_SLEEP)
        results.append((ds_id, zone, r))
        n_done += 1
        if n_done % 25 == 0:
            print(f"  {n_done}/{len(verifiable)}", file=sys.stderr)
            CACHE.write_text(json.dumps(cache, indent=2))

    CACHE.write_text(json.dumps(cache, indent=2))

    # Build report
    by_status = {"exists": [], "missing": [], "error": []}
    by_subj = {"thyroid": [], "non_thyroid_warning": []}
    for ds_id, zone, r in results:
        if r.exists is True:
            by_status["exists"].append((ds_id, zone, r))
            if r.thyroid_match:
                by_subj["thyroid"].append((ds_id, zone, r))
            else:
                by_subj["non_thyroid_warning"].append((ds_id, zone, r))
        elif r.exists is False:
            by_status["missing"].append((ds_id, zone, r))
        else:
            by_status["error"].append((ds_id, zone, r))

    lines = ["# Accession verification report", ""]
    lines.append(f"_Run at {time.strftime('%Y-%m-%d %H:%M:%S')} UTC_")
    lines.append("")
    lines.append("| Status | Count |")
    lines.append("|---|---|")
    lines.append(f"| Verified exists | **{len(by_status['exists'])}** |")
    lines.append(f"|   - thyroid-match in title/summary | {len(by_subj['thyroid'])} |")
    lines.append(f"|   - **non-thyroid warning** (review!) | {len(by_subj['non_thyroid_warning'])} |")
    lines.append(f"| Missing (no UID) | {len(by_status['missing'])} |")
    lines.append(f"| Verification error (network/parse) | {len(by_status['error'])} |")
    lines.append(f"| **Total checked** | {len(results)} |")
    lines.append("")

    if by_subj["non_thyroid_warning"]:
        lines.append("## ⚠ Non-thyroid warnings (verified, but title/summary doesn't match thyroid)")
        lines.append("")
        lines.append("These accessions exist but their public titles do not contain thyroid-related terms. They may have been mis-listed in master_dataset_catalog.csv as thyroid datasets when they are not. Manual review required.")
        lines.append("")
        for ds_id, zone, r in sorted(by_subj["non_thyroid_warning"]):
            lines.append(f"- `{ds_id}` [{zone}] **{r.accession}** ({r.kind}) — title: _{r.title or '?'}_")
        lines.append("")

    if by_status["missing"]:
        lines.append("## Missing accessions (API returned no UID)")
        lines.append("")
        lines.append("These accessions could not be found via the public API. They may be wrong, withdrawn, or in restricted-access archives.")
        lines.append("")
        for ds_id, zone, r in sorted(by_status["missing"]):
            lines.append(f"- `{ds_id}` [{zone}] **{r.accession}** ({r.kind}) — {r.note}")
        lines.append("")

    if by_status["error"]:
        lines.append("## Verification errors (network / parse failures)")
        lines.append("")
        for ds_id, zone, r in sorted(by_status["error"]):
            lines.append(f"- `{ds_id}` [{zone}] **{r.accession}** ({r.kind}) — {r.note}")
        lines.append("")

    lines.append("## Verified thyroid datasets")
    lines.append("")
    lines.append(f"({len(by_subj['thyroid'])} entries — full list)")
    lines.append("")
    for ds_id, zone, r in sorted(by_subj["thyroid"]):
        ns = f"n={r.n_samples}" if r.n_samples else ""
        plat = f"GPL{r.platform}" if r.platform else ""
        meta = " | ".join(x for x in [ns, plat] if x)
        suffix = f" — {meta}" if meta else ""
        lines.append(f"- `{ds_id}` [{zone}] **{r.accession}** ({r.kind}){suffix} — _{r.title or '?'}_")
    lines.append("")

    REPORT.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {REPORT.relative_to(REPO)} ({len(results)} entries; {len(by_subj['non_thyroid_warning'])} non-thyroid warnings)")


if __name__ == "__main__":
    main()
