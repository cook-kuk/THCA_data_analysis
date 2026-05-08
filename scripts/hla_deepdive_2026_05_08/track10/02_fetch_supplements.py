#!/usr/bin/env python3
"""
Track 10 — Step 2: Fetch full-text / supplements for the 5 resolved papers.

Strategy: try unpaywall (best_oa_location); fall back to PMC; fall back to
publisher landing page. Document each result. NEVER use Sci-Hub.
"""
from __future__ import annotations
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path
import pandas as pd

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track10_korean_lit")
PDFS = OUT / "source_pdfs"
PDFS.mkdir(parents=True, exist_ok=True)
TBL = OUT / "tables"

# Resolved citations (manually confirmed from candidate JSON)
RESOLVED = {
    "Shin_2019": dict(
        pmid="31091281",
        doi="10.1371/journal.pone.0216941",
        title="HLA alleles, especially amino-acid signatures of HLA-DPB1, might contribute to the molecular pathogenesis of early-onset autoimmune thyroid disease",
        journal="PLoS ONE",
        year="2019",
        volume="14",
        pages="e0216941",
        authors="Shin DH, Baek IC, Kim HJ, Choi EJ, Ahn M, Jung MH, Suh BK, Cho WK, Kim TG",
        n_case=71, n_control=142,
        cohort_label="Korean pediatric AITD (Graves' + Hashimoto)",
        is_open=True,
    ),
    "Cho_2011": dict(
        pmid="21952423",
        doi="10.1159/000331134",
        title="Association of HLA alleles with autoimmune thyroid disease in Korean children",
        journal="Hormone Research in Paediatrics",
        year="2011",
        volume="76",
        pages="328-334",
        authors="Cho WK, Jung MH, Choi EJ, Choi HB, Kim TG, Suh BK",
        n_case=None, n_control=None,
        cohort_label="Korean pediatric AITD",
        is_open=False,
    ),
    "Park_2005": dict(
        pmid="15993720",
        doi="10.1016/j.humimm.2005.03.001",
        title="Association of HLA-DR and -DQ genes with Graves disease in Koreans",
        journal="Human Immunology",
        year="2005",
        volume="66",
        pages="741-747",
        authors="Park MH, Park YJ, Song EY, Park H, Kim TY, Park DJ, Park KS, Cho BY",
        n_case=88, n_control=104,
        cohort_label="Korean adult Graves' disease",
        is_open=False,
    ),
    "Jang_2011": dict(
        pmid="21062236",
        doi="10.3109/08820139.2010.525571",
        title="Identification of HLA-DRB1 alleles associated with Graves' disease in Koreans by sequence-based typing",
        journal="Immunological Investigations",
        year="2011",
        volume="40",
        pages="172-182",
        authors="Jang HW, Shin HW, Cho HJ, Kim HK, Lee JI, Kim SW, Kim JW, Chung JH",
        n_case=None, n_control=None,
        cohort_label="Korean adult Graves' disease (DRB1 SBT)",
        is_open=False,
    ),
    "Baek_2021": dict(
        pmid="",  # not on PubMed
        doi="10.1111/tan.14134",
        title="Distributions of HLA-A, -B, and -DRB1 alleles typed by amplicon-based next-generation sequencing in Korean preschool children",
        journal="HLA",
        year="2021",
        volume="97",
        pages="six-digit deep typing reference cohort",
        authors="Baek IC, Choi EJ, Shin DH, Kim HJ, Choi H, Kim TG",
        n_case=None, n_control=None,
        cohort_label="Korean class-II NGS reference population (no AITD case-control; population baseline)",
        is_open=False,
    ),
}


def http_json(url: str) -> dict | None:
    try:
        req = urllib.request.Request(
            url, headers={"User-Agent": "track10-supplement-fetcher/1.0 (mailto:kukshomr@gmail.com)"})
        with urllib.request.urlopen(req, timeout=30) as r:
            return json.loads(r.read().decode("utf-8"))
    except Exception as e:
        print(f"   ! json fail: {e}")
        return None


def http_get(url: str, dest: Path, max_mb: int = 25) -> bool:
    try:
        req = urllib.request.Request(
            url, headers={
                "User-Agent": "Mozilla/5.0 (compatible; track10-fetcher/1.0; +mailto:kukshomr@gmail.com)",
                "Accept": "*/*",
            })
        with urllib.request.urlopen(req, timeout=60) as r:
            data = r.read(max_mb * 1024 * 1024)
            if not data:
                return False
            dest.write_bytes(data)
            return True
    except Exception as e:
        print(f"   ! get fail: {url[:80]} -> {e}")
        return False


def try_unpaywall(doi: str, dest_dir: Path, name: str) -> dict:
    if not doi:
        return dict(method="unpaywall", success=False, reason="no_doi")
    url = f"https://api.unpaywall.org/v2/{doi}?email=kukshomr@gmail.com"
    data = http_json(url)
    if not data:
        return dict(method="unpaywall", success=False, reason="api_fail")
    is_oa = data.get("is_oa", False)
    if not is_oa:
        return dict(method="unpaywall", success=False, reason="not_oa", license=data.get("oa_status", ""))
    loc = data.get("best_oa_location") or {}
    pdf_url = loc.get("url_for_pdf") or loc.get("url")
    landing = loc.get("url_for_landing_page", "")
    if not pdf_url:
        return dict(method="unpaywall", success=False, reason="no_pdf_url", landing=landing)
    out_pdf = dest_dir / f"{name}.pdf"
    ok = http_get(pdf_url, out_pdf)
    return dict(
        method="unpaywall",
        success=ok and out_pdf.stat().st_size > 5_000,
        path=str(out_pdf) if ok else "",
        bytes=out_pdf.stat().st_size if ok else 0,
        pdf_url=pdf_url,
        landing=landing,
        license=data.get("oa_status", ""),
    )


def try_pmc(pmid: str, dest_dir: Path, name: str) -> dict:
    if not pmid:
        return dict(method="pmc", success=False, reason="no_pmid")
    # Convert PMID -> PMCID
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=pmc&retmode=json&id={pmid}"
    data = http_json(url)
    if not data:
        return dict(method="pmc", success=False, reason="elink_fail")
    pmcids = []
    for ls in data.get("linksets", []):
        for ldb in ls.get("linksetdbs", []):
            if ldb.get("dbto") == "pmc":
                pmcids.extend(ldb.get("links", []))
    if not pmcids:
        return dict(method="pmc", success=False, reason="not_in_pmc")
    pmcid = pmcids[0]
    # PMC OA service: https://pmc.ncbi.nlm.nih.gov/articles/PMC<id>/pdf/...
    pdf_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmcid}/pdf/"
    out_pdf = dest_dir / f"{name}_PMC{pmcid}.pdf"
    ok = http_get(pdf_url, out_pdf)
    return dict(
        method="pmc",
        success=ok and out_pdf.stat().st_size > 5_000,
        path=str(out_pdf) if ok else "",
        bytes=out_pdf.stat().st_size if ok else 0,
        pmcid=pmcid,
    )


def try_plos_supplements(doi: str, dest_dir: Path, name: str) -> dict:
    """PLoS ONE supplements: https://journals.plos.org/plosone/article/file?id=...&type=supplementary"""
    if "10.1371" not in doi:
        return dict(method="plos_supp", success=False, reason="not_plos")
    # Article XML page lists supplementary info; we need to scrape it.
    article_url = f"https://journals.plos.org/plosone/article?id={doi}"
    try:
        req = urllib.request.Request(article_url, headers={"User-Agent": "Mozilla/5.0 track10"})
        with urllib.request.urlopen(req, timeout=30) as r:
            html = r.read().decode("utf-8", errors="ignore")
    except Exception as e:
        return dict(method="plos_supp", success=False, reason=f"page_fail:{e}")
    # supplemental files are referenced as ?id=info:doi/<doi>.s001&type=supplementary
    import re
    matches = re.findall(r"info:doi/[^&\"\']+\.s\d+", html)
    matches = list(dict.fromkeys(matches))
    files = []
    for m in matches:
        url = f"https://journals.plos.org/plosone/article/file?id={m}&type=supplementary"
        # Determine extension by HEAD request
        suffix = m.split(".")[-1]
        out = dest_dir / f"{name}_{suffix}.bin"
        if http_get(url, out, max_mb=15):
            # detect extension by magic
            head = out.read_bytes()[:4]
            new_suffix = ".pdf" if head.startswith(b"%PDF") else (
                ".docx" if head.startswith(b"PK") else ".bin")
            new = out.with_suffix(new_suffix)
            out.rename(new)
            files.append(dict(name=str(new.name), bytes=new.stat().st_size, src_id=m))
    return dict(
        method="plos_supp",
        success=len(files) > 0,
        files=files,
        article_url=article_url,
    )


# ----- run -----
results = {}
for name, meta in RESOLVED.items():
    print(f"\n[{name}] PMID={meta['pmid']}  DOI={meta['doi']}  ({meta['journal']} {meta['year']})")
    rec = dict(meta=meta, attempts=[])
    # 1) PMC if PMID
    if meta["pmid"]:
        r = try_pmc(meta["pmid"], PDFS, name)
        rec["attempts"].append(r)
        print(f"   PMC: success={r['success']} reason={r.get('reason','')}")
    # 2) Unpaywall
    r = try_unpaywall(meta["doi"], PDFS, name)
    rec["attempts"].append(r)
    print(f"   Unpaywall: success={r['success']} reason={r.get('reason','')} license={r.get('license','')}")
    # 3) PLoS supplements (Shin 2019)
    if "10.1371" in meta["doi"]:
        r = try_plos_supplements(meta["doi"], PDFS, name)
        rec["attempts"].append(r)
        print(f"   PLoS supp: success={r['success']} files={len(r.get('files',[]))}")
    results[name] = rec
    time.sleep(0.5)

with open(OUT / "supplement_fetch_log.json", "w") as f:
    json.dump(results, f, indent=2)
print(f"\n[saved] {OUT/'supplement_fetch_log.json'}")
print("\nDownloaded files:")
for f in sorted(PDFS.iterdir()):
    print(f"   {f.name}  {f.stat().st_size:,} B")

# Build resolution table with availability flag
rows = []
for name, rec in results.items():
    m = rec["meta"]
    fetched_main = ""
    fetched_supp = []
    for a in rec["attempts"]:
        if a.get("success"):
            if a["method"] in ("pmc", "unpaywall"):
                fetched_main = a.get("path", "")
            elif a["method"] == "plos_supp":
                fetched_supp = a.get("files", [])
    rows.append(dict(
        paper_id=name,
        pmid=m["pmid"], doi=m["doi"], title=m["title"],
        journal=m["journal"], year=m["year"], volume=m["volume"], pages=m["pages"],
        authors=m["authors"],
        n_case=m["n_case"] or "", n_control=m["n_control"] or "",
        cohort_label=m["cohort_label"],
        full_text_path=fetched_main,
        n_supplements=len(fetched_supp),
        supp_files=";".join(s["name"] for s in fetched_supp),
        url_doi=f"https://doi.org/{m['doi']}",
        url_pubmed=f"https://pubmed.ncbi.nlm.nih.gov/{m['pmid']}/" if m["pmid"] else "",
    ))
df = pd.DataFrame(rows)
df.to_csv(TBL / "T02_resolved_citations.tsv", sep="\t", index=False)
print(f"[saved] {TBL/'T02_resolved_citations.tsv'}")
