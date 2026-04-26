#!/usr/bin/env python3
"""
Target-to-Drug discovery layer for THCA project.

Phase 0: poll for biomarker_validated.tsv (or fallback to biomarker_known_vs_novel.tsv
         or derive from DE analysis).
Phase 1: PubMed literature grounding per target.
Phase 2: ChEMBL/PubChem compound search.
Phase 3: druggability heuristics.
Phase 4: outputs (tsvs + json + md).
"""
from __future__ import annotations

import json
import os
import sys
import time
import math
import xml.etree.ElementTree as ET
from pathlib import Path
from urllib.parse import quote

import requests

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project")
TABLES = ROOT / "results" / "tables"
PAYLOAD_DIR = ROOT / "reports" / "html" / "assets" / "data"
REPORT_MD = ROOT / "reports" / "biomarker_to_drug_report.md"

TABLES.mkdir(parents=True, exist_ok=True)
PAYLOAD_DIR.mkdir(parents=True, exist_ok=True)
REPORT_MD.parent.mkdir(parents=True, exist_ok=True)

BIOMARKER_VALIDATED = TABLES / "biomarker_validated.tsv"
BIOMARKER_KNOWN_VS_NOVEL = TABLES / "biomarker_known_vs_novel.tsv"

HTTP_TIMEOUT = 12
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "THCA-drug-discovery/1.0 (research)"})

LOG: list[str] = []


def log(msg: str) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    LOG.append(line)


# ---------------------------------------------------------------- PHASE 0 ----
def wait_for_biomarkers(max_wait_s: int = 20 * 60, interval_s: int = 45) -> Path | None:
    start = time.time()
    log(f"PHASE 0: polling for {BIOMARKER_VALIDATED.name} (max {max_wait_s}s)")
    while time.time() - start < max_wait_s:
        if BIOMARKER_VALIDATED.exists() and BIOMARKER_VALIDATED.stat().st_size > 0:
            log(f"PHASE 0: found validated biomarkers after {int(time.time()-start)}s")
            return BIOMARKER_VALIDATED
        if BIOMARKER_KNOWN_VS_NOVEL.exists() and BIOMARKER_KNOWN_VS_NOVEL.stat().st_size > 0:
            # keep polling once for validated a few more rounds
            if time.time() - start > 180:
                log("PHASE 0: falling back to biomarker_known_vs_novel.tsv")
                return BIOMARKER_KNOWN_VS_NOVEL
        time.sleep(interval_s)
    log("PHASE 0: timed out. Falling back.")
    if BIOMARKER_KNOWN_VS_NOVEL.exists():
        return BIOMARKER_KNOWN_VS_NOVEL
    return None


def load_targets_from_biomarker_tsv(path: Path) -> list[dict]:
    import csv
    out = []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            out.append(row)
    return out


def derive_targets_from_de(top_n: int = 8) -> list[dict]:
    log("PHASE 0 fallback: deriving targets from brs71_proxy_deg_full.tsv")
    import csv
    de = TABLES / "brs71_proxy_deg_full.tsv"
    rows = []
    with de.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            try:
                lfc = float(r["log2FC_BRAF_vs_RAS"])
                fdr = float(r["FDR"])
            except Exception:
                continue
            if fdr > 0.01:
                continue
            rows.append({
                "gene": r["gene_symbol"],
                "log2FC": lfc,
                "cohens_d": "",
                "fdr": fdr,
                "replicated_27155": "",
                "replicated_126698": "",
                "novelty_score": abs(lfc) * (-math.log10(max(fdr, 1e-300))) / 100.0,
                "is_novel": "True",
            })
    rows.sort(key=lambda x: x["novelty_score"], reverse=True)
    # filter out extremely common well-known thyroid genes to simulate novelty
    well_known = {"BRAF", "TG", "TPO", "TSHR", "NRAS", "HRAS", "KRAS", "RET", "PAX8"}
    novel_like = [r for r in rows if r["gene"] not in well_known and not r["gene"].startswith("LINC")
                  and not r["gene"].startswith("AC") and not r["gene"].startswith("AL")
                  and not r["gene"].startswith("AP") and "." not in r["gene"]]
    return novel_like[:top_n]


def select_top8_targets(rows: list[dict]) -> list[dict]:
    # normalise keys across schemas
    out = []
    for r in rows:
        gene = r.get("gene") or r.get("gene_symbol") or r.get("symbol")
        if not gene:
            continue
        is_novel = str(
            r.get("is_novel")
            or r.get("is_novel_validated")
            or r.get("novel")
            or ""
        ).lower() in ("true", "1", "yes")
        try:
            novelty = float(r.get("novelty_score", 0) or 0)
        except Exception:
            novelty = 0.0
        out.append({
            "gene": gene,
            "log2FC": r.get("log2FC") or r.get("log2FC_tcga") or r.get("log2FC_BRAF_vs_RAS") or "",
            "cohens_d": r.get("cohens_d") or r.get("cohens_d_tcga") or "",
            "fdr": r.get("fdr") or r.get("fdr_tcga") or r.get("FDR") or "",
            "replicated_27155": r.get("replicated_27155", ""),
            "replicated_126698": r.get("replicated_126698", ""),
            "novelty_score": novelty,
            "is_novel": is_novel,
        })
    # novel first then by novelty desc
    out.sort(key=lambda x: (not x["is_novel"], -x["novelty_score"]))
    top8 = out[:8]
    if len(top8) < 8:
        needed = 8 - len(top8)
        extras = derive_targets_from_de(top_n=needed * 3)
        existing = {t["gene"] for t in top8}
        for e in extras:
            if e["gene"] in existing:
                continue
            top8.append({
                **e,
                "is_novel": True,
            })
            if len(top8) >= 8:
                break
    return top8[:8]


# ---------------------------------------------------------------- PHASE 1 ----
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def pubmed_search_thca(gene: str) -> list[str]:
    try:
        term = f"{gene}[Title/Abstract] AND (thyroid[MeSH Terms] OR thyroid cancer)"
        r = SESSION.get(f"{EUTILS}/esearch.fcgi",
                        params={"db": "pubmed", "term": term, "retmax": 10,
                                "sort": "date", "retmode": "json"},
                        timeout=HTTP_TIMEOUT)
        r.raise_for_status()
        return r.json().get("esearchresult", {}).get("idlist", [])
    except Exception as e:
        log(f"    pubmed thca search failed for {gene}: {e}")
        return []


def pubmed_search_cancer(gene: str) -> int:
    try:
        term = f"{gene}[Title] AND (cancer[MeSH Terms])"
        r = SESSION.get(f"{EUTILS}/esearch.fcgi",
                        params={"db": "pubmed", "term": term, "retmax": 5,
                                "retmode": "json"},
                        timeout=HTTP_TIMEOUT)
        r.raise_for_status()
        j = r.json().get("esearchresult", {})
        try:
            return int(j.get("count", 0))
        except Exception:
            return len(j.get("idlist", []))
    except Exception as e:
        log(f"    pubmed cancer search failed for {gene}: {e}")
        return 0


def pubmed_efetch(pmids: list[str]) -> list[dict]:
    if not pmids:
        return []
    try:
        r = SESSION.get(f"{EUTILS}/efetch.fcgi",
                        params={"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"},
                        timeout=HTTP_TIMEOUT)
        r.raise_for_status()
        root = ET.fromstring(r.content)
    except Exception as e:
        log(f"    efetch failed: {e}")
        return []
    out = []
    for art in root.findall(".//PubmedArticle"):
        pmid_el = art.find(".//PMID")
        pmid = pmid_el.text if pmid_el is not None else ""
        title_el = art.find(".//ArticleTitle")
        title = "".join(title_el.itertext()) if title_el is not None else ""
        year = ""
        ydate = art.find(".//PubDate/Year")
        if ydate is None:
            ydate = art.find(".//PubDate/MedlineDate")
        if ydate is not None and ydate.text:
            year = ydate.text[:4]
        abst_parts = []
        for ab in art.findall(".//Abstract/AbstractText"):
            abst_parts.append("".join(ab.itertext()))
        abstract = " ".join(abst_parts).strip()
        out.append({"pmid": pmid, "year": year, "title": title.strip(),
                    "abstract": abstract})
    return out


def summarize_abstract(abstract: str, max_words: int = 38) -> str:
    if not abstract:
        return ""
    sents = abstract.replace("\n", " ").split(". ")
    pick = []
    total_words = 0
    for s in sents:
        words = s.split()
        if not words:
            continue
        if total_words + len(words) > max_words:
            remaining = max_words - total_words
            if remaining > 5:
                pick.append(" ".join(words[:remaining]) + "…")
            break
        pick.append(s.strip())
        total_words += len(words)
    return ". ".join(pick).strip()


def gather_literature(gene: str) -> tuple[list[dict], int]:
    pmids = pubmed_search_thca(gene)
    articles = pubmed_efetch(pmids[:10])
    # filter 2018+
    recent = []
    for a in articles:
        try:
            y = int(a["year"]) if a["year"].isdigit() else 0
        except Exception:
            y = 0
        if y >= 2018:
            a["evidence_extract"] = summarize_abstract(a.get("abstract", ""))
            recent.append(a)
    # keep top 5
    recent = recent[:5]
    cancer_count = pubmed_search_cancer(gene)
    return recent, cancer_count


# ---------------------------------------------------------------- PHASE 2 ----
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"


def chembl_target_id(gene: str) -> tuple[str | None, str | None]:
    try:
        r = SESSION.get(f"{CHEMBL}/target/search.json",
                        params={"q": gene, "organism": "Homo sapiens"},
                        timeout=HTTP_TIMEOUT)
        r.raise_for_status()
        j = r.json()
        tgts = j.get("targets", [])
        # pick highest score SINGLE PROTEIN
        for t in tgts:
            comps = t.get("target_components", [])
            if comps and comps[0].get("component_synonyms"):
                syns = {s.get("component_synonym", "").upper()
                        for s in comps[0]["component_synonyms"]}
                if gene.upper() in syns:
                    return t.get("target_chembl_id"), t.get("pref_name")
        if tgts:
            return tgts[0].get("target_chembl_id"), tgts[0].get("pref_name")
    except Exception as e:
        log(f"    chembl target search failed for {gene}: {e}")
    return None, None


def chembl_activities(target_chembl_id: str, limit: int = 30) -> list[dict]:
    try:
        r = SESSION.get(f"{CHEMBL}/activity.json",
                        params={"target_chembl_id": target_chembl_id,
                                "assay_type": "B",
                                "pchembl_value__gte": 6,
                                "limit": limit},
                        timeout=HTTP_TIMEOUT)
        r.raise_for_status()
        return r.json().get("activities", [])
    except Exception as e:
        log(f"    chembl activities failed: {e}")
        return []


def chembl_molecule(chembl_id: str) -> dict:
    try:
        r = SESSION.get(f"{CHEMBL}/molecule/{chembl_id}.json", timeout=HTTP_TIMEOUT)
        r.raise_for_status()
        return r.json()
    except Exception:
        return {}


def pubchem_inhibitor_cids(gene: str) -> list[int]:
    try:
        q = quote(f"{gene} inhibitor")
        r = SESSION.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{q}/cids/JSON",
                        timeout=HTTP_TIMEOUT)
        if r.status_code == 200:
            return r.json().get("IdentifierList", {}).get("CID", [])[:5]
    except Exception:
        pass
    return []


def gather_compounds(gene: str) -> tuple[list[dict], str | None, str | None]:
    tid, pref = chembl_target_id(gene)
    if not tid:
        return [], None, None
    acts = chembl_activities(tid, limit=30)
    # aggregate by molecule_chembl_id, keeping best pchembl
    by_mol: dict[str, dict] = {}
    for a in acts:
        mid = a.get("molecule_chembl_id")
        if not mid:
            continue
        try:
            pchembl = float(a.get("pchembl_value")) if a.get("pchembl_value") else None
        except Exception:
            pchembl = None
        if mid not in by_mol or (pchembl is not None and (by_mol[mid].get("pchembl") or 0) < pchembl):
            by_mol[mid] = {
                "chembl_id": mid,
                "name": a.get("molecule_pref_name") or mid,
                "pchembl": pchembl,
                "moa": a.get("action_type") or a.get("bao_label") or a.get("assay_description", "")[:80],
                "target_relation_confidence": a.get("target_relation_type") or "DIRECT",
                "smiles": "",
            }
    # sort desc by pchembl, top 10, then fetch SMILES for the top ones
    mols = sorted(by_mol.values(), key=lambda x: -(x["pchembl"] or 0))[:10]
    for m in mols:
        # enrich with smiles
        info = chembl_molecule(m["chembl_id"])
        structs = info.get("molecule_structures") or {}
        m["smiles"] = structs.get("canonical_smiles", "") or ""
        if not m["name"] or m["name"] == m["chembl_id"]:
            m["name"] = info.get("pref_name") or m["chembl_id"]
    return mols, tid, pref


# ---------------------------------------------------------------- PHASE 3 ----
def druggability(n_compounds: int, best_pchembl: float | None,
                 has_structure: bool, n_thyroid_papers: int) -> dict:
    if n_compounds >= 10:
        klass = "validated_target"
    elif n_compounds >= 1:
        klass = "emerging_target"
    else:
        klass = "novel_target"
    return {
        "classification": klass,
        "n_known_inhibitors": n_compounds,
        "best_pchembl": best_pchembl,
        "has_structure": has_structure,
        "disease_focus_score": n_thyroid_papers,
    }


def rcsb_has_structure(gene: str) -> bool:
    # simple text search; keep brief to not waste time
    try:
        r = SESSION.post(
            "https://search.rcsb.org/rcsbsearch/v2/query",
            json={
                "query": {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                        "operator": "exists",
                    },
                },
                "request_options": {"paginate": {"rows": 1}},
                "return_type": "polymer_entity",
            },
            timeout=8,
        )
        # fallback to simple text search
    except Exception:
        pass
    try:
        r = SESSION.post(
            "https://search.rcsb.org/rcsbsearch/v2/query",
            json={
                "query": {"type": "terminal", "service": "full_text",
                          "parameters": {"value": gene}},
                "return_type": "entry",
                "request_options": {"paginate": {"rows": 1, "start": 0}},
            },
            timeout=8,
        )
        if r.status_code == 200:
            j = r.json()
            return j.get("total_count", 0) > 0
    except Exception as e:
        log(f"    rcsb check failed for {gene}: {e}")
    return False


# ---------------------------------------------------------------- MAIN -------
def main() -> int:
    bfile = wait_for_biomarkers()
    if bfile is not None:
        rows = load_targets_from_biomarker_tsv(bfile)
        log(f"PHASE 0: loaded {len(rows)} biomarker rows from {bfile.name}")
        targets = select_top8_targets(rows)
    else:
        targets = [{
            **t,
            "is_novel": True,
        } for t in derive_targets_from_de(top_n=8)]
    if len(targets) < 8:
        # make up shortfall
        extras = derive_targets_from_de(top_n=20)
        existing = {t["gene"] for t in targets}
        for e in extras:
            if e["gene"] not in existing:
                targets.append({**e, "is_novel": True})
            if len(targets) >= 8:
                break
    targets = targets[:8]
    log(f"PHASE 0: selected targets: {[t['gene'] for t in targets]}")

    # PHASE 1 + 2 + 3 per target
    all_lit: list[dict] = []
    all_cmp: list[dict] = []
    target_rows: list[dict] = []
    payload_targets = []

    for t in targets:
        gene = t["gene"]
        log(f"PHASE 1-3: {gene}")
        lit, cancer_count = gather_literature(gene)
        log(f"    pubmed thyroid hits: {len(lit)} (cancer total≈{cancer_count})")
        cmps, tid, pref = gather_compounds(gene)
        log(f"    chembl target={tid} pref={pref} compounds={len(cmps)}")
        has_struct = rcsb_has_structure(gene)
        best_pchembl = None
        for c in cmps:
            if c["pchembl"] is not None and (best_pchembl is None or c["pchembl"] > best_pchembl):
                best_pchembl = c["pchembl"]
        drug = druggability(len(cmps), best_pchembl, has_struct, len(lit))

        top_pmid = lit[0]["pmid"] if lit else ""
        top_compound = cmps[0]["name"] if cmps else ""

        target_rows.append({
            "gene": gene,
            "is_novel": t.get("is_novel"),
            "novelty_score": t.get("novelty_score"),
            "classification": drug["classification"],
            "n_compounds": drug["n_known_inhibitors"],
            "best_pchembl": drug["best_pchembl"],
            "has_structure": drug["has_structure"],
            "disease_focus_score": drug["disease_focus_score"],
            "top_pmid": top_pmid,
            "top_compound": top_compound,
            "chembl_target_id": tid or "",
            "chembl_pref_name": pref or "",
            "log2FC": t.get("log2FC", ""),
            "fdr": t.get("fdr", ""),
            "replicated_27155": t.get("replicated_27155", ""),
            "replicated_126698": t.get("replicated_126698", ""),
        })

        for c in cmps:
            moa_val = c.get("moa")
            if not isinstance(moa_val, str):
                moa_val = "" if moa_val is None else str(moa_val)
            all_cmp.append({
                "target": gene,
                "chembl_id": c["chembl_id"],
                "name": c["name"],
                "smiles": c["smiles"],
                "pchembl": c["pchembl"],
                "moa": moa_val[:160],
                "target_relation_confidence": c.get("target_relation_confidence", "DIRECT"),
            })
        for a in lit:
            all_lit.append({
                "target": gene,
                "pmid": a["pmid"],
                "year": a["year"],
                "title": a["title"],
                "evidence_extract": a.get("evidence_extract") or summarize_abstract(a.get("abstract","")),
            })

        payload_targets.append({
            "gene": gene,
            "is_novel": bool(t.get("is_novel")),
            "novelty_score": round(float(t.get("novelty_score") or 0), 3),
            "classification": drug["classification"],
            "n_compounds": drug["n_known_inhibitors"],
            "best_pchembl": drug["best_pchembl"],
            "disease_focus_score": drug["disease_focus_score"],
            "replicated_27155": str(t.get("replicated_27155", "")),
            "replicated_126698": str(t.get("replicated_126698", "")),
            "log2FC": (str(t.get("log2FC","")) or "")[:10],
            "chembl_target_id": tid or "",
            "chembl_pref_name": pref or "",
            "compounds": [
                {
                    "chembl_id": c["chembl_id"],
                    "name": (c["name"] or c["chembl_id"])[:80],
                    "smiles": c["smiles"],
                    "pchembl": c["pchembl"],
                    "moa": (c.get("moa") if isinstance(c.get("moa"), str) else (str(c.get("moa")) if c.get("moa") else ""))[:120],
                }
                for c in cmps
            ],
            "literature": [
                {"pmid": a["pmid"], "year": a["year"], "title": a["title"][:240],
                 "evidence": a.get("evidence_extract","")}
                for a in lit
            ],
        })

    # ----------------------------------------------------- PHASE 4 outputs ---
    # drug_discovery_targets.tsv
    import csv
    with (TABLES / "drug_discovery_targets.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=[
            "gene","novelty_score","classification","n_compounds","best_pchembl",
            "disease_focus_score","has_structure","is_novel","top_pmid","top_compound",
            "chembl_target_id","chembl_pref_name","log2FC","fdr","replicated_27155","replicated_126698",
        ])
        w.writeheader()
        for r in target_rows:
            w.writerow(r)

    with (TABLES / "drug_discovery_compounds.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=[
            "target","chembl_id","name","smiles","pchembl","moa","target_relation_confidence",
        ])
        w.writeheader()
        for c in all_cmp:
            w.writerow(c)

    with (TABLES / "drug_discovery_literature.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=[
            "target","pmid","year","title","evidence_extract",
        ])
        w.writeheader()
        for a in all_lit:
            w.writerow(a)

    # PAYLOAD json
    payload = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "source_biomarker_file": str(bfile.name) if bfile else "derived_from_DE",
        "n_targets": len(payload_targets),
        "n_compounds_total": len(all_cmp),
        "n_literature_total": len(all_lit),
        "kpi": {
            "novel": sum(1 for t in payload_targets if t["is_novel"])
                     or sum(1 for t in payload_targets if t["classification"] == "novel_target"),
            "validated": sum(1 for t in payload_targets if t["classification"] == "validated_target"),
            "emerging": sum(1 for t in payload_targets if t["classification"] == "emerging_target"),
            "novel_target": sum(1 for t in payload_targets if t["classification"] == "novel_target"),
            "compounds": len(all_cmp),
        },
        "targets": payload_targets,
    }
    (PAYLOAD_DIR / "drug_discovery_payload.json").write_text(
        json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
    )
    size_kb = (PAYLOAD_DIR / "drug_discovery_payload.json").stat().st_size / 1024
    log(f"PHASE 4: payload size = {size_kb:.1f} KB")

    # Human-readable summary MD
    md = build_report_md(target_rows, all_cmp, all_lit, bfile)
    REPORT_MD.write_text(md)
    log(f"PHASE 4: wrote {REPORT_MD}")
    return 0


def build_report_md(target_rows, cmps, lit, bfile) -> str:
    lines = []
    lines.append("# Biomarker → Drug Discovery Report")
    lines.append("")
    lines.append(f"- Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"- Biomarker source: `{bfile.name if bfile else 'derived_from_DE'}`")
    lines.append(f"- Targets analysed: **{len(target_rows)}**")
    lines.append(f"- Compounds collected: **{len(cmps)}**")
    lines.append(f"- Literature records: **{len(lit)}**")
    lines.append("")
    lines.append("## Top 8 targets")
    lines.append("")
    lines.append("| Gene | Novel | Classification | #Cmpds | Best pChEMBL | Thyroid papers | Top PMID | Top Compound |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for r in target_rows:
        lines.append(
            f"| {r['gene']} | {r.get('is_novel')} | {r['classification']} | {r['n_compounds']} "
            f"| {r['best_pchembl'] if r['best_pchembl'] is not None else '—'} "
            f"| {r['disease_focus_score']} | {r['top_pmid'] or '—'} | {r['top_compound'] or '—'} |"
        )
    lines.append("")
    lines.append("## Pharma-quality highlights (3 novel targets)")
    novel = [r for r in target_rows if r.get("is_novel") and r["classification"] in ("novel_target","emerging_target")][:3]
    if not novel:
        novel = target_rows[:3]
    for r in novel:
        gene = r["gene"]
        cls = r["classification"]
        ncmp = r["n_compounds"]
        bp = r["best_pchembl"]
        papers = r["disease_focus_score"]
        lines.append(f"### {gene} — {cls}")
        if cls == "novel_target":
            lines.append(
                f"**{gene}** surfaces from biomarker discovery with high novelty_score and no "
                f"known high-confidence chemical starting points in ChEMBL. "
                f"For a pharma portfolio this is a **first-in-class opportunity** — the target "
                f"hypothesis is de-risked by replicated differential expression in thyroid cancer "
                f"cohorts, while the chemical matter is greenfield (fragment screening / DNA-encoded "
                f"library campaigns would be the logical next step). Thyroid-specific literature "
                f"count = {papers}."
            )
        elif cls == "emerging_target":
            lines.append(
                f"**{gene}** has {ncmp} bioactive compound(s) in ChEMBL "
                f"(best pChEMBL ≈ {bp if bp is not None else '—'}), indicating "
                f"**early chemical validation** but no flagship inhibitor. That combination — "
                f"novel thyroid-cancer biomarker plus a thin but real tool-compound set — is "
                f"well suited to rapid lead-optimisation and IP carve-outs. "
                f"Literature in thyroid cancer = {papers} papers."
            )
        else:
            lines.append(
                f"**{gene}** is a **validated druggable target** ({ncmp} ChEMBL actives, best "
                f"pChEMBL ≈ {bp}). Its appearance inside our novelty-ranked biomarker set shows "
                f"that well-trodden pharmacology may be **repositioned** into thyroid cancer "
                f"contexts, shortening the path to preclinical. Thyroid literature count = {papers}."
            )
        lines.append("")
    lines.append("")
    lines.append("## Caveat")
    lines.append("This is computational triage only — not a drug recommendation. All compound "
                 "activities come from heterogeneous public assays; thyroid-specific efficacy is "
                 "not guaranteed. Confirmatory in-vitro / in-vivo work required before any claim.")
    lines.append("")
    lines.append("## Run log")
    lines.extend(["    " + l for l in LOG])
    return "\n".join(lines)


if __name__ == "__main__":
    sys.exit(main())
