#!/usr/bin/env python3
"""GEO drill-down on shortlisted candidates."""
from __future__ import annotations
from pathlib import Path
import urllib.request, urllib.parse, json, time, re
import xml.etree.ElementTree as ET

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT/"project/results/geo_search_2026_05_08"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

CANDIDATES = ["GSE301163","GSE231952","GSE231953","GSE308145","GSE305978",
              "GSE145958","GSE231952","GSE273927","GSE305186"]

def fetch_geo_text(gse):
    """Fetch GEO accession plain-text view."""
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&targ=self&form=text&view=quick"
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            return r.read().decode("utf-8", errors="ignore")
    except Exception as e:
        return f"ERROR: {e}"

def parse_kv(text):
    out = {}
    for line in text.splitlines():
        if line.startswith("!"):
            kv = line.split("=", 1)
            if len(kv) == 2:
                k = kv[0].strip("!").strip()
                v = kv[1].strip()
                out.setdefault(k, []).append(v)
    return out

import csv
fout = OUT/"geo_candidates_drilldown.tsv"
with open(fout, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["gse","title","summary_short","organism","platforms","n_samples",
                "type","supp_keywords","pubmed","status","decision_hint"])
    for gse in CANDIDATES:
        time.sleep(0.5)
        text = fetch_geo_text(gse)
        if text.startswith("ERROR"):
            print(f"  {gse}: {text}")
            continue
        kv = parse_kv(text)
        title = "; ".join(kv.get("Series_title", []))[:200]
        summary = "; ".join(kv.get("Series_summary", []))[:380]
        organism = "; ".join(kv.get("Series_sample_organism", kv.get("Sample_organism_ch1", [])))
        platforms = "; ".join(kv.get("Series_platform_id", []))
        n_samples = kv.get("Series_sample_id", [])
        n_samples_n = len(n_samples)
        stype = "; ".join(kv.get("Series_type", []))[:120]
        supp = "; ".join(kv.get("Series_supplementary_file", []))[:240]
        pubmed = "; ".join(kv.get("Series_pubmed_id", []))
        status = "; ".join(kv.get("Series_status", []))[:80]
        # Decision hint
        blob = (title + " " + summary + " " + supp + " " + organism + " " + stype).lower()
        is_human = "homo sapiens" in blob or "human" in blob
        is_mouse = "mus musculus" in blob or "mouse" in blob
        is_thyroid_ca = any(k in blob for k in ["papillary thyroid","follicular thyroid","anaplastic thyroid","thyroid carcinoma","ptc","thyroid cancer","thyroid lesion","thyroid tumor"])
        is_visium = "visium" in blob or "10x genomics" in blob or "spatial" in blob
        is_ffpe = "ffpe" in blob or "fixed" in blob
        hint = []
        if is_human: hint.append("HUMAN")
        if is_mouse: hint.append("MOUSE")
        if is_thyroid_ca: hint.append("THYROID_CA")
        if is_visium: hint.append("SPATIAL")
        if is_ffpe: hint.append("FFPE")
        decision = ",".join(hint) or "AMBIG"
        w.writerow([gse, title, summary, organism, platforms, n_samples_n,
                    stype, supp, pubmed, status, decision])
        print(f"\n{'='*78}\n{gse}  ({n_samples_n} samples)  → {decision}")
        print(f"  Title:    {title}")
        print(f"  Type:     {stype}")
        print(f"  Organism: {organism}")
        print(f"  Platform: {platforms}")
        print(f"  Status:   {status}")
        print(f"  Pubmed:   {pubmed}")
        print(f"  Summary:  {summary[:280]}...")
        print(f"  Supp:     {supp[:200]}")

print(f"\n[output] {fout}")
