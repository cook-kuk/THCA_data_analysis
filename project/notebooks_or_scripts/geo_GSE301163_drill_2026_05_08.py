#!/usr/bin/env python3
"""Detailed drill-down on GSE301163 — DICER1/DGCR8 thyroid lesions spatial transcriptomics."""
import urllib.request, time, re
from pathlib import Path

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/geo_search_2026_05_08")

def fetch(gse, view="brief", targ="all"):
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&targ={targ}&form=text&view={view}"
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read().decode("utf-8", errors="ignore")

print("[fetch] GSE301163 series + samples")
text_brief = fetch("GSE301163", view="brief", targ="self")
print("--- SERIES (self/brief) ---")
print(text_brief[:5000])

print("\n--- SAMPLE LIST (gsm) ---")
text_samples = fetch("GSE301163", view="brief", targ="gsm")
gsm_blocks = re.split(r"(?=^\^SAMPLE = )", text_samples, flags=re.M)
print(f"Total sample blocks: {len(gsm_blocks)}")
sample_recs = []
for b in gsm_blocks:
    if not b.strip(): continue
    gsm = re.search(r"^\^SAMPLE = (GSM\d+)", b, flags=re.M)
    title = re.search(r"^!Sample_title\s*=\s*(.+)", b, flags=re.M)
    src = re.search(r"^!Sample_source_name_ch1\s*=\s*(.+)", b, flags=re.M)
    org = re.search(r"^!Sample_organism_ch1\s*=\s*(.+)", b, flags=re.M)
    char = re.findall(r"^!Sample_characteristics_ch1\s*=\s*(.+)", b, flags=re.M)
    plat = re.search(r"^!Sample_platform_id\s*=\s*(.+)", b, flags=re.M)
    desc = re.search(r"^!Sample_description\s*=\s*(.+)", b, flags=re.M)
    if gsm:
        sample_recs.append({
            "gsm": gsm.group(1),
            "title": (title.group(1) if title else "").strip(),
            "source": (src.group(1) if src else "").strip(),
            "organism": (org.group(1) if org else "").strip(),
            "platform": (plat.group(1) if plat else "").strip(),
            "characteristics": " | ".join([c.strip() for c in char]),
            "description": (desc.group(1) if desc else "").strip(),
        })

print(f"\nParsed: {len(sample_recs)} samples")
print("\n=== Platform breakdown ===")
plats = {}
for s in sample_recs:
    plats[s["platform"]] = plats.get(s["platform"], 0) + 1
for p, n in sorted(plats.items(), key=lambda x: -x[1]):
    print(f"  {p:14s}  n={n}")

print("\n=== First 30 samples (title / source / chars) ===")
for s in sample_recs[:30]:
    print(f"  {s['gsm']:12s}  [{s['platform']}]  {s['title'][:70]:70s}  | {s['characteristics'][:120]}")

# Save TSV
import csv
fout = OUT/"GSE301163_samples.tsv"
with open(fout, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    keys = ["gsm","title","source","organism","platform","characteristics","description"]
    w.writerow(keys)
    for s in sample_recs:
        w.writerow([s[k] for k in keys])
print(f"\n[output] {fout}")

# Identify which are spatial vs others
spatial_keys = ["spatial","visium","st-seq","st_seq","stseq"]
spatial_samples = [s for s in sample_recs if any(k in (s["title"]+s["source"]+s["description"]+s["characteristics"]).lower() for k in spatial_keys)]
print(f"\n=== Spatial-related samples in GSE301163: {len(spatial_samples)} ===")
for s in spatial_samples[:40]:
    print(f"  {s['gsm']:12s}  [{s['platform']}]  {s['title'][:80]}")

# Disease/lesion type breakdown
print(f"\n=== Spatial samples — characteristics breakdown ===")
char_counter = {}
for s in spatial_samples:
    for c in s["characteristics"].split("|"):
        c = c.strip()
        if c:
            char_counter[c] = char_counter.get(c, 0) + 1
for k, v in sorted(char_counter.items(), key=lambda x: -x[1])[:30]:
    print(f"  ({v:3d})  {k[:120]}")
