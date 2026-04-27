#!/usr/bin/env python3
"""v17 TERT promoter recovery — final attempt from 3 public sources.

Sources tried (in order):
  1. cBioPortal datahub-public · thca_tcga_pan_can_atlas_2018/data_mutations.txt
  2. Liu et al. 2017 JCO supplementary (TERT promoter status table)
  3. ICGC dcc.icgc.org · projects/THCA-US

For each: HTTP attempt + size + parsed-row count are logged.
If at least one source returns TERT-bearing rows, cross-tab with
existing v17 DM1/DM2 cluster assignment from
results/v17/tables/sample_master_v17_tert.tsv.
"""
import csv, json, re, sys, time
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

ROOT = Path("/opt/thyroid-dash/project")
RES  = ROOT / "results/v17_tert_recovery"
RES.mkdir(parents=True, exist_ok=True)

UA = {"User-Agent": "thyrai-v17-tert-recovery (kukshomr@gmail.com)"}

def fetch(url, max_bytes=80_000_000):
    """Return (status, bytes) — empty bytes on failure."""
    try:
        req = Request(url, headers=UA)
        with urlopen(req, timeout=60) as r:
            data = r.read(max_bytes)
            return r.status, data
    except HTTPError as e:
        return e.code, b""
    except URLError as e:
        return -1, str(e).encode()
    except Exception as e:
        return -2, str(e).encode()

# ---------- Source 1: cBioPortal datahub ----------
def try_cbio():
    print("[S1] cBioPortal datahub thca_tcga_pan_can_atlas_2018")
    urls = [
        "https://github.com/cBioPortal/datahub-public/raw/master/public/thca_tcga_pan_can_atlas_2018/data_mutations.txt",
        "https://raw.githubusercontent.com/cBioPortal/datahub/master/public/thca_tcga_pan_can_atlas_2018/data_mutations.txt",
        # Older filename pattern
        "https://github.com/cBioPortal/datahub/raw/master/public/thca_tcga_pan_can_atlas_2018/data_mutations.txt",
    ]
    for u in urls:
        print(f"  ▸ trying {u[:90]}…")
        status, body = fetch(u)
        print(f"    → HTTP {status}, {len(body)} B")
        if status == 200 and len(body) > 1000:
            (RES / "cbio_data_mutations_raw.txt").write_bytes(body)
            return body, u
    return b"", None

def parse_maf_for_tert(body):
    """Parse MAF / data_mutations.txt and extract TERT entries."""
    if not body: return []
    text = body.decode("utf-8", errors="ignore")
    lines = text.splitlines()
    # header line: skip leading '#' comments
    header = None; data_start = 0
    for i, l in enumerate(lines):
        if l.startswith("#"): continue
        header = l.split("\t"); data_start = i + 1; break
    if not header: return []
    cols_lc = [c.lower() for c in header]
    def col(name):
        for variant in (name, name.lower()):
            try: return cols_lc.index(variant.lower())
            except ValueError: continue
        return -1
    iH = col("Hugo_Symbol")
    iV = col("Variant_Classification")
    iC = col("Chromosome")
    iS = col("Start_Position") if col("Start_Position") >= 0 else col("Start_position")
    iSb = col("Tumor_Sample_Barcode")
    iRef = col("Reference_Allele")
    iAlt = col("Tumor_Seq_Allele2") if col("Tumor_Seq_Allele2") >= 0 else col("Variant_Allele")
    out = []
    for l in lines[data_start:]:
        if not l or l.startswith("#"): continue
        f = l.split("\t")
        if len(f) <= max(iH, iV, iSb): continue
        if iH < 0 or f[iH] != "TERT": continue
        out.append({
            "sample_barcode": f[iSb] if iSb >= 0 and iSb < len(f) else "",
            "variant_class": f[iV] if iV >= 0 and iV < len(f) else "",
            "chrom": f[iC] if iC >= 0 and iC < len(f) else "",
            "start": f[iS] if iS >= 0 and iS < len(f) else "",
            "ref": f[iRef] if iRef >= 0 and iRef < len(f) else "",
            "alt": f[iAlt] if iAlt >= 0 and iAlt < len(f) else "",
        })
    return out

# ---------- Source 2: Liu 2017 JCO supplementary ----------
def try_liu_jco():
    print("[S2] Liu et al. 2017 JCO supplementary")
    # The supplementary file is typically a docx / xlsx behind a doi-suppl URL.
    # We try several documented routes.
    urls = [
        "https://ascopubs.org/doi/suppl/10.1200/JCO.2016.71.5654/suppl_file/jco_2016_71_5654_TERT.xlsx",
        "https://ascopubs.org/doi/suppl/10.1200/JCO.2016.71.5654/suppl_file/Supplementary_Table.xlsx",
        # Most ASCO papers actually serve ds_001/002 etc.
        "https://ascopubs.org/doi/suppl/10.1200/JCO.2016.71.5654/suppl_file/ds_001.xlsx",
        "https://ascopubs.org/doi/suppl/10.1200/JCO.2016.71.5654/suppl_file/ds_002.xlsx",
    ]
    hits = []
    for u in urls:
        status, body = fetch(u, max_bytes=20_000_000)
        ok = status == 200 and len(body) > 8000
        print(f"  ▸ {u[-50:]:>50} → HTTP {status}, {len(body)} B {'OK' if ok else ''}")
        if ok:
            (RES / f"liu2017_{u.rsplit('/',1)[-1]}").write_bytes(body)
            hits.append(u)
    return hits

# ---------- Source 3: ICGC THCA-US ----------
def try_icgc():
    print("[S3] ICGC THCA-US")
    urls = [
        "https://dcc.icgc.org/api/v1/projects/THCA-US",
        # Mutations download endpoint (typically gated; try anyway)
        "https://dcc.icgc.org/api/v1/projects/THCA-US/mutations?size=10",
    ]
    out = []
    for u in urls:
        status, body = fetch(u, max_bytes=2_000_000)
        ok = status == 200 and len(body) > 200
        print(f"  ▸ {u[-50:]:>50} → HTTP {status}, {len(body)} B {'OK' if ok else ''}")
        if ok:
            try:
                j = json.loads(body.decode())
                out.append((u, j))
            except json.JSONDecodeError:
                out.append((u, body.decode("utf-8", errors="ignore")[:1500]))
    return out

# ---------- run all 3 ----------
def main():
    t0 = time.time()
    log = []

    cbio_body, cbio_url = try_cbio()
    cbio_rows = parse_maf_for_tert(cbio_body)
    log.append({"source":"cBioPortal datahub", "url": cbio_url or "all 3 URLs failed",
                "bytes": len(cbio_body), "tert_rows": len(cbio_rows),
                "ok": len(cbio_rows) > 0})

    liu_hits = try_liu_jco()
    log.append({"source":"Liu 2017 JCO suppl", "url": ", ".join(liu_hits) or "all candidate URLs failed",
                "bytes": sum((RES / f"liu2017_{u.rsplit('/',1)[-1]}").stat().st_size for u in liu_hits) if liu_hits else 0,
                "tert_rows": "n/a (xlsx — manual parse step needed)" if liu_hits else 0,
                "ok": bool(liu_hits)})

    icgc_hits = try_icgc()
    log.append({"source":"ICGC THCA-US", "url": ", ".join(u for u,_ in icgc_hits) or "all candidate URLs failed",
                "bytes": "n/a", "tert_rows": 0, "ok": bool(icgc_hits)})

    # Save attempt TSV
    with (RES / "tert_recovery_attempt.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["source","url","bytes","tert_rows","ok"], delimiter="\t")
        w.writeheader()
        for r in log: w.writerow(r)

    # If we got cbio TERT rows → save them + crosstab
    crosstab_path = RES / "tert_dm_crosstab.tsv"
    if cbio_rows:
        # Write raw TERT rows
        with (RES / "tert_mutations_cbio.tsv").open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(cbio_rows[0].keys()), delimiter="\t")
            w.writeheader()
            for r in cbio_rows: w.writerow(r)

        # Promoter-region filter (chr5 1294497-1295162 is canonical TERT promoter; C228T = 1295228, C250T = 1295250 in hg19; pan-cancer atlas typically uses GRCh37)
        promoter_classes = {"5'Flank", "5'UTR", "Promoter", "5_prime_flanking_variant"}
        promoter_rows = [r for r in cbio_rows
                         if (r["variant_class"] in promoter_classes)
                         or (r["chrom"] in ("5", "chr5") and r["start"].isdigit() and 1294000 <= int(r["start"]) <= 1296000)]
        # Build TCGA short barcode (TCGA-XX-XXXX)
        def short(b):
            parts = b.split("-")
            return "-".join(parts[:3]) if len(parts) >= 3 else b
        tert_pos_cases = {short(r["sample_barcode"]) for r in promoter_rows}

        # Cross-tab with v17 DM1/DM2 cluster
        sm = ROOT / "results/v17/tables/sample_master_v17_tert.tsv"
        cross = {("DM1","TERT+"): 0, ("DM1","TERT-"): 0,
                 ("DM2","TERT+"): 0, ("DM2","TERT-"): 0,
                 ("other","TERT+"): 0, ("other","TERT-"): 0}
        n_total = 0
        if sm.exists():
            with sm.open() as f:
                r = csv.DictReader(f, delimiter="\t")
                for row in r:
                    sid = row.get("sample_id","")
                    case = "-".join(sid.split("-")[:3])
                    dm = row.get("dm_like","") or "other"
                    if dm not in ("DM1_like","DM2_like"): cluster = "other"
                    else: cluster = "DM1" if dm == "DM1_like" else "DM2"
                    tert = "TERT+" if case in tert_pos_cases else "TERT-"
                    cross[(cluster, tert)] += 1
                    n_total += 1
        with crosstab_path.open("w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["cluster","TERT_status","n"])
            for (c,t),n in cross.items():
                w.writerow([c, t, n])
        cross_summary = {f"{c}_{t}": n for (c,t), n in cross.items()}
        cross_summary["n_total_cross_tabbed"] = n_total
        cross_summary["n_unique_tert_promoter_cases"] = len(tert_pos_cases)
        cross_summary["n_raw_tert_rows"] = len(cbio_rows)
        cross_summary["n_promoter_filtered_rows"] = len(promoter_rows)
    else:
        cross_summary = {"note": "no TERT rows recovered → cross-tab skipped"}

    # Status markdown
    lines = ["# v17 TERT promoter recovery — status", "",
             f"_Run: {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}_  ",
             f"_Wall: {round(time.time()-t0,1)} s_", "",
             "## Source-by-source result", ""]
    for r in log:
        lines.append(f"### {r['source']}")
        lines.append(f"- URL(s): `{r['url'][:200]}`")
        lines.append(f"- Bytes: {r['bytes']}")
        lines.append(f"- TERT rows: {r['tert_rows']}")
        lines.append(f"- OK: **{'✅' if r['ok'] else '❌'}**")
        lines.append("")
    lines.append("## Cross-tab with v17 DM1/DM2")
    lines.append("```json")
    lines.append(json.dumps(cross_summary, indent=2, ensure_ascii=False))
    lines.append("```")
    lines.append("")
    lines.append("## Verdict")
    if cbio_rows:
        lines.append("**SUCCESS** — cBioPortal datahub returned a TERT mutation MAF. Cross-tab written to `tert_dm_crosstab.tsv`. Next step: integrate `tert_mutations_cbio.tsv` into `sample_master_v17_tert.tsv` (column `tert_promoter`) and refresh downstream survival / pathway analyses.")
    elif liu_hits:
        lines.append("**PARTIAL** — Liu 2017 supplementary file fetched but is xlsx; manual openpyxl parse step required to extract TERT promoter status. Cross-tab not yet computed.")
    elif any(r["ok"] for r in log):
        lines.append("**PARTIAL** — at least one source responded but no TERT rows were extractable in this run. See logs above.")
    else:
        lines.append("**FAILURE** — all 3 sources unreachable / returned no parseable TERT data. Recovery attempt is exhausted; downstream analyses must remain TERT-status-blinded for now.")
    (RES / "recovery_status.md").write_text("\n".join(lines), encoding="utf-8")

    idx = {"built_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "sources_tried": 3,
           "sources_succeeded": sum(1 for r in log if r["ok"]),
           "tert_rows_recovered": len(cbio_rows),
           "crosstab_written": cbio_rows and crosstab_path.exists(),
           "elapsed_seconds": round(time.time()-t0, 1),
           "log": log,
           "cross_summary": cross_summary}
    (RES / "v17_tert_index.json").write_text(json.dumps(idx, indent=2, ensure_ascii=False))
    print("[done]", json.dumps({k:v for k,v in idx.items() if k!="log"}, indent=2))

if __name__ == "__main__": main()
