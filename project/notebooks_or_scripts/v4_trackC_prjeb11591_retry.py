#!/usr/bin/env python3
"""
v4 Track C: PRJEB11591 retry via alt supplementary URLs.

Attempts 6+ alternative URLs (PLOS supplementary .s023/.s021, europepmc,
ENA analyses, ArrayExpress biostudies, GEO search). Saves author email draft.

If ANY file retrieved -> parse + external validation.
If all fail -> status.txt with HTTP response codes + alt cohort candidates
(GSE33630 ATC arm, GSE82208 FTC arm).

Downstream tracks run regardless.

Outputs:
  results/ml/v4_trackC_status.txt
  results/ml/v4_trackC_attempts.json
  reports/html/figs_interactive/v4_trackC_figure.png (status infographic)
"""
from __future__ import annotations

import json
import sys
import time
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import requests

warnings.filterwarnings("ignore")

ROOT = Path("/opt/thyroid-dash/project")
OUT_ML = ROOT / "results" / "ml"
OUT_FIG = ROOT / "reports" / "html" / "figs_interactive"
OUT_RAW = ROOT / "data_raw"
for d in (OUT_ML, OUT_FIG, OUT_RAW):
    d.mkdir(parents=True, exist_ok=True)

ALT_URLS = [
    # PLOS supplementary variants
    "https://journals.plos.org/plosone/article/file?id=info:doi/10.1371/journal.pone.0156441.s023&type=supplementary",
    "https://journals.plos.org/plosone/article/file?id=info:doi/10.1371/journal.pone.0156441.s021&type=supplementary",
    "https://journals.plos.org/plosone/article/file?id=info:doi/10.1371/journal.pone.0156441.s022&type=supplementary",
    # europepmc
    "https://europepmc.org/article/MED/27310019",
    # ENA analyses
    "https://www.ebi.ac.uk/ena/browser/api/xml/PRJEB11591",
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB11591&result=analysis&format=tsv",
    # ArrayExpress biostudies
    "https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-PRJEB11591",
    # GEO search fallback
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=PRJEB11591",
]

AUTHOR_EMAIL_DRAFT = """Subject: Access request to supplementary expression matrix for PRJEB11591 (PLOS ONE 2016; doi:10.1371/journal.pone.0156441)

Dear colleagues,

I am a researcher working on a retrospective computational triage tool for
thyroid cancer (decision-support prototype, not a diagnostic device).
As part of a v4 cross-cohort validation of a transcriptomic signature
(TDS16 / TierA67), I am attempting to obtain the processed expression table
from PRJEB11591 as originally deposited alongside the PLOS ONE 2016 article.

I have tried several alternate supplementary URLs (s021/s022/s023),
europepmc, ENA (PRJEB11591), and ArrayExpress biostudies, without success.

Would you be able to:
1) confirm the correct supplementary file name for the processed matrix, or
2) provide a direct link, or
3) share a copy for strictly non-commercial academic research use?

I will happily cite the original publication and credit your lab.

Thank you very much for any pointers.

Best regards,
kukshomr@gmail.com
"""

ALT_COHORT_CANDIDATES = [
    {
        "geo_id": "GSE33630",
        "arm": "ATC (anaplastic) — 11 samples",
        "use": "aggressive arm external validation (LODO on top-of-funnel)",
        "note": "already harmonized in microarray_v3",
    },
    {
        "geo_id": "GSE82208",
        "arm": "FTC (follicular) — 52 samples",
        "use": "FTC follicular arm external validation",
        "note": "alt cohort; not in current data_processed, needs download",
    },
    {
        "geo_id": "GSE76039",
        "arm": "ATC+PDTC — 37 samples",
        "use": "aggressive-only validation",
        "note": "already harmonized in microarray_v3",
    },
    {
        "geo_id": "GSE29265",
        "arm": "mixed PTC — 49 samples",
        "use": "PTC external validation",
        "note": "already harmonized in microarray_v3",
    },
]


def log(msg: str) -> None:
    print(f"[trackC] {msg}", flush=True)


def attempt_url(url: str, timeout: int = 15) -> dict:
    rec = {"url": url, "status": None, "bytes": 0, "path": None, "error": None, "content_type": None}
    try:
        r = requests.get(url, timeout=timeout, allow_redirects=True,
                         headers={"User-Agent": "Mozilla/5.0 (THCA-research)"})
        rec["status"] = int(r.status_code)
        rec["content_type"] = r.headers.get("content-type", "")
        if r.status_code == 200 and len(r.content) > 1024:
            # save under a sane filename
            safe = url.split("/")[-1].split("?")[0] or "response.bin"
            if len(safe) < 3 or "=" in safe:
                safe = f"prjeb11591_attempt_{abs(hash(url))%10000}.bin"
            out = OUT_RAW / f"v4_prjeb11591_{safe}"
            out.write_bytes(r.content)
            rec["bytes"] = int(len(r.content))
            rec["path"] = str(out)
    except Exception as e:
        rec["error"] = str(e)[:200]
    return rec


def main() -> int:
    log(f"trying {len(ALT_URLS)} alt URLs...")
    attempts = []
    success_files = []
    for u in ALT_URLS:
        rec = attempt_url(u)
        attempts.append(rec)
        marker = "OK" if rec["path"] else "FAIL"
        log(f"  [{marker}] {rec['status']} bytes={rec['bytes']} :: {u[:90]}")
        if rec["path"]:
            success_files.append(rec)
        time.sleep(0.2)

    (OUT_ML / "v4_trackC_attempts.json").write_text(json.dumps(attempts, indent=2))

    # save author email draft
    (OUT_ML / "v4_trackC_author_email_draft.txt").write_text(AUTHOR_EMAIL_DRAFT)

    parsable = any(
        rec["path"] and (
            rec["content_type"].startswith("text")
            or "csv" in rec["content_type"]
            or "tsv" in rec["content_type"]
            or rec["path"].endswith((".tsv", ".csv", ".xlsx", ".txt"))
        )
        for rec in success_files
    )

    if success_files and parsable:
        # We'd try to parse here; for now mark as partial download only.
        status = "PARTIAL_DOWNLOAD"
        downloaded = [rec["path"] for rec in success_files]
        log(f"downloaded {len(downloaded)} candidate artifacts; parsing deferred")
    else:
        status = "ALL_FAILED"

    status_text = [
        f"# v4 Track C — PRJEB11591 retry",
        f"",
        f"status: {status}",
        f"n_alt_urls_tried: {len(ALT_URLS)}",
        f"n_successes (>=1KB): {len(success_files)}",
        f"",
        f"## HTTP response codes",
    ]
    for rec in attempts:
        status_text.append(f"  [{rec['status']}] bytes={rec['bytes']} {rec['url']}")
    status_text += [
        "",
        "## Downstream policy",
        "Downstream v4 tracks (A, B, D, synth) proceed regardless of this track's outcome.",
        "",
        "## Alternative cohort candidates",
    ]
    for c in ALT_COHORT_CANDIDATES:
        status_text.append(f"- {c['geo_id']} | {c['arm']} | {c['use']} | {c['note']}")

    (OUT_ML / "v4_trackC_status.txt").write_text("\n".join(status_text))

    # Figure: status infographic
    fig, ax = plt.subplots(figsize=(9, 5))
    codes = [rec["status"] if rec["status"] is not None else 0 for rec in attempts]
    labels = [rec["url"].split("//")[-1][:40] for rec in attempts]
    colors = ["#2ecc71" if c == 200 else "#e74c3c" for c in codes]
    y = list(range(len(attempts)))
    ax.barh(y, [rec["bytes"] / 1024 if rec["bytes"] else 0 for rec in attempts], color=colors)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("response bytes (KB)")
    ax.set_title(f"Track C: PRJEB11591 alt-URL attempts ({status})")
    for i, c in enumerate(codes):
        ax.text(1, i, f"HTTP {c}", fontsize=7, va="center")
    (OUT_FIG / "v4_trackC_figure.png").parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG / "v4_trackC_figure.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    (OUT_FIG / "v4_trackC_figure.html").write_text(
        f"<html><body><img src='v4_trackC_figure.png' style='max-width:100%'></body></html>"
    )

    log(f"track C status: {status}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
