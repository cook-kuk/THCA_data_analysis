#!/usr/bin/env python3
"""Catalog gap analysis.

For each claim in data_registry/manifests/claim_matrix.csv: count how many of
the claim's allowed_datasets are VERIFIED-OK, VERIFIED-OK-XCANCER,
VERIFIED-MANUAL-NEEDED, or simply unverified. Flag claims that depend on
unverified-only datasets.

Also reports zone-level verification rates.

Output: data_registry/reports/gap_analysis_2026_05_08.md
"""
from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
MASTER = REPO / "data_registry" / "manifests" / "master_dataset_catalog.csv"
CLAIMS = REPO / "data_registry" / "manifests" / "claim_matrix.csv"
OUT = REPO / "data_registry" / "reports" / "gap_analysis_2026_05_08.md"


def status_of(notes: str) -> str:
    if not notes:
        return "unverified"
    if notes.startswith("[VERIFIED-OK-XCANCER"):
        return "ok-xcancer"
    if notes.startswith("[VERIFIED-OK]"):
        return "ok"
    if notes.startswith("[VERIFIED-MANUAL-NEEDED]"):
        return "manual-needed"
    if notes.startswith("[VERIFIED-WRONG"):
        return "wrong"
    if notes.startswith("[VERIFIED-MISSING]"):
        return "missing"
    return "unverified"


def split_dataset_list(s: str) -> list[str]:
    if not s or s.strip() in ("", "-"):
        return []
    parts = [p.strip() for p in re.split(r"[;,]", s) if p.strip()]
    return parts


def main() -> None:
    master = list(csv.DictReader(MASTER.open(encoding="utf-8")))
    by_id = {r["dataset_id"]: r for r in master}
    status_by_id = {r["dataset_id"]: status_of(r.get("notes", "")) for r in master}

    # Zone-level verification rates
    zone_status: dict[str, Counter] = defaultdict(Counter)
    for r in master:
        zone_status[r["zone"]][status_of(r.get("notes", ""))] += 1

    # Claim coverage
    claims = list(csv.DictReader(CLAIMS.open(encoding="utf-8")))
    claim_rows = []
    for c in claims:
        allowed = split_dataset_list(c["allowed_datasets"])
        # Strip trailing notes like "(when used for HLA-cancer claims)" inside dataset names
        allowed = [a.split(" (")[0] for a in allowed]
        # Resolve glob patterns like "PAPER1_CANCER_CORE/*"
        allowed = [a for a in allowed if not a.endswith("/*")]
        statuses = [status_by_id.get(a, "missing-from-master") for a in allowed]
        ct = Counter(statuses)
        claim_rows.append({
            "claim_id": c["claim_id"],
            "claim_text": c["claim_text"],
            "paper": c["paper"],
            "n_total": len(allowed),
            "n_ok": ct["ok"],
            "n_ok_xcancer": ct["ok-xcancer"],
            "n_manual": ct["manual-needed"],
            "n_unverified": ct["unverified"],
            "n_missing_from_master": ct["missing-from-master"],
            "claim_strength": c["claim_strength"],
            "needs_validation": c["needs_validation"],
            "weakly_supported": ct["ok"] + ct["ok-xcancer"] == 0,
            "datasets": allowed,
            "statuses": statuses,
        })

    # Build report
    lines = ["# Catalog gap analysis — 2026-05-08", ""]
    lines.append(f"**Master:** {MASTER.relative_to(REPO)} ({len(master)} rows)")
    lines.append(f"**Claim matrix:** {CLAIMS.relative_to(REPO)} ({len(claims)} claims)")
    lines.append("")

    # Zone verification rates
    lines.append("## 1. Zone-level verification rates")
    lines.append("")
    lines.append("| Zone | Total | OK | OK-XCancer | Manual | Unverified | OK% |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for zone in sorted(zone_status):
        ct = zone_status[zone]
        total = sum(ct.values())
        ok = ct["ok"]
        ok_x = ct["ok-xcancer"]
        manual = ct["manual-needed"]
        unverified = ct["unverified"]
        pct = (ok + ok_x) / total * 100 if total else 0
        lines.append(f"| `{zone}` | {total} | {ok} | {ok_x} | {manual} | {unverified} | {pct:.1f}% |")
    lines.append("")

    # Claim → support
    lines.append("## 2. Claim-by-claim support")
    lines.append("")
    lines.append("Each claim's `allowed_datasets` list is enumerated against master verification status. **WEAK** = no VERIFIED-OK or VERIFIED-OK-XCANCER datasets back the claim.")
    lines.append("")
    lines.append("| Claim | Paper | Strength | Total | OK | OK-XC | Manual | Unverified | Missing-from-master | Status |")
    lines.append("|---|---|---|---:|---:|---:|---:|---:|---:|---|")
    for r in claim_rows:
        flag = "⚠ WEAK" if r["weakly_supported"] else "✓"
        lines.append(
            f"| `{r['claim_id']}` | {r['paper']} | {r['claim_strength']} | "
            f"{r['n_total']} | {r['n_ok']} | {r['n_ok_xcancer']} | {r['n_manual']} | "
            f"{r['n_unverified']} | {r['n_missing_from_master']} | {flag} |"
        )
    lines.append("")

    # Detail per claim
    lines.append("## 3. Per-claim dataset support detail")
    lines.append("")
    for r in claim_rows:
        lines.append(f"### `{r['claim_id']}` — {r['claim_text']}")
        lines.append(f"_{r['paper']} · strength={r['claim_strength']} · validation={r['needs_validation']}_")
        lines.append("")
        if r["weakly_supported"]:
            lines.append("⚠ **WEAK** — no VERIFIED-OK datasets currently back this claim. All supporting datasets are unverified or only manual-lookup-needed.")
            lines.append("")
        for ds, st in zip(r["datasets"], r["statuses"]):
            sym = {"ok": "✅", "ok-xcancer": "✅(x-cancer)", "manual-needed": "⚪ manual", "unverified": "○ unverified", "missing-from-master": "❌ missing-from-master", "wrong": "❌ wrong", "missing": "❌ missing"}.get(st, st)
            lines.append(f"- {sym} `{ds}`")
        lines.append("")

    # Weakly-supported summary
    weak = [r for r in claim_rows if r["weakly_supported"]]
    if weak:
        lines.append("## 4. ⚠ WEAKLY SUPPORTED CLAIMS (priority fill list)")
        lines.append("")
        lines.append("These claims need at least one VERIFIED-OK supporting dataset before they should ride into Paper 1/2/3/4 prose.")
        lines.append("")
        for r in weak:
            lines.append(f"- `{r['claim_id']}` ({r['paper']}, strength={r['claim_strength']}): {r['claim_text']}")
        lines.append("")

    # MANUAL-NEEDED dataset roll-up (which manual-needed datasets are most cited?)
    cited_manual = Counter()
    for r in claim_rows:
        for ds, st in zip(r["datasets"], r["statuses"]):
            if st == "manual-needed":
                cited_manual[ds] += 1
    if cited_manual:
        lines.append("## 5. Top MANUAL-NEEDED datasets by claim-citation count")
        lines.append("")
        lines.append("These manual-needed entries underpin multiple claims; resolving their PMIDs / citations would unlock the most claim-support.")
        lines.append("")
        lines.append("| Dataset | Cited by # claims |")
        lines.append("|---|---:|")
        for ds, n in cited_manual.most_common(20):
            lines.append(f"| `{ds}` | {n} |")
        lines.append("")

    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {OUT.relative_to(REPO)}")
    print(f"  {len(claim_rows)} claims analyzed")
    print(f"  {len(weak)} weakly-supported")


if __name__ == "__main__":
    main()
