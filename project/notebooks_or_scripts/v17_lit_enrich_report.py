"""Generate consolidated paste-ready report from deep_analysis CSVs.

Inputs:  project/manuscript_p2_brief/lit_enrich_2026_05_02/deep_analysis/*.csv
Output:  project/manuscript_p2_brief/lit_enrich_2026_05_02/08_deep_analysis_REPORT.md
"""

from datetime import datetime
from pathlib import Path

import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
DA = ROOT / "project/manuscript_p2_brief/lit_enrich_2026_05_02/deep_analysis"
OUT = ROOT / "project/manuscript_p2_brief/lit_enrich_2026_05_02/08_deep_analysis_REPORT.md"

verify = pd.read_csv(DA / "01_verified_refs_triangulation.csv")
top5 = pd.read_csv(DA / "02_related_works_top5_per_claim.csv")
gap = pd.read_csv(DA / "06_gap_candidates_top3_per_claim.csv")
trials_all = pd.read_csv(DA / "03_trials_structured.csv")
trials_drug = pd.read_csv(DA / "03_trials_targeted_drugs.csv")
hla_w = pd.read_csv(DA / "04_hla_pop_weighted.csv")
hla_long = pd.read_csv(DA / "04_hla_long.csv")
up = pd.read_csv(DA / "05_unpaywall_summary.csv")
overlap = pd.read_csv(DA / "06_overlap_with_bib.csv")


def md_table(df: pd.DataFrame, max_col_width: int = 80) -> str:
    df = df.copy().fillna("")
    for c in df.columns:
        df[c] = df[c].astype(str).str.replace("\n", " ").str.replace("|", "/").str.slice(0, max_col_width)
    return df.to_markdown(index=False)


lines = [
    "# 08 lit_enrich DEEP analysis report",
    f"",
    f"_Generated: {datetime.now():%Y-%m-%d %H:%M}_",
    f"",
    "Builds on raw `data/*.json` (5.9 MB) — does dedupe + cross-source merge + ranking,",
    "outputs paper-ready tables. Does NOT replace `01-07_*.md` (those stay as the",
    "narrative summary); this is the **structured drill-down**.",
    "",
    "---",
    "",
    "## 1. Reference verification — triangulation matrix",
    "",
    f"_All 17 entries × 4 sources (CrossRef / PubMed / OpenAlex / EuropePMC). Sorted by_",
    f"_(incomplete-first, then n_confirm-desc). **{int(verify['incomplete'].sum())} incomplete entries**_",
    f"_remain; **{int((verify['n_confirm']>=1).sum())}/17 have ≥1 confirming source**._",
    "",
    md_table(verify[["key", "incomplete", "n_confirm", "crossref", "pubmed", "openalex", "europepmc", "crossref_title"]]),
    "",
    "---",
    "",
    "## 2. Competitive landscape — Top-5 per claim (cross-source dedupe)",
    "",
    f"_Mined {len(top5.groupby('claim'))} claims × 5 sources (OpenAlex/SemScholar/E-PMC/bioRxiv/arXiv);_",
    f"_dedupe collapsed 279 raw rows → {len(top5)} top-5 papers across {top5['claim'].nunique()} claims._",
    f"_Composite rank = `n_sources*100 + min(cited_by,500) + influential*5 + (year-2000)`._",
    "",
]

for claim, sub in top5.groupby("claim"):
    lines += [
        f"### {claim}",
        "",
        md_table(sub[["n_sources", "year", "first_author", "title", "cited_by", "influential", "doi"]]),
        "",
    ]

lines += [
    "---",
    "",
    "## 3. Bibliography GAP candidates — top-3 per claim NOT in our bib",
    "",
    "_These are high-rank papers from the dedupe that don't appear in our 17-entry bib._",
    "_Triage candidates for Discussion §3 expansion._",
    "",
]
for claim, sub in gap.groupby("claim"):
    lines += [
        f"### {claim}",
        "",
        md_table(sub[["n_sources", "year", "first_author", "title", "cited_by", "doi"]]),
        "",
    ]

lines += [
    "---",
    "",
    "## 4. HLA AFND — sample-weighted East Asian matrix (Pillar 1 input)",
    "",
    f"_From {len(hla_long)} (allele × population × study) raw records, computed sample-_",
    f"_weighted mean frequency per (allele, country) using AFND `sample_size` as weights._",
    "",
    "### Weighted matrix",
    "",
    md_table(hla_w),
    "",
    "### Pillar 1 forest input — DPB1*05:01 East Asian baseline",
    "",
    "| Country | Weighted freq | Total n | Studies | vs Korean PTC pool 53.2% |",
    "|---|---|---|---|---|",
]
dpb = hla_w[hla_w["allele"] == "DPB1*05:01"]
for _, row in dpb.iterrows():
    delta = 0.532 - row["weighted_freq"]
    lines.append(
        f"| {row['population']} | {row['weighted_freq']:.3f} ({row['weighted_freq']*100:.1f}%) | "
        f"{int(row['total_n'])} | {int(row['n_studies'])} | Δ = +{delta:.3f} |"
    )

lines += [
    "",
    "_**Interpretation:** Korean PTC pool 53.2% > Korean baseline 36.7% (Δ=+16.5pp,_",
    f"_pooled n={int(dpb[dpb['population']=='Korea']['total_n'].iloc[0])}). Forest meta input ready._",
    "",
    "---",
    "",
    "## 5. Clinical trials — structured drill-down",
    "",
    f"_Total {len(trials_all)} active/recruiting trials across 5 queries; **{len(trials_drug)} use targeted_",
    f"_drugs** (RET/BRAF/IO/TKI). Phase + intervention + sponsor extracted from CT.gov v2._",
    "",
    "### Phase distribution (all 65 trials)",
    "",
    md_table(trials_all["phase"].value_counts().reset_index().rename(
        columns={"phase": "phase", "count": "n_trials"})),
    "",
    "### Sponsor class (all 65 trials)",
    "",
    md_table(trials_all["sponsor_class"].value_counts().reset_index().rename(
        columns={"sponsor_class": "class", "count": "n_trials"})),
    "",
    "### Top targeted-drug trials (sorted: phase desc, recruiting first)",
    "",
]
sorted_drug = trials_drug.sort_values(
    ["status", "phase", "enrollment"], ascending=[True, False, False]
).head(20)
lines += [md_table(sorted_drug[["nct_id", "phase", "status", "enrollment", "intervention_names", "lead_sponsor"]]), ""]

lines += [
    "---",
    "",
    "## 6. Open-access full-text discovery (Unpaywall)",
    "",
    f"_{int(up['is_oa'].sum())}/{len(up)} entries have OA copy. PMC IDs + PDF URLs ready for Methods/Discussion citations._",
    "",
    md_table(up[["key", "is_oa", "oa_status", "pmcid", "pdf_url"]]),
    "",
    "---",
    "",
    "## 7. Cross-reference: dedup ∩ verified_bib",
    "",
    f"_{len(overlap)} papers appear in BOTH our bib AND the mined top-ranked list — these are_",
    f"_**self-validated citations** (independently surfaced as top by 5-source mining)._",
    "",
    md_table(overlap),
    "",
    "---",
    "",
    "## 8. Action items for manuscript_v8",
    "",
    "**Immediate (paste-ready):**",
    "1. Pillar 1 forest: replace placeholder Korean baseline with **36.7% (n=680, 3 studies)** — §4.",
    "2. Pan-Asian gradient panel: DPB1*05:01 Korea 36.7% / China 35.5% / Japan 38.6% / Taiwan 52.1%",
    "3. Discussion §3.x translational outlook: cite NCT06458036 (RAISE selpercatinib pre-RAI),",
    "   NCT06475989 (Phase 3 targeted vs chemo), NCT04675710 (pembro+dabra+trame neoadjuvant).",
    "",
    "**Triage (review gap_candidates_top3 per claim):**",
    "- 24 high-rank papers not yet in bib — accept/reject per claim.",
    "",
    "**Pending (next pass):**",
    "- Fetch full abstracts for 10 OA Unpaywall PDFs (currently link-only).",
    "- Sem Scholar 0-hit on 6/8 claims (rate-limited) — retry with backoff.",
    "",
]

OUT.write_text("\n".join(lines))
print(f"Wrote {OUT} ({OUT.stat().st_size/1024:.1f} KB, {len(lines)} lines)")
