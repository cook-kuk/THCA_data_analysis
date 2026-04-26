#!/usr/bin/env python3
"""v12 SNUBH validation shortlist —
Cross-references top-20 novel candidates × per-target clinical translatability
× per-target PubMed depth to produce a prospective-validation shortlist for
the SNUBH thyroid-cancer cohort collaboration (유형원 교수).

Output:
  results/v12_literature/snubh_validation/snubh_shortlist.tsv
  results/v12_literature/snubh_validation/snubh_shortlist_rationale.md
"""
import csv, json, math
from pathlib import Path

ROOT = Path("/opt/thyroid-dash/project")
RES  = ROOT / "results/v12_literature"
OUT  = RES / "snubh_validation"
OUT.mkdir(parents=True, exist_ok=True)

# ---- inputs ----------------------------------------------------------------
top20 = list(csv.DictReader((RES / "crosscheck/top_20_novel_deep_dive.tsv").open(), delimiter="\t"))
trans = {r["target"]: r for r in csv.DictReader(
            (RES / "crosscheck/clinical_translatability.tsv").open(), delimiter="\t")}
targs = {r["target"]: r for r in csv.DictReader(
            (RES / "pubmed_evidence/targets_literature_summary.tsv").open(), delimiter="\t")}
biom_full = {r["gene"]: r for r in csv.DictReader(
            (ROOT / "results/tables/biomarker_validated.tsv").open(), delimiter="\t")}

# 8 druggable target list — these have v7 chemistry support
DRUGGABLE_8 = {"CYP1B1","TACSTD2","TMPRSS4","PLEKHA6","LDLR","GABRB2","B3GNT3","PTPRE"}

# ---- score function --------------------------------------------------------
def as_float(x, default=0.0):
    try: return float(x)
    except: return default

def replication_fraction(x):
    """Normalize external replication count to 0..1 across 2 cohorts."""
    return min(max(as_float(x) / 2.0, 0.0), 1.0)

def format_fdr(x):
    try:
        return f"{float(x):.1e}"
    except Exception:
        return str(x) if x is not None else ""

def snubh_score(row):
    """Composite priority for prospective validation in a clinical cohort.
    Weights chosen to reward (a) effect-size strength, (b) replication,
    (c) novelty (under-investigated = more cohort value), and (d) druggability."""
    g = row["gene"]
    novelty_pubmed = 1.0 if row["classification"] == "novel" else (0.5 if row["classification"]=="emerging" else 0.0)
    cohens_d = abs(as_float(biom_full.get(g,{}).get("cohens_d_tcga",0)))
    rep_rate = replication_fraction(biom_full.get(g,{}).get("replication_rate",0))
    druggable = 1.0 if g in DRUGGABLE_8 else 0.0
    braf_pmids = as_float(row.get("pubmed_braf",0))
    # composite
    s = (0.30 * min(cohens_d/2.5, 1.0)
       + 0.20 * rep_rate
       + 0.20 * novelty_pubmed
       + 0.20 * druggable
       + 0.10 * min(braf_pmids/5, 1.0))
    return round(s, 3)

# ---- build rows ------------------------------------------------------------
rows = []
for r in top20:
    g = r["gene"]
    bio = biom_full.get(g, {})
    cli = trans.get(g, {})
    tgt = targs.get(g, {})
    score = snubh_score(r)
    rows.append({
        "rank_input":      r["rank"],
        "gene":            g,
        "snubh_priority_score": score,
        "log2FC_tcga":     r.get("log2FC_tcga",""),
        "cohens_d_tcga":   bio.get("cohens_d_tcga",""),
        "fdr_tcga":        bio.get("fdr_tcga",""),
        "replication_rate":bio.get("replication_rate",""),
        "pubmed_thyroid":  r.get("pubmed_thyroid",""),
        "pubmed_braf":     r.get("pubmed_braf",""),
        "novelty_class":   r.get("classification",""),
        "is_v7_druggable": "yes" if g in DRUGGABLE_8 else "no",
        "ctgov_thyroid_trials": cli.get("thyroid_trials",""),
        "ctgov_any_trials":     cli.get("any_cancer_trials",""),
        "translatability_tier": cli.get("priority",""),
        "rationale":       "",  # filled below
    })

# Per-gene one-line rationale
def rationale(r):
    parts = []
    rep_count = int(as_float(r["replication_rate"]))
    if r["is_v7_druggable"]=="yes":   parts.append("v7-druggable")
    if r["novelty_class"]=="novel":   parts.append("PubMed-novel (<5 thyroid)")
    elif r["novelty_class"]=="emerging": parts.append("emerging (5-24 thyroid)")
    if r["translatability_tier"]=="high": parts.append("high translatability")
    if int(as_float(r["pubmed_braf"])) > 0:  parts.append(f"{int(as_float(r['pubmed_braf']))} BRAF-context paper(s)")
    if abs(as_float(r["cohens_d_tcga"])) > 1.5: parts.append(f"large effect (|d|={abs(as_float(r['cohens_d_tcga'])):.2f})")
    if rep_count >= 2:
        parts.append("fully replicated (2/2 cohorts)")
    elif rep_count == 1:
        parts.append("partially replicated (1/2 cohorts)")
    return " · ".join(parts) if parts else "candidate"

for r in rows:
    r["rationale"] = rationale(r)

# Sort by snubh score desc
rows.sort(key=lambda r: -r["snubh_priority_score"])

# ---- write TSV -------------------------------------------------------------
tsv_path = OUT / "snubh_shortlist.tsv"
with tsv_path.open("w", newline="") as out_f:
    w = csv.DictWriter(out_f, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    for r in rows: w.writerow(r)

# ---- write rationale MD ----------------------------------------------------
md = ["# SNUBH prospective validation shortlist",
      "",
      "_Build: 2026-04-25 · 20 candidates from v12 top-20 novel deep dive,",
      "ranked by composite SNUBH validation priority score._",
      "",
      "## Scoring",
      "",
      "`snubh_priority_score = 0.30·min(|d|/2.5,1) + 0.20·rep_rate + 0.20·novelty + 0.20·druggable + 0.10·min(BRAF_pp/5,1)`",
      "",
      "- `|d|` = absolute Cohen's d in TCGA THCA (BRAF vs RAS).",
      "- `rep_rate` = fraction of external cohorts (GSE27155, GSE126698) replicating direction.",
      "- `novelty` = 1.0 (PubMed-novel <5 thyroid papers), 0.5 (emerging 5-24), 0 (established ≥25).",
      "- `druggable` = 1.0 if in the v7 8-target druggable shortlist, else 0.",
      "- `BRAF_pp` = PubMed papers co-mentioning the gene with BRAF (any cancer).",
      "",
      "## Top tier (score ≥ 0.65) — first wave for SNUBH validation",
      "",
      "| Rank | Gene | Score | log2FC | Cohen's d | FDR | Replicated | Novelty | Drug? | Rationale |",
      "|---:|---|---:|---:|---:|---|---:|---|---|---|"]

for i, r in enumerate(rows, 1):
    if r["snubh_priority_score"] < 0.65: continue
    fdr = r["fdr_tcga"]
    fdr_s = format_fdr(fdr) if fdr else ""
    md.append(f"| {i} | **{r['gene']}** | {r['snubh_priority_score']:.3f} | "
              f"{as_float(r['log2FC_tcga']):+.2f} | {as_float(r['cohens_d_tcga']):+.2f} | {fdr_s} | "
              f"{int(as_float(r['replication_rate']))} | {r['novelty_class']} | "
              f"{r['is_v7_druggable']} | {r['rationale']} |")

md += ["",
       "## Second tier (0.50 ≤ score < 0.65)",
       "",
       "| Rank | Gene | Score | Novelty | Drug? | Rationale |",
       "|---:|---|---:|---|---|---|"]
for i, r in enumerate(rows, 1):
    if not (0.50 <= r["snubh_priority_score"] < 0.65): continue
    md.append(f"| {i} | {r['gene']} | {r['snubh_priority_score']:.3f} | "
              f"{r['novelty_class']} | {r['is_v7_druggable']} | {r['rationale']} |")

md += ["",
       "## Third tier (score < 0.50)",
       "",
       "| Rank | Gene | Score | Notes |",
       "|---:|---|---:|---|"]
for i, r in enumerate(rows, 1):
    if r["snubh_priority_score"] >= 0.50: continue
    md.append(f"| {i} | {r['gene']} | {r['snubh_priority_score']:.3f} | {r['rationale']} |")

# Per-gene one-paragraph synopsis for the top 5
md += ["", "## Top-5 detailed synopses for clinical-collaborator briefing", ""]
for i, r in enumerate(rows[:5], 1):
    g = r["gene"]; bio = biom_full.get(g,{}); cli = trans.get(g,{}); tgt = targs.get(g,{})
    cd = abs(as_float(bio.get("cohens_d_tcga",0)))
    md += [f"### {i}. {g}  (SNUBH score {r['snubh_priority_score']:.3f})", ""]
    if g in DRUGGABLE_8:
        md.append(f"- **Druggable**: yes (v7 chemistry support; "
                  f"target literature: {tgt.get('n_thyroid_papers','?')} thyroid / "
                  f"{tgt.get('n_cancer_papers','?')} cancer / "
                  f"{tgt.get('n_braf_papers','?')} BRAF papers; "
                  f"druggability flag = {tgt.get('druggable_evidence','?')}).")
    else:
        md.append(f"- **Druggable**: no v7 chemistry hit, but biomarker-grade signal.")
    md.append(f"- **TCGA effect**: log2FC = {as_float(bio.get('log2FC_tcga',0)):+.2f}, "
              f"Cohen's d = {as_float(bio.get('cohens_d_tcga',0)):+.2f}, "
              f"FDR = {format_fdr(bio.get('fdr_tcga',''))}; replicated in "
              f"{int(as_float(bio.get('replication_rate',0)))} of 2 external cohorts.")
    md.append(f"- **PubMed depth**: {r['pubmed_thyroid']} thyroid · "
              f"{r['pubmed_braf']} BRAF-context — classification: {r['novelty_class']}.")
    if cli:
        md.append(f"- **Clinical**: {cli.get('thyroid_trials','0')} thyroid trial(s), "
                  f"{cli.get('any_cancer_trials','0')} any-cancer trial(s); "
                  f"translatability tier = {cli.get('priority','—')}.")
    # Validation modality recommendation
    if g in DRUGGABLE_8 and cli.get("priority")=="high":
        md.append("- **Recommended SNUBH validation modality**: IHC + targeted RT-qPCR on "
                  "SNUBH BRAF V600E PTC archive; cross-reference with available targeted-therapy "
                  "trial inclusion criteria.")
    elif g in DRUGGABLE_8:
        md.append("- **Recommended SNUBH validation modality**: targeted RT-qPCR + IHC pilot "
                  "(n≈30 BRAF V600E vs 30 RAS-like) to confirm subtype-specific over-expression "
                  "before any therapeutic-modality discussion.")
    elif r["novelty_class"]=="novel":
        md.append("- **Recommended SNUBH validation modality**: RT-qPCR pilot only — biomarker "
                  "panel addition candidate, not therapeutic target. Confirm subtype "
                  "discrimination performance vs the BRS panel.")
    else:
        md.append("- **Recommended SNUBH validation modality**: include in expanded biomarker panel; "
                  "no standalone validation needed.")
    md.append("")

(OUT / "snubh_shortlist_rationale.md").write_text("\n".join(md) + "\n")
print(f"[done] {tsv_path}")
print(f"       {OUT}/snubh_shortlist_rationale.md")
print(f"       {len(rows)} candidates ranked")
n_top = sum(1 for r in rows if r["snubh_priority_score"]>=0.65)
n_mid = sum(1 for r in rows if 0.50 <= r["snubh_priority_score"] < 0.65)
n_low = sum(1 for r in rows if r["snubh_priority_score"] < 0.50)
print(f"       tier counts: top={n_top}  mid={n_mid}  low={n_low}")
print(json.dumps({"top": [r["gene"] for r in rows if r["snubh_priority_score"]>=0.65]},
                 indent=2))
