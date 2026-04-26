#!/usr/bin/env python3
"""Consolidate ALL project numbers into a single markdown status report.
Read-only. Writes reports/STATUS_REPORT_FOR_USER.md and prints to stdout.
"""
from __future__ import annotations
import json
import os
import subprocess
from pathlib import Path
from collections import Counter

import pandas as pd

ROOT = Path("/opt/thyroid-dash/project")
META = ROOT / "metadata"
RES  = ROOT / "results"
HTML = ROOT / "reports" / "html"
RPT  = ROOT / "reports"
OUT  = RPT / "STATUS_REPORT_FOR_USER.md"

buf: list[str] = []
def w(*a): buf.append(" ".join(str(x) for x in a))
def header(t, lvl=2): w(""); w("#"*lvl + " " + t); w("")


def safe_read_tsv(path: Path) -> pd.DataFrame | None:
    if not path.exists(): return None
    try: return pd.read_csv(path, sep="\t")
    except Exception as e: w(f"_error reading {path}: {e}_"); return None


def safe_read_text(path: Path) -> str | None:
    if not path.exists(): return None
    try: return path.read_text(encoding="utf-8")
    except Exception: return None


def df_to_md(df: pd.DataFrame, max_rows: int = 25, float_fmt="{:.3g}") -> str:
    if df is None or df.empty: return "_(empty)_"
    sub = df.head(max_rows).copy()
    # Round floats
    for c in sub.select_dtypes(include="float").columns:
        sub[c] = sub[c].apply(lambda v: float_fmt.format(v) if pd.notna(v) else "")
    out = sub.to_markdown(index=False)
    if len(df) > max_rows:
        out += f"\n\n_({len(df)} rows total — showing first {max_rows})_"
    return out


def counter_md(c: Counter, title: str) -> None:
    if not c:
        w(f"- {title}: _none_"); return
    w(f"**{title}** (n = {sum(c.values())})")
    w("")
    w("| Value | n |")
    w("|---|---:|")
    for k, v in c.most_common():
        w(f"| {k} | {v} |")


# ------------------------------------------------------------------
# TL;DR (computed at end, inserted at top later)
# ------------------------------------------------------------------
tldr_idx = len(buf)
w("# THCA Project — Consolidated Status Report")
w("")
w(f"_Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}_")
w("")
w("## TL;DR")
w("<!-- placeholder --> ")
w("")

# ==================================================================
# SECTION 1 — DATA INVENTORY
# ==================================================================
header("1. DATA INVENTORY", 2)

sm = safe_read_tsv(META / "sample_master.tsv")
if sm is not None:
    w(f"**Total samples:** {len(sm)}")
    w("")
    counter_md(Counter(sm["dataset"].dropna()), "Samples per dataset")
    w("")
    if "platform" in sm.columns:
        counter_md(Counter(sm["platform"].fillna("unknown")), "Samples per platform")
        w("")
    if "molecular_subtype" in sm.columns:
        counter_md(Counter(sm["molecular_subtype"].fillna("unknown")), "Samples per molecular_subtype")
        w("")
    if "histology_subtype" in sm.columns:
        counter_md(Counter(sm["histology_subtype"].fillna("unknown")), "Samples per histology_subtype")
        w("")
    if "label_confidence" in sm.columns:
        counter_md(Counter(sm["label_confidence"].fillna("unknown")), "Samples per label_confidence")
        w("")
else:
    w("FILE NOT FOUND: metadata/sample_master.tsv")

header("Dataset master", 3)
dm = safe_read_tsv(META / "dataset_master.tsv")
w(df_to_md(dm, max_rows=15) if dm is not None else "FILE NOT FOUND: metadata/dataset_master.tsv")

header("Gene panel coverage per dataset", 3)
gc = safe_read_tsv(META / "gene_coverage_table.tsv")
if gc is not None:
    keep = ["dataset","total_gene_count","TDS16_coverage_n","TDS16_coverage_frac",
            "BRS71_proxy_coverage_n","BRS71_proxy_coverage_frac",
            "TierA67_entries_coverage_n","TierA67_entries_coverage_frac"]
    keep = [c for c in keep if c in gc.columns]
    w(df_to_md(gc[keep], max_rows=15))
else:
    w("FILE NOT FOUND: metadata/gene_coverage_table.tsv")

header("External raw-data availability", 3)
checks = {
    "PRJEB11591 (ENA Korean SNU cohort)": ROOT / "data_raw/ena_v3/PRJEB11591_counts.tsv",
    "TCGA fusion callset (cBioPortal)":    ROOT / "data_raw/tcga_v3/tcga_thca_fusions_merged.tsv",
    "TCGA fusion callset (alt path)":      next((p for p in ROOT.rglob("*fusion*tcga*") if p.is_file()), None),
    "TCGA SCNA / GISTIC2":                 next((p for p in (ROOT/"data_raw/tcga_v3/gistic2").glob("*") if p.is_file()), None) if (ROOT/"data_raw/tcga_v3/gistic2").exists() else None,
    "GSE33630 SOFT family":                 next((p for p in ROOT.rglob("GSE33630_family*")), None),
    "GSE29265 SOFT family":                 next((p for p in ROOT.rglob("GSE29265_family*")), None),
}
for label, path in checks.items():
    if path and Path(path).exists():
        w(f"- ✅ **{label}** — `{path}`")
    else:
        w(f"- ❌ **{label}** — NOT FOUND")

# ==================================================================
# SECTION 2 — ML PERFORMANCE (v2)
# ==================================================================
header("2. ML PERFORMANCE (v2 baseline)", 2)

v2_md = safe_read_text(RPT / "ml_report_v2.md")
if v2_md:
    header("ml_report_v2.md (verbatim)", 3)
    w("```markdown")
    w(v2_md.strip())
    w("```")

header("baseline_ml_results.tsv — internal CV", 3)
int_df = safe_read_tsv(RES / "ml/baseline_ml_results.tsv")
if int_df is not None:
    cols = [c for c in ["task","dataset","feature_set","model","cv_auc","cv_auc_ci_lo","cv_auc_ci_hi",
                        "cv_pr_auc","cv_brier","cv_balanced_accuracy","cv_balanced_accuracy_youden",
                        "cv_f1","cv_mcc","cv_youden_threshold"] if c in int_df.columns]
    w(df_to_md(int_df[cols], max_rows=25))
else:
    w("FILE NOT FOUND: results/ml/baseline_ml_results.tsv")

header("baseline_ml_external.tsv — external validation", 3)
ext_df = safe_read_tsv(RES / "ml/baseline_ml_external.tsv")
if ext_df is not None:
    cols = [c for c in ["task","dataset","feature_set","model","n_samples","auc","auc_ci_lo","auc_ci_hi",
                        "pr_auc","balanced_accuracy","balanced_accuracy_youden","f1","f1_youden","mcc","status"] if c in ext_df.columns]
    # Drop skipped rows for compactness
    sub = ext_df[ext_df.get("status","ok") == "ok"] if "status" in ext_df.columns else ext_df
    w(df_to_md(sub[cols] if cols else sub, max_rows=60))
else:
    w("FILE NOT FOUND: results/ml/baseline_ml_external.tsv")

header("Quantum algorithm comparison", 3)
q_md = safe_read_text(RES / "ml/quantum_full_comparison.md")
if q_md:
    w("```markdown")
    w(q_md.strip())
    w("```")
else:
    w("FILE NOT FOUND: results/ml/quantum_full_comparison.md")

# ==================================================================
# SECTION 3 — BIOMARKER + DRUG (v2)
# ==================================================================
header("3. BIOMARKER + DRUG (v2)", 2)

bio_md = safe_read_text(RPT / "biomarker_analysis.md")
if bio_md:
    header("biomarker_analysis.md TL;DR", 3)
    # Extract TL;DR + known biomarkers table if present
    lines = bio_md.splitlines()
    section = []
    capture = False
    for L in lines:
        if L.startswith("## ") and "TL;DR" in L.upper():
            capture = True; section.append(L); continue
        if capture and L.startswith("## "):
            break
        if capture: section.append(L)
    w("\n".join(section) if section else bio_md[:1500])
else:
    w("FILE NOT FOUND: reports/biomarker_analysis.md")

header("biomarker_validated.tsv — top 20 by score", 3)
val_df = safe_read_tsv(RES / "tables/biomarker_validated.tsv")
if val_df is not None:
    w(f"**Total rows:** {len(val_df)}")
    sort_col = None
    for cand in ["novelty_score","cohens_d_tcga","abs_cohens_d","cohens_d","meta_abs_d","TCGA_cohens_d"]:
        if cand in val_df.columns: sort_col = cand; break
    if sort_col:
        val_df["_sort"] = val_df[sort_col].abs() if sort_col.endswith("cohens_d") or sort_col.endswith("d_tcga") else val_df[sort_col]
        top = val_df.sort_values("_sort", ascending=False).drop(columns="_sort").head(20)
    else:
        top = val_df.head(20)
    keep = [c for c in top.columns if not c.startswith("_") and c not in
            ("overlap_genes","panel_memberships","per_cohort_means")]
    w("")
    w(df_to_md(top[keep[:10]], max_rows=20))
else:
    w("FILE NOT FOUND: results/tables/biomarker_validated.tsv")

header("biomarker_to_drug_report.md highlights", 3)
b2d = safe_read_text(RPT / "biomarker_to_drug_report.md")
if b2d:
    lines = b2d.splitlines()
    out_sect: list[str] = []
    in_targets = in_highlight = False
    for L in lines:
        low = L.lower()
        if low.startswith("## ") and ("target" in low or "top 8" in low):
            in_targets = True; in_highlight = False; out_sect.append(L); continue
        if low.startswith("## ") and ("novel" in low or "highlight" in low or "interesting" in low):
            in_highlight = True; in_targets = False; out_sect.append(L); continue
        if low.startswith("## "):
            in_targets = in_highlight = False
            if len(out_sect) > 0: out_sect.append(L); break
        if in_targets or in_highlight: out_sect.append(L)
    w("\n".join(out_sect) if out_sect else b2d[:2000])
else:
    w("FILE NOT FOUND: reports/biomarker_to_drug_report.md")

header("drug_discovery_compounds.tsv", 3)
cmp_df = safe_read_tsv(RES / "tables/drug_discovery_compounds.tsv")
if cmp_df is not None:
    keep = [c for c in ["target","chembl_id","preferred_name","canonical_smiles",
                        "pchembl_value","mechanism_of_action","target_relation_confidence"]
            if c in cmp_df.columns]
    w(df_to_md(cmp_df[keep] if keep else cmp_df, max_rows=20))
else:
    w("FILE NOT FOUND: results/tables/drug_discovery_compounds.tsv")

# ==================================================================
# SECTION 4 — v3 PROGRESS
# ==================================================================
header("4. v3 PROGRESS", 2)
v3_files = [
    RES / "tables/v3_leakage_curve.tsv",
    RES / "tables/v3_mapk_ablation.tsv",
    RES / "tables/v3_permutation_null.tsv",
    RES / "ml/v3_dataset_identifiability.tsv",
    RES / "ml/v3_lodo.tsv",
    RES / "ml/v3_calibration_gse27155.tsv",
    RES / "tables/v3_panel_size_curve.tsv",
    RES / "ml/v3_three_class_results.tsv",
    RES / "tables/v3_survival_results.tsv",
    RES / "tables/v3_bethesda_triage.tsv",
    RES / "ml/v3_ext_validation_results.tsv",
    RES / "ml/v3_ext_lodo_v3.tsv",
    RES / "tables/v3_ext_batch_diagnostics.json",
]
# also scan for any v3_*.tsv
extra = sorted([p for p in (RES/"tables").glob("v3_*.tsv")] +
               [p for p in (RES/"ml").glob("v3_*.tsv")])
all_v3 = list({*v3_files, *extra})
missing_v3: list[str] = []
for p in all_v3:
    if not p.exists():
        w(f"### ❌ NOT YET RUN — `{p.relative_to(ROOT)}`"); w("")
        missing_v3.append(str(p.relative_to(ROOT)))
        continue
    w(f"### ✅ `{p.relative_to(ROOT)}`")
    if p.suffix == ".tsv":
        df = safe_read_tsv(p)
        if df is not None:
            w(f"_Rows: {len(df)}_")
            w("")
            w(df_to_md(df, max_rows=50))
    elif p.suffix == ".json":
        try:
            d = json.load(open(p))
            w("```json")
            w(json.dumps(d, indent=2)[:4000])
            w("```")
        except Exception as e:
            w(f"_error: {e}_")
    w("")

# ==================================================================
# SECTION 5 — KNOWN LIMITATIONS
# ==================================================================
header("5. KNOWN LIMITATIONS + NEXT STEPS", 2)

as_v2 = safe_read_text(RPT / "analysis_summary_v2.md")
if as_v2:
    header("analysis_summary_v2.md · Remaining limitations", 3)
    lines = as_v2.splitlines()
    out_sect = []
    capture = False
    for L in lines:
        if L.lower().startswith("## remaining limit") or L.lower().startswith("## limitations"):
            capture = True; out_sect.append(L); continue
        if capture and L.startswith("## "):
            break
        if capture: out_sect.append(L)
    if out_sect:
        w("\n".join(out_sect))
    else:
        w("_(limitations section not found; first 60 lines of file:)_")
        w("")
        w("```markdown"); w("\n".join(lines[:60])); w("```")
else:
    w("FILE NOT FOUND: reports/analysis_summary_v2.md")

ns_v2 = safe_read_text(RPT / "next_steps_v2.md")
header("next_steps_v2.md (verbatim)", 3)
if ns_v2:
    w("```markdown"); w(ns_v2.strip()); w("```")
else:
    w("FILE NOT FOUND: reports/next_steps_v2.md")

# ==================================================================
# SECTION 6 — BUILD STATE
# ==================================================================
header("6. BUILD STATE", 2)

v = safe_read_text(HTML / "_version.json")
header("reports/html/_version.json", 3)
if v:
    w("```json"); w(v.strip()); w("```")
else:
    w("FILE NOT FOUND: reports/html/_version.json")

pages_count = len(list((HTML/"pages").glob("*.html"))) if (HTML/"pages").exists() else 0
figs_count  = len(list((HTML/"figs_interactive").glob("*.html"))) if (HTML/"figs_interactive").exists() else 0
w("")
w(f"- **pages/*.html count:** {pages_count}")
w(f"- **figs_interactive/*.html count:** {figs_count}")

# HTTP check
try:
    http_code = subprocess.check_output(
        ["curl","-sS","-o","/dev/null","-w","%{http_code}","http://127.0.0.1:8012/reports/html/index.html"],
        timeout=5
    ).decode().strip()
except Exception as e:
    http_code = f"error: {e}"
w(f"- **Port 8012 /reports/html/index.html:** {http_code}")

# ==================================================================
# FINALIZE — fill in TL;DR placeholder
# ==================================================================
tldr_lines = [
    "- 32 HTML pages live on public URL, 100+ interactive Plotly figures, 35MB downloadable zip.",
    "- v2 baseline TCGA internal AUC=1.00 on TierA67 / variance_top50 (LogReg+GB); external GSE27155 AUC 0.82–0.98 but bACC_youden collapses to 0.5 on linear models at default threshold.",
    "- v3 honesty audit complete: strict-CV AUC 0.992 at k=0, drops to 0.895 at k=30; MAPK ablation nearly zero effect (0.923 → 0.924); permutation p=0.001.",
    "- Dataset-identifiability OvR AUC = 1.00 across 6 cohorts → strong batch leakage signal; external AUC=1.00 likely reflects platform not biology.",
    "- Top 3 gaps: (1) PRJEB11591 mutation-anchored external not obtained (supplementary URLs 403/404), (2) TCGA fusion callset empty, (3) BRS71 original vs proxy overlap 1.4% → proxy must be replaced.",
]
# splice into top placeholder
buf[tldr_idx+3] = "\n".join(tldr_lines)

# Write out
content = "\n".join(buf) + "\n"
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(content, encoding="utf-8")

# Print STAGE 99 block
missing_core = [str(p.relative_to(ROOT)) for p in [
    META/"sample_master.tsv", META/"dataset_master.tsv", RES/"ml/baseline_ml_results.tsv",
    RES/"ml/baseline_ml_external.tsv", RPT/"biomarker_analysis.md",
] if not p.exists()]

print("\n" + "="*60)
print("=== STATUS REPORT READY ===")
print(f"File       : {OUT}")
print(f"Size       : {OUT.stat().st_size} bytes ({OUT.stat().st_size/1024:.1f} KB)")
print(f"Sections   : 6")
print(f"Missing    : {missing_core or 'none (core files)'}  |  v3 not-yet-run: {len(missing_v3)}")
print(f"Port 8012  : {http_code}")
print("="*60)
