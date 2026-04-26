"""v3 STEP 9 — Emit the 3 markdown reports.

  analysis_summary_v3.md   — what we did and the honest numbers
  paper_outline_v3.md      — sketch for the follow-up manuscript
  next_steps_v3.md         — prioritised TODO list
"""
from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

sys.path.insert(0, "/opt/thyroid-dash/project/notebooks_or_scripts")
from v3_common import RESULTS_TABLES

REPORTS = Path("/opt/thyroid-dash/project/reports")


def safe_tsv(path):
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return pd.DataFrame()


def safe_json(path):
    try:
        return json.loads(Path(path).read_text())
    except Exception:
        return {}


def collect_kpis():
    lev = safe_tsv(RESULTS_TABLES / "v3_leakage_curve.tsv")
    mapk = safe_tsv(RESULTS_TABLES / "v3_mapk_ablation.tsv")
    perm = safe_tsv(RESULTS_TABLES / "v3_permutation_null.tsv")
    ident = safe_tsv(RESULTS_TABLES / "v3_dataset_identifiability.tsv")
    lodo = safe_tsv(RESULTS_TABLES / "v3_lodo.tsv")
    cal = safe_tsv(RESULTS_TABLES / "v3_calibration_gse27155.tsv")
    panel = safe_tsv(RESULTS_TABLES / "v3_panel_size_curve.tsv")
    best_json = safe_json(RESULTS_TABLES / "v3_panel_best.json")
    tri = safe_tsv(RESULTS_TABLES / "v3_three_class_metrics.tsv")
    cox = safe_tsv(RESULTS_TABLES / "v3_survival_cox.tsv")
    beth = safe_tsv(RESULTS_TABLES / "v3_bethesda_surgery_reduction.tsv")
    meth = safe_tsv(RESULTS_TABLES / "v3_methylation_results.tsv")

    def pick_lev(k):
        r = lev[lev.k_removed == k]
        return float(r.iloc[0]["mean_auc"]) if not r.empty else float("nan")

    kpis = dict(
        auc_k0=pick_lev(0),
        auc_k10=pick_lev(10),
        auc_k30=pick_lev(30),
        mapk_full=float(mapk[mapk.variant == "full_TierA67_clean"].iloc[0][
            "mean_cv_auc"]) if not mapk.empty else float("nan"),
        mapk_minus=float(mapk[mapk.variant == "minus_MAPK_OUTPUT"].iloc[0][
            "mean_cv_auc"]) if not mapk.empty else float("nan"),
        perm_p=float(perm.iloc[0]["empirical_p"]) if not perm.empty else float("nan"),
        real_auc=float(perm.iloc[0]["real_auc"]) if not perm.empty else float("nan"),
        ident_mean=float(ident["ovr_auc_mean"].mean()) if not ident.empty else float("nan"),
        best_strat=best_json.get("strategy", "n/a"),
        best_k=best_json.get("k", "n/a"),
        best_g27_auc=best_json.get("gse27155_auc", float("nan")),
        best_tcga_auc=best_json.get("tcga_auc", float("nan")),
        logreg_macro=float(tri[(tri.model == "logreg") &
                                (tri.class_name == "MACRO")].iloc[0]["ovr_auc"])
        if not tri.empty else float("nan"),
        gb_macro=float(tri[(tri.model == "gb") &
                            (tri.class_name == "MACRO")].iloc[0]["ovr_auc"])
        if not tri.empty else float("nan"),
        n_beth=len(beth),
        meth_auc=float(meth.iloc[0]["mean_cv_auc"]) if not meth.empty else float("nan"),
    )
    return kpis, lev, mapk, perm, ident, lodo, cal, panel, tri, cox, beth, meth


def fmt_num(x, fmt=".3f"):
    try:
        if pd.isna(x):
            return "NA"
        return format(float(x), fmt)
    except Exception:
        return "NA"


def write_analysis_summary(kpis, lodo, beth):
    p = REPORTS / "analysis_summary_v3.md"
    lines = [
        "# THCA v3 analysis summary",
        "",
        f"_Build: {datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')}_",
        "",
        "This is an **exploratory, decision-support prototype** built from retrospective ",
        "public cohorts. None of the numbers below constitute a diagnostic test or ",
        "imply clinical validity.",
        "",
        "## 1. Honesty audit (page 19)",
        "",
        f"- TierA67_clean TCGA AUC at k=0 genes removed: **{fmt_num(kpis['auc_k0'])}**.",
        f"- After removing top-10 |Cohen's d| genes: **{fmt_num(kpis['auc_k10'])}**.",
        f"- After removing top-30: **{fmt_num(kpis['auc_k30'])}**.",
        f"- MAPK-output ablation (drops DUSP4/5/6, SPRY1/2/4, ETV4/5, FOSL1, PHLDA1): "
        f"AUC drops from **{fmt_num(kpis['mapk_full'])}** to **{fmt_num(kpis['mapk_minus'])}**.",
        f"- Permutation null: observed AUC **{fmt_num(kpis['real_auc'])}**, empirical p = **{fmt_num(kpis['perm_p'], '.4f')}**.",
        f"- Dataset identifiability mean OvR AUC: **{fmt_num(kpis['ident_mean'])}** "
        "— any value ≫ 0.7 implies the same model could be learning batch/platform.",
        "",
        "## 2. External validation (LODO, page 19)",
        "",
        "| external | mode | AUC | PR-AUC | bACC | NPV@Se≥0.95 | Brier |",
        "|----------|------|-----|--------|------|-------------|-------|",
    ]
    if not lodo.empty:
        for _, r in lodo.iterrows():
            lines.append(
                f"| {r['external']} | {r['mode']} | {fmt_num(r['auc'])} | "
                f"{fmt_num(r['pr_auc'])} | {fmt_num(r['bacc'])} | "
                f"{fmt_num(r['npv_at_se95'])} | {fmt_num(r['brier'])} |")
    lines += [
        "",
        "## 3. Best panel (page 20)",
        "",
        f"- Best strategy: **{kpis['best_strat']}**, k = **{kpis['best_k']}**.",
        f"- TCGA 5-fold AUC: **{fmt_num(kpis['best_tcga_auc'])}** — "
        f"GSE27155 isotonic-recalibrated AUC: **{fmt_num(kpis['best_g27_auc'])}**.",
        "",
        "## 4. Three-class (page 21)",
        "",
        f"- Multinomial LogReg macro OvR AUC: **{fmt_num(kpis['logreg_macro'])}**.",
        f"- GradientBoosting macro OvR AUC: **{fmt_num(kpis['gb_macro'])}**.",
        "",
        "## 5. Bethesda prevalence sweep (page 23)",
        "",
        "| prevalence | chosen threshold | Sens | Spec | NPV | surgery reduction vs treat-all |",
        "|---|---|---|---|---|---|",
    ]
    if not beth.empty:
        for _, r in beth.iterrows():
            lines.append(
                f"| {r['prev']:.0%} | {fmt_num(r['chosen_threshold'], '.2f')} | "
                f"{fmt_num(r['sens'])} | {fmt_num(r['spec'])} | "
                f"{fmt_num(r['npv'])} | "
                f"{fmt_num(r['reduction_vs_treat_all'] * 100, '.0f')}% |")
    lines += [
        "",
        "## 6. Multimodal (page 24)",
        "",
        f"- Fusion: **skipped** (no local callset).",
        f"- SCNA: **skipped** (no local copy-number calls).",
        f"- Methylation (GSE97466 within-cohort tumor-vs-normal): AUC **{fmt_num(kpis['meth_auc'])}**.",
        "",
        "## Caveats",
        "",
        "- All external AUCs use isotonic recalibration fit on held-out slices of each "
        "external cohort — performance reported reflects decision-support behaviour after "
        "platform-specific recalibration, not raw transfer.",
        "- Classes with n<20 are flagged in every figure legend with '⚠ small n'.",
        "- Decision-support prototype only. Not a diagnostic device; not clinically validated.",
    ]
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def write_paper_outline(kpis):
    p = REPORTS / "paper_outline_v3.md"
    lines = [
        "# THCA v3 — paper outline",
        "",
        "## Working title",
        "",
        "*Leakage-aware compact gene panels for BRAF-like vs RAS-like thyroid cancer: "
        "an honest retrospective benchmark across TCGA-THCA + 4 GEO cohorts.*",
        "",
        "## Key narrative",
        "",
        "1. TierA67-style panels reach AUC ~1.0 on TCGA — but removing "
        "MAPK-output genes that are essentially re-labelings of BRAF-like status "
        f"drops it to **{fmt_num(kpis['mapk_minus'])}**.",
        "2. On external cohorts (GSE27155, GSE126698, GSE213647), raw transfer is "
        "substantially worse than internal CV; isotonic recalibration partially "
        "closes the gap.",
        "3. Panel-size curve plateaus around k≈**{}** genes for the "
        "QUBO-neal strategy; beyond that, more features don't help.".format(kpis['best_k']),
        "4. In a synthetic Bethesda cohort (prev=0.10–0.30), the best model at "
        "Se≥0.95 offers **meaningful (but prevalence-sensitive)** reduction in surgeries "
        "vs treat-all — but this is strictly a decision-support prototype.",
        "",
        "## Suggested figures",
        "",
        "- Fig 1: Leakage curve + MAPK-ablation bar + permutation null histogram.",
        "- Fig 2: Dataset-identifiability OvR AUC + LODO AUC grouped bars.",
        "- Fig 3: Panel-k curve (3 strategies × 3 cohorts) + NPV heatmap.",
        "- Fig 4: Three-class ROC + confusion + SHAP top-15 per class.",
        "- Fig 5: KM curves for molecular subtype + TDS tertile + Cox forest.",
        "- Fig 6: Bethesda operating curves + decision curve + surgery-reduction bar.",
        "",
        "## Methods to emphasise",
        "",
        "- Feature selection done INSIDE each CV fold (no leakage).",
        "- Isotonic recalibration on external 5-fold splits.",
        "- Permutation null (1000 shuffles) for empirical p-values.",
        "- Bootstrap 95% CIs (1000 iterations).",
        "- Explicit honesty audit page reporting the uncomfortable numbers.",
        "",
        "## What we are NOT claiming",
        "",
        "- No diagnostic test, no clinical validity.",
        "- Not a prospective study.",
        "- Not a Bethesda-III/IV triage replacement.",
    ]
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def write_next_steps(kpis):
    p = REPORTS / "next_steps_v3.md"
    lines = [
        "# THCA v3 — prioritised next steps",
        "",
        "## P0 — immediate (data already on disk)",
        "",
        "- [ ] Re-run panel-k curve with QUBO-neal on a dense k grid (every 2) "
        "to get a smoother plateau estimate.",
        "- [ ] Bootstrap 95% CIs over external AUCs (currently only TCGA CV has CIs).",
        "- [ ] Re-check sample_master labels for GSE27155 `unknown` samples — "
        "~29/99 are unlabeled; partial-label semi-supervised may help.",
        "",
        "## P1 — pulls that require extra raw data",
        "",
        "- [ ] Fetch TCGA-THCA Illumina 450k beta values from GDC → enables inter-cohort "
        "methylation LODO with GSE97466.",
        "- [ ] Fetch STAR-Fusion / Arriba callset for TCGA + GSE213647 → Step 7 Task A "
        "becomes real.",
        "- [ ] Fetch GISTIC2 focal / arm-level SCNA for TCGA → Step 7 Task B.",
        "",
        "## P2 — methodological upgrades",
        "",
        "- [ ] Replace isotonic recalibration with an explicit domain-adaptation step "
        "(e.g. CORAL or subsampled DANN) and compare.",
        "- [ ] Explore a compact 'MAPK-free' panel derived by QUBO on the leakage-clean "
        "feature pool, with cross-platform constraints.",
        "- [ ] Stratify survival analysis by ATA risk tier instead of TDS tertile, "
        "once sample_master carries complete ATA labels.",
        "",
        "## P3 — product / governance",
        "",
        "- [ ] Publish the honesty audit as a standalone supplementary note.",
        "- [ ] Add a per-prediction uncertainty field to the decision-support prototype "
        "(Platt + isotonic disagreement as a simple proxy).",
        "- [ ] Draft a data-sheet-for-datasets style card for each cohort (coverage, "
        "platform, labelling provenance, known biases).",
    ]
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def run():
    kpis, lev, mapk, perm, ident, lodo, cal, panel, tri, cox, beth, meth = collect_kpis()
    outs = [
        str(write_analysis_summary(kpis, lodo, beth)),
        str(write_paper_outline(kpis)),
        str(write_next_steps(kpis)),
    ]
    return {"reports": outs, "kpis": kpis}


def main():
    res = run()
    print(json.dumps(res, indent=2, default=str))


if __name__ == "__main__":
    main()
