#!/usr/bin/env python3
"""AUDIT 12 — functional decomposition of the 8-gene panel into 3 sub-modules.

Extends (does not redo) `tcga_rai_structural_disease_2026_08_06.py` /
`scripts/audit/03_rare_event_structural_outcomes.py`, which established (Firth
penalised logistic regression, TCGA-THCA, RAI-treated patients only):

  new tumour event after initial treatment (9 events / 145)   OR = 0.221 [0.062,0.784] P=0.019
  persistent disease within 3 months of surgery (13 / 72)     OR = 0.224 [0.061,0.826] P=0.025

for the aggregate RAI_8 score (mean cohort z-score across all 8 panel genes).

This script asks a genuinely new question: does the aggregate association come from all
three functional layers of the panel equally, or is it driven by one? The 8 genes are
decomposed into 3 biologically distinct sub-modules:

  LINEAGE_4  = mean cohort-z(PAX8, NKX2-1, FOXE1, TSHR)     thyrocyte identity / master TFs + TSH receptor
  UPTAKE_1   =      cohort-z(SLC5A5)                        iodide trapping (NIS) — the single gene the
                                                             radioiodine mechanism most directly depends on
  HORMONE_3  = mean cohort-z(TG, TPO, DIO1)                 hormone synthesis/organification machinery

Same cohort, same endpoints, same rare-event corrections (Firth score-only, Firth+stage,
Firth+purity, Firth+stage+purity exploratory, leave-one-event-out, patient bootstrap) as the
originals, so the sub-module results are directly comparable to the aggregate numbers above.

Cohort-z base: the same 513-sample TCGA-PTC GDC study cohort that RAI_8 itself was
z-scored against (`panel_expression_thpa_tcga_gdc.tsv`), NOT the 145/72-patient
RAI-treated subset — matching the original scoring convention exactly
(`project/results/ncomm_push_2026_05_08/_common.py::score_panel`).

Multiplicity: 3 sub-modules x 2 endpoints = 6 primary (score-only) tests where the original
analysis ran 1. Bonferroni alpha = 0.05 / 6 = 0.00833 is applied and reported alongside the
raw P-values; nothing here is presented as an independent discovery at the nominal 0.05 level.

Outputs (new files only):
  rai-response-genomics-atlas/results/tables/audit12_module_scores_2026_08_19.tsv
  rai-response-genomics-atlas/results/tables/audit12_module_firth_main_2026_08_19.tsv
  rai-response-genomics-atlas/results/tables/audit12_module_stability_2026_08_19.tsv
  audit/rai_integration_20260819/B_functional_decomposition_tcga.md
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]          # rai-response-genomics-atlas/
REPO_ROOT = ROOT.parent                              # THCA_data_analysis/
TAB = ROOT / "results" / "tables"
AUDIT_OUT = REPO_ROOT / "audit" / "rai_integration_20260819"
TAB.mkdir(parents=True, exist_ok=True)
AUDIT_OUT.mkdir(parents=True, exist_ok=True)

STAMP = "2026_08_19"
SEED = 20260806  # reuse the original seed for exact bootstrap comparability

PANEL_FULL = (REPO_ROOT / "project" / "results" / "ncomm_push_2026_05_08" / "cbioportal_sweep"
              / "panel_expression_thpa_tcga_gdc.tsv")

# ---- import the original modules directly so every helper (cohens_d, boot_ci_d, auc_mw,
#      firth_logit, the file paths CBIO/PURITY/SURV/RADIO) is byte-identical to the audit ----
_spec_base = importlib.util.spec_from_file_location(
    "rai_base", ROOT / "scripts" / "tcga_rai_best_response_2026_08_06.py")
_base = importlib.util.module_from_spec(_spec_base)
_spec_base.loader.exec_module(_base)

_spec_firth = importlib.util.spec_from_file_location(
    "rai_firth", Path(__file__).resolve().parent / "03_rare_event_structural_outcomes.py")
# executing that module runs its own main() only under __main__ guard, so import is side-effect free
_firth_mod = importlib.util.module_from_spec(_spec_firth)
_spec_firth.loader.exec_module(_firth_mod)

firth_logit = _firth_mod.firth_logit
cohens_d = _base.cohens_d
boot_ci_d = _base.boot_ci_d
auc_mw = _base.auc_mw
CBIO = _base.CBIO
PURITY = _base.PURITY
SURV = _base.SURV
RADIO = _base.RADIO

MODULES = {
    "LINEAGE_4": ["PAX8", "NKX2-1", "FOXE1", "TSHR"],
    "UPTAKE_1": ["SLC5A5"],
    "HORMONE_3": ["TG", "TPO", "DIO1"],
}
ALL8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]


def build_module_scores() -> pd.DataFrame:
    """Cohort-z each of the 8 genes over the SAME 513-sample base cohort RAI_8 used, then
    average within each sub-module. Reproduces `_common.py::score_panel`'s exact formula
    (per-gene cohort z = (x - cohort_mean) / cohort_std, then mean across genes) but applied
    to gene subsets instead of all 8."""
    raw = pd.read_csv(PANEL_FULL, sep="\t")
    missing = [g for g in ALL8 if g not in raw.columns]
    if missing:
        raise RuntimeError(f"panel_expression file missing genes: {missing}")
    sub = raw[["sampleId"] + ALL8].copy()
    z = (sub[ALL8] - sub[ALL8].mean(axis=0)) / sub[ALL8].std(axis=0).replace(0, 1)

    out = pd.DataFrame({"sampleId": sub["sampleId"]})
    for mod_name, genes in MODULES.items():
        out[mod_name] = z[genes].mean(axis=1)
    out["RAI_8_recomputed"] = z[ALL8].mean(axis=1)  # sanity check against the stored RAI_8

    # sanity check: recomputed aggregate must match the stored RAI_8 to numerical precision
    check = raw[["sampleId", "RAI_8"]].merge(out[["sampleId", "RAI_8_recomputed"]], on="sampleId")
    max_abs_diff = float((check["RAI_8"] - check["RAI_8_recomputed"]).abs().max())
    print(f"sanity check: max|RAI_8 - recomputed| across {len(check)} samples = {max_abs_diff:.2e}")
    if max_abs_diff > 1e-6:
        raise RuntimeError(
            f"module-score recomputation does not reproduce stored RAI_8 (max diff {max_abs_diff}); "
            "cohort-z base or gene set does not match the original scoring convention — STOP.")

    # weighted-mean identity check: RAI_8 = (4*LINEAGE_4 + 1*UPTAKE_1 + 3*HORMONE_3) / 8
    recon = (4 * out["LINEAGE_4"] + 1 * out["UPTAKE_1"] + 3 * out["HORMONE_3"]) / 8
    max_recon_diff = float((out["RAI_8_recomputed"] - recon).abs().max())
    print(f"module weighted-sum identity check: max diff = {max_recon_diff:.2e}")

    out["patient"] = out["sampleId"].str.slice(0, 12)
    out = out.drop_duplicates("patient")
    return out


def build_cohort(mod_scores: pd.DataFrame) -> pd.DataFrame:
    """Reproduce `tcga_rai_best_response_2026_08_06.py::build()` exactly, substituting the
    3 sub-module scores (+ recomputed RAI_8 as a sanity anchor) for the panel file it reads."""
    rad = pd.read_csv(RADIO, sep="\t", skiprows=[1, 2], low_memory=False)
    rad["dose"] = pd.to_numeric(rad["radiation_total_dose"], errors="coerce")
    rai = rad[rad["radiation_adjuvant_units"].astype(str).str.lower().eq("mci")].copy()

    agg = rai.groupby("bcr_patient_barcode").agg(
        rai_courses=("bcr_radiation_barcode", "count"),
        rai_dose_mci=("dose", "sum"),
    ).reset_index().rename(columns={"bcr_patient_barcode": "patient"})

    df = mod_scores.merge(agg, on="patient", how="inner")

    pur = pd.read_csv(PURITY, sep="\t")
    pcol = next((c for c in pur.columns if "purity" in c.lower()), None) or "leuko_frac"
    pur["patient"] = pur["sample_id"].astype(str).str.slice(0, 12)
    pur[pcol] = pd.to_numeric(pur[pcol], errors="coerce")
    pur = pur.dropna(subset=[pcol]).drop_duplicates("patient")
    df = df.merge(pur[["patient", pcol]].rename(columns={pcol: "purity"}), on="patient", how="left")

    cb = pd.read_csv(CBIO, sep="\t", low_memory=False)
    keep = ["patientId", "AGE", "AJCC_PATHOLOGIC_TUMOR_STAGE", "TUMOR_STATUS",
            "CLINICAL_STATUS_WITHIN_3_MTHS_SURGERY", "NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT"]
    cb = cb[[c for c in keep if c in cb.columns]].rename(columns={"patientId": "patient"})
    df = df.merge(cb, on="patient", how="left")
    return df


def fit_report(df: pd.DataFrame, ycol: str, score_col: str, terms_extra: list[str], label: str):
    terms = [score_col] + terms_extra
    d = df[[ycol] + terms].apply(pd.to_numeric, errors="coerce").dropna()
    if d[ycol].sum() < 2 or (1 - d[ycol]).sum() < 2:
        return None
    X = np.column_stack([np.ones(len(d))] + [d[t].values for t in terms])
    y = d[ycol].values
    b, se, pv = firth_logit(X, y)
    i = 1  # first term after intercept is the sub-module score
    return dict(model=label, n=len(d), events=int(y.sum()),
                epv=round(y.sum() / len(terms), 2),
                or_=float(np.exp(b[i])),
                ci_low=float(np.exp(b[i] - 1.96 * se[i])),
                ci_high=float(np.exp(b[i] + 1.96 * se[i])),
                p=float(pv[i]), method="Firth")


def main():
    mod_scores = build_module_scores()
    mod_scores.to_csv(TAB / f"audit12_module_scores_{STAMP}.tsv", sep="\t", index=False)

    df = build_cohort(mod_scores)

    nte = df["NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT"].astype(str).str.upper()
    cs = df["CLINICAL_STATUS_WITHIN_3_MTHS_SURGERY"].astype(str).str.lower()
    df["stage_high"] = df["AJCC_PATHOLOGIC_TUMOR_STAGE"].astype(str).str.contains(
        "III|IV", regex=True).astype(int)
    df["y_new"] = np.where(nte.eq("YES"), 1, np.where(nte.eq("NO"), 0, np.nan))
    df["y_pers"] = np.where(cs.str.contains("persistent"), 1,
                            np.where(cs.str.contains("no imaging evidence|no evidence"), 0, np.nan))

    print(f"cohort size after merge (matches original build(), RAI-treated with panel score): {len(df)}")

    endpoints = [("y_new", "new tumour event after initial treatment"),
                 ("y_pers", "persistent disease within 3 months of surgery")]
    score_cols = ["RAI_8_recomputed", "LINEAGE_4", "UPTAKE_1", "HORMONE_3"]

    rows, stab = [], []
    rng = np.random.default_rng(SEED)

    for ycol, nice in endpoints:
        sub = df[df[ycol].notna()].copy()
        for score_col in score_cols:
            a = sub.loc[sub[ycol] == 1, score_col].dropna().values
            b = sub.loc[sub[ycol] == 0, score_col].dropna().values
            if len(a) >= 2 and len(b) >= 2:
                d_val = cohens_d(a, b)
                lo, hi = boot_ci_d(a, b)
                rows.append(dict(endpoint=nice, score=score_col,
                                 model="Tier 1 . Cohen's d (patient bootstrap)",
                                 n=len(a) + len(b), events=len(a), epv=np.nan,
                                 or_=float(d_val), ci_low=float(lo), ci_high=float(hi),
                                 p=float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue),
                                 method="nonparametric"))

            for terms_extra, lab in [([], "Tier 2 . Firth, score only"),
                                     (["stage_high"], "Tier 3 . Firth + stage"),
                                     (["purity"], "Tier 3 . Firth + purity"),
                                     (["stage_high", "purity"],
                                      "Tier 3 . Firth + stage + purity (exploratory)")]:
                r = fit_report(sub, ycol, score_col, terms_extra, lab)
                if r:
                    r["endpoint"] = nice
                    r["score"] = score_col
                    rows.append(r)
                    print(f"{nice[:38]:38s} {score_col:16s} {lab:38s} OR={r['or_']:.3f} "
                          f"[{r['ci_low']:.3f},{r['ci_high']:.3f}] P={r['p']:.4f} EPV={r['epv']}")

            # stability: LOO + bootstrap, univariable Firth model only (matches audit03 scope)
            m = sub[[ycol, score_col]].apply(pd.to_numeric, errors="coerce").dropna()
            ev_idx = m.index[m[ycol] == 1].tolist()
            loo = []
            for e in ev_idx:
                mm = m.drop(index=e)
                X = np.column_stack([np.ones(len(mm)), mm[score_col].values])
                bb, _, _ = firth_logit(X, mm[ycol].values)
                loo.append(np.exp(bb[1]))
            boot = []
            for _ in range(2000):
                idx = rng.choice(m.index.values, len(m), replace=True)
                mm = m.loc[idx]
                if mm[ycol].sum() < 2 or (1 - mm[ycol]).sum() < 2:
                    continue
                X = np.column_stack([np.ones(len(mm)), mm[score_col].values])
                bb, _, _ = firth_logit(X, mm[ycol].values)
                boot.append(np.exp(bb[1]))
            boot = np.asarray(boot)
            if len(loo) and len(boot):
                stab.append(dict(endpoint=nice, score=score_col, n_events=len(ev_idx),
                                 loo_or_min=round(float(np.min(loo)), 3),
                                 loo_or_max=round(float(np.max(loo)), 3),
                                 loo_all_below_1=bool(np.all(np.asarray(loo) < 1)),
                                 boot_median_or=round(float(np.median(boot)), 3),
                                 boot_ci_low=round(float(np.percentile(boot, 2.5)), 3),
                                 boot_ci_high=round(float(np.percentile(boot, 97.5)), 3),
                                 boot_frac_or_below_1=round(float((boot < 1).mean()), 3),
                                 n_boot_used=len(boot)))
                print(f"   stability [{score_col}]: LOO OR {np.min(loo):.3f}-{np.max(loo):.3f}; "
                      f"boot median {np.median(boot):.3f}, {100*(boot<1).mean():.1f}% below 1")

    res = pd.DataFrame(rows)

    # Benjamini-Hochberg AND Bonferroni across the 6 primary (score-only, univariable) tests:
    # 3 sub-modules x 2 endpoints. This is the multiplicity family that matters because it is
    # the family actually being newly tested here (the original tested 1 aggregate score).
    primary_mask = res["model"].eq("Tier 2 . Firth, score only") & res["score"].isin(
        ["LINEAGE_4", "UPTAKE_1", "HORMONE_3"])
    n_primary = int(primary_mask.sum())
    res["bonferroni_alpha"] = np.nan
    res["p_bonferroni_sig"] = np.nan
    res.loc[primary_mask, "bonferroni_alpha"] = 0.05 / n_primary
    res.loc[primary_mask, "p_bonferroni_sig"] = res.loc[primary_mask, "p"] < (0.05 / n_primary)
    ranked = res.loc[primary_mask, "p"].rank(method="first")
    res.loc[primary_mask, "q_bh"] = (res.loc[primary_mask, "p"] * n_primary / ranked).clip(upper=1.0)

    res.to_csv(TAB / f"audit12_module_firth_main_{STAMP}.tsv", sep="\t", index=False)
    st = pd.DataFrame(stab)
    st.to_csv(TAB / f"audit12_module_stability_{STAMP}.tsv", sep="\t", index=False)

    print(f"\nn_primary tests (score-only, 3 modules x 2 endpoints) = {n_primary}, "
          f"Bonferroni alpha = {0.05 / n_primary:.5f}")
    print(res.loc[primary_mask, ["endpoint", "score", "n", "events", "or_", "ci_low", "ci_high",
                                  "p", "q_bh", "p_bonferroni_sig"]].to_string(index=False))

    return res, st, df


if __name__ == "__main__":
    main()
