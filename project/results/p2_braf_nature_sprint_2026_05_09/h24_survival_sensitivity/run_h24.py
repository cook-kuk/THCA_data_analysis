"""
H24 — Sensitivity battery for the H6 finding (DM2-vs-not_DM PFI HR=5.91 in BRAF-cPTC).

Adds external clinical fields not in h6_merged_clinical.tsv:
  - cBioPortal thca_tcga: EXTRATHYROIDAL_EXTENSION, AJCC_PATHOLOGIC_TUMOR_STAGE,
    AJCC_TUMOR_PATHOLOGIC_PT, AJCC_NODES_PATHOLOGIC_PN, AJCC_METASTASIS_PATHOLOGIC_PM,
    INITIAL_PATHOLOGIC_DX_YEAR, TISSUE_SOURCE_SITE, NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT
  - GDC: tumor_focality (multifocality)

Sensitivity slices on the H6 BRAF-cPTC stratum:
  1. Stage-stratified Cox (strata=tumor_stage)        - per-stage Cox
  2. T-stage stratified (T1 vs T2-4) + interaction
  3. ETE yes/no, N0 vs N1, M0 vs M1, multifocality - subgroup HR + interaction p
  4. Era split (<=2010 vs >2010) - subgroup HR + interaction p
  5. DM1 vs DM2 same battery (protective sensitivity)
  6. Robust SE + cluster by TISSUE_SOURCE_SITE
  7. Forest plot data table

Output dir:
  project/results/p2_braf_nature_sprint_2026_05_09/h24_survival_sensitivity/
"""

from __future__ import annotations
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import urllib.request
import urllib.parse

from lifelines import CoxPHFitter, KaplanMeierFitter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
np.random.seed(0)

OUT = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/p2_braf_nature_sprint_2026_05_09/h24_survival_sensitivity"
)
OUT.mkdir(parents=True, exist_ok=True)
H6 = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/p2_braf_nature_sprint_2026_05_09/h6_survival/h6_merged_clinical.tsv"
)
CACHE = OUT / "cache"
CACHE.mkdir(parents=True, exist_ok=True)


# ----------------------------------------------------------------------
# 1. Load H6 merged clinical (DM, TERT, OS/PFI/DFI/DSS, age, stage, sex)
# ----------------------------------------------------------------------
def load_h6():
    df = pd.read_csv(H6, sep="\t")
    df["patient12"] = df["sample_id"].str[:12]
    print(f"[H6] loaded {len(df)} samples; DM dist: {df['dm'].value_counts().to_dict()}")
    return df


# ----------------------------------------------------------------------
# 2. Pull cBioPortal thca_tcga patient clinical (cached)
# ----------------------------------------------------------------------
def fetch_cbio_thca_tcga():
    cache = CACHE / "cbio_thca_tcga_patient_clinical.tsv"
    if cache.exists():
        print(f"[cBio] using cached {cache}")
        return pd.read_csv(cache, sep="\t", dtype=str)

    url = "https://www.cbioportal.org/api/studies/thca_tcga/clinical-data?clinicalDataType=PATIENT"
    print(f"[cBio] fetching {url}")
    with urllib.request.urlopen(url, timeout=120) as resp:
        data = json.loads(resp.read())
    pats: dict[str, dict] = {}
    for r in data:
        pid = r["patientId"]
        pats.setdefault(pid, {"patientId": pid})
        pats[pid][r["clinicalAttributeId"]] = r["value"]
    df = pd.DataFrame(pats.values())
    df.to_csv(cache, sep="\t", index=False)
    print(f"[cBio] cached {len(df)} patients @ {cache}")
    return df


# ----------------------------------------------------------------------
# 3. Pull GDC tumor_focality + year_of_diagnosis (cached)
# ----------------------------------------------------------------------
def fetch_gdc_thca():
    cache = CACHE / "gdc_thca_focality.tsv"
    if cache.exists():
        print(f"[GDC] using cached {cache}")
        return pd.read_csv(cache, sep="\t", dtype=str)

    fields = [
        "submitter_id",
        "diagnoses.tumor_focality",
        "diagnoses.year_of_diagnosis",
        "diagnoses.ajcc_pathologic_t",
        "diagnoses.ajcc_pathologic_n",
        "diagnoses.ajcc_pathologic_m",
        "diagnoses.ajcc_pathologic_stage",
        "tissue_source_site.bcr_id",
    ]
    filters = {
        "op": "and",
        "content": [{
            "op": "in",
            "content": {"field": "project.project_id", "value": ["TCGA-THCA"]},
        }],
    }
    url = (
        "https://api.gdc.cancer.gov/cases"
        f"?filters={urllib.parse.quote(json.dumps(filters))}"
        f"&fields={','.join(fields)}"
        "&size=600&format=json"
    )
    print(f"[GDC] fetching")
    with urllib.request.urlopen(url, timeout=120) as resp:
        d = json.loads(resp.read())
    rows = []
    for h in d["data"]["hits"]:
        diag = (h.get("diagnoses") or [{}])[0]
        rows.append({
            "patient12": h["submitter_id"],
            "tumor_focality": diag.get("tumor_focality"),
            "year_of_diagnosis": diag.get("year_of_diagnosis"),
            "ajcc_pathologic_t": diag.get("ajcc_pathologic_t"),
            "ajcc_pathologic_n": diag.get("ajcc_pathologic_n"),
            "ajcc_pathologic_m": diag.get("ajcc_pathologic_m"),
            "ajcc_pathologic_stage": diag.get("ajcc_pathologic_stage"),
            "bcr_tss": (h.get("tissue_source_site") or {}).get("bcr_id"),
        })
    df = pd.DataFrame(rows)
    df.to_csv(cache, sep="\t", index=False)
    print(f"[GDC] cached {len(df)} cases @ {cache}")
    return df


# ----------------------------------------------------------------------
# 4. Build full merged clinical
# ----------------------------------------------------------------------
def build_merged(h6, cbio, gdc):
    cbio = cbio.rename(columns={"patientId": "patient12"})
    keep = [
        "patient12",
        "EXTRATHYROIDAL_EXTENSION",
        "AJCC_PATHOLOGIC_TUMOR_STAGE",
        "AJCC_TUMOR_PATHOLOGIC_PT",
        "AJCC_NODES_PATHOLOGIC_PN",
        "AJCC_METASTASIS_PATHOLOGIC_PM",
        "INITIAL_PATHOLOGIC_DX_YEAR",
        "TISSUE_SOURCE_SITE",
        "NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT",
        "HISTOLOGICAL_DIAGNOSIS",
    ]
    cbio = cbio[[c for c in keep if c in cbio.columns]]

    df = h6.merge(cbio, on="patient12", how="left").merge(gdc, on="patient12", how="left")
    # h6 already has molecular_subtype + histology_subtype

    # ---- coerce / harmonize ----
    # Stage canonical (I/II/III/IV). Prefer cBio AJCC then GDC.
    def stage_norm(s):
        if pd.isna(s) or s in ("", "[Not Available]", "[Discrepancy]", "[Not Applicable]"):
            return np.nan
        s = str(s).upper().replace("STAGE ", "").strip()
        # strip A/B/C suffix
        for suf in ("A", "B", "C"):
            if s.endswith(suf) and len(s) > 1:
                s = s[:-1]
        if s in ("I", "II", "III", "IV"):
            return s
        return np.nan
    df["stage_canon"] = df["AJCC_PATHOLOGIC_TUMOR_STAGE"].map(stage_norm)
    df["stage_canon"] = df["stage_canon"].fillna(df["ajcc_pathologic_stage"].map(stage_norm))

    # T-stage: collapse T1a/T1b->T1, T2->T2, T3/T4 -> T3+
    def t_norm(s):
        if pd.isna(s) or s in ("", "TX", "T0", "[Not Available]", "[Discrepancy]"):
            return np.nan
        s = str(s).upper()
        if s.startswith("T1"): return "T1"
        if s.startswith("T2"): return "T2"
        if s.startswith("T3"): return "T3"
        if s.startswith("T4"): return "T4"
        return np.nan
    # Prefer cBio
    df["t_canon"] = df["AJCC_TUMOR_PATHOLOGIC_PT"].map(t_norm)
    df["t_canon"] = df["t_canon"].fillna(df["ajcc_pathologic_t"].map(t_norm))
    df["t_small"] = df["t_canon"].map({"T1": 1, "T2": 0, "T3": 0, "T4": 0})

    # N: N0 vs N1
    def n_norm(s):
        if pd.isna(s) or s in ("", "NX", "[Not Available]"):
            return np.nan
        s = str(s).upper()
        if s.startswith("N0"): return 0
        if s.startswith("N1"): return 1
        return np.nan
    df["n_pos"] = df["AJCC_NODES_PATHOLOGIC_PN"].map(n_norm)
    df["n_pos"] = df["n_pos"].fillna(df["ajcc_pathologic_n"].map(n_norm))

    # M: M0 vs M1
    def m_norm(s):
        if pd.isna(s) or s in ("", "MX", "[Not Available]"):
            return np.nan
        s = str(s).upper()
        if s.startswith("M0"): return 0
        if s.startswith("M1"): return 1
        return np.nan
    df["m_pos"] = df["AJCC_METASTASIS_PATHOLOGIC_PM"].map(m_norm)
    df["m_pos"] = df["m_pos"].fillna(df["ajcc_pathologic_m"].map(m_norm))

    # ETE — cBio thca_tcga only records ETE-positive cases (Minimal/Moderate/Advanced),
    # so NaN = no ETE recorded. We dichotomise as "advanced ETE" (T4-coded) vs others.
    # Derive ete_any: cBio non-null → 1, NaN AND t_canon in {T1,T2} → 0, NaN AND T3/T4 → drop.
    # Derive ete_advanced: Moderate/Very Advanced (T4) → 1, Minimal/None → 0.
    def ete_any(row):
        v = row["EXTRATHYROIDAL_EXTENSION"]
        if pd.notna(v) and str(v).strip() != "":
            return 1
        if row["t_canon"] in ("T1", "T2"):
            return 0
        return np.nan
    df["ete_pos"] = df.apply(ete_any, axis=1)

    def ete_adv(s):
        if pd.isna(s): return np.nan
        s = str(s).strip().lower()
        if "moderate" in s or "very advanced" in s or "t4" in s: return 1
        if "minimal" in s or "t3" in s: return 0
        return np.nan
    df["ete_advanced"] = df["EXTRATHYROIDAL_EXTENSION"].map(ete_adv)

    # Multifocality
    def focal_norm(s):
        if pd.isna(s): return np.nan
        s = str(s).strip().lower()
        if s in ("multifocal",): return 1
        if s in ("unifocal",): return 0
        return np.nan
    df["multifocal"] = df["tumor_focality"].map(focal_norm)

    # Year of diagnosis: prefer cBio INITIAL_PATHOLOGIC_DX_YEAR else GDC year_of_diagnosis
    def yr_num(s):
        try:
            v = float(str(s))
            if 1990 <= v <= 2025: return v
        except Exception:
            pass
        return np.nan
    df["dx_year"] = df["INITIAL_PATHOLOGIC_DX_YEAR"].map(yr_num)
    df["dx_year"] = df["dx_year"].fillna(df["year_of_diagnosis"].map(yr_num))
    df["era_post2010"] = (df["dx_year"] > 2010).astype("Int64")
    df.loc[df["dx_year"].isna(), "era_post2010"] = pd.NA

    # Tissue source site for clustering
    df["tss"] = df["TISSUE_SOURCE_SITE"].fillna(df["bcr_tss"]).astype(str)

    print(
        f"[merge] n={len(df)}; "
        f"stage_canon non-null={df['stage_canon'].notna().sum()}; "
        f"t_canon={df['t_canon'].notna().sum()}; "
        f"ete_pos={df['ete_pos'].notna().sum()}; "
        f"multifocal={df['multifocal'].notna().sum()}; "
        f"dx_year={df['dx_year'].notna().sum()}; "
        f"tss unique={df['tss'].nunique()}"
    )
    return df


# ----------------------------------------------------------------------
# 5. Cox helpers
# ----------------------------------------------------------------------
def cox_one(
    df: pd.DataFrame,
    duration: str,
    event: str,
    main_var: str,
    covars: list[str] | None = None,
    strata: list[str] | None = None,
    cluster_col: str | None = None,
    label: str = "",
):
    """Fit a Cox model, return dict of summary for `main_var`. Falls back to
    bumped penalty / no robust if convergence fails. Returns row dict ready
    for table append.
    """
    covars = covars or []
    cols = [duration, event, main_var] + [c for c in covars if c != main_var]
    if strata:
        cols += list(strata)
    if cluster_col:
        cols += [cluster_col]
    use = df.dropna(subset=cols).copy()
    use = use[use[duration] > 0]
    n = len(use)
    ev = int(use[event].sum())
    base = {
        "label": label, "main_var": main_var, "n": n, "events": ev,
        "covars": ",".join(covars), "strata": ",".join(strata or []),
        "cluster": cluster_col or "",
        "HR": np.nan, "CI_low": np.nan, "CI_high": np.nan, "p": np.nan,
        "note": "",
    }
    if n < 8 or ev < 3:
        base["note"] = "insufficient n/events"
        return base

    fit_ok = False
    last_err = ""
    for pen in (0.0, 0.05, 0.2):
        cph = CoxPHFitter(penalizer=pen)
        try:
            kwargs = dict(duration_col=duration, event_col=event, robust=True)
            if strata:
                kwargs["strata"] = list(strata)
            if cluster_col:
                kwargs["cluster_col"] = cluster_col
            fit_df = use[[duration, event, main_var] + [c for c in covars if c != main_var]
                         + (list(strata) if strata else [])
                         + ([cluster_col] if cluster_col else [])]
            cph.fit(fit_df, **kwargs)
            fit_ok = True
            break
        except Exception as e:
            last_err = f"{type(e).__name__}: {str(e)[:80]}"
            continue
    if not fit_ok:
        base["note"] = f"fit_failed: {last_err}"
        return base
    s = cph.summary
    if main_var not in s.index:
        base["note"] = "main_var dropped"
        return base
    HR = float(s.loc[main_var, "exp(coef)"])
    CI_lo = float(s.loc[main_var, "exp(coef) lower 95%"])
    CI_hi = float(s.loc[main_var, "exp(coef) upper 95%"])
    base["HR"] = HR
    base["CI_low"] = CI_lo
    base["CI_high"] = CI_hi
    base["p"] = float(s.loc[main_var, "p"])
    notes = []
    # cell-size flag: count main_var=1 events, often the limiting cell
    if main_var in use.columns:
        cell1_n = int((use[main_var] == 1).sum())
        cell1_ev = int(use.loc[use[main_var] == 1, event].sum())
        if cell1_n < 5 or cell1_ev < 2:
            notes.append(f"cell main_var=1 sparse (n={cell1_n}, ev={cell1_ev})")
    if ev < 10:
        notes.append("n_events<10 (CI non-deployable)")
    # fit-pathology flag (HR essentially 0 or huge with extreme CI ratio)
    if (HR < 1e-3 or HR > 1e3 or (CI_hi / max(CI_lo, 1e-12) > 1e6)):
        notes.append("fit-pathology (extreme HR/CI - drop)")
    base["note"] = "; ".join(notes)
    return base


# ----------------------------------------------------------------------
# 6. Main analysis battery
# ----------------------------------------------------------------------
def main():
    h6 = load_h6()
    cbio = fetch_cbio_thca_tcga()
    gdc = fetch_gdc_thca()
    df = build_merged(h6, cbio, gdc)

    # save merged table
    df.to_csv(OUT / "h24_merged_with_external_clinical.tsv", sep="\t", index=False)

    # BRAF-cPTC mask (mirrors H6 S1)
    S1 = (df["molecular_subtype"] == "BRAF_like") & (df["histology_subtype"] == "cPTC")
    print(f"[S1 BRAF-cPTC] n={int(S1.sum())}")

    # DM2 vs not_DM contrast frame
    dm2_mask = df["dm"].isin(["DM2", "not_DM"])
    df["is_DM2"] = (df["dm"] == "DM2").astype(int)
    # DM1 vs DM2 (DM-positive only)
    df["is_DM1_vs_DM2"] = np.where(
        df["dm"] == "DM1", 1, np.where(df["dm"] == "DM2", 0, np.nan)
    )

    rows: list[dict] = []
    inter_rows: list[dict] = []

    BASE_COVARS = ["age_yr", "sex_m"]            # main effect contrast = is_DM2 (vs not_DM)
    BASE_COVARS_STAGE = BASE_COVARS + ["adv_stage"]  # H6 main spec

    # ---------- (A) main-spec replication of H6 PFI HR=5.91 ----------
    sub = df[S1 & dm2_mask]
    rows.append(cox_one(
        sub, "PFI.time", "PFI", "is_DM2",
        covars=BASE_COVARS_STAGE,
        label="A1_main_DM2vNot_PFI_replicate",
    ))

    # ---------- (B) Stage-stratified Cox (PFI) ----------
    rows.append(cox_one(
        sub, "PFI.time", "PFI", "is_DM2",
        covars=BASE_COVARS,                      # stage in strata, not as covariate
        strata=["stage_canon"],
        label="B1_stage_strata_DM2vNot_PFI",
    ))
    rows.append(cox_one(
        sub, "PFI.time", "PFI", "is_DM2",
        covars=BASE_COVARS,
        strata=["t_canon"],
        label="B2_T_strata_DM2vNot_PFI",
    ))
    # ----- per-stage Cox slices -----
    for st in ["I", "II", "III", "IV"]:
        s_st = sub[sub["stage_canon"] == st]
        rows.append(cox_one(
            s_st, "PFI.time", "PFI", "is_DM2",
            covars=BASE_COVARS,
            label=f"B3_within_stage_{st}_DM2vNot_PFI",
        ))

    # ---------- (C) T-stage strata + interaction ----------
    sub_t1 = sub[sub["t_canon"] == "T1"]
    sub_t234 = sub[sub["t_canon"].isin(["T2", "T3", "T4"])]
    rows.append(cox_one(
        sub_t1, "PFI.time", "PFI", "is_DM2",
        covars=BASE_COVARS_STAGE,
        label="C1_T1only_DM2vNot_PFI",
    ))
    rows.append(cox_one(
        sub_t234, "PFI.time", "PFI", "is_DM2",
        covars=BASE_COVARS_STAGE,
        label="C2_T2to4_DM2vNot_PFI",
    ))
    # Interaction: is_DM2 * t_small (1 if T1, 0 else)
    sub_int = sub.dropna(subset=["t_small"]).copy()
    sub_int["dm2_x_t1"] = sub_int["is_DM2"] * sub_int["t_small"]
    inter = cox_one(
        sub_int, "PFI.time", "PFI", "dm2_x_t1",
        covars=BASE_COVARS_STAGE + ["is_DM2", "t_small"],
        label="C3_interaction_DM2_x_T1",
    )
    rows.append(inter)
    inter_rows.append({"feature": "T_small_(T1)", **inter})

    # ---------- (D) Anatomic features × DM ----------
    for fname, fcol, levels in [
        ("ETE", "ete_pos", [(0, "neg"), (1, "pos")]),
        ("N", "n_pos", [(0, "N0"), (1, "N1")]),
        ("M", "m_pos", [(0, "M0"), (1, "M1")]),
        ("multifocal", "multifocal", [(0, "uni"), (1, "multi")]),
    ]:
        for v, lab in levels:
            ssub = sub[sub[fcol] == v]
            rows.append(cox_one(
                ssub, "PFI.time", "PFI", "is_DM2",
                covars=BASE_COVARS_STAGE,
                label=f"D_{fname}_{lab}_DM2vNot_PFI",
            ))
        # interaction term
        sint = sub.dropna(subset=[fcol]).copy()
        sint[f"dm2_x_{fname}"] = sint["is_DM2"] * sint[fcol]
        inter = cox_one(
            sint, "PFI.time", "PFI", f"dm2_x_{fname}",
            covars=BASE_COVARS_STAGE + ["is_DM2", fcol],
            label=f"D_interaction_DM2_x_{fname}",
        )
        rows.append(inter)
        inter_rows.append({"feature": fname, **inter})

    # ---------- (E) Era split ----------
    for era_v, lab in [(0, "le2010"), (1, "post2010")]:
        ssub = sub[sub["era_post2010"] == era_v]
        rows.append(cox_one(
            ssub, "PFI.time", "PFI", "is_DM2",
            covars=BASE_COVARS_STAGE,
            label=f"E_era_{lab}_DM2vNot_PFI",
        ))
    sint = sub.dropna(subset=["era_post2010"]).copy()
    sint["era_post2010"] = sint["era_post2010"].astype(int)
    sint["dm2_x_era"] = sint["is_DM2"] * sint["era_post2010"]
    inter = cox_one(
        sint, "PFI.time", "PFI", "dm2_x_era",
        covars=BASE_COVARS_STAGE + ["is_DM2", "era_post2010"],
        label="E_interaction_DM2_x_post2010",
    )
    rows.append(inter)
    inter_rows.append({"feature": "era_post2010", **inter})

    # ---------- (F) DM1 vs DM2 protective sensitivity ----------
    sub_pos = df[S1 & df["dm"].isin(["DM1", "DM2"])].copy()
    sub_pos["is_DM1"] = (sub_pos["dm"] == "DM1").astype(int)
    rows.append(cox_one(
        sub_pos, "PFI.time", "PFI", "is_DM1",
        covars=BASE_COVARS_STAGE,
        label="F1_main_DM1vDM2_PFI_replicate",
    ))
    rows.append(cox_one(
        sub_pos, "PFI.time", "PFI", "is_DM1",
        covars=BASE_COVARS,
        strata=["stage_canon"],
        label="F2_stage_strata_DM1vDM2_PFI",
    ))
    rows.append(cox_one(
        sub_pos, "PFI.time", "PFI", "is_DM1",
        covars=BASE_COVARS,
        strata=["t_canon"],
        label="F3_T_strata_DM1vDM2_PFI",
    ))
    for era_v, lab in [(0, "le2010"), (1, "post2010")]:
        ssub = sub_pos[sub_pos["era_post2010"] == era_v]
        rows.append(cox_one(
            ssub, "PFI.time", "PFI", "is_DM1",
            covars=BASE_COVARS_STAGE,
            label=f"F4_era_{lab}_DM1vDM2_PFI",
        ))
    # interaction era for DM1 vs DM2
    sint = sub_pos.dropna(subset=["era_post2010"]).copy()
    sint["era_post2010"] = sint["era_post2010"].astype(int)
    sint["dm1_x_era"] = sint["is_DM1"] * sint["era_post2010"]
    inter = cox_one(
        sint, "PFI.time", "PFI", "dm1_x_era",
        covars=BASE_COVARS_STAGE + ["is_DM1", "era_post2010"],
        label="F_interaction_DM1_x_post2010",
    )
    rows.append(inter)
    inter_rows.append({"feature": "DM1_x_era_post2010", **inter})

    # ---------- (G) Robust SE / cluster by tissue source site ----------
    # cluster column needs >=2 levels and at least one event per cluster ideally
    sub_cl = sub.copy()
    sub_cl["tss"] = sub_cl["tss"].fillna("UNK")
    rows.append(cox_one(
        sub_cl, "PFI.time", "PFI", "is_DM2",
        covars=BASE_COVARS_STAGE,
        cluster_col="tss",
        label="G1_cluster_TSS_DM2vNot_PFI",
    ))

    # ---------- table ----------
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "h24_sensitivity_table.tsv", sep="\t", index=False)
    print(f"[ok] wrote {OUT/'h24_sensitivity_table.tsv'} (rows={len(tab)})")

    inter_tab = pd.DataFrame(inter_rows)
    inter_tab["interaction_p"] = inter_tab["p"]
    inter_tab[["feature", "label", "n", "events", "HR", "CI_low", "CI_high", "interaction_p", "note"]].to_csv(
        OUT / "h24_interactions.tsv", sep="\t", index=False
    )
    print(f"[ok] wrote {OUT/'h24_interactions.tsv'}")

    # ---------- forest plot data ----------
    # Curate one row per subgroup for the DM2-vs-not_DM PFI contrast.
    SUBG_LABELS = [
        ("Overall (BRAF-cPTC, H6 main)", "A1_main_DM2vNot_PFI_replicate"),
        ("Stage-stratified", "B1_stage_strata_DM2vNot_PFI"),
        ("T-stratified", "B2_T_strata_DM2vNot_PFI"),
        ("Stage I", "B3_within_stage_I_DM2vNot_PFI"),
        ("Stage II", "B3_within_stage_II_DM2vNot_PFI"),
        ("Stage III", "B3_within_stage_III_DM2vNot_PFI"),
        ("Stage IV", "B3_within_stage_IV_DM2vNot_PFI"),
        ("T1", "C1_T1only_DM2vNot_PFI"),
        ("T2-T4", "C2_T2to4_DM2vNot_PFI"),
        ("ETE neg", "D_ETE_neg_DM2vNot_PFI"),
        ("ETE pos", "D_ETE_pos_DM2vNot_PFI"),
        ("N0", "D_N_N0_DM2vNot_PFI"),
        ("N1 (LN+)", "D_N_N1_DM2vNot_PFI"),
        ("M0", "D_M_M0_DM2vNot_PFI"),
        ("M1 (mets+)", "D_M_M1_DM2vNot_PFI"),
        ("Unifocal", "D_multifocal_uni_DM2vNot_PFI"),
        ("Multifocal", "D_multifocal_multi_DM2vNot_PFI"),
        ("Dx <=2010", "E_era_le2010_DM2vNot_PFI"),
        ("Dx >2010", "E_era_post2010_DM2vNot_PFI"),
        ("Cluster by TSS (robust)", "G1_cluster_TSS_DM2vNot_PFI"),
    ]
    forest = []
    for nice, lab in SUBG_LABELS:
        r = tab[tab["label"] == lab]
        if len(r) == 0:
            continue
        r0 = r.iloc[0]
        forest.append({
            "subgroup": nice,
            "label_id": lab,
            "n": int(r0["n"]),
            "events": int(r0["events"]),
            "HR": r0["HR"],
            "CI_low": r0["CI_low"],
            "CI_high": r0["CI_high"],
            "p": r0["p"],
            "note": r0["note"],
        })
    f_df = pd.DataFrame(forest)
    f_df.to_csv(OUT / "h24_forest_data.tsv", sep="\t", index=False)
    print(f"[ok] wrote {OUT/'h24_forest_data.tsv'} (rows={len(f_df)})")

    # ---------- forest plot rendering ----------
    plot_forest(f_df, OUT / "h24_forest_plot.png")

    # ---------- sentinel: how many subgroups keep HR>1 and lower-CI>1?
    # We track two tiers:
    #   (i) ALL evaluable subgroups (HR fit, not fit-pathology) — broadest count
    #   (ii) DEPLOYABLE subgroups: HR fit, n_events>=10, no cell-sparsity flag,
    #        no fit-pathology — the strict count for the headline verdict.
    valid = tab[tab["label"].isin([s[1] for s in SUBG_LABELS])]
    valid = valid.dropna(subset=["HR", "CI_low"]).copy()
    valid["is_pathology"] = valid["note"].astype(str).str.contains("fit-pathology", na=False)
    valid["is_sparse"] = valid["note"].astype(str).str.contains("sparse", na=False)
    valid["is_underpowered"] = valid["note"].astype(str).str.contains("non-deployable", na=False)

    evaluable = valid[~valid["is_pathology"]]
    deployable = valid[(~valid["is_pathology"]) & (~valid["is_sparse"]) & (~valid["is_underpowered"])]

    n_total = len(evaluable)
    n_hr_gt_1 = int((evaluable["HR"] > 1).sum())
    n_ci_gt_1 = int((evaluable["CI_low"] > 1).sum())
    n_ci_crosses_1 = int((evaluable["CI_low"] <= 1).sum())

    n_dep_total = len(deployable)
    n_dep_hr_gt_1 = int((deployable["HR"] > 1).sum())
    n_dep_ci_gt_1 = int((deployable["CI_low"] > 1).sum())

    summary = {
        "n_subgroups_evaluable": n_total,
        "n_HR_gt_1": n_hr_gt_1,
        "n_CIlow_gt_1_strict": n_ci_gt_1,
        "n_CI_crosses_1": n_ci_crosses_1,
        "DEPLOYABLE_n_total": n_dep_total,
        "DEPLOYABLE_n_HR_gt_1": n_dep_hr_gt_1,
        "DEPLOYABLE_n_CIlow_gt_1": n_dep_ci_gt_1,
        "DEPLOYABLE_subgroups": deployable["subgroup"].tolist() if "subgroup" in deployable.columns else deployable["label"].tolist(),
    }
    with open(OUT / "h24_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=lambda x: None if pd.isna(x) else float(x))
    print(f"[ok] summary: {summary}")

    return tab, f_df, inter_tab, summary


def plot_forest(df: pd.DataFrame, fpath: Path):
    if len(df) == 0:
        return
    df = df.copy()
    df = df[~df["note"].astype(str).str.contains("fit-pathology", na=False)]
    df = df.dropna(subset=["HR", "CI_low", "CI_high"])
    if len(df) == 0:
        return
    # Cap CI for plot only
    df["CI_low_cap"] = df["CI_low"].clip(lower=0.05)
    df["CI_high_cap"] = df["CI_high"].clip(upper=200)
    fig, ax = plt.subplots(figsize=(8.5, 0.42 * len(df) + 1.2))
    y = np.arange(len(df))[::-1]
    ax.errorbar(
        df["HR"], y,
        xerr=[(df["HR"] - df["CI_low_cap"]).clip(lower=0),
              (df["CI_high_cap"] - df["HR"]).clip(lower=0)],
        fmt="s", color="#222", ecolor="#888", capsize=3, markersize=5,
    )
    ax.set_xscale("log")
    ax.axvline(1.0, color="r", ls="--", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(df["subgroup"].tolist(), fontsize=9)
    ax.set_xlabel("HR (DM2 vs not_DM, PFI) [log scale]")
    ax.set_xlim(0.05, 250)
    # annotate n/events
    for yi, (_, r) in zip(y, df.iterrows()):
        ax.text(
            220, yi,
            f"n={int(r['n'])}, ev={int(r['events'])}",
            va="center", fontsize=8, color="#444",
        )
    ax.set_title("H24 — DM2-vs-not_DM PFI HR sensitivity (BRAF-cPTC)")
    plt.tight_layout()
    plt.savefig(fpath, dpi=180)
    plt.close()
    print(f"[ok] wrote {fpath}")


if __name__ == "__main__":
    main()
