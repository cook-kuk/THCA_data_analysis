#!/usr/bin/env python3
"""H9 — MAPK-inhibitor sensitivity in DM1 thyroid cell lines (Paper 1+2 BRAF Nature sprint).

Inputs (DepMap Public 24Q2, on-disk):
  /data/thca/repo_results/p3_p9_full_execution/paper9/raw/Model.csv
  /data/thca/repo_results/p3_p9_full_execution/paper9/raw/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
  /data/thca/repo_results/p3_p9_full_execution/paper9/raw/CRISPRGeneEffect.csv
  /data/thca/repo_results/p3_p9_full_execution/paper9/raw/Repurposing_Public_24Q2_LFC_COLLAPSED.csv
  /data/thca/repo_results/p3_p9_full_execution/paper9/raw/Repurposing_Public_24Q2_Treatment_Meta_Data.csv

DM1 score (per cell line):
  HT-13 + FA-12 panel z-score sum, where:
    HT-13 = HLA-DRA, HLA-DRB1, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1,
            CD79A, CD79B, MS4A1, AICDA, CXCL13, CCR6, IFNG
    FA-12 = FASN, ACACA, ACLY, SCD, FADS1, FADS2, ELOVL6, ACOX1,
            CPT1A, HMGCS2, HADH, ACADM
  Per-gene z-score across thyroid panel; sum HT z-scores, sum FA z-scores,
  DM1_score = HT_z_sum  -  FA_z_sum  (HT-up + FA-down direction = DM1).
  (Matches H2/H4 sprint convention: HT-axis up + FA/lipid-output down = DM1.)

Outputs:
  h9_thyroid_dm1_score.tsv          # per-line DM1 score
  h9_depmap_essentiality_top.tsv    # CRISPR Spearman + DM1-high vs DM1-low t-test, all genes
  h9_prism_drug_ranking.tsv         # PRISM Spearman + DM1-high vs DM1-low t-test, all drugs + MAPK class flag
  h9_mapk_drug_focus.tsv            # focused MEK/BRAF/ERK inhibitor table
  H9_REPORT.md                      # <400 word report
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

OUT = Path(__file__).resolve().parent
RAW = Path("/data/thca/repo_results/p3_p9_full_execution/paper9/raw")

MODEL_CSV = RAW / "Model.csv"
EXPR_CSV  = RAW / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
CRISPR_CSV = RAW / "CRISPRGeneEffect.csv"
PRISM_LFC = RAW / "Repurposing_Public_24Q2_LFC_COLLAPSED.csv"
PRISM_META = RAW / "Repurposing_Public_24Q2_Treatment_Meta_Data.csv"

HT13 = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
        "CD79A", "CD79B", "MS4A1", "AICDA", "CXCL13", "CCR6", "IFNG"]
FA12 = ["FASN", "ACACA", "ACLY", "SCD", "FADS1", "FADS2", "ELOVL6", "ACOX1",
        "CPT1A", "HMGCS2", "HADH", "ACADM"]

# Canonical MEK/BRAF/ERK inhibitor list (we filter case-insensitive on `name`)
MEK = {"selumetinib", "trametinib", "azd-0364", "azd0364", "cobimetinib",
       "binimetinib", "pimasertib", "refametinib", "pd0325901", "pd-0325901",
       "u0126", "u-0126", "ci-1040", "mirdametinib", "pd-184352",
       "pd0184352", "tak-733"}
BRAF = {"dabrafenib", "vemurafenib", "encorafenib", "raf265", "lgx-818",
        "plx4720", "plx-4720", "plx-4032", "sb590885", "sb-590885",
        "raf709", "az-628", "az628", "lifirafenib", "tovorafenib",
        "bgb-3245", "bgb3245", "regorafenib", "sorafenib"}
ERK = {"sch772984", "sch-772984", "ulixertinib", "bvd-523", "bvd523",
       "ly3214996", "ravoxertinib", "gdc-0994", "gdc0994", "vx-11e",
       "vx11e", "ml786", "cc-90003", "cc90003"}


def log(msg: str) -> None:
    print(f"[H9] {msg}", flush=True)


def load_thyroid_models() -> pd.DataFrame:
    log("loading Model.csv (thyroid filter)")
    m = pd.read_csv(MODEL_CSV, low_memory=False)
    th = m[m["OncotreeLineage"] == "Thyroid"].copy()
    log(f"  thyroid cell lines: {len(th)}")
    return th


def parse_expr_columns() -> dict[str, str]:
    """Return mapping HUGO -> exact column name in expression file.

    Columns are: ['Unnamed: 0', 'SequencingID', 'ModelConditionID', 'ModelID',
    'IsDefaultEntryForMC', '<gene1> (<entrez1>)', ...].
    """
    log("parsing expression header")
    hdr = pd.read_csv(EXPR_CSV, nrows=0)
    cols = list(hdr.columns)
    META_COLS = {"Unnamed: 0", "SequencingID", "ModelConditionID",
                 "ModelID", "IsDefaultEntryForMC"}
    mapping: dict[str, str] = {}
    for c in cols:
        if c in META_COLS:
            continue
        gene = c.split(" (")[0]
        mapping[gene] = c
    log(f"  parsed {len(mapping):,} gene columns")
    return mapping


def load_panel_expression(thyroid_ids: list[str], gene_map: dict[str, str]
                          ) -> pd.DataFrame:
    """Load HT13 + FA12 expression for thyroid models only."""
    panel = HT13 + FA12
    found = [g for g in panel if g in gene_map]
    missing = [g for g in panel if g not in gene_map]
    log(f"  panel genes resolved: {len(found)}/{len(panel)} (missing: {missing})")
    use_cols = ["ModelID", "IsDefaultEntryForMC"] + [gene_map[g] for g in found]
    log("  reading panel-only columns from expression file (~6 MB)")
    df = pd.read_csv(EXPR_CSV, usecols=use_cols)
    df = df[df["ModelID"].isin(thyroid_ids)].copy()
    # If multiple sequencing entries per ModelID, keep IsDefault
    if "IsDefaultEntryForMC" in df.columns:
        if df["IsDefaultEntryForMC"].notna().any():
            d_def = df[df["IsDefaultEntryForMC"] == True]
            if len(d_def) >= 1:
                df = d_def
        df = df.drop_duplicates(subset="ModelID", keep="first")
        df = df.drop(columns=["IsDefaultEntryForMC"])
    # Rename gene cols back to HUGO
    inv = {gene_map[g]: g for g in found}
    df = df.rename(columns=inv)
    log(f"  thyroid lines with expression: {len(df)}")
    return df, found, missing


def compute_dm1_score(expr: pd.DataFrame, found: list[str]) -> pd.DataFrame:
    """DM1 score = z(HT genes).sum() - z(FA genes).sum() across thyroid panel."""
    ht_used = [g for g in HT13 if g in found]
    fa_used = [g for g in FA12 if g in found]
    log(f"  z-score across {len(expr)} thyroid lines: HT={len(ht_used)}, FA={len(fa_used)}")

    z = expr.copy()
    for g in ht_used + fa_used:
        v = z[g].astype(float)
        z[g] = (v - v.mean()) / (v.std(ddof=0) if v.std(ddof=0) > 0 else 1.0)
    z["HT_z_sum"] = z[ht_used].sum(axis=1)
    z["FA_z_sum"] = z[fa_used].sum(axis=1)
    z["DM1_score"] = z["HT_z_sum"] - z["FA_z_sum"]
    return z[["ModelID", "HT_z_sum", "FA_z_sum", "DM1_score"]]


# ----------------------------------------------------------------------
# CRISPR essentiality
# ----------------------------------------------------------------------
def crispr_essentiality(dm1: pd.DataFrame) -> pd.DataFrame:
    """Spearman per gene + DM1-high vs DM1-low t-test on CRISPR effect."""
    thyroid_ids = set(dm1["ModelID"])
    log(f"loading CRISPR effect (filtered to {len(thyroid_ids)} thyroid lines)")

    # The CSV has unnamed first col with ModelIDs + many gene columns. 421 MB -> stream.
    hdr = pd.read_csv(CRISPR_CSV, nrows=0).columns.tolist()
    first_col = hdr[0]  # 'Unnamed: 0'
    # Iterate in chunks on rows.
    chunks = []
    for chunk in pd.read_csv(CRISPR_CSV, chunksize=200, low_memory=False):
        sub = chunk[chunk[first_col].isin(thyroid_ids)]
        if len(sub):
            chunks.append(sub)
    crispr = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
    crispr = crispr.rename(columns={first_col: "ModelID"})
    # Strip "(entrez)" from gene column names for downstream readability
    rename_map = {c: c.split(" (")[0]
                  for c in crispr.columns if c != "ModelID"}
    crispr = crispr.rename(columns=rename_map)
    log(f"  CRISPR thyroid rows: {len(crispr)}")

    merged = dm1.merge(crispr, on="ModelID", how="inner")
    log(f"  merged with DM1 score: n={len(merged)}")

    gene_cols = [c for c in merged.columns
                 if c not in ("ModelID", "HT_z_sum", "FA_z_sum", "DM1_score")]
    log(f"  CRISPR genes to test: {len(gene_cols):,}")

    dm1_vec = merged["DM1_score"].to_numpy()
    median = float(np.median(dm1_vec))
    high_mask = dm1_vec > median
    low_mask = dm1_vec <= median
    n_high = int(high_mask.sum())
    n_low = int(low_mask.sum())
    log(f"  DM1-high n={n_high}, DM1-low n={n_low} (median split @ {median:.3f})")

    Y = merged[gene_cols].to_numpy()  # (n_lines, n_genes)
    # Spearman per gene (vectorized via rank)
    n = Y.shape[0]
    valid_mask = ~np.isnan(Y)
    rec = []
    # Pre-rank DM1
    dm1_rank = stats.rankdata(dm1_vec)
    dm1_rank_centered = dm1_rank - dm1_rank.mean()
    for j, gene in enumerate(gene_cols):
        y = Y[:, j]
        ok = ~np.isnan(y)
        if ok.sum() < 8:
            continue
        # Spearman with available mask (re-rank within valid)
        if ok.sum() == n:
            yr = stats.rankdata(y) - stats.rankdata(y).mean()
            num = (dm1_rank_centered * yr).sum()
            den = np.sqrt((dm1_rank_centered ** 2).sum() * (yr ** 2).sum())
            r = num / den if den > 0 else np.nan
            # Two-sided p
            df_ = n - 2
            t = r * np.sqrt(df_ / max(1e-12, 1 - r ** 2))
            p_sp = 2 * stats.t.sf(abs(t), df_)
        else:
            r, p_sp = stats.spearmanr(dm1_vec[ok], y[ok])

        # high vs low t-test
        h = y[ok & high_mask]
        l = y[ok & low_mask]
        if len(h) >= 3 and len(l) >= 3:
            t_stat, p_t = stats.ttest_ind(h, l, equal_var=False)
            d = ((h.mean() - l.mean()) /
                 np.sqrt((h.var(ddof=1) + l.var(ddof=1)) / 2))
        else:
            t_stat, p_t, d = np.nan, np.nan, np.nan

        rec.append((gene, ok.sum(), len(h), len(l),
                    h.mean() if len(h) else np.nan,
                    l.mean() if len(l) else np.nan,
                    r, p_sp, t_stat, p_t, d))

    out = pd.DataFrame(rec, columns=["gene", "n", "n_high", "n_low",
                                     "mean_eff_high", "mean_eff_low",
                                     "spearman_r", "spearman_p",
                                     "t_stat", "t_p", "cohens_d"])
    # FDR per test
    out["spearman_fdr"] = multipletests(out["spearman_p"].fillna(1.0),
                                        method="fdr_bh")[1]
    out["t_fdr"] = multipletests(out["t_p"].fillna(1.0), method="fdr_bh")[1]
    out = out.sort_values("spearman_r")  # most-essential-in-high first (negative effect = essential)
    return out


# ----------------------------------------------------------------------
# PRISM drug response
# ----------------------------------------------------------------------
def prism_drug_response(dm1: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-drug Spearman + DM1-high vs DM1-low; flag MEK/BRAF/ERK class.

    PRISM LFC is in long format (1.5M rows): row_id contains ACH-id::compound::...
    plus a `broad_id` and `LFC` column. We average LFC per (ModelID, broad_id)
    across doses/screens, then test against DM1 score.
    """
    thyroid_ids = set(dm1["ModelID"])
    log("loading PRISM Repurposing 24Q2 LFC (long format)")
    prism = pd.read_csv(PRISM_LFC, usecols=["row_id", "broad_id", "dose", "LFC"])
    prism["ModelID"] = prism["row_id"].str.split("::").str[0]
    log(f"  PRISM rows total: {len(prism):,}, unique cells: {prism.ModelID.nunique()}, "
        f"unique drugs: {prism.broad_id.nunique()}")
    prism = prism[prism["ModelID"].isin(thyroid_ids)].copy()
    log(f"  thyroid-line rows: {len(prism):,}, thyroid lines present: "
        f"{prism.ModelID.nunique()}")
    # QC filter
    prism = prism[~prism["broad_id"].astype(str).str.contains("QC Failure",
                                                              na=False)]

    # Average LFC per (ModelID, broad_id) over doses & screens
    cd = (prism.groupby(["ModelID", "broad_id"])["LFC"]
                .mean().reset_index())
    log(f"  cell × drug rows: {len(cd):,}")

    # Wide pivot: ModelID rows, drug cols
    wide = cd.pivot(index="ModelID", columns="broad_id", values="LFC")
    wide = wide.reset_index()
    merged = dm1.merge(wide, on="ModelID", how="inner")
    log(f"  merged with DM1 score: n={len(merged)} thyroid lines, "
        f"{wide.shape[1]-1} drugs")

    drug_cols = [c for c in merged.columns
                 if c not in ("ModelID", "HT_z_sum", "FA_z_sum", "DM1_score")]

    dm1_vec = merged["DM1_score"].to_numpy()
    median = float(np.median(dm1_vec))
    high_mask = dm1_vec > median
    low_mask = dm1_vec <= median
    log(f"  DM1-high n={int(high_mask.sum())}, DM1-low n={int(low_mask.sum())} "
        f"(median split @ {median:.3f})")
    n = len(dm1_vec)
    dm1_rank = stats.rankdata(dm1_vec)
    dm1_rc = dm1_rank - dm1_rank.mean()

    Y = merged[drug_cols].to_numpy(dtype=float)
    rec = []
    for j, drug in enumerate(drug_cols):
        y = Y[:, j]
        ok = ~np.isnan(y)
        if ok.sum() < 6:
            continue
        if ok.sum() == n:
            yr = stats.rankdata(y) - stats.rankdata(y).mean()
            num = (dm1_rc * yr).sum()
            den = np.sqrt((dm1_rc ** 2).sum() * (yr ** 2).sum())
            r = num / den if den > 0 else np.nan
            df_ = n - 2
            t = r * np.sqrt(df_ / max(1e-12, 1 - r ** 2))
            p_sp = 2 * stats.t.sf(abs(t), df_)
        else:
            r, p_sp = stats.spearmanr(dm1_vec[ok], y[ok])

        h = y[ok & high_mask]
        l = y[ok & low_mask]
        if len(h) >= 3 and len(l) >= 3:
            t_stat, p_t = stats.ttest_ind(h, l, equal_var=False)
            pooled_var = (h.var(ddof=1) + l.var(ddof=1)) / 2
            d = ((h.mean() - l.mean()) /
                 np.sqrt(pooled_var)) if pooled_var > 0 else np.nan
        else:
            t_stat, p_t, d = np.nan, np.nan, np.nan
        rec.append((drug, ok.sum(), len(h), len(l),
                    h.mean() if len(h) else np.nan,
                    l.mean() if len(l) else np.nan,
                    r, p_sp, t_stat, p_t, d))

    out = pd.DataFrame(rec, columns=["broad_id", "n", "n_high", "n_low",
                                     "mean_LFC_high", "mean_LFC_low",
                                     "spearman_r", "spearman_p",
                                     "t_stat", "t_p", "cohens_d"])

    log("  joining drug metadata")
    meta = pd.read_csv(PRISM_META, low_memory=False)
    log(f"  PRISM meta cols: {meta.columns.tolist()}")
    # Pick name + moa + target (Repurposing 24Q2 schema may not have moa here;
    # check)
    bid_col = "broad_id"
    name_col = "name" if "name" in meta.columns else None
    moa_col = next((c for c in meta.columns
                    if "moa" in c.lower() or "mechanism" in c.lower()), None)
    target_col = next((c for c in meta.columns
                       if c.lower() in ("target", "targets", "gene_target")), None)
    keep = [bid_col] + [c for c in (name_col, moa_col, target_col) if c]
    meta_small = meta[keep].drop_duplicates(subset=[bid_col])
    rename = {}
    if name_col: rename[name_col] = "name"
    if moa_col: rename[moa_col] = "moa"
    if target_col: rename[target_col] = "target"
    meta_small = meta_small.rename(columns=rename)
    out = out.merge(meta_small, on="broad_id", how="left")
    if "moa" not in out.columns: out["moa"] = ""
    if "target" not in out.columns: out["target"] = ""
    if "name" not in out.columns: out["name"] = ""

    # FDR
    out["spearman_fdr"] = multipletests(out["spearman_p"].fillna(1.0),
                                        method="fdr_bh")[1]
    out["t_fdr"] = multipletests(out["t_p"].fillna(1.0), method="fdr_bh")[1]

    # Class flag
    def classify(row) -> str:
        nm = str(row.get("name", "")).lower()
        moa = str(row.get("moa", "")).lower()
        tgt = str(row.get("target", "")).lower()
        if any(k in nm for k in MEK) or "mek" in moa or "map2k" in tgt:
            return "MEK"
        if any(k in nm for k in BRAF) or ("braf" in moa and "inhibitor" in moa) \
                or "braf" in tgt:
            return "BRAF/RAF"
        if any(k in nm for k in ERK) or "erk" in moa or "mapk1" in tgt or "mapk3" in tgt:
            return "ERK"
        if "mapk" in moa:
            return "MAPK_other"
        return ""

    out["mapk_class"] = out.apply(classify, axis=1)

    # Sort: most DM1-high-selective at top (cohens_d most negative = LFC lower in high = more sensitive)
    out = out.sort_values("cohens_d")

    # Focused table: MEK/BRAF/ERK subset
    focus = out[out["mapk_class"].isin(["MEK", "BRAF/RAF", "ERK"])].copy()
    focus = focus.sort_values("cohens_d")
    return out, focus


# ----------------------------------------------------------------------
# GDSC + CTRP — clinical MAPK monotherapy (dabrafenib, trametinib, SCH772984, VX-11E)
# ----------------------------------------------------------------------
GDSC_CTRP = {
    "GDSC2": {
        "auc":  RAW / "GDSC2AUCMatrix.csv",
        "ic50": RAW / "GDSC2Log2IC50Matrix.csv",
        "viab": RAW / "GDSC2Log2ViabilityCollapsedConditions.csv",
    },
    "CTRP": {
        "auc":  RAW / "CTRPAUCMatrix.csv",
        "ic50": RAW / "CTRPLog2IC50Matrix.csv",
        "viab": RAW / "CTRPLog2ViabilityCollapsedConditions.csv",
    },
    "GDSC1": {
        "auc":  RAW / "GDSC1AUCMatrix.csv",
        "ic50": RAW / "GDSC1Log2IC50Matrix.csv",
        "viab": RAW / "GDSC1Log2ViabilityCollapsedConditions.csv",
    },
}
CLINICAL_MAPK = ["dabrafenib", "trametinib", "vemurafenib", "selumetinib",
                 "cobimetinib", "binimetinib", "encorafenib", "regorafenib",
                 "sorafenib", "sch772984", "vx-11e", "vx11e"]


def gdsc_ctrp_clinical(dm1: pd.DataFrame) -> pd.DataFrame | None:
    """For each clinical MAPK drug present in GDSC1/2/CTRP, compute Spearman
    of AUC vs DM1 score and DM1-high vs DM1-low t-test on AUC.

    AUC: lower = more sensitive; so cohens_d < 0 means DM1-high more sensitive.
    """
    rec = []
    th_ids = set(dm1.ModelID)
    dm1_scores = dict(zip(dm1.ModelID, dm1.DM1_score))
    for src, paths in GDSC_CTRP.items():
        viab = pd.read_csv(paths["viab"], usecols=["CompoundName", "CompoundID"])
        viab = viab.drop_duplicates()
        # find clinical drugs
        viab["nm_l"] = viab.CompoundName.str.lower()
        # Strict: only monotherapies (no ":" in name)
        viab_mono = viab[~viab.CompoundName.str.contains(":", na=False)].copy()
        keep = viab_mono[viab_mono.nm_l.isin(CLINICAL_MAPK)]
        keep = keep.drop_duplicates(["CompoundName", "CompoundID"])
        if not len(keep):
            continue
        log(f"  {src}: clinical MAPK drugs found: "
            f"{keep.CompoundName.tolist()}")

        auc = pd.read_csv(paths["auc"], low_memory=False)
        auc = auc.rename(columns={auc.columns[0]: "ModelID"})
        for _, row in keep.iterrows():
            dpc = row["CompoundID"]; nm = row["CompoundName"]
            if dpc not in auc.columns:
                continue
            sub = auc[["ModelID", dpc]].rename(columns={dpc: "AUC"})
            sub = sub[sub.ModelID.isin(th_ids)].copy()
            sub["DM1"] = sub.ModelID.map(dm1_scores)
            sub = sub.dropna(subset=["AUC", "DM1"])
            if len(sub) < 5:
                rec.append({"source": src, "drug": nm, "compound_id": dpc,
                            "n": len(sub), "spearman_r": None,
                            "spearman_p": None, "cohens_d": None,
                            "mean_auc_high": None, "mean_auc_low": None,
                            "median_auc": None,
                            "note": "n<5 thyroid lines with AUC"})
                continue
            r, p = stats.spearmanr(sub.DM1, sub.AUC)
            med = float(sub.DM1.median())
            high = sub.loc[sub.DM1 > med, "AUC"]
            low = sub.loc[sub.DM1 <= med, "AUC"]
            pv = (high.var(ddof=1) + low.var(ddof=1)) / 2
            d = (high.mean() - low.mean()) / np.sqrt(pv) \
                if pv > 0 and len(high) >= 2 and len(low) >= 2 else None
            t_stat, t_p = (stats.ttest_ind(high, low, equal_var=False)
                           if len(high) >= 3 and len(low) >= 3 else (None, None))
            rec.append({"source": src, "drug": nm, "compound_id": dpc,
                        "n": len(sub), "n_high": len(high), "n_low": len(low),
                        "median_auc": float(sub.AUC.median()),
                        "mean_auc_high": float(high.mean()),
                        "mean_auc_low": float(low.mean()),
                        "spearman_r": float(r), "spearman_p": float(p),
                        "t_stat": float(t_stat) if t_stat is not None else None,
                        "t_p": float(t_p) if t_p is not None else None,
                        "cohens_d": float(d) if d is not None else None,
                        "note": ""})
    out = pd.DataFrame(rec)
    return out


def main() -> None:
    th = load_thyroid_models()
    th_ids = th["ModelID"].tolist()

    gene_map = parse_expr_columns()
    expr, found, missing = load_panel_expression(th_ids, gene_map)
    dm1 = compute_dm1_score(expr, found)

    # Annotate with name + oncotree
    dm1 = dm1.merge(th[["ModelID", "StrippedCellLineName", "OncotreeCode",
                        "OncotreePrimaryDisease"]],
                    on="ModelID", how="left")
    dm1 = dm1.sort_values("DM1_score", ascending=False)
    dm1.to_csv(OUT / "h9_thyroid_dm1_score.tsv", sep="\t", index=False)
    log(f"saved DM1 score n={len(dm1)} thyroid lines (DM1-high top)")

    log("=== CRISPR essentiality ===")
    crispr_out = crispr_essentiality(dm1[["ModelID", "DM1_score"]])
    crispr_out.to_csv(OUT / "h9_depmap_essentiality_top.tsv", sep="\t", index=False)
    log(f"saved CRISPR table: {len(crispr_out):,} genes")

    log("=== PRISM drug response ===")
    prism_out, prism_focus = prism_drug_response(dm1[["ModelID", "DM1_score"]])
    prism_out.to_csv(OUT / "h9_prism_drug_ranking.tsv", sep="\t", index=False)
    prism_focus.to_csv(OUT / "h9_mapk_drug_focus.tsv", sep="\t", index=False)
    log(f"saved PRISM table: {len(prism_out):,} drugs (focus={len(prism_focus)})")

    log("=== GDSC2 / CTRP clinical MAPK monotherapy AUC ===")
    gdsc_out = gdsc_ctrp_clinical(dm1[["ModelID", "DM1_score"]])
    if gdsc_out is not None and len(gdsc_out):
        gdsc_out.to_csv(OUT / "h9_gdsc_ctrp_clinical_mapk.tsv",
                        sep="\t", index=False)
        log(f"saved GDSC/CTRP clinical: {len(gdsc_out)} (drug, source) rows")

    # Tiny summary json for the report
    summary = {
        "n_thyroid_lines": int(len(dm1)),
        "n_dm1_high": int((dm1["DM1_score"] > dm1["DM1_score"].median()).sum()),
        "n_dm1_low":  int((dm1["DM1_score"] <= dm1["DM1_score"].median()).sum()),
        "panel_resolved": found,
        "panel_missing": missing,
        "dm1_score_range": [float(dm1["DM1_score"].min()),
                            float(dm1["DM1_score"].max())],
        "top10_crispr_more_essential_in_dm1_high": (
            crispr_out.head(10)
            [["gene", "spearman_r", "spearman_p", "spearman_fdr", "cohens_d",
              "mean_eff_high", "mean_eff_low"]]
            .to_dict(orient="records")
        ),
        "top10_prism_dm1_high_selective": (
            prism_out.head(10)
            [["broad_id", "name", "moa", "target", "mapk_class",
              "spearman_r", "cohens_d", "t_p", "t_fdr"]]
            .to_dict(orient="records")
        ),
        "mapk_class_focus_summary": {
            cls: {
                "n_drugs": int((prism_focus["mapk_class"] == cls).sum()),
                "median_cohens_d": float(
                    prism_focus.loc[prism_focus["mapk_class"] == cls,
                                    "cohens_d"].median()),
                "frac_negative_d": float(
                    (prism_focus.loc[prism_focus["mapk_class"] == cls,
                                     "cohens_d"] < 0).mean()),
            }
            for cls in ["MEK", "BRAF/RAF", "ERK"]
            if (prism_focus["mapk_class"] == cls).any()
        },
    }
    (OUT / "h9_summary.json").write_text(json.dumps(summary, indent=2,
                                                    default=str))
    log("done.")


if __name__ == "__main__":
    sys.exit(main() or 0)
