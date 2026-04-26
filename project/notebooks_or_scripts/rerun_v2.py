#!/usr/bin/env python3
from __future__ import annotations

import base64
import gzip
import io
import json
import math
import re
import shutil
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib
import warnings

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from scipy.stats import ttest_ind
from sklearn.base import clone
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings(
    "ignore",
    message="This pattern is interpreted as a regular expression, and has match groups",
    category=UserWarning,
)
warnings.filterwarnings(
    "ignore",
    message="'penalty' was deprecated in version 1.8",
    category=FutureWarning,
)

ROOT = Path(__file__).resolve().parents[2]
PROJECT = ROOT / "project"
META = PROJECT / "metadata"
PROCESSED = PROJECT / "data_processed"
RESULTS = PROJECT / "results"
REPORTS = PROJECT / "reports"
LOGS = PROJECT / "logs"
RAW = PROJECT / "data_raw"

LOG_PATH = LOGS / f"rerun_v2_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
LABEL_LOG = LOGS / "label_correction_GSE213647.log"


def log(msg: str) -> None:
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    with LOG_PATH.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def ensure_dirs() -> None:
    for p in [
        RESULTS / "tables",
        RESULTS / "figs",
        RESULTS / "ml",
        REPORTS,
        META,
        LOGS,
    ]:
        p.mkdir(parents=True, exist_ok=True)


def bh_fdr(pvalues: pd.Series) -> pd.Series:
    p = pvalues.fillna(1.0).to_numpy(dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return pd.Series(out, index=pvalues.index)


def normalize_hgnc_symbol(gene: str) -> str:
    return str(gene).strip().upper().replace(" ", "")


TDS16 = [
    "DIO1",
    "DIO2",
    "DUOX1",
    "DUOX2",
    "FOXE1",
    "GLIS3",
    "NKX2-1",
    "PAX8",
    "SLC26A4",
    "SLC5A5",
    "SLC5A8",
    "TG",
    "THRA",
    "THRB",
    "TPO",
    "TSHR",
]

TIERA67_CATEGORIES = {
    "TDS_core": [
        "DIO1",
        "DIO2",
        "DUOX1",
        "DUOX2",
        "FOXE1",
        "GLIS3",
        "NKX2-1",
        "PAX8",
        "SLC26A4",
        "SLC5A5",
        "SLC5A8",
        "TG",
        "THRA",
        "THRB",
        "TPO",
        "TSHR",
    ],
    "MAPK_output_ERK": [
        "DUSP4",
        "DUSP5",
        "DUSP6",
        "SPRY1",
        "SPRY2",
        "SPRY4",
        "ETV4",
        "ETV5",
        "PHLDA1",
        "FOSL1",
    ],
    "Driver_anchor": [
        "BRAF",
        "NRAS",
        "HRAS",
        "KRAS",
        "RET",
        "NTRK1",
        "NTRK3",
        "ALK",
        "PAX8",
        "PPARG",
        "TERT",
        "EIF1AX",
    ],
    "Aggressive_marker": [
        "TP53",
        "CDKN2A",
        "CDKN2B",
        "PIK3CA",
        "AKT1",
        "PTEN",
        "ATM",
        "CTNNB1",
        "APC",
        "MSH2",
    ],
    "Dediff_invasion": [
        "VIM",
        "ZEB1",
        "ZEB2",
        "SNAI1",
        "SNAI2",
        "TWIST1",
        "CDH1",
        "CDH2",
        "MMP9",
        "LOX",
    ],
    "Immune_stromal_light": ["CD274", "CD8A", "FOXP3", "IDO1", "HLA-DRA"],
    "Thyroid_lineage_extra": ["IYD", "THADA", "MET", "KLK10"],
}
TIERA67_ENTRIES = [normalize_hgnc_symbol(g) for genes in TIERA67_CATEGORIES.values() for g in genes]
TIERA67_UNIQUE = list(dict.fromkeys(TIERA67_ENTRIES))

# Driver-free variant: removes the Driver_anchor category whose genes (BRAF/NRAS/KRAS/RET/...) are
# also used to assign BRAF_like vs RAS_like labels in TCGA. Used as a leakage-clean comparator.
TIERA67_CLEAN_CATEGORIES = {k: v for k, v in TIERA67_CATEGORIES.items() if k != "Driver_anchor"}
TIERA67_CLEAN_ENTRIES = [normalize_hgnc_symbol(g) for genes in TIERA67_CLEAN_CATEGORIES.values() for g in genes]
TIERA67_CLEAN_UNIQUE = list(dict.fromkeys(TIERA67_CLEAN_ENTRIES))


def write_fixed_gene_lists(brs71: list[str], brs71_md: str) -> None:
    (META / "tds16_genes.txt").write_text(
        "# source: TCGA PTC Cell 2014 (PMID 25417114), Table S5\n"
        + "\n".join(TDS16)
        + "\n",
        encoding="utf-8",
    )
    tier_lines = ["# hardcoded composite TierA67 entries", "# note: PAX8 appears in two categories by design"]
    for cat, genes in TIERA67_CATEGORIES.items():
        tier_lines.append(f"# [{cat}]")
        tier_lines.extend(genes)
    (META / "tierA67_genes.txt").write_text("\n".join(tier_lines) + "\n", encoding="utf-8")
    (META / "brs71_genes.txt").write_text(
        "# BRS71 proxy derived internally from TCGA-THCA, not original Chakravarty 2011\n"
        + "\n".join(brs71)
        + "\n",
        encoding="utf-8",
    )
    (META / "brs71_derivation.md").write_text(brs71_md.strip() + "\n", encoding="utf-8")


def read_expr(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    gene_col = df.columns[0]
    df = df.rename(columns={gene_col: "gene_symbol"}).set_index("gene_symbol")
    df.index = df.index.map(normalize_hgnc_symbol)
    df = df[~df.index.duplicated(keep="first")]
    return df


def backup_sample_master() -> Path:
    cur = META / "sample_master.tsv"
    backup = META / "sample_master_v1.tsv"
    if cur.exists() and not backup.exists():
        cur.rename(backup)
        log(f"Backed up sample_master.tsv -> {backup.name}")
    elif backup.exists():
        log("sample_master_v1.tsv already exists; keeping existing backup")
    return backup if backup.exists() else cur


def load_sample_master_for_update() -> pd.DataFrame:
    backup = META / "sample_master_v1.tsv"
    src = backup if backup.exists() else META / "sample_master.tsv"
    return pd.read_csv(src, sep="\t")


def parse_tcga_mutation_groups(tcga_samples: set[str]) -> pd.DataFrame:
    mut_dir = RAW / "gdc" / "TCGA-THCA" / "mutation"
    rows = []
    for path in sorted(mut_dir.glob("*.maf.gz")):
        try:
            maf = pd.read_csv(
                path,
                sep="\t",
                comment="#",
                usecols=["Hugo_Symbol", "Tumor_Sample_Barcode", "HGVSp_Short", "Variant_Classification"],
                low_memory=False,
            )
        except Exception as exc:
            log(f"skip maf parse {path.name}: {exc}")
            continue
        maf["sample_id"] = maf["Tumor_Sample_Barcode"].astype(str).str.extract(r"^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-\d{2}[A-Z])")
        maf = maf[maf["sample_id"].isin(tcga_samples)].copy()
        if maf.empty:
            continue
        nonsilent = {
            "Missense_Mutation",
            "Nonsense_Mutation",
            "Frame_Shift_Del",
            "Frame_Shift_Ins",
            "In_Frame_Del",
            "In_Frame_Ins",
            "Splice_Site",
            "Translation_Start_Site",
            "Nonstop_Mutation",
        }
        maf = maf[maf["Variant_Classification"].isin(nonsilent)].copy()
        rows.append(maf)
    if not rows:
        log("No TCGA MAF rows available for mutation anchor derivation")
        return pd.DataFrame(columns=["sample_id", "mutation_group", "driver_anchor"])
    maf = pd.concat(rows, ignore_index=True)
    out = []
    for sample_id, sub in maf.groupby("sample_id"):
        genes = set(sub["Hugo_Symbol"].astype(str))
        hgvsp = set(sub["HGVSp_Short"].astype(str))
        braf_v600e = (
            (sub["Hugo_Symbol"] == "BRAF")
            & (sub["HGVSp_Short"].astype(str).str.contains(r"V(?:600|640)E", regex=True, na=False))
        ).any()
        ras_mut = sub["Hugo_Symbol"].isin(["NRAS", "HRAS", "KRAS"]).any()
        if braf_v600e and not ras_mut:
            group = "BRAF_like"
            driver = "BRAF"
        elif ras_mut and not braf_v600e:
            group = "RAS_like"
            driver = "RAS"
        else:
            group = None
            driver = None
        out.append(
            {
                "sample_id": sample_id,
                "mutation_group": group,
                "driver_anchor": driver,
                "has_braf_v600e": braf_v600e,
                "has_ras_mut": ras_mut,
                "observed_genes": ",".join(sorted(genes & {"BRAF", "NRAS", "HRAS", "KRAS", "RET", "EIF1AX", "TP53", "TERT"})),
            }
        )
    mut_df = pd.DataFrame(out)
    mut_df.to_csv(RESULTS / "tables" / "tcga_thca_mutation_groups.tsv", sep="\t", index=False)
    log(f"TCGA mutation anchor groups available for {mut_df['mutation_group'].notna().sum()} samples")
    return mut_df


def derive_brs71_proxy(tcga_expr: pd.DataFrame, mut_df: pd.DataFrame) -> tuple[list[str], pd.DataFrame, str]:
    eligible = mut_df.dropna(subset=["mutation_group"]).copy()
    common = [s for s in eligible["sample_id"] if s in tcga_expr.columns]
    eligible = eligible[eligible["sample_id"].isin(common)].copy()
    if eligible["mutation_group"].nunique() < 2:
        log("BRS71 proxy derivation skipped: insufficient mutation-anchor classes")
        return [], pd.DataFrame(), "BRS71 proxy derivation failed due to insufficient mutation-anchor classes."

    braf_samples = eligible.loc[eligible["mutation_group"] == "BRAF_like", "sample_id"].tolist()
    ras_samples = eligible.loc[eligible["mutation_group"] == "RAS_like", "sample_id"].tolist()
    Xb = tcga_expr.loc[:, braf_samples].astype(float)
    Xr = tcga_expr.loc[:, ras_samples].astype(float)
    log2fc = Xb.mean(axis=1) - Xr.mean(axis=1)
    stat = ttest_ind(Xb.to_numpy(), Xr.to_numpy(), axis=1, equal_var=False, nan_policy="omit")
    deg = pd.DataFrame(
        {
            "gene_symbol": tcga_expr.index,
            "log2FC_BRAF_vs_RAS": log2fc.to_numpy(),
            "pvalue": stat.pvalue,
        }
    ).dropna(subset=["pvalue"]).sort_values("pvalue")
    deg["FDR"] = bh_fdr(deg["pvalue"])
    sig = deg[(deg["FDR"] < 0.05) & (deg["log2FC_BRAF_vs_RAS"].abs() > 1)].copy()
    sig["rank_score"] = sig["FDR"] * -np.log10(np.clip(sig["pvalue"], 1e-300, 1))
    brs71 = sig.sort_values(["FDR", "pvalue", "log2FC_BRAF_vs_RAS"], ascending=[True, True, False]).head(71)["gene_symbol"].tolist()
    deg.to_csv(RESULTS / "tables" / "brs71_proxy_deg_full.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_symbol": brs71}).to_csv(RESULTS / "tables" / "brs71_proxy_genes.tsv", sep="\t", index=False)
    md = f"""
# BRS71 proxy derivation

- Status: proxy
- Source cohort: TCGA-THCA
- Source expression: `project/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv`
- Mutation anchor source: open GDC masked somatic MAF files available locally
- Positive class: BRAF V600E only
- Negative class: NRAS/HRAS/KRAS mutant
- Exclusions: samples with both anchors, neither anchor, or unavailable MAF
- Differential expression implementation: Welch t-test on log2 normalized expression
- Threshold: abs(log2FC) > 1 and FDR < 0.05
- Selected genes: top 71 by significance among threshold-passing genes
- Important note: this is a BRS71 proxy derived internally from TCGA-THCA, not original Chakravarty 2011
- Mutation-anchor counts used: {dict(Counter(eligible['mutation_group']))}
"""
    log(f"Derived BRS71 proxy with {len(brs71)} genes")
    return brs71, deg, md


def load_gse213647_supplementary() -> pd.DataFrame | None:
    xlsx = META / "GSE213647_supplementary_data.xlsx"
    if not xlsx.exists():
        url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-45366-0/MediaObjects/41467_2024_45366_MOESM5_ESM.xlsx"
        try:
            r = requests.get(url, timeout=300)
            r.raise_for_status()
            xlsx.write_bytes(r.content)
            log("Downloaded GSE213647 supplementary data workbook")
        except Exception as exc:
            log(f"Failed to download GSE213647 supplementary workbook: {exc}")
            return None
    try:
        df = pd.read_excel(xlsx, sheet_name="Supple data13. RNAseq")
        return df
    except Exception as exc:
        log(f"Failed to parse GSE213647 supplementary workbook: {exc}")
        return None


def tissue_to_histology(tissue: str) -> str:
    mapping = {
        "Normal": "normal",
        "PTC": "cPTC",
        "PDTC": "PDTC",
        "ATC": "ATC",
    }
    return mapping.get(str(tissue), "unknown")


def relabel_gse213647(sample_master: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    supp = load_gse213647_supplementary()
    out = sample_master.copy()
    changes = []
    if supp is None:
        mask = out["dataset"] == "GSE213647"
        bad = mask & out["normal_vs_tumor"].eq("normal") & out["histology_subtype"].eq("FTC")
        out.loc[bad, "normal_vs_tumor"] = "unknown_due_to_mislabel"
        out.loc[bad, "label_confidence"] = "medium→excluded_relabeling"
        for _, row in out.loc[bad, ["sample_id", "clinical_subtype_tag"]].iterrows():
            changes.append(f"{row['sample_id']}\t{row['clinical_subtype_tag']}\tfallback_excluded")
        LABEL_LOG.write_text(
            "Supplementary relabel failed; fallback exclusion applied to suspicious normal/FN-labeled samples.\n"
            + "\n".join(changes)
            + "\n",
            encoding="utf-8",
        )
        return out, "excluded"

    supp["sample_id"] = supp["GEO ID"].astype(str)
    supp["histology_subtype_new"] = supp["Tissue type"].map(tissue_to_histology)
    supp["normal_vs_tumor_new"] = np.where(supp["Tissue type"].eq("Normal"), "normal", "tumor")
    supp = supp.set_index("sample_id")

    mask = out["dataset"] == "GSE213647"
    for idx, row in out.loc[mask].iterrows():
        sid = row["sample_id"]
        if sid not in supp.index:
            out.loc[idx, "normal_vs_tumor"] = "unknown_due_to_mislabel"
            out.loc[idx, "label_confidence"] = "medium→excluded_relabeling"
            changes.append(f"{sid}\t{row['histology_subtype']}\tmissing_in_supplementary\texcluded")
            continue
        new_tissue = supp.at[sid, "Tissue type"]
        new_hist = supp.at[sid, "histology_subtype_new"]
        new_nvst = supp.at[sid, "normal_vs_tumor_new"]
        old_hist = row["histology_subtype"]
        old_nvst = row["normal_vs_tumor"]
        out.loc[idx, "normal_vs_tumor"] = new_nvst
        out.loc[idx, "histology_subtype"] = new_hist
        out.loc[idx, "clinical_subtype_tag"] = f"{row['clinical_subtype_tag']} | supplementary_tissue:{new_tissue}"
        out.loc[idx, "source_note"] = f"{row['source_note']} | supplementary_tissue_type={new_tissue}"
        out.loc[idx, "label_confidence"] = "medium→high_corrected"
        changes.append(f"{sid}\t{old_nvst}/{old_hist}\t{new_nvst}/{new_hist}\t{new_tissue}")

    LABEL_LOG.write_text(
        "source: Nature Communications 2024 supplementary data13 RNAseq\n"
        "columns used: GEO ID, Tissue type\n"
        + "\n".join(changes)
        + "\n",
        encoding="utf-8",
    )
    log(f"GSE213647 supplementary relabel applied to {mask.sum()} samples")
    return out, "fixed"


def compute_coverage(dataset_expr: dict[str, pd.DataFrame], brs71: list[str]) -> pd.DataFrame:
    rows = []
    panels = {
        "TDS16": TDS16,
        "BRS71_proxy": brs71,
        "TierA67_entries": TIERA67_ENTRIES,
    }
    for dataset, expr in dataset_expr.items():
        genes = set(expr.index.map(normalize_hgnc_symbol))
        row = {"dataset": dataset, "total_gene_count": len(genes)}
        for panel_name, panel in panels.items():
            present = [g for g in panel if normalize_hgnc_symbol(g) in genes]
            row[f"{panel_name}_coverage_n"] = len(present)
            row[f"{panel_name}_coverage_frac"] = len(present) / len(panel) if panel else np.nan
            row[f"{panel_name}_covered_genes"] = ",".join(sorted(set(present)))
        rows.append(row)
    cov = pd.DataFrame(rows)
    cov.to_csv(RESULTS / "tables" / "gene_coverage_table.tsv", sep="\t", index=False)
    cov.to_csv(META / "gene_coverage_table.tsv", sep="\t", index=False)
    return cov


def get_tcga_training_labels(sample_master: pd.DataFrame, mut_df: pd.DataFrame, tcga_expr: pd.DataFrame) -> pd.DataFrame:
    tcga_meta = sample_master[sample_master["dataset"] == "TCGA-THCA"].copy()
    tcga_meta = tcga_meta.merge(mut_df[["sample_id", "mutation_group", "driver_anchor"]], on="sample_id", how="left")
    tcga_meta = tcga_meta[tcga_meta["sample_id"].isin(tcga_expr.columns)].copy()
    tcga_meta = tcga_meta[tcga_meta["mutation_group"].isin(["BRAF_like", "RAS_like"])].copy()
    tcga_meta["label_bin"] = np.where(tcga_meta["mutation_group"] == "BRAF_like", 1, 0)
    return tcga_meta


def classify_external_datasets(sample_master: pd.DataFrame) -> dict[str, pd.DataFrame]:
    candidates = {}

    g27155 = sample_master[sample_master["dataset"] == "GSE27155"].copy()
    labels = []
    reasons = []
    for _, row in g27155.iterrows():
        note = json.loads(row["source_note"])
        tissue = note.get("tissue", "")
        morph = note.get("morphology of papillary carcinomas", "")
        if tissue == "Medullary Thyroid Carcinoma":
            labels.append(np.nan)
            reasons.append("exclude_MTC")
        elif tissue == "Papillary Thyroid Carcinoma" and morph in {
            "classical type of papillary carcinoma",
            "tall cell variant of papillary carcinoma",
        }:
            labels.append("BRAF_like")
            reasons.append("histology_proxy")
        elif tissue == "Papillary Thyroid Carcinoma" and morph == "follicular variant of papillary carcinoma":
            labels.append("RAS_like")
            reasons.append("histology_proxy")
        elif tissue in {"Follicular Thyroid Carcinoma", "Oncocytic Thyroid Carcinoma"}:
            labels.append("RAS_like")
            reasons.append("histology_proxy")
        else:
            labels.append(np.nan)
            reasons.append("exclude_noncomparable")
    g27155["external_label"] = labels
    g27155["external_reason"] = reasons
    candidates["GSE27155"] = g27155

    g126698 = sample_master[sample_master["dataset"] == "GSE126698"].copy()
    g126698["expr_sample_id"] = g126698["clinical_subtype_tag"].str.replace(" [RNA-Seq]", "", regex=False)
    g126698["external_label"] = pd.Series(index=g126698.index, dtype="object")
    g126698.loc[g126698["histology_subtype"].eq("cPTC"), "external_label"] = "BRAF_like"
    g126698.loc[g126698["histology_subtype"].eq("FTC"), "external_label"] = "RAS_like"
    g126698["external_reason"] = np.where(g126698["external_label"].notna(), "histology_proxy", "exclude_noncomparable")
    candidates["GSE126698"] = g126698

    g213647 = sample_master[sample_master["dataset"] == "GSE213647"].copy()
    g213647["external_label"] = np.nan
    g213647["external_reason"] = "exclude_no_ras_like_tumor_class_after_correction"
    candidates["GSE213647"] = g213647

    return candidates


def select_feature_genes(feature_set: str, X_train: pd.DataFrame, brs71: list[str]) -> list[str]:
    if feature_set == "TDS16":
        genes = TDS16
    elif feature_set == "TierA67":
        genes = TIERA67_UNIQUE
    elif feature_set == "TierA67_clean":
        genes = TIERA67_CLEAN_UNIQUE
    elif feature_set == "BRS71_proxy":
        genes = brs71
    elif feature_set == "variance_top50":
        return X_train.var(axis=0).sort_values(ascending=False).head(50).index.tolist()
    else:
        raise ValueError(feature_set)
    return [g for g in genes if g in X_train.columns]


def bootstrap_auc_ci(y_true, scores, n_boot: int = 1000, seed: int = 1) -> tuple[float, float]:
    """Percentile bootstrap 95% CI for ROC-AUC. Returns (lo, hi) or (nan, nan)."""
    y_true = np.asarray(y_true)
    scores = np.asarray(scores, dtype=float)
    n = len(y_true)
    if n == 0 or len(np.unique(y_true)) < 2:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        yb = y_true[idx]
        if len(np.unique(yb)) < 2:
            continue
        vals.append(roc_auc_score(yb, scores[idx]))
    if not vals:
        return np.nan, np.nan
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def youden_threshold(y_true, scores) -> float:
    """Threshold that maximizes TPR - FPR (Youden's J). Defaults to 0.5 when undefined."""
    y_true = np.asarray(y_true)
    scores = np.asarray(scores, dtype=float)
    if len(np.unique(y_true)) < 2:
        return 0.5
    fpr, tpr, thr = roc_curve(y_true, scores)
    j = tpr - fpr
    idx = int(np.argmax(j))
    return float(thr[idx]) if np.isfinite(thr[idx]) else 0.5


def threshold_curve_payload(y_true, scores, n_points: int = 101) -> dict:
    """Per-threshold precision/recall/balanced_accuracy/F1 for interactive threshold slider."""
    y_true = np.asarray(y_true, dtype=int)
    scores = np.asarray(scores, dtype=float)
    if len(scores) == 0 or len(np.unique(y_true)) < 2:
        return {"thresholds": [], "precision": [], "recall": [], "f1": [], "balanced_accuracy": []}
    thrs = np.linspace(np.min(scores), np.max(scores), n_points)
    prec, rec, f1s, bacc = [], [], [], []
    for t in thrs:
        pred = (scores >= t).astype(int)
        tp = int(((pred == 1) & (y_true == 1)).sum())
        fp = int(((pred == 1) & (y_true == 0)).sum())
        fn = int(((pred == 0) & (y_true == 1)).sum())
        tn = int(((pred == 0) & (y_true == 0)).sum())
        p = tp / (tp + fp) if (tp + fp) else 0.0
        r = tp / (tp + fn) if (tp + fn) else 0.0
        f = 2 * p * r / (p + r) if (p + r) else 0.0
        tpr = r
        tnr = tn / (tn + fp) if (tn + fp) else 0.0
        prec.append(p)
        rec.append(r)
        f1s.append(f)
        bacc.append((tpr + tnr) / 2.0)
    return {
        "thresholds": [float(x) for x in thrs],
        "precision": prec,
        "recall": rec,
        "f1": f1s,
        "balanced_accuracy": bacc,
    }


def model_score(estimator, X):
    if hasattr(estimator, "predict_proba"):
        return estimator.predict_proba(X)[:, 1]
    if hasattr(estimator, "decision_function"):
        return estimator.decision_function(X)
    pred = estimator.predict(X)
    return pred.astype(float)


def make_roc_html(y_true: np.ndarray, scores: np.ndarray, title: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    html = [f"<html><body><h2>{title}</h2>"]
    if len(np.unique(y_true)) < 2:
        html.append("<p>ROC unavailable: only one class present in evaluation set.</p></body></html>")
        out_path.write_text("".join(html), encoding="utf-8")
        return
    fpr, tpr, _ = roc_curve(y_true, scores)
    auc = roc_auc_score(y_true, scores)
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.plot(fpr, tpr, label=f"AUC={auc:.3f}")
    ax.plot([0, 1], [0, 1], linestyle="--", color="grey")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title(title)
    ax.legend(loc="lower right")
    buf = io.BytesIO()
    fig.tight_layout()
    fig.savefig(buf, format="png", dpi=180)
    plt.close(fig)
    b64 = base64.b64encode(buf.getvalue()).decode("ascii")
    html.append(f'<img src="data:image/png;base64,{b64}"/>')
    html.append("</body></html>")
    out_path.write_text("".join(html), encoding="utf-8")


def run_ml_v2(sample_master: pd.DataFrame, tcga_expr: pd.DataFrame, brs71: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    candidates = classify_external_datasets(sample_master)
    train_meta = get_tcga_training_labels(sample_master, parse_tcga_mutation_groups(set(tcga_expr.columns)), tcga_expr)
    X_all = tcga_expr[train_meta["sample_id"].tolist()].T
    X_all.columns = X_all.columns.map(normalize_hgnc_symbol)
    y_all = train_meta.set_index("sample_id").loc[X_all.index, "label_bin"].astype(int)
    log(
        f"ML training matrix prepared: {X_all.shape[0]} samples x {X_all.shape[1]} genes; "
        f"class counts={dict(Counter(y_all))}"
    )

    model_defs = {
        "LogReg_l2": LogisticRegression(max_iter=5000, class_weight="balanced"),
        "LogReg_elasticnet": LogisticRegression(
            max_iter=5000, class_weight="balanced", solver="saga", penalty="elasticnet", l1_ratio=0.5
        ),
        "RandomForest": RandomForestClassifier(n_estimators=500, random_state=1, class_weight="balanced"),
        "GradientBoosting": GradientBoostingClassifier(random_state=1),
    }
    try:
        from xgboost import XGBClassifier

        model_defs["XGBoost"] = XGBClassifier(
            n_estimators=400,
            max_depth=4,
            learning_rate=0.05,
            subsample=0.8,
            colsample_bytree=0.8,
            eval_metric="logloss",
            random_state=1,
        )
        log("xgboost available; included in model set")
    except Exception as exc:
        log(f"xgboost unavailable; skip ({exc})")

    internal_rows = []
    external_rows = []
    best_external = {"dataset": "none", "auc": np.nan, "feature_set": None, "model": None}
    # Interactive-dashboard payload: per-(feature_set, model, dataset) threshold curves + bootstrap CI
    curves_payload: dict[str, dict] = {}
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=1)

    for feature_set in ["TDS16", "TierA67", "TierA67_clean", "variance_top50"]:
        for model_name, base_model in model_defs.items():
            log(f"ML start: feature_set={feature_set}, model={model_name}")
            scores = pd.Series(index=X_all.index, dtype=float)
            preds = pd.Series(index=X_all.index, dtype=float)
            per_fold_genes = []
            for fold_i, (tr, te) in enumerate(cv.split(X_all, y_all), start=1):
                Xtr = X_all.iloc[tr]
                Xte = X_all.iloc[te]
                ytr = y_all.iloc[tr]
                genes = select_feature_genes(feature_set, Xtr, brs71)
                per_fold_genes.append(len(genes))
                if len(genes) < 3:
                    continue
                pipe = Pipeline(
                    [
                        ("imputer", SimpleImputer(strategy="median")),
                        ("scaler", StandardScaler()),
                        ("model", clone(base_model)),
                    ]
                )
                pipe.fit(Xtr[genes], ytr)
                fold_scores = model_score(pipe, Xte[genes])
                scores.iloc[te] = fold_scores
                preds.iloc[te] = (fold_scores >= 0.5).astype(int)

            valid = scores.dropna().index
            y_valid = y_all.loc[valid]
            score_valid = scores.loc[valid]
            pred_valid = preds.loc[valid].astype(int)
            two_class = len(np.unique(y_valid)) == 2
            auc = roc_auc_score(y_valid, score_valid) if two_class else np.nan
            pr_auc = average_precision_score(y_valid, score_valid) if two_class else np.nan
            ci_lo, ci_hi = bootstrap_auc_ci(y_valid.to_numpy(), score_valid.to_numpy())
            youden_thr = youden_threshold(y_valid.to_numpy(), score_valid.to_numpy())
            # Reclassify using Youden threshold for a more honest balanced-accuracy
            pred_youden = (score_valid.to_numpy() >= youden_thr).astype(int)
            bacc = balanced_accuracy_score(y_valid, pred_valid) if len(valid) else np.nan
            bacc_youden = balanced_accuracy_score(y_valid, pred_youden) if len(valid) else np.nan
            f1 = f1_score(y_valid, pred_valid) if two_class else np.nan
            mcc = matthews_corrcoef(y_valid, pred_valid) if len(valid) else np.nan
            brier = brier_score_loss(y_valid, score_valid) if two_class and np.all((score_valid >= 0) & (score_valid <= 1)) else np.nan
            cm = confusion_matrix(y_valid, pred_valid, labels=[1, 0])
            pd.DataFrame(cm, index=["true_BRAF_like", "true_RAS_like"], columns=["pred_BRAF_like", "pred_RAS_like"]).to_csv(
                RESULTS / "ml" / f"confusion_v2_{feature_set}_{model_name}.tsv", sep="\t"
            )
            internal_rows.append(
                {
                    "task": "BRAF_like_vs_RAS_like",
                    "dataset": "TCGA-THCA",
                    "feature_set": feature_set,
                    "model": model_name,
                    "n_samples": len(valid),
                    "cv_auc": auc,
                    "cv_pr_auc": pr_auc,
                    "cv_auc_ci_lo": ci_lo,
                    "cv_auc_ci_hi": ci_hi,
                    "cv_brier": brier,
                    "cv_balanced_accuracy": bacc,
                    "cv_balanced_accuracy_youden": bacc_youden,
                    "cv_f1": f1,
                    "cv_mcc": mcc,
                    "cv_youden_threshold": youden_thr,
                    "median_feature_count": int(np.median(per_fold_genes)) if per_fold_genes else 0,
                    "status": "ok" if len(valid) else "skipped",
                }
            )
            # Interactive dashboard data (ROC + PR + per-threshold curves)
            if two_class:
                fpr, tpr, _ = roc_curve(y_valid, score_valid)
                prec_arr, rec_arr, _ = precision_recall_curve(y_valid, score_valid)
                curves_payload.setdefault(feature_set, {}).setdefault(model_name, {})["TCGA-THCA_internal"] = {
                    "n_samples": int(len(valid)),
                    "auc": float(auc),
                    "pr_auc": float(pr_auc),
                    "auc_ci_lo": float(ci_lo) if not pd.isna(ci_lo) else None,
                    "auc_ci_hi": float(ci_hi) if not pd.isna(ci_hi) else None,
                    "youden_threshold": float(youden_thr),
                    "brier": float(brier) if not pd.isna(brier) else None,
                    "roc": {"fpr": [float(x) for x in fpr], "tpr": [float(x) for x in tpr]},
                    "pr": {"precision": [float(x) for x in prec_arr], "recall": [float(x) for x in rec_arr]},
                    "scores": [float(x) for x in score_valid.to_numpy()],
                    "y_true": [int(x) for x in y_valid.to_numpy()],
                    "threshold_curve": threshold_curve_payload(y_valid.to_numpy(), score_valid.to_numpy()),
                }
            log(
                f"ML internal done: feature_set={feature_set}, model={model_name}, "
                f"auc={auc:.3f} (95% CI {ci_lo:.3f}-{ci_hi:.3f}) pr_auc={pr_auc:.3f} "
                f"bacc@0.5={bacc:.3f} bacc@youden={bacc_youden:.3f} mcc={mcc:.3f}"
            )

            full_genes = select_feature_genes(feature_set, X_all, brs71)
            if len(full_genes) >= 3:
                full_pipe = Pipeline(
                    [
                        ("imputer", SimpleImputer(strategy="median")),
                        ("scaler", StandardScaler()),
                        ("model", clone(base_model)),
                    ]
                )
                full_pipe.fit(X_all[full_genes], y_all)
            else:
                full_pipe = None

            for dataset, meta in candidates.items():
                use = meta[meta["external_label"].isin(["BRAF_like", "RAS_like"])].copy()
                expr_sample_col = "sample_id"
                if dataset == "GSE126698":
                    expr_sample_col = "expr_sample_id"
                    expr = read_expr(PROCESSED / "bulk_rnaseq" / "GSE126698_rnaseq_expression_log2.tsv")
                elif dataset == "GSE27155":
                    expr = read_expr(PROCESSED / "microarray" / "GSE27155_microarray_expression_log2.tsv")
                elif dataset == "GSE213647":
                    expr = read_expr(PROCESSED / "bulk_rnaseq" / "GSE213647_rnaseq_expression_log2.tsv")
                else:
                    continue
                use = use[use[expr_sample_col].isin(expr.columns)].copy()
                if use.empty or full_pipe is None:
                    external_rows.append(
                        {
                            "task": "BRAF_like_vs_RAS_like",
                            "dataset": dataset,
                            "feature_set": feature_set,
                            "model": model_name,
                            "n_samples": 0,
                            "auc": np.nan,
                            "balanced_accuracy": np.nan,
                            "f1": np.nan,
                            "mcc": np.nan,
                            "status": "skipped_no_eligible_samples_or_features",
                        }
                    )
                    continue
                Xext = expr[use[expr_sample_col].tolist()].T
                Xext.columns = Xext.columns.map(normalize_hgnc_symbol)
                Xext = Xext.reindex(columns=full_genes)
                yext = (use.set_index(expr_sample_col).loc[Xext.index, "external_label"] == "BRAF_like").astype(int)
                ext_scores = model_score(full_pipe, Xext)
                ext_pred_05 = (ext_scores >= 0.5).astype(int)
                # Apply TCGA-CV-derived Youden threshold to the external cohort (domain-honest)
                ext_pred_youden = (np.asarray(ext_scores) >= youden_thr).astype(int)
                ext_two_class = len(np.unique(yext)) == 2
                ext_auc = roc_auc_score(yext, ext_scores) if ext_two_class else np.nan
                ext_pr_auc = average_precision_score(yext, ext_scores) if ext_two_class else np.nan
                ext_ci_lo, ext_ci_hi = bootstrap_auc_ci(yext.to_numpy(), np.asarray(ext_scores))
                ext_bacc = balanced_accuracy_score(yext, ext_pred_05) if len(yext) else np.nan
                ext_bacc_youden = balanced_accuracy_score(yext, ext_pred_youden) if len(yext) else np.nan
                ext_f1 = f1_score(yext, ext_pred_05) if ext_two_class else np.nan
                ext_f1_youden = f1_score(yext, ext_pred_youden) if ext_two_class else np.nan
                ext_mcc = matthews_corrcoef(yext, ext_pred_05) if len(yext) else np.nan
                external_rows.append(
                    {
                        "task": "BRAF_like_vs_RAS_like",
                        "dataset": dataset,
                        "feature_set": feature_set,
                        "model": model_name,
                        "n_samples": len(yext),
                        "auc": ext_auc,
                        "pr_auc": ext_pr_auc,
                        "auc_ci_lo": ext_ci_lo,
                        "auc_ci_hi": ext_ci_hi,
                        "balanced_accuracy": ext_bacc,
                        "balanced_accuracy_youden": ext_bacc_youden,
                        "f1": ext_f1,
                        "f1_youden": ext_f1_youden,
                        "mcc": ext_mcc,
                        "threshold_applied_youden": float(youden_thr),
                        "status": "ok" if len(yext) else "skipped",
                    }
                )
                if ext_two_class:
                    fpr_e, tpr_e, _ = roc_curve(yext, ext_scores)
                    prec_e, rec_e, _ = precision_recall_curve(yext, ext_scores)
                    curves_payload.setdefault(feature_set, {}).setdefault(model_name, {})[dataset] = {
                        "n_samples": int(len(yext)),
                        "auc": float(ext_auc),
                        "pr_auc": float(ext_pr_auc),
                        "auc_ci_lo": float(ext_ci_lo) if not pd.isna(ext_ci_lo) else None,
                        "auc_ci_hi": float(ext_ci_hi) if not pd.isna(ext_ci_hi) else None,
                        "youden_threshold_from_tcga": float(youden_thr),
                        "roc": {"fpr": [float(x) for x in fpr_e], "tpr": [float(x) for x in tpr_e]},
                        "pr": {"precision": [float(x) for x in prec_e], "recall": [float(x) for x in rec_e]},
                        "scores": [float(x) for x in np.asarray(ext_scores)],
                        "y_true": [int(x) for x in yext.to_numpy()],
                        "threshold_curve": threshold_curve_payload(yext.to_numpy(), np.asarray(ext_scores)),
                    }
                log(
                    f"ML external done: dataset={dataset}, feature_set={feature_set}, model={model_name}, "
                    f"n={len(yext)} auc={ext_auc if not pd.isna(ext_auc) else 'nan'}"
                )
                if not pd.isna(ext_auc) and (pd.isna(best_external["auc"]) or ext_auc > best_external["auc"]):
                    best_external = {
                        "dataset": dataset,
                        "auc": float(ext_auc),
                        "feature_set": feature_set,
                        "model": model_name,
                    }
                    make_roc_html(
                        yext.to_numpy(),
                        np.asarray(ext_scores, dtype=float),
                        f"BRAF_like vs RAS_like ROC - {dataset}",
                        RESULTS / "figs" / f"roc_braf_like_vs_ras_like_{dataset}.html",
                    )

            make_roc_html(
                y_valid.to_numpy(),
                score_valid.to_numpy(dtype=float),
                f"BRAF_like vs RAS_like ROC - TCGA-THCA CV - {feature_set} - {model_name}",
                RESULTS / "figs" / f"roc_braf_like_vs_ras_like_TCGA-THCA_{feature_set}_{model_name}.html",
            )

    internal_df = pd.DataFrame(internal_rows)
    external_df = pd.DataFrame(external_rows)
    internal_df.to_csv(RESULTS / "ml" / "baseline_ml_results.tsv", sep="\t", index=False)
    external_df.to_csv(RESULTS / "ml" / "baseline_ml_external.tsv", sep="\t", index=False)
    # Write compact interactive-curves JSON for the v2 interactive dashboard.
    # Uses a DIFFERENT filename than the legacy model_curves.json produced by
    # build_html_reports.py so both can coexist without either overwriting the other.
    curves_out_primary = RESULTS / "ml" / "model_curves_v2.json"
    curves_out_dashboard = REPORTS / "html" / "assets" / "data" / "model_curves_v2.json"
    curves_out_dashboard.parent.mkdir(parents=True, exist_ok=True)
    payload_text = json.dumps(curves_payload, separators=(",", ":"))
    curves_out_primary.write_text(payload_text, encoding="utf-8")
    curves_out_dashboard.write_text(payload_text, encoding="utf-8")
    log(f"model_curves_v2.json written ({len(payload_text)/1024:.1f} KB)")
    best_msg = (
        f"{best_external['dataset']} AUC={best_external['auc']:.3f}"
        if not pd.isna(best_external["auc"])
        else "none"
    )
    return internal_df, external_df, best_msg


def build_ml_report_v2(internal_df: pd.DataFrame, external_df: pd.DataFrame, brs71: list[str], gse213647_status: str) -> None:
    best_cv = internal_df.sort_values("cv_auc", ascending=False).head(10)
    best_ext = external_df.sort_values("auc", ascending=False, na_position="last").head(10)
    report = f"""
# ml_report_v2

## Setup

- Training cohort: TCGA-THCA
- Task: BRAF_like vs RAS_like
- Label source: open GDC masked somatic MAF anchor
- BRAF_like definition: BRAF V600E only
- RAS_like definition: KRAS/HRAS/NRAS mutant
- Feature sets: TDS16, TierA67, variance_top50
- BRS71 handling: proxy list derived internally from TCGA-THCA ({len(brs71)} genes); not used as a direct model feature in v2 because the required feature sets were fixed by task
- GSE213647 relabel status: {gse213647_status}

## Internal CV

```tsv
{best_cv.to_csv(sep='\t', index=False).strip()}
```

## External validation

```tsv
{best_ext.to_csv(sep='\t', index=False).strip()}
```

## Interpretation

- External validation based on GSE27155 is still label-inferred, not mutation-verified.
- GSE213647 is no longer treated as a valid BRAF_like vs RAS_like external set because the supplementary sample table resolves it as PTC/PDTC/ATC/Normal without a comparable RAS-like differentiated tumor class.
- GSE126698 can be used only as a very small exploratory external set after restricting to PTC vs FTC.
"""
    (REPORTS / "ml_report_v2.md").write_text(report.strip() + "\n", encoding="utf-8")


def write_analysis_summary_v2(
    sample_master: pd.DataFrame,
    coverage: pd.DataFrame,
    internal_df: pd.DataFrame,
    external_df: pd.DataFrame,
    brs71: list[str],
    gse213647_status: str,
) -> None:
    best_ext = external_df.sort_values("auc", ascending=False, na_position="last").head(1)
    best_ext_line = (
        f"{best_ext.iloc[0]['dataset']} AUC={best_ext.iloc[0]['auc']:.3f}"
        if not best_ext.empty and pd.notna(best_ext.iloc[0]["auc"])
        else "none"
    )
    limitations = [
        "BRS71 is a proxy derived from TCGA-THCA, not the original published list.",
        "TierA67 hardcoded composite contains 67 category entries but 66 unique HGNC-normalized symbols because PAX8 appears in two categories.",
        "GSE27155 external labels are phenotype proxies inferred from histology/morphology, not direct molecular subtypes.",
        "GSE126698 external subset is very small after restricting to PTC vs FTC.",
        "TCGA mutation anchors depend on locally available open MAFs, not a controlled-access harmonized mutation callset.",
    ]
    txt = f"""
# analysis_summary_v2

## What changed from v1

- TDS16 is now fixed to the official Cell 2014 list.
- BRS71 is now handled explicitly as a proxy derived internally from TCGA-THCA.
- TierA67 is now a hardcoded composite panel with category annotations.
- GSE213647 labels were revisited against the Nature Communications 2024 supplementary RNA cohort sheet.
- The ML task was redefined to TCGA mutation-anchor BRAF_like vs RAS_like instead of the previous weak histology-only internal task.

## Gene panels

- TDS16 confirmed from Cell 2014 (PMID 25417114).
- BRS71 proxy gene count: {len(brs71)}.
- TierA67 hardcoded entries: 67.

## GSE213647 relabeling

- Status: {gse213647_status}.
- Evidence source: supplementary RNA cohort sheet from the 2024 Nature Communications SHMT2 paper.
- Outcome: the cohort resolves to Normal/PTC/PDTC/ATC. It does not provide a directly comparable RAS-like differentiated tumor class for 2-class external validation.

## Coverage

```tsv
{coverage.to_csv(sep='\t', index=False).strip()}
```

## ML performance

- Best external validation candidate: {best_ext_line}

### Internal CV

```tsv
{internal_df.sort_values('cv_auc', ascending=False).head(10).to_csv(sep='\t', index=False).strip()}
```

### External validation

```tsv
{external_df.sort_values('auc', ascending=False, na_position='last').head(10).to_csv(sep='\t', index=False).strip()}
```

## Remaining limitations

""" + "\n".join(f"- {x}" for x in limitations)
    (REPORTS / "analysis_summary_v2.md").write_text(txt.strip() + "\n", encoding="utf-8")


def write_next_steps_v2() -> None:
    txt = """
# next_steps_v2

1. Replace the BRS71 proxy with the original source-verified supplementary list once the exact official table is pinned down.
2. Add a better external RNA-seq validation cohort with explicit follicular-patterned tumor labels rather than generic PTC-only labels.
3. Reprocess Affymetrix CEL files in R/Bioconductor for a cleaner microarray validation layer.
4. Expand TCGA anchors beyond BRAF V600E and RAS to fusion-aware subtype definitions when a consistent public callset is available.
5. Prepare a small Phase 1 pilot model package using only TDS16 and TierA67 on clearly labeled differentiated thyroid tumors.
"""
    (REPORTS / "next_steps_v2.md").write_text(txt.strip() + "\n", encoding="utf-8")


def main() -> None:
    ensure_dirs()
    backup_sample_master()
    sample_master = load_sample_master_for_update()

    tcga_expr = read_expr(PROCESSED / "bulk_rnaseq" / "TCGA-THCA_rnaseq_expression_log2.tsv")
    g213647_expr = read_expr(PROCESSED / "bulk_rnaseq" / "GSE213647_rnaseq_expression_log2.tsv")
    g126698_expr = read_expr(PROCESSED / "bulk_rnaseq" / "GSE126698_rnaseq_expression_log2.tsv")
    g27155_expr = read_expr(PROCESSED / "microarray" / "GSE27155_microarray_expression_log2.tsv")
    g76039_expr = read_expr(PROCESSED / "microarray" / "GSE76039_microarray_expression_log2.tsv")

    mut_df = parse_tcga_mutation_groups(set(tcga_expr.columns))
    brs71, brs_deg, brs_md = derive_brs71_proxy(tcga_expr, mut_df)
    write_fixed_gene_lists(brs71, brs_md)
    log(f"Gene panels written: TDS16={len(TDS16)}, BRS71_proxy={len(brs71)}, TierA67_entries={len(TIERA67_ENTRIES)}")

    sample_master, gse213647_status = relabel_gse213647(sample_master)
    sample_master["driver_anchor"] = sample_master["driver_anchor"].astype("object")
    sample_master["molecular_subtype"] = sample_master["molecular_subtype"].astype("object")
    tcga_updates = mut_df.dropna(subset=["mutation_group"]).set_index("sample_id")
    tcga_mask = sample_master["dataset"].eq("TCGA-THCA") & sample_master["sample_id"].isin(tcga_updates.index)
    sample_master.loc[tcga_mask, "driver_anchor"] = sample_master.loc[tcga_mask, "sample_id"].map(tcga_updates["driver_anchor"])
    sample_master.loc[tcga_mask, "molecular_subtype"] = sample_master.loc[tcga_mask, "sample_id"].map(tcga_updates["mutation_group"])

    sample_master.to_csv(META / "sample_master.tsv", sep="\t", index=False)
    log("Updated sample_master.tsv")

    coverage = compute_coverage(
        {
            "TCGA-THCA": tcga_expr,
            "GSE213647": g213647_expr,
            "GSE126698": g126698_expr,
            "GSE27155": g27155_expr,
            "GSE76039": g76039_expr,
        },
        brs71,
    )
    log("Gene coverage table updated")
    internal_df, external_df, best_ext_line = run_ml_v2(sample_master, tcga_expr, brs71)
    build_ml_report_v2(internal_df, external_df, brs71, gse213647_status)
    write_analysis_summary_v2(sample_master, coverage, internal_df, external_df, brs71, gse213647_status)
    write_next_steps_v2()
    log("v2 reports written")

    crit = [
        "BRS71 remains proxy, not the original published list",
        "GSE27155 external labels are inferred proxies",
        "GSE126698 external subset is small",
    ]
    print("=== RERUN COMPLETE ===")
    print("TDS16 genes: confirmed (source: Cell 2014)")
    print(f"BRS71 genes: {len(brs71)} (proxy derived from TCGA-THCA)")
    print("TierA67 genes: 67 (hardcoded composite)")
    print(f"GSE213647 labels: {gse213647_status}")
    print(f"Best external validation: {best_ext_line}")
    print(f"Critical limitations still present: {', '.join(crit)}")
    print("Reports: project/reports/analysis_summary_v2.md")
    print("=====================")


if __name__ == "__main__":
    main()
