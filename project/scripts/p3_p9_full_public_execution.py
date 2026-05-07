#!/usr/bin/env python3
"""Full-public Paper 3/Paper 9 execution.

This script intentionally writes only into project/results/p3_p9_full_execution
and project/reports/p3_p9_full_execution. It does not edit manuscript prose.
"""

from __future__ import annotations

import gzip
import json
import math
import re
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from scipy.stats import pearsonr, spearmanr
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[2]
RESULT_ROOT = ROOT / "project/results/p3_p9_full_execution"
REPORT_ROOT = ROOT / "project/reports/p3_p9_full_execution"
P3 = RESULT_ROOT / "paper3"
P9 = RESULT_ROOT / "paper9"
RAW = P9 / "raw"
SHARED = RESULT_ROOT / "shared"
MAX_SINGLE_DOWNLOAD = 10 * 1024**3
MAX_TOTAL_DOWNLOAD = 50 * 1024**3
DEPMAP_API = "https://depmap.org/portal/data_page/api/data"
DEPMAP_BASE = "https://depmap.org"


SIGNATURES: dict[str, list[str]] = {
    "IFNG_TIS": ["IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA", "GZMB", "PRF1"],
    "CYTOLYTIC": ["GZMA", "GZMB", "PRF1", "GNLY", "NKG7"],
    "MHC1_CORE": ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "NLRC5"],
    "MHC2_CORE": ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "CIITA"],
    "ANTIGEN_PROCESSING": ["TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "PSMB10", "ERAP1", "ERAP2"],
    "TLS_CABRITA": ["CXCL13", "CCL19", "CCL21", "LTB", "CD79A", "MS4A1", "BANK1"],
    "CXCL13_AXIS": ["CXCL13", "CXCR5", "ICOS", "PDCD1", "TOX", "BCL6"],
    "T_CELL_EXHAUSTION": ["PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX", "CXCL13"],
    "TREG": ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18"],
    "M2_TAM": ["CD163", "MRC1", "MSR1", "IL10", "VSIG4", "MARCO"],
    "MDSC_LIKE": ["S100A8", "S100A9", "ARG1", "IL1B", "CXCR2", "FCGR3B"],
    "STROMAL_EXCLUSION": ["TGFB1", "TGFBI", "COL1A1", "COL3A1", "ACTA2", "FAP", "CXCL12"],
    "IFN_GAMMA_RESPONSE": ["STAT1", "IRF1", "GBP1", "GBP2", "CXCL9", "CXCL10", "IDO1"],
    "DEDIFFERENTIATED_EPI": ["VIM", "FN1", "SNAI2", "TWIST1", "ZEB1", "ITGA5", "MMP2"],
}

LINEAGE_SIGNATURES: dict[str, list[str]] = {
    "RAI_8_score": ["SLC5A5", "TPO", "TG", "TSHR", "DUOX2", "DUOXA2", "DIO1", "SLC26A7"],
    "THYROID_NONOVERLAP_score": ["TG", "TPO", "TSHR", "SLC5A5", "IYD", "SLC26A4", "DUOX2", "DUOXA2"],
    "TDS_like_score": ["TG", "TPO", "TSHR", "SLC5A5", "DIO1", "DIO2", "FOXE1", "PAX8", "NKX2-1"],
    "TF_collapse_score": ["PAX8", "NKX2-1", "FOXE1", "HHEX"],
    "STAT3_AP1_DNMT_score": ["STAT3", "FOS", "JUN", "JUNB", "FOSL1", "DNMT1", "DNMT3B"],
}

CANDIDATES: dict[str, str] = {
    "NAMPT": "NAD_salvage",
    "STAT3": "JAK_STAT",
    "JAK1": "JAK_STAT",
    "JAK2": "JAK_STAT",
    "IL6R": "JAK_STAT",
    "OSMR": "JAK_STAT",
    "DNMT1": "DNMT",
    "DNMT3B": "DNMT",
    "LYN": "SRC_family",
    "SRC": "SRC_family",
    "FYN": "SRC_family",
    "KCNN4": "ion_channel",
    "ATR": "DDR",
    "CHEK2": "DDR",
    "MYC": "MYC",
    "GLS": "metabolism",
    "LDHA": "metabolism",
    "BCL2L1": "apoptosis",
    "MCL1": "apoptosis",
    "BCL2": "apoptosis",
    "AURKA": "cell_cycle",
    "CDK7": "transcription_cycle",
    "CDK9": "transcription_cycle",
    "HDAC1": "chromatin",
    "HDAC2": "chromatin",
    "EZH2": "chromatin",
}

VERDICT_PRIORITY = {"HIGH_PRIORITY": 0, "MEDIUM": 1, "NEEDS_WETLAB": 2, "LOW": 3, "DROP": 4}

DRUG_REGEX: dict[str, str] = {
    "NAMPT": r"nampt|fk866|ap866|gmx1778|kpt-9274|ot-82",
    "STAT3": r"stat3|stattic|napabucasin|c188-9",
    "JAK1": r"jak|ruxolitinib|tofacitinib|baricitinib|momelotinib|fedratinib",
    "JAK2": r"jak|ruxolitinib|momelotinib|fedratinib",
    "IL6R": r"il-?6|tocilizumab|sarilumab",
    "OSMR": r"osm|jak",
    "DNMT1": r"azacitidine|decitabine|dnmt|guadecitabine",
    "DNMT3B": r"azacitidine|decitabine|dnmt|guadecitabine",
    "LYN": r"dasatinib|bosutinib|saracatinib|src|lyn",
    "SRC": r"dasatinib|bosutinib|saracatinib|src",
    "FYN": r"dasatinib|bosutinib|saracatinib|fyn|src",
    "ATR": r"atr|ceralasertib|berzosertib|ve-821|azd6738|m6620",
    "CHEK2": r"chk|chek|prexasertib|azd7762",
    "CDK7": r"cdk7|samuraciclib|thz1|sns-032",
    "CDK9": r"cdk9|flavopiridol|dinaciclib|atuveciclib|sns-032",
    "GLS": r"glutaminase|telaglenastat|cb-839",
    "MYC": r"myc|brd4|jq1|otx015|birabresib|proteasome|bortezomib",
    "HDAC1": r"hdac|vorinostat|panobinostat|romidepsin|belinostat",
    "HDAC2": r"hdac|vorinostat|panobinostat|romidepsin|belinostat",
    "EZH2": r"ezh2|tazemetostat|gsk126",
    "AURKA": r"aurora|alisertib|tozasertib",
    "BCL2L1": r"navitoclax|abt-263|bcl",
    "MCL1": r"mcl1|s63845|azd5991",
    "BCL2": r"venetoclax|abt-199|bcl2",
}


def ensure_dirs() -> None:
    for path in [RESULT_ROOT, REPORT_ROOT, P3, P9, RAW, SHARED]:
        path.mkdir(parents=True, exist_ok=True)


def read_table(path: Path, **kwargs) -> pd.DataFrame:
    sep = kwargs.pop("sep", "\t" if path.suffix in {".tsv", ".gz"} and ".csv" not in path.name else ",")
    return pd.read_csv(path, sep=sep, **kwargs)


def zscore(s: pd.Series) -> pd.Series:
    s = pd.to_numeric(s, errors="coerce")
    sd = s.std(skipna=True)
    if not np.isfinite(sd) or sd == 0:
        return s * np.nan
    return (s - s.mean(skipna=True)) / sd


def clean_gene_name(name: str) -> str:
    x = str(name).strip()
    x = re.sub(r"\s*\(\d+\)$", "", x)
    x = re.sub(r"\s*\[\S+\]$", "", x)
    return x.upper()


def numeric_corr(x: pd.Series, y: pd.Series, method: str = "spearman") -> tuple[float, float, int]:
    tmp = pd.concat([pd.to_numeric(x, errors="coerce"), pd.to_numeric(y, errors="coerce")], axis=1).dropna()
    if len(tmp) < 4 or tmp.iloc[:, 0].nunique() < 2 or tmp.iloc[:, 1].nunique() < 2:
        return np.nan, np.nan, len(tmp)
    stat = spearmanr(tmp.iloc[:, 0], tmp.iloc[:, 1]) if method == "spearman" else pearsonr(tmp.iloc[:, 0], tmp.iloc[:, 1])
    return float(stat.statistic), float(stat.pvalue), len(tmp)


def fdr_bh(pvals: Iterable[float]) -> list[float]:
    p = np.asarray([1.0 if pd.isna(v) else float(v) for v in pvals])
    n = len(p)
    order = np.argsort(p)
    ranked = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        idx = order[i]
        val = min(prev, p[idx] * n / (i + 1))
        ranked[idx] = val
        prev = val
    return ranked.tolist()


def score_signatures(expr_gene_by_sample: pd.DataFrame, signatures: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = expr_gene_by_sample.copy()
    expr.index = [clean_gene_name(g) for g in expr.index]
    expr = expr.groupby(expr.index).mean(numeric_only=True)
    gene_z = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1).replace(0, np.nan), axis=0)
    scores = {}
    coverage = []
    for name, genes in signatures.items():
        wanted = [clean_gene_name(g) for g in genes]
        present = [g for g in wanted if g in gene_z.index]
        scores[name] = gene_z.loc[present].mean(axis=0) if present else pd.Series(index=gene_z.columns, dtype=float)
        coverage.append(
            {
                "signature": name,
                "n_genes": len(wanted),
                "n_present": len(present),
                "genes_present": ";".join(present),
                "genes_missing": ";".join([g for g in wanted if g not in present]),
            }
        )
    return pd.DataFrame(scores).reset_index(names="sample_id"), pd.DataFrame(coverage)


def load_expression_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    gene_col = df.columns[0]
    df = df.rename(columns={gene_col: "gene_symbol"})
    df = df.set_index("gene_symbol")
    return df.apply(pd.to_numeric, errors="coerce")


def build_inventory_report() -> None:
    inventory_path = REPORT_ROOT / "local_inventory_files.txt"
    files = [line.strip() for line in inventory_path.read_text().splitlines() if line.strip()] if inventory_path.exists() else []
    buckets = {
        "A. Paper 1 lineage-state files": ["rai", "dm1", "tds", "lineage", "external", "gpl570", "gse76039", "landa"],
        "B. TCGA expression / clinical / mutation files": ["tcga", "thca", "mutation", "clinical", "pancan_geneexp", "phenotype"],
        "C. external expression validation files": ["external_expression", "gpl570", "gse", "expression_gene_log", "landa"],
        "D. immune / HLA / TLS / BCR files": ["immune", "hla", "tls", "bcr"],
        "E. scRNA files": ["scrna", "h5ad", "loom", "seurat", "single_cell"],
        "F. dependency / drug-response files": ["depmap", "ccle", "prism", "gdsc", "ctrp", "paper9", "dependency", "drug"],
    }
    lines = [
        "# Local Inventory Report",
        "",
        f"- Inventory file: `{inventory_path}`",
        f"- Matching files found: {len(files)}",
        "- Protected/raw WES/BAM/FASTQ are not executed from this inventory.",
        "",
    ]
    lower = [(f, f.lower()) for f in files]
    assigned = set()
    for title, keys in buckets.items():
        hits = [f for f, fl in lower if any(k in fl for k in keys)]
        assigned.update(hits)
        lines.append(f"## {title}")
        lines.append(f"- Count: {len(hits)}")
        for f in hits[:40]:
            lines.append(f"- `{f}`")
        if len(hits) > 40:
            lines.append(f"- ... {len(hits) - 40} more listed in `local_inventory_files.txt`")
        lines.append("")
    missing = [
        "Protected HLA LOH inputs: local BAM/WES plus matched normal metadata not identified for execution.",
        "Neoantigen execution inputs: local HLA genotypes plus VCF/peptide binding outputs not identified as a complete runnable set.",
        "Full public DepMap/PRISM/GDSC/CTRP files: checked in Paper 9 acquisition step and downloaded if absent.",
    ]
    lines.append("## G. missing files")
    for m in missing:
        lines.append(f"- {m}")
    (REPORT_ROOT / "01_local_inventory_report.md").write_text("\n".join(lines) + "\n")


def build_lineage_master() -> pd.DataFrame:
    rows = []
    wanted = [
        "dataset",
        "sample_id",
        "cohort",
        "cancer_type",
        "histology",
        "RAI_8_score",
        "DM1_like_score",
        "THYROID_NONOVERLAP_score",
        "TDS_like_score",
        "TF_collapse_score",
        "STAT3_AP1_DNMT_score",
        "source_file",
        "notes",
    ]

    def add_from(df: pd.DataFrame, source: Path, dataset_col: str | None = None) -> None:
        tmp = pd.DataFrame()
        tmp["dataset"] = df[dataset_col].astype(str) if dataset_col and dataset_col in df else source.parent.name
        tmp["sample_id"] = df["sample_id"].astype(str) if "sample_id" in df else df.get("sample", pd.Series(index=df.index, dtype=str)).astype(str)
        tmp["cohort"] = tmp["dataset"]
        tmp["cancer_type"] = "THCA/thyroid"
        tmp["histology"] = df.get("histology_clean", df.get("histology", df.get("disease_group", "")))
        for col in wanted[5:11]:
            tmp[col] = pd.to_numeric(df[col], errors="coerce") if col in df else np.nan
        tmp["source_file"] = str(source)
        tmp["notes"] = "already-local processed lineage scores; within-dataset scale retained"
        rows.append(tmp[wanted])

    source = ROOT / "project/results/p_external_expression_validation/external_sample_scores.tsv.gz"
    if source.exists():
        add_from(pd.read_csv(source, sep="\t"), source, "dataset")
    for source in [
        ROOT / "project/results/p_gpl570_validation/GSE29265_scores.tsv",
        ROOT / "project/results/p_gpl570_validation/GSE33630_scores.tsv",
        ROOT / "project/results/p_gpl570_validation/GSE65144_scores.tsv",
        ROOT / "project/results/p_landa_2016/scores.tsv",
    ]:
        if source.exists():
            df = pd.read_csv(source, sep="\t")
            add_from(df, source, "dataset" if "dataset" in df else None)

    source = ROOT / "project/results/v17p3/tables/A2_dm_score_full_cohort.tsv"
    if source.exists():
        df = pd.read_csv(source, sep="\t")
        tmp = pd.DataFrame()
        tmp["dataset"] = df.get("dataset", "TCGA-THCA").astype(str)
        tmp["sample_id"] = df["sample_id"].astype(str)
        tmp["cohort"] = "TCGA-THCA"
        tmp["cancer_type"] = "THCA"
        tmp["histology"] = df.get("histology_subtype", "")
        tmp["RAI_8_score"] = zscore(df.get("rai_score_recalc", pd.Series(np.nan, index=df.index)))
        tmp["DM1_like_score"] = zscore(df.get("dedifferentiation_proxy_score", pd.Series(np.nan, index=df.index)))
        tmp["THYROID_NONOVERLAP_score"] = np.nan
        tmp["TDS_like_score"] = zscore(df.get("tds_score", pd.Series(np.nan, index=df.index)))
        tmp["TF_collapse_score"] = np.nan
        tmp["STAT3_AP1_DNMT_score"] = np.nan
        tmp["source_file"] = str(source)
        tmp["notes"] = "TCGA local processed cohort; scores z-scored within TCGA for shared state table"
        rows.append(tmp[wanted])

    master = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=wanted)
    for score in wanted[5:11]:
        if score in master:
            master[score] = master.groupby("dataset", group_keys=False)[score].transform(lambda x: zscore(x) if x.notna().sum() > 2 else x)
    master = master.drop_duplicates(["dataset", "sample_id"], keep="first")
    out = SHARED / "lineage_state_master.tsv"
    master.to_csv(out, sep="\t", index=False)

    lines = [
        "# Lineage State Report",
        "",
        f"- Output: `{out}`",
        f"- Rows: {len(master)}",
        f"- Datasets: {master['dataset'].nunique() if len(master) else 0}",
        "- Platform pooling rule: no absolute cross-platform pooling; available score columns were retained or z-scored within dataset.",
        "- Claim boundary: these are research-use lineage-state variables, not clinical treatment variables.",
        "",
        "## Dataset Counts",
    ]
    for dataset, n in master["dataset"].value_counts().items():
        lines.append(f"- {dataset}: {n}")
    lines.extend(
        [
            "",
            "## Notes",
            "- External processed score files were prioritized where present.",
            "- TCGA lineage state was derived from already-local processed TCGA table columns, not from new protected data.",
            "- Missing score columns are left as unavailable rather than imputed.",
        ]
    )
    (REPORT_ROOT / "02_lineage_state_report.md").write_text("\n".join(lines) + "\n")
    return master


def write_signature_registry() -> pd.DataFrame:
    rows = []
    expected = {
        "DEDIFFERENTIATED_EPI": "higher_with_lineage_silencing",
        "STROMAL_EXCLUSION": "higher_with_lineage_silencing",
        "TREG": "context_dependent_suppressive",
        "M2_TAM": "context_dependent_suppressive",
        "MDSC_LIKE": "context_dependent_suppressive",
    }
    for name, genes in SIGNATURES.items():
        rows.append(
            {
                "signature": name,
                "signature_group": name,
                "genes": ";".join(genes),
                "n_genes": len(genes),
                "expected_direction": expected.get(name, "immune_readiness_axis"),
                "source": "canonical/provisional",
                "claim_guard": "signature score only; not an ICI response predictor",
            }
        )
    reg = pd.DataFrame(rows)
    reg.to_csv(P3 / "signature_registry.tsv", sep="\t", index=False)
    return reg


def run_paper3(master: pd.DataFrame) -> dict[str, object]:
    registry = write_signature_registry()
    expr_files = sorted((ROOT / "project/results/p_external_expression_validation").glob("*_expression_gene_log.tsv.gz"))
    landa = ROOT / "project/results/p_landa_2016/expression_matrix_log_gene.tsv.gz"
    if landa.exists():
        expr_files.append(landa)

    score_frames = []
    cov_frames = []
    for path in expr_files:
        dataset = path.name.replace("_expression_gene_log.tsv.gz", "").replace("expression_matrix_log_gene.tsv.gz", "GSE76039")
        expr = load_expression_matrix(path)
        scores, cov = score_signatures(expr, SIGNATURES)
        scores.insert(0, "dataset", dataset)
        scores["source_file"] = str(path)
        score_frames.append(scores)
        cov.insert(0, "dataset", dataset)
        cov["source_file"] = str(path)
        cov_frames.append(cov)

    immune = pd.concat(score_frames, ignore_index=True) if score_frames else pd.DataFrame()
    if not immune.empty:
        score_cols = [c for c in SIGNATURES if c in immune]
        immune[score_cols] = immune.groupby("dataset", group_keys=False)[score_cols].transform(zscore)

    tcga_immune_path = ROOT / "project/results/tables/immune_signatures_tcga.tsv"
    if tcga_immune_path.exists():
        tcga = pd.read_csv(tcga_immune_path, sep="\t")
        mapped = pd.DataFrame(
            {
                "dataset": "TCGA-THCA",
                "sample_id": tcga["sample_id"].astype(str),
                "IFNG_TIS": pd.to_numeric(tcga.get("IFN_gamma_response"), errors="coerce"),
                "CYTOLYTIC": pd.to_numeric(tcga.get("CD8_cytotoxic"), errors="coerce"),
                "TREG": pd.to_numeric(tcga.get("Treg"), errors="coerce"),
                "M2_TAM": pd.to_numeric(tcga.get("Macrophage_M2"), errors="coerce"),
                "T_CELL_EXHAUSTION": pd.to_numeric(tcga.get("Checkpoint_exhaustion"), errors="coerce"),
                "source_file": str(tcga_immune_path),
            }
        )
        for c in SIGNATURES:
            if c not in mapped:
                mapped[c] = np.nan
        mapped[list(SIGNATURES)] = mapped[list(SIGNATURES)].apply(zscore)
        immune = pd.concat([immune, mapped[["dataset", "sample_id", *SIGNATURES.keys(), "source_file"]]], ignore_index=True)

    immune.to_csv(P3 / "p3_immune_scores.tsv", sep="\t", index=False)
    if cov_frames:
        pd.concat(cov_frames, ignore_index=True).to_csv(P3 / "p3_signature_coverage.tsv", sep="\t", index=False)

    join = immune.merge(master, on=["dataset", "sample_id"], how="left", suffixes=("", "_lineage"))
    join.to_csv(P3 / "p3_lineage_immune_join.tsv", sep="\t", index=False)

    audit_rows = []
    for sig in SIGNATURES:
        r_dm, p_dm, n_dm = numeric_corr(join.get("DM1_like_score", pd.Series(dtype=float)), join[sig])
        r_mhc1, p_mhc1, n_mhc1 = numeric_corr(join[sig], join.get("MHC1_CORE", pd.Series(dtype=float)))
        r_mhc2, p_mhc2, n_mhc2 = numeric_corr(join[sig], join.get("MHC2_CORE", pd.Series(dtype=float)))
        r_tls, p_tls, n_tls = numeric_corr(join[sig], join.get("TLS_CABRITA", pd.Series(dtype=float)))
        dyn = float(pd.to_numeric(join[sig], errors="coerce").quantile(0.9) - pd.to_numeric(join[sig], errors="coerce").quantile(0.1))
        present_n = int(pd.to_numeric(join[sig], errors="coerce").notna().sum())
        exp = registry.loc[registry.signature == sig, "expected_direction"].iloc[0]
        verdict = "PASS"
        if present_n < 10 or not np.isfinite(dyn) or dyn < 0.5:
            verdict = "COLLAPSE"
        elif exp == "immune_readiness_axis" and np.isfinite(r_dm) and r_dm < -0.25 and p_dm < 0.05:
            verdict = "FLIP"
        elif exp.startswith("higher") and np.isfinite(r_dm) and r_dm < -0.25 and p_dm < 0.05:
            verdict = "FLIP"
        elif pd.isna(r_dm):
            verdict = "AMBIGUOUS"
        audit_rows.append(
            {
                "signature": sig,
                "expected_direction": exp,
                "n_scored_samples": present_n,
                "dynamic_range_p90_p10": dyn,
                "spearman_DM1_like": r_dm,
                "p_DM1_like": p_dm,
                "n_DM1_like": n_dm,
                "spearman_MHC1": r_mhc1,
                "spearman_MHC2": r_mhc2,
                "spearman_TLS_CABRITA": r_tls,
                "verdict": verdict,
            }
        )
    audit = pd.DataFrame(audit_rows)
    audit["q_DM1_like"] = fdr_bh(audit["p_DM1_like"])
    audit.to_csv(P3 / "p3_dial_audit.tsv", sep="\t", index=False)
    pass_sigs = set(audit.loc[audit.verdict == "PASS", "signature"])

    readiness_terms = {
        "IFNG_TIS": 1,
        "CYTOLYTIC": 1,
        "MHC1_CORE": 1,
        "MHC2_CORE": 1,
        "TLS_CABRITA": 1,
        "CXCL13_AXIS": 1,
        "TREG": -1,
        "M2_TAM": -1,
        "MDSC_LIKE": -1,
        "STROMAL_EXCLUSION": -1,
    }
    ready = join[["dataset", "sample_id", "DM1_like_score", "RAI_8_score", "histology"]].copy()
    used_terms = [t for t in readiness_terms if t in pass_sigs and t in join]
    comp = []
    for t in used_terms:
        comp.append(readiness_terms[t] * pd.to_numeric(join[t], errors="coerce"))
        ready[t] = join[t]
    ready["ICI_READINESS"] = pd.concat(comp, axis=1).mean(axis=1) if comp else np.nan
    ready["HLA_LOH_status"] = "unavailable"
    ready["neoantigen_presentability"] = "unavailable"
    ready["score_terms_used"] = ";".join(used_terms)
    ready["claim_guard"] = "ICI readiness score only; no response prediction or treatment recommendation"
    ready.to_csv(P3 / "p3_ici_readiness_scores.tsv", sep="\t", index=False)

    cluster_cols = sorted(set(used_terms + ["DM1_like_score", "RAI_8_score"]))
    clust = join[["dataset", "sample_id"] + cluster_cols].copy()
    X = clust[cluster_cols].apply(pd.to_numeric, errors="coerce")
    valid = X.dropna(thresh=max(3, int(len(cluster_cols) * 0.5))).copy()
    labels = pd.Series("unassigned", index=clust.index)
    if len(valid) >= 6:
        Ximp = valid.fillna(valid.mean())
        Xs = StandardScaler().fit_transform(Ximp)
        k = 4 if len(valid) >= 20 else 3
        labels.loc[valid.index] = [f"ecotype_{i+1}" for i in KMeans(n_clusters=k, random_state=7, n_init=20).fit_predict(Xs)]
    ecotypes = clust[["dataset", "sample_id"]].copy()
    ecotypes["immune_ecotype"] = labels.values
    ecotypes = ecotypes.merge(ready[["dataset", "sample_id", "ICI_READINESS"]], on=["dataset", "sample_id"], how="left")
    ecotypes.to_csv(P3 / "p3_immune_ecotypes.tsv", sep="\t", index=False)

    sns.set_theme(style="whitegrid")
    heat_cols = [c for c in used_terms if c in join]
    if heat_cols and len(join) > 2:
        heat = join[heat_cols].apply(pd.to_numeric, errors="coerce").dropna(how="all")
        heat = heat.fillna(heat.mean()).clip(-3, 3)
        if len(heat) > 80:
            heat = heat.sample(80, random_state=7)
        plt.figure(figsize=(max(7, len(heat_cols) * 0.55), 8))
        sns.heatmap(heat, cmap="vlag", center=0, yticklabels=False)
        plt.tight_layout()
        plt.savefig(P3 / "fig_p3_immune_heatmap.png", dpi=200)
        plt.close()
    else:
        placeholder_plot(P3 / "fig_p3_immune_heatmap.png", "Insufficient PASS immune signatures")

    plt.figure(figsize=(6, 4.5))
    tmp = ready.dropna(subset=["DM1_like_score", "ICI_READINESS"])
    if len(tmp) > 2:
        sns.scatterplot(data=tmp, x="DM1_like_score", y="ICI_READINESS", hue="dataset", s=28, linewidth=0)
    else:
        plt.text(0.5, 0.5, "Insufficient matched lineage/readiness data", ha="center", va="center")
    plt.tight_layout()
    plt.savefig(P3 / "fig_p3_lineage_vs_ici_readiness.png", dpi=200)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    order = audit.sort_values("spearman_DM1_like")
    sns.barplot(data=order, x="spearman_DM1_like", y="signature", hue="verdict", dodge=False)
    plt.axvline(0, color="black", lw=0.8)
    plt.tight_layout()
    plt.savefig(P3 / "fig_p3_dial_audit.png", dpi=200)
    plt.close()

    write_paper3_blockers()
    write_paper3_report(registry, immune, join, audit, ready, ecotypes, used_terms)
    return {
        "immune_rows": len(immune),
        "joined_rows": len(join),
        "pass_signatures": sorted(pass_sigs),
        "used_readiness_terms": used_terms,
        "top_readiness": ready.sort_values("ICI_READINESS", ascending=False).head(5)[["dataset", "sample_id", "ICI_READINESS"]].to_dict("records")
        if "ICI_READINESS" in ready
        else [],
    }


def placeholder_plot(path: Path, text: str) -> None:
    plt.figure(figsize=(6, 4))
    plt.text(0.5, 0.5, text, ha="center", va="center")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


def write_paper3_blockers() -> None:
    patterns = ["*.vcf", "*.vcf.gz", "*hla*", "*HLA*", "*lohhla*", "*LOHHLA*", "*neoantigen*", "*netmhc*", "*wes*", "*WES*"]
    hits = []
    for pat in patterns:
        for p in (ROOT / "project/results").rglob(pat):
            if p.is_file() and "p3_p9_full_execution" not in str(p):
                hits.append(p)
    complete = any("lohhla" in p.name.lower() for p in hits) and any("neoantigen" in p.name.lower() for p in hits)
    lines = [
        "# Paper 3 Protected Access Blockers",
        "",
        "- Protected data download: not attempted.",
        "- LOHHLA execution: not run.",
        "- NetMHCpan/neoantigen calling: not run.",
        f"- Complete local HLA LOH/neoantigen runnable set detected: {'yes' if complete else 'no'}",
        "",
        "## Local Candidate Files Observed",
    ]
    if hits:
        for p in sorted(set(hits))[:80]:
            lines.append(f"- `{p.relative_to(ROOT)}`")
        if len(hits) > 80:
            lines.append(f"- ... {len(hits) - 80} additional candidate paths omitted from this blocker report.")
    else:
        lines.append("- None found.")
    lines.extend(
        [
            "",
            "## Blocker Decision",
            "- HLA_LOH_status remains `unavailable`.",
            "- neoantigen_presentability remains `unavailable`.",
            "- No imputation was performed.",
        ]
    )
    (REPORT_ROOT / "03_paper3_protected_access_blockers.md").write_text("\n".join(lines) + "\n")


def write_paper3_report(registry: pd.DataFrame, immune: pd.DataFrame, join: pd.DataFrame, audit: pd.DataFrame, ready: pd.DataFrame, ecotypes: pd.DataFrame, used_terms: list[str]) -> None:
    lines = [
        "# Paper 3 ICI Full-Public Report",
        "",
        "## 1. Data used",
        f"- Local/public expression signature score rows: {len(immune)}",
        f"- Lineage-immune joined rows: {len(join)}",
        "- TCGA immune scores used only from already-local processed immune signature table when available.",
        "- HLA LOH and neoantigen data were not imputed.",
        "",
        "## 2. Signature registry",
        f"- Registry output: `{P3 / 'signature_registry.tsv'}`",
        f"- Signatures registered: {len(registry)}",
        "",
        "## 3. DIAL audit result",
        f"- PASS: {(audit.verdict == 'PASS').sum()}",
        f"- FLIP: {(audit.verdict == 'FLIP').sum()}",
        f"- COLLAPSE: {(audit.verdict == 'COLLAPSE').sum()}",
        f"- AMBIGUOUS: {(audit.verdict == 'AMBIGUOUS').sum()}",
        "",
        "## 4. ICI readiness score",
        f"- Terms used after PASS filter: {', '.join(used_terms) if used_terms else 'none'}",
        f"- Scored samples: {ready['ICI_READINESS'].notna().sum() if 'ICI_READINESS' in ready else 0}",
        "",
        "## 5. Immune ecotypes",
        f"- Ecotype rows: {len(ecotypes)}",
    ]
    for label, n in ecotypes["immune_ecotype"].value_counts().items():
        lines.append(f"- {label}: {n}")
    lines.extend(
        [
            "",
            "## 6. HLA/neoantigen blockers",
            "- See `03_paper3_protected_access_blockers.md`.",
            "",
            "## 7. Claim boundary",
            "- Allowed: ICI readiness, vulnerability, and research prioritization language.",
            "- Forbidden: ICI response predictor, treatment recommendation, patient treatment selection.",
            "",
            "## 8. Figure list",
            f"- `{P3 / 'fig_p3_immune_heatmap.png'}`",
            f"- `{P3 / 'fig_p3_lineage_vs_ici_readiness.png'}`",
            f"- `{P3 / 'fig_p3_dial_audit.png'}`",
            "",
            "## 9. Next gates",
            "- Add protected HLA LOH or neoantigen analyses only if local/credentialed inputs are explicitly provided.",
            "- Treat weak or negative immune associations as honest negative results.",
        ]
    )
    (REPORT_ROOT / "04_paper3_ici_full_public_report.md").write_text("\n".join(lines) + "\n")


def parse_size(size: str | None) -> int | None:
    if not size:
        return None
    m = re.search(r"([\d.]+)\s*([KMGTP]?B)", size.upper())
    if not m:
        return None
    val = float(m.group(1))
    unit = m.group(2)
    mult = {"KB": 1024, "MB": 1024**2, "GB": 1024**3, "TB": 1024**4, "B": 1}[unit]
    return int(val * mult)


@dataclass
class DownloadSpec:
    release: str
    file_name: str
    role: str
    required: bool = True


def depmap_rows() -> list[dict]:
    data = requests.get(DEPMAP_API, timeout=120, headers={"User-Agent": "Mozilla/5.0"}).json()
    rows = []
    for key in ["currentRelease", "allDownloads", "table"]:
        if isinstance(data.get(key), list):
            rows.extend(data[key])
    return rows


def select_downloads(rows: list[dict]) -> list[dict]:
    specs = [
        DownloadSpec("DepMap Public 26Q1", "Model.csv", "model_metadata"),
        DownloadSpec("DepMap Public 26Q1", "CRISPRGeneEffect.csv", "crispr_gene_effect"),
        DownloadSpec("DepMap Public 26Q1", "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", "ccle_expression"),
        DownloadSpec("DepMap Public 26Q1", "AchillesCommonEssentialControls.csv", "common_essential"),
        DownloadSpec("DepMap Public 26Q1", "AchillesNonessentialControls.csv", "nonessential_controls"),
        DownloadSpec("PRISM Primary Repurposing DepMap Public 24Q2", "Repurposing_Public_24Q2_LFC_COLLAPSED.csv", "prism_lfc_collapsed", False),
        DownloadSpec("PRISM Primary Repurposing DepMap Public 24Q2", "Repurposing_Public_24Q2_Treatment_Meta_Data.csv", "prism_treatment_metadata", False),
        DownloadSpec("PRISM Primary Repurposing DepMap Public 24Q2", "Repurposing_Public_24Q2_Cell_Line_Meta_Data.csv", "prism_cell_metadata", False),
        DownloadSpec("Harmonized GDSC 25Q2", "GDSC1AUCMatrix.csv", "gdsc1_auc", False),
        DownloadSpec("Harmonized GDSC 25Q2", "GDSC1Log2IC50Matrix.csv", "gdsc1_ic50", False),
        DownloadSpec("Harmonized GDSC 25Q2", "GDSC1Log2ViabilityCollapsedConditions.csv", "gdsc1_conditions", False),
        DownloadSpec("Harmonized GDSC 25Q2", "GDSC2AUCMatrix.csv", "gdsc2_auc", False),
        DownloadSpec("Harmonized GDSC 25Q2", "GDSC2Log2IC50Matrix.csv", "gdsc2_ic50", False),
        DownloadSpec("Harmonized GDSC 25Q2", "GDSC2Log2ViabilityCollapsedConditions.csv", "gdsc2_conditions", False),
        DownloadSpec("Harmonized CTD^2 25Q2", "CTRPAUCMatrix.csv", "ctrp_auc", False),
        DownloadSpec("Harmonized CTD^2 25Q2", "CTRPLog2IC50Matrix.csv", "ctrp_ic50", False),
        DownloadSpec("Harmonized CTD^2 25Q2", "CTRPLog2ViabilityCollapsedConditions.csv", "ctrp_conditions", False),
    ]
    selected = []
    for spec in specs:
        match = next((r for r in rows if r.get("releaseName") == spec.release and r.get("fileName") == spec.file_name), None)
        if match is None:
            if spec.required:
                raise RuntimeError(f"Required DepMap file not found in API: {spec.release} / {spec.file_name}")
            continue
        item = dict(match)
        item["role"] = spec.role
        item["expected_bytes"] = parse_size(item.get("size"))
        selected.append(item)
    return selected


def download_public_files() -> pd.DataFrame:
    rows = depmap_rows()
    selected = select_downloads(rows)
    total = 0
    manifest = []
    for item in selected:
        expected = item.get("expected_bytes")
        if expected and expected > MAX_SINGLE_DOWNLOAD:
            raise RuntimeError(f"Single download exceeds 10 GB guard: {item['fileName']} {item.get('size')}")
        total += expected or 0
    if total > MAX_TOTAL_DOWNLOAD:
        raise RuntimeError(f"Selected downloads exceed 50 GB guard: {total / 1024**3:.1f} GB")

    for item in selected:
        out = RAW / item["fileName"]
        url = item["downloadUrl"]
        full_url = url if url.startswith("http") else DEPMAP_BASE + url
        status = "existing"
        if not out.exists() or out.stat().st_size < 100:
            status = "downloaded"
            with requests.get(full_url, stream=True, timeout=300, headers={"User-Agent": "Mozilla/5.0"}) as r:
                r.raise_for_status()
                content_length = int(r.headers.get("content-length") or 0)
                if content_length > MAX_SINGLE_DOWNLOAD:
                    raise RuntimeError(f"Single download content-length exceeds 10 GB guard: {item['fileName']}")
                tmp = out.with_suffix(out.suffix + ".part")
                with tmp.open("wb") as handle:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            handle.write(chunk)
                tmp.replace(out)
        manifest.append(
            {
                "role": item["role"],
                "release": item.get("releaseName"),
                "file_name": item["fileName"],
                "size_label": item.get("size"),
                "expected_bytes": item.get("expected_bytes"),
                "actual_bytes": out.stat().st_size if out.exists() else np.nan,
                "status": status,
                "source_url": full_url,
                "local_path": str(out),
            }
        )
    mdf = pd.DataFrame(manifest)
    mdf.to_csv(P9 / "p9_public_download_manifest.tsv", sep="\t", index=False)
    return mdf


def load_depmap_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "ModelID" in df.columns:
        model_ids = df["ModelID"].astype(str)
        metadata_cols = {
            "Unnamed: 0",
            "SequencingID",
            "ModelConditionID",
            "ModelID",
            "IsDefaultEntryForMC",
            "IsDefaultEntryForModel",
        }
        df = df.drop(columns=[c for c in metadata_cols if c in df.columns])
        df.index = model_ids
    else:
        first = df.columns[0]
        df = df.rename(columns={first: "ModelID"}).set_index("ModelID")
    df.columns = [clean_gene_name(c) for c in df.columns]
    df = df.groupby(level=0).first()
    df = df.T.groupby(level=0).mean(numeric_only=True).T
    return df.apply(pd.to_numeric, errors="coerce")


def find_thyroid_models(model: pd.DataFrame) -> pd.Series:
    text_cols = [c for c in model.columns if model[c].dtype == object]
    text = model[text_cols].fillna("").astype(str).agg(" ".join, axis=1).str.lower()
    return text.str.contains("thyroid|thca|papillary thyroid|anaplastic thyroid|follicular thyroid", regex=True)


def run_paper9() -> dict[str, object]:
    manifest = download_public_files()
    model = pd.read_csv(RAW / "Model.csv")
    model_id_col = "ModelID" if "ModelID" in model.columns else model.columns[0]
    model = model.rename(columns={model_id_col: "ModelID"}).set_index("ModelID", drop=False)
    model["is_thyroid"] = find_thyroid_models(model)
    expr = load_depmap_matrix(RAW / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
    gene_effect = load_depmap_matrix(RAW / "CRISPRGeneEffect.csv")

    common = set()
    common_path = RAW / "AchillesCommonEssentialControls.csv"
    if common_path.exists():
        cdf = pd.read_csv(common_path)
        vals = cdf.iloc[:, 0].dropna().astype(str).map(clean_gene_name)
        common = set(vals)

    expr_samples_by_gene = expr.copy()
    expr_gene_by_sample = expr_samples_by_gene.T
    lineage, lin_cov = score_signatures(expr_gene_by_sample, LINEAGE_SIGNATURES)
    lineage = lineage.rename(columns={"sample_id": "ModelID"})
    lineage["DM1_like_score"] = -pd.to_numeric(lineage["RAI_8_score"], errors="coerce")
    lineage = lineage.merge(model.reset_index(drop=True), on="ModelID", how="left")
    lineage["is_thyroid_model"] = lineage["is_thyroid"].fillna(False)
    lineage["cancer_type"] = lineage.get("OncotreeLineage", lineage.get("OncotreePrimaryDisease", "unknown"))
    lineage["lineage_silenced_score"] = zscore(lineage["DM1_like_score"])
    lineage.to_csv(P9 / "p9_cellline_lineage_scores.tsv", sep="\t", index=False)
    lin_cov.to_csv(P9 / "p9_cellline_lineage_signature_coverage.tsv", sep="\t", index=False)

    dep = gene_effect.merge(lineage[["ModelID", "lineage_silenced_score", "is_thyroid_model", "cancer_type"]], left_index=True, right_on="ModelID", how="inner")
    rows = []
    for gene, klass in CANDIDATES.items():
        col = clean_gene_name(gene)
        if col not in dep.columns:
            rows.append({"target": gene, "class": klass, "in_depmap": False})
            continue
        y = pd.to_numeric(dep[col], errors="coerce")
        x = dep["lineage_silenced_score"]
        r, p, n = numeric_corr(x, y)
        td = dep[dep["is_thyroid_model"] == True]
        rt, pt, nt = numeric_corr(td["lineage_silenced_score"], pd.to_numeric(td[col], errors="coerce"))
        centered_y = y - dep.groupby("cancer_type")[col].transform("mean")
        centered_x = x - dep.groupby("cancer_type")["lineage_silenced_score"].transform("mean")
        rr, pr, nr = numeric_corr(centered_x, centered_y)
        rows.append(
            {
                "target": gene,
                "class": klass,
                "in_depmap": True,
                "spearman_lineage_dependency": r,
                "p_lineage_dependency": p,
                "n_models": n,
                "spearman_lineage_dependency_thyroid": rt,
                "p_lineage_dependency_thyroid": pt,
                "n_thyroid_models": nt,
                "pan_cancer_lineage_centered_spearman": rr,
                "pan_cancer_lineage_centered_p": pr,
                "common_essential_control": clean_gene_name(gene) in common,
                "median_gene_effect": float(y.median(skipna=True)),
                "artifact_risk": "HIGH_common_essential_or_panlethal" if clean_gene_name(gene) in common or y.median(skipna=True) < -0.7 else "moderate_or_low",
                "dependency_interpretation": "negative correlation means higher lineage silencing associates with stronger dependency",
            }
        )
    dep_assoc = pd.DataFrame(rows)
    dep_assoc["q_lineage_dependency"] = fdr_bh(dep_assoc.get("p_lineage_dependency", pd.Series(dtype=float)))
    dep_assoc.to_csv(P9 / "p9_candidate_dependency_associations.tsv", sep="\t", index=False)

    drug = run_drug_associations(lineage)
    drug.to_csv(P9 / "p9_drug_response_associations.tsv", sep="\t", index=False)

    expr_support_rows = []
    for gene, klass in CANDIDATES.items():
        col = clean_gene_name(gene)
        if col in expr.columns:
            e = expr[[col]].rename(columns={col: "target_expression"}).merge(lineage[["ModelID", "lineage_silenced_score", "is_thyroid_model"]], left_index=True, right_on="ModelID", how="inner")
            r, p, n = numeric_corr(e["lineage_silenced_score"], e["target_expression"])
            rt, pt, nt = numeric_corr(e.loc[e.is_thyroid_model == True, "lineage_silenced_score"], e.loc[e.is_thyroid_model == True, "target_expression"])
        else:
            r = p = rt = pt = np.nan
            n = nt = 0
        expr_support_rows.append({"target": gene, "target_expression_spearman": r, "target_expression_p": p, "target_expression_n": n, "thyroid_target_expression_spearman": rt, "thyroid_target_expression_p": pt, "thyroid_expression_n": nt})
    expr_support = pd.DataFrame(expr_support_rows)

    ranking = dep_assoc.merge(expr_support, on="target", how="left")
    drug_summary = (
        drug.groupby("target")
        .agg(
            drug_response_support=("p_lineage_sensitivity", lambda x: int((pd.to_numeric(x, errors="coerce") < 0.05).sum())),
            best_drug_p=("p_lineage_sensitivity", "min"),
            best_drug_abs_r=("spearman_lineage_sensitivity_metric", lambda x: np.nanmax(np.abs(pd.to_numeric(x, errors="coerce"))) if len(x) else np.nan),
        )
        .reset_index()
        if not drug.empty and "target" in drug
        else pd.DataFrame({"target": list(CANDIDATES), "drug_response_support": 0, "best_drug_p": np.nan, "best_drug_abs_r": np.nan})
    )
    ranking = ranking.merge(drug_summary, on="target", how="left")
    ranking["drug_response_support"] = ranking["drug_response_support"].fillna(0)
    ranking["dependency_association"] = ranking["spearman_lineage_dependency"]
    ranking["target_expression_support"] = ranking["target_expression_spearman"]
    ranking["thyroid_model_support"] = ranking.apply(
        lambda r: "nominal" if pd.notna(r.get("p_lineage_dependency_thyroid")) and r.get("p_lineage_dependency_thyroid") < 0.1 else ("underpowered" if r.get("n_thyroid_models", 0) < 8 else "none"),
        axis=1,
    )
    ranking["pan_cancer_specificity"] = ranking["pan_cancer_lineage_centered_spearman"].apply(lambda x: "lineage_centered_signal" if pd.notna(x) and abs(x) >= 0.1 else "weak_or_none")
    ranking["druggability"] = ranking["target"].map(lambda g: "drugged_or_tool_compound" if g in DRUG_REGEX else "target_or_pathway_only")

    def tier(row: pd.Series) -> tuple[str, str]:
        dep_ok = pd.notna(row.get("p_lineage_dependency")) and row["p_lineage_dependency"] < 0.05 and row.get("spearman_lineage_dependency", 0) < 0
        drug_ok = row.get("drug_response_support", 0) > 0
        expr_ok = pd.notna(row.get("target_expression_p")) and row["target_expression_p"] < 0.05
        artifact_high = str(row.get("artifact_risk", "")).startswith("HIGH")
        if artifact_high:
            return "artifact_flagged", "DROP" if not drug_ok else "LOW"
        if dep_ok and drug_ok:
            return "dependency_plus_drug_nominal", "HIGH_PRIORITY"
        if dep_ok and (expr_ok or row.get("thyroid_model_support") == "nominal"):
            return "dependency_plus_expression_or_thyroid", "MEDIUM"
        if dep_ok:
            return "dependency_only_nominal", "NEEDS_WETLAB"
        if drug_ok:
            return "drug_only_nominal", "LOW"
        return "weak_or_negative", "LOW"

    tv = ranking.apply(tier, axis=1, result_type="expand")
    ranking["evidence_tier"] = tv[0]
    ranking["verdict"] = tv[1]
    ranking[
        [
            "target",
            "class",
            "dependency_association",
            "drug_response_support",
            "target_expression_support",
            "thyroid_model_support",
            "pan_cancer_specificity",
            "druggability",
            "artifact_risk",
            "evidence_tier",
            "verdict",
            "p_lineage_dependency",
            "q_lineage_dependency",
            "n_models",
            "n_thyroid_models",
        ]
    ].to_csv(P9 / "p9_vulnerability_ranked_targets.tsv", sep="\t", index=False)

    make_paper9_figures(dep_assoc, ranking, drug)
    write_paper9_report(manifest, lineage, dep_assoc, drug, ranking)
    top_ranked = (
        ranking.assign(_priority=ranking["verdict"].map(VERDICT_PRIORITY).fillna(9))
        .sort_values(["_priority", "q_lineage_dependency", "p_lineage_dependency"], na_position="last")
        .head(8)
    )
    return {
        "downloaded_files": manifest.to_dict("records"),
        "thyroid_models": int(lineage["is_thyroid_model"].sum()),
        "dep_rows": len(dep_assoc),
        "drug_rows": len(drug),
        "top_targets": top_ranked[["target", "verdict", "evidence_tier", "dependency_association"]].to_dict("records"),
    }


def read_drug_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    first = df.columns[0]
    return df.rename(columns={first: "ModelID"})


def metadata_label_map(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    df = pd.read_csv(path)
    cols = df.columns.tolist()
    id_cols = [c for c in cols if re.search(r"broad|compound|treatment|drug|condition|id", c, re.I)]
    label_cols = [c for c in cols if re.search(r"name|target|moa|mechanism|perturbagen|compound", c, re.I)]
    mapping = {}
    for _, r in df.iterrows():
        labels = " ".join(str(r[c]) for c in label_cols if c in df and pd.notna(r[c]))
        for c in id_cols:
            if c in df and pd.notna(r[c]):
                mapping[str(r[c])] = labels
    return mapping


def run_drug_associations(lineage: pd.DataFrame) -> pd.DataFrame:
    rows = []
    matrices = [
        ("PRISM24Q2_LFC_COLLAPSED", RAW / "Repurposing_Public_24Q2_LFC_COLLAPSED.csv", metadata_label_map(RAW / "Repurposing_Public_24Q2_Treatment_Meta_Data.csv"), "lower_more_sensitive"),
        ("GDSC1_AUC", RAW / "GDSC1AUCMatrix.csv", metadata_label_map(RAW / "GDSC1Log2ViabilityCollapsedConditions.csv"), "lower_more_sensitive"),
        ("GDSC2_AUC", RAW / "GDSC2AUCMatrix.csv", metadata_label_map(RAW / "GDSC2Log2ViabilityCollapsedConditions.csv"), "lower_more_sensitive"),
        ("CTRP_AUC", RAW / "CTRPAUCMatrix.csv", metadata_label_map(RAW / "CTRPLog2ViabilityCollapsedConditions.csv"), "lower_more_sensitive"),
    ]
    lineage_small = lineage[["ModelID", "lineage_silenced_score", "is_thyroid_model"]].dropna(subset=["lineage_silenced_score"])
    for source, path, label_map, direction in matrices:
        if not path.exists():
            continue
        try:
            mat = read_drug_matrix(path)
        except Exception as exc:
            rows.append({"source": source, "target": "LOAD_FAILED", "drug_id": path.name, "notes": str(exc)})
            continue
        joined = mat.merge(lineage_small, on="ModelID", how="inner")
        if joined.empty:
            continue
        numeric_cols = [c for c in mat.columns if c != "ModelID" and pd.api.types.is_numeric_dtype(joined[c])]
        for target, pattern in DRUG_REGEX.items():
            rx = re.compile(pattern, re.I)
            selected = []
            for c in numeric_cols:
                label = f"{c} {label_map.get(str(c), '')}"
                if rx.search(label):
                    selected.append(c)
            for c in selected[:20]:
                metric = pd.to_numeric(joined[c], errors="coerce")
                sensitivity_metric = -metric if direction == "lower_more_sensitive" else metric
                r, p, n = numeric_corr(joined["lineage_silenced_score"], sensitivity_metric)
                jt = joined[joined["is_thyroid_model"] == True]
                rt, pt, nt = numeric_corr(jt["lineage_silenced_score"], -pd.to_numeric(jt[c], errors="coerce"))
                rows.append(
                    {
                        "source": source,
                        "target": target,
                        "drug_id": c,
                        "drug_label": label_map.get(str(c), ""),
                        "spearman_lineage_sensitivity_metric": r,
                        "p_lineage_sensitivity": p,
                        "n_models": n,
                        "spearman_lineage_sensitivity_thyroid": rt,
                        "p_lineage_sensitivity_thyroid": pt,
                        "n_thyroid_models": nt,
                        "metric_direction": direction,
                        "notes": "positive correlation means higher lineage silencing associates with greater sensitivity after direction normalization",
                    }
                )
    return pd.DataFrame(rows)


def make_paper9_figures(dep_assoc: pd.DataFrame, ranking: pd.DataFrame, drug: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid")
    vals = dep_assoc.set_index("target")[["spearman_lineage_dependency", "spearman_lineage_dependency_thyroid", "pan_cancer_lineage_centered_spearman"]]
    plt.figure(figsize=(6.5, max(5, len(vals) * 0.22)))
    sns.heatmap(vals.apply(pd.to_numeric, errors="coerce"), cmap="vlag", center=0, annot=False)
    plt.tight_layout()
    plt.savefig(P9 / "fig_p9_lineage_dependency_heatmap.png", dpi=200)
    plt.close()

    plot_rank = ranking.copy()
    plot_rank["score"] = -pd.to_numeric(plot_rank["dependency_association"], errors="coerce")
    plot_rank = plot_rank.sort_values("score", ascending=False).head(20)
    plt.figure(figsize=(7, 5))
    sns.barplot(data=plot_rank, x="score", y="target", hue="verdict", dodge=False)
    plt.xlabel("Lineage-silencing dependency association (-Spearman)")
    plt.tight_layout()
    plt.savefig(P9 / "fig_p9_candidate_target_rank.png", dpi=200)
    plt.close()

    if not drug.empty and "target" in drug:
        dsum = drug.pivot_table(index="target", columns="source", values="spearman_lineage_sensitivity_metric", aggfunc="max")
        plt.figure(figsize=(7, max(4, len(dsum) * 0.25)))
        sns.heatmap(dsum, cmap="vlag", center=0)
        plt.tight_layout()
        plt.savefig(P9 / "fig_p9_drug_response_map.png", dpi=200)
        plt.close()
    else:
        placeholder_plot(P9 / "fig_p9_drug_response_map.png", "No matched drug response associations")

    tier = pd.crosstab(ranking["class"], ranking["verdict"])
    plt.figure(figsize=(7, max(4, len(tier) * 0.35)))
    sns.heatmap(tier, cmap="YlGnBu", annot=True, fmt="d")
    plt.tight_layout()
    plt.savefig(P9 / "fig_p9_evidence_tier_matrix.png", dpi=200)
    plt.close()


def write_paper9_report(manifest: pd.DataFrame, lineage: pd.DataFrame, dep_assoc: pd.DataFrame, drug: pd.DataFrame, ranking: pd.DataFrame) -> None:
    top = (
        ranking.assign(_priority=ranking["verdict"].map(VERDICT_PRIORITY).fillna(9))
        .sort_values(["_priority", "q_lineage_dependency", "p_lineage_dependency"], na_position="last")
        .head(10)
    )
    lines = [
        "# Paper 9 Synthetic Lethality Full-Public Report",
        "",
        "## 1. Data downloaded/used",
    ]
    for _, r in manifest.iterrows():
        lines.append(f"- {r['role']}: `{r['file_name']}` ({r['release']}, {r['size_label']}, {r['status']})")
    lines.extend(
        [
            "",
            "## 2. Model/cell-line coverage",
            f"- DepMap expression-scored models: {len(lineage)}",
            f"- Thyroid-flagged models: {int(lineage['is_thyroid_model'].sum())}",
            "",
            "## 3. Lineage-state scoring result",
            "- RAI/thyroid differentiation and TF-collapse scores were computed from public processed CCLE/DepMap expression.",
            "- DM1-like/lineage-silenced score is the inverse of RAI_8 within the DepMap expression context.",
            "",
            "## 4. Dependency association result",
            f"- Candidate targets tested: {len(dep_assoc)}",
            f"- Nominal pan-cancer dependency associations with p < 0.05: {(pd.to_numeric(dep_assoc['p_lineage_dependency'], errors='coerce') < 0.05).sum()}",
            "- Negative dependency correlations indicate stronger dependency as lineage silencing increases.",
            "",
            "## 5. Drug-response result",
            f"- Drug association rows: {len(drug)}",
            "- Drug metrics are direction-normalized where lower viability/AUC/LFC indicates greater sensitivity.",
            "",
            "## 6. Top candidate targets",
        ]
    )
    for _, r in top.iterrows():
        lines.append(f"- {r['target']}: {r['verdict']} ({r['evidence_tier']}), dependency association={r['dependency_association']:.3g}")
    lines.extend(
        [
            "",
            "## 7. What failed",
            "- Thyroid-only dependency associations are underpowered if few thyroid models overlap the CRISPR matrix.",
            "- Drug metadata matching is keyword-based and therefore candidate-generating, not confirmatory.",
            "- Protected patient-level data and raw sequencing were not used.",
            "",
            "## 8. Claim boundary",
            "- Allowed: candidate dependency nomination and research-use vulnerability prioritization.",
            "- Forbidden: clinical treatment selection, proven drug response, or patient treatment recommendation.",
            "",
            "## 9. Wet-lab validation priority",
            "- Prioritize targets with negative lineage-dependency association, non-high artifact risk, expression support, and drug/dependency concordance.",
            "- Treat common-essential or proliferation-linked candidates as artifact-risk until validated with orthogonal assays.",
        ]
    )
    (REPORT_ROOT / "05_paper9_synthetic_lethality_full_public_report.md").write_text("\n".join(lines) + "\n")


def write_boundary_summary(p3_summary: dict[str, object], p9_summary: dict[str, object]) -> None:
    p3_terms = set(p3_summary.get("used_readiness_terms", []))
    p9_targets = set(CANDIDATES)
    overlap = sorted(p3_terms.intersection(p9_targets))
    lines = [
        "# P3/P9 Execution Boundary Summary",
        "",
        "## Paper 3 result summary",
        f"- Immune score rows: {p3_summary.get('immune_rows')}",
        f"- Joined lineage/immune rows: {p3_summary.get('joined_rows')}",
        f"- PASS signatures: {', '.join(p3_summary.get('pass_signatures', []))}",
        "- Output is an ICI-readiness atlas, not an ICI response predictor.",
        "",
        "## Paper 9 result summary",
        f"- Thyroid DepMap models: {p9_summary.get('thyroid_models')}",
        f"- Candidate dependency rows: {p9_summary.get('dep_rows')}",
        f"- Drug association rows: {p9_summary.get('drug_rows')}",
        "- Output is a cell-intrinsic vulnerability map, not clinical utility evidence.",
        "",
        "## Overlap genes/pathways",
        f"- Direct signature-target overlap: {', '.join(overlap) if overlap else 'none'}",
        "- Shared pathway themes: interferon/JAK-STAT context, dedifferentiation/stromal programs, chromatin and survival dependencies.",
        "",
        "## Difference",
        "- Paper 3 = immune readiness and tumor microenvironment state.",
        "- Paper 9 = cell-intrinsic dependency and drug-response nomination from cell-line data.",
        "",
        "## Claim guard",
        "- No treatment recommendation.",
        "- No patient treatment selection.",
        "- No ICI response predictor claim.",
        "- No TROP2 synthetic-lethality claim.",
        "",
        "## Next gates",
        "- Protected-access planning only if local credentials/files are explicitly supplied.",
        "- Wet-lab validation planning for Paper 9 candidates with non-high artifact risk.",
        "- Manuscript prose remains untouched.",
    ]
    (REPORT_ROOT / "06_p3_p9_execution_boundary_summary.md").write_text("\n".join(lines) + "\n")


def write_commit_proposal() -> None:
    lines = [
        "# Commit Proposal Only",
        "",
        "No automatic commit was made.",
        "",
        "## Proposed groups",
        "- P3-1 scripts and signature registry: `project/scripts/p3_p9_full_public_execution.py`, `project/results/p3_p9_full_execution/paper3/signature_registry.tsv`",
        "- P3-2 results and figures: Paper 3 TSV and PNG outputs under `project/results/p3_p9_full_execution/paper3/`",
        "- P3-3 report: `project/reports/p3_p9_full_execution/04_paper3_ici_full_public_report.md` and blocker report",
        "- P9-1 scripts and data manifest: execution script and `p9_public_download_manifest.tsv`",
        "- P9-2 processed scores/results: Paper 9 processed TSV outputs under `project/results/p3_p9_full_execution/paper9/` excluding raw",
        "- P9-3 figures: Paper 9 PNG outputs under `project/results/p3_p9_full_execution/paper9/`",
        "- P9-4 report: `project/reports/p3_p9_full_execution/05_paper9_synthetic_lethality_full_public_report.md`",
        "",
        "## Exclude",
        "- `project/results/p3_p9_full_execution/paper9/raw/`",
        "- Huge files",
        "- Protected data",
        "- Paper 1 manuscript files",
        "- Paper 2/4 files",
    ]
    (REPORT_ROOT / "07_commit_proposal.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    ensure_dirs()
    build_inventory_report()
    master = build_lineage_master()
    p3_summary = run_paper3(master)
    p9_summary = run_paper9()
    write_boundary_summary(p3_summary, p9_summary)
    write_commit_proposal()
    final = {
        "completed_at": datetime.now().isoformat(timespec="seconds"),
        "p3": p3_summary,
        "p9": {
            "thyroid_models": p9_summary.get("thyroid_models"),
            "dep_rows": p9_summary.get("dep_rows"),
            "drug_rows": p9_summary.get("drug_rows"),
            "top_targets": p9_summary.get("top_targets"),
        },
        "claim_guards": [
            "No protected access used",
            "No manuscript prose written",
            "No ICI response predictor claim",
            "No clinical treatment recommendation",
        ],
    }
    (REPORT_ROOT / "execution_summary.json").write_text(json.dumps(final, indent=2) + "\n")
    print(json.dumps(final, indent=2))


if __name__ == "__main__":
    main()
