"""
External expression validation — Paper 1 driver-orthogonal RAI / DM1 axis.

Steps 3-7 in one driver:
  - Parse GPL570 SOFT annotation (probe -> gene symbol)
  - Parse Series Matrix (metadata + probe x sample expression) for each dataset
  - Build sample_metadata.tsv (dataset, sample_id, histology_raw, histology_clean,
    disease_group, platform, modality, source_file, notes)
  - Collapse probe-level expression to gene-level (mean of probes per symbol)
  - Z-score within dataset, compute panel scores, run Mann-Whitney + Spearman

No raw FASTQ/CEL ingestion. Processed Series Matrix only. No cross-dataset pooling.
"""
from __future__ import annotations

import gzip
import io
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "raw"
OUT = ROOT
QC = ROOT / "qc"

# ---------------------------------------------------------------------------
# Gene panels (from spec)
# ---------------------------------------------------------------------------
RAI_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
THYROID_NONOVERLAP = [
    "SLC26A4", "IYD", "DUOX1", "DUOX2", "TFF3", "HHEX", "GLIS3", "DIO2",
]
TDS_TF = ["FOXE1", "NKX2-1", "PAX8", "HHEX"]
MECHANISM = ["DNMT1", "DNMT3B", "STAT3", "FOSL1", "JUNB", "TACSTD2"]
OPTIONAL_TARGETS = [
    "NAMPT", "KCNN4", "LYN", "OSMR", "IL6R", "SRC", "FYN",
    "ATR", "CHEK2", "MYC", "GLS", "LDHA",
]
ALL_GENES = sorted(set(RAI_8 + THYROID_NONOVERLAP + TDS_TF + MECHANISM + OPTIONAL_TARGETS))

# Affymetrix HG-U133+2 sometimes records NKX2-1 as TITF1; PAX8 has standard symbol;
# FOXE1 was historically TTF2 / FKHL15. Capture aliases for robust mapping.
GENE_ALIASES = {
    "NKX2-1": ["NKX2-1", "TITF1", "NKX2.1", "NKX2A"],
    "FOXE1": ["FOXE1", "TTF2", "FKHL15"],
    "SLC5A5": ["SLC5A5", "NIS"],
    "DUOX1": ["DUOX1"],
    "DUOX2": ["DUOX2"],
}

# ---------------------------------------------------------------------------
# Dataset registry
# ---------------------------------------------------------------------------
DATASETS = {
    "GSE33630":  {"platform": "GPL570",   "modality": "Affymetrix HG-U133+2 array"},
    "GSE29265":  {"platform": "GPL570",   "modality": "Affymetrix HG-U133+2 array"},
    "GSE65144":  {"platform": "GPL570",   "modality": "Affymetrix HG-U133+2 array"},
    "GSE53157":  {"platform": "GPL570",   "modality": "Affymetrix HG-U133+2 array"},
}

# ---------------------------------------------------------------------------
# GPL570 probe -> gene symbol map
# ---------------------------------------------------------------------------
def load_gpl570_probe2gene(path: Path) -> dict[str, list[str]]:
    """Return probe_id -> list of HGNC symbols (split on '///')."""
    probe2gene: dict[str, list[str]] = {}
    in_table = False
    sym_idx = None
    id_idx = None
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if not in_table:
                continue
            cols = line.rstrip("\n").split("\t")
            if id_idx is None:
                id_idx = cols.index("ID")
                sym_idx = cols.index("Gene symbol")
                continue
            probe = cols[id_idx]
            sym_field = cols[sym_idx] if sym_idx < len(cols) else ""
            if not sym_field:
                continue
            symbols = [s.strip() for s in sym_field.split("///") if s.strip()]
            if symbols:
                probe2gene[probe] = symbols
    return probe2gene


# ---------------------------------------------------------------------------
# Series Matrix parsing
# ---------------------------------------------------------------------------
def _quoted_cells(line: str) -> list[str]:
    # Tokens may be quoted with embedded tabs only between fields. Split on \t,
    # then strip surrounding quotes.
    parts = line.rstrip("\n").split("\t")
    return [p.strip().strip('"') for p in parts]


def parse_series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (metadata_df indexed by GSM, expression probe x sample DataFrame)."""
    meta_lines: list[list[str]] = []
    expr_lines: list[str] = []
    in_table = False
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                in_table = False
                continue
            if in_table:
                expr_lines.append(line)
            elif line.startswith("!Sample_"):
                meta_lines.append(_quoted_cells(line))

    # ---- metadata ----
    meta_dict: dict[str, list[str]] = {}
    char_counter: dict[str, int] = {}
    geo_acc = None
    for cells in meta_lines:
        key = cells[0].lstrip("!")
        vals = cells[1:]
        # multi-row characteristics_ch1 -> rename uniquely
        if key == "Sample_characteristics_ch1":
            char_counter[key] = char_counter.get(key, 0) + 1
            uniq = f"{key}_{char_counter[key]}"
            meta_dict[uniq] = vals
        else:
            meta_dict[key] = vals
    meta_df = pd.DataFrame(meta_dict)
    meta_df.index = meta_df["Sample_geo_accession"]
    meta_df.index.name = "sample_id"

    # ---- expression ----
    text = "".join(expr_lines)
    # First line is header: "ID_REF"\t"GSMxxx"\t...
    expr_df = pd.read_csv(
        io.StringIO(text),
        sep="\t",
        index_col=0,
        na_values=["", "null", "NA"],
        low_memory=False,
    )
    # Drop probe rows that are all-NA (microarray controls etc) or end-of-file marker
    expr_df = expr_df.dropna(how="all")
    # Remove quotes from column headers if present
    expr_df.columns = [c.strip('"') for c in expr_df.columns]
    # Strip quotes from probe IDs
    expr_df.index = expr_df.index.astype(str).str.strip().str.strip('"')
    return meta_df, expr_df


# ---------------------------------------------------------------------------
# Histology label cleaning per dataset
# ---------------------------------------------------------------------------
def _clean_label(raw: str) -> str:
    return raw.strip().lower()


def label_GSE33630(meta_df: pd.DataFrame) -> pd.DataFrame:
    src = meta_df.get("Sample_source_name_ch1", pd.Series(["", ] * len(meta_df), index=meta_df.index))
    char1 = meta_df.get("Sample_characteristics_ch1_1", pd.Series([""] * len(meta_df), index=meta_df.index))
    title = meta_df.get("Sample_title", pd.Series([""] * len(meta_df), index=meta_df.index))
    out = []
    for s, c, t in zip(src, char1, title):
        c_l = _clean_label(c); t_l = _clean_label(t); s_l = _clean_label(s)
        is_non_tumor = ("non-tumor" in s_l) or ("non-tumoral" in s_l) or ("non-tumor" in c_l)
        if is_non_tumor:
            out.append("normal")
        elif "anaplastic" in c_l or t_l.startswith("atc"):
            out.append("ATC")
        elif "papillary" in c_l or s_l == "tumour" or s_l == "tumor":
            out.append("PTC")
        elif "normal" in c_l or "normal" in s_l or "normal" in t_l:
            out.append("normal")
        else:
            out.append("UNKNOWN")
    return pd.Series(out, index=meta_df.index, name="histology_clean")


def label_GSE29265(meta_df: pd.DataFrame) -> pd.DataFrame:
    title = meta_df["Sample_title"].fillna("")
    src = meta_df.get("Sample_source_name_ch1", pd.Series([""] * len(meta_df), index=meta_df.index))
    out = []
    for t, s in zip(title, src):
        tl = _clean_label(t); sl = _clean_label(s)
        if "anaplastic" in tl or tl.startswith("atc"):
            out.append("ATC")
        elif "papillary" in tl or tl.startswith("ptc") or "ptc" in tl:
            # check it's not a "non-tumor" control
            if "non-tumor" in sl or "nontumor" in sl or "patient-matched" in sl:
                out.append("normal")
            else:
                out.append("PTC")
        elif "non-tumor" in tl or "nontumor" in tl or "patient-matched" in sl:
            out.append("normal")
        else:
            out.append("UNKNOWN")
    return pd.Series(out, index=meta_df.index, name="histology_clean")


def label_GSE65144(meta_df: pd.DataFrame) -> pd.DataFrame:
    src = meta_df["Sample_source_name_ch1"].fillna("")
    out = []
    for s in src:
        sl = _clean_label(s)
        if "anaplastic" in sl or "atc" in sl:
            out.append("ATC")
        elif "normal" in sl:
            out.append("normal")
        else:
            out.append("UNKNOWN")
    return pd.Series(out, index=meta_df.index, name="histology_clean")


def label_GSE53157(meta_df: pd.DataFrame) -> pd.DataFrame:
    src = meta_df["Sample_source_name_ch1"].fillna("")
    title = meta_df["Sample_title"].fillna("")
    out = []
    for s, t in zip(src, title):
        sl = _clean_label(s); tl = _clean_label(t)
        # Drop commercial pool
        if "commercial" in sl or "rna pool" in tl or "pool" in tl:
            out.append("DROP_pool")
        elif "pdtc" in sl or tl.startswith("pdtc"):
            out.append("PDTC")
        elif "fvptc" in sl or tl.startswith("fvptc"):
            out.append("FVPTC")
        elif "ftc" in sl or tl.startswith("ftc"):
            out.append("FTC")
        elif "ptc" in sl or "papillary" in sl or tl.startswith("ptc"):
            out.append("PTC")
        elif "normal" in sl or "normal" in tl:
            out.append("normal")
        else:
            out.append("UNKNOWN")
    return pd.Series(out, index=meta_df.index, name="histology_clean")


LABEL_FUNCS = {
    "GSE33630": label_GSE33630,
    "GSE29265": label_GSE29265,
    "GSE65144": label_GSE65144,
    "GSE53157": label_GSE53157,
}


# ---------------------------------------------------------------------------
# Probe -> gene collapse
# ---------------------------------------------------------------------------
def probe_to_gene(expr_probe: pd.DataFrame, probe2gene: dict[str, list[str]]) -> pd.DataFrame:
    """Collapse probe x sample matrix to gene x sample using mean across probes per symbol."""
    # Detect log-scale: GPL570 series matrices already RMA log2, but verify.
    sample_min = float(np.nanmin(expr_probe.values))
    sample_max = float(np.nanmax(expr_probe.values))
    if sample_max > 50:  # likely raw scale
        expr_probe = np.log2(expr_probe.clip(lower=1.0))
    rows: dict[str, list[np.ndarray]] = {}
    for probe, vals in expr_probe.iterrows():
        symbols = probe2gene.get(probe)
        if not symbols:
            continue
        v = vals.to_numpy(dtype=float)
        for sym in symbols:
            rows.setdefault(sym, []).append(v)
    gene_data = {sym: np.nanmean(np.vstack(stack), axis=0) for sym, stack in rows.items()}
    gene_df = pd.DataFrame(gene_data, index=expr_probe.columns).T
    gene_df.index.name = "gene_symbol"
    return gene_df


def lookup_gene(gene_df: pd.DataFrame, target: str) -> pd.Series | None:
    """Try canonical symbol then aliases."""
    candidates = [target] + [a for a in GENE_ALIASES.get(target, []) if a != target]
    for cand in candidates:
        if cand in gene_df.index:
            return gene_df.loc[cand]
    return None


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------
def zscore_within_dataset(gene_df: pd.DataFrame) -> pd.DataFrame:
    arr = gene_df.to_numpy(dtype=float)
    mu = np.nanmean(arr, axis=1, keepdims=True)
    sd = np.nanstd(arr, axis=1, keepdims=True, ddof=0)
    sd[sd == 0] = np.nan
    z = (arr - mu) / sd
    return pd.DataFrame(z, index=gene_df.index, columns=gene_df.columns)


def panel_score(z_df: pd.DataFrame, panel: list[str]) -> tuple[pd.Series, list[str], list[str]]:
    found = []
    missing = []
    series_list = []
    for g in panel:
        s = lookup_gene(z_df, g)
        if s is None:
            missing.append(g)
        else:
            found.append(g)
            series_list.append(s.values)
    if not series_list:
        return pd.Series(np.nan, index=z_df.columns), found, missing
    score = np.nanmean(np.vstack(series_list), axis=0)
    return pd.Series(score, index=z_df.columns), found, missing


def cohen_d(a: np.ndarray, b: np.ndarray) -> float:
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    ma, mb = a.mean(), b.mean()
    va, vb = a.var(ddof=1), b.var(ddof=1)
    sp = np.sqrt(((len(a) - 1) * va + (len(b) - 1) * vb) / (len(a) + len(b) - 2))
    if sp == 0:
        return float("nan")
    return float((ma - mb) / sp)


def mw_test(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan"), float("nan")
    res = mannwhitneyu(a, b, alternative="two-sided")
    return float(res.statistic), float(res.pvalue)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    QC.mkdir(parents=True, exist_ok=True)

    print("[1/4] loading GPL570 annotation...")
    probe2gene = load_gpl570_probe2gene(RAW / "GPL570" / "GPL570.annot.gz")
    print(f"      {len(probe2gene):,} probes mapped")

    all_meta_rows: list[dict] = []
    coverage_rows: list[dict] = []
    score_rows: list[dict] = []
    test_rows: list[dict] = []
    direction_rows: list[dict] = []

    panel_definitions = {
        "RAI_8": RAI_8,
        "THYROID_NONOVERLAP": THYROID_NONOVERLAP,
        "TDS_TF": TDS_TF,
        "TF_collapse": TDS_TF,
        "STAT3_AP1_DNMT": ["STAT3", "FOSL1", "JUNB", "DNMT1", "DNMT3B"],
    }

    sample_score_dump: list[pd.DataFrame] = []

    for acc, info in DATASETS.items():
        print(f"\n=== {acc} ===")
        sm_path = RAW / acc / f"{acc}_series_matrix.txt.gz"
        meta_df, probe_expr = parse_series_matrix(sm_path)
        print(f"  meta n={len(meta_df)}  probe matrix shape={probe_expr.shape}")

        # ---- histology labels ----
        clean = LABEL_FUNCS[acc](meta_df)
        # Source field for histology_raw
        src_field = meta_df.get("Sample_source_name_ch1", meta_df["Sample_title"]).astype(str)
        title_field = meta_df["Sample_title"].astype(str)

        # ---- gene-level collapse ----
        gene_df = probe_to_gene(probe_expr, probe2gene)
        print(f"  gene matrix shape={gene_df.shape}")
        # write gz expression
        out_path = OUT / f"{acc}_expression_gene_log.tsv.gz"
        gene_df.to_csv(out_path, sep="\t", compression="gzip")

        # ---- coverage ----
        for panel_name, panel_genes in {
            "RAI_8": RAI_8,
            "THYROID_NONOVERLAP": THYROID_NONOVERLAP,
            "TDS_TF": TDS_TF,
            "MECHANISM": MECHANISM,
            "OPTIONAL_TARGETS": OPTIONAL_TARGETS,
        }.items():
            for g in panel_genes:
                row = lookup_gene(gene_df, g)
                coverage_rows.append({
                    "dataset": acc,
                    "panel": panel_name,
                    "gene": g,
                    "available": int(row is not None),
                    "matched_symbol": (
                        next((cand for cand in [g] + GENE_ALIASES.get(g, []) if cand in gene_df.index), "")
                        if row is not None else ""
                    ),
                })

        # ---- z-scores ----
        z_df = zscore_within_dataset(gene_df)

        # ---- panel scores ----
        rai_score, rai_found, rai_missing = panel_score(z_df, RAI_8)
        non_score, non_found, non_missing = panel_score(z_df, THYROID_NONOVERLAP)
        tds_score, tds_found, tds_missing = panel_score(z_df, TDS_TF)
        tf_score = tds_score
        sad_score, sad_found, sad_missing = panel_score(z_df, ["STAT3", "FOSL1", "JUNB", "DNMT1", "DNMT3B"])
        tacstd2 = lookup_gene(z_df, "TACSTD2")
        tacstd2_z = pd.Series(tacstd2.values if tacstd2 is not None else [np.nan] * z_df.shape[1],
                              index=z_df.columns)

        dm1_like = -rai_score

        sample_score = pd.DataFrame({
            "dataset": acc,
            "sample_id": z_df.columns,
            "histology_clean": clean.reindex(z_df.columns).values,
            "histology_raw": src_field.reindex(z_df.columns).values,
            "title": title_field.reindex(z_df.columns).values,
            "RAI_8_score": rai_score.values,
            "DM1_like_score": dm1_like.values,
            "THYROID_NONOVERLAP_score": non_score.values,
            "TDS_like_score": tds_score.values,
            "TF_collapse_score": tf_score.values,
            "STAT3_AP1_DNMT_score": sad_score.values,
            "TACSTD2_z": tacstd2_z.values,
        })
        sample_score_dump.append(sample_score)
        print(f"  RAI_8 found={len(rai_found)}/8 missing={rai_missing}")
        print(f"  NONOV found={len(non_found)}/8 missing={non_missing}")

        # ---- metadata table ----
        for sid in z_df.columns:
            cl = clean.get(sid, "UNKNOWN")
            disease_group = cl
            if cl in {"DROP_pool", "UNKNOWN"}:
                disease_group = cl
            all_meta_rows.append({
                "dataset": acc,
                "sample_id": sid,
                "histology_raw": src_field.get(sid, ""),
                "title": title_field.get(sid, ""),
                "histology_clean": cl,
                "disease_group": disease_group,
                "platform": info["platform"],
                "modality": info["modality"],
                "source_file": str(sm_path.relative_to(ROOT.parent)),
                "notes": "",
            })

        # ---- contrasts ----
        groups_present = sample_score.groupby("histology_clean")
        avail = {k: g for k, g in groups_present if k not in {"DROP_pool", "UNKNOWN"}}
        contrasts = []
        if "ATC" in avail and "PTC" in avail: contrasts.append(("ATC", "PTC"))
        if "ATC" in avail and "normal" in avail: contrasts.append(("ATC", "normal"))
        if "ATC" in avail and "PDTC" in avail: contrasts.append(("ATC", "PDTC"))
        if "PDTC" in avail and "PTC" in avail: contrasts.append(("PDTC", "PTC"))
        if "PTC" in avail and "normal" in avail: contrasts.append(("PTC", "normal"))
        if "PDTC" in avail and "normal" in avail: contrasts.append(("PDTC", "normal"))
        if "FVPTC" in avail and "PTC" in avail: contrasts.append(("FVPTC", "PTC"))
        if "FTC" in avail and "normal" in avail: contrasts.append(("FTC", "normal"))

        score_columns = [
            "RAI_8_score", "DM1_like_score", "THYROID_NONOVERLAP_score",
            "TDS_like_score", "TF_collapse_score", "STAT3_AP1_DNMT_score", "TACSTD2_z",
        ]
        # store summary row per (dataset, contrast)
        for g1, g2 in contrasts:
            a_df = avail[g1]; b_df = avail[g2]
            row = {"dataset": acc, "contrast": f"{g1}_vs_{g2}",
                   "n_group1": len(a_df), "n_group2": len(b_df)}
            for col in score_columns:
                a = a_df[col].to_numpy(dtype=float)
                b = b_df[col].to_numpy(dtype=float)
                d = cohen_d(a, b)
                stat, p = mw_test(a, b)
                row[f"{col}_d"] = d
                row[f"{col}_p"] = p
            test_rows.append(row)

        # ---- direction consistency: signs of advanced-vs-DTC RAI delta ----
        # Define "advanced" >= ATC || PDTC; "DTC" PTC, FVPTC, FTC; "normal" reference
        adv = pd.concat([avail.get(k, pd.DataFrame()) for k in ("ATC", "PDTC")])
        dtc = pd.concat([avail.get(k, pd.DataFrame()) for k in ("PTC", "FVPTC", "FTC")])
        norm = avail.get("normal", pd.DataFrame())
        if len(adv) >= 2 and len(dtc) >= 2:
            for col in score_columns:
                d = cohen_d(adv[col].to_numpy(dtype=float), dtc[col].to_numpy(dtype=float))
                _, p = mw_test(adv[col].to_numpy(dtype=float), dtc[col].to_numpy(dtype=float))
                direction_rows.append({"dataset": acc, "contrast": "advanced_vs_DTC",
                                       "score": col, "cohen_d": d, "p": p,
                                       "n1": len(adv), "n2": len(dtc)})
        if len(adv) >= 2 and len(norm) >= 2:
            for col in score_columns:
                d = cohen_d(adv[col].to_numpy(dtype=float), norm[col].to_numpy(dtype=float))
                _, p = mw_test(adv[col].to_numpy(dtype=float), norm[col].to_numpy(dtype=float))
                direction_rows.append({"dataset": acc, "contrast": "advanced_vs_normal",
                                       "score": col, "cohen_d": d, "p": p,
                                       "n1": len(adv), "n2": len(norm)})

        # ---- spearman correlations ----
        keep = sample_score[~sample_score["histology_clean"].isin(["DROP_pool", "UNKNOWN"])]
        if len(keep) >= 5:
            for x, y in [
                ("DM1_like_score", "THYROID_NONOVERLAP_score"),
                ("TF_collapse_score", "DM1_like_score"),
                ("STAT3_AP1_DNMT_score", "DM1_like_score"),
                ("TACSTD2_z", "DM1_like_score"),
            ]:
                if keep[x].notna().sum() >= 3 and keep[y].notna().sum() >= 3:
                    rho, p = spearmanr(keep[x], keep[y], nan_policy="omit")
                    score_rows.append({
                        "dataset": acc,
                        "x": x, "y": y,
                        "rho": float(rho), "p": float(p), "n": int(keep[x].notna().sum()),
                    })

    # ---- write tables ----
    sample_metadata_df = pd.DataFrame(all_meta_rows)
    sample_metadata_df.to_csv(OUT / "sample_metadata.tsv", sep="\t", index=False)

    coverage_df = pd.DataFrame(coverage_rows)
    coverage_df.to_csv(OUT / "external_gene_coverage.tsv", sep="\t", index=False)

    test_df = pd.DataFrame(test_rows)
    test_df.to_csv(OUT / "external_score_tests.tsv", sep="\t", index=False)

    direction_df = pd.DataFrame(direction_rows)
    direction_df.to_csv(OUT / "external_direction_consistency.tsv", sep="\t", index=False)

    spearman_df = pd.DataFrame(score_rows)
    spearman_df.to_csv(OUT / "external_spearman.tsv", sep="\t", index=False)

    full_scores = pd.concat(sample_score_dump, ignore_index=True)
    full_scores.to_csv(OUT / "external_sample_scores.tsv.gz", sep="\t", index=False, compression="gzip")

    print("\n[done] tables written:")
    for f in [
        "sample_metadata.tsv", "external_gene_coverage.tsv",
        "external_score_tests.tsv", "external_direction_consistency.tsv",
        "external_spearman.tsv", "external_sample_scores.tsv.gz",
    ]:
        p = OUT / f
        print(f"  {p.relative_to(ROOT.parent)}  ({p.stat().st_size:,} B)")


if __name__ == "__main__":
    main()
