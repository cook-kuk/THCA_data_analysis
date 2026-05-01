#!/usr/bin/env python
"""v17 Q20 — Mutation co-occurrence + mutual exclusivity meta-analysis.

Combines DNA-only thyroid cancer cohorts:
  - Wang 2023 Chinese (n=458)
  - Pozdeyev 2018 advanced thyroid (n=779)
  - Landa 2016 MSK PDTC/ATC (n=117)
Total ~1,354 patients.

Outputs co-occurrence tables, 4-phenotype classification, and 3 plotly figures.
Author: Seungho Cook
"""
from __future__ import annotations

import gzip
import json
import logging
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.stats import fisher_exact

RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)

PROJECT_ROOT = Path("/opt/thyroid-dash/project")
WANG_FILE = PROJECT_ROOT / "results/v17_wang2023/Wang2023_per_patient_flags.tsv"
POZ_FILE = PROJECT_ROOT / "results/v17_pozdeyev/Pozdeyev2018_per_patient.tsv"
MSK_MUT = PROJECT_ROOT / "results/v17_tert_recovery/v3/v3_thyroid_mskcc_2016_mutations.tsv.gz"
MSK_SAMPLE = PROJECT_ROOT / "results/v17_tert_recovery/v3/v3_thyroid_mskcc_2016_clinical_sample.tsv"

OUT_DIR = PROJECT_ROOT / "results/v17_dna_meta"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR = PROJECT_ROOT / "results/v17_dna_meta"
LOG_DIR = PROJECT_ROOT / "logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)
LOG_FILE = LOG_DIR / "v17_q20_dna_meta.log"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[logging.FileHandler(LOG_FILE, mode="w"), logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

DARK_BG = "#0b0e12"

GENES = [
    "BRAF_V600E",
    "RAS_any",
    "TERT_promoter",
    "TP53",
    "CDKN2A_loss",
    "CDKN2A_mut",
    "PIK3CA",
    "PTEN",
    "EIF1AX",
]

# Histology harmonization ----------------------------------------------------
HISTOLOGY_MAP_POZ = {
    "Papillary Thyroid Cancer": "PTC_advanced",
    "Pediatric Papillary Thyroid Cancer": "PTC_advanced",
    "Follicular Thyroid Cancer": "FTC",
    "Hurthle Cell Thyroid Cancer": "FTC",  # Hurthle (oncocytic) grouped with FTC
    "Anaplastic Thyroid Cancer": "ATC",
}
HISTOLOGY_MAP_WANG = {
    "PTC": "PTC",  # generic Chinese PTC
    "FTC": "FTC",
    "ATC": "ATC",
    "PDTC": "PDTC",
}
HISTOLOGY_MAP_MSK = {
    "Anaplastic Thyroid Carcinoma": "ATC",
    "Poorly Differentiated Thyroid Carcinoma": "PDTC",
}


def load_wang() -> pd.DataFrame:
    df = pd.read_csv(WANG_FILE, sep="\t")
    log.info("Wang raw n=%d cols=%s", len(df), list(df.columns))
    out = pd.DataFrame()
    out["patient_id"] = df["Patient ID"]
    out["histology"] = df["Histology"].map(HISTOLOGY_MAP_WANG).fillna(df["Histology"])
    out["BRAF_V600E"] = df["BRAF_V600E"].astype(int)
    out["RAS_any"] = df["RAS"].astype(int)
    out["TERT_promoter"] = df["TERTp"].astype(int)
    out["TP53"] = df["TP53"].astype(int)
    # Wang panel did not report these flags in the cached table -> NaN (excluded from per-pair tests)
    for g in ["CDKN2A_loss", "CDKN2A_mut", "PIK3CA", "PTEN", "EIF1AX"]:
        out[g] = np.nan
    out["cohort"] = "Wang2023"
    return out


def load_pozdeyev() -> pd.DataFrame:
    df = pd.read_csv(POZ_FILE, sep="\t")
    log.info("Pozdeyev raw n=%d", len(df))
    out = pd.DataFrame()
    out["patient_id"] = df["Specimen_Id"]
    out["histology"] = df["Tumor_Type"].map(HISTOLOGY_MAP_POZ).fillna(df["Tumor_Type"])
    out["BRAF_V600E"] = pd.to_numeric(df["BRAF_V600E"], errors="coerce").fillna(0).astype(int)
    out["RAS_any"] = pd.to_numeric(df["RAS_any"], errors="coerce").fillna(0).astype(int)
    out["TERT_promoter"] = pd.to_numeric(df["TERT_promoter"], errors="coerce").fillna(0).astype(int)
    out["TP53"] = pd.to_numeric(df["TP53"], errors="coerce").fillna(0).astype(int)
    out["CDKN2A_loss"] = pd.to_numeric(df["CDKN2A_loss"], errors="coerce").fillna(0).astype(int)
    # Pozdeyev table didn't expose these as flags -> NaN
    for g in ["CDKN2A_mut", "PIK3CA", "PTEN", "EIF1AX"]:
        out[g] = np.nan
    out["cohort"] = "Pozdeyev2018"
    return out


def load_landa_msk() -> pd.DataFrame:
    """Re-derive per-patient mutation flags from raw MSK 2016 MAF."""
    sample = pd.read_csv(MSK_SAMPLE, sep="\t")
    sample = sample[sample["CANCER_TYPE_DETAILED"].isin(HISTOLOGY_MAP_MSK)].copy()
    sample["histology"] = sample["CANCER_TYPE_DETAILED"].map(HISTOLOGY_MAP_MSK)
    keep = sample[["SAMPLE_ID", "histology"]].rename(columns={"SAMPLE_ID": "patient_id"})

    muts = pd.read_csv(MSK_MUT, sep="\t", low_memory=False)
    log.info("MSK mutations rows=%d, samples=%d", len(muts), len(keep))

    def has_braf_v600e(g: pd.DataFrame) -> int:
        sub = g[g["Hugo_Symbol"] == "BRAF"]
        if "HGVSp_Short" in sub.columns:
            return int(sub["HGVSp_Short"].astype(str).str.contains("V600E", na=False).any())
        return 0

    def has_ras(g: pd.DataFrame) -> int:
        sub = g[g["Hugo_Symbol"].isin(["NRAS", "HRAS", "KRAS"])]
        if sub.empty:
            return 0
        # treat any non-silent missense at codon 12/13/61 as RAS hot
        codons = sub["HGVSp_Short"].astype(str)
        is_hot = codons.str.contains(r"(?:G12|G13|Q61)", na=False, regex=True)
        return int(is_hot.any())

    def has_tert_promoter(g: pd.DataFrame) -> int:
        sub = g[(g["Hugo_Symbol"] == "TERT") & (g["Variant_Classification"] == "5'Flank")]
        return int(len(sub) > 0)

    def has_gene_nonsilent(g: pd.DataFrame, gene: str) -> int:
        sub = g[(g["Hugo_Symbol"] == gene) & (g["Variant_Classification"] != "Silent")]
        return int(len(sub) > 0)

    rows = []
    grouped = muts.groupby("Tumor_Sample_Barcode")
    for sid, hist in keep.itertuples(index=False):
        if sid in grouped.groups:
            sub = grouped.get_group(sid)
        else:
            sub = muts.iloc[0:0]
        rows.append(
            dict(
                patient_id=sid,
                histology=hist,
                BRAF_V600E=has_braf_v600e(sub),
                RAS_any=has_ras(sub),
                TERT_promoter=has_tert_promoter(sub),
                TP53=has_gene_nonsilent(sub, "TP53"),
                CDKN2A_loss=np.nan,  # CNA, not in MAF
                CDKN2A_mut=has_gene_nonsilent(sub, "CDKN2A"),
                PIK3CA=has_gene_nonsilent(sub, "PIK3CA"),
                PTEN=has_gene_nonsilent(sub, "PTEN"),
                EIF1AX=has_gene_nonsilent(sub, "EIF1AX"),
                cohort="Landa2016_MSK",
            )
        )
    return pd.DataFrame(rows)


def fisher_or(a: np.ndarray, b: np.ndarray):
    """Return (OR, p, n11, n10, n01, n00). Uses Haldane 0.5 correction for OR."""
    mask = (~pd.isna(a)) & (~pd.isna(b))
    a_, b_ = a[mask].astype(int), b[mask].astype(int)
    n11 = int(((a_ == 1) & (b_ == 1)).sum())
    n10 = int(((a_ == 1) & (b_ == 0)).sum())
    n01 = int(((a_ == 0) & (b_ == 1)).sum())
    n00 = int(((a_ == 0) & (b_ == 0)).sum())
    table = [[n11, n10], [n01, n00]]
    try:
        _, p = fisher_exact(table, alternative="two-sided")
    except ValueError:
        p = 1.0
    # Haldane-corrected OR for log-OR display
    a11, a10, a01, a00 = n11 + 0.5, n10 + 0.5, n01 + 0.5, n00 + 0.5
    or_ = (a11 * a00) / (a10 * a01)
    return or_, p, n11, n10, n01, n00


def pairwise_table(df: pd.DataFrame, genes: list[str], stratum_label: str) -> pd.DataFrame:
    rows = []
    for g1, g2 in combinations(genes, 2):
        if g1 not in df.columns or g2 not in df.columns:
            continue
        a = df[g1].to_numpy()
        b = df[g2].to_numpy()
        mask = (~pd.isna(a)) & (~pd.isna(b))
        if mask.sum() < 5:
            continue
        or_, p, n11, n10, n01, n00 = fisher_or(a, b)
        rows.append(
            dict(
                stratum=stratum_label,
                gene_a=g1,
                gene_b=g2,
                n_eval=int(mask.sum()),
                n11=n11,
                n10=n10,
                n01=n01,
                n00=n00,
                odds_ratio=or_,
                log2_or=np.log2(or_) if or_ > 0 else np.nan,
                fisher_p=p,
                relation="co_occur" if or_ > 1 else "mutually_exclusive",
            )
        )
    return pd.DataFrame(rows)


def classify_phenotype(row) -> str:
    flags = [int(row.get(g, 0) or 0) for g in ["BRAF_V600E", "RAS_any", "TERT_promoter", "TP53"]]
    n = sum(flags)
    if n == 0:
        return "Triple-negative"
    if n >= 2:
        return "Multi-hit"
    if flags[0] == 1:
        return "BRAF-only"
    if flags[1] == 1:
        return "RAS-only"
    if flags[2] == 1:
        return "TERT-only"
    if flags[3] == 1:
        return "TP53-only"
    return "Other"


def make_cooccurrence_heatmap(pairs: pd.DataFrame, path: Path) -> None:
    overall = pairs[pairs["stratum"] == "ALL"].copy()
    genes = sorted(set(overall["gene_a"]) | set(overall["gene_b"]))
    mat = pd.DataFrame(np.nan, index=genes, columns=genes, dtype=float)
    pmat = pd.DataFrame(np.nan, index=genes, columns=genes, dtype=float)
    for _, r in overall.iterrows():
        mat.loc[r.gene_a, r.gene_b] = r.log2_or
        mat.loc[r.gene_b, r.gene_a] = r.log2_or
        pmat.loc[r.gene_a, r.gene_b] = r.fisher_p
        pmat.loc[r.gene_b, r.gene_a] = r.fisher_p
    np.fill_diagonal(mat.values, 0.0)
    text = pmat.applymap(lambda p: "*" if pd.notna(p) and p < 0.001 else "")
    fig = go.Figure(
        data=go.Heatmap(
            z=mat.values,
            x=list(mat.columns),
            y=list(mat.index),
            colorscale="RdBu_r",
            zmid=0,
            zmin=-3,
            zmax=3,
            text=text.values,
            texttemplate="%{text}",
            colorbar=dict(title="log2(OR)"),
        )
    )
    fig.update_layout(
        title="Pairwise gene co-occurrence (combined n=1,354) — * Fisher p<0.001",
        template="plotly_dark",
        paper_bgcolor=DARK_BG,
        plot_bgcolor=DARK_BG,
        width=820,
        height=720,
    )
    fig.write_html(path)
    log.info("wrote %s", path)


def make_phenotype_bar(phen: pd.DataFrame, path: Path) -> None:
    grp = phen.groupby(["cohort", "histology", "phenotype"]).size().rename("n").reset_index()
    grp["total"] = grp.groupby(["cohort", "histology"])["n"].transform("sum")
    grp["pct"] = 100 * grp["n"] / grp["total"]
    grp["xlabel"] = grp["cohort"] + " | " + grp["histology"] + " (n=" + grp["total"].astype(str) + ")"
    order = [
        "BRAF-only",
        "RAS-only",
        "TERT-only",
        "TP53-only",
        "Multi-hit",
        "Triple-negative",
    ]
    palette = {
        "BRAF-only": "#ef553b",
        "RAS-only": "#00cc96",
        "TERT-only": "#ab63fa",
        "TP53-only": "#ffa15a",
        "Multi-hit": "#19d3f3",
        "Triple-negative": "#7f7f7f",
    }
    fig = go.Figure()
    xorder = sorted(grp["xlabel"].unique())
    for ph in order:
        sub = grp[grp["phenotype"] == ph]
        if sub.empty:
            continue
        fig.add_bar(
            name=ph,
            x=sub["xlabel"],
            y=sub["pct"],
            marker_color=palette[ph],
            customdata=np.stack([sub["n"], sub["total"]], axis=-1),
            hovertemplate="%{x}<br>" + ph + ": %{y:.1f}% (%{customdata[0]}/%{customdata[1]})<extra></extra>",
        )
    fig.update_layout(
        barmode="stack",
        title="4-group molecular phenotype by histology × cohort",
        template="plotly_dark",
        paper_bgcolor=DARK_BG,
        plot_bgcolor=DARK_BG,
        yaxis_title="% of patients",
        xaxis=dict(categoryorder="array", categoryarray=xorder, tickangle=-30),
        width=1100,
        height=620,
    )
    fig.write_html(path)
    log.info("wrote %s", path)


def make_tert_acquisition(combined: pd.DataFrame, path: Path) -> None:
    """Sankey: histology stage -> background driver -> TERT status."""
    df = combined.copy()
    # Use only cohorts with all 3 flags
    df = df.dropna(subset=["BRAF_V600E", "RAS_any", "TERT_promoter"]).copy()
    df["background"] = np.where(
        df["BRAF_V600E"] == 1,
        "BRAF",
        np.where(df["RAS_any"] == 1, "RAS", "Driver-neg"),
    )
    df["tert_state"] = np.where(df["TERT_promoter"] == 1, "TERT+", "TERT-")
    # Stage groups
    stage_map = {
        "PTC": "Indolent",
        "FTC": "Indolent",
        "PTC_advanced": "Advanced",
        "PDTC": "Advanced",
        "ATC": "Anaplastic",
    }
    df["stage"] = df["histology"].map(stage_map).fillna("Other")

    stage_levels = ["Indolent", "Advanced", "Anaplastic"]
    bg_levels = ["BRAF", "RAS", "Driver-neg"]
    tert_levels = ["TERT+", "TERT-"]

    nodes = stage_levels + [f"{s}|{b}" for s in stage_levels for b in bg_levels] + tert_levels
    node_idx = {n: i for i, n in enumerate(nodes)}

    src, tgt, val, color = [], [], [], []
    bg_color = {"BRAF": "#ef553b", "RAS": "#00cc96", "Driver-neg": "#7f7f7f"}
    tert_color = {"TERT+": "rgba(171,99,250,0.55)", "TERT-": "rgba(120,120,120,0.35)"}

    for s in stage_levels:
        for b in bg_levels:
            n = int(((df["stage"] == s) & (df["background"] == b)).sum())
            if n == 0:
                continue
            src.append(node_idx[s])
            tgt.append(node_idx[f"{s}|{b}"])
            val.append(n)
            color.append(bg_color[b].replace(")", ",0.45)").replace("rgb", "rgba") if bg_color[b].startswith("rgb") else "rgba(150,150,150,0.35)")
            for t in tert_levels:
                nt = int(
                    (
                        (df["stage"] == s)
                        & (df["background"] == b)
                        & (df["tert_state"] == t)
                    ).sum()
                )
                if nt == 0:
                    continue
                src.append(node_idx[f"{s}|{b}"])
                tgt.append(node_idx[t])
                val.append(nt)
                color.append(tert_color[t])

    node_colors = (
        ["#19d3f3"] * len(stage_levels)
        + [bg_color[b] for _ in stage_levels for b in bg_levels]
        + ["#ab63fa", "#7f7f7f"]
    )

    fig = go.Figure(
        data=[
            go.Sankey(
                node=dict(
                    label=nodes,
                    color=node_colors,
                    pad=18,
                    thickness=18,
                ),
                link=dict(source=src, target=tgt, value=val, color=color),
            )
        ]
    )
    fig.update_layout(
        title="TERT acquisition path: stage → background driver → TERT promoter status",
        template="plotly_dark",
        paper_bgcolor=DARK_BG,
        plot_bgcolor=DARK_BG,
        font=dict(color="#e6e6e6"),
        width=1080,
        height=640,
    )
    fig.write_html(path)
    log.info("wrote %s", path)


def main() -> None:
    log.info("=== v17 Q20 DNA meta-analysis (random_state=%d) ===", RANDOM_STATE)

    wang = load_wang()
    poz = load_pozdeyev()
    msk = load_landa_msk()
    log.info("Loaded: Wang n=%d, Pozdeyev n=%d, Landa-MSK n=%d", len(wang), len(poz), len(msk))

    combined = pd.concat([wang, poz, msk], ignore_index=True)
    combined["n"] = 1
    log.info("Combined n=%d", len(combined))

    combined.to_csv(OUT_DIR / "combined_3cohort_n1354.tsv", sep="\t", index=False)

    # Pairwise co-occurrence
    pair_frames = [pairwise_table(combined, GENES, "ALL")]
    for hist, sub in combined.groupby("histology"):
        if len(sub) < 10:
            continue
        pair_frames.append(pairwise_table(sub, GENES, f"hist={hist}"))
    for cohort, sub in combined.groupby("cohort"):
        pair_frames.append(pairwise_table(sub, GENES, f"cohort={cohort}"))
    pairs = pd.concat(pair_frames, ignore_index=True)
    pairs.to_csv(OUT_DIR / "q20_pairwise_cooccurrence.tsv", sep="\t", index=False)
    log.info("pairwise rows=%d", len(pairs))

    # 4-group phenotype
    combined["phenotype"] = combined.apply(classify_phenotype, axis=1)
    phen_summary = (
        combined.groupby(["cohort", "histology", "phenotype"]).size().rename("n").reset_index()
    )
    totals = combined.groupby(["cohort", "histology"]).size().rename("total").reset_index()
    phen_summary = phen_summary.merge(totals, on=["cohort", "histology"])
    phen_summary["pct"] = 100 * phen_summary["n"] / phen_summary["total"]
    phen_summary.to_csv(OUT_DIR / "q20_phenotype_by_histology.tsv", sep="\t", index=False)

    # Figures
    make_cooccurrence_heatmap(pairs, FIG_DIR / "v17_q20_cooccurrence_matrix.html")
    make_phenotype_bar(combined, FIG_DIR / "v17_q20_phenotype_by_histology.html")
    make_tert_acquisition(combined, FIG_DIR / "v17_q20_TERT_acquisition_path.html")

    # Key pair extraction for the JSON summary
    overall = pairs[pairs["stratum"] == "ALL"].set_index(["gene_a", "gene_b"])

    def lookup(g1, g2):
        if (g1, g2) in overall.index:
            r = overall.loc[(g1, g2)]
        elif (g2, g1) in overall.index:
            r = overall.loc[(g2, g1)]
        else:
            return None
        return dict(
            n_eval=int(r["n_eval"]),
            n11=int(r["n11"]),
            n10=int(r["n10"]),
            n01=int(r["n01"]),
            n00=int(r["n00"]),
            odds_ratio=float(r["odds_ratio"]),
            log2_or=float(r["log2_or"]),
            fisher_p=float(r["fisher_p"]),
            relation=str(r["relation"]),
        )

    key_pairs = {
        f"{a}__{b}": lookup(a, b)
        for a, b in [
            ("BRAF_V600E", "TERT_promoter"),
            ("BRAF_V600E", "RAS_any"),
            ("BRAF_V600E", "TP53"),
            ("TERT_promoter", "CDKN2A_loss"),
            ("RAS_any", "TERT_promoter"),
            ("RAS_any", "TP53"),
            ("RAS_any", "EIF1AX"),
            ("BRAF_V600E", "EIF1AX"),
            ("TP53", "TERT_promoter"),
            ("PIK3CA", "TP53"),
        ]
    }

    # Look for "novel" significant pairs (p<0.001) outside the focus list
    focus = {
        ("BRAF_V600E", "TERT_promoter"),
        ("BRAF_V600E", "RAS_any"),
        ("BRAF_V600E", "TP53"),
        ("TERT_promoter", "CDKN2A_loss"),
        ("RAS_any", "TERT_promoter"),
    }
    overall_reset = pairs[pairs["stratum"] == "ALL"].copy()
    sig = overall_reset[overall_reset["fisher_p"] < 0.001].copy()
    sig["pair_key"] = sig.apply(lambda r: tuple(sorted([r.gene_a, r.gene_b])), axis=1)
    novel = sig[~sig["pair_key"].apply(lambda x: x in {tuple(sorted(p)) for p in focus})]
    novel_records = novel.drop(columns=["pair_key"]).to_dict(orient="records")

    summary = dict(
        random_state=RANDOM_STATE,
        cohorts=dict(
            Wang2023=int((combined["cohort"] == "Wang2023").sum()),
            Pozdeyev2018=int((combined["cohort"] == "Pozdeyev2018").sum()),
            Landa2016_MSK=int((combined["cohort"] == "Landa2016_MSK").sum()),
            total=int(len(combined)),
        ),
        histology_counts={
            k: int(v) for k, v in combined["histology"].value_counts().items()
        },
        phenotype_counts={
            k: int(v) for k, v in combined["phenotype"].value_counts().items()
        },
        key_pairs=key_pairs,
        novel_significant_pairs=novel_records,
        files=dict(
            combined=str(OUT_DIR / "combined_3cohort_n1354.tsv"),
            pairwise=str(OUT_DIR / "q20_pairwise_cooccurrence.tsv"),
            phenotype=str(OUT_DIR / "q20_phenotype_by_histology.tsv"),
            fig_cooccurrence=str(FIG_DIR / "v17_q20_cooccurrence_matrix.html"),
            fig_phenotype=str(FIG_DIR / "v17_q20_phenotype_by_histology.html"),
            fig_tert_path=str(FIG_DIR / "v17_q20_TERT_acquisition_path.html"),
        ),
    )

    with open(OUT_DIR / "q20_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    log.info("=== done ===")


if __name__ == "__main__":
    main()
