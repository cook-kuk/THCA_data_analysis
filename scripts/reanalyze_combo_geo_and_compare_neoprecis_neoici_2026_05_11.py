#!/usr/bin/env python3
from __future__ import annotations

import gzip
import json
import math
import re
import shutil
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import scanpy as sc


ROOT = Path(__file__).resolve().parents[1]
TMP = Path("/data/thca/_tmp/neoici_combo_reanalysis")
RAW = TMP / "raw"
OUT = ROOT / "project" / "results" / "p_neoprecis_vs_neoici_combo_reanalysis_2026_05_11"
HUB = ROOT / "project" / "papers_hub_2026_05_04"
WWW = Path("/var/www/papers/papers_hub_2026_05_04")

GSE = {
    "GSE222011": {
        "kind": "GEO/SRA",
        "processed": "open processed files; SRA raw available",
        "title": "autogene cevumeran + atezolizumab + mFOLFIRINOX",
        "cancer": "PDAC",
        "claim": "positive-control longitudinal immune-state and TCR-expansion lane",
    },
    "GSE255830": {
        "kind": "GEO/dbGaP",
        "processed": "open processed GEO files; raw data withheld to dbGaP",
        "title": "GNOS-PV02 + INO-9012 + pembrolizumab",
        "cancer": "HCC",
        "claim": "post-vaccine TCR clonotype and phenotype validation lane",
    },
}

CYTOTOXIC = ["NKG7", "GNLY", "PRF1", "GZMB", "CCL5", "IFNG"]
EXHAUST = ["PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "CTLA4"]
MEMORY = ["IL7R", "LTB", "CCR7", "SELL", "TCF7"]
ACTIVATION = ["CD69", "IL2RA", "HLA-DRA", "ICOS", "TNFRSF9", "TNFRSF4"]
PROLIF = ["MKI67", "TOP2A", "STMN1", "TYMS"]


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def ftp_text(url: str) -> str:
    return requests.get(url, timeout=120).text


def download(url: str, dest: Path) -> None:
    ensure_dir(dest.parent)
    if dest.exists() and dest.stat().st_size > 0:
        return
    with requests.get(url, stream=True, timeout=300) as r:
        r.raise_for_status()
        with dest.open("wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                if chunk:
                    f.write(chunk)


def parse_filelist(acc: str) -> pd.DataFrame:
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/filelist.txt"
    txt = ftp_text(url)
    rows = []
    for line in txt.splitlines():
        if not line.startswith("File\t"):
            continue
        _, name, time, size, typ = line.split("\t")
        rows.append({"name": name, "time": time, "size": int(size), "type": typ})
    return pd.DataFrame(rows)


def sample_title_map(acc: str) -> dict[str, str]:
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/matrix/{acc}_series_matrix.txt.gz"
    txt = gzip.decompress(requests.get(url, timeout=120).content).decode("utf-8", errors="replace")
    lines = txt.splitlines()
    titles_line = next(l for l in lines if l.startswith("!Sample_title"))
    access_line = next(l for l in lines if l.startswith("!Sample_geo_accession"))
    titles = [t.strip('"') for t in titles_line.split("\t")[1:]]
    accessions = [a.strip('"') for a in access_line.split("\t")[1:]]
    return dict(zip(accessions, titles))


def split_group(title: str) -> str:
    m = re.search(r"_(mFOL1|mFOL2|followup4|followup6|followup8|presurgery|atezo|pcv1|pcv5)", title)
    if m:
        return m.group(1)
    return title


def load_10x_h5(path: Path) -> sc.AnnData:
    adata = sc.read_10x_h5(path)
    adata.var_names_make_unique()
    return adata


def load_10x_mtx(path: Path) -> sc.AnnData:
    adata = sc.read_10x_mtx(path, var_names="gene_symbols", make_unique=True)
    adata.var_names_make_unique()
    return adata


def score_module(adata: sc.AnnData, name: str, genes: Iterable[str]) -> None:
    genes = [g for g in genes if g in adata.var_names]
    if len(genes) >= 3:
        sc.tl.score_genes(adata, gene_list=genes, score_name=name, use_raw=False)
    else:
        adata.obs[name] = np.nan


def entropy(values: Iterable[str]) -> float:
    c = Counter(values)
    n = sum(c.values())
    if n <= 1:
        return 0.0
    probs = np.array([v / n for v in c.values()], dtype=float)
    return float(-(probs * np.log2(probs)).sum())


def gini(values: Iterable[int]) -> float:
    arr = np.array(list(values), dtype=float)
    if len(arr) == 0:
        return np.nan
    arr = np.sort(arr)
    n = len(arr)
    if arr.sum() == 0:
        return 0.0
    cum = np.cumsum(arr)
    return float((n + 1 - 2 * (cum / cum[-1]).sum()) / n)


def tcr_summary(contig_path: Path) -> dict[str, float]:
    df = pd.read_csv(contig_path)
    if "barcode" not in df.columns:
        return {"tcr_cells": 0, "tcr_clonotypes": 0, "top_clone_frac": np.nan, "tcr_entropy": np.nan, "tcr_gini": np.nan}
    if "productive" in df.columns:
        prod_mask = df["productive"].astype(str).str.lower().isin({"true", "1", "yes"})
        productive = df[prod_mask].copy()
    else:
        productive = df.copy()
    if productive.empty:
        return {"tcr_cells": 0, "tcr_clonotypes": 0, "top_clone_frac": np.nan, "tcr_entropy": np.nan, "tcr_gini": np.nan}
    clonotype_col = None
    for col in ["raw_clonotype_id", "clonotype_id", "clone_id", "cdr3"]:
        if col in productive.columns:
            clonotype_col = col
            break
    if clonotype_col is None:
        clonotype_col = "barcode"
    counts = productive[clonotype_col].astype(str).value_counts()
    return {
        "tcr_cells": float(productive.shape[0]),
        "tcr_clonotypes": float(counts.shape[0]),
        "top_clone_frac": float(counts.iloc[0] / counts.sum()),
        "tcr_entropy": float(entropy(productive[clonotype_col].astype(str))),
        "tcr_gini": float(gini(counts.values)),
    }


def download_dataset(acc: str) -> dict[str, dict[str, Path | str]]:
    fl = parse_filelist(acc)
    out: dict[str, dict[str, Path | str]] = {}
    for _, row in fl.iterrows():
        if row["type"] not in {"H5", "CSV", "MTX", "TSV"}:
            continue
        name = row["name"]
        if acc == "GSE222011" and not (name.endswith(".h5") or name.endswith(".csv.gz")):
            continue
        if acc == "GSE255830" and "readme" in name.lower():
            continue
        gsm = name.split("_", 1)[0]
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/{name}"
        label = sample_label_from_title(sample_title_map(acc).get(gsm, gsm))
        if acc == "GSE255830" and (name.endswith(".mtx.gz") or name.endswith("_barcodes.tsv.gz") or name.endswith("_features.tsv.gz")):
            if name.endswith("_matrix.mtx.gz"):
                canon = "matrix.mtx.gz"
            elif name.endswith("_barcodes.tsv.gz"):
                canon = "barcodes.tsv.gz"
            elif name.endswith("_features.tsv.gz"):
                canon = "features.tsv.gz"
            else:
                canon = name
            dest = RAW / acc / label / canon
            out.setdefault(label, {})["gex_dir"] = dest.parent
        elif acc == "GSE255830" and name.endswith(".csv.gz"):
            dest = RAW / acc / label / name
            out.setdefault(label, {})["tcr"] = dest
        elif acc == "GSE222011" and name.endswith(".h5"):
            dest = RAW / acc / label / name
            out.setdefault(label, {})["gex"] = dest
        elif acc == "GSE222011" and name.endswith(".csv.gz"):
            dest = RAW / acc / label / name
            out.setdefault(label, {})["tcr"] = dest
        else:
            continue
        download(url, dest)
    return out


def normalize_title(title: str) -> str:
    title = title.strip()
    for suffix in ["-GEX", "-TCR", ", scRNA", ", TCR-seq", ", scRNAseq", ", scRNA-seq"]:
        if title.endswith(suffix):
            title = title[: -len(suffix)]
    return title.strip()


def sample_label_from_title(title: str) -> str:
    return normalize_title(title)


def analyze_gse222011(paths: dict[str, dict[str, Path | str]]) -> pd.DataFrame:
    rows = []
    for label, d in paths.items():
        if "gex" not in d or "tcr" not in d:
            continue
        adata = load_10x_h5(Path(d["gex"]))
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        for nm, genes in [("cytotoxic", CYTOTOXIC), ("exhaustion", EXHAUST), ("memory", MEMORY), ("activation", ACTIVATION), ("proliferation", PROLIF)]:
            score_module(adata, nm, genes)
        df = adata.obs.copy()
        tcr = pd.read_csv(Path(d["tcr"]))
        tcr["barcode"] = tcr["barcode"].astype(str).str.replace(r"-1$", "", regex=True)
        # map TCR clonotypes to barcodes from the same library
        tcr_stats = tcr_summary(Path(d["tcr"]))
        rows.append({
            "dataset": "GSE222011",
            "sample": label,
            "title": label,
            "n_cells": int(adata.n_obs),
            "cytotoxic_score": float(np.nanmean(df["cytotoxic"])),
            "exhaustion_score": float(np.nanmean(df["exhaustion"])),
            "memory_score": float(np.nanmean(df["memory"])),
            "activation_score": float(np.nanmean(df["activation"])),
            "proliferation_score": float(np.nanmean(df["proliferation"])),
            "response_index": float(np.nanmean(df["cytotoxic"]) - np.nanmean(df["exhaustion"]) + 0.5 * np.nanmean(df["activation"])),
            **tcr_stats,
        })
    return pd.DataFrame(rows)


def analyze_gse255830(paths: dict[str, dict[str, Path | str]]) -> pd.DataFrame:
    rows = []
    for label, d in paths.items():
        if (("gex" not in d) and ("gex_dir" not in d)) or "tcr" not in d:
            continue
        if "gex_dir" in d:
            adata = load_10x_mtx(Path(d["gex_dir"]))
        else:
            adata = load_10x_h5(Path(d["gex"]))
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        for nm, genes in [("cytotoxic", CYTOTOXIC), ("exhaustion", EXHAUST), ("memory", MEMORY), ("activation", ACTIVATION), ("proliferation", PROLIF)]:
            score_module(adata, nm, genes)
        df = adata.obs.copy()
        tcr_stats = tcr_summary(Path(d["tcr"]))
        rows.append({
            "dataset": "GSE255830",
            "sample": label,
            "title": label,
            "n_cells": int(adata.n_obs),
            "cytotoxic_score": float(np.nanmean(df["cytotoxic"])),
            "exhaustion_score": float(np.nanmean(df["exhaustion"])),
            "memory_score": float(np.nanmean(df["memory"])),
            "activation_score": float(np.nanmean(df["activation"])),
            "proliferation_score": float(np.nanmean(df["proliferation"])),
            "response_index": float(np.nanmean(df["cytotoxic"]) - np.nanmean(df["exhaustion"]) + 0.5 * np.nanmean(df["activation"])),
            **tcr_stats,
        })
    return pd.DataFrame(rows)


def plot_panels(df: pd.DataFrame) -> list[Path]:
    paths = []
    for metric in ["response_index", "cytotoxic_score", "exhaustion_score", "top_clone_frac"]:
        fig, ax = plt.subplots(figsize=(9, 4), dpi=160)
        plot_df = df.copy()
        plot_df["label"] = plot_df["sample"]
        if metric == "top_clone_frac":
            plot_df[metric] = plot_df[metric].fillna(0)
        for ds, sub in plot_df.groupby("dataset"):
            ax.plot(sub["sample"], sub[metric], marker="o", label=ds)
        ax.set_xticklabels(plot_df["sample"], rotation=35, ha="right")
        ax.set_title(metric)
        ax.legend(frameon=False)
        fig.tight_layout()
        out = OUT / "figures" / f"Fig_combo_{metric}.png"
        fig.savefig(out)
        plt.close(fig)
        paths.append(out)
    return paths


def compare_contracts() -> pd.DataFrame:
    neoici = pd.read_csv(ROOT / "project" / "results" / "p_neo_ici_combo_2026_05_11" / "neoici_algorithm_gauntlet_metrics.tsv", sep="\t")
    rows = [
        {
            "system": "NeoPrecis",
            "type": "published contract",
            "best_claim": "Immunogenicity + clonality-aware landscape improves ICI response prediction and has CEDAR/NCI/melanoma/NSCLC benchmarks.",
            "current_status": "published; not locally paired to the combo cohorts in this workspace",
        },
        {
            "system": "NeoICI_GA_RL_combo_v1",
            "type": "local combo prioritizer",
            "best_claim": "experiment-priority combo-therapy layer",
            "current_status": f"academic_validation_like AUPRC {neoici.query('split == \"academic_validation_like\" and algorithm == \"NeoICI_GA_RL_combo_v1\"')['AUPRC'].iloc[0]:.3f}; industrial_locked_v0 AUPRC {neoici.query('split == \"industrial_locked_v0\" and algorithm == \"NeoICI_GA_RL_combo_v1\"')['AUPRC'].iloc[0]:.3f}",
        },
        {
            "system": "KG_GA_evolved",
            "type": "local champion",
            "best_claim": "retrospective validation-like champion",
            "current_status": f"academic_validation_like AUPRC {neoici.query('split == \"academic_validation_like\" and algorithm == \"KG_GA_evolved\"')['AUPRC'].iloc[0]:.3f}",
        },
        {
            "system": "BioDarwin_mode_bank_v4",
            "type": "industrial anchor",
            "best_claim": "locked public mini-set anchor",
            "current_status": f"industrial_locked_v0 AUPRC {neoici.query('split == \"industrial_locked_v0\" and algorithm == \"BioDarwin_mode_bank_v4\"')['AUPRC'].iloc[0]:.3f}",
        },
    ]
    return pd.DataFrame(rows)


def main() -> None:
    ensure_dir(RAW)
    ensure_dir(OUT / "figures")
    ensure_dir(HUB)
    ensure_dir(WWW)

    g222 = download_dataset("GSE222011")
    g255 = download_dataset("GSE255830")

    df222 = analyze_gse222011(g222)
    df255 = analyze_gse255830(g255)
    combined = pd.concat([df222, df255], axis=0, ignore_index=True)
    combined.to_csv(OUT / "combo_geo_reanalysis.tsv", sep="\t", index=False)

    contract = compare_contracts()
    contract.to_csv(OUT / "neoprecis_neoici_contract.tsv", sep="\t", index=False)

    figs = plot_panels(combined)

    summary = [
        "# Open combo GEO reanalysis",
        "",
        "## Datasets",
        "",
        f"- GSE222011 rows: {df222.shape[0]}",
        f"- GSE255830 rows: {df255.shape[0]}",
        "",
        "## What was measured",
        "",
        "- Per-sample GEX module scores: cytotoxicity, exhaustion, memory, activation, proliferation",
        "- TCR contig summary: productive TCR cell count, clonotype count, top clone fraction, entropy, gini",
        "",
        "## Early readout",
        "",
        f"- Combined response-index median: {combined['response_index'].median():.3f}",
        f"- Combined top-clone-fraction median: {combined['top_clone_frac'].median():.3f}",
        "",
        "## Contract comparison",
        "",
        "- NeoPrecis remains the published benchmark framework.",
        "- NeoICI remains the broader combo-therapy prioritization layer, but no same-dataset paired score table is available yet.",
    ]
    (OUT / "COMBO_GEO_REANALYSIS_SUMMARY.md").write_text("\n".join(summary) + "\n", encoding="utf-8")

    html = [
        "<!doctype html><html lang='en'><head><meta charset='utf-8'/>",
        "<meta name='viewport' content='width=device-width, initial-scale=1'/>",
        "<title>Open Combo GEO Reanalysis</title>",
        "<style>body{margin:0;background:#0b1020;color:#e5e7eb;font-family:system-ui,sans-serif}.w{max-width:1200px;margin:0 auto;padding:28px 22px 40px}table{width:100%;border-collapse:collapse;font-size:14px}th,td{border:1px solid #243047;padding:8px;vertical-align:top}th{background:#111827}.g{display:grid;grid-template-columns:1.2fr .8fr;gap:18px}.p{border:1px solid #243047;background:#0f172a;padding:14px}.lead{color:#cbd5e1;line-height:1.5}.fig{margin-top:14px}.fig img{max-width:100%;display:block;border:1px solid #243047}</style></head><body><div class='w'>",
        "<h1>Open Combo GEO Reanalysis</h1>",
        "<p class='lead'>This page re-analyzes the open combo cohorts GSE222011 and GSE255830 with sample-level T-cell state and clonotype summaries, while keeping the NeoPrecis vs NeoICI comparison contract-level rather than claiming a paired score rerun.</p>",
        "<div class='g'><div class='p'><h2>Combined sample table</h2>",
        combined.to_html(index=False, escape=False),
        "</div><div class='p'><h2>Contract table</h2>",
        contract.to_html(index=False, escape=False),
        "<div class='fig'><h2>Figures</h2>",
    ]
    for p in figs:
        html.append(f"<img src='assets/open_combo_geo/{p.name}' alt='{p.name}'/>")
    html += ["</div></div></div></body></html>"]
    html_text = "".join(html)
    html_path = HUB / "open_combo_geo_reanalysis.html"
    html_path.write_text(html_text, encoding="utf-8")
    (WWW / "open_combo_geo_reanalysis.html").write_text(html_text, encoding="utf-8")
    asset_dir = HUB / "assets" / "open_combo_geo"
    asset_www = WWW / "assets" / "open_combo_geo"
    ensure_dir(asset_dir)
    ensure_dir(asset_www)
    for p in figs:
        shutil.copy2(p, asset_dir / p.name)
        shutil.copy2(p, asset_www / p.name)

    meta = {
        "generated": "2026-05-11",
        "combined_rows": int(combined.shape[0]),
        "html": str(html_path),
        "summary": str(OUT / "COMBO_GEO_REANALYSIS_SUMMARY.md"),
    }
    (OUT / "combo_geo_reanalysis.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
