#!/usr/bin/env python3
"""Parse GSE250521 Visium GEO supplementary into per-sample AnnData h5ad."""
from __future__ import annotations
import argparse, gzip, json, logging, re, shutil, sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("parse")

STAGE_FROM_LABEL = {"N": "PT", "PTC": "PTC", "LPTC": "LPTC", "ATC": "ATC"}
SAMPLE_RE = re.compile(r"^(GSM\d+)_([A-Z]+)-(\d+)_visium_(.+)$")


def find_samples(extract_dir: Path) -> dict[str, dict]:
    """Group GEO files by sample id."""
    samples: dict[str, dict] = {}
    for f in sorted(extract_dir.iterdir()):
        m = SAMPLE_RE.match(f.name)
        if not m:
            log.warning("skip unrecognized: %s", f.name)
            continue
        gsm, label, idx, suffix = m.groups()
        sid = f"{gsm}_{label}-{idx}"
        s = samples.setdefault(sid, {"gsm": gsm, "label": label, "idx": int(idx),
                                     "stage_inferred": STAGE_FROM_LABEL.get(label, "UNKNOWN"),
                                     "stage_metadata_raw": label})
        if "barcodes" in suffix: s["barcodes"] = f
        elif "features" in suffix: s["features"] = f
        elif "matrix" in suffix: s["matrix"] = f
        elif "scalefactors" in suffix: s["scalefactors"] = f
        elif "tissue_hires_image" in suffix: s["image"] = f
        elif "tissue_positions_list" in suffix: s["positions"] = f
        else: log.warning("unknown suffix: %s", f.name)
    return samples


def read_gz(p: Path):
    with gzip.open(p, "rt") as fh:
        return fh.read()


def build_anndata(sid: str, info: dict, processed_dir: Path) -> ad.AnnData:
    required = ["barcodes", "features", "matrix", "positions", "scalefactors", "image"]
    missing = [k for k in required if k not in info]
    if missing:
        raise FileNotFoundError(f"{sid}: missing {missing}")

    log.info("[%s] reading mtx", sid)
    with gzip.open(info["matrix"], "rb") as fh:
        X = mmread(fh).tocsr()  # genes x spots
    barcodes = pd.read_csv(info["barcodes"], header=None, sep="\t", compression="gzip")[0].tolist()
    features = pd.read_csv(info["features"], header=None, sep="\t", compression="gzip")
    if features.shape[1] >= 2:
        features.columns = ["ensembl_id", "gene_symbol"] + [f"c{i}" for i in range(features.shape[1] - 2)]
    else:
        features.columns = ["gene_symbol"]
        features["ensembl_id"] = features["gene_symbol"]

    if X.shape[0] == len(features) and X.shape[1] == len(barcodes):
        X = X.T.tocsr()
    elif X.shape[0] == len(barcodes) and X.shape[1] == len(features):
        pass
    else:
        raise ValueError(f"{sid}: mtx shape {X.shape} vs features {len(features)} barcodes {len(barcodes)}")

    obs = pd.DataFrame(index=pd.Index(barcodes, name="spot_id"))
    var = features.set_index("gene_symbol")
    var.index = var.index.astype(str)
    var.index.name = None
    if not var.index.is_unique:
        idx = pd.Series(var.index)
        dup_mask = idx.duplicated(keep=False)
        if dup_mask.any():
            counts = idx.groupby(idx).cumcount()
            idx = idx.where(~dup_mask, idx + "-" + counts.astype(str))
            var.index = idx.values
            log.warning("[%s] %d duplicate gene symbols disambiguated with -N suffix",
                        sid, int(dup_mask.sum()))

    pos = pd.read_csv(info["positions"], header=None, compression="gzip")
    pos_cols = ["barcode", "in_tissue", "array_row", "array_col",
                "pxl_row_in_fullres", "pxl_col_in_fullres"]
    if pos.shape[1] == 6:
        pos.columns = pos_cols
    elif pos.shape[1] == 7:
        pos.columns = pos_cols + ["extra"]
    else:
        raise ValueError(f"{sid}: positions has {pos.shape[1]} cols")
    pos = pos.set_index("barcode")
    obs = obs.join(pos[["in_tissue", "array_row", "array_col",
                        "pxl_row_in_fullres", "pxl_col_in_fullres"]], how="left")

    sf = json.loads(read_gz(info["scalefactors"]))

    img_dst = processed_dir / sid / "tissue_hires_image.png"
    img_dst.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(info["image"], "rb") as src, open(img_dst, "wb") as dst:
        shutil.copyfileobj(src, dst)

    adata = ad.AnnData(X=csr_matrix(X), obs=obs, var=var)
    adata.obs["sample_id"] = sid
    adata.obs["stage"] = info["stage_inferred"]
    adata.obs["stage_metadata_raw"] = info["stage_metadata_raw"]
    adata.obs["gsm"] = info["gsm"]
    adata.obsm["spatial"] = obs[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
    adata.uns["spatial"] = {sid: {
        "scalefactors": sf,
        "images": {"hires_path": str(img_dst)},
    }}
    log.info("[%s] anndata: %d spots × %d genes (in_tissue=%d)",
             sid, adata.n_obs, adata.n_vars, int(obs["in_tissue"].fillna(0).sum()))
    return adata


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--extract-dir", default="project/data/raw/GSE250521/extracted")
    ap.add_argument("--processed-dir", default="project/data/processed/GSE250521")
    ap.add_argument("--meta-out", default="project/data/processed/GSE250521/sample_metadata.tsv")
    args = ap.parse_args()

    extract_dir = Path(args.extract_dir)
    processed_dir = Path(args.processed_dir)
    processed_dir.mkdir(parents=True, exist_ok=True)

    samples = find_samples(extract_dir)
    log.info("found %d samples", len(samples))
    rows = []
    for sid, info in samples.items():
        try:
            adata = build_anndata(sid, info, processed_dir)
        except Exception as e:
            log.error("[%s] build failed: %s", sid, e)
            rows.append({"sample_id": sid, "status": "fail", "error": str(e)})
            continue
        out = processed_dir / sid / f"{sid}.raw.h5ad"
        out.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(out, compression="gzip")
        rows.append({"sample_id": sid, "gsm": info["gsm"],
                     "stage_metadata_raw": info["stage_metadata_raw"],
                     "stage_inferred": info["stage_inferred"],
                     "n_spots": adata.n_obs, "n_genes": adata.n_vars,
                     "n_in_tissue": int(adata.obs["in_tissue"].fillna(0).sum()),
                     "h5ad": str(out), "image": str(processed_dir / sid / "tissue_hires_image.png"),
                     "status": "ok"})
    meta = pd.DataFrame(rows)
    Path(args.meta_out).parent.mkdir(parents=True, exist_ok=True)
    meta.to_csv(args.meta_out, sep="\t", index=False)
    log.info("metadata → %s", args.meta_out)
    print(meta.to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
