#!/usr/bin/env python3
"""GSE248205 (8 Visium samples: 2 control, 3 HT, 3 GD; Martinez-Hernandez et al.) parser.
Input: project_external_st/data/raw/GSE248205/extracted/Processed files/{label}/{label}_*
Output: project_external_st/data/processed/GSE248205/{sample_id}/*.raw.h5ad
"""
from __future__ import annotations
import argparse, gzip, json, logging, sys, shutil
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("parse_gse248205")

LABEL_TO_CONDITION = {"C1":"CONTROL","C2":"CONTROL",
                      "HT1":"HT","HT2":"HT","HT3":"HT",
                      "GD1":"GD","GD2":"GD","GD3":"GD"}
LABEL_TO_GSM = {"C1":"GSM7908359","C2":"GSM7908360",
                "HT1":"GSM7908361","HT2":"GSM7908362","HT3":"GSM7908363",
                "GD1":"GSM7908364","GD2":"GSM7908365","GD3":"GSM7908366"}


def maybe_gz_open(p: Path):
    return gzip.open(p, "rb") if str(p).endswith(".gz") else open(p, "rb")


def build(label: str, sample_dir: Path, processed: Path) -> ad.AnnData:
    sid = f"{LABEL_TO_GSM[label]}_{label}"
    files = {f.name: f for f in sample_dir.iterdir()}
    needed = {"barcodes": f"{label}_barcodes.tsv.gz",
              "features": f"{label}_features.tsv.gz",
              "matrix":   f"{label}_matrix.mtx.gz",
              "positions":f"{label}_tissue_positions_list.csv",
              "scalefactors": f"{label}_scalefactors_json.json",
              "image":    f"{label}_tissue_hires_image.png"}
    paths = {k: files.get(v) for k, v in needed.items()}
    miss = [k for k, v in paths.items() if v is None]
    if miss: raise FileNotFoundError(f"{sid}: missing {miss}")
    log.info("[%s] reading mtx", sid)
    with gzip.open(paths["matrix"], "rb") as fh:
        X = mmread(fh).tocsr()
    barcodes = pd.read_csv(paths["barcodes"], header=None, sep="\t", compression="gzip")[0].tolist()
    feats = pd.read_csv(paths["features"], header=None, sep="\t", compression="gzip")
    if feats.shape[1] >= 2:
        feats.columns = ["ensembl_id","gene_symbol"] + [f"c{i}" for i in range(feats.shape[1]-2)]
    else:
        feats.columns = ["gene_symbol"]; feats["ensembl_id"] = feats["gene_symbol"]
    if X.shape[0] == len(feats) and X.shape[1] == len(barcodes):
        X = X.T.tocsr()
    elif X.shape[0] == len(barcodes) and X.shape[1] == len(feats):
        pass
    else:
        raise ValueError(f"{sid}: mtx shape {X.shape} vs feats {len(feats)} bc {len(barcodes)}")
    obs = pd.DataFrame(index=pd.Index(barcodes, name="spot_id"))
    var = feats.set_index("gene_symbol")
    var.index = var.index.astype(str); var.index.name = None
    if not var.index.is_unique:
        idx = pd.Series(var.index)
        dup = idx.duplicated(keep=False)
        if dup.any():
            counts = idx.groupby(idx).cumcount()
            idx = idx.where(~dup, idx + "-" + counts.astype(str))
            var.index = idx.values
            log.warning("[%s] %d duplicate gene symbols disambiguated", sid, int(dup.sum()))
    pos = pd.read_csv(paths["positions"], header=None)
    cols6 = ["barcode","in_tissue","array_row","array_col","pxl_row_in_fullres","pxl_col_in_fullres"]
    if pos.shape[1] == 6: pos.columns = cols6
    elif pos.shape[1] == 7: pos.columns = cols6 + ["extra"]
    else: raise ValueError(f"{sid}: positions {pos.shape[1]} cols")
    pos = pos.set_index("barcode")
    obs = obs.join(pos[["in_tissue","array_row","array_col","pxl_row_in_fullres","pxl_col_in_fullres"]],
                   how="left")
    sf = json.loads(paths["scalefactors"].read_text())
    img_dst = processed / sid / "tissue_hires_image.png"
    img_dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(paths["image"], img_dst)
    a = ad.AnnData(X=csr_matrix(X), obs=obs, var=var)
    a.obs["sample_id"] = sid
    a.obs["dataset"] = "GSE248205"
    a.obs["source_accession"] = LABEL_TO_GSM[label]
    a.obs["condition_raw"] = label
    a.obs["condition_inferred"] = LABEL_TO_CONDITION[label]
    a.obs["disease_axis"] = LABEL_TO_CONDITION[label]
    a.obsm["spatial"] = obs[["pxl_col_in_fullres","pxl_row_in_fullres"]].to_numpy()
    a.uns["spatial"] = {sid: {"scalefactors": sf, "images": {"hires_path": str(img_dst)}}}
    log.info("[%s] %d spots × %d genes (in_tissue=%d)",
             sid, a.n_obs, a.n_vars, int(obs["in_tissue"].fillna(0).sum()))
    return a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--extract", default="project_external_st/data/raw/GSE248205/extracted/Processed files")
    ap.add_argument("--processed", default="project_external_st/data/processed/GSE248205")
    ap.add_argument("--meta", default="project_external_st/data/processed/GSE248205/sample_metadata.tsv")
    args = ap.parse_args()
    extract = Path(args.extract); processed = Path(args.processed)
    processed.mkdir(parents=True, exist_ok=True)
    rows = []
    for label_dir in sorted(extract.iterdir()):
        if not label_dir.is_dir(): continue
        label = label_dir.name
        if label not in LABEL_TO_CONDITION:
            log.warning("skip non-sample dir: %s", label); continue
        sid = f"{LABEL_TO_GSM[label]}_{label}"
        try:
            a = build(label, label_dir, processed)
        except Exception as e:
            log.error("[%s] %s", sid, e); rows.append({"sample_id": sid, "status": "fail", "error": str(e)}); continue
        out = processed / sid / f"{sid}.raw.h5ad"
        a.write_h5ad(out, compression="gzip")
        rows.append({"dataset": "GSE248205", "sample_id": sid,
                     "source_accession": LABEL_TO_GSM[label],
                     "condition_raw": label, "condition_inferred": LABEL_TO_CONDITION[label],
                     "disease_axis": LABEL_TO_CONDITION[label],
                     "has_HE_image": True, "has_spatial_coords": True,
                     "n_spots_raw": a.n_obs, "n_genes": a.n_vars,
                     "n_in_tissue": int(a.obs["in_tissue"].fillna(0).sum()),
                     "h5ad": str(out), "image": str(processed / sid / "tissue_hires_image.png"),
                     "notes": "Martinez-Hernandez 2024 PMID 39003267; full Visium output (incl scalefactors+hires png)",
                     "status": "ok"})
    meta = pd.DataFrame(rows)
    Path(args.meta).parent.mkdir(parents=True, exist_ok=True)
    meta.to_csv(args.meta, sep="\t", index=False)
    print(meta.to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
