#!/usr/bin/env python3
"""GSE230424 (4 PTC+HT Visium samples, Yan PMC) parser.
Format per sample: {GSM}_{Pn}_{barcodes,features,matrix,tissue_positions_list,HE}
NO scalefactors_json or tissue_hires_image.png — store HE.jpg path only.
"""
from __future__ import annotations
import argparse, gzip, logging, re, shutil, sys
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("parse_gse230424")

# user spec: "4 samples: Hashimoto thyroiditis + PTC combined with Hashimoto thyroiditis"
# GEO sample_titles say only "P1, ST-seq" .. "P4, ST-seq" with no per-sample disease state.
# Default → all 4 as PTC_HT, mark ambiguity in registry notes.
SAMPLE_RE = re.compile(r"^(GSM\d+)_(P\d+)_(.+)$")
DEFAULT_CONDITION = "PTC_HT"
DEFAULT_NOTE = "GEO sample_title=P[1-4],ST-seq; per-patient HT vs PTC+HT split not in metadata; user spec=PTC+HT all"


def find_samples(extract: Path) -> dict[str, dict]:
    samples: dict[str, dict] = {}
    for f in sorted(extract.iterdir()):
        m = SAMPLE_RE.match(f.name)
        if not m:
            log.warning("skip unrecognized: %s", f.name); continue
        gsm, label, suffix = m.groups()
        sid = f"{gsm}_{label}"
        s = samples.setdefault(sid, {"gsm": gsm, "label": label})
        if suffix.startswith("barcodes"): s["barcodes"] = f
        elif suffix.startswith("features"): s["features"] = f
        elif suffix.startswith("matrix"): s["matrix"] = f
        elif suffix.startswith("tissue_positions"): s["positions"] = f
        elif suffix.startswith("HE"): s["image"] = f
        else: log.warning("unknown suffix: %s", f.name)
    return samples


def build(sid: str, info: dict, processed: Path) -> ad.AnnData:
    for k in ("barcodes","features","matrix","positions"):
        if k not in info: raise FileNotFoundError(f"{sid}: missing {k}")
    log.info("[%s] reading mtx", sid)
    with gzip.open(info["matrix"], "rb") as fh:
        X = mmread(fh).tocsr()
    barcodes = pd.read_csv(info["barcodes"], header=None, sep="\t", compression="gzip")[0].tolist()
    feats = pd.read_csv(info["features"], header=None, sep="\t", compression="gzip")
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
    pos = pd.read_csv(info["positions"], header=None, compression="gzip")
    cols6 = ["barcode","in_tissue","array_row","array_col","pxl_row_in_fullres","pxl_col_in_fullres"]
    if pos.shape[1] == 6: pos.columns = cols6
    elif pos.shape[1] == 7: pos.columns = cols6 + ["extra"]
    else: raise ValueError(f"{sid}: positions {pos.shape[1]} cols")
    pos = pos.set_index("barcode")
    obs = obs.join(pos[["in_tissue","array_row","array_col","pxl_row_in_fullres","pxl_col_in_fullres"]],
                   how="left")
    img_dst = processed / sid / "HE.jpg"
    img_dst.parent.mkdir(parents=True, exist_ok=True)
    if "image" in info:
        with gzip.open(info["image"], "rb") as src, open(img_dst, "wb") as dst:
            shutil.copyfileobj(src, dst)
    a = ad.AnnData(X=csr_matrix(X), obs=obs, var=var)
    a.obs["sample_id"] = sid
    a.obs["dataset"] = "GSE230424"
    a.obs["source_accession"] = info["gsm"]
    a.obs["condition_raw"] = info["label"]
    a.obs["condition_inferred"] = DEFAULT_CONDITION
    a.obs["disease_axis"] = DEFAULT_CONDITION
    a.obsm["spatial"] = obs[["pxl_col_in_fullres","pxl_row_in_fullres"]].to_numpy()
    a.uns["spatial"] = {sid: {"images": {"HE_path": str(img_dst) if "image" in info else None},
                              "scalefactors": {}}}  # GSE230424 lacks scalefactors_json
    log.info("[%s] %d spots × %d genes (in_tissue=%d)",
             sid, a.n_obs, a.n_vars, int(obs["in_tissue"].fillna(0).sum()))
    return a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--extract", default="project_external_st/data/raw/GSE230424/extracted")
    ap.add_argument("--processed", default="project_external_st/data/processed/GSE230424")
    ap.add_argument("--meta", default="project_external_st/data/processed/GSE230424/sample_metadata.tsv")
    args = ap.parse_args()
    extract = Path(args.extract); processed = Path(args.processed)
    processed.mkdir(parents=True, exist_ok=True)
    samples = find_samples(extract)
    rows = []
    for sid, info in samples.items():
        try:
            a = build(sid, info, processed)
        except Exception as e:
            log.error("[%s] %s", sid, e); rows.append({"sample_id": sid, "status": "fail", "error": str(e)}); continue
        out = processed / sid / f"{sid}.raw.h5ad"
        a.write_h5ad(out, compression="gzip")
        rows.append({"dataset": "GSE230424", "sample_id": sid, "source_accession": info["gsm"],
                     "condition_raw": info["label"], "condition_inferred": DEFAULT_CONDITION,
                     "disease_axis": DEFAULT_CONDITION,
                     "has_HE_image": "image" in info, "has_spatial_coords": True,
                     "n_spots_raw": a.n_obs, "n_genes": a.n_vars,
                     "n_in_tissue": int(a.obs["in_tissue"].fillna(0).sum()),
                     "h5ad": str(out), "image": str(processed / sid / "HE.jpg") if "image" in info else "",
                     "notes": DEFAULT_NOTE, "status": "ok"})
    meta = pd.DataFrame(rows)
    Path(args.meta).parent.mkdir(parents=True, exist_ok=True)
    meta.to_csv(args.meta, sep="\t", index=False)
    print(meta.to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
