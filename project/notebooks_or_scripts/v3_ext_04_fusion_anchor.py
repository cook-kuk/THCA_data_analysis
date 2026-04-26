#!/usr/bin/env python
"""v3_ext_04_fusion_anchor.py

Build 6-class v3_anchor_6class from merged fusion callset + existing mutation calls:
  BRAF_V600E / RAS_mutant / RET_fusion / NTRK_fusion / PAX8_PPARG / other

Output:
  metadata/v3_fusion_anchor_tcga.tsv
  metadata/v3_fusion_anchor_prjeb11591.tsv
  per-class counts; flag classes with n<15 as small_n
"""
import json
import logging
from pathlib import Path
import pandas as pd

ROOT = Path("/opt/thyroid-dash/project")
META = ROOT / "metadata"
RAW = ROOT / "data_raw" / "v3_ext"
LOGDIR = ROOT / "logs"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[logging.FileHandler(LOGDIR / "v3_ext_04_fusion_anchor.log", mode="w"),
              logging.StreamHandler()],
)
log = logging.getLogger("v3_ext_04")


FUSION_REGEX = {
    "RET_fusion": ["RET-", "-RET", "CCDC6-RET", "NCOA4-RET", "KIF5B-RET"],
    "NTRK_fusion": ["NTRK1-", "-NTRK1", "NTRK3-", "-NTRK3", "ETV6-NTRK3", "TPM3-NTRK1"],
    "PAX8_PPARG": ["PAX8-PPARG", "PPARG-PAX8"],
}


def parse_cbioportal_sv(path: Path):
    if not path.exists():
        return pd.DataFrame()
    try:
        data = json.loads(path.read_text())
        if not data:
            return pd.DataFrame()
        rows = []
        for r in data:
            s1 = r.get("site1HugoSymbol", "") or ""
            s2 = r.get("site2HugoSymbol", "") or ""
            sid = r.get("sampleId", "")
            rows.append({"sample_id": sid, "fusion_pair": f"{s1}-{s2}"})
        return pd.DataFrame(rows)
    except Exception as e:
        log.warning(f"cbio sv parse: {e}")
        return pd.DataFrame()


def parse_tumorfusions(path: Path):
    if not path.exists():
        return pd.DataFrame()
    try:
        df = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
        if "Cancer" in df.columns:
            df = df[df["Cancer"].astype(str).str.contains("THCA", case=False, na=False)]
        gene_cols = [c for c in df.columns if c.lower() in ("gene_a", "gene_b", "fusion_name", "fusionname")]
        sid_cols = [c for c in df.columns if c.lower() in ("sample", "sample_id", "tumor_sample_barcode")]
        if not gene_cols or not sid_cols:
            return pd.DataFrame()
        out = df[[sid_cols[0]] + gene_cols].copy()
        out.columns = ["sample_id"] + gene_cols
        if "Fusion_name" in df.columns or "fusion_name" in df.columns:
            col = "Fusion_name" if "Fusion_name" in df.columns else "fusion_name"
            out["fusion_pair"] = df[col]
        else:
            out["fusion_pair"] = out[gene_cols[0]].astype(str) + "-" + out[gene_cols[-1]].astype(str)
        return out[["sample_id", "fusion_pair"]]
    except Exception as e:
        log.warning(f"tumorfusions parse: {e}")
        return pd.DataFrame()


def classify_fusion(pair: str) -> str:
    up = (pair or "").upper()
    for cls, pats in FUSION_REGEX.items():
        for p in pats:
            if p in up:
                return cls
    return ""


def build_tcga():
    mut_f = ROOT / "results/tables/tcga_thca_mutation_groups.tsv"
    mut = pd.read_csv(mut_f, sep="\t") if mut_f.exists() else pd.DataFrame()
    sv1 = parse_cbioportal_sv(RAW / "tcga_fusions" / "cbio_structural_variants.json")
    sv2 = parse_tumorfusions(RAW / "tcga_fusions" / "pancanfus.txt.gz")
    sv = pd.concat([sv1, sv2], ignore_index=True) if not (sv1.empty and sv2.empty) else pd.DataFrame()
    log.info(f"TCGA SV rows: {len(sv)}")
    if not sv.empty:
        sv["fusion_class"] = sv["fusion_pair"].apply(classify_fusion)
        # keep most-informative per sample
        sv = sv[sv["fusion_class"] != ""].drop_duplicates("sample_id")

    # Build per-sample class
    sm = pd.read_csv(META / "sample_master_v3.tsv", sep="\t") if (META / "sample_master_v3.tsv").exists() \
        else pd.read_csv(META / "sample_master.tsv", sep="\t")
    tcga = sm[sm["dataset"].astype(str).str.contains("TCGA", case=False, na=False)].copy()

    braf = set(mut.loc[mut.get("has_braf_v600e") == True, "sample_id"]) if not mut.empty else set()
    ras = set(mut.loc[mut.get("has_ras_mut") == True, "sample_id"]) if not mut.empty else set()
    fus_map = dict(zip(sv["sample_id"], sv["fusion_class"])) if not sv.empty else {}

    def _cls(sid):
        if sid in braf:
            return "BRAF_V600E"
        if sid in ras:
            return "RAS_mutant"
        if sid in fus_map:
            return fus_map[sid]
        return "other"

    tcga["v3_anchor_6class"] = tcga["sample_id"].apply(_cls)
    counts = tcga["v3_anchor_6class"].value_counts().to_dict()
    log.info(f"TCGA 6-class counts: {counts}")
    small_n = {k: (v < 15) for k, v in counts.items()}
    tcga["small_n"] = tcga["v3_anchor_6class"].map(small_n).fillna(False)
    out = META / "v3_fusion_anchor_tcga.tsv"
    tcga[["sample_id", "dataset", "v3_anchor_6class", "small_n"]].to_csv(out, sep="\t", index=False)
    log.info(f"wrote {out}")
    return counts


def build_prjeb():
    out = META / "v3_fusion_anchor_prjeb11591.tsv"
    ena_f = RAW / "prjeb11591" / "ena_filereport.tsv"
    if not ena_f.exists():
        pd.DataFrame({"sample_id": [], "v3_anchor_6class": [], "small_n": []}).to_csv(out, sep="\t", index=False)
        log.info(f"PRJEB11591 metadata missing -> empty {out}")
        return {}
    em = pd.read_csv(ena_f, sep="\t")
    # sample_title often contains histology e.g. "PDTC", "ATC", "PTC BRAF"
    rows = []
    for _, r in em.iterrows():
        title = str(r.get("sample_title", "")).upper()
        cls = "other"
        if "BRAF" in title:
            cls = "BRAF_V600E"
        elif "RAS" in title or "NRAS" in title or "HRAS" in title or "KRAS" in title:
            cls = "RAS_mutant"
        elif "RET" in title:
            cls = "RET_fusion"
        rows.append({"sample_id": r.get("run_accession"), "v3_anchor_6class": cls, "small_n": False})
    df = pd.DataFrame(rows)
    if not df.empty:
        counts = df["v3_anchor_6class"].value_counts().to_dict()
        df["small_n"] = df["v3_anchor_6class"].map(lambda c: counts.get(c, 0) < 15)
    df.to_csv(out, sep="\t", index=False)
    log.info(f"wrote {out}  counts={df['v3_anchor_6class'].value_counts().to_dict() if not df.empty else {}}")
    return df["v3_anchor_6class"].value_counts().to_dict() if not df.empty else {}


def main():
    log.info("v3_ext_04_fusion_anchor start")
    c1 = build_tcga()
    c2 = build_prjeb()
    (ROOT / "results" / "tables" / "v3_fusion_anchor_counts.json").write_text(
        json.dumps({"tcga": c1, "prjeb11591": c2}, indent=2))
    log.info("v3_ext_04_fusion_anchor done")


if __name__ == "__main__":
    main()
