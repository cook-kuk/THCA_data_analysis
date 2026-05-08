"""
Phase 2a — TCGA-THCA diagnostic WSI selective download from GDC public API.

NOTE: TCGA-THCA diagnostic SVS is OPEN ACCESS (not controlled).
Strategy: balance DM1/DM2/not_DM. Default 60-100 slides.

Input:
  master TCGA labels: project/results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv
  (sample_id with -01A suffix; dm = DM1 / DM2 / not_DM)

Output:
  WSI files into out_dir/
  slide_manifest.tsv  (sample_id, dm, file_id, file_name, file_size, status)

Manifest first; download via gdc-client or direct HTTP.
"""
from __future__ import annotations
import argparse, json, time, os
from pathlib import Path
import pandas as pd
import requests

GDC_FILES_API = "https://api.gdc.cancer.gov/files"
GDC_DOWNLOAD_API = "https://api.gdc.cancer.gov/data/"

def query_thca_diagnostic(case_ids: list[str], data_format: str = "SVS") -> pd.DataFrame:
    """Query GDC for diagnostic-WSI files for given TCGA-THCA case IDs."""
    all_records = []
    for chunk_start in range(0, len(case_ids), 50):
        chunk = case_ids[chunk_start:chunk_start+50]
        filters = {
            "op": "and",
            "content": [
                {"op": "in", "content": {"field": "cases.case_id", "value": chunk}},
                {"op": "in", "content": {"field": "data_format", "value": [data_format]}},
                {"op": "in", "content": {"field": "experimental_strategy", "value": ["Diagnostic Slide"]}},
            ],
        }
        params = {"filters": json.dumps(filters),
                  "fields": "file_id,file_name,file_size,cases.submitter_id,cases.case_id",
                  "format": "JSON", "size": 200}
        r = requests.get(GDC_FILES_API, params=params, timeout=60)
        r.raise_for_status()
        for hit in r.json().get("data", {}).get("hits", []):
            all_records.append({
                "file_id": hit.get("file_id"),
                "file_name": hit.get("file_name"),
                "file_size": hit.get("file_size"),
                "case_id": hit.get("cases", [{}])[0].get("case_id"),
                "submitter_id": hit.get("cases", [{}])[0].get("submitter_id"),
            })
    return pd.DataFrame(all_records)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--master", type=Path,
                   default=Path("project/results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv"))
    p.add_argument("--out_dir", type=Path,
                   default=Path("project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_clam"))
    p.add_argument("--n_per_group", type=int, default=30,
                   help="DM1=30 + DM2=30 + not_DM=30 (default; 90 total)")
    p.add_argument("--manifest_only", action="store_true",
                   help="Only build manifest; do not download")
    p.add_argument("--max_size_gb_per_slide", type=float, default=2.5,
                   help="Skip slides larger than this (typical TCGA SVS = 0.5-2 GB)")
    args = p.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    wsi_dir = args.out_dir / "wsi"
    wsi_dir.mkdir(exist_ok=True)

    # --- load master + balance ---
    master = pd.read_csv(args.master, sep="\t")
    master["short"] = master["sample_id"].str[:12]
    print(f"[load] master n={len(master)} ; dm distribution:")
    print(master["dm"].value_counts())

    # GDC needs case UUID, not submitter_id (TCGA-XX-XXXX).
    # Map TCGA-XX-XXXX → case_id via GDC cases API
    print("[map] TCGA submitter_id → case_id ...")
    submitter_to_case = {}
    submitters = sorted(master["short"].unique())
    for chunk_start in range(0, len(submitters), 100):
        chunk = submitters[chunk_start:chunk_start+100]
        f = {"op": "in", "content": {"field": "submitter_id", "value": chunk}}
        r = requests.get("https://api.gdc.cancer.gov/cases",
                         params={"filters": json.dumps(f),
                                 "fields": "case_id,submitter_id",
                                 "format": "JSON", "size": 200}, timeout=60)
        r.raise_for_status()
        for hit in r.json().get("data", {}).get("hits", []):
            submitter_to_case[hit["submitter_id"]] = hit["case_id"]
    print(f"  mapped {len(submitter_to_case)} cases")

    # --- balance sampling ---
    rng = pd.Series(range(len(master)))
    selected = []
    for grp in ["DM1", "DM2", "not_DM"]:
        subset = master[master["dm"] == grp].sample(n=min(args.n_per_group, (master["dm"]==grp).sum()),
                                                      random_state=42)
        selected.append(subset)
    sel = pd.concat(selected)
    sel["case_id"] = sel["short"].map(submitter_to_case)
    sel = sel.dropna(subset=["case_id"])
    print(f"[balance] selected {len(sel)} cases (DM1={sum(sel['dm']=='DM1')}, DM2={sum(sel['dm']=='DM2')}, not_DM={sum(sel['dm']=='not_DM')})")

    # --- query GDC for diagnostic SVS files ---
    print("[query] GDC for diagnostic WSI files ...")
    files_df = query_thca_diagnostic(sel["case_id"].tolist(), "SVS")
    print(f"  {len(files_df)} SVS files found")

    # join dm label
    case2dm = dict(zip(sel["case_id"], sel["dm"]))
    files_df["dm"] = files_df["case_id"].map(case2dm)

    # filter by size
    files_df["file_size_gb"] = files_df["file_size"] / 1e9
    files_df = files_df[files_df["file_size_gb"] <= args.max_size_gb_per_slide].copy()
    print(f"  after size filter: {len(files_df)} files")

    # one slide per case (pick smallest if multiple)
    files_df = files_df.sort_values("file_size").groupby("case_id").head(1)

    # save manifest
    manifest_path = args.out_dir / "slide_manifest.tsv"
    files_df.to_csv(manifest_path, sep="\t", index=False)
    print(f"[manifest] saved {len(files_df)} slides to {manifest_path}")

    if args.manifest_only:
        print("[stop] --manifest_only requested; not downloading")
        return

    # --- download ---
    total_gb = files_df["file_size_gb"].sum()
    print(f"[download] total {total_gb:.1f} GB ; into {wsi_dir}")
    for i, row in files_df.reset_index().iterrows():
        out_path = wsi_dir / f"{row['file_id']}_{row['file_name']}"
        if out_path.exists() and out_path.stat().st_size == row["file_size"]:
            print(f"  [{i+1}/{len(files_df)}] skip exist: {out_path.name}")
            continue
        url = GDC_DOWNLOAD_API + row["file_id"]
        try:
            with requests.get(url, stream=True, timeout=300) as r:
                r.raise_for_status()
                with out_path.open("wb") as f:
                    for chunk in r.iter_content(chunk_size=1<<20):
                        f.write(chunk)
            print(f"  [{i+1}/{len(files_df)}] {out_path.name}  ({out_path.stat().st_size/1e9:.2f} GB)")
        except Exception as e:
            print(f"  [{i+1}/{len(files_df)}] FAIL: {e}")
    print("[done]")

if __name__ == "__main__":
    main()
