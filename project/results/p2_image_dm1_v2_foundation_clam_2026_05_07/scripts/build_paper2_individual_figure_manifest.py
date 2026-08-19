#!/usr/bin/env python3
"""Publish separated Paper 2 figure assets and a browser manifest."""
from __future__ import annotations

import json
import re
import shutil
from pathlib import Path


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
P2 = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET_DIR_NAME = "paper2_individual_figures"
LOCAL_ASSET = HUB / "assets" / ASSET_DIR_NAME
LIVE_ASSET = LIVE / "assets" / ASSET_DIR_NAME

GROUP_ORDER = {
    "GSE250521 slide overlays": 10,
    "Paper 2 main figures": 20,
    "Path2Space spatial reanalysis": 30,
    "GSE230424 external thyroid": 40,
    "Classifier audit controls": 50,
    "Supplementary model checks": 60,
    "GSE248205 AITD axis": 70,
    "Pathology exploratory axes": 80,
    "Phase 3 integration": 90,
}


def image_paths() -> list[Path]:
    paths: list[Path] = []
    for pattern in [
        "phase1_gse250521/*.png",
        "figures/*.png",
        "analysis_supp/figures/*.png",
        "analysis_supp/audit_ras_auc100/figures/*.png",
        "analysis_supp/audit_uni_loto/*.png",
        "analysis_supp/path2space_*/*.png",
        "analysis_supp/path2space_*/*/*.png",
        "analysis_supp/gse230424_pathology_thyroid_axis_2026_05_09/*.png",
        "analysis_supp/gse230424_pathology_thyroid_axis_2026_05_09/*/*.png",
        "analysis_supp/gse248205_pathology_aitd_axis_2026_05_09/*.png",
        "analysis_supp/pathology_*/*.png",
        "phase3_integration/*.png",
    ]:
        paths.extend(P2.glob(pattern))

    seen: set[Path] = set()
    unique: list[Path] = []
    for path in paths:
        if not path.is_file():
            continue
        if "paper2_cv2_visual_summary_2026_05_10" in path.parts:
            continue
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        unique.append(path)
    return sorted(unique, key=lambda p: (GROUP_ORDER.get(group_for(p), 999), p.relative_to(P2).as_posix()))


def group_for(path: Path) -> str:
    rel = path.relative_to(P2).as_posix()
    if rel.startswith("phase1_gse250521/"):
        return "GSE250521 slide overlays"
    if rel.startswith("figures/"):
        return "Paper 2 main figures"
    if rel.startswith("phase3_integration/"):
        return "Phase 3 integration"
    if "gse230424_pathology_thyroid_axis" in rel:
        return "GSE230424 external thyroid"
    if "gse248205_pathology_aitd_axis" in rel:
        return "GSE248205 AITD axis"
    if "/audit_" in f"/{rel}":
        return "Classifier audit controls"
    if "path2space" in rel:
        return "Path2Space spatial reanalysis"
    if rel.startswith("analysis_supp/figures/"):
        return "Supplementary model checks"
    if "pathology_" in rel:
        return "Pathology exploratory axes"
    return "Supplementary model checks"


def slug_for(path: Path) -> str:
    rel = path.relative_to(P2).as_posix()
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "__", rel)
    slug = slug.replace("/", "__")
    return slug


def title_for(path: Path) -> str:
    stem = path.stem
    if stem.startswith("uni_dm1_overlay_"):
        return stem.replace("uni_dm1_overlay_", "GSE250521 ").replace("_", " ")
    label = re.sub(r"^(fig|figure)[_-]?", "", stem, flags=re.IGNORECASE)
    label = label.replace("_", " ").replace("-", " ")
    label = re.sub(r"\s+", " ", label).strip()
    if not label:
        label = stem
    return label[:1].upper() + label[1:]


def caption_for(path: Path) -> str:
    rel = path.relative_to(P2).as_posix()
    group = group_for(path)
    if group == "GSE250521 slide overlays":
        return "Single-slide observed/predicted UNI-DM1 spatial overlay; kept separate from the 16-slide mosaic."
    if group == "GSE230424 external thyroid":
        return "External thyroid Visium/H&E analysis figure from the GSE230424 support layer."
    if group == "Path2Space spatial reanalysis":
        return "Path2Space-inspired spatial RNA projection or control figure."
    if group == "Classifier audit controls":
        return "Image-DM1 classifier audit/control figure."
    return f"Separated original figure from {rel}."


def write_manifest(entries: list[dict[str, str]], asset_dir: Path) -> None:
    payload = json.dumps(entries, indent=2, ensure_ascii=False)
    (asset_dir / "manifest.json").write_text(payload + "\n", encoding="utf-8")
    (asset_dir / "manifest.js").write_text(
        "window.paper2IndividualFigures = " + payload + ";\n",
        encoding="utf-8",
    )


def main() -> None:
    LOCAL_ASSET.mkdir(parents=True, exist_ok=True)
    LIVE_ASSET.mkdir(parents=True, exist_ok=True)

    entries: list[dict[str, str]] = []
    for index, path in enumerate(image_paths(), start=1):
        dest_name = slug_for(path)
        local_dest = LOCAL_ASSET / dest_name
        live_dest = LIVE_ASSET / dest_name
        shutil.copy2(path, local_dest)
        shutil.copy2(path, live_dest)
        entries.append(
            {
                "id": f"p2fig{index:03d}",
                "group": group_for(path),
                "title": title_for(path),
                "caption": caption_for(path),
                "src": f"assets/{ASSET_DIR_NAME}/{dest_name}",
                "sourcePath": path.relative_to(ROOT).as_posix(),
            }
        )

    write_manifest(entries, LOCAL_ASSET)
    write_manifest(entries, LIVE_ASSET)
    print(f"published {len(entries)} separated figures")
    print(LOCAL_ASSET)
    print(LIVE_ASSET)


if __name__ == "__main__":
    main()
