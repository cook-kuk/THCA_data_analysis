"""Consolidate structure + compound data into assets/data/chem_targets.json.

The showcase page (`pages/platform_chem_showcase.html`) consumes a single
payload describing each of the 8 shortlisted targets: biomarker stats,
structure file reference, list of compounds (with SDF paths and pChEMBL),
rationale excerpt, and attribution.
"""

from __future__ import annotations

import csv
import json
import logging
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path("/opt/thyroid-dash/project")
HTML_ROOT = PROJECT_ROOT / "reports" / "html"
TARGETS_TSV = PROJECT_ROOT / "results" / "tables" / "drug_discovery_targets.tsv"
COMPOUNDS_TSV = PROJECT_ROOT / "results" / "tables" / "drug_discovery_compounds.tsv"
CATALOG_PATH = HTML_ROOT / "assets" / "data" / "structure_catalog.json"
LIGAND_MANIFEST = HTML_ROOT / "assets" / "data" / "ligands" / "_manifest.json"
REPORT_PATH = PROJECT_ROOT / "reports" / "biomarker_to_drug_report.md"
OUT_PATH = HTML_ROOT / "assets" / "data" / "chem_targets.json"
VERSION_PATH = HTML_ROOT / "_version.json"
LOG_PATH = PROJECT_ROOT / "logs" / "chem_showcase_build.log"


def _setup_logger() -> logging.Logger:
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("chem_showcase_build")
    logger.setLevel(logging.INFO)
    for h in list(logger.handlers):
        logger.removeHandler(h)
    fmt = logging.Formatter("%(asctime)s %(levelname)s chem_showcase_build %(message)s")
    fh = logging.FileHandler(LOG_PATH, mode="a"); fh.setFormatter(fmt)
    sh = logging.StreamHandler(sys.stdout); sh.setFormatter(fmt)
    logger.addHandler(fh); logger.addHandler(sh)
    return logger


def _parse_rationales() -> dict[str, str]:
    """Pull per-target rationale paragraphs out of biomarker_to_drug_report.md."""
    if not REPORT_PATH.exists():
        return {}
    text = REPORT_PATH.read_text()
    out: dict[str, str] = {}
    # pattern matches "### GENE — classification\n\nparagraph" blocks
    for m in re.finditer(r"### (\S+) — (\S+)\n(.*?)(?=\n###|\n## |\Z)", text, flags=re.DOTALL):
        gene = m.group(1).strip()
        body = m.group(3).strip()
        # first paragraph only
        first = body.split("\n\n", 1)[0].strip()
        out[gene] = first
    return out


def _fmt_float(v: str | float | None, digits: int = 3) -> float | None:
    if v is None or v == "":
        return None
    try:
        return round(float(v), digits)
    except (TypeError, ValueError):
        return None


def main() -> int:
    logger = _setup_logger()
    logger.info("=== chem_showcase_build START ===")

    catalog = json.loads(CATALOG_PATH.read_text()) if CATALOG_PATH.exists() else {}
    logger.info("loaded structure catalog: %d entries", len(catalog))

    ligand_manifest = {}
    if LIGAND_MANIFEST.exists():
        ligand_manifest = json.loads(LIGAND_MANIFEST.read_text())
    sdf_by_chembl = {
        c["chembl_id"]: c.get("sdf")
        for c in ligand_manifest.get("compounds", [])
        if c.get("ok")
    }
    logger.info("ligand SDFs available: %d", len(sdf_by_chembl))

    rationales = _parse_rationales()

    targets_rows: list[dict] = []
    with TARGETS_TSV.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            targets_rows.append(row)
    logger.info("targets rows: %d", len(targets_rows))

    compounds_by_target: dict[str, list[dict]] = {}
    with COMPOUNDS_TSV.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            t = row.get("target") or ""
            compounds_by_target.setdefault(t, []).append(row)

    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": {
            "targets": "results/tables/drug_discovery_targets.tsv",
            "compounds": "results/tables/drug_discovery_compounds.tsv",
            "structure_catalog": "assets/data/structure_catalog.json",
            "report": "reports/biomarker_to_drug_report.md",
        },
        "attribution": {
            "pdb":         "Structural data from RCSB PDB under public terms (rcsb.org).",
            "alphafold":   "AlphaFold DB predicted models by DeepMind / EMBL-EBI, licensed CC-BY-4.0.",
            "chembl":      "Bioactivity data from ChEMBL (EMBL-EBI), CC-BY-SA 3.0.",
            "visual_note": "Visual system inspired by modern drug-discovery platform conventions. Not affiliated with any commercial AI drug-discovery company.",
        },
        "targets": [],
    }

    n_with_struct = 0
    total_compound_count = 0

    for row in targets_rows:
        gene = row.get("gene") or ""
        cat_entry = catalog.get(gene, {}) or {}
        has_struct = bool(cat_entry.get("file"))
        if has_struct:
            n_with_struct += 1
        cmpds_raw = compounds_by_target.get(gene, [])
        cmpds_raw_sorted = sorted(
            cmpds_raw,
            key=lambda r: (-(float(r.get("pchembl") or 0.0))),
        )
        compounds = []
        for c in cmpds_raw_sorted:
            chid = (c.get("chembl_id") or "").strip()
            compounds.append({
                "chembl_id": chid,
                "name": c.get("name") or chid,
                "smiles": c.get("smiles") or "",
                "pchembl": _fmt_float(c.get("pchembl"), 2),
                "moa": (c.get("moa") or "")[:160],
                "relation": c.get("target_relation_confidence") or "",
                "sdf": sdf_by_chembl.get(chid),
            })
        total_compound_count += len(compounds)
        tgt = {
            "gene": gene,
            "chembl_target_id": row.get("chembl_target_id") or None,
            "chembl_pref_name": row.get("chembl_pref_name") or None,
            "classification": row.get("classification") or None,
            "novelty_score": _fmt_float(row.get("novelty_score"), 2),
            "biomarker_stats": {
                "log2FC":          _fmt_float(row.get("log2FC"), 3),
                "fdr":             _fmt_float(row.get("fdr"), 6),
                "replicated_27155":  int(row.get("replicated_27155") or 0),
                "replicated_126698": int(row.get("replicated_126698") or 0),
                "disease_focus_score": _fmt_float(row.get("disease_focus_score"), 1),
                "best_pchembl":    _fmt_float(row.get("best_pchembl"), 2),
                "n_compounds_total": int(row.get("n_compounds") or 0),
            },
            "rationale": rationales.get(gene, ""),
            "top_pmid": row.get("top_pmid") or None,
            "top_compound": row.get("top_compound") or None,
            "structure": {
                "source_kind": cat_entry.get("source_kind"),
                "identifier":  cat_entry.get("identifier"),
                "file":        cat_entry.get("file"),
                "ligand_hetcode":  cat_entry.get("ligand_hetcode"),
                "ligand_sdf":      cat_entry.get("ligand_sdf"),
                "binding_site_residues": cat_entry.get("binding_site_residues") or [],
                "source_url":  cat_entry.get("source_url"),
                "attribution": cat_entry.get("attribution"),
            },
            "compounds": compounds,
        }
        payload["targets"].append(tgt)

    payload["summary"] = {
        "n_targets": len(payload["targets"]),
        "n_targets_with_structure": n_with_struct,
        "n_compounds_total": total_compound_count,
        "n_ligand_sdfs": len(sdf_by_chembl),
    }

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(payload, indent=2))
    logger.info("wrote %s (%d targets, %d compounds, %d with structure)",
                OUT_PATH, len(payload["targets"]), total_compound_count, n_with_struct)

    # bump _version.json
    ver = {}
    if VERSION_PATH.exists():
        try:
            ver = json.loads(VERSION_PATH.read_text())
        except json.JSONDecodeError:
            ver = {}
    ver["chem_showcase_build_time"] = datetime.now(timezone.utc).isoformat()
    ver["chem_showcase_targets"] = len(payload["targets"])
    ver["chem_showcase_compounds"] = total_compound_count
    ver["chem_showcase_structures"] = n_with_struct
    VERSION_PATH.write_text(json.dumps(ver, indent=2))
    logger.info("bumped _version.json")

    logger.info("=== chem_showcase_build DONE ===")
    print(json.dumps(payload["summary"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
