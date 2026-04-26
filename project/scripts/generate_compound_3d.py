"""Generate 3D SDF conformers for every compound in drug_discovery_compounds.tsv.

For each row:
  1. parse SMILES via RDKit (fetch from ChEMBL API if the TSV cell is empty)
  2. add explicit hydrogens, embed a single ETKDGv3 conformer, UFF-minimize
  3. write an SDF to assets/data/ligands/{CHEMBL_ID}.sdf

Per-compound status lines are appended to logs/compound_3d.log.
"""

from __future__ import annotations

import csv
import json
import logging
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")  # silence the RDKit C++ noise

PROJECT_ROOT = Path("/opt/thyroid-dash/project")
TSV_PATH = PROJECT_ROOT / "results" / "tables" / "drug_discovery_compounds.tsv"
LIGAND_DIR = PROJECT_ROOT / "reports" / "html" / "assets" / "data" / "ligands"
LOG_PATH = PROJECT_ROOT / "logs" / "compound_3d.log"
BUILD_LOG = PROJECT_ROOT / "logs" / "chem_showcase_build.log"


def _setup_logger() -> logging.Logger:
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("compound_3d")
    logger.setLevel(logging.INFO)
    for h in list(logger.handlers):
        logger.removeHandler(h)
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    fh = logging.FileHandler(LOG_PATH, mode="a")
    fh.setFormatter(fmt)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(sh)
    # also copy INFO lines into the unified build log
    bh = logging.FileHandler(BUILD_LOG, mode="a")
    bh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s compound_3d %(message)s"))
    logger.addHandler(bh)
    return logger


def _fetch_chembl_smiles(chembl_id: str, logger: logging.Logger, timeout: int = 20) -> str | None:
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "thyrai-chem-showcase/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            payload = json.loads(resp.read())
    except Exception as exc:
        logger.warning("%s: ChEMBL fetch failed %s", chembl_id, exc)
        return None
    structs = payload.get("molecule_structures") or {}
    smi = structs.get("canonical_smiles")
    if smi:
        logger.info("%s: pulled SMILES from ChEMBL API", chembl_id)
    return smi


def _embed(smiles: str) -> Chem.Mol | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xC0FFEE
    if AllChem.EmbedMolecule(mol, params) != 0:
        # fall-back: try basic embedding if ETKDGv3 fails (rare rings, macrocycles)
        if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
            return None
    try:
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass
    return mol


def main() -> int:
    logger = _setup_logger()
    logger.info("=== generate_compound_3d START ===")
    LIGAND_DIR.mkdir(parents=True, exist_ok=True)

    rows = []
    with TSV_PATH.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)

    n_total = len(rows)
    n_ok = 0
    n_fail = 0
    per_compound: list[dict] = []
    for row in rows:
        chembl_id = (row.get("chembl_id") or "").strip()
        if not chembl_id:
            continue
        smiles = (row.get("smiles") or "").strip()
        if not smiles:
            smiles = _fetch_chembl_smiles(chembl_id, logger) or ""
            time.sleep(0.25)
        status = {"chembl_id": chembl_id, "target": row.get("target"), "name": row.get("name"),
                  "smiles": smiles, "sdf": None, "ok": False, "error": None}
        if not smiles:
            status["error"] = "no_smiles"
            logger.warning("%s: no SMILES available", chembl_id)
            per_compound.append(status); n_fail += 1
            continue
        mol = _embed(smiles)
        if mol is None:
            status["error"] = "embed_failed"
            logger.warning("%s: RDKit embed failed for %s", chembl_id, smiles)
            per_compound.append(status); n_fail += 1
            continue
        sdf_path = LIGAND_DIR / f"{chembl_id}.sdf"
        writer = Chem.SDWriter(str(sdf_path))
        mol.SetProp("_Name", row.get("name") or chembl_id)
        mol.SetProp("target", row.get("target") or "")
        mol.SetProp("pchembl", str(row.get("pchembl") or ""))
        writer.write(mol)
        writer.close()
        status["sdf"] = f"assets/data/ligands/{sdf_path.name}"
        status["ok"] = True
        per_compound.append(status)
        n_ok += 1
        logger.info("%s: SDF written (%s, pChEMBL=%s)", chembl_id, row.get("target"), row.get("pchembl"))

    # write a compact manifest so the showcase build step can consume it
    manifest_path = LIGAND_DIR / "_manifest.json"
    manifest_path.write_text(json.dumps({
        "total": n_total,
        "success": n_ok,
        "failed": n_fail,
        "compounds": per_compound,
    }, indent=2))
    logger.info("manifest written: %s", manifest_path)
    logger.info("=== generate_compound_3d DONE  ok=%d fail=%d total=%d ===", n_ok, n_fail, n_total)
    print(json.dumps({"ok": n_ok, "fail": n_fail, "total": n_total}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
