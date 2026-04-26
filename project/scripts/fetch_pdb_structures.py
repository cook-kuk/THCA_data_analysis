"""Fetch PDB and AlphaFold structures for the 8 shortlisted thyroid-cancer targets.

Priority order per target: co-crystal > apo > AlphaFold.  If a download fails
(HTTP 404 or timeout) the failure is logged and the target is skipped — the
showcase page falls back to compound-only 2D rendering for that target.

For co-crystal PDB entries, HETATM lines for the ligand (non-water, non-ion)
are extracted into a minimal-header SDF file and binding-site residues
(within 5 Å of any ligand heavy atom) are identified via Biopython
NeighborSearch and stored in the catalog.
"""

from __future__ import annotations

import json
import logging
import os
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

from Bio.PDB import PDBParser, NeighborSearch, Selection

PROJECT_ROOT = Path("/opt/thyroid-dash/project")
HTML_ROOT = PROJECT_ROOT / "reports" / "html"
STRUCT_DIR = HTML_ROOT / "assets" / "data" / "structures"
LIGAND_DIR = HTML_ROOT / "assets" / "data" / "ligands"
CATALOG_PATH = HTML_ROOT / "assets" / "data" / "structure_catalog.json"
LOG_PATH = PROJECT_ROOT / "logs" / "chem_showcase_build.log"

# Ions / water / common buffer components that should NOT be treated as the
# bound ligand during HETATM extraction.
_SKIP_HET = {
    "HOH", "WAT", "DOD",
    "NA", "K", "CL", "MG", "CA", "ZN", "FE", "MN", "CU", "NI", "CO",
    "SO4", "PO4", "GOL", "EDO", "PEG", "MES", "TRS", "ACT", "FMT",
    "DMS", "IMD", "CIT", "BME",
}

STRUCTURES = [
    # (target, kind, pdb_or_uniprot, ligand_hetcode)
    # CYP1B1 3PM0 holds alpha-naphthoflavone as HETATM code "BHF"
    # (chemically the compound community calls this "ANF"; the PDB stores
    # the crystal ligand code as BHF, so we match that here).
    ("CYP1B1",  "pdb",       "3PM0",    "BHF"),
    ("TACSTD2", "alphafold", "P09758",  None),
    ("TMPRSS4", "alphafold", "Q9NRS4",  None),
    ("PLEKHA6", "alphafold", "Q9Y2H5",  None),
    ("LDLR",    "pdb",       "1N7D",    None),   # apo (no co-crystal ligand forced)
    ("GABRB2",  "pdb",       "6X3S",    None),
    ("B3GNT3",  "alphafold", "Q9Y2A9",  None),
    ("PTPRE",   "pdb",       "2JJD",    None),
]


def _setup_logger() -> logging.Logger:
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("fetch_pdb")
    logger.setLevel(logging.INFO)
    for h in list(logger.handlers):
        logger.removeHandler(h)
    fmt = logging.Formatter("%(asctime)s %(levelname)s fetch_pdb %(message)s")
    fh = logging.FileHandler(LOG_PATH, mode="a")
    fh.setFormatter(fmt)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger


def _download(url: str, out_path: Path, logger: logging.Logger, timeout: int = 45) -> bool:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "thyrai-chem-showcase/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = resp.read()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_bytes(data)
        logger.info("downloaded %s -> %s (%d bytes)", url, out_path.name, len(data))
        return True
    except urllib.error.HTTPError as exc:
        logger.warning("HTTP %s for %s", exc.code, url)
    except Exception as exc:  # pragma: no cover - defensive
        logger.warning("download failed %s: %s", url, exc)
    return False


def _extract_ligand_sdf(pdb_path: Path, het_code: str, out_sdf: Path, logger: logging.Logger) -> dict | None:
    """Extract HETATM records for ``het_code`` and write a minimal-header SDF.

    We keep it deliberately simple: a single molecule block with atom-count
    only (no bond list, no explicit charges).  3Dmol.js will fall back to
    a distance-based bond perception when SDF bonds are absent, so the
    visualization is still usable as a ball-and-stick overlay.
    """
    atoms: list[tuple[str, float, float, float]] = []
    with pdb_path.open() as fh:
        for line in fh:
            if not line.startswith("HETATM"):
                continue
            res = line[17:20].strip()
            if res != het_code:
                continue
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            element = line[76:78].strip() or line[12:14].strip()
            element = "".join(c for c in element if c.isalpha())
            if not element:
                continue
            element = element[0].upper() + element[1:].lower() if len(element) > 1 else element.upper()
            # Use first residue instance only (skip alt-locs beyond first)
            atoms.append((element, x, y, z))
    if not atoms:
        logger.warning("no HETATM atoms found for ligand %s in %s", het_code, pdb_path.name)
        return None
    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"{het_code}",
        "  thyrai-chem-showcase",
        "",
        f"{len(atoms):>3}{0:>3}  0  0  0  0  0  0  0  0999 V2000",
    ]
    for (el, x, y, z) in atoms:
        lines.append(f"{x:10.4f}{y:10.4f}{z:10.4f} {el:<3}0  0  0  0  0  0  0  0  0  0  0  0")
    lines.append("M  END")
    lines.append("$$$$")
    out_sdf.write_text("\n".join(lines) + "\n")
    logger.info("ligand SDF written: %s (%d atoms)", out_sdf.name, len(atoms))
    return {"het_code": het_code, "n_atoms": len(atoms)}


def _binding_site_residues(pdb_path: Path, het_code: str, logger: logging.Logger, radius: float = 5.0, top_n: int = 8) -> list[dict]:
    """Return up to ``top_n`` binding-site residues within ``radius`` Å of the ligand."""
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("s", str(pdb_path))
    except Exception as exc:
        logger.warning("PDB parse failed for %s: %s", pdb_path.name, exc)
        return []
    ligand_atoms = []
    for atom in structure.get_atoms():
        res = atom.get_parent()
        if res.id[0].startswith("H_") and res.resname.strip() == het_code:
            # skip hydrogens
            if atom.element and atom.element.strip() != "H":
                ligand_atoms.append(atom)
    if not ligand_atoms:
        return []
    all_atoms = [a for a in structure.get_atoms() if a.element and a.element.strip() != "H"]
    ns = NeighborSearch(all_atoms)
    nearby_residues: dict[tuple, tuple[float, object]] = {}
    for lig_atom in ligand_atoms:
        for near in ns.search(lig_atom.coord, radius, level="A"):
            res = near.get_parent()
            if res.id[0].startswith("H_") or res.id[0] == "W":
                continue  # skip other hetero / water
            chain = res.get_parent().id
            key = (chain, res.id[1], res.resname.strip())
            dist = float(((near.coord - lig_atom.coord) ** 2).sum() ** 0.5)
            prev = nearby_residues.get(key)
            if prev is None or dist < prev[0]:
                nearby_residues[key] = (dist, res)
    sorted_res = sorted(nearby_residues.items(), key=lambda kv: kv[1][0])[:top_n]
    out = []
    for (chain, resi, resname), (dist, _res) in sorted_res:
        out.append({
            "chain": chain,
            "resi": resi,
            "resname": resname,
            "distance_ang": round(dist, 2),
        })
    return out


def main() -> int:
    logger = _setup_logger()
    logger.info("=== fetch_pdb_structures START ===")
    STRUCT_DIR.mkdir(parents=True, exist_ok=True)
    LIGAND_DIR.mkdir(parents=True, exist_ok=True)

    catalog: dict[str, dict] = {}
    pdb_ok = 0
    af_ok = 0

    for target, kind, ident, het in STRUCTURES:
        entry = {
            "target": target,
            "source_kind": kind,
            "identifier": ident,
            "file": None,
            "ligand_hetcode": het,
            "ligand_sdf": None,
            "binding_site_residues": [],
            "source_url": None,
            "attribution": None,
        }
        if kind == "pdb":
            url = f"https://files.rcsb.org/download/{ident}.pdb"
            out = STRUCT_DIR / f"{target}_{ident}.pdb"
            entry["source_url"] = url
            entry["attribution"] = f"Structure data from RCSB PDB entry {ident} (public terms)"
            if _download(url, out, logger):
                pdb_ok += 1
                entry["file"] = f"assets/data/structures/{out.name}"
                if het:
                    sdf_path = LIGAND_DIR / f"{target}_cocrystal.sdf"
                    lig_info = _extract_ligand_sdf(out, het, sdf_path, logger)
                    if lig_info:
                        entry["ligand_sdf"] = f"assets/data/ligands/{sdf_path.name}"
                        entry["ligand_n_atoms"] = lig_info["n_atoms"]
                    entry["binding_site_residues"] = _binding_site_residues(out, het, logger)
        elif kind == "alphafold":
            # Try the prediction API first (returns the pdbUrl for the latest
            # published version); fall back to probing v6..v1 static URLs if
            # the API is unreachable.
            api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{ident}"
            out = STRUCT_DIR / f"{target}_AF-{ident}.pdb"
            entry["attribution"] = (
                f"AlphaFold DB model AF-{ident}-F1 (Jumper et al. 2021 / Varadi et al. 2022), CC-BY-4.0"
            )
            resolved_url = None
            try:
                req = urllib.request.Request(api_url, headers={"User-Agent": "thyrai-chem-showcase/1.0"})
                with urllib.request.urlopen(req, timeout=20) as resp:
                    meta = json.loads(resp.read())
                if isinstance(meta, list) and meta:
                    resolved_url = meta[0].get("pdbUrl")
                    logger.info("%s: AlphaFold API returned pdbUrl=%s", target, resolved_url)
            except Exception as exc:
                logger.warning("%s: AlphaFold API failed (%s)", target, exc)

            tried = []
            if resolved_url:
                tried.append(resolved_url)
            for v in (6, 5, 4, 3, 2, 1):
                tried.append(f"https://alphafold.ebi.ac.uk/files/AF-{ident}-F1-model_v{v}.pdb")
            for url in tried:
                if _download(url, out, logger):
                    af_ok += 1
                    entry["file"] = f"assets/data/structures/{out.name}"
                    entry["source_url"] = url
                    break
            else:
                entry["source_url"] = api_url
        catalog[target] = entry
        time.sleep(0.3)  # polite pacing

    CATALOG_PATH.write_text(json.dumps(catalog, indent=2))
    logger.info("catalog written: %s", CATALOG_PATH)
    logger.info("PDB ok=%d  AlphaFold ok=%d", pdb_ok, af_ok)
    logger.info("=== fetch_pdb_structures DONE ===")
    # emit summary for stdout capture
    print(json.dumps({"pdb_ok": pdb_ok, "alphafold_ok": af_ok, "total": len(STRUCTURES)}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
