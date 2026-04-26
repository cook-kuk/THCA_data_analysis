"""
v13 Task 1 — Fetch 3D structures for 8 THCA targets.

Strategy:
  1. Query RCSB PDB for each UniProt ID → download best-resolution X-ray/cryo-EM if available.
  2. Fallback to AlphaFold DB precomputed AF2 model.
  3. Record provenance in structure_sources.tsv.

We explicitly do NOT call AlphaFold3 API (gated access, per prompt).
"""
from __future__ import annotations
import json, time, os, sys, io
from pathlib import Path
import requests
import pandas as pd

PROJECT = Path("/opt/thyroid-dash/project")
RES = PROJECT / "results" / "v13_drug_discovery"
STRUCT = RES / "structures"
STRUCT.mkdir(parents=True, exist_ok=True)

TARGETS = {
    # gene : uniprot accession (canonical human)
    "CYP1B1":  "Q16678",
    "LDLR":    "P01130",
    "TACSTD2": "P09758",   # TROP2
    "TMPRSS4": "Q9NRS4",
    "GABRB2":  "P47870",
    "PLEKHA6": "Q9Y2H5",
    "PTPRE":   "P23469",
    "B3GNT3":  "Q9Y2A9",
}

# Curated starting-point PDBs (from prompt + literature). Override empty for API discovery.
CURATED_PDB = {
    "CYP1B1":  "3PM0",   # CYP1B1–alpha-naphthoflavone 2.7 Å
    "LDLR":    "1N7D",   # EGF domain
    "TACSTD2": "",       # no crystal — AF2 only
    "TMPRSS4": "",
    "GABRB2":  "6X3X",   # cryo-EM GABA-A α1β2γ2
    "PLEKHA6": "",
    "PTPRE":   "2JJD",   # PTPRE PTP domain
    "B3GNT3":  "",
}

UA = {"User-Agent": "thca-v13-pipeline (kukshomr@gmail.com)"}


def fetch(url: str, timeout: int = 30, binary: bool = False):
    try:
        r = requests.get(url, headers=UA, timeout=timeout)
        if r.status_code == 200:
            return r.content if binary else r.text
    except Exception as e:
        print(f"  fetch {url} → {type(e).__name__}: {e}")
    return None


def discover_pdb_for_uniprot(up: str) -> list[dict]:
    """Query RCSB polymer-entity endpoint by UniProt."""
    url = (
        "https://search.rcsb.org/rcsbsearch/v2/query?json="
        + requests.utils.quote(json.dumps({
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                    "operator": "exact_match",
                    "value": up,
                }
            },
            "return_type": "polymer_entity",
            "request_options": {"paginate": {"rows": 10}},
        }))
    )
    txt = fetch(url)
    if not txt:
        return []
    try:
        js = json.loads(txt)
        hits = [h["identifier"] for h in js.get("result_set", [])]
        return hits[:5]
    except Exception:
        return []


def download_pdb(pdb_id: str) -> bytes | None:
    pdb_id = pdb_id.split("_")[0].lower()
    return fetch(f"https://files.rcsb.org/download/{pdb_id}.pdb", binary=True)


def download_alphafold(up: str) -> tuple[bytes | None, dict]:
    meta_url = f"https://alphafold.ebi.ac.uk/api/prediction/{up}"
    meta_txt = fetch(meta_url)
    if not meta_txt:
        return None, {}
    try:
        meta = json.loads(meta_txt)
        if not meta:
            return None, {}
        entry = meta[0]
        pdb_url = entry.get("pdbUrl")
        if not pdb_url:
            return None, {}
        pdb_bytes = fetch(pdb_url, binary=True)
        return pdb_bytes, entry
    except Exception as e:
        print(f"  AF parse fail {up}: {e}")
        return None, {}


def main():
    rows = []
    for gene, up in TARGETS.items():
        print(f"[{gene}] UniProt={up}")
        source = None
        outfile = STRUCT / f"{gene}.pdb"
        meta = {"gene": gene, "uniprot": up}

        def is_long_enough(b: bytes) -> bool:
            if not b:
                return False
            # count unique residues in the FIRST model only
            try:
                in_first = True
                seen = set()
                for ln in b.decode("latin-1").splitlines():
                    if ln.startswith("MODEL") and ln.split()[1] not in ("", "1"):
                        in_first = False
                        continue
                    if ln.startswith("ENDMDL"):
                        break
                    if in_first and ln.startswith("ATOM") and ln[12:16].strip() == "CA":
                        seen.add((ln[21], ln[22:27].strip()))
                return len(seen) >= 100
            except Exception:
                return False

        # 1) curated PDB
        pdb_id = CURATED_PDB.get(gene, "")
        pdb_bytes = None
        if pdb_id:
            cand = download_pdb(pdb_id)
            if is_long_enough(cand):
                pdb_bytes = cand
                source = f"PDB:{pdb_id}"
                meta["pdb_id"] = pdb_id

        # 2) PDB discovery via search API
        if not pdb_bytes:
            hits = discover_pdb_for_uniprot(up)
            for h in hits:
                cand = download_pdb(h)
                if is_long_enough(cand):
                    pdb_bytes = cand
                    source = f"PDB:{h}"
                    meta["pdb_id"] = h
                    break

        # 3) AlphaFold DB fallback
        used_af = False
        if not pdb_bytes:
            af_bytes, af_meta = download_alphafold(up)
            if af_bytes:
                pdb_bytes = af_bytes
                source = f"AF2:{af_meta.get('entryId','?')}"
                meta["af_entry"] = af_meta.get("entryId", "")
                meta["af_version"] = af_meta.get("latestVersion", "")
                meta["af_model_created"] = af_meta.get("modelCreatedDate", "")
                meta["af_plddt_mean"] = ""  # we'll compute below
                used_af = True

        if not pdb_bytes:
            print(f"  WARN no structure obtained for {gene}")
            meta["source"] = "NONE"
            rows.append(meta)
            continue

        outfile.write_bytes(pdb_bytes)
        size_kb = len(pdb_bytes) // 1024
        meta["source"] = source
        meta["file_kb"] = size_kb

        # quick stats — residues + mean B-factor (= plDDT in AF models)
        try:
            from Bio.PDB import PDBParser
            parser = PDBParser(QUIET=True)
            st = parser.get_structure(gene, str(outfile))
            res_ct = 0
            bfs = []
            for model in st:
                for chain in model:
                    for residue in chain:
                        if residue.id[0] != " ":
                            continue
                        res_ct += 1
                        for atom in residue:
                            bfs.append(atom.get_bfactor())
                break  # first model only
            meta["residues"] = res_ct
            if bfs:
                import statistics
                meta["bfactor_mean"] = round(statistics.mean(bfs), 2)
            if used_af and bfs:
                meta["af_plddt_mean"] = round(statistics.mean(bfs), 2)
        except Exception as e:
            print(f"  parse fail {gene}: {e}")

        print(f"  → {source}  {size_kb} KB  residues={meta.get('residues','?')}  B/plDDT={meta.get('bfactor_mean','?')}")
        rows.append(meta)
        time.sleep(0.5)

    df = pd.DataFrame(rows)
    df.to_csv(STRUCT / "structure_sources.tsv", sep="\t", index=False)
    print(f"\nWrote {STRUCT/'structure_sources.tsv'}  ({len(df)} rows)")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
