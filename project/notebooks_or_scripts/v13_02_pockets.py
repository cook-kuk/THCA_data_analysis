"""
v13 Task 2 — Pocket detection + druggability tier.

No p2rank/fpocket available, so we use:
  (a) HETATM ligand extraction — explicit binding sites for crystal structures.
  (b) Alpha-sphere-lite: find buried hydrophobic residue clusters as candidate pockets.
  (c) Druggability tier = combination of pocket PLB proxy + subcellular location
      + amenable modality class.

Output:
  pockets/{GENE}_pockets.tsv — top 3 pockets per target with center+residues
  pockets/druggability_tiers.tsv — tier + rationale per target
  pockets/pocket_detail.json — machine-readable details (for dashboard)
"""
from __future__ import annotations
import json, math
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd

PROJECT = Path("/opt/thyroid-dash/project")
RES = PROJECT / "results" / "v13_drug_discovery"
STRUCT = RES / "structures"
POCKETS = RES / "pockets"
POCKETS.mkdir(parents=True, exist_ok=True)

HYDROPHOBIC = set("ALA VAL LEU ILE MET PHE TRP TYR PRO CYS".split())
AROMATIC = set("PHE TRP TYR HIS".split())

# Subcellular location + biology-informed modality class
TARGET_CONTEXT = {
    # gene: (subcellular_location, target_class, preferred_modalities)
    "CYP1B1":  ("endoplasmic_reticulum_membrane",
                "heme_monooxygenase",
                ["SM"]),
    "LDLR":    ("plasma_membrane_extracellular",
                "internalizing_receptor",
                ["mAb", "ADC", "siRNA", "mRNA_replacement"]),
    "TACSTD2": ("plasma_membrane_extracellular",
                "adhesion_internalizing_receptor",
                ["ADC", "mAb", "bispecific", "CAR_T", "radioligand"]),
    "TMPRSS4": ("plasma_membrane_transmembrane_serine_protease",
                "serine_protease",
                ["SM", "mAb"]),
    "GABRB2":  ("plasma_membrane_ion_channel",
                "pentameric_ligand_gated_ion_channel",
                ["SM"]),
    "PLEKHA6": ("cytoplasm_junctional",
                "PH_domain_scaffold",
                ["PROTAC", "siRNA", "ASO"]),
    "PTPRE":   ("plasma_membrane_intracellular_catalytic",
                "receptor_tyrosine_phosphatase",
                ["SM", "PROTAC", "siRNA"]),
    "B3GNT3":  ("golgi_membrane",
                "glycosyltransferase",
                ["SM", "siRNA", "ASO"]),
}


def parse_pdb_first_model(path: Path):
    """Return list of (atom_name, res_name, chain, res_seq, icode, xyz, hetatm?)
       from the first model only."""
    atoms = []
    in_first = True
    for ln in path.read_text(errors="replace").splitlines():
        if ln.startswith("MODEL"):
            parts = ln.split()
            if len(parts) > 1 and parts[1] not in ("", "1"):
                in_first = False
            continue
        if ln.startswith("ENDMDL"):
            break
        if not in_first:
            continue
        if ln.startswith("ATOM") or ln.startswith("HETATM"):
            try:
                atoms.append({
                    "record": ln[:6].strip(),
                    "atom": ln[12:16].strip(),
                    "element": ln[76:78].strip() or ln[12:14].strip(),
                    "res": ln[17:20].strip(),
                    "chain": ln[21].strip() or "A",
                    "resseq": ln[22:26].strip(),
                    "icode": ln[26].strip(),
                    "x": float(ln[30:38]),
                    "y": float(ln[38:46]),
                    "z": float(ln[46:54]),
                })
            except Exception:
                continue
    return atoms


def extract_hetatm_ligands(atoms):
    """Group HETATM into ligands (exclude water, ions, common buffers)."""
    IGNORE = {"HOH", "WAT", "DOD", "NA", "CL", "K", "ZN", "MG", "CA", "FE",
              "SO4", "PO4", "EDO", "GOL", "PEG", "TRS", "ACT", "DMS"}
    ligs = defaultdict(list)
    for a in atoms:
        if a["record"] != "HETATM":
            continue
        if a["res"] in IGNORE:
            continue
        key = (a["res"], a["chain"], a["resseq"])
        ligs[key].append(a)
    return ligs


def residues_within(atoms, center_xyz, cutoff=5.0):
    """Return set of (chain, resseq) protein residues within cutoff of center."""
    cx, cy, cz = center_xyz
    hit = set()
    for a in atoms:
        if a["record"] != "ATOM":
            continue
        d = (a["x"]-cx)**2 + (a["y"]-cy)**2 + (a["z"]-cz)**2
        if d <= cutoff*cutoff:
            hit.add((a["chain"], a["resseq"], a["res"]))
    return hit


def detect_pockets_alpha(atoms, top_n=3):
    """
    Alpha-sphere-lite:
      1. For each heavy-atom pair within 4-9 Å, compute midpoint.
      2. Keep midpoints whose neighborhood (3.5Å) contains 0 atoms
         but whose wider neighborhood (8Å) contains >=15 atoms (buried).
      3. Cluster surviving midpoints (6Å radius, greedy).
      4. Score each cluster by:
           - size (α-sphere count)
           - hydrophobic fraction of lining residues
           - aromatic contribution
           - buriedness (mean neighbor density)
    """
    prot = [a for a in atoms if a["record"] == "ATOM" and a["element"] != "H"]
    if len(prot) < 200:
        return []
    coords = np.array([[a["x"], a["y"], a["z"]] for a in prot])

    # Subsample for speed on very large structures (e.g., 6X3X 2121 residues)
    N = len(coords)
    if N > 8000:
        # limit pairs: use KD-tree
        pass

    # Build simple KD-tree via scipy? Use numpy broadcasting in chunks.
    try:
        from scipy.spatial import cKDTree  # scipy is a standard ML dep
        kd = cKDTree(coords)
    except Exception:
        kd = None

    alpha_centers = []
    # Find close pairs (4-9 Å)
    if kd is not None:
        pairs = kd.query_pairs(r=9.0, output_type="ndarray")
        # filter to 4-9 Å
        if len(pairs) > 0:
            d = np.linalg.norm(coords[pairs[:, 0]] - coords[pairs[:, 1]], axis=1)
            keep = (d >= 4.0) & (d <= 9.0)
            pairs = pairs[keep]
            # limit
            if len(pairs) > 40000:
                rng = np.random.default_rng(0)
                idx = rng.choice(len(pairs), 40000, replace=False)
                pairs = pairs[idx]
            mids = 0.5 * (coords[pairs[:, 0]] + coords[pairs[:, 1]])
            # alpha sphere: no atom within 3.0Å, >=15 atoms within 8Å
            near_cts = kd.query_ball_point(mids, r=3.0, return_length=True)
            far_cts = kd.query_ball_point(mids, r=8.0, return_length=True)
            mask = (near_cts == 0) & (far_cts >= 15)
            alpha_centers = mids[mask]

    if len(alpha_centers) == 0:
        return []

    # Greedy cluster alpha centers (6 Å)
    clusters = []
    used = np.zeros(len(alpha_centers), dtype=bool)
    if kd is not None:
        ckd = cKDTree(alpha_centers)
        order = np.argsort(-alpha_centers[:, 0])  # arbitrary ordering
        for i in order:
            if used[i]:
                continue
            members = ckd.query_ball_point(alpha_centers[i], r=6.0)
            members = [m for m in members if not used[m]]
            if len(members) < 8:
                continue
            for m in members:
                used[m] = True
            pts = alpha_centers[members]
            center = pts.mean(axis=0)
            # volume estimate via bbox
            bbox = pts.max(axis=0) - pts.min(axis=0)
            vol = float(np.prod(bbox + 1e-3))
            clusters.append({"center": center.tolist(),
                             "n_alpha": int(len(members)),
                             "volume_bbox": round(vol, 2)})

    # Score each cluster
    for c in clusters:
        near = residues_within(atoms, c["center"], cutoff=5.0)
        c["lining_residues"] = len(near)
        if not near:
            c["hydrophobic_frac"] = 0.0
            c["aromatic_frac"] = 0.0
            c["residue_list"] = []
            c["plb"] = 0.0
            continue
        hydrop = sum(1 for _, _, r in near if r in HYDROPHOBIC) / len(near)
        arom = sum(1 for _, _, r in near if r in AROMATIC) / len(near)
        # proxy PLB (0..1)
        size_term = min(1.0, c["n_alpha"] / 40.0)
        plb = 0.40*size_term + 0.35*hydrop + 0.15*arom + 0.10*min(1.0, c["lining_residues"]/20.0)
        c["hydrophobic_frac"] = round(hydrop, 3)
        c["aromatic_frac"] = round(arom, 3)
        c["residue_list"] = [f"{ch}:{r}{seq}" for ch, seq, r in sorted(near, key=lambda x: (x[0], int(x[1]) if x[1].isdigit() else 0))]
        c["plb"] = round(plb, 3)

    clusters.sort(key=lambda c: -c["plb"])
    return clusters[:top_n]


def main():
    pocket_rows = []
    tier_rows = []
    detail = {}

    for gene in TARGET_CONTEXT:
        pdb = STRUCT / f"{gene}.pdb"
        if not pdb.exists():
            print(f"[{gene}] structure missing, skip")
            continue
        print(f"[{gene}] parsing {pdb.name}")
        atoms = parse_pdb_first_model(pdb)
        prot_atoms = [a for a in atoms if a["record"] == "ATOM"]
        # For crystallographic lattices with many chains, restrict to chain A to
        # keep pocket detection tractable and avoid inter-chain artifacts.
        chains = sorted({a["chain"] for a in prot_atoms})
        if len(prot_atoms) > 10000 and len(chains) > 2:
            atoms = [a for a in atoms if a["chain"] == chains[0]]
            print(f"  large multi-chain structure: restricting to chain {chains[0]}")
            prot_atoms = [a for a in atoms if a["record"] == "ATOM"]
        print(f"  atoms: {len(atoms)} total, {len(prot_atoms)} ATOM, chains={chains}")

        # Bound ligands (known binding sites)
        ligs = extract_hetatm_ligands(atoms)
        bound = []
        for key, la in ligs.items():
            if len(la) < 6:  # skip tiny
                continue
            cx = np.mean([a["x"] for a in la])
            cy = np.mean([a["y"] for a in la])
            cz = np.mean([a["z"] for a in la])
            near = residues_within(atoms, (cx, cy, cz), cutoff=5.0)
            if not near:
                continue
            hydrop = sum(1 for _, _, r in near if r in HYDROPHOBIC) / len(near)
            arom = sum(1 for _, _, r in near if r in AROMATIC) / len(near)
            # PLB proxy for known ligand site is high (0.7+)
            plb = round(0.70 + 0.15*hydrop + 0.10*arom, 3)
            bound.append({
                "source": f"HETATM:{key[0]}",
                "center": [round(cx, 2), round(cy, 2), round(cz, 2)],
                "volume_bbox": None,
                "lining_residues": len(near),
                "hydrophobic_frac": round(hydrop, 3),
                "aromatic_frac": round(arom, 3),
                "residue_list": [f"{ch}:{r}{seq}" for ch, seq, r in sorted(near, key=lambda x: (x[0], int(x[1]) if x[1].isdigit() else 0))],
                "plb": plb,
                "n_alpha": None,
                "ligand_hetcode": key[0],
            })
        bound.sort(key=lambda c: -c["plb"])
        print(f"  HETATM ligand sites: {len(bound)}")

        # De-novo pockets
        denovo = detect_pockets_alpha(atoms, top_n=5)
        print(f"  De-novo pockets detected: {len(denovo)}")

        # Merge, dedupe by center distance (<=6Å)
        merged = list(bound)
        for p in denovo:
            too_close = False
            for b in merged:
                d = math.dist(p["center"], b["center"])
                if d <= 6.0:
                    too_close = True
                    break
            if not too_close:
                merged.append(p)
        top3 = merged[:3]

        for rank, p in enumerate(top3, 1):
            pocket_rows.append({
                "gene": gene,
                "rank": rank,
                "plb": p["plb"],
                "volume_bbox": p.get("volume_bbox"),
                "n_alpha": p.get("n_alpha"),
                "lining_residues": p["lining_residues"],
                "hydrophobic_frac": p["hydrophobic_frac"],
                "aromatic_frac": p["aromatic_frac"],
                "center_x": p["center"][0],
                "center_y": p["center"][1],
                "center_z": p["center"][2],
                "ligand_hetcode": p.get("ligand_hetcode", ""),
                "residues_top10": ";".join(p["residue_list"][:10]),
            })

        detail[gene] = {"pockets": top3}

        # Druggability tier
        loc, cls, mods = TARGET_CONTEXT[gene]
        best_plb = top3[0]["plb"] if top3 else 0.0
        if best_plb >= 0.70:
            tier = "High"
        elif best_plb >= 0.50:
            tier = "Medium"
        else:
            tier = "Low"
        rationale = (
            f"Best pocket PLB={best_plb}. Target class={cls}. "
            f"Location={loc}. Preferred modalities={','.join(mods)}."
        )
        tier_rows.append({
            "gene": gene,
            "uniprot_class": cls,
            "subcellular": loc,
            "best_plb": best_plb,
            "pocket_count": len(top3),
            "has_crystal_ligand": any(p.get("ligand_hetcode") for p in top3),
            "druggability_tier": tier,
            "preferred_modalities": ",".join(mods),
            "rationale": rationale,
        })

    pd.DataFrame(pocket_rows).to_csv(POCKETS / "all_pockets.tsv", sep="\t", index=False)
    pd.DataFrame(tier_rows).to_csv(POCKETS / "druggability_tiers.tsv", sep="\t", index=False)
    # Per-gene TSVs too
    for gene in TARGET_CONTEXT:
        rows = [r for r in pocket_rows if r["gene"] == gene]
        if rows:
            pd.DataFrame(rows).to_csv(POCKETS / f"{gene}_pockets.tsv", sep="\t", index=False)
    (POCKETS / "pocket_detail.json").write_text(json.dumps(detail, indent=2, default=str))

    print("\n=== Druggability tiers ===")
    print(pd.DataFrame(tier_rows).to_string(index=False))


if __name__ == "__main__":
    main()
