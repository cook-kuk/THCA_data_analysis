#!/usr/bin/env python3
"""Shared helpers for CROSS-Neo MD audit scripts."""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
FIG = OUT / "figures"
QC = OUT / "qc"
CONTACTS = OUT / "contacts"
COUNTER = OUT / "counterfactual"
REPL = OUT / "replicates"
INTEGRATED = OUT / "integrated"

EXPECTED_PEPTIDES = {
    "HMTEVVRHC": "HMTEVVRHC",
    "GADGVGKSAL": "GADGVGKSAL",
}

AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "MSE": "M",
    "CYX": "C",
    "CYM": "C",
}

WATER_ION_RESNAMES = {
    "HOH",
    "WAT",
    "SOL",
    "NA",
    "CL",
    "K",
    "MG",
    "CA",
    "ZN",
    "SOD",
    "CLA",
}


@dataclass(frozen=True)
class RunRecord:
    run_id: str
    candidate: str
    peptide: str
    hla: str
    source: str
    run_dir: Path
    topology: Path | None
    trajectory: Path | None
    state: Path | None
    metadata: Path | None
    log: Path | None
    checkpoint: Path | None
    final_pdb: Path | None
    condition: str
    target_ns: float | None
    timestep_fs: float | None
    is_primary: bool
    is_replicate: bool
    is_smoke: bool


def ensure_dirs() -> None:
    for d in [OUT, FIG, QC, CONTACTS, COUNTER, REPL, INTEGRATED]:
        d.mkdir(parents=True, exist_ok=True)


def read_json(path: Path | None) -> dict:
    if path and path.exists():
        try:
            return json.loads(path.read_text())
        except Exception:
            return {}
    return {}


def infer_peptide(text: str) -> str:
    upper = text.upper()
    if "6VRN" in upper or "6VRM" in upper or "6VQO" in upper or "7RM4" in upper:
        return "HMTEVVRHC"
    if "6UON" in upper or "6ULK" in upper or "6ULI" in upper or "6ULR" in upper or "6ULN" in upper:
        return "GADGVGKSAL"
    for pep in EXPECTED_PEPTIDES:
        if pep in text:
            return pep
    return ""


def infer_hla(text: str) -> str:
    upper = text.upper()
    if "6VRN" in upper or "6VRM" in upper or "6VQO" in upper or "7RM4" in upper:
        return "HLA-A*02:01"
    if "6UON" in upper or "6ULK" in upper or "6ULI" in upper or "6ULR" in upper or "6ULN" in upper:
        return "HLA-C*08:02"
    if "HLA-A0201" in text or "A0201" in text:
        return "HLA-A*02:01"
    if "HLA-C0802" in text or "C0802" in text:
        return "HLA-C*08:02"
    m = re.search(r"HLA-([ABC])\*?(\d{2}):?(\d{2})", text)
    if m:
        return f"HLA-{m.group(1)}*{m.group(2)}:{m.group(3)}"
    return ""


def safe_id(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text)).strip("_")


def first_existing(paths: Iterable[Path]) -> Path | None:
    for p in paths:
        if p.exists():
            return p
    return None


def discover_runs() -> list[RunRecord]:
    ensure_dirs()
    roots = [
        OUT / "remote_sync",
        REPO / "project/results/cross_neo_v2_sota_sprint_2026_05_09/tcr_extension/md_escalation/p0_structures/openmm_pilot_10ns_package",
        REPO / "project/results/cross_neo_v2_sota_sprint_2026_05_09/tcr_extension/md_escalation/p0_structures/openmm_p0_replicate_screen_package",
    ]
    run_dirs: set[Path] = set()
    markers = {"state.tsv", "trajectory.dcd", "run_metadata.json", "prepared_start.pdb"}
    for root in roots:
        if not root.exists():
            continue
        for p in root.rglob("*"):
            if p.is_file() and p.name in markers:
                run_dirs.add(p.parent)

    records: list[RunRecord] = []
    for d in sorted(run_dirs):
        name = d.name
        text = str(d)
        peptide = infer_peptide(text)
        hla = infer_hla(text)
        candidate = f"{peptide}/{hla}" if peptide or hla else name
        metadata = first_existing([d / "run_metadata.json"])
        meta = read_json(metadata)
        topology = first_existing([d / "prepared_start.pdb", d / "minimized_start.pdb", d / "final.pdb"])
        trajectory = first_existing([d / "trajectory.dcd", d / "trajectory.xtc"])
        state = first_existing([d / "state.tsv"])
        log = first_existing([d / f"{name}.nohup.log", d.parent / f"{name}.nohup.log"])
        checkpoint = first_existing([d / "final.chk"])
        final_pdb = first_existing([d / "final.pdb"])
        source = "local"
        if "pod1_thca_neo_bayesian_aux" in text:
            source = "runpod_pod1_thca_neo_bayesian_aux"
        elif "pod2_thca_img_dm1_a6000" in text:
            source = "runpod_pod2_thca_img_dm1_a6000"
        elif "/smoke/" in text:
            source = "local_cpu_smoke"
        condition = "explicit_cuda" if source.startswith("runpod") else "local_smoke_or_prepared"
        if "smoke" in name or "/smoke/" in text:
            condition = "vacuum_smoke"
        if "prod_10ns" in name:
            condition = "explicit_cuda_10ns"
        elif "screens" in text:
            condition = "explicit_cuda_replicate_screen"
        target_ns = meta.get("ns_requested")
        if target_ns is None:
            if "10ns" in name:
                target_ns = 10.0
            elif "1p0ns" in name:
                target_ns = 1.0
            elif "0p5ns" in name:
                target_ns = 0.5
        timestep_fs = meta.get("timestep_fs")
        records.append(
            RunRecord(
                run_id=safe_id(name),
                candidate=candidate,
                peptide=peptide,
                hla=hla,
                source=source,
                run_dir=d,
                topology=topology,
                trajectory=trajectory,
                state=state,
                metadata=metadata,
                log=log,
                checkpoint=checkpoint,
                final_pdb=final_pdb,
                condition=condition,
                target_ns=float(target_ns) if target_ns is not None else None,
                timestep_fs=float(timestep_fs) if timestep_fs is not None else None,
                is_primary="prod_10ns" in name,
                is_replicate="screens" in text or "replicate" in text,
                is_smoke="smoke" in text,
            )
        )
    return records


def parse_state(path: Path | None) -> pd.DataFrame:
    if path is None or not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        df = pd.read_csv(path, sep="\t", comment=None)
    except Exception:
        return pd.DataFrame()
    df.columns = [str(c).replace('"', "").replace("#", "").strip() for c in df.columns]
    rename = {
        "Step": "step",
        "Time (ps)": "time_ps",
        "Potential Energy (kJ/mole)": "potential_energy_kj_mol",
        "Kinetic Energy (kJ/mole)": "kinetic_energy_kj_mol",
        "Total Energy (kJ/mole)": "total_energy_kj_mol",
        "Temperature (K)": "temperature_k",
        "Box Volume (nm^3)": "volume_nm3",
        "Density (g/mL)": "density_g_ml",
        "Speed (ns/day)": "speed_ns_per_day",
    }
    df = df.rename(columns=rename)
    for c in df.columns:
        converted = pd.to_numeric(df[c], errors="coerce")
        if converted.notna().sum() == df[c].notna().sum():
            df[c] = converted
    return df


def latest_state_summary(run: RunRecord) -> dict:
    df = parse_state(run.state)
    meta = read_json(run.metadata)
    out = {
        "run_id": run.run_id,
        "candidate": run.candidate,
        "peptide": run.peptide,
        "hla": run.hla,
        "source": run.source,
        "condition": run.condition,
        "run_dir": str(run.run_dir),
        "target_ns": run.target_ns,
        "timestep_fs": run.timestep_fs or meta.get("timestep_fs"),
        "state_rows": len(df),
        "trajectory_path": str(run.trajectory) if run.trajectory else "",
        "trajectory_bytes": run.trajectory.stat().st_size if run.trajectory and run.trajectory.exists() else 0,
        "topology_path": str(run.topology) if run.topology else "",
        "metadata_path": str(run.metadata) if run.metadata else "",
        "checkpoint_path": str(run.checkpoint) if run.checkpoint else "",
        "final_pdb_path": str(run.final_pdb) if run.final_pdb else "",
        "complete_by_metadata": bool(meta.get("final_pdb") or meta.get("elapsed_seconds")),
    }
    if not df.empty:
        last = df.iloc[-1].to_dict()
        out.update({k: last.get(k) for k in ["step", "time_ps", "potential_energy_kj_mol", "kinetic_energy_kj_mol", "total_energy_kj_mol", "temperature_k", "volume_nm3", "density_g_ml", "speed_ns_per_day"] if k in df.columns})
        if run.target_ns:
            out["completion_fraction"] = min(float(out.get("time_ps", 0)) / (run.target_ns * 1000.0), 1.0)
        else:
            out["completion_fraction"] = math.nan
    else:
        out["completion_fraction"] = 0.0
    return out


def pdb_chain_sequences(path: Path) -> pd.DataFrame:
    rows = []
    if not path or not path.exists():
        return pd.DataFrame()
    chain_residues: dict[str, list[tuple[int, str, int]]] = {}
    seen: set[tuple[str, int, str]] = set()
    for line in path.open(errors="ignore"):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        resname = line[17:20].strip()
        chain = line[21].strip() or "_"
        if resname in WATER_ION_RESNAMES:
            continue
        if resname not in AA3_TO_1:
            continue
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        icode = line[26].strip()
        key = (chain, resseq, icode)
        if key in seen:
            continue
        seen.add(key)
        chain_residues.setdefault(chain, []).append((resseq, AA3_TO_1[resname], len(chain_residues.get(chain, [])) + 1))
    for chain, residues in sorted(chain_residues.items()):
        seq = "".join(r[1] for r in residues)
        rows.append(
            {
                "chain_id": chain,
                "n_residues": len(residues),
                "sequence": seq,
                "residue_start": residues[0][0] if residues else "",
                "residue_end": residues[-1][0] if residues else "",
            }
        )
    return pd.DataFrame(rows)


def annotate_chains(topology: Path, expected_peptide: str = "") -> pd.DataFrame:
    seqs = pdb_chain_sequences(topology)
    if seqs.empty:
        return seqs
    annotations = []
    peptide_chains = set(seqs.loc[seqs["sequence"].eq(expected_peptide), "chain_id"]) if expected_peptide else set()
    if not peptide_chains:
        short = seqs[seqs["n_residues"].between(8, 15)].copy()
        peptide_chains = set(short.sort_values("n_residues").head(1)["chain_id"])
    mhc_candidates = seqs[seqs["n_residues"].ge(250)].sort_values("n_residues", ascending=False)
    mhc_chain = str(mhc_candidates.iloc[0]["chain_id"]) if not mhc_candidates.empty else ""
    b2m_candidates = seqs[seqs["n_residues"].between(80, 120)]
    b2m_chain = str(b2m_candidates.iloc[0]["chain_id"]) if not b2m_candidates.empty else ""
    for _, r in seqs.iterrows():
        chain = str(r["chain_id"])
        role = "other_protein"
        if chain in peptide_chains:
            role = "peptide"
        elif chain == mhc_chain:
            role = "mhc_heavy_chain"
        elif chain == b2m_chain:
            role = "beta2m"
        elif 150 <= int(r["n_residues"]) <= 260:
            role = "tcr_candidate"
        annotations.append({**r.to_dict(), "expected_peptide": expected_peptide, "role": role, "peptide_sequence_match": bool(expected_peptide and r["sequence"] == expected_peptide)})
    return pd.DataFrame(annotations)


def chain_ids(annotation: pd.DataFrame, role: str) -> list[str]:
    if annotation.empty or "role" not in annotation:
        return []
    return annotation.loc[annotation["role"].eq(role), "chain_id"].astype(str).tolist()


def mdtraj_chain_atom_indices(topology, chain_ids_: Iterable[str], heavy_only: bool = False, backbone_only: bool = False) -> list[int]:
    ids = set(chain_ids_)
    out = []
    for atom in topology.atoms:
        cid = str(atom.residue.chain.chain_id)
        if cid not in ids:
            continue
        if heavy_only and atom.element.symbol == "H":
            continue
        if backbone_only and atom.name not in {"N", "CA", "C", "O"}:
            continue
        out.append(atom.index)
    return out


def load_mdtraj(run: RunRecord, stride: int = 1):
    import mdtraj as md

    if not run.trajectory or not run.topology or not run.trajectory.exists() or not run.topology.exists():
        return None
    return md.load(str(run.trajectory), top=str(run.topology), stride=stride)


def write_markdown_table(df: pd.DataFrame, cols: list[str], max_rows: int = 30) -> str:
    if df.empty:
        return "_No rows._"
    sub = df[cols].head(max_rows).copy()
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for _, row in sub.iterrows():
        vals = []
        for c in cols:
            v = row.get(c, "")
            if isinstance(v, float):
                vals.append(f"{v:.4g}")
            else:
                vals.append(str(v))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)
