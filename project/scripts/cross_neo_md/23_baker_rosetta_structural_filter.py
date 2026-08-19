#!/usr/bin/env python3
"""Baker/Rosetta-style cheap structural filter for CROSS-Neo candidates.

The intended order is DL first, then cheap structural checks, then MD.
This script therefore does three things:

1. Audit local and RunPod availability for Rosetta/PyRosetta/ProteinMPNN/RF* tools.
2. Build command manifests for Rosetta, ProteinMPNN, and structure-ensemble tools.
3. Run a local fallback static-interface score from available PDB files when Rosetta
   is not installed.

The fallback score is not an energy and not immunogenicity evidence. It is a
triage layer to decide which structures deserve expensive modeling/MD.
"""

from __future__ import annotations

import json
import os
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1


REPO = Path(__file__).resolve().parents[3]
MD_OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
OUT = MD_OUT / "baker_rosetta_filter"
FIG = MD_OUT / "figures"
DL_FIRST = MD_OUT / "dl_first_funnel"
COUNTER = MD_OUT / "counterfactual_design"

RUNPOD_KEY = Path.home() / ".runpod/ssh/RunPod-Key-Go"
RUNPOD_CONFIG = Path.home() / ".runpod/config.toml"

COMMANDS = [
    "rosetta_scripts.default.linuxgccrelease",
    "rosetta_scripts",
    "score_jd2.default.linuxgccrelease",
    "score_jd2",
    "relax.default.linuxgccrelease",
    "relax",
    "InterfaceAnalyzer.default.linuxgccrelease",
    "InterfaceAnalyzer",
    "protein_mpnn_run.py",
    "boltz",
    "chai-lab",
    "chai",
    "python",
    "python3",
]
MODULES = ["pyrosetta", "Bio", "MDAnalysis", "openmm", "prody", "pennylane", "qiskit", "torch"]


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def safe_float(x: object, default: float = 0.0) -> float:
    try:
        if pd.isna(x):
            return default
        return float(x)
    except Exception:
        return default


def norm_hla(x: object) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip().upper()
    return s if s.startswith("HLA-") else f"HLA-{s}"


def key(peptide: object, hla: object) -> str:
    return f"{str(peptide).upper().strip()}|{norm_hla(hla)}"


def run(cmd: list[str], timeout: int = 15) -> tuple[int, str]:
    try:
        p = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=timeout, check=False)
        return p.returncode, p.stdout.strip()
    except Exception as e:
        return 999, f"{type(e).__name__}: {e}"


def audit_local() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for c in COMMANDS:
        path = shutil.which(c)
        rows.append(
            {
                "environment": "local",
                "host": os.uname().nodename,
                "tool_type": "command",
                "tool": c,
                "available": bool(path),
                "path_or_version": path or "",
                "notes": "",
            }
        )
    for m in MODULES:
        code = f"import {m}; print(getattr({m}, '__version__', 'OK'))"
        rc, out = run(["python", "-c", code], timeout=10)
        rows.append(
            {
                "environment": "local",
                "host": os.uname().nodename,
                "tool_type": "python_module",
                "tool": m,
                "available": rc == 0,
                "path_or_version": out.splitlines()[-1] if out else "",
                "notes": "" if rc == 0 else out[:200],
            }
        )
    return rows


def runpod_api_key() -> str:
    if not RUNPOD_CONFIG.exists():
        return ""
    for line in RUNPOD_CONFIG.read_text().splitlines():
        if line.strip().startswith("apikey"):
            return line.split('"')[1] if '"' in line else line.split("=", 1)[1].strip()
    return ""


def list_runpod_pods() -> pd.DataFrame:
    api_key = runpod_api_key()
    if not api_key:
        return pd.DataFrame()
    query = "query { myself { pods { id name desiredStatus runtime { ports { ip publicPort privatePort type } } } } }"
    rc, out = run(
        [
            "curl",
            "-s",
            "-X",
            "POST",
            "https://api.runpod.io/graphql",
            "-H",
            f"Authorization: Bearer {api_key}",
            "-H",
            "Content-Type: application/json",
            "-d",
            json.dumps({"query": query}),
        ],
        timeout=20,
    )
    if rc != 0 or not out:
        return pd.DataFrame([{"pod_query_status": "failed", "raw": out}])
    try:
        payload = json.loads(out)
    except Exception:
        return pd.DataFrame([{"pod_query_status": "json_parse_failed", "raw": out[:500]}])
    rows = []
    for p in payload.get("data", {}).get("myself", {}).get("pods", []):
        ssh_ip = ""
        ssh_port = ""
        for port in (p.get("runtime") or {}).get("ports", []) or []:
            if port.get("privatePort") == 22 and port.get("type") == "tcp":
                ssh_ip = port.get("ip", "")
                ssh_port = port.get("publicPort", "")
        rows.append(
            {
                "pod_id": p.get("id", ""),
                "pod_name": p.get("name", ""),
                "desired_status": p.get("desiredStatus", ""),
                "ssh_ip": ssh_ip,
                "ssh_port": ssh_port,
                "ssh_available": bool(ssh_ip and ssh_port),
            }
        )
    return pd.DataFrame(rows)


def audit_runpod(pods: pd.DataFrame) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    if pods.empty or not RUNPOD_KEY.exists():
        return rows
    remote_cmd = (
        "hostname; "
        "nvidia-smi --query-gpu=name,memory.total,memory.used --format=csv,noheader 2>/dev/null || true; "
        "for c in "
        + " ".join(shlex.quote(c) for c in COMMANDS)
        + "; do p=$(command -v $c 2>/dev/null || true); [ -n \"$p\" ] && echo CMD:$c:$p || echo CMD:$c:MISS; done; "
        "python3 - <<'PY'\n"
        f"mods={MODULES!r}\n"
        "for m in mods:\n"
        "    try:\n"
        "        mod=__import__(m); print('MOD:%s:OK:%s' % (m, getattr(mod, '__version__', '')))\n"
        "    except Exception as e:\n"
        "        print('MOD:%s:MISS:%s' % (m, type(e).__name__))\n"
        "PY"
    )
    for _, p in pods.iterrows():
        if not p.get("ssh_available", False):
            continue
        cmd = [
            "ssh",
            "-o",
            "StrictHostKeyChecking=no",
            "-o",
            "BatchMode=yes",
            "-o",
            "ConnectTimeout=8",
            "-i",
            str(RUNPOD_KEY),
            "-p",
            str(p["ssh_port"]),
            f"root@{p['ssh_ip']}",
            remote_cmd,
        ]
        rc, out = run(cmd, timeout=35)
        lines = out.splitlines()
        gpu = "|".join([line for line in lines if "NVIDIA" in line])[:250]
        if rc != 0:
            rows.append(
                {
                    "environment": "runpod",
                    "host": p["pod_name"],
                    "tool_type": "ssh",
                    "tool": "ssh_connect",
                    "available": False,
                    "path_or_version": "",
                    "notes": out[:300],
                    "gpu": gpu,
                }
            )
            continue
        for line in lines:
            if line.startswith("CMD:"):
                _, tool, path = line.split(":", 2)
                rows.append(
                    {
                        "environment": "runpod",
                        "host": p["pod_name"],
                        "tool_type": "command",
                        "tool": tool,
                        "available": path != "MISS",
                        "path_or_version": "" if path == "MISS" else path,
                        "notes": "",
                        "gpu": gpu,
                    }
                )
            elif line.startswith("MOD:"):
                _, tool, status, version = (line.split(":", 3) + [""])[:4]
                rows.append(
                    {
                        "environment": "runpod",
                        "host": p["pod_name"],
                        "tool_type": "python_module",
                        "tool": tool,
                        "available": status == "OK",
                        "path_or_version": version if status == "OK" else "",
                        "notes": "" if status == "OK" else version,
                        "gpu": gpu,
                    }
                )
    return rows


def chain_lookup_from_annotations() -> dict[str, dict[str, str]]:
    ann = read_tsv(MD_OUT / "complex_chain_annotation.tsv")
    lookup: dict[str, dict[str, str]] = {}
    for _, r in ann.iterrows():
        pdb = str(r.get("topology", ""))
        if pdb:
            lookup[str(Path(pdb).resolve())] = {
                "peptide_chain": str(r.get("peptide_chains", "")),
                "mhc_chain": str(r.get("mhc_heavy_chain", "")),
                "b2m_chain": str(r.get("beta2m_chain", "")),
                "tcr_chains": str(r.get("tcr_candidate_chains", "")),
            }
    pos = read_tsv(COUNTER / "positive_control_pdb_extraction.tsv")
    for _, r in pos.iterrows():
        pdb = str(r.get("extracted_pdb", ""))
        if pdb:
            lookup[str(Path(pdb).resolve())] = {
                "peptide_chain": str(r.get("peptide_chain", "")),
                "mhc_chain": str(r.get("mhc_chain", "")),
                "b2m_chain": str(r.get("b2m_chain", "")),
                "tcr_chains": str(r.get("tcr_chains", "")),
            }
    return lookup


def infer_chains(pdb: Path, expected_peptide: str, lookup: dict[str, dict[str, str]]) -> dict[str, str]:
    rp = str(pdb.resolve())
    if rp in lookup:
        return lookup[rp]
    parser = PDBParser(QUIET=True)
    try:
        model = next(parser.get_structure(pdb.stem, str(pdb)).get_models())
    except Exception:
        name = pdb.name.upper()
        if "6VRN" in name or "HMTEVVRHC" in name:
            return {"peptide_chain": "P", "mhc_chain": "A", "b2m_chain": "B", "tcr_chains": "D|E"}
        if "6UON" in name or "GADGVGKSAL" in name:
            return {"peptide_chain": "F", "mhc_chain": "D", "b2m_chain": "E", "tcr_chains": "G|H"}
        return {"peptide_chain": "", "mhc_chain": "", "b2m_chain": "", "tcr_chains": ""}
    seqs = {}
    for ch in model:
        residues = [res for res in ch if res.id[0] == " "]
        aa = []
        for res in residues:
            try:
                aa.append(seq1(res.resname))
            except Exception:
                aa.append("X")
        seqs[ch.id] = "".join(aa)
    peptide_chain = next((c for c, s in seqs.items() if s == expected_peptide), "")
    if not peptide_chain and expected_peptide:
        peptide_chain = next((c for c, s in seqs.items() if expected_peptide in s), "")
    if not peptide_chain:
        peptide_chain = min(seqs, key=lambda c: len(seqs[c])) if seqs else ""
    rest = [c for c in seqs if c != peptide_chain]
    mhc_chain = max(rest, key=lambda c: len(seqs[c]), default="")
    b2m_chain = next((c for c in rest if c != mhc_chain and 75 <= len(seqs[c]) <= 130), "")
    tcr_chains = "|".join([c for c in rest if c not in {mhc_chain, b2m_chain}])
    return {"peptide_chain": peptide_chain, "mhc_chain": mhc_chain, "b2m_chain": b2m_chain, "tcr_chains": tcr_chains}


def build_manifest() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    dl = read_tsv(DL_FIRST / "dl_first_md_escalation_candidates.tsv")
    openmm = read_tsv(COUNTER / "openmm_ready_structure_jobs.tsv")
    lookup = chain_lookup_from_annotations()
    if openmm.empty:
        return pd.DataFrame()
    wanted = set(dl["row_id"].astype(str)) if not dl.empty else set(openmm["row_id"].astype(str))
    openmm = openmm[openmm["row_id"].astype(str).isin(wanted)].copy()
    for i, r in openmm.iterrows():
        pdb = Path(str(r.get("input_template_pdb", "")))
        chains = infer_chains(pdb, str(r.get("sequence", r.get("peptide", ""))).upper(), lookup) if pdb.exists() else {}
        rows.append(
            {
                "filter_id": f"BAKERFILTER_{len(rows):04d}",
                "prep_id": r.get("prep_id", ""),
                "row_id": r.get("row_id", ""),
                "peptide": r.get("peptide", ""),
                "hla_4digit": r.get("hla_4digit", ""),
                "candidate_key": key(r.get("peptide", ""), r.get("hla_4digit", "")),
                "control_type": r.get("control_type", ""),
                "complex_kind": r.get("complex_kind", ""),
                "sequence": r.get("sequence", ""),
                "sequence_status": r.get("sequence_status", ""),
                "input_pdb": str(pdb),
                "pdb_exists": pdb.exists(),
                "peptide_chain": chains.get("peptide_chain", ""),
                "mhc_chain": chains.get("mhc_chain", ""),
                "b2m_chain": chains.get("b2m_chain", ""),
                "tcr_chains": chains.get("tcr_chains", ""),
                "dl_first_source": "md_escalation_after_cheap_gates",
                "claim_status": "diagnostic_only_structural_filter_not_immunogenicity_evidence",
            }
        )
    return pd.DataFrame(rows)


@dataclass
class AtomSet:
    coords: np.ndarray
    residue_keys: list[str]


def atomset(model, chain_ids: list[str]) -> AtomSet:
    coords = []
    keys = []
    for chain_id in chain_ids:
        if not chain_id:
            continue
        try:
            chain = model[chain_id]
        except Exception:
            continue
        for res in chain:
            if res.id[0] != " ":
                continue
            rkey = f"{chain_id}:{res.resname}{res.id[1]}"
            for atom in res:
                element = (atom.element or atom.name[0]).upper()
                if element == "H":
                    continue
                coords.append(atom.coord.astype(float))
                keys.append(rkey)
    if not coords:
        return AtomSet(np.zeros((0, 3), dtype=float), [])
    return AtomSet(np.vstack(coords), keys)


def contact_metrics(a: AtomSet, b: AtomSet, contact_cutoff: float = 4.0, clash_cutoff: float = 2.0) -> dict[str, float]:
    if len(a.coords) == 0 or len(b.coords) == 0:
        return {
            "heavy_atom_contacts_4A": 0,
            "residue_pair_contacts_4A": 0,
            "chain_clashes_2A": 0,
            "min_distance_A": np.nan,
        }
    d = np.linalg.norm(a.coords[:, None, :] - b.coords[None, :, :], axis=2)
    contact_idx = np.argwhere(d <= contact_cutoff)
    residue_pairs = {f"{a.residue_keys[i]}--{b.residue_keys[j]}" for i, j in contact_idx}
    return {
        "heavy_atom_contacts_4A": int(contact_idx.shape[0]),
        "residue_pair_contacts_4A": int(len(residue_pairs)),
        "chain_clashes_2A": int((d <= clash_cutoff).sum()),
        "min_distance_A": float(np.min(d)),
    }


def score_manifest(manifest: pd.DataFrame) -> pd.DataFrame:
    parser = PDBParser(QUIET=True)
    rows = []
    for _, r in manifest.iterrows():
        pdb = Path(str(r.get("input_pdb", "")))
        row = r.to_dict()
        if not pdb.exists():
            row.update({"fallback_status": "missing_pdb", "fallback_structural_score": 0.0})
            rows.append(row)
            continue
        try:
            model = next(parser.get_structure(pdb.stem, str(pdb)).get_models())
        except Exception as e:
            row.update({"fallback_status": f"parse_failed:{type(e).__name__}", "fallback_structural_score": 0.0})
            rows.append(row)
            continue
        pep = atomset(model, [str(r.get("peptide_chain", ""))])
        mhc = atomset(model, [str(r.get("mhc_chain", ""))])
        inferred_tcr_ids = [x for x in str(r.get("tcr_chains", "")).split("|") if x and x.lower() != "nan"]
        score_as_tcr_pmhc = str(r.get("complex_kind", "")) == "TCR-pMHC"
        tcr_ids = inferred_tcr_ids if score_as_tcr_pmhc else []
        tcr = atomset(model, tcr_ids)
        pmhc = contact_metrics(pep, mhc)
        tcrpep = contact_metrics(pep, tcr)
        tcrmhc = contact_metrics(tcr, mhc)

        pep_mhc_res = pmhc["residue_pair_contacts_4A"]
        tcr_pep_res = tcrpep["residue_pair_contacts_4A"]
        clashes = pmhc["chain_clashes_2A"] + tcrpep["chain_clashes_2A"] + tcrmhc["chain_clashes_2A"]
        pmhc_proxy = min(1.0, pep_mhc_res / max(8.0, 2.0 * len(str(r.get("sequence", "")))))
        tcr_proxy = min(1.0, tcr_pep_res / 10.0) if score_as_tcr_pmhc and tcr_ids else 0.0
        clash_penalty = min(0.35, clashes / 500.0)
        if score_as_tcr_pmhc:
            score = 0.60 * pmhc_proxy + 0.35 * tcr_proxy + 0.05 * (1.0 if pdb.exists() else 0.0) - clash_penalty
        else:
            score = 0.90 * pmhc_proxy + 0.05 * (1.0 if pdb.exists() else 0.0) - clash_penalty
        score = float(np.clip(score, 0, 1))

        row.update(
            {
                "fallback_status": "scored",
                "pmhc_heavy_atom_contacts_4A": pmhc["heavy_atom_contacts_4A"],
                "pmhc_residue_pair_contacts_4A": pep_mhc_res,
                "pmhc_min_distance_A": pmhc["min_distance_A"],
                "tcr_peptide_heavy_atom_contacts_4A": tcrpep["heavy_atom_contacts_4A"],
                "tcr_peptide_residue_pair_contacts_4A": tcr_pep_res,
                "tcr_peptide_min_distance_A": tcrpep["min_distance_A"],
                "tcr_mhc_residue_pair_contacts_4A": tcrmhc["residue_pair_contacts_4A"],
                "cross_chain_clashes_2A": clashes,
                "pmhc_static_proxy_score": pmhc_proxy,
                "tcr_static_proxy_score": tcr_proxy,
                "fallback_structural_score": score,
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def write_command_templates(manifest: pd.DataFrame) -> None:
    cmd_dir = OUT / "commands"
    cmd_dir.mkdir(parents=True, exist_ok=True)
    score_lines = ["#!/usr/bin/env bash", "set -euo pipefail", "mkdir -p baker_rosetta_score_out"]
    relax_lines = ["#!/usr/bin/env bash", "set -euo pipefail", "mkdir -p baker_rosetta_relax_out"]
    interface_lines = ["#!/usr/bin/env bash", "set -euo pipefail", "mkdir -p baker_rosetta_interface_out"]
    mpnn_lines = ["#!/usr/bin/env bash", "set -euo pipefail", "mkdir -p proteinmpnn_out"]
    ensemble_lines = ["#!/usr/bin/env bash", "set -euo pipefail", "mkdir -p structure_ensemble_out"]
    runpod_lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "# Non-destructive template. Edit POD_HOST/POD_PORT if RunPod ports change.",
        "POD_HOST=${POD_HOST:-135.84.176.142}",
        "POD_PORT=${POD_PORT:-20878}",
        "KEY=${KEY:-$HOME/.runpod/ssh/RunPod-Key-Go}",
        "REMOTE=${REMOTE:-/workspace/cross_neo_baker_filter}",
        "ssh -i \"$KEY\" -p \"$POD_PORT\" -o StrictHostKeyChecking=no root@\"$POD_HOST\" \"mkdir -p $REMOTE\"",
        f"scp -i \"$KEY\" -P \"$POD_PORT\" -o StrictHostKeyChecking=no {shlex.quote(str(OUT / 'baker_rosetta_filter_manifest.tsv'))} root@\"$POD_HOST\":\"$REMOTE/\"",
        "# Copy PDBs separately if needed; large trajectory files are intentionally not copied here.",
    ]
    for _, r in manifest.iterrows():
        pdb = str(r.get("input_pdb", ""))
        if not pdb or not Path(pdb).exists():
            continue
        q = shlex.quote(pdb)
        tag = str(r.get("filter_id", "job"))
        score_lines.append(f"score_jd2.default.linuxgccrelease -s {q} -out:file:scorefile baker_rosetta_score_out/{tag}.sc || true")
        relax_lines.append(
            f"relax.default.linuxgccrelease -s {q} -relax:fast -nstruct 3 -out:path:all baker_rosetta_relax_out/{tag} || true"
        )
        interface_lines.append(
            "# InterfaceAnalyzer chain grouping may need manual adjustment for multi-chain TCR-pMHC complexes."
        )
        interface_lines.append(
            f"InterfaceAnalyzer.default.linuxgccrelease -s {q} -out:file:scorefile baker_rosetta_interface_out/{tag}.sc || true"
        )
        mpnn_lines.append(f"python protein_mpnn_run.py --pdb_path {q} --out_folder proteinmpnn_out/{tag} --num_seq_per_target 16 --sampling_temp 0.1 || true")
        ensemble_lines.append(f"# boltz predict {q} --out_dir structure_ensemble_out/{tag}_boltz")
        ensemble_lines.append(f"# chai-lab fold --input {q} --output-dir structure_ensemble_out/{tag}_chai")
    for name, lines in [
        ("run_rosetta_score_templates.sh", score_lines),
        ("run_rosetta_relax_templates.sh", relax_lines),
        ("run_rosetta_interface_templates.sh", interface_lines),
        ("run_proteinmpnn_templates.sh", mpnn_lines),
        ("run_structure_ensemble_templates.sh", ensemble_lines),
        ("runpod_baker_filter_sync_template.sh", runpod_lines),
    ]:
        path = cmd_dir / name
        path.write_text("\n".join(lines) + "\n")
        path.chmod(0o755)


def plot_outputs(scores: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    if scores.empty:
        return
    show = scores.sort_values("fallback_structural_score", ascending=False).head(20).copy()
    label = show["sequence"].astype(str) + "\n" + show["complex_kind"].astype(str) + "/" + show["control_type"].astype(str)
    fig, ax = plt.subplots(figsize=(11.5, 6.2))
    ax.barh(label, show["fallback_structural_score"], color="#4d96ff")
    ax.set_xlim(0, 1)
    ax.invert_yaxis()
    ax.set_xlabel("Fallback static structural score")
    ax.set_title("Baker/Rosetta cheap-filter fallback scores")
    ax.grid(axis="x", alpha=0.25)
    fig.savefig(FIG / "fig_md20_baker_rosetta_fallback_scores.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md20_baker_rosetta_fallback_scores.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.8, 5.6))
    ax.scatter(
        scores["pmhc_residue_pair_contacts_4A"],
        scores["tcr_peptide_residue_pair_contacts_4A"],
        s=70,
        c=scores["fallback_structural_score"],
        cmap="viridis",
        edgecolor="#111",
    )
    for _, r in scores.head(12).iterrows():
        ax.annotate(str(r.get("sequence", ""))[:12], (r["pmhc_residue_pair_contacts_4A"], r["tcr_peptide_residue_pair_contacts_4A"]), fontsize=7)
    ax.set_xlabel("pMHC residue-pair contacts <=4 A")
    ax.set_ylabel("TCR-peptide residue-pair contacts <=4 A")
    ax.set_title("Static interface contact map for structural triage")
    ax.grid(alpha=0.25)
    fig.savefig(FIG / "fig_md21_baker_contact_triage_map.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md21_baker_contact_triage_map.pdf", bbox_inches="tight")
    plt.close(fig)


def write_reports(tool_audit: pd.DataFrame, pods: pd.DataFrame, manifest: pd.DataFrame, scores: pd.DataFrame) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    tool_audit.to_csv(OUT / "baker_rosetta_tool_audit.tsv", sep="\t", index=False)
    pods.to_csv(OUT / "runpod_pod_status.tsv", sep="\t", index=False)
    manifest.to_csv(OUT / "baker_rosetta_filter_manifest.tsv", sep="\t", index=False)
    scores.to_csv(OUT / "baker_rosetta_fallback_interface_scores.tsv", sep="\t", index=False)
    available = tool_audit[tool_audit["available"].astype(bool)]
    top = scores.sort_values("fallback_structural_score", ascending=False).head(20) if not scores.empty else pd.DataFrame()
    lines = [
        "# Baker/Rosetta Cheap Structural Filter",
        "",
        "## Purpose",
        "",
        "This is the structural layer after the DL-first screen and before expensive MD. It prepares Rosetta/ProteinMPNN/RF-style jobs and runs a static contact fallback when Rosetta is not installed.",
        "",
        "## Tool Availability",
        "",
        available[["environment", "host", "tool_type", "tool", "path_or_version"]].to_markdown(index=False) if not available.empty else "No structural-design tools detected.",
        "",
        "## RunPod Status",
        "",
        pods.to_markdown(index=False) if not pods.empty else "RunPod status unavailable.",
        "",
        "## Manifest Summary",
        "",
        manifest.groupby(["candidate_key", "control_type", "complex_kind"]).size().reset_index(name="n").to_markdown(index=False) if not manifest.empty else "No manifest rows.",
        "",
        "## Fallback Static Interface Scores",
        "",
        top[
            [
                "filter_id",
                "sequence",
                "candidate_key",
                "control_type",
                "complex_kind",
                "fallback_structural_score",
                "pmhc_residue_pair_contacts_4A",
                "tcr_peptide_residue_pair_contacts_4A",
                "cross_chain_clashes_2A",
            ]
        ].to_markdown(index=False)
        if not top.empty
        else "No fallback scores.",
        "",
        "## Interpretation Boundary",
        "",
        "- Rosetta/RF/ProteinMPNN scores are structural triage signals, not immunogenicity proof.",
        "- Static fallback contacts are weaker than Rosetta energies and much weaker than replicated MD.",
        "- TCR-discordant main-model positives should be treated as false-positive audit/control candidates.",
        "- MD escalation should remain limited to candidates that survived DL uncertainty and TCR concordance.",
    ]
    (OUT / "baker_rosetta_filter_report.md").write_text("\n".join(lines) + "\n")
    (OUT / "baker_rosetta_summary.json").write_text(
        json.dumps(
            {
                "n_manifest_rows": int(len(manifest)),
                "n_scored_rows": int((scores.get("fallback_status", pd.Series(dtype=str)) == "scored").sum()) if not scores.empty else 0,
                "available_tools": available[["environment", "host", "tool"]].to_dict("records"),
                "top_scores": top[["filter_id", "sequence", "candidate_key", "control_type", "complex_kind", "fallback_structural_score"]].head(10).to_dict("records") if not top.empty else [],
            },
            indent=2,
        )
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    pods = list_runpod_pods()
    tool_audit = pd.DataFrame(audit_local() + audit_runpod(pods))
    manifest = build_manifest()
    write_command_templates(manifest)
    scores = score_manifest(manifest)
    plot_outputs(scores)
    write_reports(tool_audit, pods, manifest, scores)
    print("[baker-rosetta-filter]", (OUT / "baker_rosetta_summary.json").read_text())


if __name__ == "__main__":
    main()
