#!/usr/bin/env python3
"""Prepare a diverse OpenMM launch pack for P0 neoantigen MD follow-up.

The pack separates:
- short template-diversity screens,
- full 10 ns replicate/control jobs that are ready now,
- jobs that should wait for HMTEVVRHC completion/sync,
- WT/decoy rows that remain blocked until sequence/structure curation.

It does not launch GPU work.
"""

from __future__ import annotations

import hashlib
import json
import shlex
import shutil
from pathlib import Path

import pandas as pd

from common_md import OUT, REPO, annotate_chains, safe_id, write_markdown_table


DESIGN = OUT / "counterfactual_design"
DIVERSE = OUT / "diverse_simulation_plan"
INPUTS = DIVERSE / "inputs"
LOGS = DIVERSE / "logs"

P0_PACKAGE = (
    REPO
    / "project/results/cross_neo_v2_sota_sprint_2026_05_09/tcr_extension/md_escalation/p0_structures/openmm_p0_replicate_screen_package"
)
RUNNER_SRC = P0_PACKAGE / "run_openmm_pilot.py"
ANALYZER_SRC = P0_PACKAGE / "analyze_openmm_pilot.py"
SCREEN_MANIFEST = P0_PACKAGE / "openmm_p0_replicate_screen_manifest.tsv"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def has_value(value: object) -> bool:
    if pd.isna(value):
        return False
    return str(value).strip().lower() not in {"", "nan", "none", "na"}


def ns_label(ns: float) -> str:
    return str(ns).replace(".", "p") + "ns"


def file_hash(path: Path) -> str:
    h = hashlib.sha1()
    h.update(str(path).encode())
    h.update(str(path.stat().st_size if path.exists() else 0).encode())
    return h.hexdigest()[:8]


def copy_input(path: Path) -> str:
    if not path.exists():
        return ""
    INPUTS.mkdir(parents=True, exist_ok=True)
    out = INPUTS / f"{path.stem}_{file_hash(path)}{path.suffix}"
    if not out.exists() or out.stat().st_size != path.stat().st_size:
        shutil.copy2(path, out)
    return f"inputs/{out.name}"


def status_index() -> dict[str, dict]:
    status = read_tsv(OUT / "md_status_summary.tsv")
    out: dict[str, dict] = {}
    if status.empty:
        return out
    for (peptide, hla), sub in status.groupby(["peptide", "hla"], dropna=False):
        key = f"{peptide}/{hla}"
        complete_10ns = sub[
            pd.to_numeric(sub.get("target_ns", pd.Series(dtype=float)), errors="coerce").ge(9.0)
            & pd.to_numeric(sub.get("completion_fraction", pd.Series(dtype=float)), errors="coerce").ge(0.99)
        ]
        best = sub.copy()
        best["_completion"] = pd.to_numeric(best.get("completion_fraction", pd.Series(dtype=float)), errors="coerce").fillna(0)
        best = best.sort_values("_completion", ascending=False)
        top = best.iloc[0].to_dict() if not best.empty else {}
        out[key] = {
            "complete_10ns_count": int(len(complete_10ns)),
            "best_completion_fraction": float(top.get("completion_fraction", 0) or 0),
            "best_run_id": top.get("run_id", ""),
            "best_time_ps": float(top.get("time_ps", 0) or 0),
            "best_source": top.get("source", ""),
        }
    return out


def completed_run_ids() -> set[str]:
    status = read_tsv(OUT / "md_status_summary.tsv")
    if status.empty:
        return set()
    completion = pd.to_numeric(status.get("completion_fraction", pd.Series(dtype=float)), errors="coerce").fillna(0)
    return set(status.loc[completion.ge(0.99), "run_id"].astype(str))


def chain_summary(local_input: str, expected_peptide: str) -> tuple[dict, str, str, str, str]:
    path = DIVERSE / local_input
    ann = annotate_chains(path, expected_peptide)
    if ann.empty:
        return (
            {
                "local_input": local_input,
                "expected_peptide": expected_peptide,
                "peptide_chain": "",
                "mhc_chain": "",
                "b2m_chain": "",
                "tcr_chains": "",
                "quality_flag": "no_parseable_protein_chains",
            },
            "",
            "",
            "",
            "no_parseable_protein_chains",
        )
    peptide = ann.loc[ann["role"].eq("peptide"), "chain_id"].astype(str).tolist()
    mhc = ann.loc[ann["role"].eq("mhc_heavy_chain"), "chain_id"].astype(str).tolist()
    b2m = ann.loc[ann["role"].eq("beta2m"), "chain_id"].astype(str).tolist()
    tcr = ann.loc[ann["role"].eq("tcr_candidate"), "chain_id"].astype(str).tolist()
    flags: list[str] = []
    if expected_peptide and not bool(ann.get("peptide_sequence_match", pd.Series(dtype=bool)).any()):
        flags.append("expected_peptide_not_exact")
    if not peptide:
        flags.append("missing_peptide_chain")
    if not mhc:
        flags.append("missing_mhc_chain")
    if not b2m:
        flags.append("missing_beta2m")
    if len(tcr) < 2:
        flags.append("tcr_chain_count_lt_2")
    flag = ";".join(flags) if flags else "ok"
    row = {
        "local_input": local_input,
        "expected_peptide": expected_peptide,
        "peptide_chain": "|".join(peptide),
        "mhc_chain": "|".join(mhc[:1]),
        "b2m_chain": "|".join(b2m[:1]),
        "tcr_chains": "|".join(tcr[:2]),
        "n_chains": int(len(ann)),
        "quality_flag": flag,
    }
    return row, "|".join(peptide[:1]), "|".join(mhc[:1]), "|".join(tcr[:2]), flag


def make_seeded_runner() -> None:
    """Copy the existing runner and add explicit seed support in the launch pack."""
    text = RUNNER_SRC.read_text()
    if "--random-seed" not in text:
        text = text.replace(
            '    ap.add_argument("--minimize-iterations", type=int, default=500)\n',
            '    ap.add_argument("--minimize-iterations", type=int, default=500)\n'
            '    ap.add_argument("--random-seed", type=int, default=None)\n',
        )
        text = text.replace(
            "    if args.mode == \"explicit\":\n"
            "        system.addForce(MonteCarloBarostat(1.0 * unit.atmosphere, args.temperature_k * unit.kelvin, 25))\n\n"
            "    integrator = LangevinMiddleIntegrator(args.temperature_k * unit.kelvin, 1 / unit.picosecond, dt)\n",
            "    if args.mode == \"explicit\":\n"
            "        barostat = MonteCarloBarostat(1.0 * unit.atmosphere, args.temperature_k * unit.kelvin, 25)\n"
            "        if args.random_seed is not None:\n"
            "            barostat.setRandomNumberSeed(int(args.random_seed) + 100_000)\n"
            "        system.addForce(barostat)\n\n"
            "    integrator = LangevinMiddleIntegrator(args.temperature_k * unit.kelvin, 1 / unit.picosecond, dt)\n"
            "    if args.random_seed is not None:\n"
            "        integrator.setRandomNumberSeed(int(args.random_seed))\n",
        )
        text = text.replace(
            "    with minimized_pdb.open(\"w\") as handle:\n"
            "        PDBFile.writeFile(modeller.topology, state.getPositions(), handle)\n\n"
            "    simulation.reporters.append(DCDReporter(str(args.outdir / \"trajectory.dcd\"), args.report_steps))\n",
            "    with minimized_pdb.open(\"w\") as handle:\n"
            "        PDBFile.writeFile(modeller.topology, state.getPositions(), handle)\n\n"
            "    if args.random_seed is not None:\n"
            "        simulation.context.setVelocitiesToTemperature(args.temperature_k * unit.kelvin, int(args.random_seed))\n"
            "    else:\n"
            "        simulation.context.setVelocitiesToTemperature(args.temperature_k * unit.kelvin)\n\n"
            "    simulation.reporters.append(DCDReporter(str(args.outdir / \"trajectory.dcd\"), args.report_steps))\n",
        )
        text = text.replace(
            '        "temperature_k": args.temperature_k,\n',
            '        "temperature_k": args.temperature_k,\n'
            '        "random_seed": args.random_seed,\n',
        )
    (DIVERSE / "run_openmm_pilot.py").write_text(text)
    shutil.copy2(ANALYZER_SRC, DIVERSE / "analyze_openmm_pilot.py")


def add_job(
    rows: list[dict],
    audits: list[dict],
    *,
    tier: str,
    launch_group: str,
    source_kind: str,
    peptide: str,
    hla: str,
    sequence: str,
    control_type: str,
    complex_scope: str,
    source_pdb: Path,
    ns: float,
    timestep_fs: float,
    report_steps: int,
    replicate_idx: int,
    seed: int,
    immediate_launch: bool,
    hold_reason: str,
    rationale: str,
    analysis_targets: str,
    priority: int,
    peptide_chain_hint: str = "",
    mhc_chain_hint: str = "",
    tcr_chains_hint: str = "",
) -> None:
    local_input = copy_input(source_pdb)
    if not local_input:
        return
    audit, peptide_chain, mhc_chain, tcr_chains, flag = chain_summary(local_input, sequence or peptide)
    audits.append({**audit, "source_pdb": str(source_pdb), "source_kind": source_kind})
    peptide_chain = peptide_chain_hint or peptide_chain
    mhc_chain = mhc_chain_hint or mhc_chain
    tcr_chains = tcr_chains_hint or tcr_chains
    base = safe_id(f"{tier}_{peptide}_{hla}_{control_type}_{complex_scope}_{source_pdb.stem}")
    job_id = f"{base}_rep{replicate_idx:02d}_{ns_label(ns)}"
    command = (
        "python run_openmm_pilot.py "
        f"--input {shlex.quote(local_input)} "
        f"--outdir {shlex.quote('runs/' + job_id)} "
        "--mode explicit --platform CUDA "
        f"--ns {ns:g} --report-steps {int(report_steps)} "
        "--minimize-iterations 1000 "
        f"--timestep-fs {timestep_fs:g} --temperature-k 300 "
        f"--random-seed {int(seed)}"
    )
    analysis_command = (
        "python analyze_openmm_pilot.py "
        f"--top {shlex.quote('runs/' + job_id + '/prepared_start.pdb')} "
        f"--traj {shlex.quote('runs/' + job_id + '/trajectory.dcd')} "
        f"--peptide-chain {shlex.quote(peptide_chain)} "
        f"--mhc-chain {shlex.quote(mhc_chain)} "
        f"--tcr-chains {shlex.quote(tcr_chains)} "
        f"--peptide-sequence {shlex.quote(sequence or peptide)} "
        f"--out {shlex.quote('runs/' + job_id + '/analysis.json')}"
    )
    can_launch = flag != "no_parseable_protein_chains" and "missing_peptide_chain" not in flag and "missing_mhc_chain" not in flag
    if complex_scope == "TCR-pMHC" and "tcr_chain_count_lt_2" in flag:
        can_launch = False
        hold_reason = hold_reason or "tcr_chain_count_lt_2"
    rows.append(
        {
            "priority": priority,
            "tier": tier,
            "launch_group": launch_group,
            "job_id": job_id,
            "source_kind": source_kind,
            "peptide": peptide,
            "hla_4digit": hla,
            "sequence": sequence,
            "control_type": control_type,
            "complex_scope": complex_scope,
            "source_pdb": str(source_pdb),
            "local_input": local_input,
            "ns": ns,
            "timestep_fs": timestep_fs,
            "report_steps": report_steps,
            "replicate_idx": replicate_idx,
            "random_seed": seed,
            "immediate_launch": bool(immediate_launch and can_launch),
            "can_launch": bool(can_launch),
            "hold_reason": hold_reason if not immediate_launch or not can_launch else "",
            "chain_quality_flag": flag,
            "peptide_chain": peptide_chain,
            "mhc_chain": mhc_chain,
            "tcr_chains": tcr_chains,
            "analysis_targets": analysis_targets,
            "rationale": rationale,
            "run_command": command,
            "analysis_command": analysis_command,
            "claim_status": "simulation_prioritization_only_not_immunogenicity_proof",
        }
    )


def build_jobs() -> tuple[pd.DataFrame, pd.DataFrame]:
    DIVERSE.mkdir(parents=True, exist_ok=True)
    INPUTS.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)

    ready = read_tsv(DESIGN / "openmm_ready_structure_jobs.tsv")
    blocked = read_tsv(DESIGN / "blocked_structure_jobs.tsv")
    screens = read_tsv(SCREEN_MANIFEST)
    md_status = status_index()
    already_complete = completed_run_ids()

    rows: list[dict] = []
    audits: list[dict] = []

    if not screens.empty:
        for rank, (_, s) in enumerate(screens.iterrows(), start=1):
            if str(s.get("screen_priority", "")) != "run_short_replicate_screen":
                continue
            input_pdb = P0_PACKAGE / str(s["input_pdb"])
            peptide = str(s["target_peptide"])
            hla = str(s["hla_4digit"])
            screen_ns = float(s.get("screen_ns", 0.5))
            expected_run_id = f"{s.get('pilot_id', '')}_{ns_label(screen_ns)}"
            if expected_run_id in already_complete:
                continue
            immediate = peptide == "GADGVGKSAL"
            hold = "" if immediate else "finish_or_sync_current_HMTEVVRHC_10ns_before_more_HMTEVVRHC_GPU_work"
            add_job(
                rows,
                audits,
                tier="tier0_short_template_diversity_screen",
                launch_group="tier0_short_screens",
                source_kind="existing_p0_replicate_screen_manifest",
                peptide=peptide,
                hla=hla,
                sequence=peptide,
                control_type="mutant_alt_template",
                complex_scope="TCR-pMHC",
                source_pdb=input_pdb,
                ns=screen_ns,
                timestep_fs=float(s.get("timestep_fs", 1.0)),
                report_steps=int(s.get("report_steps", 50000)),
                replicate_idx=1,
                seed=10_000 + rank,
                immediate_launch=immediate,
                hold_reason=hold,
                rationale="short alternate-template instability screen before committing long GPU time",
                analysis_targets=str(s.get("pilot_id", "")),
                priority=rank,
                peptide_chain_hint=str(s.get("peptide_chain", "")),
                mhc_chain_hint=str(s.get("mhc_chain", "")),
                tcr_chains_hint=str(s.get("tcr_chains", "")),
            )

    if not ready.empty:
        grouped = ready.groupby(["peptide", "hla_4digit", "control_type", "sequence", "input_template_pdb"], dropna=False)
        for rank, (key, sub) in enumerate(grouped, start=101):
            peptide, hla, control_type, sequence, input_template = key
            peptide = str(peptide)
            hla = str(hla)
            control_type = str(control_type)
            sequence = str(sequence)
            source_pdb = Path(str(input_template))
            complex_kinds = sorted(set(sub["complex_kind"].astype(str)))
            scope = "TCR-pMHC" if "TCR-pMHC" in complex_kinds else "pMHC"
            status = md_status.get(f"{peptide}/{hla}", {})
            completed_10ns = int(status.get("complete_10ns_count", 0))
            target_reps = 3
            new_reps = target_reps
            if control_type == "mutant" and completed_10ns:
                new_reps = max(0, target_reps - completed_10ns)
            if peptide == "HMTEVVRHC" and control_type == "mutant" and completed_10ns == 0:
                new_reps = 2
            if new_reps == 0:
                continue

            if peptide == "GADGVGKSAL":
                immediate = control_type in {"mutant", "same_or_similar_hla_positive_control"}
                launch_group = "tier1_gad_full_replicates_controls"
                hold = ""
            else:
                immediate = False
                launch_group = "after_hmtevv_sync"
                hold = "finish_or_sync_current_HMTEVVRHC_10ns_before_duplicate_or_control_batch"

            if control_type == "same_or_similar_hla_positive_control" and peptide == "GADGVGKSAL" and scope == "TCR-pMHC":
                # The current 6ULR extraction is pMHC-only; keep it as pMHC sanity control.
                scope = "pMHC"

            for rep in range(1, new_reps + 1):
                add_job(
                    rows,
                    audits,
                    tier="tier1_full_10ns_replicate_or_control",
                    launch_group=launch_group,
                    source_kind="counterfactual_openmm_ready_structure_job",
                    peptide=peptide,
                    hla=hla,
                    sequence=sequence,
                    control_type=control_type,
                    complex_scope=scope,
                    source_pdb=source_pdb,
                    ns=10.0,
                    timestep_fs=2.0,
                    report_steps=50000,
                    replicate_idx=rep,
                    seed=20_000 + rank * 10 + rep,
                    immediate_launch=immediate,
                    hold_reason=hold,
                    rationale="3x10ns replicate/control batch for specificity and robustness triage",
                    analysis_targets=";".join(sub["batch_id"].astype(str).tolist()),
                    priority=rank * 10 + rep,
                )

    jobs = pd.DataFrame(rows).sort_values(["immediate_launch", "priority"], ascending=[False, True])
    audit = pd.DataFrame(audits).drop_duplicates(["local_input", "expected_peptide", "source_pdb"])

    if not blocked.empty:
        blocked_summary = (
            blocked.groupby(["peptide", "hla_4digit", "control_type", "complex_kind", "blocked_by"], dropna=False)
            .size()
            .reset_index(name="n")
        )
    else:
        blocked_summary = pd.DataFrame()
    blocked_summary.to_csv(DIVERSE / "blocked_wt_decoy_structure_work.tsv", sep="\t", index=False)
    return jobs, audit


def shell_script(name: str, jobs: pd.DataFrame, header: str) -> None:
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        'SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
        'cd "$SCRIPT_DIR"',
        "mkdir -p runs logs",
        "",
        f"echo {shlex.quote(header)}",
        "",
    ]
    for _, r in jobs.iterrows():
        job_id = str(r["job_id"])
        lines.extend(
            [
                f"echo '[job] {job_id}'",
                f"if [ -f {shlex.quote('runs/' + job_id + '/run_metadata.json')} ]; then",
                f"  echo '[skip] {job_id} already has run_metadata.json'",
                "else",
                f"  {r['run_command']} 2>&1 | tee {shlex.quote('logs/' + job_id + '.openmm.log')}",
                "fi",
                f"if [ -f {shlex.quote('runs/' + job_id + '/trajectory.dcd')} ]; then",
                f"  {r['analysis_command']} 2>&1 | tee {shlex.quote('logs/' + job_id + '.analysis.log')}",
                "fi",
                "",
            ]
        )
    out = DIVERSE / name
    out.write_text("\n".join(lines) + "\n")
    out.chmod(0o755)


def write_scripts(jobs: pd.DataFrame) -> None:
    make_seeded_runner()
    runnable = jobs[jobs["can_launch"]].copy()
    shell_script(
        "run_tier0_short_screens.sh",
        runnable[runnable["launch_group"].eq("tier0_short_screens")],
        "tier0 short template-diversity screens",
    )
    shell_script(
        "run_tier1_gad_full_replicates_controls.sh",
        runnable[runnable["launch_group"].eq("tier1_gad_full_replicates_controls")],
        "tier1 GADGVGKSAL full 10ns replicate/control batch",
    )
    shell_script(
        "run_after_hmtevv_sync.sh",
        runnable[runnable["launch_group"].eq("after_hmtevv_sync")],
        "HMTEVVRHC follow-up batch after current 10ns completion/sync",
    )
    immediate = runnable[runnable["immediate_launch"]].copy()
    shell_script("run_immediate_ready_batch.sh", immediate, "immediate ready diverse MD batch")

    upload = """#!/usr/bin/env bash
set -euo pipefail
POD_HOST=${POD_HOST:-135.84.176.142}
POD_PORT=${POD_PORT:-20878}
KEY=${KEY:-$HOME/.runpod/ssh/RunPod-Key-Go}
REMOTE=${REMOTE:-/runpod-volume/cross_neo_md_diverse_2026_05_10}
SCRIPT=${SCRIPT:-run_immediate_ready_batch.sh}

ssh -i "$KEY" -p "$POD_PORT" -o StrictHostKeyChecking=no root@"$POD_HOST" "mkdir -p '$REMOTE'"
rsync -az --info=progress2 -e "ssh -i $KEY -p $POD_PORT -o StrictHostKeyChecking=no" ./ root@"$POD_HOST":"$REMOTE"/
ssh -i "$KEY" -p "$POD_PORT" -o StrictHostKeyChecking=no root@"$POD_HOST" "cd '$REMOTE' && nohup bash '$SCRIPT' > '${SCRIPT%.sh}.nohup.log' 2>&1 & echo \\$! > '${SCRIPT%.sh}.pid'"
"""
    out = DIVERSE / "runpod_upload_and_start_template.sh"
    out.write_text(upload)
    out.chmod(0o755)

    queue = """#!/usr/bin/env bash
set -euo pipefail
WAIT_PID=${1:?usage: queue_after_pid.sh WAIT_PID [SCRIPT]}
SCRIPT=${2:-run_immediate_ready_batch.sh}
POLL_SECONDS=${POLL_SECONDS:-300}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "[queue] waiting for PID ${WAIT_PID} before ${SCRIPT}"
while kill -0 "$WAIT_PID" 2>/dev/null; do
  date -Is
  sleep "$POLL_SECONDS"
done
echo "[queue] PID ${WAIT_PID} finished; starting ${SCRIPT}"
bash "$SCRIPT"
"""
    out = DIVERSE / "queue_after_pid.sh"
    out.write_text(queue)
    out.chmod(0o755)


def write_reports(jobs: pd.DataFrame, audit: pd.DataFrame) -> None:
    jobs.to_csv(DIVERSE / "diverse_simulation_job_plan.tsv", sep="\t", index=False)
    jobs[jobs["immediate_launch"]].to_csv(DIVERSE / "immediate_launch_jobs.tsv", sep="\t", index=False)
    jobs[~jobs["immediate_launch"]].to_csv(DIVERSE / "hold_or_after_sync_jobs.tsv", sep="\t", index=False)
    audit.to_csv(DIVERSE / "input_pdb_chain_audit.tsv", sep="\t", index=False)

    summary = {
        "total_jobs": int(len(jobs)),
        "immediate_launch_jobs": int(jobs["immediate_launch"].sum()) if not jobs.empty else 0,
        "can_launch_jobs": int(jobs["can_launch"].sum()) if not jobs.empty else 0,
        "jobs_by_launch_group": jobs.groupby("launch_group").size().to_dict() if not jobs.empty else {},
        "jobs_by_tier": jobs.groupby("tier").size().to_dict() if not jobs.empty else {},
        "total_requested_ns_immediate": float(jobs.loc[jobs["immediate_launch"], "ns"].sum()) if not jobs.empty else 0.0,
        "boundary": "MD launch plan only; not immunogenicity, vaccine efficacy, or clinical validation evidence.",
    }
    (DIVERSE / "launch_pack_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    immediate = jobs[jobs["immediate_launch"]].copy()
    hold = jobs[~jobs["immediate_launch"]].copy()
    lines = [
        "# Diverse MD Simulation Launch Pack",
        "",
        "## Boundary",
        "",
        "This is a GPU launch scaffold for structural robustness testing. It does not prove immunogenicity, cancer specificity, vaccine efficacy, or SOTA superiority.",
        "",
        "## Summary",
        "",
        f"- total planned jobs: {summary['total_jobs']}",
        f"- immediate ready jobs: {summary['immediate_launch_jobs']}",
        f"- runnable total jobs: {summary['can_launch_jobs']}",
        f"- immediate requested GPU trajectory length: {summary['total_requested_ns_immediate']:.1f} ns",
        "",
        "## Immediate Jobs",
        "",
        write_markdown_table(
            immediate,
            ["job_id", "peptide", "hla_4digit", "control_type", "complex_scope", "ns", "random_seed", "chain_quality_flag"],
            30,
        ),
        "",
        "## Hold / After-Sync Jobs",
        "",
        write_markdown_table(
            hold,
            ["job_id", "peptide", "hla_4digit", "control_type", "complex_scope", "ns", "hold_reason", "chain_quality_flag"],
            30,
        ),
        "",
        "## Blocked Work",
        "",
        "WT and anchor-preserved decoy simulations remain structure/sequence-curation work. See `blocked_wt_decoy_structure_work.tsv`.",
        "",
        "## Run Commands",
        "",
        "On a CUDA/OpenMM RunPod after upload:",
        "",
        "```bash",
        "bash run_immediate_ready_batch.sh",
        "bash run_after_hmtevv_sync.sh",
        "```",
        "",
        "For upload/start, edit environment variables if the pod endpoint changes, then run:",
        "",
        "```bash",
        "bash runpod_upload_and_start_template.sh",
        "```",
    ]
    (DIVERSE / "README.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    jobs, audit = build_jobs()
    write_scripts(jobs)
    write_reports(jobs, audit)
    print(f"[diverse-md-pack] jobs={len(jobs)} immediate={int(jobs['immediate_launch'].sum()) if not jobs.empty else 0}")
    print(DIVERSE)


if __name__ == "__main__":
    main()
