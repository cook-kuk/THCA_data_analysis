#!/usr/bin/env python3
"""Discover and prioritize additional TCR-recognition models for CROSS-Neo-TCR."""

from __future__ import annotations

import importlib.util
import shutil
from pathlib import Path

import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from common import OUT, REPO


TCR_OUT = OUT / "tcr_extension"
DISC_OUT = TCR_OUT / "model_discovery"
NETTCR_LOCAL = Path("/data/thca/_tmp/relocated_2026_05_09/NetTCR-2.2")
MODEL_REPOS = Path("/data/thca/tcr_model_repos")


def module_available(name: str) -> bool:
    return importlib.util.find_spec(name) is not None


def exists(path: Path) -> bool:
    return path.exists()


def local_status(model: str) -> str:
    if model == "NetTCR-2.2":
        if exists(NETTCR_LOCAL / "models/nettcr_2_2_pretrained") and exists(NETTCR_LOCAL / "src/make_webserver_prediction.py"):
            return f"local_ready:{NETTCR_LOCAL}"
    if model == "tFold-TCR":
        return "python_module_available" if module_available("tfold") else "not_found"
    if model == "AlphaFold-Multimer/ColabFold":
        return f"local_ready:{shutil.which('colabfold_batch')}" if shutil.which("colabfold_batch") else "not_found"
    repo_checks = {
        "ERGO-II": (MODEL_REPOS / "ERGO-II", "Predict.py"),
        "pMTnet": (MODEL_REPOS / "pMTnet", "pMTnet.py"),
        "TSpred": (MODEL_REPOS / "TSpred", "example_run2/predict.py"),
        "PanPep": (MODEL_REPOS / "PanPep", "PanPep.py"),
        "TEPCAM": (MODEL_REPOS / "TEPCAM", "ckpts/tepcam_test.pt"),
        "TCR-H": (MODEL_REPOS / "TCR-H", ""),
    }
    if model in repo_checks:
        root, marker = repo_checks[model]
        if root.exists() and (not marker or (root / marker).exists()):
            return f"local_ready:{root}"
    if model == "TCRdock":
        return "local_binary_found" if shutil.which("tcrdock") or shutil.which("TCRdock") else "not_found"
    if model == "Boltz/Boltz-2":
        return "python_module_available" if module_available("boltz") else "not_found"
    if model == "Chai-1":
        return "python_module_available" if module_available("chai_lab") else "not_found"
    return "not_found"


MODELS: list[dict[str, object]] = [
    {
        "model": "NetTCR-2.2",
        "family": "sequence paired-chain CNN",
        "input": "peptide + CDR1/2/3 alpha + CDR1/2/3 beta",
        "hla_modeled": "no explicit HLA",
        "paired_alpha_beta": "required/preferred",
        "code_url": "https://github.com/mnielLab/NetTCR-2.2",
        "source_url": "https://doi.org/10.7554/eLife.93934.3",
        "weights": "local pretrained TFLite weights found",
        "throughput": "high",
        "priority": "P0",
        "claim_use": "TCR-available subset baseline/expert only",
        "blocker": "CROSS-Neo registry lacks CDR1/CDR2; need V-gene CDR reconstruction or full-chain TCR import",
    },
    {
        "model": "ERGO-II",
        "family": "sequence LSTM/attention TCR-peptide",
        "input": "peptide + CDR3 beta; optional alpha, V/J, MHC, T-cell type",
        "hla_modeled": "optional MHC feature",
        "paired_alpha_beta": "not required",
        "code_url": "https://github.com/IdoSpringer/ERGO-II",
        "source_url": "https://github.com/IdoSpringer/ERGO-II",
        "weights": "local pretrained checkpoint and autoencoder weights found",
        "throughput": "high",
        "priority": "P0",
        "claim_use": "strong next practical model because CDR3-only path matches current registry",
        "blocker": "CROSS-Neo input adapter and strict split evaluation",
    },
    {
        "model": "pMTnet",
        "family": "sequence beta-chain peptide-HLA model",
        "input": "CDR3 beta + peptide + HLA",
        "hla_modeled": "yes",
        "paired_alpha_beta": "no, beta-only",
        "code_url": "https://github.com/tianshilu/pMTnet",
        "source_url": "https://github.com/tianshilu/pMTnet",
        "weights": "local trained model/library run successfully",
        "throughput": "high after legacy setup",
        "priority": "P0/P1",
        "claim_use": "diagnostic peptide-HLA-TCR score; useful for beta-only rows",
        "blocker": "strict leakage-filtered benchmark and score calibration",
    },
    {
        "model": "TSpred",
        "family": "paired-chain CNN + reciprocal attention",
        "input": "paired alpha/beta TCR sequence + epitope",
        "hla_modeled": "no explicit HLA in headline task",
        "paired_alpha_beta": "yes",
        "code_url": "https://github.com/ha01994/TSpred",
        "source_url": "https://academic.oup.com/bioinformatics/article/40/8/btae472/7721043",
        "weights": "local example checkpoint run successfully",
        "throughput": "high if pretrained model is available",
        "priority": "P1",
        "claim_use": "paired TCR subset benchmark; seen/unseen epitope split baseline",
        "blocker": "requires CDR1/CDR2 reconstruction for registry-scale use",
    },
    {
        "model": "PanPep",
        "family": "meta-learning peptide-TCR",
        "input": "peptide + TCR sequence, commonly CDR3 beta in public usage; extensions reported for alpha/alpha-beta",
        "hla_modeled": "no explicit HLA in original core task",
        "paired_alpha_beta": "not required in original core task",
        "code_url": "https://github.com/bm2-lab/PanPep",
        "source_url": "https://www.nature.com/articles/s42256-023-00619-3",
        "weights": "local model and memory files run successfully",
        "throughput": "high/medium",
        "priority": "P1",
        "claim_use": "unseen-peptide stress test and negative-control/diagnostic baseline",
        "blocker": "weak zero-shot pilot; try few-shot only if compatible support labels exist",
    },
    {
        "model": "TEPCAM",
        "family": "cross-attention + multi-channel convolution",
        "input": "TCR sequence + peptide",
        "hla_modeled": "no explicit HLA",
        "paired_alpha_beta": "verify code; likely sequence-level TCR/epitope",
        "code_url": "https://github.com/Chenjw99/TEPCAM",
        "source_url": "https://pubmed.ncbi.nlm.nih.gov/37983648/",
        "weights": "local checkpoint run successfully",
        "throughput": "high",
        "priority": "P1/P2",
        "claim_use": "sequence-only TCRbeta-peptide expert; compare to pMTnet",
        "blocker": "no HLA input; needs strict source/epitope-heldout validation",
    },
    {
        "model": "TEINet",
        "family": "deep TCR-epitope sequence model",
        "input": "CDR3 beta + epitope",
        "hla_modeled": "no explicit HLA",
        "paired_alpha_beta": "no",
        "code_url": "https://github.com/jiangdada1221/TEINet",
        "source_url": "https://academic.oup.com/bib/article/doi/10.1093/bib/bbad086/7076118",
        "weights": "check repo",
        "throughput": "high",
        "priority": "P2",
        "claim_use": "legacy/secondary beta-only baseline",
        "blocker": "install and compare against ERGO/pMTnet first",
    },
    {
        "model": "TCR-H",
        "family": "SVM physicochemical CDR3beta-epitope",
        "input": "CDR3 beta + epitope",
        "hla_modeled": "no explicit HLA",
        "paired_alpha_beta": "no",
        "code_url": "https://github.com/rajitha-tatikonda/TCR-H",
        "source_url": "https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1426173/full",
        "weights": "classical model/code",
        "throughput": "very high",
        "priority": "P1/P2",
        "claim_use": "explainable low-compute hard-split baseline",
        "blocker": "clone and map feature script",
    },
    {
        "model": "MixTCRpred",
        "family": "dual-alpha-aware TCR-epitope model",
        "input": "paired/dual alpha TCR context + beta + epitope depending dataset",
        "hla_modeled": "not central",
        "paired_alpha_beta": "yes/optional depending mode",
        "code_url": "https://github.com/viragbioinfo/IMMREP_2022_TCRSpecificity",
        "source_url": "https://www.nature.com/articles/s41467-024-47461-8",
        "weights": "check associated release",
        "throughput": "medium",
        "priority": "P2",
        "claim_use": "dual-alpha case-study model; secondary benchmark",
        "blocker": "repo/tool packaging verification",
    },
    {
        "model": "TCRen",
        "family": "structure/statistical potential",
        "input": "TCR-pMHC model/structure or homology model + peptide context",
        "hla_modeled": "yes via pMHC structure",
        "paired_alpha_beta": "yes for TCR complex",
        "code_url": "https://github.com/antigenomics/tcren-ms",
        "source_url": "https://www.nature.com/articles/s43588-024-00653-0",
        "weights": "code/data; statistical potential",
        "throughput": "medium; structure generation bottleneck",
        "priority": "P1/P2",
        "claim_use": "structure-guided case diagnostics, especially unseen epitopes",
        "blocker": "requires modeled/experimental TCR-pMHC structures",
    },
    {
        "model": "TCRLens",
        "family": "structure-aware EGNN",
        "input": "TCR-pMHC-I structure/interface graph",
        "hla_modeled": "yes, class I focus",
        "paired_alpha_beta": "yes",
        "code_url": "https://github.com/paopitsiri/TCRLens",
        "source_url": "https://www.lifescience.net/publications/1944865/tcrlens-structure-aware-equivariant-graph-learning/",
        "weights": "check repo",
        "throughput": "medium; structure generation bottleneck",
        "priority": "P2",
        "claim_use": "future structure branch benchmark",
        "blocker": "requires actual TCR-pMHC-I structures and installation",
    },
    {
        "model": "tFold-TCR",
        "family": "structure prediction",
        "input": "full TCR alpha/beta + MHC + beta2M + peptide sequences",
        "hla_modeled": "yes via structure",
        "paired_alpha_beta": "yes",
        "code_url": "https://github.com/TencentAI4S/tfold",
        "source_url": "https://github.com/TencentAI4S/tfold",
        "weights": "downloaded via package functions if installed",
        "throughput": "medium/high on GPU",
        "priority": "P1",
        "claim_use": "structure feature generator, not recognition classifier alone",
        "blocker": "not installed; full-chain sequences missing for many jobs",
    },
    {
        "model": "TCRdock",
        "family": "TCR-pMHC structure/docking geometry",
        "input": "peptide-HLA + TCR alpha/beta V/J/CDR3 or structures",
        "hla_modeled": "yes",
        "paired_alpha_beta": "yes",
        "code_url": "https://github.com/phbradley/TCRdock",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9859041/",
        "weights": "AlphaFold/TCRdock workflow assets",
        "throughput": "medium",
        "priority": "P1",
        "claim_use": "geometry/contact diagnostics",
        "blocker": "not installed locally",
    },
    {
        "model": "AlphaFold-Multimer/ColabFold",
        "family": "general complex structure prediction",
        "input": "full chain FASTA: TCR alpha/beta, MHC, beta2M, peptide",
        "hla_modeled": "yes",
        "paired_alpha_beta": "yes for TCR-pMHC",
        "code_url": "https://github.com/google-deepmind/alphafold",
        "source_url": "https://github.com/sokrypton/ColabFold",
        "weights": "local colabfold wrapper available",
        "throughput": "medium",
        "priority": "P1",
        "claim_use": "structure feature generator and parser validation",
        "blocker": "full-chain sequences missing for most registry rows",
    },
    {
        "model": "UniPMT",
        "family": "unified peptide-MHC-TCR graph/multitask",
        "input": "peptide + MHC + TCR features depending task",
        "hla_modeled": "yes",
        "paired_alpha_beta": "verify, often CDR3 beta heavy",
        "code_url": "https://github.com/ethanmock/UniPMT",
        "source_url": "https://www.nature.com/articles/s42256-025-01002-0",
        "weights": "check repo",
        "throughput": "high/medium",
        "priority": "P1/P2",
        "claim_use": "unified pMHC/TCR expert benchmark",
        "blocker": "install and verify checkpoint/data format",
    },
    {
        "model": "UnifyImmun",
        "family": "unified cross-attention pHLA/pTCR",
        "input": "peptide + HLA and/or peptide + TCR",
        "hla_modeled": "yes for pHLA branch",
        "paired_alpha_beta": "verify",
        "code_url": "https://github.com/hliulab/UnifyImmun",
        "source_url": "https://www.nature.com/articles/s42256-024-00973-w",
        "weights": "repo trained model folder reported",
        "throughput": "high/medium",
        "priority": "P2",
        "claim_use": "cross-attention baseline, pHLA+pTCR multitask comparison",
        "blocker": "install and checkpoint verification",
    },
    {
        "model": "ImRex",
        "family": "interaction-map CNN",
        "input": "CDR3 beta + epitope",
        "hla_modeled": "no explicit HLA",
        "paired_alpha_beta": "no",
        "code_url": "https://github.com/pmoris/ImRex",
        "source_url": "https://github.com/pmoris/ImRex",
        "weights": "small pretrained models reported",
        "throughput": "high",
        "priority": "P2",
        "claim_use": "legacy beta-only baseline",
        "blocker": "lower priority than ERGO/pMTnet/PanPep",
    },
    {
        "model": "TITAN",
        "family": "bimodal attention",
        "input": "TCR sequence + epitope",
        "hla_modeled": "no explicit HLA",
        "paired_alpha_beta": "not central",
        "code_url": "https://github.com/PaccMann/TITAN",
        "source_url": "https://arxiv.org/abs/2105.03323",
        "weights": "trained_model directory in repo",
        "throughput": "high/medium",
        "priority": "P2",
        "claim_use": "legacy benchmark",
        "blocker": "not a first next model",
    },
]


def net_tcr_sanity_metrics() -> dict[str, object]:
    path = DISC_OUT / "nettcr_sanity/nettcr_small_example_predictions.csv"
    if not path.exists():
        return {"status": "not_run"}
    df = pd.read_csv(path)
    if {"binder", "prediction"}.issubset(df.columns):
        return {
            "status": "ran",
            "n": int(len(df)),
            "n_pos": int(pd.to_numeric(df["binder"], errors="coerce").fillna(0).sum()),
            "auprc": float(average_precision_score(df["binder"].astype(int), df["prediction"].astype(float))),
            "auroc": float(roc_auc_score(df["binder"].astype(int), df["prediction"].astype(float))),
            "path": str(path),
        }
    return {"status": "output_missing_required_columns", "path": str(path)}


def make_matrix() -> pd.DataFrame:
    rows = []
    for i, rec in enumerate(MODELS, start=1):
        row = rec.copy()
        row["rank"] = i
        row["local_status"] = local_status(str(row["model"]))
        row["runnable_now"] = row["local_status"].astype(str) if hasattr(row["local_status"], "astype") else (
            "yes" if str(row["local_status"]).startswith("local_ready") else "no"
        )
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    DISC_OUT.mkdir(parents=True, exist_ok=True)
    matrix = make_matrix()
    matrix.to_csv(DISC_OUT / "tcr_model_discovery_matrix.tsv", sep="\t", index=False, na_rep="NA")

    sanity = net_tcr_sanity_metrics()
    next_actions = pd.DataFrame(
        [
            {
                "order": 1,
                "action": "Promote pMTnet and TEPCAM into the next strict external-model benchmark.",
                "owner_script_to_create": "15_train_tcr_aware_models.py extension",
                "reason": "Both produced strong same-pilot shuffled-decoy signal without CDR1/CDR2 reconstruction.",
            },
            {
                "order": 2,
                "action": "Build ERGO-II adapter for rows with peptide, MHC, CDR3 alpha/beta, V/J, and T-cell type when available.",
                "owner_script_to_create": "20_run_ergo_cross_neo_adapter.py",
                "reason": "Runtime now works; adapter quality is the next blocker.",
            },
            {
                "order": 3,
                "action": "Add CDR1/CDR2 reconstruction from V genes for NetTCR-2.2/TSpred.",
                "owner_script_to_create": "21_prepare_full_cdr_inputs_for_nettcr.py",
                "reason": "NetTCR-2.2 and TSpred are local and runnable, but need all six CDRs.",
            },
            {
                "order": 4,
                "action": "Keep PanPep as a diagnostic or negative-control unless few-shot mode improves.",
                "owner_script_to_create": "19_compare_external_tcr_models.py",
                "reason": "Zero-shot same-pilot result was near random.",
            },
            {
                "order": 5,
                "action": "Keep tFold/TCRdock/TCRen/TCRLens as structure branch after full-chain sequence recovery.",
                "owner_script_to_create": "13/14 structure scripts extension",
                "reason": "Current structure manifest is blocked by missing full TCR/MHC sequences.",
            },
        ]
    )
    next_actions.to_csv(DISC_OUT / "tcr_model_next_actions.tsv", sep="\t", index=False)

    p0 = matrix[matrix["priority"].astype(str).str.startswith("P0")]
    p1 = matrix[matrix["priority"].astype(str).str.startswith("P1")]
    lines = [
        "# CROSS-Neo-TCR Model Discovery Sprint",
        "",
        "Date: 2026-05-09",
        "",
        "## Bottom Line",
        "",
        "We found enough additional model surface to keep pushing. The immediate model queue should not be structure-first; it should be **pMTnet + TEPCAM first**, then **ERGO-II adapter**, then **NetTCR-2.2/TSpred after CDR1/CDR2 reconstruction**.",
        "",
        "## Local Runtime Discovery",
        "",
        f"- NetTCR-2.2 local repo: `{NETTCR_LOCAL}`",
        f"- NetTCR-2.2 pretrained TFLite weights found: `{(NETTCR_LOCAL / 'models/nettcr_2_2_pretrained').exists()}`",
        f"- ColabFold batch executable: `{shutil.which('colabfold_batch') or 'not_found'}`",
        f"- TensorFlow available: `{module_available('tensorflow')}`",
        f"- PyTorch available: `{module_available('torch')}`",
        f"- Transformers available: `{module_available('transformers')}`",
        f"- ESM available: `{module_available('esm')}`",
        f"- Stitchr/ANARCI/tcrdist3 available: `stitchr={module_available('stitchr')}`, `anarci={module_available('anarci')}`, `tcrdist3={module_available('tcrdist3')}`",
        "",
        "## NetTCR-2.2 Sanity Run",
        "",
    ]
    if sanity.get("status") == "ran":
        lines.extend(
            [
                f"- Small example predictions completed: n={sanity['n']}, positives={sanity['n_pos']}",
                f"- Example AUPRC={sanity['auprc']:.3f}, AUROC={sanity['auroc']:.3f}",
                f"- Output: `{sanity['path']}`",
                "- This is only a runtime sanity check on repository example data, not a CROSS-Neo validation result.",
            ]
        )
    else:
        lines.append(f"- Sanity status: {sanity}")
    comp_path = DISC_OUT / "external_tcr_model_pilot_compare/external_tcr_model_pilot_comparison.tsv"
    if comp_path.exists():
        comp = pd.read_csv(comp_path, sep="\t")
        lines.extend(
            [
                "",
                "## Same-Pilot External Model Comparison",
                "",
                "| Model | Inputs | n | AUPRC | AUROC | Read |",
                "|---|---|---:|---:|---:|---|",
            ]
        )
        for r in comp.itertuples(index=False):
            lines.append(f"| {r.model} | {r.inputs} | {r.n} | {r.auprc:.4f} | {r.auroc:.4f} | {r.claim_boundary} |")
        lines.append("")
        lines.append("These are fast shuffled-decoy pilot numbers, not clean external benchmark claims.")
    bench_path = DISC_OUT / "external_tcr_expert_benchmark/external_tcr_expert_benchmark_metrics.tsv"
    if bench_path.exists():
        bench = pd.read_csv(bench_path, sep="\t")
        pooled = bench[bench["panel"].eq("pooled")].sort_values("auprc", ascending=False)
        lines.extend(
            [
                "",
                "## Fast Multi-Panel Challenge Benchmark",
                "",
                "| Model | n | AUPRC | AUROC |",
                "|---|---:|---:|---:|",
            ]
        )
        for r in pooled.itertuples(index=False):
            lines.append(f"| {r.model} | {r.n} | {r.auprc:.4f} | {r.auroc:.4f} |")
        lines.append("")
        lines.append("The no-training mean pMTnet/TEPCAM ensemble is currently the best practical diagnostic score.")
    lines.extend(
        [
            "",
            "## P0 Models",
            "",
            "| Model | Why P0 | Immediate blocker | Local status |",
            "|---|---|---|---|",
        ]
    )
    for r in p0.itertuples(index=False):
        lines.append(f"| {r.model} | {r.claim_use} | {r.blocker} | {r.local_status} |")
    lines.extend(
        [
            "",
            "## P1 Models",
            "",
            "| Model | Role | Immediate blocker |",
            "|---|---|---|",
        ]
    )
    for r in p1.itertuples(index=False):
        lines.append(f"| {r.model} | {r.claim_use} | {r.blocker} |")
    lines.extend(
        [
            "",
            "## Expanded Model Inventory",
            "",
            "| Model | Family | Input | HLA modeled | Paired alpha/beta | Priority | Local status |",
            "|---|---|---|---|---|---|---|",
        ]
    )
    for r in matrix.itertuples(index=False):
        lines.append(f"| {r.model} | {r.family} | {r.input} | {r.hla_modeled} | {r.paired_alpha_beta} | {r.priority} | {r.local_status} |")
    lines.extend(
        [
            "",
            "## Hard Decision",
            "",
            "1. **Run pMTnet and TEPCAM before structure-first work**, because current rows already have peptide and CDR3 beta, and pMTnet additionally uses HLA.",
            "2. **Move ERGO-II next**, because runtime is now fixed but CROSS-Neo input mapping still needs a clean adapter.",
            "3. **Do not throw away NetTCR-2.2/TSpred**. They are concrete local assets, but need a CDR1/CDR2 reconstruction layer from V genes or full-chain TCR imports.",
            "4. **Do not make structure the next bottleneck**. tFold/TCRdock/TCRen/TCRLens become valuable after full-chain sequence recovery; today they are mostly blocked.",
            "5. **Use all external models as optional TCR expert scores**, not as replacements for the pMHC ranker.",
            "",
            "## Sources",
            "",
            "- NetTCR-2.2: https://github.com/mnielLab/NetTCR-2.2",
            "- TSpred: https://academic.oup.com/bioinformatics/article/40/8/btae472/7721043",
            "- PanPep: https://www.nature.com/articles/s42256-023-00619-3",
            "- TEPCAM: https://pubmed.ncbi.nlm.nih.gov/37983648/",
            "- TCR-H: https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1426173/full",
            "- TCRen: https://www.nature.com/articles/s43588-024-00653-0",
            "- TCRLens: https://www.lifescience.net/publications/1944865/tcrlens-structure-aware-equivariant-graph-learning/",
            "- tFold-TCR: https://github.com/TencentAI4S/tfold",
            "- Comparative model dependency warning: https://pmc.ncbi.nlm.nih.gov/articles/PMC10152969/",
        ]
    )
    report = "\n".join(lines) + "\n"
    (DISC_OUT / "01_tcr_model_discovery.md").write_text(report)
    print(f"[tcr-model-discovery] models={len(matrix)} p0={len(p0)} p1={len(p1)} sanity={sanity.get('status')}")


if __name__ == "__main__":
    main()
