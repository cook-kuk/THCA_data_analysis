#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from neoimmune_common import ensure_run_dir, write_tsv


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)
    rows = [
        ["Structure_LR", "local", "structure_proxy/logistic", "clean baseline", True, "internal artifact", "peptide,HLA,structure proxies", "immunogenicity score", True, True],
        ["Wave8_TCR_SelfSim_full", "local", "TCR/self-similarity", "immunogenicity branch", True, "internal artifact", "peptide,TCR/self-similarity features", "TCR-visible score", True, True],
        ["ESM2_Bayesian", "local", "protein language embedding Bayesian", "sequence representation branch", True, "local embeddings/artifact", "peptide/HLA embeddings", "immunogenicity score", True, True],
        ["GP_quantum", "local", "Gaussian process/quantum features", "experimental branch", True, "internal artifact", "quantum kernel features", "immunogenicity score", True, True],
        ["VQC", "local", "variational quantum classifier", "experimental branch", False, "internal/pending rerun", "quantum features", "immunogenicity score", True, True],
        ["W7A_QK_only", "local", "quantum-kernel-only", "clean algorithm branch", True, "internal artifact", "quantum kernel matrix/features", "immunogenicity score", True, True],
        ["W7A_full", "local", "multi-feature clean model", "clean algorithm branch", True, "internal artifact", "clean peptide/HLA/structure/TCR features", "immunogenicity score", True, True],
        ["W7B_stacked", "local", "stacked local ensemble", "clean stack if trained leakage-safe", True, "internal artifact", "local scores only", "stack score", True, True],
        ["quantum_kernel_no_anchor_gamma1.0", "local", "quantum kernel", "must survive source-heldout", True, "internal artifact", "peptide/HLA quantum kernel", "kernel score", True, True],
        ["quantum_kernel_compact_gamma0.5", "local", "quantum kernel", "compact comparator", True, "internal artifact", "peptide/HLA quantum kernel", "kernel score", True, True],
        ["Stack_mean", "local", "score stack", "simple clean ensemble", True, "internal artifact", "local scores only", "mean score", True, True],
        ["Stack_median", "local", "score stack", "robust clean ensemble", True, "internal artifact", "local scores only", "median score", True, True],
        ["Stack_LR", "local", "logistic stack", "clean stack if leakage-safe", True, "internal artifact", "local scores only", "calibrated stack", True, True],
        ["BAR-Neo", "local", "patient-gated score", "business/practical comparator", True, "internal artifact", "candidate context/local scores", "candidate score", False, True],
        ["BAR-Neo-X", "local", "reviewer-aware ensemble", "production support", True, "internal artifact", "candidate context/local+support", "candidate score", False, True],
        ["KG_GA_evolved_controller", "local_architecture_search", "knowledge-graph genetic/evolutionary controller", "production/experiment-priority branch", True, "internal retrospective GA artifact; includes public predictor components", "candidate_id + local/public/component scores", "kg_ga_evolved_score", False, True],
        ["Darwin_RL_blueprint", "local_architecture_search", "RL/GA method-mining blueprint", "future search policy blueprint", True, "internal blueprint artifact", "method graph/action space", "architecture policy", False, False],
        ["NetMHCpan_4.1_EL", "external", "HLA ligand elution predictor", "presentation comparator/frozen feature", False, "license/manual install", "peptide,HLA", "EL rank/score", False, True],
        ["NetMHCpan_4.1_BA", "external", "binding affinity predictor", "binding comparator/frozen feature", False, "license/manual install", "peptide,HLA", "BA rank/nM", False, True],
        ["MHCflurry_2.0_presentation", "external", "presentation predictor", "presentation comparator/frozen feature", True, "public package/artifact", "peptide,HLA", "presentation score", False, True],
        ["MHCflurry_2.0_affinity", "external", "binding predictor", "binding comparator/frozen feature", True, "public package/artifact", "peptide,HLA", "affinity score", False, True],
        ["BigMHC_EL", "external", "EL predictor", "presentation comparator/frozen feature", False, "external model/artifact pending", "peptide,HLA", "EL score", False, True],
        ["BigMHC_IM", "external", "immunogenicity predictor", "public comparator/frozen feature", True, "public model/artifact", "peptide,HLA", "IM score", False, True],
        ["PRIME", "external", "immunogenicity/presentation predictor", "public comparator/frozen feature", True, "license sensitive/artifact", "peptide,HLA", "PRIME score", False, True],
        ["PRIME2.1", "external", "updated PRIME family", "pending comparator", False, "license/install check required", "peptide,HLA", "PRIME2.1 score", False, True],
        ["HLApollo", "external", "HLA/presentation model", "pending comparator", False, "runnability unknown", "peptide,HLA", "presentation score", False, True],
        ["pVACtools_parser", "external_parser", "pipeline output parser", "ingest external clinical pipeline output", True, "parser only", "pVACtools TSV", "frozen features", False, True],
        ["NeoDisc_parser", "external_parser", "pipeline output parser", "ingest external pipeline output", False, "parser stub", "NeoDisc output", "frozen features", False, True],
        ["LLM_Rationale_Auditor", "analysis_infra", "latest LLM optional", "report/rationale QA only", False, "provider API/env needed", "ranked candidates + evidence", "narrative audit", False, False],
    ]
    df = pd.DataFrame(
        rows,
        columns=[
            "model_name",
            "source",
            "type",
            "role",
            "runnable",
            "license_note",
            "input_required",
            "output_score",
            "clean_track_allowed",
            "production_track_allowed",
        ],
    )
    write_tsv(df, outdir / "registry" / "model_registry.tsv")


if __name__ == "__main__":
    main()
