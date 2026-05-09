#!/usr/bin/env python3
"""Write reviewer-safe decision addendum for RunPod ESM2/QK gate."""

from __future__ import annotations

import pandas as pd

from common import OUT, ensure_dirs


def main() -> None:
    ensure_dirs()
    gate = pd.read_csv(OUT / "metrics/fast_esm2_qk_gate_primary_comparison.tsv", sep="\t")
    esm = pd.read_csv(OUT / "metrics/runpod_esm2_primary_comparison.tsv", sep="\t")
    selected = pd.read_csv(OUT / "metrics/fast_esm2_qk_gate_selected.tsv", sep="\t")
    exact = gate[gate["split_name"].eq("exact_peptide_hla_holdout")].iloc[0]
    near = gate[gate["split_name"].eq("near_peptide_cluster_holdout")].iloc[0]
    verdict = "HOLD_FOR_SOTA__USE_AS_EXACT_SPLIT_RESCUE_EXPERT"
    lines = [
        "# CROSS-Neo v2 RunPod ESM2/QK Gate Decision Addendum",
        "",
        f"Decision: **{verdict}**",
        "",
        "## What Changed",
        "",
        "- RunPod A6000 generated frozen ESM2-35M, ESM2-150M, and ESM2-650M features.",
        "- ESM2 standalone did not beat the v1/v2 locked leaders by AUPRC.",
        "- A fold-safe ESM2+QK gate was then tested using same-split OOF predictions from outer-train rows only.",
        "",
        "## Main Result",
        "",
        gate.to_markdown(index=False),
        "",
        "## Interpretation",
        "",
        f"- Exact peptide-HLA improved from AUPRC {exact['pre_gate_AUPRC']:.3f} to {exact['gate_AUPRC']:.3f}; this is a real locked-split signal.",
        f"- Near-peptide holdout fell from AUPRC {near['pre_gate_AUPRC']:.3f} to {near['gate_AUPRC']:.3f}; this blocks promotion to headline model.",
        "- HLA group gained top10 but lost AUPRC; HLA supertype also lost AUPRC.",
        "- Therefore the gate is a bounded rescue expert, not the main model.",
        "",
        "## Selected Experts",
        "",
        selected.to_markdown(index=False),
        "",
        "## Claim Boundary",
        "",
        "- Allowed: ESM2/QK gate improves exact peptide-HLA locked AUPRC in an internal split-safe test.",
        "- Allowed: ESM2 contributes source/top-k rescue diagnostics, especially TESLA_mmc4/mc7 descriptive recovery.",
        "- Forbidden: SOTA predictor claim, external validation claim, quantum advantage claim, or public-comparator clean claim.",
        "",
        "## Next Move",
        "",
        "Do not keep adding arbitrary models. The next reviewer-relevant step is a stability-filtered gate that abstains from ESM2 fusion on near-peptide/HLA splits unless inner-fold stability criteria are met, plus completion of public training-corpus overlap.",
    ]
    (OUT / "CROSS_Neo_v2_RunPod_ESM2_QK_gate_decision_addendum.md").write_text("\n".join(lines) + "\n")
    print(f"[gate-addendum] {verdict}")


if __name__ == "__main__":
    main()
