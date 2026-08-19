#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd

from neoimmune_common import (
    EXTERNAL_MODELS,
    HUB_ROOT,
    LOCAL_MODELS,
    REPO_ROOT,
    ensure_run_dir,
    find_files,
    rel,
    write_md,
    write_tsv,
)


KEYWORDS = [
    "neoantigen_master.csv",
    "esm2_embeddings.npy",
    "active_learning_top50.tsv",
    "kras_g12_predictions.tsv",
    "leakage_audit.json",
    "all_analyses.json",
    "qx_analyses.json",
    "kkpp_analyses.json",
    "esm2_full_analyses.json",
    "cross_neo_v1",
    "cross_neo_v2",
    "clean_neobench",
    "barneo",
    "strict_class",
    "quantum_kernel",
    "Structure_LR",
    "Wave8",
    "BigMHC",
    "MHCflurry",
    "PRIME",
    "TESLA",
    "IMPROVE",
    "CEDAR",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)
    roots = [REPO_ROOT / "project", REPO_ROOT / "scripts", HUB_ROOT]
    rows = []
    for root in roots:
        if not root.exists():
            continue
        for p in root.rglob("*"):
            if not p.is_file():
                continue
            name = str(p)
            hit = [k for k in KEYWORDS if k.lower() in name.lower()]
            if hit:
                rows.append(
                    {
                        "path": rel(p),
                        "size_mb": round(p.stat().st_size / 1024 / 1024, 3),
                        "matched_keywords": ";".join(hit),
                    }
                )
    files = pd.DataFrame(rows).sort_values(["matched_keywords", "path"]) if rows else pd.DataFrame()
    write_tsv(files, outdir / "audit" / "repo_audit_files.tsv")

    dep_rows = []
    for exe in ["netMHCpan", "netmhcpan", "mhcflurry-predict", "python", "Rscript", "pvacseq", "pvacbind"]:
        dep_rows.append({"dependency": exe, "available_on_path": bool(shutil.which(exe)), "path": shutil.which(exe) or ""})
    deps = pd.DataFrame(dep_rows)
    write_tsv(deps, outdir / "audit" / "dependency_check.tsv")

    reusable = files[files["path"].str.contains("cross_neo|clean_neobench|barneo|wave11|neoantigen_master|leakage", case=False, na=False)]
    md = [
        "# NeoImmune-Stack / CLEAN-Neo++ repo audit",
        "",
        "## Decision",
        "Build an integration layer, not a new foundation model. Existing local algorithms are reusable and should be separated from frozen public predictors.",
        "",
        "## Existing local algorithms to integrate",
    ]
    for m in LOCAL_MODELS:
        present = files["path"].str.contains(m.replace(".", r"\."), case=False, na=False).any() if not files.empty else False
        md.append(f"- {m}: {'found or referenced' if present else 'registry-required; direct artifact not yet found'}")
    md += [
        "",
        "## Existing public/external comparator artifacts",
    ]
    for m in EXTERNAL_MODELS:
        token = m.split("_")[0].replace("2.0", "").replace("4.1", "")
        present = files["path"].str.contains(token, case=False, na=False).any() if not files.empty else False
        md.append(f"- {m}: {'existing output or raw source found' if present else 'pending adapter/status stub'}")
    md += [
        "",
        "## High-value reusable artifacts",
    ]
    for p in reusable["path"].head(80).tolist() if not reusable.empty else []:
        md.append(f"- `{p}`")
    md += [
        "",
        "## Missing dependencies / execution risk",
    ]
    for _, r in deps.iterrows():
        md.append(f"- {r['dependency']}: {'available' if r['available_on_path'] else 'not on PATH'}")
    md += [
        "",
        "## High-risk assumptions",
        "- Public predictor scores may overlap with public training data; they are barred from the clean science track.",
        "- Presentation labels and immunogenicity labels must remain separated.",
        "- Patient-level top-N metrics are only meaningful when patient IDs and positive labels exist.",
        "- Existing strong random-split signals are not accepted as mechanistic or deployable evidence without source/patient/peptide leakage controls.",
        "- LLM usage is restricted to narrative audit, rationale summarization, and report drafting; it is not a vaccine efficacy predictor.",
        "",
        "## Integration map",
        "- Canonical candidates: `clean_neobench_master.tsv`, `cross_neo_v0/master_table.tsv`, `/data/neoantigen_vaccine_hub/data_processed/neoantigen_master.csv`.",
        "- Local scores: `clean_neobench_method_scores.tsv`, BAR-Neo score files, cross-neo score files, wave11 local model predictions.",
        "- External scores: wave11 predictions for BigMHC/MHCflurry/PRIME/NetMHCpan plus adapter stubs for unavailable tools.",
        "- Leakage evidence: clean_neobench overlap flags, cross_neo retrieval audit, neoantigen hub leakage JSON.",
    ]
    write_md("\n".join(md) + "\n", outdir / "audit" / "repo_audit.md")


if __name__ == "__main__":
    main()

