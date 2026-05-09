#!/usr/bin/env python3
"""Build public-overlap audit manifest and run local audits where possible."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from cross_neo_v1_lockdown_common import OUT, REPO, ensure_dirs, load_master, peptide_cluster


SOURCES = {
    "MHCflurry": {
        "search": "MHCflurry training data presentation affinity IEDB MS ligands",
        "columns": "peptide,HLA/allele,measurement/source",
        "target": "project/data/external/mhcflurry_training/",
    },
    "NetMHCpan": {
        "search": "NetMHCpan 4.1 training data BA EL peptide HLA",
        "columns": "peptide,HLA,BA/EL label",
        "target": "project/data/external/netmhcpan_training/",
    },
    "BigMHC": {
        "search": "BigMHC Mendeley el_train im_train peptide HLA immunogenicity",
        "columns": "peptide,HLA,label,split",
        "target": "project/data/external/bigmhc/",
    },
    "PRIME": {
        "search": "PRIME neoantigen immunogenicity training data peptide HLA",
        "columns": "peptide,HLA,immunogenicity",
        "target": "project/data/external/prime/",
    },
    "MixMHCpred": {
        "search": "MixMHCpred training ligands peptide HLA allele",
        "columns": "peptide,HLA/source allele",
        "target": "project/data/external/mixmhcpred/",
    },
    "NetMHCstabpan": {
        "search": "NetMHCstabpan training data peptide HLA stability",
        "columns": "peptide,HLA,stability",
        "target": "project/data/external/netmhcstabpan_training/",
    },
    "IEDB": {
        "search": "IEDB T cell assay export peptide MHC class I immunogenicity",
        "columns": "peptide,HLA,assay,result,reference",
        "target": "project/data/external/iedb/",
    },
    "CEDAR": {
        "search": "local CROSS-Neo master table CEDAR",
        "columns": "peptide_mut,hla,label,study",
        "target": "project/results/cross_neo_v0/master_table.tsv",
    },
    "NEPdb": {
        "search": "local CROSS-Neo master table NEPdb",
        "columns": "peptide_mut,hla,label,study",
        "target": "project/results/cross_neo_v0/master_table.tsv",
    },
    "TESLA": {
        "search": "local CROSS-Neo master table TESLA mmc4 mmc7",
        "columns": "peptide_mut,hla,label,study",
        "target": "project/results/cross_neo_v0/master_table.tsv",
    },
}


def find_local(source: str, target: str) -> list[str]:
    p = REPO / target
    if p.exists():
        return [str(p.relative_to(REPO))]
    patterns = [source.lower(), source.replace("MHC", "mhc").lower()]
    hits = []
    roots = [REPO / "project/data", REPO / "project/results", REPO / "external"]
    for root in roots:
        if not root.exists():
            continue
        for f in root.rglob("*"):
            if not f.is_file():
                continue
            name = f.name.lower()
            if any(pat in name for pat in patterns):
                hits.append(str(f.relative_to(REPO)))
                if len(hits) >= 6:
                    return hits
    return hits


def local_source_rows(master: pd.DataFrame, source: str) -> pd.DataFrame:
    if source == "TESLA":
        return master[master["study"].astype(str).str.contains("TESLA", na=False)].copy()
    if source in {"CEDAR", "NEPdb"}:
        return master[master["study"].eq(source)].copy()
    if source == "IEDB":
        return master[master["study"].isin(["CEDAR", "NEPdb"])].copy()
    return pd.DataFrame()


def main() -> None:
    ensure_dirs()
    master = load_master()
    strict = master[master["strict_set_flag"].astype(bool)].copy()
    strict_pep = set(strict["peptide_mut"].dropna().astype(str))
    strict_pair = set(zip(strict["peptide_mut"].astype(str), strict["hla"].astype(str)))
    strict_cluster = set(strict["peptide_mut"].map(peptide_cluster))

    manifest_rows = []
    audit_rows = []
    for source, meta in SOURCES.items():
        hits = find_local(source, meta["target"])
        local_df = local_source_rows(master, source)
        possible = len(local_df) > 0
        status = "local_audit_available" if possible else ("local_file_found_unparsed" if hits else "missing_unresolved")
        manifest_rows.append(
            {
                "source": source,
                "status": status,
                "local_file_found_or_missing": ";".join(hits) if hits else "missing",
                "expected_url_or_search_term": meta["search"],
                "expected_columns": meta["columns"],
                "local_target_path": meta["target"],
                "exact_peptide_audit_possible": bool(possible),
                "peptide_hla_audit_possible": bool(possible and "hla" in local_df.columns),
                "near_peptide_audit_possible": bool(possible),
                "source_protein_audit_possible": bool(possible and local_df.get("source_protein", pd.Series(dtype=object)).notna().any()),
                "action_needed": "download/parse official training corpus" if not possible else "row-level audit completed for local rows",
            }
        )
        if possible:
            local_pep = set(local_df["peptide_mut"].dropna().astype(str))
            local_pair = set(zip(local_df["peptide_mut"].astype(str), local_df["hla"].astype(str)))
            local_cluster = set(local_df["peptide_mut"].map(peptide_cluster))
            source_protein_overlap = 0
            if "source_protein" in local_df.columns and "source_protein" in strict.columns:
                sp = set(local_df["source_protein"].dropna().astype(str))
                source_protein_overlap = len(sp & set(strict["source_protein"].dropna().astype(str)))
            audit_rows.append(
                {
                    "source": source,
                    "local_n": len(local_df),
                    "strict_n": len(strict),
                    "exact_peptide_overlap_n": len(strict_pep & local_pep),
                    "exact_peptide_hla_overlap_n": len(strict_pair & local_pair),
                    "near_peptide_cluster_overlap_n": len(strict_cluster & local_cluster),
                    "hla_normalized_overlap_n": len(set(strict["hla"].dropna().astype(str)) & set(local_df["hla"].dropna().astype(str))),
                    "source_protein_overlap_n": source_protein_overlap,
                    "audit_scope": "local master rows only; public comparator training corpus still separate unless source is local",
                }
            )

    manifest = pd.DataFrame(manifest_rows)
    audit = pd.DataFrame(audit_rows)
    manifest.to_csv(OUT / "public_overlap_manifest.tsv", sep="\t", index=False)
    audit.to_csv(OUT / "public_overlap_audit_if_available.tsv", sep="\t", index=False)

    unresolved = manifest[manifest["status"].ne("local_audit_available")]
    lines = [
        "# CROSS-Neo v1 Public Overlap Manifest",
        "",
        "This manifest turns unresolved public-pretrained comparator overlap into concrete next actions. Public predictor scores remain forbidden as model features.",
        "",
        "## Manifest",
        "",
        manifest.to_markdown(index=False),
        "",
        "## Local Audit If Available",
        "",
        audit.to_markdown(index=False) if len(audit) else "No local audit rows.",
        "",
        "## Unresolved Downloads",
        "",
        unresolved[["source", "expected_url_or_search_term", "local_target_path", "action_needed"]].to_markdown(index=False) if len(unresolved) else "All sources locally auditable.",
    ]
    (OUT / "public_overlap_report.md").write_text("\n".join(lines) + "\n")
    print(f"[overlap-manifest] manifest={len(manifest)} audit={len(audit)} unresolved={len(unresolved)}")


if __name__ == "__main__":
    main()
