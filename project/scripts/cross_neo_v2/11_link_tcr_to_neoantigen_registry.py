#!/usr/bin/env python3
"""Link the CROSS-Neo-TCR registry to the main CROSS-Neo pMHC registry.

The output is an evidence annotation layer. It does not transfer TCR labels to
the main neoantigen task unless a downstream script explicitly verifies label
compatibility and leakage boundaries.
"""

from __future__ import annotations

import re
from collections import defaultdict

import numpy as np
import pandas as pd

from common import OUT, hla_supertype, levenshtein, peptide_cluster, safe_parquet


TCR_OUT = OUT / "tcr_extension"
RECURRENT_GENES = {"KRAS", "NRAS", "TP53", "BRAF"}


def write_tsv(df: pd.DataFrame, path) -> None:
    view = df.copy()
    view = view.where(view.notna(), "NA").replace("", "NA")
    view.to_csv(path, sep="\t", index=False, na_rep="NA")


def clean(value: object) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    text = str(value).strip()
    return "" if text.lower() in {"", "nan", "none", "na", "<na>"} else text


def norm_peptide(value: object) -> str:
    text = clean(value).upper().replace(" ", "")
    return text if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", text or "") else ""


def get_gene_text(row: pd.Series) -> str:
    parts = []
    for col in ["source_protein", "antigen_source", "mutation_key", "source_window", "source_window_15aa", "source_window_30aa"]:
        if col in row:
            parts.append(clean(row.get(col, "")))
    return " ".join(parts).upper()


def row_recurrent_genes(row: pd.Series) -> set[str]:
    text = get_gene_text(row)
    return {g for g in RECURRENT_GENES if g in text}


def summarize_hits(hits: pd.DataFrame) -> dict[str, object]:
    if hits.empty:
        return {
            "tcr_evidence_count": 0,
            "tcr_evidence_sources": "",
            "paired_tcr_evidence_count": 0,
            "beta_only_tcr_evidence_count": 0,
            "peptide_hla_tcr_label_count": 0,
            "cancer_context_evidence_count": 0,
            "pathogen_context_evidence_count": 0,
            "tcr_binding_positive_count": 0,
            "tcr_binding_negative_count": 0,
            "structure_evidence_count": 0,
            "example_tcr_registry_rows": "",
            "example_tcr_cdr3b": "",
            "example_tcr_cdr3a": "",
            "example_tcr_antigen_sources": "",
        }
    labels = pd.to_numeric(hits.get("binding_label_binary", pd.Series([], dtype=float)), errors="coerce")
    def examples(col: str, n: int = 5) -> str:
        if col not in hits:
            return ""
        vals = hits[col].fillna("").astype(str)
        vals = vals[vals.ne("")]
        return "|".join(vals.drop_duplicates().head(n))

    return {
        "tcr_evidence_count": int(len(hits)),
        "tcr_evidence_sources": "|".join(hits["source_dataset"].dropna().astype(str).drop_duplicates().head(8)),
        "paired_tcr_evidence_count": int(hits.get("paired_tcr_available", False).fillna(False).astype(bool).sum()),
        "beta_only_tcr_evidence_count": int(hits.get("beta_only_tcr_available", False).fillna(False).astype(bool).sum()),
        "peptide_hla_tcr_label_count": int((hits.get("peptide_hla_available", False).fillna(False).astype(bool) & labels.notna()).sum()),
        "cancer_context_evidence_count": int(hits.get("is_cancer_context", False).fillna(False).astype(bool).sum()),
        "pathogen_context_evidence_count": int(hits.get("is_pathogen_context", False).fillna(False).astype(bool).sum()),
        "tcr_binding_positive_count": int((labels == 1).sum()),
        "tcr_binding_negative_count": int((labels == 0).sum()),
        "structure_evidence_count": int(hits.get("structure_pdb_id", "").fillna("").astype(str).ne("").sum()),
        "example_tcr_registry_rows": "|".join(hits["row_id"].astype(str).head(8)),
        "example_tcr_cdr3b": examples("cdr3_beta"),
        "example_tcr_cdr3a": examples("cdr3_alpha"),
        "example_tcr_antigen_sources": examples("antigen_source"),
    }


def main() -> None:
    TCR_OUT.mkdir(parents=True, exist_ok=True)
    neo = pd.read_csv(OUT / "canonical_registry.tsv", sep="\t", low_memory=False)
    tcr = pd.read_csv(TCR_OUT / "tcr_registry.tsv", sep="\t", low_memory=False)

    neo = neo.copy()
    neo["peptide_norm_for_tcr_link"] = neo["peptide"].map(norm_peptide)
    neo["hla_supertype_for_tcr_link"] = neo["hla_4digit"].map(hla_supertype)
    neo["near_peptide_cluster_for_tcr_link"] = neo["peptide_norm_for_tcr_link"].map(peptide_cluster)
    neo["tcr_exact_peptide_hla_key"] = neo["peptide_norm_for_tcr_link"].astype(str) + "|" + neo["hla_4digit"].fillna("").astype(str)
    neo["tcr_near_cluster_supertype_key"] = neo["near_peptide_cluster_for_tcr_link"].astype(str) + "|" + neo["hla_supertype_for_tcr_link"].astype(str)

    tcr = tcr.copy()
    tcr["peptide"] = tcr["peptide"].map(norm_peptide)
    tcr["exact_peptide_hla_key"] = tcr["peptide"].astype(str) + "|" + tcr["hla_4digit"].fillna("").astype(str)
    tcr["peptide_only_key"] = tcr["peptide"].astype(str)
    tcr["near_peptide_cluster_key"] = tcr["peptide"].map(peptide_cluster).astype(str) + "|" + tcr["hla_supertype"].fillna("").astype(str)
    tcr["mutation_gene_norm"] = tcr.get("mutation_gene", "").fillna("").astype(str).str.upper()

    exact_idx: dict[str, list[int]] = defaultdict(list)
    peptide_idx: dict[str, np.ndarray] = defaultdict(list)
    near_idx: dict[str, np.ndarray] = defaultdict(list)
    gene_idx: dict[str, np.ndarray] = defaultdict(list)
    for idx, row in tcr.iterrows():
        pep = clean(row.get("peptide", ""))
        if pep:
            exact_key = clean(row.get("exact_peptide_hla_key", ""))
            if exact_key:
                exact_idx[exact_key].append(idx)
            peptide_idx[pep].append(idx)
            near_idx[clean(row.get("near_peptide_cluster_key", ""))].append(idx)
        gene = clean(row.get("mutation_gene_norm", ""))
        if gene in RECURRENT_GENES:
            gene_idx[gene].append(idx)
        else:
            text = f"{clean(row.get('antigen_source',''))} {clean(row.get('mutation_key',''))}".upper()
            for g in RECURRENT_GENES:
                if g in text:
                    gene_idx[g].append(idx)
    exact_idx = {k: np.asarray(v, dtype=int) for k, v in exact_idx.items()}
    peptide_idx = {k: np.asarray(v, dtype=int) for k, v in peptide_idx.items()}
    near_idx = {k: np.asarray(v, dtype=int) for k, v in near_idx.items()}
    gene_idx = {k: np.asarray(sorted(set(v)), dtype=int) for k, v in gene_idx.items()}

    rows: list[dict[str, object]] = []
    evidence_rows: list[dict[str, object]] = []
    for _, r in neo.iterrows():
        pep = clean(r.get("peptide_norm_for_tcr_link", ""))
        exact_key = clean(r.get("tcr_exact_peptide_hla_key", ""))
        near_key = clean(r.get("tcr_near_cluster_supertype_key", ""))
        genes = row_recurrent_genes(r)

        exact_hits = tcr.loc[exact_idx.get(exact_key, np.asarray([], dtype=int))]
        peptide_hits = tcr.loc[peptide_idx.get(pep, np.asarray([], dtype=int))]
        near_hits = tcr.loc[near_idx.get(near_key, np.asarray([], dtype=int))]
        gene_hits = pd.concat([tcr.loc[gene_idx.get(g, np.asarray([], dtype=int))] for g in genes], ignore_index=False) if genes else pd.DataFrame(columns=tcr.columns)

        if len(exact_hits) and exact_hits.get("paired_tcr_available", False).fillna(False).astype(bool).any():
            category = "exact_tcr_pmhc_match"
            chosen = exact_hits
            evidence_basis = "exact peptide+HLA with at least one paired alpha/beta TCR"
        elif len(exact_hits):
            category = "peptide_hla_match"
            chosen = exact_hits
            evidence_basis = "exact normalized peptide+HLA TCR-resource match; paired TCR may be absent"
        elif len(peptide_hits):
            category = "peptide_only_match"
            chosen = peptide_hits
            evidence_basis = "exact peptide-only match; HLA/source may differ"
        elif len(near_hits):
            category = "near_peptide_match"
            chosen = near_hits
            evidence_basis = "peptide cluster + HLA-supertype match"
        elif len(gene_hits):
            category = "recurrent_cancer_gene_context_match"
            chosen = gene_hits
            evidence_basis = "shared recurrent cancer gene annotation; not peptide-level evidence"
        else:
            category = "no_tcr_match"
            chosen = pd.DataFrame(columns=tcr.columns)
            evidence_basis = "no local TCR-resource evidence"

        rec = r.to_dict()
        rec.update(summarize_hits(chosen))
        rec["tcr_link_category"] = category
        rec["tcr_link_evidence_basis"] = evidence_basis
        rec["tcr_label_transfer_allowed"] = False
        rec["tcr_link_claim_status"] = (
            "diagnostic_only_exact_paired" if category == "exact_tcr_pmhc_match"
            else "diagnostic_only_nondefinitive" if category != "no_tcr_match"
            else "no_tcr_evidence"
        )
        rec["tcr_link_warning"] = (
            "Do not treat peptide-only, near-peptide, pathogen-derived, or recurrent-gene matches as neoantigen ground truth."
            if category not in {"no_tcr_match", "exact_tcr_pmhc_match"}
            else "Exact paired TCR-pMHC evidence still requires source/assay compatibility before model-label use."
            if category == "exact_tcr_pmhc_match"
            else ""
        )
        rows.append(rec)

        for _, hit in chosen.head(20).iterrows():
            evidence_rows.append({
                "neo_row_id": r.get("row_id", ""),
                "neo_peptide": r.get("peptide", ""),
                "neo_hla": r.get("hla_4digit", ""),
                "tcr_row_id": hit.get("row_id", ""),
                "tcr_source_dataset": hit.get("source_dataset", ""),
                "tcr_peptide": hit.get("peptide", ""),
                "tcr_hla": hit.get("hla_4digit", ""),
                "tcr_cdr3_alpha": hit.get("cdr3_alpha", ""),
                "tcr_cdr3_beta": hit.get("cdr3_beta", ""),
                "tcr_paired_available": hit.get("paired_tcr_available", ""),
                "tcr_binding_label_binary": hit.get("binding_label_binary", ""),
                "tcr_cancer_context": hit.get("is_cancer_context", ""),
                "tcr_pathogen_context": hit.get("is_pathogen_context", ""),
                "tcr_antigen_source": hit.get("antigen_source", ""),
                "link_category": category,
                "link_evidence_basis": evidence_basis,
                "peptide_levenshtein": levenshtein(r.get("peptide", ""), hit.get("peptide", ""), max_cutoff=3),
            })

    linked = pd.DataFrame(rows)
    evidence = pd.DataFrame(evidence_rows)
    safe_parquet(linked, TCR_OUT / "tcr_neo_linked_registry.parquet")
    write_tsv(linked, TCR_OUT / "tcr_neo_linked_registry.tsv")
    write_tsv(evidence, TCR_OUT / "tcr_neo_linkage_evidence_examples.tsv")

    summary = (
        linked.groupby("tcr_link_category", dropna=False)
        .agg(
            n_neo_rows=("row_id", "size"),
            mean_tcr_evidence_count=("tcr_evidence_count", "mean"),
            total_tcr_evidence_count=("tcr_evidence_count", "sum"),
            paired_tcr_evidence_count=("paired_tcr_evidence_count", "sum"),
            cancer_context_evidence_count=("cancer_context_evidence_count", "sum"),
            pathogen_context_evidence_count=("pathogen_context_evidence_count", "sum"),
        )
        .reset_index()
        .sort_values("n_neo_rows", ascending=False)
    )
    write_tsv(summary, TCR_OUT / "tcr_neo_linkage_summary.tsv")

    exact = int((linked["tcr_link_category"] == "exact_tcr_pmhc_match").sum())
    pep_hla = int((linked["tcr_link_category"] == "peptide_hla_match").sum())
    pep_only = int((linked["tcr_link_category"] == "peptide_only_match").sum())
    near = int((linked["tcr_link_category"] == "near_peptide_match").sum())
    no = int((linked["tcr_link_category"] == "no_tcr_match").sum())
    recurrent = int((linked["tcr_link_category"] == "recurrent_cancer_gene_context_match").sum())
    lines = [
        "# CROSS-Neo-TCR Linkage Report",
        "",
        "## Claim Boundary",
        "",
        "TCR evidence is added as a diagnostic annotation layer only. The linker does not transfer labels from public TCR resources into the main CROSS-Neo immunogenicity task.",
        "",
        "## Linkage Counts",
        "",
        f"- Main CROSS-Neo rows: {len(linked):,}",
        f"- `exact_tcr_pmhc_match`: {exact:,}",
        f"- `peptide_hla_match`: {pep_hla:,}",
        f"- `peptide_only_match`: {pep_only:,}",
        f"- `near_peptide_match`: {near:,}",
        f"- `recurrent_cancer_gene_context_match`: {recurrent:,}",
        f"- `no_tcr_match`: {no:,}",
        "",
        "## Summary Table",
        "",
        summary.to_markdown(index=False),
        "",
        "## Interpretation Rules",
        "",
        "- Exact peptide+HLA with paired TCR is the strongest local evidence, but still diagnostic until source, assay, and split compatibility are checked.",
        "- Peptide-only matches are not definitive because MHC restriction and antigen source can differ.",
        "- Pathogen-derived matches are not transferred to cancer neoantigens without clear cancer annotation.",
        "- Recurrent-gene matches are context only; they do not imply TCR recognition of the same neoepitope.",
        "",
        "## Outputs",
        "",
        "- `tcr_neo_linked_registry.parquet` / `tcr_neo_linked_registry.tsv`",
        "- `tcr_neo_linkage_summary.tsv`",
        "- `tcr_neo_linkage_evidence_examples.tsv`",
    ]
    (TCR_OUT / "tcr_neo_linkage_report.md").write_text("\n".join(lines) + "\n")
    print(f"[tcr-link] neo_rows={len(linked)} exact_paired={exact} peptide_hla={pep_hla} peptide_only={pep_only} near={near} out={TCR_OUT}")


if __name__ == "__main__":
    main()
