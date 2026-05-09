#!/usr/bin/env python3
"""Build a canonical TCR-pMHC registry for the CROSS-Neo-TCR extension.

The registry is deliberately permissive: rows with missing peptide, HLA, or TCR
chains are retained and flagged rather than dropped. This keeps the missingness
reality visible for downstream claim-boundary decisions.
"""

from __future__ import annotations

import argparse
import ast
import json
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from common import OUT, REPO, hla_supertype, peptide_cluster, safe_parquet


TCR_OUT = OUT / "tcr_extension"
AA20 = set("ACDEFGHIKLMNPQRSTVWY")
AA_EXT = set("ACDEFGHIKLMNPQRSTVWYX")
MISSING = {"", "NA", "NAN", "NONE", "NULL", "N/A", "NOT_AVAILABLE", "NOT AVAILABLE", "UNSPECIFIED"}
KEYWORDS = (
    "vdjdb", "mcpas", "iedb", "pird", "10x", "tettcr", "tet-tcr", "stcrdab",
    "tcr3d", "tcr", "cdr3", "tra", "trb", "epitope", "peptide", "hla", "mhc",
)
RECURRENT_GENES = {"KRAS", "NRAS", "TP53", "BRAF"}

CANONICAL_COLS = [
    "row_id",
    "source_dataset",
    "study_id",
    "assay_type",
    "species",
    "disease_context",
    "cancer_type",
    "antigen_source",
    "peptide",
    "peptide_length",
    "hla_raw",
    "hla_4digit",
    "mhc_class",
    "mhc_gene",
    "tcr_alpha_v",
    "tcr_alpha_j",
    "cdr3_alpha",
    "tcr_beta_v",
    "tcr_beta_j",
    "cdr3_beta",
    "paired_tcr_available",
    "clone_count",
    "clone_frequency",
    "binding_label_raw",
    "binding_label_binary",
    "activation_label",
    "tetramer_label",
    "cytokine_label",
    "structure_pdb_id",
    "train_test_original",
    "public_overlap_flags",
]

EXTRA_COLS = [
    "source_record_id",
    "source_path",
    "source_row_count",
    "peptide_raw",
    "hla_supertype",
    "near_peptide_cluster",
    "exact_peptide_hla_key",
    "peptide_only_key",
    "near_peptide_cluster_key",
    "beta_only_tcr_available",
    "alpha_only_tcr_available",
    "any_tcr_available",
    "peptide_hla_available",
    "mutation_gene",
    "mutation_key",
    "is_cancer_context",
    "is_pathogen_context",
    "parser_notes",
]


def clean_text(value: object) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    text = str(value).strip()
    if text.upper() in MISSING:
        return ""
    return text


def first_nonempty(*values: object) -> str:
    for value in values:
        text = clean_text(value)
        if text:
            return text
    return ""


def normalize_peptide(value: object) -> str:
    text = clean_text(value).upper().replace(" ", "")
    if not text:
        return ""
    # Preserve non-peptide rows as missing rather than manufacturing false AA strings.
    if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", text):
        return ""
    return text


def normalize_cdr3(value: object) -> str:
    text = clean_text(value).upper().replace(" ", "")
    if not text:
        return ""
    if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWYX]+", text):
        return ""
    return text


def normalize_gene(value: object) -> str:
    text = clean_text(value).upper().replace(" ", "")
    return text


def normalize_hla_any(value: object) -> str:
    raw = clean_text(value)
    if not raw:
        return ""
    h = raw.upper().replace(" ", "").replace("_", "*")
    h = re.sub(r"^HLA-", "", h)
    h = h.replace("HLA", "")
    if h in {"B2M", "BETA2MICROGLOBULIN"}:
        return ""
    if h.startswith(("H-2", "H2-", "MAMU", "PATR", "RT1", "HUMANCD1", "CD1")):
        return raw.upper().replace(" ", "")

    m = re.match(r"^([ABC])\*?0?(\d{1,2})(?::?(\d{2}))?(?::\d+)?$", h)
    if m:
        gene, fam, allele = m.group(1), int(m.group(2)), m.group(3) or "01"
        return f"HLA-{gene}*{fam:02d}:{allele}"

    m = re.match(r"^(DRB[1-5]|DQA1|DQB1|DPA1|DPB1)\*?(\d{2})(?::?(\d{2}))?(?::\d+)?$", h)
    if m:
        gene, fam, allele = m.group(1), m.group(2), m.group(3) or "01"
        return f"HLA-{gene}*{fam}:{allele}"

    m = re.match(r"^DR(\d{1,2})$", h)
    if m:
        return f"HLA-DRB1*{int(m.group(1)):02d}:01"

    if re.match(r"^(DQ|DP)\d+$", h):
        return f"HLA-{h}"
    return raw.upper() if raw.upper().startswith("HLA-") else f"HLA-{raw.upper()}"


def mhc_gene_from_hla(hla: object) -> str:
    h = normalize_hla_any(hla)
    if not h.startswith("HLA-"):
        return ""
    body = h.replace("HLA-", "", 1)
    if "*" in body:
        return body.split("*", 1)[0]
    m = re.match(r"([A-Z]+)", body)
    return m.group(1) if m else ""


def normalize_mhc_class(value: object, hla: object = "") -> str:
    text = clean_text(value).upper().replace(" ", "")
    if text in {"MHCI", "MHC-I", "I", "CLASSI", "CLASSI."}:
        return "I"
    if text in {"MHCII", "MHC-II", "II", "CLASSII", "CLASSII."}:
        return "II"
    gene = mhc_gene_from_hla(hla)
    if gene in {"A", "B", "C"}:
        return "I"
    if gene.startswith(("DR", "DQ", "DP")):
        return "II"
    return ""


def binary_label(value: object) -> object:
    text = clean_text(value).upper()
    if not text:
        return pd.NA
    if text in {"1", "TRUE", "T", "P", "POS", "POSITIVE", "POSITIVE-HIGH", "POSITIVE-LOW"}:
        return 1
    if text in {"0", "FALSE", "F", "N", "NEG", "NEGATIVE"}:
        return 0
    if "POSITIVE" in text and "NEGATIVE" not in text:
        return 1
    if "NEGATIVE" in text:
        return 0
    return pd.NA


def parse_jsonish(value: object) -> dict[str, object]:
    text = clean_text(value)
    if not text:
        return {}
    try:
        obj = json.loads(text)
        return obj if isinstance(obj, dict) else {}
    except Exception:
        try:
            obj = ast.literal_eval(text)
            return obj if isinstance(obj, dict) else {}
        except Exception:
            return {}


def context_flags(*values: object) -> tuple[bool, bool]:
    text = " ".join(clean_text(v) for v in values).lower()
    cancer = bool(re.search(r"cancer|tumou?r|melanoma|carcinoma|leukemia|lymphoma|myeloma|glioma|sarcoma|neoantigen", text))
    pathogen = bool(re.search(r"hiv|cmv|ebv|influenza|sars|cov|virus|viral|bacteria|tuberculosis|parasite|plasmodium|hepatitis", text))
    return cancer, pathogen


def assay_flags(*values: object) -> tuple[str, str, str]:
    text = " ".join(clean_text(v) for v in values).lower()
    activation = "activation_assay" if re.search(r"activation|elispot|ifn|tnf|il-2|cytotoxic|proliferation|cd69|cd137|killing", text) else ""
    tetramer = "tetramer_or_multimer" if re.search(r"tetramer|multimer|dextramer|streptamer", text) else ""
    cytokine = "cytokine_readout" if re.search(r"ifn|interferon|tnf|il-2|il2|cytokine", text) else ""
    return activation, tetramer, cytokine


def finalize_records(records: list[dict[str, object]]) -> pd.DataFrame:
    cols = CANONICAL_COLS + EXTRA_COLS
    if not records:
        return pd.DataFrame(columns=cols)
    reg = pd.DataFrame(records)
    for col in cols:
        if col not in reg.columns:
            reg[col] = pd.NA
    reg = reg[cols].copy()
    reg["peptide"] = reg["peptide"].map(normalize_peptide)
    reg["peptide_length"] = reg["peptide"].map(lambda s: len(s) if s else pd.NA)
    reg["hla_4digit"] = reg["hla_raw"].map(normalize_hla_any)
    reg["mhc_class"] = [
        normalize_mhc_class(c, h) for c, h in zip(reg["mhc_class"], reg["hla_4digit"])
    ]
    reg["mhc_gene"] = reg["hla_4digit"].map(mhc_gene_from_hla)
    for col in ["cdr3_alpha", "cdr3_beta"]:
        reg[col] = reg[col].map(normalize_cdr3)
    for col in ["tcr_alpha_v", "tcr_alpha_j", "tcr_beta_v", "tcr_beta_j", "mutation_gene"]:
        reg[col] = reg[col].map(normalize_gene)
    reg["paired_tcr_available"] = reg["cdr3_alpha"].astype(str).ne("") & reg["cdr3_beta"].astype(str).ne("")
    reg["beta_only_tcr_available"] = reg["cdr3_beta"].astype(str).ne("") & reg["cdr3_alpha"].astype(str).eq("")
    reg["alpha_only_tcr_available"] = reg["cdr3_alpha"].astype(str).ne("") & reg["cdr3_beta"].astype(str).eq("")
    reg["any_tcr_available"] = reg[["cdr3_alpha", "cdr3_beta"]].astype(str).ne("").any(axis=1)
    reg["peptide_hla_available"] = reg["peptide"].astype(str).ne("") & reg["hla_4digit"].astype(str).ne("")
    reg["hla_supertype"] = reg["hla_4digit"].map(hla_supertype)
    reg["near_peptide_cluster"] = reg["peptide"].map(peptide_cluster)
    reg["exact_peptide_hla_key"] = reg["peptide"].astype(str) + "|" + reg["hla_4digit"].astype(str)
    reg["peptide_only_key"] = reg["peptide"].astype(str)
    reg["near_peptide_cluster_key"] = reg["near_peptide_cluster"].astype(str) + "|" + reg["hla_supertype"].astype(str)
    reg["binding_label_binary"] = reg["binding_label_binary"].map(binary_label).astype("Int64")
    reg["row_id"] = [f"TCRREG_{i:07d}" for i in range(len(reg))]
    return reg


def parquet_compatible(df: pd.DataFrame) -> pd.DataFrame:
    """Make the heterogeneous public-resource registry stable for pyarrow."""
    out = df.copy()
    for col in out.columns:
        if pd.api.types.is_object_dtype(out[col]) or isinstance(out[col].dtype, pd.StringDtype):
            out[col] = out[col].where(out[col].notna(), "").astype(str)
    return out


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    view = df.copy()
    view = view.where(view.notna(), "NA").replace("", "NA")
    view.to_csv(path, sep="\t", index=False, na_rep="NA")


def join_values(values: Iterable[object], max_items: int = 6) -> str:
    seen: list[str] = []
    for value in values:
        text = clean_text(value)
        if text and text not in seen:
            seen.append(text)
        if len(seen) >= max_items:
            break
    return "|".join(seen)


def parse_vdjdb(path: Path) -> list[dict[str, object]]:
    if not path.exists():
        return []
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False, encoding_errors="replace")
    df["_idx"] = np.arange(len(df)).astype(str)
    complex_id = df.get("complex.id", pd.Series("", index=df.index)).fillna("").astype(str)
    valid_complex = complex_id.str.fullmatch(r"[1-9]\d*")
    df["_group"] = np.where(valid_complex, "complex:" + complex_id, "row:" + df["_idx"])
    records: list[dict[str, object]] = []
    for source_group, group in df.groupby("_group", sort=False):
        first = group.iloc[0]
        alpha = group[group.get("gene", "").astype(str).eq("TRA")]
        beta = group[group.get("gene", "").astype(str).eq("TRB")]
        method = parse_jsonish(first.get("method", ""))
        meta = parse_jsonish(first.get("meta", ""))
        assay = join_values([
            method.get("identification", ""),
            method.get("verification", ""),
            first.get("web.method", ""),
            first.get("web.method.seq", ""),
        ])
        activation, tetramer, cytokine = assay_flags(assay, first.get("method", ""))
        antigen_source = join_values([first.get("antigen.species", ""), first.get("antigen.gene", "")])
        disease_context = join_values([meta.get("subject.cohort", ""), meta.get("tissue", "")])
        cancer, pathogen = context_flags(disease_context, antigen_source)
        structure_id = first_nonempty(meta.get("structure.id", ""), meta.get("meta.structure.id", ""))
        records.append({
            "row_id": "",
            "source_dataset": "VDJdb",
            "source_record_id": source_group,
            "study_id": first.get("reference.id", ""),
            "assay_type": assay,
            "species": first.get("species", ""),
            "disease_context": disease_context,
            "cancer_type": disease_context if cancer else "",
            "antigen_source": antigen_source,
            "peptide": first.get("antigen.epitope", ""),
            "peptide_raw": first.get("antigen.epitope", ""),
            "hla_raw": first.get("mhc.a", ""),
            "mhc_class": first.get("mhc.class", ""),
            "tcr_alpha_v": join_values(alpha.get("v.segm", [])),
            "tcr_alpha_j": join_values(alpha.get("j.segm", [])),
            "cdr3_alpha": join_values(alpha.get("cdr3", [])),
            "tcr_beta_v": join_values(beta.get("v.segm", [])),
            "tcr_beta_j": join_values(beta.get("j.segm", [])),
            "cdr3_beta": join_values(beta.get("cdr3", [])),
            "clone_count": int(len(group)),
            "clone_frequency": first_nonempty(method.get("frequency", ""), meta.get("meta.subset.frequency", "")),
            "binding_label_raw": first_nonempty(first.get("vdjdb.score", ""), "VDJdb specificity record"),
            "binding_label_binary": 1,
            "activation_label": activation,
            "tetramer_label": tetramer,
            "cytokine_label": cytokine,
            "structure_pdb_id": structure_id,
            "train_test_original": "",
            "public_overlap_flags": "public_tcr_specificity_database;diagnostic_only",
            "source_path": str(path),
            "source_row_count": int(len(group)),
            "mutation_gene": first.get("antigen.gene", "") if clean_text(first.get("antigen.gene", "")).upper() in RECURRENT_GENES else "",
            "is_cancer_context": cancer,
            "is_pathogen_context": pathogen,
            "parser_notes": "VDJdb chain rows collapsed by nonzero complex.id; complex.id=0 retained as single-chain rows.",
        })
    return records


def parse_vdjdb_chunk(path: Path, source_name: str) -> list[dict[str, object]]:
    if not path.exists():
        return []
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False, encoding_errors="replace")
    records: list[dict[str, object]] = []
    for i, row in enumerate(df.to_dict("records")):
        assay = join_values([row.get("method.identification", ""), row.get("method.verification", ""), row.get("method.sequencing", "")])
        activation, tetramer, cytokine = assay_flags(assay)
        antigen_source = join_values([row.get("antigen.species", ""), row.get("antigen.gene", "")])
        disease_context = join_values([row.get("meta.subject.cohort", ""), row.get("meta.tissue", "")])
        cancer, pathogen = context_flags(disease_context, antigen_source)
        records.append({
            "row_id": "",
            "source_dataset": source_name,
            "source_record_id": f"{path.stem}_{i}",
            "study_id": row.get("reference.id", ""),
            "assay_type": assay,
            "species": row.get("species", ""),
            "disease_context": disease_context,
            "cancer_type": disease_context if cancer else "",
            "antigen_source": antigen_source,
            "peptide": row.get("antigen.epitope", ""),
            "peptide_raw": row.get("antigen.epitope", ""),
            "hla_raw": row.get("mhc.a", ""),
            "mhc_class": row.get("mhc.class", ""),
            "tcr_alpha_v": row.get("v.alpha", ""),
            "tcr_alpha_j": row.get("j.alpha", ""),
            "cdr3_alpha": row.get("cdr3.alpha", ""),
            "tcr_beta_v": row.get("v.beta", ""),
            "tcr_beta_j": row.get("j.beta", ""),
            "cdr3_beta": row.get("cdr3.beta", ""),
            "clone_count": 1,
            "clone_frequency": row.get("method.frequency", ""),
            "binding_label_raw": "VDJdb chunk specificity record",
            "binding_label_binary": 1,
            "activation_label": activation,
            "tetramer_label": tetramer,
            "cytokine_label": cytokine,
            "structure_pdb_id": row.get("meta.structure.id", ""),
            "public_overlap_flags": "public_tcr_specificity_database;diagnostic_only",
            "source_path": str(path),
            "source_row_count": 1,
            "mutation_gene": row.get("antigen.gene", "") if clean_text(row.get("antigen.gene", "")).upper() in RECURRENT_GENES else "",
            "is_cancer_context": cancer,
            "is_pathogen_context": pathogen,
            "parser_notes": f"{source_name} VDJdb chunk parsed as already-paired records.",
        })
    return records


def parse_mcpas(path: Path) -> list[dict[str, object]]:
    if not path.exists():
        return []
    df = pd.read_csv(path, dtype=str, low_memory=False, encoding_errors="replace")
    records: list[dict[str, object]] = []
    for i, row in enumerate(df.to_dict("records")):
        assay = join_values([
            row.get("Antigen.identification.method", ""),
            row.get("T.Cell.Type", ""),
            row.get("T.cell.characteristics", ""),
            row.get("Single.cell", ""),
            row.get("NGS", ""),
        ])
        activation, tetramer, cytokine = assay_flags(assay, row.get("Remarks", ""))
        disease = row.get("Pathology", "")
        antigen_source = join_values([row.get("Antigen.protein", ""), row.get("Protein.ID", "")])
        cancer, pathogen = context_flags(row.get("Category", ""), disease, antigen_source)
        records.append({
            "row_id": "",
            "source_dataset": "McPAS-TCR",
            "source_record_id": f"McPAS_{i:07d}",
            "study_id": row.get("PubMed.ID", ""),
            "assay_type": assay,
            "species": row.get("Species", ""),
            "disease_context": disease,
            "cancer_type": disease if cancer else "",
            "antigen_source": antigen_source,
            "peptide": row.get("Epitope.peptide", ""),
            "peptide_raw": row.get("Epitope.peptide", ""),
            "hla_raw": row.get("MHC", ""),
            "mhc_class": "",
            "tcr_alpha_v": row.get("TRAV", ""),
            "tcr_alpha_j": row.get("TRAJ", ""),
            "cdr3_alpha": row.get("CDR3.alpha.aa", ""),
            "tcr_beta_v": row.get("TRBV", ""),
            "tcr_beta_j": row.get("TRBJ", ""),
            "cdr3_beta": row.get("CDR3.beta.aa", ""),
            "clone_count": 1,
            "clone_frequency": "",
            "binding_label_raw": "McPAS curated specificity record",
            "binding_label_binary": 1,
            "activation_label": activation,
            "tetramer_label": tetramer,
            "cytokine_label": cytokine,
            "structure_pdb_id": "",
            "public_overlap_flags": "public_tcr_specificity_database;diagnostic_only",
            "source_path": str(path),
            "source_row_count": 1,
            "mutation_gene": row.get("Antigen.protein", "") if clean_text(row.get("Antigen.protein", "")).upper() in RECURRENT_GENES else "",
            "is_cancer_context": cancer,
            "is_pathogen_context": pathogen,
            "parser_notes": "McPAS rows are curated positive TCR/antigen associations; many are beta-only or pathology-only.",
        })
    return records


def parse_nepdb(path: Path) -> list[dict[str, object]]:
    if not path.exists():
        return []
    df = pd.read_csv(path, dtype=str, low_memory=False, encoding_errors="replace")
    records: list[dict[str, object]] = []
    for i, row in enumerate(df.to_dict("records")):
        hla_raw = first_nonempty(row.get("alleleA", ""), row.get("alleleB", ""))
        if hla_raw and not hla_raw.upper().startswith("HLA-"):
            hla_raw = "HLA-" + hla_raw
        assay = join_values([row.get("assay", ""), row.get("Tcell_source", ""), row.get("antigen_type", "")])
        activation, tetramer, cytokine = assay_flags(assay)
        cancer_type = row.get("Tumor Type", "")
        gene = row.get("genesymbol", "")
        mutation_key = join_values([gene, row.get("mut_aa_pos", ""), row.get("wt_aa", ""), row.get("mut_aa", "")])
        records.append({
            "row_id": "",
            "source_dataset": "NEPdb",
            "source_record_id": first_nonempty(row.get("id", ""), f"NEPdb_{i:07d}"),
            "study_id": row.get("PMID", ""),
            "assay_type": assay,
            "species": "Human",
            "disease_context": cancer_type,
            "cancer_type": cancer_type,
            "antigen_source": gene,
            "peptide": first_nonempty(row.get("minimal_peptide", ""), row.get("mut_peptide", "")),
            "peptide_raw": first_nonempty(row.get("minimal_peptide", ""), row.get("mut_peptide", "")),
            "hla_raw": hla_raw,
            "mhc_class": "I",
            "tcr_alpha_v": row.get("TRAV", ""),
            "tcr_alpha_j": row.get("TRAJ", ""),
            "cdr3_alpha": row.get("CDR3A", ""),
            "tcr_beta_v": row.get("TRBV", ""),
            "tcr_beta_j": row.get("TRBJ", ""),
            "cdr3_beta": row.get("CDR3B", ""),
            "clone_count": 1,
            "clone_frequency": "",
            "binding_label_raw": row.get("response", ""),
            "binding_label_binary": binary_label(row.get("response", "")),
            "activation_label": activation,
            "tetramer_label": tetramer,
            "cytokine_label": cytokine,
            "structure_pdb_id": "",
            "public_overlap_flags": "neoantigen_database;diagnostic_label_compatibility_required",
            "source_path": str(path),
            "source_row_count": 1,
            "mutation_gene": gene,
            "mutation_key": mutation_key,
            "is_cancer_context": True,
            "is_pathogen_context": False,
            "parser_notes": "NEPdb rows preserve positive and negative neoantigen response labels; labels are not transferred to CROSS-Neo rows.",
        })
    return records


def assign_iedb_chain(row: dict[str, object], chain_idx: int) -> tuple[str, str]:
    ctype = clean_text(row.get(f"receptor_chain{chain_idx}_types", "")).lower()
    cdr3 = row.get(f"receptor_chain{chain_idx}_cdr3_seqs", "")
    if "alpha" in ctype or "tra" in ctype:
        return "alpha", clean_text(cdr3)
    if "beta" in ctype or "trb" in ctype:
        return "beta", clean_text(cdr3)
    return "", clean_text(cdr3)


def parse_iedb_api(path: Path) -> list[dict[str, object]]:
    if not path.exists():
        return []
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False, encoding_errors="replace")
    records: list[dict[str, object]] = []
    for i, row in enumerate(df.to_dict("records")):
        chain1_type, chain1_cdr3 = assign_iedb_chain(row, 1)
        chain2_type, chain2_cdr3 = assign_iedb_chain(row, 2)
        cdr3_alpha = chain1_cdr3 if chain1_type == "alpha" else (chain2_cdr3 if chain2_type == "alpha" else "")
        cdr3_beta = chain1_cdr3 if chain1_type == "beta" else (chain2_cdr3 if chain2_type == "beta" else "")
        assay = join_values([row.get("assay_names", ""), row.get("assay_description", ""), row.get("qualitative_measure", "")])
        activation, tetramer, cytokine = assay_flags(assay)
        disease = row.get("disease_names", "")
        antigen_source = join_values([
            row.get("curated_source_antigen.name", ""),
            row.get("parent_source_antigen_name", ""),
            row.get("source_organism_name", ""),
        ])
        cancer, pathogen = context_flags(disease, antigen_source)
        records.append({
            "row_id": "",
            "source_dataset": "IEDB_tcell_api",
            "source_record_id": first_nonempty(row.get("tcell_id", ""), f"IEDB_API_{i:07d}"),
            "study_id": row.get("pubmed_id", ""),
            "assay_type": assay,
            "species": row.get("host_organism_name", ""),
            "disease_context": disease,
            "cancer_type": disease if cancer else "",
            "antigen_source": antigen_source,
            "peptide": row.get("linear_sequence", ""),
            "peptide_raw": row.get("linear_sequence", ""),
            "hla_raw": first_nonempty(row.get("mhc_allele_name", ""), row.get("mhc_restriction", "")),
            "mhc_class": row.get("mhc_class", ""),
            "tcr_alpha_v": "",
            "tcr_alpha_j": "",
            "cdr3_alpha": cdr3_alpha,
            "tcr_beta_v": "",
            "tcr_beta_j": "",
            "cdr3_beta": cdr3_beta,
            "clone_count": 1,
            "clone_frequency": "",
            "binding_label_raw": row.get("qualitative_measure", ""),
            "binding_label_binary": binary_label(row.get("qualitative_measure", "")),
            "activation_label": activation,
            "tetramer_label": tetramer,
            "cytokine_label": cytokine,
            "structure_pdb_id": row.get("pdb_id", ""),
            "public_overlap_flags": "IEDB_public_tcell_assay;training_overlap_risk;diagnostic_only",
            "source_path": str(path),
            "source_row_count": 1,
            "is_cancer_context": cancer,
            "is_pathogen_context": pathogen,
            "parser_notes": "IEDB IQ-API export; non-peptidic epitopes retained with empty normalized peptide.",
        })
    return records


def parse_iedb_legacy(path: Path, max_rows: int | None = None) -> list[dict[str, object]]:
    if not path.exists():
        return []
    keep = {
        "Assay_ID_IEDB_IRI", "Reference_PMID", "Epitope_Object_Type", "Epitope_Name",
        "Epitope_Source_Molecule", "Epitope_Source_Organism", "Host_Name",
        "1st_in_vivo_Process_Disease", "Assay_Method", "Assay_Response_measured",
        "Assay_Qualitative_Measurement", "Assay_Number_of_Subjects_Tested",
        "Assay_Number_of_Subjects_Positive", "Effector_Cell_TCR_Name", "Complex_PDB_ID",
        "MHC_Restriction_Name", "MHC_Restriction_Class", "Assay_Antigen_Source_Molecule",
        "Assay_Antigen_Source_Organism",
    }
    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        low_memory=False,
        usecols=lambda col: col in keep,
        nrows=max_rows,
        encoding_errors="replace",
    )
    records: list[dict[str, object]] = []
    for i, row in enumerate(df.to_dict("records")):
        assay = join_values([row.get("Assay_Method", ""), row.get("Assay_Response_measured", ""), row.get("Assay_Qualitative_Measurement", "")])
        activation, tetramer, cytokine = assay_flags(assay)
        disease = row.get("1st_in_vivo_Process_Disease", "")
        antigen_source = join_values([
            row.get("Epitope_Source_Molecule", ""),
            row.get("Assay_Antigen_Source_Molecule", ""),
            row.get("Epitope_Source_Organism", ""),
            row.get("Assay_Antigen_Source_Organism", ""),
        ])
        cancer, pathogen = context_flags(disease, antigen_source)
        records.append({
            "row_id": "",
            "source_dataset": "IEDB_tcell_full_export",
            "source_record_id": first_nonempty(row.get("Assay_ID_IEDB_IRI", ""), f"IEDB_FULL_{i:07d}"),
            "study_id": row.get("Reference_PMID", ""),
            "assay_type": assay,
            "species": row.get("Host_Name", ""),
            "disease_context": disease,
            "cancer_type": disease if cancer else "",
            "antigen_source": antigen_source,
            "peptide": row.get("Epitope_Name", "") if row.get("Epitope_Object_Type", "") == "Linear peptide" else "",
            "peptide_raw": row.get("Epitope_Name", ""),
            "hla_raw": row.get("MHC_Restriction_Name", ""),
            "mhc_class": row.get("MHC_Restriction_Class", ""),
            "tcr_alpha_v": "",
            "tcr_alpha_j": "",
            "cdr3_alpha": "",
            "tcr_beta_v": "",
            "tcr_beta_j": "",
            "cdr3_beta": "",
            "clone_count": row.get("Assay_Number_of_Subjects_Positive", ""),
            "clone_frequency": "",
            "binding_label_raw": row.get("Assay_Qualitative_Measurement", ""),
            "binding_label_binary": binary_label(row.get("Assay_Qualitative_Measurement", "")),
            "activation_label": activation,
            "tetramer_label": tetramer,
            "cytokine_label": cytokine,
            "structure_pdb_id": row.get("Complex_PDB_ID", ""),
            "public_overlap_flags": "IEDB_public_tcell_assay;training_overlap_risk;diagnostic_only",
            "source_path": str(path),
            "source_row_count": 1,
            "is_cancer_context": cancer,
            "is_pathogen_context": pathogen,
            "parser_notes": "IEDB legacy export has peptide-HLA/T-cell assay labels but usually no receptor sequence.",
        })
    return records


def discover_candidate_files() -> pd.DataFrame:
    roots = [REPO / "project", Path("/data/neoantigen_vaccine_hub")]
    rows: list[dict[str, object]] = []
    allowed = {".tsv", ".csv", ".txt", ".json", ".parquet", ".pkl", ".html", ".zip"}
    seen: set[str] = set()
    for root in roots:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            if ".git" in path.parts:
                continue
            if path.suffix.lower() not in allowed:
                continue
            spath = str(path)
            if spath in seen:
                continue
            seen.add(spath)
            hay = spath.lower()
            path_hits = [kw for kw in KEYWORDS if kw in hay]
            content_hits: list[str] = []
            if path.suffix.lower() in {".tsv", ".csv", ".txt", ".json", ".md"} and path.stat().st_size <= 20_000_000:
                try:
                    sample = path.read_text(errors="ignore")[:65536].lower()
                    content_hits = [kw for kw in KEYWORDS if kw in sample]
                except Exception:
                    content_hits = []
            if path_hits or content_hits:
                rows.append({
                    "path": spath,
                    "size_bytes": int(path.stat().st_size),
                    "path_keyword_hits": ",".join(path_hits),
                    "content_keyword_hits": ",".join(content_hits),
                    "parsed_status": "candidate_unparsed",
                    "parser_note": "",
                })
    inv = pd.DataFrame(rows).sort_values(["path"]).reset_index(drop=True)
    return inv


def update_inventory_status(inv: pd.DataFrame, parsed_paths: dict[str, str]) -> pd.DataFrame:
    if inv.empty:
        return inv
    inv = inv.copy()
    for path, note in parsed_paths.items():
        mask = inv["path"].eq(path)
        inv.loc[mask, "parsed_status"] = "parsed"
        inv.loc[mask, "parser_note"] = note
    duplicate_patterns = ["vdjdb.slim.txt", "vdjdb_full.txt", "tcell_full_v3.csv", "mcpas_input.csv", "mcpas_im.csv"]
    for pat in duplicate_patterns:
        mask = inv["path"].str.contains(re.escape(pat), case=False, regex=True, na=False) & inv["parsed_status"].eq("candidate_unparsed")
        inv.loc[mask, "parsed_status"] = "represented_or_derived_duplicate"
        inv.loc[mask, "parser_note"] = "Not parsed directly to avoid duplicate rows; represented by primary source or prediction-bundle derivative."
    return inv


def write_reports(reg: pd.DataFrame, inv: pd.DataFrame) -> None:
    TCR_OUT.mkdir(parents=True, exist_ok=True)
    source_counts = (
        reg.groupby("source_dataset", dropna=False)
        .agg(
            n=("row_id", "size"),
            paired_tcr=("paired_tcr_available", "sum"),
            beta_only_tcr=("beta_only_tcr_available", "sum"),
            alpha_only_tcr=("alpha_only_tcr_available", "sum"),
            any_tcr=("any_tcr_available", "sum"),
            peptide_hla=("peptide_hla_available", "sum"),
            positives=("binding_label_binary", "sum"),
            prevalence=("binding_label_binary", "mean"),
            cancer_context=("is_cancer_context", "sum"),
            pathogen_context=("is_pathogen_context", "sum"),
            structures=("structure_pdb_id", lambda s: int(s.fillna("").astype(str).ne("").sum())),
        )
        .reset_index()
    )
    write_tsv(source_counts, TCR_OUT / "tcr_registry_source_counts.tsv")

    miss = []
    for col in CANONICAL_COLS:
        missing = reg[col].isna() | reg[col].astype(str).isin(["", "nan", "None", "<NA>"])
        miss.append({"column": col, "missing_n": int(missing.sum()), "missing_frac": float(missing.mean()) if len(reg) else np.nan})
    write_tsv(pd.DataFrame(miss), TCR_OUT / "tcr_registry_missingness.tsv")

    label_balance = reg.groupby(["source_dataset", "binding_label_binary"], dropna=False).size().reset_index(name="n")
    write_tsv(label_balance, TCR_OUT / "tcr_registry_label_balance.tsv")
    write_tsv(inv, TCR_OUT / "tcr_registry_source_inventory.tsv")

    partition = pd.DataFrame([
        {"partition": "paired_alpha_beta", "n": int(reg["paired_tcr_available"].sum())},
        {"partition": "beta_only", "n": int(reg["beta_only_tcr_available"].sum())},
        {"partition": "alpha_only", "n": int(reg["alpha_only_tcr_available"].sum())},
        {"partition": "no_tcr_sequence", "n": int((~reg["any_tcr_available"]).sum())},
        {"partition": "peptide_hla_available", "n": int(reg["peptide_hla_available"].sum())},
        {"partition": "structure_pdb_id_available", "n": int(reg["structure_pdb_id"].fillna("").astype(str).ne("").sum())},
    ])
    write_tsv(partition, TCR_OUT / "tcr_registry_partition_counts.tsv")

    top_sources = source_counts.sort_values("n", ascending=False).head(12).to_markdown(index=False)
    lines = [
        "# CROSS-Neo-TCR Registry Report",
        "",
        "## Scope",
        "",
        "This registry is an optional TCR-aware evidence layer. It is not a replacement for the main CROSS-Neo pMHC ranking registry because most neoantigen rows do not contain paired TCR alpha/beta chains.",
        "",
        "## Headline Counts",
        "",
        f"- Registry rows retained: {len(reg):,}",
        f"- Rows with paired TCR alpha/beta CDR3: {int(reg['paired_tcr_available'].sum()):,}",
        f"- Rows with beta-only TCR CDR3: {int(reg['beta_only_tcr_available'].sum()):,}",
        f"- Rows with peptide + normalized HLA: {int(reg['peptide_hla_available'].sum()):,}",
        f"- Rows with cancer-context annotation: {int(reg['is_cancer_context'].sum()):,}",
        f"- Rows with pathogen-context annotation: {int(reg['is_pathogen_context'].sum()):,}",
        f"- Rows with PDB/structure identifier: {int(reg['structure_pdb_id'].fillna('').astype(str).ne('').sum()):,}",
        "",
        "## Source Balance",
        "",
        top_sources,
        "",
        "## Missingness Policy",
        "",
        "- No parsed source rows are silently discarded; non-peptidic IEDB rows are retained with `peptide_raw` populated and normalized `peptide` left empty.",
        "- Paired alpha/beta, beta-only, alpha-only, and no-TCR rows are explicitly flagged.",
        "- Public TCR evidence remains diagnostic unless labels and leakage boundaries are compatible with the downstream task.",
        "",
        "## Parsed Primary Local Sources",
        "",
        "- `/data/neoantigen_vaccine_hub/data_raw/tcr/vdjdb.txt`",
        "- `/data/neoantigen_vaccine_hub/data_raw/tcr/mcpas/McPAS-TCR.csv`",
        "- `/data/neoantigen_vaccine_hub/data_raw/nepdb/nepdb_all.csv`",
        "- `/data/neoantigen_vaccine_hub/data_raw/iedb/iedb_tcell_human_positive.tsv`",
        "- `/data/neoantigen_vaccine_hub/data_raw/iedb/iedb_tcell_human_cancer.tsv`",
        "- VDJdb 10x and PDB chunks are parsed as paired/structure supplements when present.",
        "",
        "## Outputs",
        "",
        "- `tcr_registry.parquet` / `tcr_registry.tsv`",
        "- `tcr_registry_source_counts.tsv`",
        "- `tcr_registry_missingness.tsv`",
        "- `tcr_registry_label_balance.tsv`",
        "- `tcr_registry_source_inventory.tsv`",
        "- `tcr_registry_partition_counts.tsv`",
    ]
    (TCR_OUT / "tcr_registry_report.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--max-iedb-legacy-rows",
        type=int,
        default=None,
        help="Optional cap for the large IEDB legacy export; default parses all rows.",
    )
    args = parser.parse_args()

    TCR_OUT.mkdir(parents=True, exist_ok=True)
    inv = discover_candidate_files()
    records: list[dict[str, object]] = []
    parsed_paths: dict[str, str] = {}

    source_paths = {
        "vdjdb": Path("/data/neoantigen_vaccine_hub/data_raw/tcr/vdjdb.txt"),
        "mcpas": Path("/data/neoantigen_vaccine_hub/data_raw/tcr/mcpas/McPAS-TCR.csv"),
        "nepdb": Path("/data/neoantigen_vaccine_hub/data_raw/nepdb/nepdb_all.csv"),
        "iedb_api": Path("/data/neoantigen_vaccine_hub/data_raw/iedb/iedb_tcell_human_positive.tsv"),
        "iedb_legacy": Path("/data/neoantigen_vaccine_hub/data_raw/iedb/iedb_tcell_human_cancer.tsv"),
        "vdjdb_10x": Path("/data/neoantigen_vaccine_hub/data_raw/tcr/github_antigenomics_vdjdb/chunks/10xgenomics-2019-07-09.txt"),
        "vdjdb_pdb": Path("/data/neoantigen_vaccine_hub/data_raw/tcr/github_antigenomics_vdjdb/chunks/PDB_Database.txt"),
    }

    for path, parser_fn, note in [
        (source_paths["vdjdb"], parse_vdjdb, "primary VDJdb full export"),
        (source_paths["mcpas"], parse_mcpas, "primary McPAS-TCR export"),
        (source_paths["nepdb"], parse_nepdb, "NEPdb neoantigen/TCR fields"),
        (source_paths["iedb_api"], parse_iedb_api, "IEDB IQ-API tcell export with receptor columns"),
    ]:
        if path.exists():
            rows = parser_fn(path)
            records.extend(rows)
            parsed_paths[str(path)] = note

    if source_paths["iedb_legacy"].exists():
        rows = parse_iedb_legacy(source_paths["iedb_legacy"], max_rows=args.max_iedb_legacy_rows)
        records.extend(rows)
        parsed_paths[str(source_paths["iedb_legacy"])] = "IEDB full tcell export with peptide-HLA assay labels"

    # Add paired chunk-level supplements that are not faithfully recoverable from chain-row-only records.
    for path, source_name in [
        (source_paths["vdjdb_10x"], "VDJdb_10x_chunk"),
        (source_paths["vdjdb_pdb"], "VDJdb_PDB_chunk"),
    ]:
        if path.exists():
            rows = parse_vdjdb_chunk(path, source_name)
            records.extend(rows)
            parsed_paths[str(path)] = f"{source_name} paired supplement"

    reg = finalize_records(records)
    inv = update_inventory_status(inv, parsed_paths)

    safe_parquet(parquet_compatible(reg), TCR_OUT / "tcr_registry.parquet")
    write_tsv(reg, TCR_OUT / "tcr_registry.tsv")
    write_reports(reg, inv)
    print(
        "[tcr-registry] "
        f"rows={len(reg)} paired={int(reg['paired_tcr_available'].sum())} "
        f"beta_only={int(reg['beta_only_tcr_available'].sum())} "
        f"peptide_hla={int(reg['peptide_hla_available'].sum())} out={TCR_OUT}"
    )


if __name__ == "__main__":
    main()
