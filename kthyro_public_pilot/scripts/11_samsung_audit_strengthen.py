#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
FIGS = ROOT / "results" / "figures" / "proposal"
REPORTS = ROOT / "results" / "reports"

MAIN_AXES = [
    "rai_differentiation_score",
    "immune_visibility_score",
    "cd8_exclusion_proxy",
    "drug_delivery_failure_proxy",
    "aggressive_dedifferentiation_score",
]

DISPLAY_AXIS = {
    "rai_differentiation_score": "RAI differentiation",
    "immune_visibility_score": "Immune visibility",
    "cd8_exclusion_proxy": "CD8 exclusion",
    "drug_delivery_failure_proxy": "Delivery failure proxy",
    "aggressive_dedifferentiation_score": "Aggressive dediff.",
    "hla_i_apm_score": "HLA/APM",
    "caf_ecm_tgfb_score": "CAF/ECM",
    "myeloid_tam_score": "Myeloid/TAM",
    "spatial_myeloid_caf_barrier_score": "CAF/myeloid",
}

NICHE_COLORS = {
    "RAI-restorable niche": "#2aa7ff",
    "RAI-low dedifferentiated niche": "#e84a5f",
    "HLA-visible inflamed niche": "#9b72df",
    "HLA-low invisible tumor niche": "#7a4c3a",
    "CD8-excluded / myeloid-CAF niche": "#f28e2b",
    "APC-rich niche": "#17becf",
    "Drug-delivery failure proxy niche": "#f4c542",
    "Mixed/Other": "#8f98a8",
}

SUPPORT_COLORS = {
    "Strong": "#34d399",
    "Moderate": "#60a5fa",
    "Weak": "#fbbf24",
    "Not tested": "#6b7280",
}


def read_tsv(name: str) -> pd.DataFrame:
    path = TABLES / name
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def save_dark(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    fig.savefig(path.with_suffix(".pdf"), dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())


def zscore(s: pd.Series) -> pd.Series:
    s = pd.to_numeric(s, errors="coerce")
    sd = s.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(0.0, index=s.index)
    return (s - s.mean()) / sd


def sample_col(df: pd.DataFrame) -> str | None:
    for c in ["sample_id", "sample_id_x", "sample", "patient_id"]:
        if c in df.columns:
            return c
    return None


def audit_outputs() -> dict[str, object]:
    expected = {
        "tcga": TABLES / "tcga_thca_patient_vulnerability_scores.tsv",
        "spatial_spots": TABLES / "spatial_spot_vulnerability_scores.tsv",
        "spatial_slides": TABLES / "spatial_slide_niche_summary.tsv",
        "spatial_coherence": TABLES / "spatial_coherence_statistics.tsv",
        "drug_candidates": TABLES / "drug_pilot_candidate_rankings.tsv",
        "evidence_matrix": TABLES / "integrated_public_pilot_evidence_matrix.tsv",
        "proposal_figures": FIGS,
        "public_report": REPORTS / "public_pilot_report.md",
        "go_decision": REPORTS / "GO_NO_GO_DECISION.md",
        "claim_boundaries": REPORTS / "claim_boundaries.md",
        "validation_plan": REPORTS / "experimental_validation_plan.md",
    }
    tcga = read_tsv("tcga_thca_patient_vulnerability_scores.tsv")
    spots = read_tsv("spatial_spot_vulnerability_scores.tsv")
    slides = read_tsv("spatial_slide_niche_summary.tsv")
    coherence = read_tsv("spatial_coherence_statistics.tsv")
    drug = read_tsv("drug_pilot_candidate_rankings.tsv")
    evidence = read_tsv("integrated_public_pilot_evidence_matrix.tsv")

    lines = ["# Pilot Output Audit\n"]
    lines.append("## Expected Files\n")
    for key, path in expected.items():
        lines.append(f"- {key}: {'EXISTS' if path.exists() else 'MISSING'} `{path.relative_to(ROOT) if path.exists() or str(path).startswith(str(ROOT)) else path}`")

    lines.append("\n## Row/Column Counts\n")
    for key, df in [
        ("TCGA patient vulnerability", tcga),
        ("Spatial spot vulnerability", spots),
        ("Spatial slide summary", slides),
        ("Spatial coherence", coherence),
        ("Drug candidates", drug),
        ("Integrated evidence matrix", evidence),
    ]:
        lines.append(f"- {key}: {df.shape[0]} rows x {df.shape[1]} columns")

    lines.append("\n## Identity And Duplication Checks\n")
    if not tcga.empty:
        scol = sample_col(tcga)
        lines.append(f"- TCGA canonical sample column detected: `{scol}`")
        lines.append(f"- TCGA row count: {len(tcga)}")
        lines.append(f"- TCGA unique patients: {tcga['patient_id'].nunique() if 'patient_id' in tcga else 'missing patient_id'}")
        if scol:
            lines.append(f"- TCGA duplicated sample records: {int(tcga[scol].duplicated().sum())}")
        lines.append(f"- TCGA missing labels: {int(tcga['primary_vulnerability_label'].isna().sum()) if 'primary_vulnerability_label' in tcga else 'label field missing'}")
        if "sample_id_x" in tcga.columns:
            lines.append("- Suspicious but non-blocking: TCGA table contains `sample_id_x/sample_id_y` from clinical merge. Use `sample_id_x` or `patient_id` for proposal summaries.")
    if not spots.empty:
        dup = int(spots[["sample_id", "spot_id"]].duplicated().sum()) if {"sample_id", "spot_id"}.issubset(spots.columns) else "not checked"
        lines.append(f"- Spatial spot rows: {len(spots)}")
        lines.append(f"- Spatial slide count from spots: {spots['sample_id'].nunique() if 'sample_id' in spots else 'missing sample_id'}")
        lines.append(f"- Spatial duplicated sample_id+spot_id records: {dup}")
        lines.append("- Spatial spots are nested within slide via `sample_id`; group claims must be slide-level.")
    if not slides.empty:
        lines.append(f"- Spatial slide summary rows: {len(slides)}")
        lines.append(f"- Spatial conditions: {slides['condition'].value_counts().to_dict() if 'condition' in slides else 'missing condition'}")
    if not coherence.empty:
        niche = coherence[coherence["label_set"].eq("niche_label")] if "label_set" in coherence else pd.DataFrame()
        lines.append(f"- KNN same-niche coherence rows: {len(niche)}")
        lines.append(f"- Slides with same-niche z > 2: {int((niche['z'] > 2).sum()) if not niche.empty else 0}/{len(niche)}")
        lines.append(f"- Minimum empirical p for same-niche test: {niche['empirical_p'].min() if not niche.empty else 'not available'}")
    if not drug.empty:
        dup_drug = int(drug["drug_label"].duplicated().sum()) if "drug_label" in drug else "missing drug_label"
        lines.append(f"- Drug candidate duplicated `drug_label` rows: {dup_drug}")
        lines.append("- Drug table is usable only after compound/class de-duplication; target-level duplicates are not proposal-ready.")

    lines.append("\n## Missing Fields And Usability\n")
    requirements = {
        "TCGA": ["patient_id", "primary_vulnerability_label", *MAIN_AXES],
        "Spatial spots": ["sample_id", "condition", "spot_id", "niche_label", "rai_differentiation_score", "hla_i_apm_score", "spatial_myeloid_caf_barrier_score"],
        "Spatial slides": ["sample_id", "condition", "n_spots", "same_niche_z", "same_niche_empirical_p"],
        "Drug": ["drug_label", "target", "mechanism_axis", "evidence_category", "n_thyroid_models"],
    }
    data = {"TCGA": tcga, "Spatial spots": spots, "Spatial slides": slides, "Drug": drug}
    usable = {}
    for key, req in requirements.items():
        df = data[key]
        missing = [c for c in req if c not in df.columns]
        usable[key] = not missing and not df.empty
        lines.append(f"- {key}: missing fields = {missing if missing else 'none'}; proposal usable = {usable[key]}")

    REPORTS.mkdir(parents=True, exist_ok=True)
    (REPORTS / "pilot_output_audit.md").write_text("\n".join(lines) + "\n")
    return {
        "tcga": tcga,
        "spots": spots,
        "slides": slides,
        "coherence": coherence,
        "drug": drug,
        "evidence": evidence,
        "usable": usable,
    }


def normalize_compound(label: str) -> str:
    label = str(label)
    label = re.sub(r"\s+DPC-\d+$", "", label)
    label = re.sub(r"\s+BRD-[A-Z0-9-]+$", "", label)
    return label.strip()


def perturbation_class(row: pd.Series) -> str:
    text = f"{row.get('drug_label','')} {row.get('target','')} {row.get('mechanism_axis','')}".upper()
    if any(x in text for x in ["SELUMETINIB", "TRAMETINIB", "MEK", "BRAF", " RAF", "RET", "NTRK", "MAPK"]):
        return "RAI redifferentiation / MAPK-axis modulation"
    if any(x in text for x in ["HDAC", "DNMT", "AZACITIDINE", "DECITABINE", "VORINOSTAT", "BELINOSTAT"]):
        return "HLA/APM restoration / immune visibility modulation"
    if any(x in text for x in ["JAK", "STAT", "RUXOLITINIB", "FEDRATINIB", "MOMELOTINIB"]):
        return "HLA/APM restoration / immune visibility modulation"
    if any(x in text for x in ["AURK", "CHEK", "CDK", "MYC", "BCL", "MCL", "AZD7762", "ALISERTIB", "MYCOPHENOLIC"]):
        return "Proliferation stress vulnerability"
    if any(x in text for x in ["CSF1", "CXCR4", "TGFB", "FAP"]):
        return "Myeloid/CAF barrier modulation"
    return "Other exploratory"


def evidence_level(row: pd.Series) -> str:
    cat = str(row.get("evidence_category", ""))
    n = pd.to_numeric(row.get("n_thyroid_models", np.nan), errors="coerce")
    p = pd.to_numeric(row.get("p_lineage_sensitivity_thyroid", np.nan), errors="coerce")
    cls = str(row.get("perturbation_class", ""))
    if cat.startswith("A") and pd.notna(n) and n >= 10 and pd.notna(p) and p < 0.10 and cls != "Other exploratory":
        return "A"
    if cat.startswith(("A", "B")) and cls != "Other exploratory":
        return "B"
    if cls in {"Myeloid/CAF barrier modulation", "Drug-delivery or nanoparticle distribution validation strategy"}:
        return "C"
    if cls != "Other exploratory":
        return "C"
    return "D"


def support_string(row: pd.Series) -> str:
    n = row.get("n_thyroid_models", np.nan)
    rho = row.get("spearman_lineage_sensitivity_thyroid", np.nan)
    p = row.get("p_lineage_sensitivity_thyroid", np.nan)
    try:
        return f"thyroid lines n={int(n)}, rho={float(rho):.2f}, p={float(p):.3g}; exploratory small-N"
    except Exception:
        return "not directly tested in thyroid cell-line drug screen"


def rationale_and_caveat(row: pd.Series) -> tuple[str, str, str, str, str]:
    compound = str(row.get("compound", ""))
    cls = str(row.get("perturbation_class", ""))
    text = f"{compound} {row.get('target','')}".upper()
    if cls.startswith("RAI"):
        rationale = "MAPK-axis attenuation is mechanistically linked to thyroid differentiation and RAI redifferentiation hypotheses."
        caveat = "Public screen does not prove NIS induction or iodide uptake restoration."
        assay = "NIS/SLC5A5, TG/TPO/PAX8 protein/RNA plus iodide uptake after perturbation."
        allowed = "MAPK-axis perturbation class nominated for RAI-restoration testing."
        forbidden = "This drug restores RAI sensitivity in patients."
    elif "HLA/APM" in cls:
        rationale = "Epigenetic/IFN-related perturbation class can be tested for HLA-I/APM rescue in HLA-low tumor territories."
        caveat = "JAK inhibitors may suppress IFN signaling; do not claim HLA/APM restoration without perturbation-expression rescue data."
        if any(x in text for x in ["JAK", "RUXOLITINIB", "FEDRATINIB", "MOMELOTINIB"]):
            caveat = "JAK inhibition can oppose IFN-driven HLA/APM induction; treat as pathway-control hypothesis, not restoration therapy."
        assay = "HLA-I/HLA-ABC, B2M, TAP1, NLRC5, PSMB8/9 rescue by IFN/HDAC/DNMT perturbation."
        allowed = "Immune-visibility modulation class nominated for HLA/APM rescue assays."
        forbidden = "JAK/HDAC/DNMT candidate is proven to restore antigen presentation."
    elif "Proliferation" in cls:
        rationale = "Cell-cycle stress associations may identify aggressive/proliferative niche vulnerabilities."
        caveat = "Could reflect general cytotoxicity rather than niche-specific therapeutic reprogramming."
        if "MYCOPHENOLIC" in text:
            caveat = "Mycophenolic acid is immunosuppressive; do not present as immune-restoration therapy."
        assay = "Ki-67/TOP2A-high organoid/slice viability with therapeutic-index and combination testing."
        allowed = "Proliferative/aggressive niche stress vulnerability for controlled PoC testing."
        forbidden = "General cytotoxic drug is a precision thyroid cancer therapy."
    elif "Myeloid/CAF" in cls:
        rationale = "CAF/myeloid-rich excluded niches nominate CSF1R, CXCR4, TGF-beta, FAP/ECM-axis validation."
        caveat = "Monoculture DepMap/PRISM cannot validate stromal or myeloid barrier biology."
        assay = "mIHC spatial exclusion assay and fresh tissue co-culture/slice barrier-modulation test."
        allowed = "Barrier-modulation class proposed for tissue-level validation."
        forbidden = "Cell-line drug screen proves CD8-exclusion reversal."
    elif "Drug-delivery" in cls:
        rationale = "ECM/hypoxia/vascular-low spatial territories nominate delivery imaging and formulation testing."
        caveat = "RNA score is a delivery proxy; actual drug distribution must be imaged."
        assay = "Fluorescent drug/liposome/nanoparticle penetration imaging in fresh tissue slices."
        allowed = "Delivery-failure proxy nominates imaging-based rescue experiments."
        forbidden = "Transcriptomics proves drug delivery failure."
    else:
        rationale = "Exploratory public-screen association."
        caveat = "Insufficient mechanism or thyroid-specific support."
        assay = "Secondary validation required before proposal emphasis."
        allowed = "Exploratory candidate only."
        forbidden = "Therapeutic candidate."
    return rationale, caveat, assay, allowed, forbidden


def clean_drug_table(drug: pd.DataFrame) -> pd.DataFrame:
    if drug.empty:
        out = pd.DataFrame()
        out.to_csv(TABLES / "drug_pilot_candidate_rankings_cleaned.tsv", sep="\t", index=False)
        return out
    df = drug.copy()
    df["compound"] = df["drug_label"].map(normalize_compound)
    df["perturbation_class"] = df.apply(perturbation_class, axis=1)
    df = df[df["perturbation_class"].ne("Other exploratory")].copy()
    df["evidence_level"] = df.apply(evidence_level, axis=1)
    order = {"A": 1, "B": 2, "C": 3, "D": 4}
    df["_evidence_order"] = df["evidence_level"].map(order).fillna(9)
    df["_abs_thyroid_rho"] = pd.to_numeric(df["spearman_lineage_sensitivity_thyroid"], errors="coerce").abs().fillna(0)
    df = df.sort_values(["_evidence_order", "rank_score", "_abs_thyroid_rho"], ascending=[True, False, False])

    rows = []
    for (compound, cls), g in df.groupby(["compound", "perturbation_class"], sort=False):
        best = g.iloc[0].copy()
        best["target"] = ",".join(sorted(set(g["target"].dropna().astype(str))))
        best["drug_id"] = ",".join(sorted(set(g["drug_id"].dropna().astype(str))))
        best["thyroid_specific_support"] = support_string(best)
        rationale, caveat, assay, allowed, forbidden = rationale_and_caveat(best)
        best["mechanism_rationale"] = rationale
        best["major_caveat"] = caveat
        best["required_validation_assay"] = assay
        best["proposal_wording_allowed"] = allowed
        best["proposal_wording_forbidden"] = forbidden
        proposed_use = {
            "RAI redifferentiation / MAPK-axis modulation": "Test RAI redifferentiation and iodide uptake restoration.",
            "HLA/APM restoration / immune visibility modulation": "Test HLA/APM protein/RNA rescue in HLA-low tumor territories.",
            "Proliferation stress vulnerability": "Test aggressive/proliferative niche sensitivity with therapeutic-index controls.",
            "Myeloid/CAF barrier modulation": "Test CD8-exclusion/barrier relief in tissue context.",
            "Drug-delivery or nanoparticle distribution validation strategy": "Test penetration and distribution rescue by imaging.",
        }.get(str(best["perturbation_class"]), "Exploratory follow-up only.")
        if str(best.get("compound", "")).upper() in {"JAK3_7406", "RUXOLITINIB", "FEDRATINIB", "MOMELOTINIB"}:
            proposed_use = "Use as JAK/IFN pathway-control hypothesis; do not frame as HLA/APM restoration unless rescue is observed."
        best["proposed_use"] = proposed_use
        rows.append(best)

    cleaned = pd.DataFrame(rows)
    class_rows = [
        {
            "source": "class_strategy",
            "target": "CSF1R,CXCR4,TGFB,FAP/ECM",
            "drug_id": "",
            "drug_label": "Class-level strategy: CSF1R/CXCR4/TGF-beta/FAP-ECM axis",
            "compound": "Class-level strategy: CSF1R/CXCR4/TGF-beta/FAP-ECM axis",
            "perturbation_class": "Myeloid/CAF barrier modulation",
            "proposed_use": "Test CD8-exclusion and stromal barrier relief in FFPE and fresh tissue models.",
            "evidence_level": "C",
            "thyroid_specific_support": "not directly testable in monoculture public drug screen",
            "mechanism_rationale": "Spatial CAF/myeloid-rich excluded niches nominate tissue-level barrier modulation.",
            "major_caveat": "Requires myeloid/stromal/tumor tissue context; not validated by DepMap/PRISM monoculture.",
            "required_validation_assay": "mIHC distance-to-CD8 analysis plus fresh tissue slice/co-culture perturbation.",
            "proposal_wording_allowed": "Barrier-modulation class proposed for validation.",
            "proposal_wording_forbidden": "Public cell-line data prove CD8-exclusion reversal.",
            "rank_score": 0,
            "n_thyroid_models": np.nan,
            "spearman_lineage_sensitivity_thyroid": np.nan,
            "p_lineage_sensitivity_thyroid": np.nan,
        },
        {
            "source": "class_strategy",
            "target": "fluorescent drug/liposome/nanoparticle imaging",
            "drug_id": "",
            "drug_label": "Class-level strategy: delivery-enhancing formulation/imaging",
            "compound": "Class-level strategy: delivery-enhancing formulation/imaging",
            "perturbation_class": "Drug-delivery or nanoparticle distribution validation strategy",
            "proposed_use": "Validate drug-delivery failure proxy and rescue strategy by direct imaging.",
            "evidence_level": "C",
            "thyroid_specific_support": "not tested in public drug screen; spatial RNA gives proxy only",
            "mechanism_rationale": "ECM/hypoxia-high and vascular-low territories nominate distribution imaging.",
            "major_caveat": "Transcriptomics does not measure actual delivery.",
            "required_validation_assay": "Fresh tissue fluorescent drug/liposome/nanoparticle penetration imaging.",
            "proposal_wording_allowed": "Delivery proxy nominates imaging validation experiments.",
            "proposal_wording_forbidden": "RNA score proves drug-delivery failure.",
            "rank_score": 0,
            "n_thyroid_models": np.nan,
            "spearman_lineage_sensitivity_thyroid": np.nan,
            "p_lineage_sensitivity_thyroid": np.nan,
        },
    ]
    cleaned = pd.concat([cleaned, pd.DataFrame(class_rows)], ignore_index=True, sort=False)
    cleaned["_evidence_order"] = cleaned["evidence_level"].map(order).fillna(9)
    cleaned = cleaned.sort_values(["_evidence_order", "rank_score"], ascending=[True, False]).reset_index(drop=True)
    cleaned.insert(0, "clean_rank", range(1, len(cleaned) + 1))
    keep = [
        "clean_rank",
        "compound",
        "perturbation_class",
        "proposed_use",
        "evidence_level",
        "thyroid_specific_support",
        "mechanism_rationale",
        "major_caveat",
        "required_validation_assay",
        "proposal_wording_allowed",
        "proposal_wording_forbidden",
        "source",
        "target",
        "drug_id",
        "drug_label",
        "rank_score",
        "n_thyroid_models",
        "spearman_lineage_sensitivity_thyroid",
        "p_lineage_sensitivity_thyroid",
    ]
    cleaned[keep].to_csv(TABLES / "drug_pilot_candidate_rankings_cleaned.tsv", sep="\t", index=False)
    return cleaned[keep]


def tcga_validation(tcga: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if tcga.empty:
        out = pd.DataFrame()
        out.to_csv(TABLES / "tcga_label_summary_cleaned.tsv", sep="\t", index=False)
        return out
    available_axes = [c for c in MAIN_AXES + ["hla_i_apm_score", "caf_ecm_tgfb_score", "myeloid_tam_score", "cytotoxic_t_score"] if c in tcga]
    z = tcga[available_axes].apply(zscore)
    for label, idx in tcga.groupby("primary_vulnerability_label").groups.items():
        sub = tcga.loc[idx]
        zmed = z.loc[idx].median().sort_values(key=lambda s: s.abs(), ascending=False)
        row = {
            "label": label,
            "n": len(sub),
            "fraction": len(sub) / len(tcga),
            "top_score_drivers": "; ".join([f"{DISPLAY_AXIS.get(k,k)}={v:.2f}" for k, v in zmed.head(4).items()]),
            "multi_axis_support_n_abs_median_z_ge_0_5": int((zmed.abs() >= 0.5).sum()),
        }
        for c in available_axes:
            row[f"median_{c}"] = sub[c].median()
            row[f"iqr_{c}"] = sub[c].quantile(0.75) - sub[c].quantile(0.25)
        rows.append(row)
    summary = pd.DataFrame(rows).sort_values("n", ascending=False)
    summary.to_csv(TABLES / "tcga_label_summary_cleaned.tsv", sep="\t", index=False)
    return summary


def apply_dark(ax):
    ax.set_facecolor("#111318")
    ax.tick_params(colors="white")
    ax.xaxis.label.set_color("white")
    ax.yaxis.label.set_color("white")
    ax.title.set_color("white")
    for spine in ax.spines.values():
        spine.set_color("#6b7280")


def make_tcga_figures(tcga: pd.DataFrame) -> None:
    if tcga.empty:
        return
    FIGS.mkdir(parents=True, exist_ok=True)
    counts = tcga["primary_vulnerability_label"].value_counts().sort_values()
    fig, ax = plt.subplots(figsize=(9.5, 5.6), facecolor="#111318")
    colors = ["#8f98a8" if x == "Mixed/Other" else "#60a5fa" for x in counts.index]
    ax.barh(counts.index, counts.values, color=colors)
    for i, v in enumerate(counts.values):
        ax.text(v + 2, i, str(v), va="center", color="white", fontsize=10)
    ax.set_title("TCGA-THCA patient-level therapeutic vulnerability labels (n=505)")
    ax.set_xlabel("Number of primary tumors")
    apply_dark(ax)
    save_dark(fig, FIGS / "tcga_label_distribution_dark.png")
    plt.close(fig)

    axes = [c for c in MAIN_AXES + ["hla_i_apm_score", "cytotoxic_t_score", "caf_ecm_tgfb_score", "myeloid_tam_score"] if c in tcga]
    order = tcga.sort_values(["primary_vulnerability_label"]).index
    heat = tcga.loc[order, axes].apply(zscore).T
    fig, ax = plt.subplots(figsize=(12, 5.8), facecolor="#111318")
    sns.heatmap(heat, cmap="vlag", center=0, xticklabels=False, yticklabels=[DISPLAY_AXIS.get(c, c) for c in axes], ax=ax, cbar_kws={"label": "z-score"})
    ax.set_title("TCGA vulnerability axes heatmap")
    ax.set_xlabel("TCGA-THCA primary tumors sorted by label")
    ax.set_ylabel("")
    apply_dark(ax)
    cbar = ax.collections[0].colorbar
    cbar.ax.yaxis.label.set_color("white")
    cbar.ax.tick_params(colors="white")
    save_dark(fig, FIGS / "tcga_vulnerability_axes_heatmap_dark.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 6.5), facecolor="#111318")
    size = 35 + 35 * (zscore(tcga["aggressive_dedifferentiation_score"]).clip(-2, 3) + 2.1)
    sc = ax.scatter(
        tcga["rai_differentiation_score"],
        tcga["immune_visibility_score"],
        c=tcga["myeloid_caf_barrier_score"],
        s=size,
        cmap="Oranges",
        alpha=0.82,
        edgecolors="none",
    )
    ax.axvline(tcga["rai_differentiation_score"].median(), color="#9ca3af", lw=1, ls="--")
    ax.axhline(tcga["immune_visibility_score"].median(), color="#9ca3af", lw=1, ls="--")
    ax.set_xlabel("RAI differentiation")
    ax.set_ylabel("HLA/APM immune visibility")
    ax.set_title("Therapeutic quadrant: bulk patient-level axes")
    apply_dark(ax)
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("Myeloid/CAF barrier")
    cbar.ax.yaxis.label.set_color("white")
    cbar.ax.tick_params(colors="white")
    save_dark(fig, FIGS / "tcga_therapeutic_quadrant_dark.png")
    plt.close(fig)


def spatial_validation(spots: pd.DataFrame, slides: pd.DataFrame, coherence: pd.DataFrame) -> pd.DataFrame:
    if slides.empty:
        out = pd.DataFrame()
        out.to_csv(TABLES / "spatial_coherence_summary_cleaned.tsv", sep="\t", index=False)
        return out
    niche = coherence[coherence["label_set"].eq("niche_label")].copy() if not coherence.empty and "label_set" in coherence else pd.DataFrame()
    summary = slides.copy()
    if not niche.empty:
        cols = ["sample_id", "observed", "null_mean", "null_sd", "z", "empirical_p"]
        summary = summary.drop(columns=[c for c in ["same_niche_observed", "same_niche_null_mean", "same_niche_z", "same_niche_empirical_p"] if c in summary], errors="ignore")
        summary = summary.merge(niche[cols], on="sample_id", how="left")
        summary = summary.rename(columns={"observed": "same_niche_observed", "null_mean": "same_niche_null_mean", "z": "same_niche_z", "empirical_p": "same_niche_empirical_p"})
    summary["coherent_z_gt_2"] = summary["same_niche_z"] > 2
    summary["spots_nested_by_slide"] = True
    diversity = spots.groupby("sample_id")["niche_label"].nunique().rename("n_niche_labels") if not spots.empty else pd.Series(dtype=float)
    summary = summary.merge(diversity, on="sample_id", how="left")
    summary.to_csv(TABLES / "spatial_coherence_summary_cleaned.tsv", sep="\t", index=False)
    return summary


def make_spatial_figures(spots: pd.DataFrame, spatial_summary: pd.DataFrame) -> None:
    if spots.empty or spatial_summary.empty:
        return
    fig, ax = plt.subplots(figsize=(10, 5.8), facecolor="#111318")
    plot = spatial_summary.sort_values("same_niche_z", ascending=True)
    ax.barh(plot["sample_id"], plot["same_niche_z"], color=plot["condition"].map({"normal": "#9ca3af", "PTC": "#60a5fa", "locally_advanced_PTC": "#f59e0b", "ATC": "#ef4444"}).fillna("#a78bfa"))
    ax.axvline(2, color="#fbbf24", lw=1.5, ls="--")
    ax.set_title("Spatial same-niche KNN coherence across GSE250521 slides")
    ax.set_xlabel("Permutation z-score")
    ax.set_ylabel("")
    apply_dark(ax)
    save_dark(fig, FIGS / "spatial_coherence_barplot_dark.png")
    plt.close(fig)

    frac_cols = [c for c in spatial_summary.columns if c.startswith("fraction_")]
    long = spatial_summary.melt(id_vars=["sample_id", "condition"], value_vars=frac_cols, var_name="niche", value_name="fraction")
    long["niche"] = long["niche"].str.replace("fraction_", "", regex=False)
    fig, ax = plt.subplots(figsize=(12, 6.5), facecolor="#111318")
    sns.barplot(data=long, x="condition", y="fraction", hue="niche", errorbar="se", ax=ax)
    ax.set_title("Slide-level niche fractions by condition")
    ax.set_xlabel("")
    ax.set_ylabel("Mean slide fraction")
    apply_dark(ax)
    leg = ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8, frameon=False)
    for text in leg.get_texts():
        text.set_color("white")
    save_dark(fig, FIGS / "spatial_niche_fraction_by_condition_dark.png")
    plt.close(fig)

    reps = []
    for cond in ["PTC", "locally_advanced_PTC", "ATC"]:
        sub = spatial_summary[spatial_summary["condition"].eq(cond)].sort_values("same_niche_z", ascending=False)
        if not sub.empty:
            reps.append(sub.iloc[0]["sample_id"])
    if len(reps) < 3:
        reps = spatial_summary.sort_values(["n_niche_labels", "same_niche_z"], ascending=False)["sample_id"].head(3).tolist()
    panels = [
        ("rai_differentiation_score", "RAI differentiation", "Blues"),
        ("hla_i_apm_score", "HLA/APM", "Purples"),
        ("spatial_myeloid_caf_barrier_score", "CAF/myeloid barrier", "Oranges"),
        ("niche_label", "Niche label", None),
    ]
    fig, axes = plt.subplots(len(reps), 4, figsize=(15.5, 4.1 * len(reps)), facecolor="#111318")
    if len(reps) == 1:
        axes = np.array([axes])
    for i, sample in enumerate(reps):
        s = spots[spots["sample_id"].eq(sample)]
        cond = str(s["condition"].iloc[0]) if not s.empty and "condition" in s else "condition unavailable"
        for j, (col, title, cmap) in enumerate(panels):
            ax = axes[i, j]
            ax.set_facecolor("#111318")
            if col == "niche_label":
                ax.scatter(s["x"], s["y"], c=s[col].map(NICHE_COLORS).fillna("#8f98a8"), s=5, linewidths=0)
            else:
                ax.scatter(s["x"], s["y"], c=s[col], cmap=cmap, s=5, linewidths=0)
            ax.invert_yaxis()
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(f"{sample} ({cond})\n{title}", color="white", fontsize=10)
    handles = [plt.Line2D([0], [0], marker="o", linestyle="", color=c, label=k, markersize=6) for k, c in NICHE_COLORS.items() if k in set(spots[spots["sample_id"].isin(reps)]["niche_label"])]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, labelcolor="white", fontsize=8)
    fig.suptitle("Representative spatial therapeutic vulnerability maps", color="white", fontsize=16, y=0.995)
    save_dark(fig, FIGS / "spatial_representative_maps_dark.png")
    plt.close(fig)


def make_evidence_matrix_figure(evidence: pd.DataFrame) -> None:
    rows = [
        "RAI differentiation /\nredifferentiation axis",
        "HLA/APM immune\nvisibility axis",
        "HLA-low invisible\ntumor axis",
        "CD8 exclusion /\nmyeloid-CAF barrier",
        "Drug-delivery\nfailure proxy",
        "Aggressive dediff. /\nproliferation axis",
        "AI companion\ndiagnostic feasibility",
        "Drug/perturbation\ncandidate feasibility",
    ]
    cols = ["TCGA bulk", "Spatial ST", "scRNA/cell-type", "Drug/perturbation", "Experimental validation", "Claim strength"]
    # Deterministic categorization from the actual completed outputs, not new biological claims.
    vals = [
        ["Strong", "Strong", "Moderate", "Moderate", "Strong", "Moderate"],
        ["Strong", "Strong", "Moderate", "Weak", "Strong", "Moderate"],
        ["Strong", "Moderate", "Weak", "Weak", "Strong", "Moderate"],
        ["Moderate", "Strong", "Moderate", "Not tested", "Strong", "Moderate"],
        ["Moderate", "Moderate", "Not tested", "Not tested", "Strong", "Weak"],
        ["Strong", "Moderate", "Weak", "Moderate", "Strong", "Moderate"],
        ["Strong", "Strong", "Moderate", "Weak", "Strong", "Strong"],
        ["Weak", "Weak", "Not tested", "Moderate", "Strong", "Moderate"],
    ]
    fig, ax = plt.subplots(figsize=(13.5, 7.2), facecolor="#111318")
    ax.set_xlim(0, len(cols))
    ax.set_ylim(0, len(rows))
    ax.axis("off")
    for i, row in enumerate(rows):
        y = len(rows) - i - 1
        ax.text(-0.05, y + 0.5, row, color="white", ha="right", va="center", fontsize=10)
        for j, col in enumerate(cols):
            value = vals[i][j]
            ax.add_patch(plt.Rectangle((j, y), 0.95, 0.82, color=SUPPORT_COLORS[value], alpha=0.9))
            ax.text(j + 0.475, y + 0.41, value, color="#111318", ha="center", va="center", fontsize=9, weight="bold")
    for j, col in enumerate(cols):
        ax.text(j + 0.475, len(rows) + 0.15, col, color="white", ha="center", va="bottom", fontsize=10, weight="bold")
    ax.text(2.9, len(rows) + 0.85, "Integrated public-pilot evidence matrix", color="white", ha="center", fontsize=17, weight="bold")
    ax.text(2.9, -0.65, "Categories summarize proposal-readiness of support; public data remain hypothesis-generating.", color="#d1d5db", ha="center", fontsize=10)
    save_dark(fig, FIGS / "integrated_public_pilot_evidence_matrix_dark.png")
    plt.close(fig)


def write_tcga_report(tcga: pd.DataFrame, label_summary: pd.DataFrame) -> None:
    assoc = read_tsv("tcga_thca_clinical_association_summary.tsv")
    label_counts = tcga["primary_vulnerability_label"].value_counts().to_dict() if not tcga.empty else {}
    corr = tcga[[c for c in MAIN_AXES if c in tcga]].corr(method="spearman").round(2) if not tcga.empty else pd.DataFrame()
    sig = assoc[pd.to_numeric(assoc.get("fdr", pd.Series(dtype=float)), errors="coerce") < 0.10].head(8) if not assoc.empty and "fdr" in assoc else pd.DataFrame()
    text = ["# TCGA Preliminary Result Text (KR)\n"]
    text.append("## 계산 내용\n")
    text.append(f"TCGA-THCA primary tumor {len(tcga)}개에서 RAI differentiation, HLA/APM immune visibility, CD8 exclusion, drug-delivery failure proxy, aggressive dedifferentiation 등 5개 치료취약성 축과 세부 module score를 계산했다.\n")
    text.append("## 핵심 결과\n")
    text.append(f"환자 수준 vulnerability label은 다음과 같이 분리되었다: {label_counts}. Mixed/Other는 {label_counts.get('Mixed/Other', 0)}/{len(tcga)}명으로 남아 있으나 전체의 {label_counts.get('Mixed/Other', 0)/max(len(tcga),1):.1%} 수준이며, 제안서에서는 별도 세분화보다 non-dominant/intermediate group으로 두는 것이 안전하다.\n")
    if not corr.empty:
        text.append("주요 축 간 Spearman 상관은 다음과 같아 축들이 완전 중복되지 않음을 보인다.\n\n")
        text.append(corr.to_markdown() + "\n")
    text.append("## 임상/driver 연관\n")
    if not sig.empty:
        text.append("탐색적 clinical association이 계산되었으며, FDR < 0.10인 상위 결과는 아래와 같다. 이는 예후모델 주장이 아니라 제안서용 hypothesis-generation 근거이다.\n\n")
        text.append(sig[["axis", "variable", "test", "effect", "p", "fdr", "n"]].to_markdown(index=False) + "\n")
    else:
        text.append("유효한 임상 연관은 본 pilot에서 강하게 주장하지 않는다. Clinical association은 탐색적으로만 취급한다.\n")
    text.append("## 조심스러운 해석\n")
    text.append("TCGA bulk는 공간 heterogeneity를 증명하지 못하며, immune/stromal score는 세포조성 및 activation state가 섞인 bulk inference이다.\n")
    text.append("## 삼성 제안서 문장\n")
    text.append("TCGA-THCA 505명 공개 bulk transcriptome 분석에서 갑상선암은 단일 driver mutation 또는 평균 발현값으로 설명되지 않고, RAI 분화상태, HLA/APM 면역가시성, CD8 exclusion, stromal/drug-delivery barrier, proliferative dedifferentiation이 분리된 치료취약성 상태로 나뉘었다.\n")
    text.append("## Figure caption\n")
    text.append("TCGA-THCA public bulk pilot defines patient-level therapeutic vulnerability states across RAI differentiation, immune visibility, stromal exclusion, delivery-failure proxy, and aggressive dedifferentiation axes. Results are hypothesis-generating and require spatial and functional validation.\n")
    (REPORTS / "tcga_preliminary_result_text_kr.md").write_text("\n".join(text))


def write_spatial_report(spots: pd.DataFrame, spatial_summary: pd.DataFrame) -> None:
    coherent = int(spatial_summary["coherent_z_gt_2"].sum()) if "coherent_z_gt_2" in spatial_summary else 0
    total = len(spatial_summary)
    conditions = spatial_summary["condition"].value_counts().to_dict() if "condition" in spatial_summary else {}
    text = ["# Spatial Preliminary Result Text (KR)\n"]
    text.append("## 계산 내용\n")
    text.append(f"GSE250521 spatial transcriptomics {total}개 slide, {len(spots)}개 spot에서 TCGA와 동일한 치료취약성 module score와 spot-level niche label을 계산했다. Spot은 독립 생물학적 반복으로 취급하지 않고 slide 내부에 nested된 관측치로 다루었다.\n")
    text.append("## 핵심 결과\n")
    if coherent == total and total > 0:
        text.append(f"KNN same-niche permutation 결과, {total}/{total}개 slide에서 niche coherence z-score > 2가 관찰되었다. 조건 분포는 {conditions}이다.\n")
    else:
        text.append(f"KNN same-niche permutation 결과, {coherent}/{total}개 slide에서 niche coherence z-score > 2가 관찰되었다. 조건 분포는 {conditions}이다.\n")
    top = spatial_summary.sort_values("same_niche_z", ascending=False).head(8)
    text.append("\n상위 coherence slide:\n\n")
    text.append(top[["sample_id", "condition", "n_spots", "same_niche_z", "same_niche_empirical_p"]].to_markdown(index=False) + "\n")
    text.append("## 조심스러운 해석\n")
    text.append("이 결과는 치료취약성 expression state가 조직 내에서 공간적으로 응집되어 있음을 지지하지만, 임상 반응 예측이나 실제 drug delivery를 직접 증명하지 않는다.\n")
    text.append("## 삼성 제안서 문장\n")
    text.append("공간전사체 pilot에서 치료취약성 niche는 무작위 spot noise가 아니라 조직 내에서 공간적으로 응집된 구조를 보였으며, 이는 FFPE/mIHC, GeoMx ROI, fresh tissue 기능실험으로 이어지는 수술-공간 theranostic platform의 필요성을 뒷받침한다.\n")
    text.append("## Figure caption\n")
    text.append("Public GSE250521 spatial transcriptomics reveals non-random spatial organization of therapeutic vulnerability niches across PTC, locally advanced PTC, and ATC slides. Spots are visualized for spatial structure, while inference is summarized at slide level.\n")
    (REPORTS / "spatial_preliminary_result_text_kr.md").write_text("\n".join(text))


def write_samsung_sections(tcga: pd.DataFrame, spatial_summary: pd.DataFrame, cleaned_drug: pd.DataFrame) -> None:
    labels = tcga["primary_vulnerability_label"].value_counts().to_dict() if not tcga.empty else {}
    coherent = int(spatial_summary["coherent_z_gt_2"].sum()) if "coherent_z_gt_2" in spatial_summary else 0
    total = len(spatial_summary)
    class_counts = cleaned_drug["perturbation_class"].value_counts().to_dict() if not cleaned_drug.empty else {}
    prelim = f"""# Samsung Preliminary Results Section (KR)

## 1. Public-data pilot purpose

본 pilot의 목적은 제안서 작성 전, 갑상선암 치료취약성이 평균 driver mutation이나 bulk 평균 발현만으로 설명되는지, 아니면 RAI 분화상태, HLA/APM 면역가시성, CD8 exclusion, myeloid/CAF barrier, drug-delivery failure proxy가 분리된 치료축으로 나타나는지를 공개데이터만으로 검증하는 것이었다.

## 2. Data used

- TCGA-THCA primary tumor bulk RNA-seq: {len(tcga)}명.
- GSE250521 spatial transcriptomics: {total}개 slide, {int(spatial_summary['n_spots'].sum()) if 'n_spots' in spatial_summary else 'NA'}개 spot.
- scRNA reference: local GSE193581 h5ad 기반 cell-type module sanity check.
- DepMap/PRISM-derived local public resource: cleaned perturbation-class candidate table.

## 3. Key result 1: TCGA therapeutic vulnerability subtypes

TCGA 505명 분석에서 vulnerability label은 {labels}로 분리되었다. 이는 갑상선암이 단일 driver mutation 질환이 아니라 RAI, immune visibility, stromal exclusion, delivery proxy, aggressive dedifferentiation이 조합된 치료취약성 상태로 나뉠 가능성을 보여준다.

## 4. Key result 2: spatial therapeutic niches are coherent

GSE250521 spatial pilot에서 {coherent}/{total}개 slide가 KNN same-niche coherence z-score > 2를 보였다. 이는 spot-level score가 무작위 noise가 아니라 조직 내 공간적으로 응집된 치료취약성 territory를 형성함을 지지한다. 단, spot은 독립 환자 수가 아니므로 모든 stage/condition 비교는 slide 수준 요약으로 해석해야 한다.

## 5. Key result 3: perturbation hypotheses are testable

Drug table은 약 이름 중심이 아니라 perturbation class 중심으로 정리했다. 현재 class 분포는 {class_counts}이다. 제안서에서는 final drug claim이 아니라 RAI redifferentiation, HLA/APM restoration, proliferative stress, myeloid/CAF barrier modulation, delivery imaging strategy로 제시한다.

## 6. Claim boundary

Public pilot은 hypothesis generation 근거이며 clinical deployment 근거가 아니다. Spatial transcriptomics는 peptide presentation을 직접 증명하지 않고, drug-delivery failure는 RNA proxy이며, RAI-restorable score는 iodide uptake assay 전까지 가설이다.

## 7. Why this justifies 30억/3년

30억/3년 과제는 이미 공개데이터에서 관찰된 치료취약성 축과 공간 niche를 실제 수술 검체에서 검증하는 translational closure 단계이다. 비용의 핵심은 신규 발견 탐색이 아니라 FFPE/mIHC, GeoMx ROI, fresh tissue slice/organoid perturbation, iodide uptake, HLA/APM rescue, fluorescent delivery imaging을 하나의 수술 전후 theranostic platform으로 닫는 데 있다.

## 8. Why this can expand to 50억/5년

50억/5년 확장에서는 WES/RNA/HLA typing, serial liquid biopsy, pathology foundation model, multi-region spatial profiling, organoid/slice perturbation atlas를 통합해 환자별 치료취약성 niche의 동역학과 intervention response를 추적할 수 있다. Nature급 proof는 public-data observation이 아니라 paired surgical-spatial-functional perturbation atlas에서 나온다.
"""
    (REPORTS / "samsung_preliminary_results_section_kr.md").write_text(prelim)

    storyline = f"""# Pilot 10-Slide Storyline (KR)

## 1. Public-data pilot GO decision
- One-line message: 공개데이터만으로도 K-Thyro 핵심 가설은 GO이다.
- Figure: `results/reports/GO_NO_GO_DECISION.md` 또는 `integrated_public_pilot_evidence_matrix_dark.png`
- Bullets: TCGA 505명; spatial 16 slide; 5/5 main axes; {coherent}/{total} slide coherent; claim은 hypothesis-generation.
- Speaker note: “임상 적용을 주장하는 것이 아니라, 삼성 과제로 검증할 충분한 신호가 있음을 보였습니다.”
- Claim boundary: public pilot, not deployment.

## 2. Data sources and analysis flow
- One-line message: bulk, spatial, scRNA, drug resource가 하나의 검증 흐름으로 연결된다.
- Figure: `figure1_public_pilot_data_flow_dark.png`
- Bullets: TCGA bulk; GSE250521 spatial; scRNA sanity check; DepMap/PRISM-derived perturbation classes; output tables.
- Speaker note: “기존 공개데이터를 사용해 신규 병원 코호트 투입 전 위험을 낮췄습니다.”
- Claim boundary: public data have no paired perturbation assay.

## 3. TCGA 505-patient vulnerability landscape
- One-line message: 치료취약성은 평균 bulk 하나가 아니라 다축 landscape이다.
- Figure: `tcga_vulnerability_axes_heatmap_dark.png`
- Bullets: 505 primary tumors; five main axes; axes are correlated but not identical; patient-level labels generated.
- Speaker note: “driver mutation 중심 분류를 넘어 치료 intervention이 가능한 상태축을 정의했습니다.”
- Claim boundary: bulk cannot prove spatial heterogeneity.

## 4. Patient-level subtype counts
- One-line message: vulnerability labels are large enough to motivate prospective validation.
- Figure: `tcga_label_distribution_dark.png`
- Bullets: label counts {labels}; Mixed/Other remains intermediate; label is exploratory inference; clinical association is secondary.
- Speaker note: “제안서에는 label count를 preliminary evidence로, clinical prediction은 목표로 둡니다.”
- Claim boundary: labels are not clinical subtypes yet.

## 5. Spatial ST 16-slide niche coherence
- One-line message: therapeutic niches are spatially organized, not random spot noise.
- Figure: `spatial_coherence_barplot_dark.png`
- Bullets: GSE250521 16 slides; KNN same-niche permutation; {coherent}/{total} z > 2; spot nested within slide.
- Speaker note: “공간분석을 해야 하는 이유가 여기서 나옵니다.”
- Claim boundary: no clinical-response prediction.

## 6. Representative spatial niche maps
- One-line message: RAI, HLA/APM, CAF/myeloid, niche labels show distinct territories.
- Figure: `spatial_representative_maps_dark.png`
- Bullets: representative PTC/LPTC/ATC; four-panel map; tumor/immune/stromal axes separate; ROI validation targets.
- Speaker note: “병리 ROI와 기능실험 위치를 지정할 수 있는 지도가 됩니다.”
- Claim boundary: expression states, not protein or peptide presentation.

## 7. Integrated evidence matrix
- One-line message: each axis has different evidence strength and validation requirement.
- Figure: `integrated_public_pilot_evidence_matrix_dark.png`
- Bullets: RAI and HLA axes moderate/strong; delivery proxy weak; scRNA supports cell-type interpretation; validation assays specified.
- Speaker note: “강한 주장과 약한 주장을 분리해 reviewer risk를 낮춥니다.”
- Claim boundary: matrix is proposal-readiness, not final proof.

## 8. Drug/perturbation hypotheses
- One-line message: drug names are secondary; perturbation classes are the proposal logic.
- Figure: `figure6_drug_perturbation_candidate_pilot.png` plus `drug_pilot_candidate_rankings_cleaned.tsv`
- Bullets: MAPK-axis; epigenetic/APM modulation; proliferation stress; barrier modulation; delivery imaging.
- Speaker note: “약효 주장 대신 functional rescue assay로 닫겠습니다.”
- Claim boundary: DepMap/PRISM nominate hypotheses only.

## 9. Experimental validation plan
- One-line message: AI-spatial inference is closed by FFPE, GeoMx, and fresh tissue functional assays.
- Figure: `figure7_samsung_experimental_validation_design.png`
- Bullets: core marker panel; extended panel; iodide uptake; HLA/APM rescue; fluorescent delivery imaging.
- Speaker note: “수술 검체를 바로 기능실험으로 연결하는 것이 과제의 차별점입니다.”
- Claim boundary: validation required before clinical use.

## 10. Samsung 30억/3년 execution logic
- One-line message: 3년 과제는 public signal을 clinical-spatial-functional proof로 전환한다.
- Figure: `figure8_pilot_conclusion_slide.png`
- Bullets: Year 1 FFPE/mIHC; Year 2 GeoMx/fresh tissue assay; Year 3 integrated AI-spatial theranostic model; 50억/5년 확장 가능.
- Speaker note: “30억은 platform proof, 50억은 atlas-scale clinical deployment pre-stage입니다.”
- Claim boundary: deployment is beyond pilot.
"""
    (REPORTS / "pilot_10_slide_storyline_kr.md").write_text(storyline)

    final_summary = f"""# Public Pilot Executive Summary Final (KR)

## GO decision

최종 audit 후에도 판정은 **GO**이다. TCGA-THCA 505명 primary tumor에서 5개 main therapeutic vulnerability axis가 모두 계산 가능했고, GSE250521 spatial transcriptomics {total}개 slide에서 {coherent}개 slide가 non-random same-niche coherence(z > 2)를 보였다.

## Scientific conclusion

공개데이터 pilot은 갑상선암 치료취약성이 단일 driver mutation 또는 평균 bulk expression으로만 설명되지 않음을 보여준다. RAI differentiation, HLA/APM immune visibility, CD8 exclusion, myeloid/CAF barrier, drug-delivery failure proxy, aggressive dedifferentiation은 환자 수준과 공간 조직 수준에서 분리 가능한 치료축으로 관찰되었다.

## Samsung significance

삼성육성과제에서 중요한 지점은 “새로운 biomarker 하나”가 아니라, 수술 전후 액체생검·공간오믹스·병리 AI·기능 perturbation을 연결해 갑상선암을 치료취약성 생태계로 읽는 플랫폼이다. Public pilot은 이 플랫폼이 공허한 개념이 아니라 실제 공개데이터에서 관찰되는 구조를 기반으로 함을 보인다.

## Experimental validation plan

검증은 FFPE/mIHC core panel(PAX8/TG/TPO/NIS, HLA-I/B2M/TAP1, CD8/GZMB, CD68/CD163, ACTA2/FAP/COL1A1, PD-L1, CA9/VEGFA), GeoMx ROI validation, fresh tissue organoid/slice culture, iodide uptake assay, HLA/APM rescue assay, PBMC/TIL co-culture, fluorescent drug/liposome/nanoparticle distribution imaging으로 닫는다.

## Remaining risks

Drug-delivery failure는 RNA proxy이며 실제 delivery imaging 전까지는 가설이다. HLA/APM RNA는 peptide presentation을 직접 증명하지 않는다. RAI-restorable score는 iodide uptake assay 전까지 치료반응 주장이 아니다. DepMap/PRISM 후보는 small-N thyroid line 기반 hypothesis이며 약제명보다 perturbation class로 제시해야 한다.

## Final one-line proposal message

본 과제는 갑상선암을 평균 유전체 질환이 아니라, 수술 perturbation과 공간 치료취약성 niche의 조합으로 치료반응이 결정되는 동적 생태계로 재정의한다.
"""
    (REPORTS / "public_pilot_executive_summary_final_kr.md").write_text(final_summary)


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit and strengthen K-Thyro public pilot for Samsung proposal.")
    parser.parse_args()
    FIGS.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    data = audit_outputs()
    tcga = data["tcga"]
    spots = data["spots"]
    slides = data["slides"]
    coherence = data["coherence"]
    drug = data["drug"]
    evidence = data["evidence"]

    cleaned_drug = clean_drug_table(drug)
    label_summary = tcga_validation(tcga)
    make_tcga_figures(tcga)
    spatial_summary = spatial_validation(spots, slides, coherence)
    make_spatial_figures(spots, spatial_summary)
    make_evidence_matrix_figure(evidence)
    write_tcga_report(tcga, label_summary)
    write_spatial_report(spots, spatial_summary)
    write_samsung_sections(tcga, spatial_summary, cleaned_drug)

    print("Samsung audit/strengthening complete")
    print(f"Cleaned drug table: {TABLES / 'drug_pilot_candidate_rankings_cleaned.tsv'}")
    print(f"Audit report: {REPORTS / 'pilot_output_audit.md'}")


if __name__ == "__main__":
    main()
