#!/usr/bin/env python3
"""Build cold method comparison with training-data provenance."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
CUR = ROOT / "curation_2026_05_09"
OUT = ROOT / "method_training_audit_2026_05_09"
STRICT = ROOT / "wave8b_strict_structure_2026_05_09" / "strict_esmfold_method_comparison.tsv"

LOCAL_TRAIN = "local train pool n=2,396: CEDAR 909; TESLA_mmc4 605; NEPdb 572; TESLA_mmc7_validation 310"
LOCAL_TEST = "ITSNdb MHC-I immunogenicity: combined n=319; no-overlap n=106; strict no TCR/self exact + HLA pseudo n=89"


def info(
    training_data: str,
    source: str,
    overlap_risk: str,
    disposition: str,
    note: str,
    urls: str = "",
    confidence: str = "local_verified",
) -> dict:
    return {
        "training_data_summary": training_data,
        "training_data_source": source,
        "training_overlap_risk": overlap_risk,
        "paper_disposition": disposition,
        "provenance_confidence": confidence,
        "training_data_note": note,
        "source_urls": urls,
    }


PUBLIC = {
    "MHCflurry": info(
        "public pretrained MHC-I presentation model; local install mhcflurry 2.2.1 model generated 2020-06-11; pretrained on MS-identified MHC-I ligands and peptide/MHC affinity measurements from IEDB plus other sources",
        "official MHCflurry docs + local model files",
        "unresolved/high: ITSNdb or close public IEDB/MS-derived measurements could overlap the pretrained corpus",
        "caveated comparator only",
        "use as public-pretrained upper-bound comparator, not clean external baseline",
        "https://openvax.github.io/mhcflurry/intro.html",
        "official_doc_plus_local",
    ),
    "NetMHCpan_4.1": info(
        "NetMHCpan 4.1 public ANN; >850,000 quantitative binding-affinity and MS-eluted-ligand peptides across MHC-I molecules",
        "DTU NetMHCpan 4.1 service docs",
        "unresolved/high: public BA/EL training; not immunogenicity-specific; overlap with IEDB-like benchmark cannot be excluded without train data audit",
        "caveated binder comparator",
        "binding/presentation proxy, not immunogenicity model",
        "https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/",
        "official_doc",
    ),
    "NetMHCstabpan": info(
        "NetMHCstabpan public ANN; >25,000 quantitative peptide-MHC-I stability measurements covering 75 HLA molecules",
        "DTU NetMHCstabpan service docs",
        "medium/high: public stability training; biological endpoint mismatched to immunogenicity",
        "negative/supplement",
        "below chance on no-overlap; keep as endpoint-mismatch control",
        "https://services.healthtech.dtu.dk/services/NetMHCstabpan-1.0/",
        "official_doc",
    ),
    "MHCnuggets_2": info(
        "public MHC class I/II binding predictor; paper describes IEDB-curated binding data and allele-specific transfer learning",
        "MHCnuggets paper / IEDB-derived training",
        "unresolved/high: IEDB-derived training; endpoint is binding, not direct immunogenicity",
        "negative/supplement",
        "near chance on no-overlap in this benchmark",
        "https://www.cs.jhu.edu/~rohit/research/mhcnuggets_paper.pdf",
        "paper_reported",
    ),
    "BigMHC_IM": info(
        "BigMHC public model; Mendeley release includes el_train/el_val/el_test for presentation and im_train/im_val/im_test for immunogenicity transfer learning, plus IEDB and MANAFEST files",
        "BigMHC Mendeley dataset + paper",
        "unresolved/high: explicit immunogenicity-transfer training and public neoepitope/IEDB resources; strong in-master inflation locally",
        "caveated comparator only",
        "good in-master, weaker no-overlap; not clean primary evidence",
        "https://data.mendeley.com/datasets/dvmz6pkzvb/4; https://www.nature.com/articles/s42256-023-00694-6",
        "official_dataset",
    ),
    "PRIME": info(
        "PRIME public immunogenicity predictor; reported curated 4,958 peptide dataset with pathogen/cancer-testis and cancer mutation immunogenicity labels, plus antigen-presentation/TCR-recognition features",
        "PRIME Cell Reports Medicine paper",
        "unresolved/high: direct immunogenicity public training; benchmark overlap cannot be excluded without row-level audit",
        "caveated comparator",
        "flat but modest locally; keep as public immunogenicity comparator",
        "https://pubmed.ncbi.nlm.nih.gov/33665637/; https://pmc.ncbi.nlm.nih.gov/articles/PMC7897774/",
        "paper_reported",
    ),
    "DeepImmuno": info(
        "DeepImmuno-CNN public immunogenicity model; >9,000 IEDB T-cell assay records filtered to 8,971 peptide-HLA instances from IEDB as of 2020-08-13; evaluated on TESLA/COVID independent sets",
        "DeepImmuno Briefings in Bioinformatics paper",
        "unresolved/high: IEDB/TESLA-style immunogenicity training/evaluation; local in-master inflation",
        "caveated comparator",
        "not clean external; inflation gap makes it supplement-only",
        "https://academic.oup.com/bib/article/22/6/bbab160/6261914",
        "paper_reported",
    ),
    "TransPHLA": info(
        "TransPHLA public pHLA binding transformer; positive pHLA binders from Anthem plus IEDB, EPIMHC, MHCBN, SYFPEITHI, MS HLA ligands, and training sets of other pHLA tools; negatives sampled from IEDB immunopeptidome source proteins",
        "TransPHLA Nature Machine Intelligence paper",
        "unresolved/high: broad public binder/MS-ligand training; binding endpoint only",
        "caveated binder comparator",
        "not direct immunogenicity; weak no-overlap here",
        "https://www.nature.com/articles/s42256-022-00459-7; https://xlab.sjtu.edu.cn/pdf/2022-NMI-TransMut.pdf",
        "paper_reported",
    ),
    "TSCAPE_TITANiAN": info(
        "T-SCAPE/TITANiAN public multidomain immunogenicity model; reported pretraining on MHC presentation, pMHC binding, TCR-pMHC interaction, source organism, T-cell activation; sources include IEDB, VDJdb, McPAS-TCR, ImmuneCODE, TBAdb, 10x Genomics, OAS, UniProt, PRIME, MHCBN, BigMHC, BioPhi, PanPep",
        "T-SCAPE Science Advances / PMC data availability",
        "unresolved/high: many public immunogenicity/TCR/presentation sources; exact training overlap must be audited",
        "caveated comparator / future benchmark target",
        "modern but still not clean until row-level overlap audit",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC12680054/; https://pubmed.ncbi.nlm.nih.gov/41348893/",
        "paper_reported",
    ),
}


def local_method(method: str, wave: str, category: str) -> dict:
    if method.startswith("Wave8"):
        if "TCR" in method and "SelfSim" in method:
            td = f"{LOCAL_TRAIN}; logistic regression on Wave8 features. Features use VDJdb + McPAS-TCR reference peptide-HLA sets and UniProt human SwissProt self-similarity."
        elif "TCR" in method:
            td = f"{LOCAL_TRAIN}; logistic regression on TCR-reference features from local VDJdb + McPAS-TCR peptide-HLA reference sets."
        else:
            td = f"{LOCAL_TRAIN}; logistic regression on UniProt human SwissProt self-similarity features."
        return info(
            td,
            "local wave8 scripts and WAVE8_REPORT",
            "medium/high: no benchmark-label training leakage for no-overlap, but exact TCR/self reference hits create artifact risk",
            "candidate retest only",
            "headline no-overlap gain collapses after exact-reference stress; use only with exact/near-match flags",
            "local: wave8/track8a_tcr_features.py; wave8/track8c_self_similarity.py",
        )
    if method == "Structure_LR":
        return info(
            f"{LOCAL_TRAIN}; local logistic regression on CPU biophysical proxies: HLA-specific anchor PWM from train positives, BLOSUM similarity to train binders, hydrophobicity/charge/helix/sheet/length/aromatic features.",
            "local wave2 compute_structure_features.py",
            "low/medium: known local training split; still trained on public CEDAR/TESLA/NEPdb labels but no hidden pretrained corpus",
            "primary clean local baseline",
            "strict-set leader among locked non-CV scores",
            "local: wave2/compute_structure_features.py",
        )
    if method in {"VQC", "GP_quantum", "W7A_QK_only", "W7A_full"}:
        return info(
            f"{LOCAL_TRAIN}; quantum/VQC or quantum-kernel models using PCA/structure-derived local features; some runs trained on capped train subsets for speed.",
            "local wave2/wave7 scripts",
            "low/medium: known local train pool; internal tuning/model-selection risk",
            "quantum branch candidate / not primary claim yet",
            "keep as active branch; requires locked external validation before paper claim",
            "local: wave2/train_vqc.py; wave2/train_quantum_kernel_svm.py; wave7/train_w7.py",
        )
    if method.startswith("Stack_") or method == "W7B_stacked":
        risk = "very high" if "inmaster" in method.lower() else "high"
        disp = "supplement only" if "inmaster" in method.lower() else "negative/supplement"
        return info(
            f"local stack/ensemble over existing model predictions; {LOCAL_TRAIN}; some variants explicitly fit on ITSNdb in-master rows.",
            "local wave4c/wave7 ensemble artifacts",
            f"{risk}: stacking inflates in-master and is not fair as a primary external method",
            disp,
            "use to demonstrate why naive ensembling is unsafe",
            "local: wave4c/ensemble_predictions.tsv; wave7/predictions_combined.tsv",
        )
    if method.startswith("MultiTask_"):
        return info(
            f"{LOCAL_TRAIN}; local deep multitask variants built from ESM/peptide-HLA features and auxiliary supervision.",
            "local wave5b scripts",
            "medium/high: strong in-master performance, below-chance no-overlap; likely dataset/source shortcut",
            "negative ablation",
            "do not headline; useful as failure mode",
            "local: wave5b/train_multitask.py",
        )
    if method in {"ESM2_Bayesian", "RF_biophys"}:
        if method == "ESM2_Bayesian":
            return info(
                f"{LOCAL_TRAIN}; ESM2-150M frozen peptide/HLA-pseudo embeddings plus Bayesian MLP/deep ensemble trained on bundle train split.",
                "local wave1 train_bayesian_neo.py",
                "medium/high: huge in-master AUROC but no-overlap collapse; source/overlap shortcut baseline",
                "negative headline baseline",
                "canonical failure baseline for leakage/generalization story",
                "local: train_bayesian_neo.py",
            )
        return info(
            f"{LOCAL_TRAIN}; local random forest/biophysical baseline from prior robustness audit.",
            "local wave11 / robustness audit",
            "medium/high: in-master inflation and external collapse after overlap audit",
            "negative/supplement",
            "do not use as clean primary baseline",
            "local: wave11 artifacts",
        )
    if method in {"kNN", "TENT"}:
        if method == "kNN":
            return info(
                f"{LOCAL_TRAIN}; retrieval over local ESM2 peptide embeddings, per-HLA train buckets with global fallback.",
                "local wave4b m5_knn.py",
                "high: retrieval/memorization behavior; in-master inflated",
                "negative/supplement",
                "use to show retrieval without leakage control is unsafe",
                "local: wave4b/m5_knn.py",
            )
        return info(
            f"{LOCAL_TRAIN}; surrogate ESM2 head plus test-time entropy minimization on ITSNdb test rows.",
            "local wave4b m6_tent.py",
            "high: adaptation directly on test distribution; no-overlap below chance",
            "negative/supplement",
            "not paper-primary",
            "local: wave4b/m6_tent.py",
        )
    if method.startswith("NetMHCIIpan") or "class-II" in method or "HLA-D" in method:
        return info(
            "not run locally for peptide immunogenicity endpoint; NetMHCIIpan 4.3 itself is trained on >650,000 BA/EL measurements across HLA-DR/DQ/DP, H-2, BoLA-DRB3",
            "local class-II status + DTU NetMHCIIpan docs",
            "not comparable to current ITSNdb MHC-I endpoint",
            "not evaluated",
            "requires separate class-II benchmark",
            "https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/",
            "official_doc",
        )
    return info(
        f"{LOCAL_TRAIN}; local model/ablation derived from bundle train pool unless marked public-pretrained.",
        "local artifacts",
        "medium: known local train split but no dedicated external/time split",
        "negative/supplement",
        "not a primary winner under no-overlap audit",
        "local artifacts",
    )


def main() -> None:
    OUT.mkdir(exist_ok=True)
    comp = pd.read_csv(CUR / "class1_class2_all_method_comparison.tsv", sep="\t")
    rows = []
    for _, r in comp.iterrows():
        method = str(r["method"])
        meta = PUBLIC.get(method, local_method(method, str(r["wave"]), str(r["category"])))
        rows.append({**r.to_dict(), **meta})
    out = pd.DataFrame(rows)

    if STRICT.exists():
        st = pd.read_csv(STRICT, sep="\t")
        st = st[["method", "auroc", "auprc", "n", "n_pos"]].rename(
            columns={
                "auroc": "strict_no_tcr_self_AUROC",
                "auprc": "strict_no_tcr_self_AUPRC",
                "n": "strict_n",
                "n_pos": "strict_n_pos",
            }
        )
        out = out.merge(st, on="method", how="left")

    out.to_csv(OUT / "method_training_audit.tsv", sep="\t", index=False)

    # Compact report.
    mhc1 = out[out["class"] == "MHC-I"].copy()
    cols = [
        "method",
        "no_overlap_AUROC",
        "no_overlap_AUPRC",
        "strict_no_tcr_self_AUROC",
        "in_master_AUROC",
        "training_data_summary",
        "training_overlap_risk",
        "paper_disposition",
    ]
    top = mhc1[cols].copy()
    top["no_overlap_AUROC"] = top["no_overlap_AUROC"].map(lambda x: "" if pd.isna(x) else f"{x:.3f}")
    top["no_overlap_AUPRC"] = top["no_overlap_AUPRC"].map(lambda x: "" if pd.isna(x) else f"{x:.3f}")
    top["strict_no_tcr_self_AUROC"] = top["strict_no_tcr_self_AUROC"].map(lambda x: "" if pd.isna(x) else f"{x:.3f}")
    top["in_master_AUROC"] = top["in_master_AUROC"].map(lambda x: "" if pd.isna(x) else f"{x:.3f}")

    lines = [
        "# Cold Method And Training-Data Audit - 2026-05-09",
        "",
        f"Current endpoint: {LOCAL_TEST}.",
        "",
        "Key rule: `no_overlap` only excludes overlap with the in-house master table. It does not prove absence from each public tool's own training corpus.",
        "",
        "## Verdict",
        "",
        "| role | method | reason |",
        "|---|---|---|",
        "| clean local anchor | Structure_LR | known local train split; strict no-reference leader among locked scores |",
        "| apparent but not clean | Wave8_TCR_SelfSim_full | best full no-overlap AUROC, but exact-reference sensitive |",
        "| public upper-bound comparator | MHCflurry / BigMHC / PRIME / DeepImmuno / T-SCAPE | training corpora unresolved or explicitly public immunogenicity/presentation data |",
        "| negative failure baseline | ESM2_Bayesian / kNN / multitask / RF_biophys | strong in-master or combined scores but weak no-overlap |",
        "| active research branch | quantum-kernel + structure | promising strict internal pilot, but not external paper evidence yet |",
        "| not comparable | class II rows | no local class-II peptide immunogenicity endpoint yet |",
        "",
        "## MHC-I Method Table",
        "",
        top.to_markdown(index=False),
        "",
        "## Training-Data Sources Checked",
        "",
        "- Local bundle: `bundle_summary.json` gives train n=2,396 from CEDAR, TESLA_mmc4, NEPdb, TESLA_mmc7_validation.",
        "- MHCflurry docs: pretrained models use MS-identified MHC-I ligands and IEDB peptide/MHC affinity measurements plus other sources.",
        "- NetMHCpan 4.1 docs: trained on >850,000 BA and MS-eluted-ligand peptides.",
        "- NetMHCpan 4.2 docs: newer version adds IEDB pathogen epitope and CEDAR neoepitope fine-tuning; do not use it as a clean comparator without a time/overlap audit.",
        "- NetMHCIIpan 4.3 docs: class-II predictor trained on >650,000 BA/EL measurements; not run locally for comparable immunogenicity AUROC.",
        "- BigMHC Mendeley data: contains `el_train/el_val/el_test` and `im_train/im_val/im_test`, IEDB, MANAFEST, and pseudosequence files.",
        "- DeepImmuno paper: IEDB 2020 T-cell assay training set, 8,971 retained peptide-HLA instances, TESLA used as independent test.",
        "- TransPHLA paper: Anthem plus IEDB/EPIMHC/MHCBN/SYFPEITHI, MS HLA ligands, and other pHLA tool training datasets.",
        "- T-SCAPE paper: multidomain public sources including IEDB, VDJdb, McPAS-TCR, ImmuneCODE, TBAdb, 10x, OAS, UniProt, PRIME, MHCBN, BigMHC, BioPhi, PanPep.",
        "",
        "## Source URLs",
        "",
        "- MHCflurry: https://openvax.github.io/mhcflurry/intro.html",
        "- NetMHCpan 4.1: https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/",
        "- NetMHCpan 4.2: https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/",
        "- NetMHCIIpan 4.3: https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/",
        "- NetMHCstabpan: https://services.healthtech.dtu.dk/services/NetMHCstabpan-1.0/",
        "- BigMHC data: https://data.mendeley.com/datasets/dvmz6pkzvb/4",
        "- DeepImmuno: https://academic.oup.com/bib/article/22/6/bbab160/6261914",
        "- PRIME: https://pubmed.ncbi.nlm.nih.gov/33665637/",
        "- TransPHLA: https://www.nature.com/articles/s42256-022-00459-7",
        "- T-SCAPE: https://pmc.ncbi.nlm.nih.gov/articles/PMC12680054/",
    ]
    (OUT / "METHOD_TRAINING_AUDIT.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
