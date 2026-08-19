#!/usr/bin/env python3
"""Build a literature-mined RL/GA blueprint for a new CROSS-Neo algorithm.

The goal is not to copy any proprietary "secret sauce". Public papers and
publicly runnable tools are decomposed into method primitives: preprocessing,
modalities, encoders, fusion rules, losses, search strategies and guardrails.
Those primitives define the search space for a new knowledge-graph-guided
GA/RL algorithm that can be locked before wetlab validation.
"""

from __future__ import annotations

import html
import json
import shutil
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "project/results/cross_neo_method_mining_rl_ga_blueprint_2026_05_10"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")


METHOD_CARDS = [
    {
        "method": "BigMHC",
        "domain": "neoantigen_immunogenicity",
        "source_type": "biomedical_ai",
        "source_url": "https://doi.org/10.1038/s42256-023-00694-6",
        "paper_anchor": "Nature Machine Intelligence 2023",
        "core_primitives": "MS_EL_pretraining;immunogenicity_transfer_learning;seven_model_ensemble;pan_allelic_peptide_HLA_encoder",
        "preprocessing_primitives": "class_I_filter;8_15mer_peptide_filter;HLA_normalization;EL_random_negative_design",
        "modalities": "peptide_sequence;HLA_allele;MS_eluted_ligand;immune_response_label",
        "search_gene": "two_stage_EL_to_IM_transfer;public_prior_score;ensemble_diversity_weight",
        "risk_guard": "training_overlap_audit;public_tool_dependency_boundary",
    },
    {
        "method": "ImmunoStruct",
        "domain": "neoantigen_immunogenicity",
        "source_type": "biomedical_ai",
        "source_url": "https://www.nature.com/articles/s42256-025-01163-y",
        "paper_anchor": "Nature Machine Intelligence 2025",
        "core_primitives": "multimodal_sequence_structure_biochemistry;interpretable_pMHC_features;multi_allele_class_I_prediction",
        "preprocessing_primitives": "pMHC_structure_standardization;biochemical_descriptor_generation;allele_stratified_split",
        "modalities": "peptide_sequence;HLA_allele;pMHC_structure;biochemical_features;immune_response_label",
        "search_gene": "structure_branch;biochemistry_branch;late_fusion_attention;structure_uncertainty_gate",
        "risk_guard": "structure_availability_bias;allele_coverage_check",
    },
    {
        "method": "DeepImmuno",
        "domain": "neoantigen_immunogenicity",
        "source_type": "biomedical_ai",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7781330/",
        "paper_anchor": "Frontiers Immunology 2021",
        "core_primitives": "CNN_immunogenicity_model;physiochemical_aware_encoding;top_k_sensitivity_benchmark",
        "preprocessing_primitives": "beta_binomial_response_confidence;HLA_peptide_pair_encoding;low_confidence_label_filter",
        "modalities": "peptide_sequence;HLA_allele;amino_acid_properties;immune_response_label",
        "search_gene": "AA_property_channel;confidence_weighted_label;small_CNN_branch;top_k_reward",
        "risk_guard": "small_data_overfit_check;threshold_sensitivity_audit",
    },
    {
        "method": "PRIME",
        "domain": "neoantigen_immunogenicity",
        "source_type": "biomedical_ai",
        "source_url": "https://www.sciencedirect.com/science/article/pii/S2666379121000057",
        "paper_anchor": "Cell Reports Medicine 2021",
        "core_primitives": "presentation_plus_TCR_recognition_propensity;TCR_recognition_determinants;immunoediting_signal",
        "preprocessing_primitives": "TCR_facing_residue_features;presentation_rank_features;neoepitope_label_harmonization",
        "modalities": "peptide_sequence;HLA_presentation_score;TCR_recognition_features;immune_response_label",
        "search_gene": "TCR_propensity_branch;TCR_facing_position_mask;immunoediting_prior",
        "risk_guard": "non_tumor_training_boundary;label_context_check",
    },
    {
        "method": "NetMHCpan_4_1",
        "domain": "antigen_presentation",
        "source_type": "biomedical_ai",
        "source_url": "https://doi.org/10.1093/nar/gkaa379",
        "paper_anchor": "Nucleic Acids Research 2020",
        "core_primitives": "motif_deconvolution;binding_affinity_plus_EL_integration;pan_specific_MHC_prediction",
        "preprocessing_primitives": "allele_resolution_normalization;peptide_length_windows;MS_EL_and_BA_label_merge",
        "modalities": "peptide_sequence;HLA_allele;binding_affinity;MS_eluted_ligand",
        "search_gene": "presentation_gate;motif_deconvolution_prior;BA_EL_dual_score",
        "risk_guard": "presentation_not_immunogenicity_boundary",
    },
    {
        "method": "pMTnet",
        "domain": "TCR_pMHC_recognition",
        "source_type": "biomedical_ai",
        "source_url": "https://www.nature.com/articles/s42256-021-00383-2",
        "paper_anchor": "Nature Machine Intelligence 2021",
        "core_primitives": "transfer_learning_TCR_pMHC;CDR3beta_peptide_HLA_input;pairing_specificity_prediction",
        "preprocessing_primitives": "CDR3beta_normalization;peptide_HLA_pairing;TCR_pair_availability_gate",
        "modalities": "TCR_CDR3beta;peptide_sequence;HLA_allele;pMHC_context",
        "search_gene": "TCR_expert_branch;TCR_missingness_policy;pairing_affinity_gate",
        "risk_guard": "TCR_sparse_proxy_guard;HLA_A02_bias_check",
    },
    {
        "method": "Harmonized_Neoantigen_ML",
        "domain": "patient_level_neoantigen_selection",
        "source_type": "biomedical_dataset_method",
        "source_url": "https://doi.org/10.1016/j.immuni.2023.09.002",
        "paper_anchor": "Immunity 2023",
        "core_primitives": "WES_RNA_reprocessing;presentation_hotspots;binding_promiscuity;oncogenicity_context",
        "preprocessing_primitives": "matched_WES_RNA_harmonization;patient_cohort_split;SNV_to_neopeptide_expansion;expression_filter",
        "modalities": "DNA_variant;RNA_expression;peptide_sequence;HLA_allele;gene_context;patient_context",
        "search_gene": "expression_gate;hotspot_feature;binding_promiscuity_feature;oncogenicity_prior",
        "risk_guard": "patient_level_leakage_guard;cohort_generalization_split",
    },
    {
        "method": "CEDAR",
        "domain": "curated_cancer_epitope_data",
        "source_type": "biomedical_database",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9825495/",
        "paper_anchor": "Nucleic Acids Research 2023",
        "core_primitives": "curated_cancer_epitope_labels;antigen_subtype_taxonomy;assay_metadata",
        "preprocessing_primitives": "neoantigen_viral_self_other_taxonomy;assay_type_filter;full_HLA_resolution_filter",
        "modalities": "peptide_sequence;HLA_allele;assay_metadata;antigen_subtype;response_label",
        "search_gene": "label_provenance_node;assay_confidence_weight;antigen_subtype_gate",
        "risk_guard": "publication_overlap_audit;assay_context_stratification",
    },
    {
        "method": "AutoML_Zero",
        "domain": "ai_algorithm_discovery",
        "source_type": "ai_method",
        "source_url": "https://arxiv.org/abs/2003.03384",
        "paper_anchor": "ICML 2020",
        "core_primitives": "evolve_learning_algorithms_from_math_ops;low_human_bias_search;emergent_regularization",
        "preprocessing_primitives": "primitive_operation_library;task_suite_curriculum;program_safety_checks",
        "modalities": "algorithm_code;math_operations;validation_tasks",
        "search_gene": "program_tree_gene;loss_function_gene;optimizer_gene;regularizer_gene",
        "risk_guard": "compute_budget_cap;invalid_program_rejection;complexity_penalty",
    },
    {
        "method": "NAS_RL",
        "domain": "ai_architecture_search",
        "source_type": "ai_method",
        "source_url": "https://research.google/pubs/neural-architecture-search-with-reinforcement-learning/",
        "paper_anchor": "ICLR 2017",
        "core_primitives": "RNN_controller_generates_architecture;validation_reward;policy_gradient_search",
        "preprocessing_primitives": "architecture_tokenization;child_model_training_budget;reward_normalization",
        "modalities": "architecture_tokens;validation_metric;training_trace",
        "search_gene": "RL_controller_action;mutation_policy_learning;reward_shaping",
        "risk_guard": "validation_overfit_guard;early_stop_fidelity_check",
    },
    {
        "method": "Regularized_Evolution",
        "domain": "ai_architecture_search",
        "source_type": "ai_method",
        "source_url": "https://arxiv.org/abs/1802.01548",
        "paper_anchor": "AAAI 2019",
        "core_primitives": "aging_evolution;tournament_selection;architecture_mutation;limited_compute_efficiency",
        "preprocessing_primitives": "population_initialization;age_tracking;fitness_cache",
        "modalities": "genotype;fitness_metric;age_metadata",
        "search_gene": "aging_selection;diversity_preserving_mutation;fitness_cache_reuse",
        "risk_guard": "random_search_baseline;search_cost_reporting",
    },
    {
        "method": "Deep_Distilling",
        "domain": "ai_algorithm_discovery",
        "source_type": "ai_method",
        "source_url": "https://www.nature.com/articles/s43588-024-00593-9",
        "paper_anchor": "Nature Computational Science 2024",
        "core_primitives": "distill_neural_solution_to_code;symbolic_essence_network;human_comprehensible_algorithm",
        "preprocessing_primitives": "task_trace_generation;symbolic_program_extraction;OOD_generalization_test",
        "modalities": "training_data;neural_parameters;distilled_code;OOD_task_suite",
        "search_gene": "distill_to_interpretable_rule;posthoc_code_simplification;OOD_reward",
        "risk_guard": "equivalence_test;overcompression_failure_check",
    },
]


ACTION_CARDS = [
    {
        "action_family": "preprocessing",
        "action": "choose_label_confidence_model",
        "choices": "binary;beta_binomial;assay_weighted;patient_weighted",
        "borrowed_from": "DeepImmuno;CEDAR",
        "rl_state_dependency": "assay metadata completeness, replicate count, source label noise",
        "guardrail": "no label-confidence feature may encode source-only positives",
    },
    {
        "action_family": "preprocessing",
        "action": "choose_split_policy",
        "choices": "source_holdout;patient_holdout;HLA_holdout;peptide_cluster_holdout;low_leakage_holdout",
        "borrowed_from": "Harmonized_Neoantigen_ML;AutoML_review",
        "rl_state_dependency": "source imbalance, HLA skew, duplicate peptide/HLA pairs",
        "guardrail": "all score claims report source and leakage split",
    },
    {
        "action_family": "modality",
        "action": "activate_modality_branches",
        "choices": "sequence;HLA;MS_EL;RNA_expression;TCR;structure;biochemistry;gene_context;assay_metadata",
        "borrowed_from": "BigMHC;ImmunoStruct;pMTnet;Harmonized_Neoantigen_ML",
        "rl_state_dependency": "candidate has modality available, missingness pattern, wetlab budget",
        "guardrail": "missing modality branch must have explicit missingness policy",
    },
    {
        "action_family": "encoder",
        "action": "choose_encoder_family",
        "choices": "tabular_GBDT;small_CNN;LSTM;transformer;protein_LM_embedding;structure_GNN;late_fusion_MLP",
        "borrowed_from": "BigMHC;DeepImmuno;ImmunoStruct;pMTnet",
        "rl_state_dependency": "sample size, modality count, HLA diversity, compute budget",
        "guardrail": "deep encoder requires nested validation and calibration check",
    },
    {
        "action_family": "fusion",
        "action": "choose_fusion_rule",
        "choices": "weighted_sum;mixture_of_experts;late_attention;rank_aggregation;claim_safe_cap;Pareto_front",
        "borrowed_from": "BigMHC ensemble;CROSS-Neo;AutoML",
        "rl_state_dependency": "comparator disagreement, failure mode, claim layer",
        "guardrail": "single-concept dominance cap",
    },
    {
        "action_family": "search",
        "action": "choose_search_operator",
        "choices": "aging_evolution;mutation;cross_over;RL_policy_mutation;Bayesian_local_search;random_baseline",
        "borrowed_from": "Regularized_Evolution;NAS_RL;AutoML_Zero",
        "rl_state_dependency": "search stagnation, diversity, compute budget",
        "guardrail": "random search and fixed ensemble baselines always reported",
    },
    {
        "action_family": "reward",
        "action": "compose_reward",
        "choices": "AUPRC;source_balanced_AUPRC;low_leakage_AUPRC;top96_precision;wetlab_hit_yield;complexity_penalty;novelty_bonus",
        "borrowed_from": "AutoML_review;CROSS-Neo v6;Deep_Distilling",
        "rl_state_dependency": "train/validation split, assay objective, product claim boundary",
        "guardrail": "reward must include at least one leakage-resistant metric",
    },
    {
        "action_family": "wetlab",
        "action": "allocate_96_well_arms",
        "choices": "clean_discovery;TCR_MD_mechanism;label_rescue;specificity_moat;positive_QC;model_boundary",
        "borrowed_from": "CROSS-Neo v6 interpreter",
        "rl_state_dependency": "model uncertainty, comparator disagreement, claim blocker",
        "guardrail": "every claimed hit needs matched control route",
    },
]


REWARD_COMPONENTS = [
    ("source_balanced_AUPRC", 0.24, "prevents CEDAR/TESLA/source dominance"),
    ("low_leakage_AUPRC", 0.18, "rewards clean generalization"),
    ("top96_precision", 0.18, "matches 96-well wetlab yield"),
    ("wetlab_interpreter_hit_yield", 0.14, "direct assay utility after v6"),
    ("external_comparator_delta", 0.10, "must beat or complement BigMHC/PRIME/NetMHCpan"),
    ("modality_completeness", 0.06, "penalizes impossible candidates"),
    ("novelty_diversity", 0.05, "prevents copying one public model"),
    ("interpretability", 0.03, "supports reviewer defense"),
    ("complexity_penalty", 0.02, "keeps deployable search"),
]


GUARDRAILS = [
    ("public_secret_boundary", "Use public papers/code/features only; no proprietary extraction."),
    ("source_leakage_guard", "Patient/source/HLA/peptide-cluster splits are separate required reports."),
    ("sparse_proxy_guard", "Sparse features such as TCR availability cannot dominate without ablation."),
    ("single_concept_cap", "No primitive family can exceed a prespecified weight cap in claim-safe mode."),
    ("random_search_baseline", "GA/RL gains must beat fixed ensemble and random search under same budget."),
    ("locked_wetlab_validation", "Prospective claims require locked score before v6 wetlab results."),
    ("modality_missingness_audit", "Every modality branch needs missingness and availability stratification."),
    ("complexity_budget", "Report search generations, population, model count and CPU/GPU budget."),
]


def split_semicolon(text: str) -> list[str]:
    return [x.strip() for x in str(text).split(";") if x.strip()]


def build_tables(out_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    methods = pd.DataFrame(METHOD_CARDS)
    actions = pd.DataFrame(ACTION_CARDS)
    reward = pd.DataFrame(REWARD_COMPONENTS, columns=["reward_component", "weight", "purpose"])
    guards = pd.DataFrame(GUARDRAILS, columns=["guardrail", "rule"])
    methods.to_csv(out_dir / "method_primitives.tsv", sep="\t", index=False)
    actions.to_csv(out_dir / "rl_ga_action_space.tsv", sep="\t", index=False)
    reward.to_csv(out_dir / "reward_components.tsv", sep="\t", index=False)
    guards.to_csv(out_dir / "guardrails.tsv", sep="\t", index=False)
    return methods, actions, reward, guards


def build_graph(methods: pd.DataFrame, actions: pd.DataFrame, reward: pd.DataFrame, guards: pd.DataFrame, out_dir: Path) -> tuple[nx.DiGraph, pd.DataFrame, pd.DataFrame]:
    graph = nx.DiGraph()
    graph.add_node("CROSS-Neo_Darwin_RL", node_type="new_algorithm", layer="controller")
    for _, row in methods.iterrows():
        graph.add_node(row["method"], node_type="source_method", domain=row["domain"], source_type=row["source_type"])
        graph.add_edge(row["method"], "CROSS-Neo_Darwin_RL", relation="contributes_primitives")
        for field, ntype in [
            ("core_primitives", "core_primitive"),
            ("preprocessing_primitives", "preprocessing_primitive"),
            ("modalities", "modality"),
            ("search_gene", "search_gene"),
            ("risk_guard", "risk_guard"),
        ]:
            for token in split_semicolon(row[field]):
                graph.add_node(token, node_type=ntype, layer=field)
                graph.add_edge(row["method"], token, relation=f"has_{ntype}")
                graph.add_edge(token, "CROSS-Neo_Darwin_RL", relation="available_to_search")
    for _, row in actions.iterrows():
        graph.add_node(row["action"], node_type="rl_action", action_family=row["action_family"])
        graph.add_edge("CROSS-Neo_Darwin_RL", row["action"], relation="controller_can_choose")
        for choice in split_semicolon(row["choices"]):
            graph.add_node(choice, node_type="action_choice", action_family=row["action_family"])
            graph.add_edge(row["action"], choice, relation="has_choice")
    for _, row in reward.iterrows():
        graph.add_node(row["reward_component"], node_type="reward", weight=float(row["weight"]))
        graph.add_edge("CROSS-Neo_Darwin_RL", row["reward_component"], relation="optimizes")
    for _, row in guards.iterrows():
        graph.add_node(row["guardrail"], node_type="guardrail")
        graph.add_edge(row["guardrail"], "CROSS-Neo_Darwin_RL", relation="constrains")

    nodes = pd.DataFrame([{"node": node, **attrs} for node, attrs in graph.nodes(data=True)])
    edges = pd.DataFrame([{"source": u, "target": v, **attrs} for u, v, attrs in graph.edges(data=True)])
    nodes.to_csv(out_dir / "method_kg_nodes.tsv", sep="\t", index=False)
    edges.to_csv(out_dir / "method_kg_edges.tsv", sep="\t", index=False)
    nx.write_graphml(graph, out_dir / "cross_neo_darwin_rl_method_graph.graphml")
    return graph, nodes, edges


def modality_matrix(methods: pd.DataFrame) -> pd.DataFrame:
    modalities = sorted({m for row in methods["modalities"] for m in split_semicolon(row)})
    mat = []
    for _, row in methods.iterrows():
        present = set(split_semicolon(row["modalities"]))
        mat.append({"method": row["method"], **{m: int(m in present) for m in modalities}})
    return pd.DataFrame(mat)


def primitive_summary(methods: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for field in ["core_primitives", "preprocessing_primitives", "modalities", "search_gene", "risk_guard"]:
        counts: dict[str, int] = {}
        for row in methods[field]:
            for token in split_semicolon(row):
                counts[token] = counts.get(token, 0) + 1
        for token, n in counts.items():
            rows.append({"primitive_type": field, "primitive": token, "n_methods": n})
    return pd.DataFrame(rows).sort_values(["primitive_type", "n_methods", "primitive"], ascending=[True, False, True])


def write_algorithm_blueprint(out_dir: Path, methods: pd.DataFrame, actions: pd.DataFrame, reward: pd.DataFrame, guards: pd.DataFrame) -> dict[str, object]:
    blueprint = {
        "name": "CROSS-Neo Darwin-RL",
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "thesis": "A literature-mined knowledge graph decomposes public AI and biomedical methods into primitives; GA evolves candidate pipelines while an RL controller learns which mutation/search actions to use under leakage, modality and wetlab constraints.",
        "genome": {
            "preprocessing_genes": [
                "label_confidence_model",
                "split_policy",
                "HLA_normalization",
                "peptide_length_filter",
                "patient_source_cluster_guard",
                "modality_missingness_policy",
            ],
            "modality_genes": [
                "peptide_sequence",
                "HLA_allele",
                "MS_EL",
                "binding_affinity",
                "RNA_expression",
                "TCR_CDR3beta",
                "pMHC_structure",
                "biochemical_features",
                "gene_context",
                "assay_metadata",
            ],
            "model_genes": [
                "tabular_GBDT",
                "CNN_branch",
                "LSTM_branch",
                "protein_LM_embedding",
                "structure_GNN",
                "late_fusion_attention",
                "mixture_of_experts",
                "claim_safe_cap",
            ],
            "learning_genes": [
                "binary_loss",
                "confidence_weighted_loss",
                "pairwise_ranking_loss",
                "top_k_surrogate_loss",
                "calibration_loss",
                "distilled_rule_head",
            ],
        },
        "rl_controller": {
            "state": [
                "source_split_metrics",
                "low_leakage_metrics",
                "HLA_coverage",
                "modality_missingness",
                "comparator_disagreement",
                "search_diversity",
                "wetlab_budget_remaining",
            ],
            "actions": actions["action"].tolist(),
            "policy": "contextual bandit first; upgrade to PPO/MCTS only after enough search traces exist",
        },
        "reward": reward.to_dict(orient="records"),
        "guardrails": guards.to_dict(orient="records"),
        "source_methods": methods[["method", "source_url", "paper_anchor"]].to_dict(orient="records"),
        "locked_validation_path": "Run GA/RL on retrospective data only, freeze architecture and thresholds, then evaluate actual v6 96-well interpreter calls.",
    }
    (out_dir / "cross_neo_darwin_rl_algorithm_blueprint.json").write_text(
        json.dumps(blueprint, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return blueprint


def plot_outputs(
    out_dir: Path,
    graph: nx.DiGraph,
    methods: pd.DataFrame,
    actions: pd.DataFrame,
    reward: pd.DataFrame,
    guards: pd.DataFrame,
) -> None:
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    mm = modality_matrix(methods)
    mm.to_csv(out_dir / "modality_matrix.tsv", sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(12, 6.2))
    data = mm.drop(columns=["method"]).to_numpy()
    ax.imshow(data, aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
    ax.set_yticks(np.arange(len(mm)))
    ax.set_yticklabels(mm["method"], fontsize=8)
    ax.set_xticks(np.arange(len(mm.columns) - 1))
    ax.set_xticklabels(mm.columns[1:], rotation=55, ha="right", fontsize=8)
    ax.set_title("Method-to-modality coverage")
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            ax.text(j, i, str(int(data[i, j])), ha="center", va="center", fontsize=7, color="#0b1820")
    fig.tight_layout()
    fig.savefig(fig_dir / "fig1_method_modality_matrix.png", dpi=190)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5.2))
    rr = reward.iloc[::-1]
    ax.barh(rr["reward_component"], rr["weight"], color="#2f6f73")
    ax.set_xlabel("Reward weight")
    ax.set_title("Darwin-RL multi-objective reward")
    for y, v in enumerate(rr["weight"]):
        ax.text(v + 0.005, y, f"{v:.2f}", va="center", fontsize=8)
    ax.set_xlim(0, max(reward["weight"]) + 0.07)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig2_reward_stack.png", dpi=190)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(12, 7.5))
    node_type_colors = {
        "new_algorithm": "#c59b3b",
        "source_method": "#69a7a2",
        "core_primitive": "#d7dde2",
        "preprocessing_primitive": "#9fb7d7",
        "modality": "#8fcf92",
        "search_gene": "#d7b56d",
        "risk_guard": "#d16b5f",
        "rl_action": "#8f9cc2",
        "action_choice": "#b8c0ce",
        "reward": "#cda65b",
        "guardrail": "#d16b5f",
    }
    pos = nx.spring_layout(graph, seed=7, k=0.85, iterations=80)
    colors = [node_type_colors.get(graph.nodes[n].get("node_type", ""), "#cccccc") for n in graph.nodes]
    sizes = [850 if n == "CROSS-Neo_Darwin_RL" else 220 if graph.nodes[n].get("node_type") == "source_method" else 70 for n in graph.nodes]
    nx.draw_networkx_nodes(graph, pos, node_color=colors, node_size=sizes, linewidths=0.4, edgecolors="#1f2a31", ax=ax)
    nx.draw_networkx_edges(graph, pos, alpha=0.22, width=0.55, edge_color="#9aa7af", arrows=False, ax=ax)
    label_nodes = [n for n in graph.nodes if graph.nodes[n].get("node_type") in {"new_algorithm", "source_method", "rl_action", "reward", "guardrail"}]
    nx.draw_networkx_labels(graph, pos, labels={n: n for n in label_nodes}, font_size=6.2, ax=ax)
    ax.set_title("Literature-mined method primitive knowledge graph")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(fig_dir / "fig3_method_knowledge_graph.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 5.8))
    ax.axis("off")
    boxes = [
        ("Web / papers", 0.07, 0.62, "#69a7a2"),
        ("Method cards", 0.25, 0.62, "#9fb7d7"),
        ("Primitive KG", 0.43, 0.62, "#c59b3b"),
        ("GA genome", 0.61, 0.62, "#d7b56d"),
        ("RL controller", 0.79, 0.62, "#8f9cc2"),
        ("Retrospective CV", 0.25, 0.24, "#d7dde2"),
        ("Guardrail audit", 0.43, 0.24, "#d16b5f"),
        ("Locked 96-well", 0.61, 0.24, "#2f6f73"),
        ("Interpreter claims", 0.79, 0.24, "#c59b3b"),
    ]
    for text, x, y, color in boxes:
        ax.add_patch(plt.Rectangle((x, y), 0.145, 0.16, facecolor=color, edgecolor="#26333c", lw=1.2, alpha=0.95))
        ax.text(x + 0.0725, y + 0.08, text, ha="center", va="center", fontsize=9, color="#101418", fontweight="bold")
    arrows = [
        ((0.215, 0.70), (0.25, 0.70)),
        ((0.395, 0.70), (0.43, 0.70)),
        ((0.575, 0.70), (0.61, 0.70)),
        ((0.755, 0.70), (0.79, 0.70)),
        ((0.862, 0.62), (0.862, 0.40)),
        ((0.79, 0.32), (0.755, 0.32)),
        ((0.61, 0.32), (0.575, 0.32)),
        ((0.43, 0.32), (0.395, 0.32)),
        ((0.322, 0.40), (0.322, 0.62)),
    ]
    for start, end in arrows:
        ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="->", color="#dce3e7", lw=1.4))
    ax.set_title("Darwin-RL method discovery loop")
    fig.tight_layout()
    fig.savefig(fig_dir / "fig4_darwin_rl_loop.png", dpi=190)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    guard_counts = [len(split_semicolon(row)) for row in methods["risk_guard"]]
    ax.bar(methods["method"], guard_counts, color="#d16b5f")
    ax.set_ylabel("Risk guards per method card")
    ax.set_title("Reviewer-risk surface captured in the primitive graph")
    ax.tick_params(axis="x", rotation=55, labelsize=8)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig5_risk_guard_surface.png", dpi=190)
    plt.close(fig)


def table_html(df: pd.DataFrame, max_rows: int = 30) -> str:
    if df.empty:
        return "<p>No rows.</p>"
    show = df.head(max_rows).copy()
    for col in show.columns:
        if pd.api.types.is_float_dtype(show[col]):
            show[col] = show[col].map(lambda x: "" if pd.isna(x) else f"{x:.3f}")
    return show.to_html(index=False, escape=True, classes="data")


def write_report(
    out_dir: Path,
    methods: pd.DataFrame,
    actions: pd.DataFrame,
    reward: pd.DataFrame,
    guards: pd.DataFrame,
    primitives: pd.DataFrame,
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    blueprint: dict[str, object],
) -> None:
    lines = [
        "# CROSS-Neo Darwin-RL method-mining blueprint",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Bottom line",
        "",
        "- Build a new algorithm from public method primitives, not from copied black-box scores.",
        "- The search space includes preprocessing, modality activation, encoders, fusion, reward composition and wetlab arm allocation.",
        "- GA evolves candidate pipelines; an RL controller learns which mutation/search actions to use from search history.",
        "- Prospective claims require freezing the selected architecture before the actual v6 96-well interpreter run.",
        "",
        "## Knowledge graph size",
        "",
        f"- Nodes: {len(nodes)}",
        f"- Edges: {len(edges)}",
        f"- Source method cards: {len(methods)}",
        f"- RL/GA action families: {actions['action_family'].nunique()}",
        "",
        "## Source method cards",
        "",
        methods.to_markdown(index=False),
        "",
        "## RL/GA action space",
        "",
        actions.to_markdown(index=False),
        "",
        "## Reward function",
        "",
        reward.to_markdown(index=False),
        "",
        "## Guardrails",
        "",
        guards.to_markdown(index=False),
        "",
        "## Primitive summary",
        "",
        primitives.to_markdown(index=False),
        "",
        "## Claim boundary",
        "",
        "This dossier is a design and search-space artifact. It does not claim that the final algorithm is prospectively superior until a locked candidate architecture is evaluated on actual wetlab v6 interpreter calls.",
        "",
    ]
    (out_dir / "CROSS_NEO_DARWIN_RL_BLUEPRINT_REPORT.md").write_text("\n".join(lines), encoding="utf-8")


def write_html(
    out_dir: Path,
    methods: pd.DataFrame,
    actions: pd.DataFrame,
    reward: pd.DataFrame,
    guards: pd.DataFrame,
    primitives: pd.DataFrame,
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
) -> dict[str, object]:
    asset_dir = HUB / "assets/cross_neo_method_mining_rl_ga_blueprint"
    asset_dir.mkdir(parents=True, exist_ok=True)
    for fig in (out_dir / "figures").glob("*.png"):
        shutil.copy2(fig, asset_dir / fig.name)

    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>CROSS-Neo Darwin-RL Blueprint</title>
  <style>
    :root {{ --bg:#101418; --panel:#151c22; --ink:#e9edf0; --muted:#9aa7af; --gold:#c59b3b; --line:#29343c; --teal:#69a7a2; --red:#d16b5f; }}
    body {{ margin:0; background:var(--bg); color:var(--ink); font:15px/1.55 system-ui, -apple-system, Segoe UI, sans-serif; }}
    header {{ padding:42px 5vw 30px; background:#0c1115; border-bottom:1px solid var(--line); }}
    .kicker {{ color:var(--gold); text-transform:uppercase; letter-spacing:.12em; font-weight:700; font-size:12px; }}
    h1 {{ margin:.3rem 0 .65rem; font-size:clamp(30px,4vw,54px); line-height:1.05; }}
    .lead {{ color:#d0d8dd; max-width:1040px; font-size:18px; }}
    .stats {{ display:grid; grid-template-columns:repeat(5,minmax(130px,1fr)); gap:10px; margin-top:22px; }}
    .stat {{ border:1px solid var(--line); background:var(--panel); padding:13px 14px; border-radius:8px; }}
    .stat b {{ display:block; font-size:23px; }}
    .stat span {{ color:var(--muted); font-size:12px; }}
    main {{ padding:30px 5vw 60px; display:grid; grid-template-columns:270px minmax(0,1fr); gap:28px; }}
    nav {{ position:sticky; top:0; align-self:start; max-height:100vh; overflow:auto; padding:14px; border:1px solid var(--line); border-radius:8px; background:#111820; }}
    nav a {{ display:block; color:#d6dee2; text-decoration:none; padding:7px 0; border-bottom:1px solid #202b32; }}
    h2 {{ border-bottom:1px solid var(--line); padding-bottom:8px; }}
    .num {{ color:var(--gold); margin-right:8px; }}
    .note {{ border-left:3px solid var(--gold); background:#171f26; padding:10px 14px; color:#dce3e7; }}
    table.data {{ width:100%; border-collapse:collapse; font-size:13px; margin:12px 0 28px; }}
    table.data th, table.data td {{ border-bottom:1px solid var(--line); padding:7px 8px; text-align:left; vertical-align:top; }}
    table.data th {{ color:#f4d891; background:#141b21; }}
    .grid {{ display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:18px; }}
    img {{ max-width:100%; border:1px solid var(--line); border-radius:8px; background:#fff; }}
    code {{ color:#f4d891; }}
    @media(max-width:900px) {{ main {{ grid-template-columns:1fr; }} nav {{ position:static; }} .stats {{ grid-template-columns:repeat(2,1fr); }} .grid {{ grid-template-columns:1fr; }} }}
  </style>
</head>
<body>
<header>
  <div class="kicker">CROSS-Neo · web-mined methods · GA/RL algorithm discovery</div>
  <h1>Darwin-RL blueprint for a new immunogenicity algorithm</h1>
  <p class="lead">Public AI and biomedical methods are decomposed into preprocessing, modality, encoder, fusion, reward and guardrail primitives. GA searches pipeline genomes; an RL controller learns which search actions to use under leakage and wetlab constraints.</p>
  <div class="stats">
    <div class="stat"><b>{len(methods)}</b><span>source method cards</span></div>
    <div class="stat"><b>{len(nodes)}</b><span>KG nodes</span></div>
    <div class="stat"><b>{len(edges)}</b><span>KG edges</span></div>
    <div class="stat"><b>{actions['action_family'].nunique()}</b><span>action families</span></div>
    <div class="stat"><b>{len(guards)}</b><span>guardrails</span></div>
  </div>
</header>
<main>
<nav>
  <a href="#summary">01 Summary</a>
  <a href="#figures">02 Figures</a>
  <a href="#methods">03 Method Cards</a>
  <a href="#actions">04 RL/GA Actions</a>
  <a href="#reward">05 Reward</a>
  <a href="#guards">06 Guardrails</a>
  <a href="#primitives">07 Primitive Library</a>
</nav>
<article>
<section id="summary">
  <h2><span class="num">01</span>Summary</h2>
  <p class="note">This is the next algorithm concept: not a weighted ensemble, but a literature-mined method graph with GA/RL search over preprocessing, modality, architecture, fusion and assay-validation choices.</p>
</section>
<section id="figures">
  <h2><span class="num">02</span>Figures</h2>
  <div class="grid">
    <img src="assets/cross_neo_method_mining_rl_ga_blueprint/fig1_method_modality_matrix.png" alt="method modality matrix">
    <img src="assets/cross_neo_method_mining_rl_ga_blueprint/fig2_reward_stack.png" alt="reward stack">
  </div>
  <div class="grid">
    <img src="assets/cross_neo_method_mining_rl_ga_blueprint/fig4_darwin_rl_loop.png" alt="Darwin RL loop">
    <img src="assets/cross_neo_method_mining_rl_ga_blueprint/fig5_risk_guard_surface.png" alt="risk guard surface">
  </div>
  <img src="assets/cross_neo_method_mining_rl_ga_blueprint/fig3_method_knowledge_graph.png" alt="method knowledge graph">
</section>
<section id="methods">
  <h2><span class="num">03</span>Source Method Cards</h2>
  {table_html(methods, 20)}
</section>
<section id="actions">
  <h2><span class="num">04</span>RL/GA Action Space</h2>
  {table_html(actions, 20)}
</section>
<section id="reward">
  <h2><span class="num">05</span>Reward</h2>
  {table_html(reward, 20)}
</section>
<section id="guards">
  <h2><span class="num">06</span>Guardrails</h2>
  {table_html(guards, 20)}
</section>
<section id="primitives">
  <h2><span class="num">07</span>Primitive Library</h2>
  {table_html(primitives, 80)}
</section>
</article>
</main>
</body>
</html>
"""
    html_path = HUB / "cross_neo_method_mining_rl_ga_blueprint.html"
    html_path.write_text(html_text, encoding="utf-8")
    live_ok = False
    warnings: list[str] = []
    try:
        live_asset_dir = LIVE_HUB / "assets/cross_neo_method_mining_rl_ga_blueprint"
        live_asset_dir.mkdir(parents=True, exist_ok=True)
        for fig in asset_dir.glob("*.png"):
            shutil.copy2(fig, live_asset_dir / fig.name)
        shutil.copy2(html_path, LIVE_HUB / html_path.name)
        live_ok = True
    except Exception as exc:
        warnings.append(str(exc))
    return {
        "html_path": str(html_path),
        "live_html_path": str(LIVE_HUB / html_path.name),
        "live_deploy_ok": live_ok,
        "live_deploy_warnings": warnings,
    }


def main() -> None:
    out_dir = DEFAULT_OUT
    out_dir.mkdir(parents=True, exist_ok=True)
    methods, actions, reward, guards = build_tables(out_dir)
    primitives = primitive_summary(methods)
    primitives.to_csv(out_dir / "primitive_library_summary.tsv", sep="\t", index=False)
    graph, nodes, edges = build_graph(methods, actions, reward, guards, out_dir)
    blueprint = write_algorithm_blueprint(out_dir, methods, actions, reward, guards)
    plot_outputs(out_dir, graph, methods, actions, reward, guards)
    write_report(out_dir, methods, actions, reward, guards, primitives, nodes, edges, blueprint)
    html_info = write_html(out_dir, methods, actions, reward, guards, primitives, nodes, edges)
    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "out_dir": str(out_dir),
        "n_method_cards": int(len(methods)),
        "n_actions": int(len(actions)),
        "n_reward_components": int(len(reward)),
        "n_guardrails": int(len(guards)),
        "n_kg_nodes": int(len(nodes)),
        "n_kg_edges": int(len(edges)),
        "algorithm_name": blueprint["name"],
        "claim_boundary": "design/search-space artifact; prospective superiority requires locked wetlab v6 validation",
        **html_info,
        "output_files": [
            str(out_dir / "method_primitives.tsv"),
            str(out_dir / "rl_ga_action_space.tsv"),
            str(out_dir / "reward_components.tsv"),
            str(out_dir / "guardrails.tsv"),
            str(out_dir / "cross_neo_darwin_rl_algorithm_blueprint.json"),
            str(out_dir / "CROSS_NEO_DARWIN_RL_BLUEPRINT_REPORT.md"),
        ],
    }
    (out_dir / "cross_neo_darwin_rl_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
