#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "project" / "results" / "hla_deepdive_2026_05_08"
OUT = ROOT / "project" / "results" / "hla_two_paper_synthesis_2026_05_09"
REPORTS = ROOT / "project" / "reports"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fmt_p(x: float | str | None) -> str:
    if x is None or x == "":
        return "NA"
    x = float(x)
    if x == 0:
        return "<1e-300"
    if abs(x) < 0.001:
        return f"{x:.1e}"
    if abs(x) < 0.01:
        return f"{x:.3g}"
    return f"{x:.3g}"


def fmt_num(x: float | str | None, nd: int = 2) -> str:
    if x is None or x == "":
        return "NA"
    x = float(x)
    return f"{x:.{nd}f}"


def find_row(rows: list[dict[str, str]], **criteria: str) -> dict[str, str]:
    for row in rows:
        if all(row.get(k) == v for k, v in criteria.items()):
            return row
    raise KeyError(criteria)


def first_json_item(items: list[dict], key: str, value: str) -> dict:
    for item in items:
        if item.get(key) == value:
            return item
    raise KeyError((key, value))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)

    track1 = load_json(RESULTS / "track1_paper4_gd_panasian" / "track1_summary.json")
    track3 = load_json(RESULTS / "track3_pan_asian_atlas" / "track3_summary.json")
    track4 = load_json(RESULTS / "track4_ht_ptc_hla2" / "track4_summary.json")
    track5 = load_json(RESULTS / "track5_dm1_hla1_module" / "track5_summary.json")
    track6 = load_json(RESULTS / "track6_pancan_hla_dm1" / "track6_summary.json")
    track12 = load_json(RESULTS / "track12_scrna_hla_celltype" / "track12_summary.json")
    track18 = load_json(RESULTS / "track18_ici_hla" / "track18_summary.json")
    track19 = load_json(RESULTS / "track19_aire_tolerance" / "track19_summary.json")
    track21 = load_json(RESULTS / "track21_hla_ai_proxy" / "track21_summary.json")
    track25 = load_json(RESULTS / "track25_hla_somatic_mut" / "track25_summary.json")
    track26 = load_json(RESULTS / "track26_ifng_hla_dm1" / "track26_summary.json")
    track27 = load_json(RESULTS / "track27_hla_g_nonclassical" / "track27_summary.json")
    track28 = load_json(RESULTS / "track28_immune_cold_thca" / "track28_summary.json")
    track29 = load_json(RESULTS / "track29_neoantigen_hla" / "track29_summary.json")
    track33 = load_json(RESULTS / "track33_spatial_tls_hla" / "track33_summary.json")

    gd_v3 = load_tsv(RESULTS / "track10_korean_lit" / "tables" / "T05_panasian_GD_DL_random_effects_v3.tsv")
    gd_ready = load_tsv(RESULTS / "track10_korean_lit" / "tables" / "T11_readiness_scoreboard_v3.tsv")
    specificity = load_tsv(RESULTS / "track31_cross_autoimmune" / "tables" / "T05_specificity_scores.tsv")
    aitd_breadth = load_tsv(RESULTS / "track31_cross_autoimmune" / "tables" / "T06b_AITD_vs_broad_summary.tsv")
    trans = load_tsv(RESULTS / "track11_finngen_panukbb" / "tables" / "T05_track1_tagsnp_transancestry.tsv")
    ht_vs_base = load_tsv(RESULTS / "track4_ht_ptc_hla2" / "tables" / "T3_HT_overlap_PTC_vs_baseline.tsv")
    ht_within = load_tsv(RESULTS / "track4_ht_ptc_hla2" / "tables" / "T4_within_GSE286332_HT_vs_HTneg.tsv")
    power = load_tsv(RESULTS / "track4_ht_ptc_hla2" / "tables" / "T9_power_calc_top3_candidates.tsv")
    rai_hla = load_tsv(RESULTS / "track34_rai_hla_decouple" / "tables" / "T02_headline_correlations.tsv")

    t1 = {r["allele"]: r for r in track1["top_alleles"]}
    v3 = {r["allele"]: r for r in gd_v3}
    ready = {r["allele"]: r for r in gd_ready}
    spec = {r["allele"]: r for r in specificity}
    breadth = {r["allele"]: r for r in aitd_breadth}
    dpb1_finngen_gd = find_row(trans, allele="HLA-DPB1*05:01", source="FinnGen_R12", phenocode="E4_GRAVES_STRICT")
    dpb1_finngen_ht = find_row(trans, allele="HLA-DPB1*05:01", source="FinnGen_R12", phenocode="E4_HYTHY_AI_STRICT")
    drb10405 = find_row(ht_vs_base, allele="DRB1*04:05")
    dpb10501_ht = find_row(ht_vs_base, allele="DPB1*05:01")
    drb10405_power = find_row(power, candidate="DRB1*04:05", or_target="1.8")
    dqb10601_within = find_row(ht_within, allele="DQB1*06:01")
    dm1_hla1 = first_json_item(track5["headline_correlations"], "comparison", "DM1 vs HLA-I module")
    dm1_hla2 = first_json_item(track5["headline_correlations"], "comparison", "DM1 vs HLA-II module")
    hla1_purity = first_json_item(track5["purity_partial"], "module", "HLA-I")
    hla2_purity = first_json_item(track5["purity_partial"], "module", "HLA-II")
    ifng_hla1 = first_json_item(track26["headline_correlations"], "comparison", "IFN-gamma Hallmark vs HLA-I module")
    ifng_hla2 = first_json_item(track26["headline_correlations"], "comparison", "IFN-gamma Hallmark vs HLA-II module")
    med_hla1 = first_json_item(track26["mediation_baron_kenny_sobel"], "label", "DM1 -> IFNG Hallmark -> HLA-I")
    med_hla2 = first_json_item(track26["mediation_baron_kenny_sobel"], "label", "DM1 -> IFNG Hallmark -> HLA-II")
    rai_dm1 = find_row(rai_hla, comparison="DM1 vs RAI module")
    rai_hla2 = find_row(rai_hla, comparison="HLA-II module vs RAI module")
    hlaga = first_json_item(track27["dm1_classical_nonclassical_corr"], "gene", "HLA-G")
    hlag_decouple = first_json_item(track27["decoupling_dm1"], "nonclassical_gene", "HLA-G")

    hypotheses: list[dict[str, str]] = [
        {
            "paper": "Paper 4 GD HLA",
            "id": "GD-H1",
            "hypothesis": "DPB1*05:01 is the most publication-ready Pan-Asian Graves/AITD HLA risk allele.",
            "status": "supported; manuscript-grade with caveats",
            "key_result": f"v3 pooled OR {fmt_num(v3['DPB1*05:01']['pooled_or'])} ({fmt_num(v3['DPB1*05:01']['ci_lo'])}-{fmt_num(v3['DPB1*05:01']['ci_hi'])}), k={v3['DPB1*05:01']['k']}, I2={fmt_num(v3['DPB1*05:01']['I2_pct'])}%; v2 PI {fmt_num(t1['DPB1*05:01']['pi_lo'])}-{fmt_num(t1['DPB1*05:01']['pi_hi'])}.",
            "source": "track10/T05_panasian_GD_DL_random_effects_v3.tsv; track1_summary.json",
            "claim_grade": ready["DPB1*05:01"]["grade_v3"],
        },
        {
            "paper": "Paper 4 GD HLA",
            "id": "GD-H2",
            "hypothesis": "B*46:01 is a real Pan-Asian GD risk allele but effect size is population-heterogeneous.",
            "status": "supported, but heterogeneity limits single-effect storytelling",
            "key_result": f"v3 pooled OR {fmt_num(v3['B*46:01']['pooled_or'])} ({fmt_num(v3['B*46:01']['ci_lo'])}-{fmt_num(v3['B*46:01']['ci_hi'])}), I2={fmt_num(v3['B*46:01']['I2_pct'])}%, PI {fmt_num(v3['B*46:01']['pi_lo'])}-{fmt_num(v3['B*46:01']['pi_hi'])}.",
            "source": "track10/T05_panasian_GD_DL_random_effects_v3.tsv",
            "claim_grade": ready["B*46:01"]["grade_v3"],
        },
        {
            "paper": "Paper 4 GD HLA",
            "id": "GD-H3",
            "hypothesis": "A*02:07 is now upgraded from partial to full Pan-Asian meta-ready after Korean literature extraction.",
            "status": "supported",
            "key_result": f"v3 pooled OR {fmt_num(v3['A*02:07']['pooled_or'])} ({fmt_num(v3['A*02:07']['ci_lo'])}-{fmt_num(v3['A*02:07']['ci_hi'])}); grade changed {ready['A*02:07']['grade_v2']} -> {ready['A*02:07']['grade_v3']}.",
            "source": "track10/T05_panasian_GD_DL_random_effects_v3.tsv; track10/T11_readiness_scoreboard_v3.tsv",
            "claim_grade": ready["A*02:07"]["grade_v3"],
        },
        {
            "paper": "Paper 4 GD HLA",
            "id": "GD-H4",
            "hypothesis": "C*01:02 is the most GD-specific allele in the cross-autoimmune panel.",
            "status": "supported as literature-curated specificity, not de novo cohort proof",
            "key_result": f"GD specificity score {fmt_num(spec['C*01:02']['GD_specificity_score'])}; non-AITD mean abs logOR {fmt_num(spec['C*01:02']['non_AITD_mean_abs_logOR'], 3)}; AITD breadth ratio {fmt_num(breadth['C*01:02']['AITD_breadth_ratio'])}.",
            "source": "track31/T05_specificity_scores.tsv; track31/T06b_AITD_vs_broad_summary.tsv",
            "claim_grade": "specificity-supporting",
        },
        {
            "paper": "Paper 4 GD HLA",
            "id": "GD-H5",
            "hypothesis": "DPB1*05:01 is AITD-broad rather than GD-specific or pan-autoimmune.",
            "status": "supported with publication-bias caveat",
            "key_result": f"GD specificity score {fmt_num(spec['DPB1*05:01']['GD_specificity_score'])}; AITD breadth ratio {fmt_num(breadth['DPB1*05:01']['AITD_breadth_ratio'])}; n strong outside GD={spec['DPB1*05:01']['n_strong_outside_GD']}.",
            "source": "track31/T05_specificity_scores.tsv; track31/T06b_AITD_vs_broad_summary.tsv",
            "claim_grade": "mechanistic framing only",
        },
        {
            "paper": "Paper 4 GD HLA",
            "id": "GD-H6",
            "hypothesis": "The DPB1*05:01 signal is not merely a B*46:01/A*02:07 haplotype passenger in Koreans.",
            "status": "supported directionally, underpowered for final LD map",
            "key_result": "K2 carrier-level LD: DPB1*05:01 vs B*46:01 r2≈0.077, p≈0.073; B*46:01 vs A*02:07 r2≈0.337.",
            "source": "track1_report.md; track1/tables/T07_LD_matrix_K2_carrier_level.tsv",
            "claim_grade": "supportive LD context",
        },
        {
            "paper": "Paper 4 GD HLA",
            "id": "GD-H7",
            "hypothesis": "European/Finnish public GWAS corroborates an MHC class-II thyroid autoimmunity signal at the DPB1*05:01 tag, but not Korean effect size.",
            "status": "supported as trans-ancestry sensitivity only",
            "key_result": f"FinnGen GD strict rs9277534 beta {fmt_num(dpb1_finngen_gd['beta'], 3)}, p={fmt_p(dpb1_finngen_gd['pval'])}; FinnGen HT-AI beta {fmt_num(dpb1_finngen_ht['beta'], 3)}, p={fmt_p(dpb1_finngen_ht['pval'])}, opposite direction.",
            "source": "track11/T05_track1_tagsnp_transancestry.tsv",
            "claim_grade": "corroborative, not transferable",
        },
        {
            "paper": "Paper 4 GD HLA",
            "id": "GD-H8",
            "hypothesis": "Healthy Korean HLA atlas supports using arcasHLA RNA-seq for population-level frequency work.",
            "status": "supported for resource/methods layer",
            "key_result": f"Healthy atlas: {track3['n_populations']} populations, {track3['n_alleles_in_master_matrix']} alleles; Korean platform concordance reported at r≈0.86-0.93; K2 normal n={track3['k2_normal_n']}, GSE213647 normal n={track3['gse213647_normal_n']}.",
            "source": "track3_summary.json; track3_report.md",
            "claim_grade": "resource support",
        },
        {
            "paper": "Paper 4 GD HLA",
            "id": "GD-H9",
            "hypothesis": "DQB1*02:01 / DRB1*07:01 are protective GD alleles, but the Korean evidence differs by allele.",
            "status": "partially supported",
            "key_result": f"DRB1*07:01 v3 pooled OR {fmt_num(v3['DRB1*07:01']['pooled_or'])} ({fmt_num(v3['DRB1*07:01']['ci_lo'])}-{fmt_num(v3['DRB1*07:01']['ci_hi'])}), k={v3['DRB1*07:01']['k']}; DQB1*02:01 remains single-source OR {fmt_num(v3['DQB1*02:01']['pooled_or'])}.",
            "source": "track10/T05_panasian_GD_DL_random_effects_v3.tsv",
            "claim_grade": "DRB1 stronger than DQB1",
        },
        {
            "paper": "Paper 4 GD HLA",
            "id": "GD-H10",
            "hypothesis": "A Korean adult GD NGS cohort is the key upgrade from strong meta-analysis to top-tier original study.",
            "status": "open; new data required",
            "key_result": "Readiness table marks DPB1*05:01/B*46:01/A*02:07 as ready-to-publish, but Korean adult NGS remains the A+ blocker for de novo Korean claim.",
            "source": "track10/T11_readiness_scoreboard_v3.tsv; track1_report.md",
            "claim_grade": "must-do validation",
        },
        {
            "paper": "Paper 2 HT-PTC HLA/Immune",
            "id": "PTC-H1",
            "hypothesis": "HT-overlap PTC has one defensible allele candidate: DRB1*04:05.",
            "status": "candidate only; underpowered and FDR-negative",
            "key_result": f"GSE286332 HT-overlap n=9: DRB1*04:05 OR {fmt_num(drb10405['or_allele'])} ({fmt_num(drb10405['ci_lo_allele'])}-{fmt_num(drb10405['ci_hi_allele'])}), Fisher p={fmt_p(drb10405['p_allele'])}, FDR={fmt_p(drb10405['fdr_bh_allele'])}; n/group for OR=1.8 ≈ {drb10405_power['n_per_group_required']}.",
            "source": "track4/T3_HT_overlap_PTC_vs_baseline.tsv; track4/T9_power_calc_top3_candidates.tsv",
            "claim_grade": "hypothesis-generating only",
        },
        {
            "paper": "Paper 2 HT-PTC HLA/Immune",
            "id": "PTC-H2",
            "hypothesis": "DPB1*05:01 is not currently an enriched HT-PTC allele despite its AITD literature role.",
            "status": "not supported in current HT-PTC allele data",
            "key_result": f"HT-overlap DPB1*05:01 OR {fmt_num(dpb10501_ht['or_allele'])} ({fmt_num(dpb10501_ht['ci_lo_allele'])}-{fmt_num(dpb10501_ht['ci_hi_allele'])}), p={fmt_p(dpb10501_ht['p_allele'])}; carrier OR {fmt_num(dpb10501_ht['or_carrier'])}, direction depletion vs Kim carrier baseline.",
            "source": "track4/T3_HT_overlap_PTC_vs_baseline.tsv",
            "claim_grade": "negative/weak candidate",
        },
        {
            "paper": "Paper 2 HT-PTC HLA/Immune",
            "id": "PTC-H3",
            "hypothesis": "Within GSE286332, HT-overlap vs HT-negative PTC allele differences are not significant.",
            "status": "not validated",
            "key_result": f"Largest listed zero-cell contrast DQB1*06:01 OR {fmt_num(dqb10601_within['or_haldane'])}, p={fmt_p(dqb10601_within['p_fisher'])}, FDR={fmt_p(dqb10601_within['fdr_bh'])}; all FDR=1.0 in table.",
            "source": "track4/T4_within_GSE286332_HT_vs_HTneg.tsv",
            "claim_grade": "no allele-level finding",
        },
        {
            "paper": "Paper 2 HT-PTC HLA/Immune",
            "id": "PTC-H4",
            "hypothesis": "HLA-I/II transcriptomic modules increase with the DM1/immune-active tumor state in TCGA-THCA.",
            "status": "strongly supported as expression-module biology",
            "key_result": f"DM1 vs HLA-I rho {fmt_num(dm1_hla1['spearman_rho'], 3)}, p={fmt_p(dm1_hla1['spearman_p'])}; DM1 vs HLA-II rho {fmt_num(dm1_hla2['spearman_rho'], 3)}, p={fmt_p(dm1_hla2['spearman_p'])}; n={dm1_hla2['n']}.",
            "source": "track5_summary.json",
            "claim_grade": "expression module, not allele genotype",
        },
        {
            "paper": "Paper 2 HT-PTC HLA/Immune",
            "id": "PTC-H5",
            "hypothesis": "The DM1-HLA relationship is not a purity/stromal artifact.",
            "status": "supported by partial correlations",
            "key_result": f"Partial immune+stromal rho: HLA-I {fmt_num(hla1_purity['rho_partial_immune+stromal'], 3)}, p={fmt_p(hla1_purity['p_partial_immune+stromal'])}; HLA-II {fmt_num(hla2_purity['rho_partial_immune+stromal'], 3)}, p={fmt_p(hla2_purity['p_partial_immune+stromal'])}.",
            "source": "track5_summary.json",
            "claim_grade": "robustness support",
        },
        {
            "paper": "Paper 2 HT-PTC HLA/Immune",
            "id": "PTC-H6",
            "hypothesis": "IFN-gamma is the dominant mediator connecting DM1 to antigen-presentation module induction.",
            "status": "supported observationally",
            "key_result": f"IFNG vs HLA-I rho {fmt_num(ifng_hla1['spearman_rho'], 3)}, p={fmt_p(ifng_hla1['spearman_p'])}; IFNG vs HLA-II rho {fmt_num(ifng_hla2['spearman_rho'], 3)}, p={fmt_p(ifng_hla2['spearman_p'])}; pct mediated Hallmark HLA-I {fmt_num(med_hla1['pct_mediated'])}%, HLA-II {fmt_num(med_hla2['pct_mediated'])}%.",
            "source": "track26_summary.json",
            "claim_grade": "observational mediation",
        },
        {
            "paper": "Paper 2 HT-PTC HLA/Immune",
            "id": "PTC-H7",
            "hypothesis": "HLA/IFN activation is inversely coupled to the thyroid RAI differentiation module.",
            "status": "supported as transcriptional decoupling",
            "key_result": f"DM1 vs RAI rho {fmt_num(rai_dm1['spearman_rho'], 3)}, p={fmt_p(rai_dm1['spearman_p'])}; HLA-II vs RAI rho {fmt_num(rai_hla2['spearman_rho'], 3)}, p={fmt_p(rai_hla2['spearman_p'])}.",
            "source": "track34/T02_headline_correlations.tsv",
            "claim_grade": "expression module",
        },
        {
            "paper": "Paper 2 HT-PTC HLA/Immune",
            "id": "PTC-H8",
            "hypothesis": "Across pan-cancer, DM1-HLA co-induction is a general lineage-portable pattern, with thyroid among the strongest lineages.",
            "status": "supported",
            "key_result": f"Pan-cancer n={track6['n_samples']}, lineages={track6['n_lineages']}; sign coherence HLA-I {fmt_num(track6['sign_coherence']['HLA-I']['fraction'] * 100)}%, HLA-II {fmt_num(track6['sign_coherence']['HLA-II']['fraction'] * 100)}%; thyroid HLA-II rho {fmt_num(track6['top5_sign_coherent_HLA2'][0]['rho'], 3)}.",
            "source": "track6_summary.json",
            "claim_grade": "contextual generalization",
        },
        {
            "paper": "Paper 2 HT-PTC HLA/Immune",
            "id": "PTC-H9",
            "hypothesis": "Single-cell/spatial data localize HLA-II to APC compartments and show tumor-stage acquisition, not just bulk contamination.",
            "status": "supported directionally",
            "key_result": f"scRNA cohorts: Pu n={track12['GSE184362_pu']['n_cells']} cells, Lu n={track12['GSE193581_lu']['n_cells']} cells; DC/Myeloid/B dominate HLA-II; Visium n={track33['n_samples_processed']} samples, niches={track33['n_total_niches']}.",
            "source": "track12_summary.json; track33_summary.json",
            "claim_grade": "cell-context support",
        },
        {
            "paper": "Paper 2 HT-PTC HLA/Immune",
            "id": "PTC-H10",
            "hypothesis": "TLS niches are HLA-II-hot in thyroid tumors, especially PTC/ATC.",
            "status": "strongly supported in spatial module data",
            "key_result": f"Niche-level global Cohen d {fmt_num(track33['niche_level_mean_cohen_d_global'], 3)}, mean delta {fmt_num(track33['niche_level_mean_diff_in_minus_out_global'], 3)}; stage d PTC {fmt_num(track33['stage_aggregate'][1]['mean_cohen_d'], 3)}, ATC {fmt_num(track33['stage_aggregate'][3]['mean_cohen_d'], 3)}.",
            "source": "track33_summary.json; track33/T6_stage_aggregate_effect.tsv",
            "claim_grade": "spatial support",
        },
        {
            "paper": "Paper 2 HT-PTC HLA/Immune",
            "id": "PTC-H11",
            "hypothesis": "Somatic loss-of-function in HLA presentation genes does not explain DM1/HLA state in TCGA-THCA.",
            "status": "supported negative",
            "key_result": f"Class-I LoF events {track25['hla_pathway_overall']['samples_any_classI_lof']}; class-I nonsilent samples {track25['hla_pathway_overall']['samples_any_classI_nonsilent']}; DM1 tertile class-I nonsilent p={fmt_p(track25['dm1_tertile_classI_nonsilent']['p_kruskal_3grp'])}.",
            "source": "track25_summary.json",
            "claim_grade": "negative mechanism",
        },
        {
            "paper": "Paper 2 HT-PTC HLA/Immune",
            "id": "PTC-H12",
            "hypothesis": "HLA-G rises with DM1 but is not a selective HLA-G immune-escape program.",
            "status": "supported",
            "key_result": f"HLA-G vs DM1 rho {fmt_num(hlaga['spearman_rho'], 3)}, p={fmt_p(hlaga['spearman_p'])}; HLA-G decoupling vs DM1 rho {fmt_num(hlag_decouple['spearman_rho_DM1_decoupling'], 3)}, p={fmt_p(hlag_decouple['spearman_p'])}.",
            "source": "track27_summary.json",
            "claim_grade": "anti-overclaim control",
        },
        {
            "paper": "Paper 2 HT-PTC HLA/Immune",
            "id": "PTC-H13",
            "hypothesis": "A BRAF-rich DM1-high/IFN-low fail-to-induce subgroup exists.",
            "status": "supported as subgroup, clinical outcome weak",
            "key_result": f"Fail-to-induce n={track28['n_fail_to_induce']}/{track28['n_total_THCA']}; BRAF OR {fmt_num(track28['fail_to_induce_top_drivers'][0]['OR'])}, FDR={fmt_p(track28['fail_to_induce_top_drivers'][0]['fdr'])}; survival not robust.",
            "source": "track28_summary.json",
            "claim_grade": "subgroup biology",
        },
        {
            "paper": "Paper 2 HT-PTC HLA/Immune",
            "id": "PTC-H14",
            "hypothesis": "DM1 has only a weak neoantigen/TMB link; HLA-II induction remains after neoantigen adjustment.",
            "status": "supported, but not central",
            "key_result": f"TMB vs DM1 rho {fmt_num(track29['tmb_dm1']['spearman_rho'], 3)}, p={fmt_p(track29['tmb_dm1']['p'])}; direct DM1->HLA-II beta {fmt_num(first_json_item(track29['mediation'], 'step', 'direct_DM1_to_HLA2')['beta_DM1'], 3)}, p={fmt_p(first_json_item(track29['mediation'], 'step', 'direct_DM1_to_HLA2')['p_DM1'])}.",
            "source": "track29_summary.json",
            "claim_grade": "supportive mechanism control",
        },
        {
            "paper": "Paper 2 HT-PTC HLA/Immune",
            "id": "PTC-H15",
            "hypothesis": "Thymic tolerance/AIRE and thyroid-antigen peptide binding provide a mechanistic hypothesis for AITD-HLA, but current data are not direct functional validation.",
            "status": "hypothesis-only",
            "key_result": f"HPA thymus used={track19['thymus_used_via_HPA']}; AIRE thymus rank={track19['AIRE_thymus_rank_HPA']}; GTEx thyroid donors={track19['n_gtex_thyroid_donors_used']}; peptidomics present={track19['track16_peptidomics_present']}.",
            "source": "track19_summary.json; track16 peptide tables",
            "claim_grade": "functional follow-up rationale",
        },
        {
            "paper": "Paper 2 HT-PTC HLA/Immune",
            "id": "PTC-H16",
            "hypothesis": "HLA-I gene-expression predicts ICI response modestly, but HLA-II does not add a strong response signal.",
            "status": "supported, external-cancer context only",
            "key_result": f"ICI cohorts n={track18['n_samples_with_response']}; HLA-I meta OR {fmt_num(track18['phase_C_v2_reconfirm']['HLA_I']['OR'])}, p={fmt_p(track18['phase_C_v2_reconfirm']['HLA_I']['p'])}; HLA-II meta OR {fmt_num(track18['new_HLA_II']['OR_meta'])}, p={fmt_p(track18['new_HLA_II']['p_meta'])}.",
            "source": "track18_summary.json",
            "claim_grade": "context only",
        },
        {
            "paper": "Paper 2 HT-PTC HLA/Immune",
            "id": "PTC-H17",
            "hypothesis": "Prospective HT-overlap PTC validation must use germline DNA NGS HLA and matched HT labels.",
            "status": "open; new data required",
            "key_result": f"Current explicit HT-overlap allele test is n={track4['n_HT_overlap_GSE286332']} vs n={track4['n_HT_negative_GSE286332']}; DRB1*04:05 OR=1.8 requires about {drb10405_power['n_per_group_required']} per group.",
            "source": "track4_summary.json; track4/T9_power_calc_top3_candidates.tsv",
            "claim_grade": "must-do validation",
        },
    ]

    write_tsv(OUT / "hypothesis_validation_matrix.tsv", hypotheses)

    priority_rows = [
        {
            "paper": "Paper 4 GD HLA",
            "submission_spine": "Pan-Asian Graves/AITD HLA architecture with Korean literature extraction, healthy Korean HLA atlas, cross-autoimmune specificity, and trans-ancestry MHC sensitivity.",
            "strongest_claim": "DPB1*05:01/B*46:01/A*02:07 are Pan-Asian GD/AITD alleles; C*01:02 is the most GD-specific focus allele.",
            "top_tier_gap": "No de novo Korean adult GD germline NGS case-control cohort yet.",
            "next_dataset": "Korean adult GD n>=250 cases + n>=500 matched controls; 8-locus NGS HLA.",
        },
        {
            "paper": "Paper 2 HT-PTC HLA/Immune",
            "submission_spine": "HT-overlap PTC as an antigen-presentation/TLS/IFN-gamma immune-context program, with HLA-II allele candidates explicitly kept exploratory.",
            "strongest_claim": "HLA-I/II gene-expression modules, IFN-gamma, spatial TLS niches, and APC compartments form a reproducible thyroid cancer immune-context axis.",
            "top_tier_gap": "Allele-level HT-PTC data are n=9 vs 9 and FDR-negative; no prospective germline NGS HT-overlap cohort.",
            "next_dataset": "PTC+HT vs PTC-only n>=600/group if testing DRB1*04:05 OR~1.8; smaller if DPB1-like high-frequency allele but current DPB1 result is null.",
        },
    ]
    write_tsv(OUT / "paper_priority_table.tsv", priority_rows)

    supported = [r for r in hypotheses if not r["status"].startswith(("open", "not validated")) and "not supported" not in r["status"]]
    open_or_negative = [r for r in hypotheses if r not in supported]

    md = f"""# Two HLA Paper NComm-Level Synthesis

**Generated:** 2026-05-09  
**Scope:** Paper 4 Graves/AITD HLA + Paper 2 HT-overlap PTC HLA/immune-context.  
**Hard rule:** cancer-cohort HLA allele genotype claims remain forbidden; HLA in cancer data means transcriptomic module only unless a new germline NGS cohort is assembled.

## Executive Verdict

There are two publishable stories, but they must be sold at different evidentiary levels.

1. **Paper 4, Graves/AITD HLA:** strongest near-manuscript core. The data support a Pan-Asian autoimmune-thyroid HLA architecture centered on DPB1*05:01, B*46:01, A*02:07 and C*01:02, with Korean literature extraction, healthy-HLA atlas support, cross-autoimmune specificity, and FinnGen/Pan-UKBB MHC sensitivity. The Nature Communications blocker is not analysis; it is the absence of a de novo Korean adult GD germline NGS case-control cohort.
2. **Paper 2, HT-overlap PTC HLA/immune:** top-tier potential only if framed as **immune-context biology**, not allele association. The strongest data are HLA-I/II expression modules, IFN-gamma mediation, TLS-HLA-II spatial niches, single-cell APC localization, and negative controls against somatic HLA-loss/HLA-G overclaim. Allele-level HT-PTC evidence is exploratory and underpowered.

## Paper 4 Claim Spine

- **Core thesis:** Pan-Asian Graves/AITD risk is shaped by a DPB1*05:01-centered class-II axis plus B*46:01/A*02:07/C*01:02 class-I linked architecture; C*01:02 is the most GD-specific focus allele, while DPB1*05:01 is AITD-broad.
- **Headline data:** DPB1*05:01 v3 pooled OR {fmt_num(v3['DPB1*05:01']['pooled_or'])} ({fmt_num(v3['DPB1*05:01']['ci_lo'])}-{fmt_num(v3['DPB1*05:01']['ci_hi'])}); B*46:01 OR {fmt_num(v3['B*46:01']['pooled_or'])} with high heterogeneity I2={fmt_num(v3['B*46:01']['I2_pct'])}%; A*02:07 upgraded to {ready['A*02:07']['grade_v3']}; C*01:02 specificity score {fmt_num(spec['C*01:02']['GD_specificity_score'])}.
- **Honest ceiling:** without Korean adult GD NGS cases and matched controls, the paper is an excellent synthesis/resource/meta-analysis, not a definitive Korean genetic association discovery.

## Paper 2 Claim Spine

- **Core thesis:** HT-overlap and immune-active PTC occupy an antigen-presentation/TLS/IFN-gamma state; HLA-II allele candidates are useful for prospective design but not currently validated.
- **Headline data:** TCGA-THCA DM1 vs HLA-II rho {fmt_num(dm1_hla2['spearman_rho'], 3)} (p={fmt_p(dm1_hla2['spearman_p'])}); IFN-gamma vs HLA-II rho {fmt_num(ifng_hla2['spearman_rho'], 3)}; spatial TLS-HLA-II global niche Cohen d {fmt_num(track33['niche_level_mean_cohen_d_global'], 3)} across {track33['n_total_niches']} niches; DRB1*04:05 HT-overlap candidate OR {fmt_num(drb10405['or_allele'])}, but FDR={fmt_p(drb10405['fdr_bh_allele'])}.
- **Honest ceiling:** the immune-context story is strong; the allele story is only a validation roadmap.

## Highest-Value Hypotheses

| Paper | ID | Status | Key Result |
|---|---|---|---|
"""
    for row in hypotheses:
        md += f"| {row['paper']} | {row['id']} | {row['status']} | {row['key_result']} |\n"

    md += f"""

## Nature Communications Upgrade Plan

| Paper | Add This Dataset | Why It Matters |
|---|---|---|
| Paper 4 GD HLA | Korean adult GD n>=250 cases + n>=500 matched controls, 8-locus germline NGS HLA | Converts meta-analysis/resource story into Korean-led original association paper; resolves DPB1*05:01 Korean-vs-Han effect and B*46:01 heterogeneity. |
| Paper 2 HT-PTC HLA/Immune | PTC+HT vs PTC-only cohort with explicit HT pathology/TPOAb/TgAb labels and germline NGS HLA | Tests DRB1*04:05 prospectively; current n=9 vs 9 is not publishable as allele association. |
| Paper 2 HT-PTC HLA/Immune | HT-PTC Visium or Xenium/CosMx with TLS markers + HLA-II module + B/T/myeloid labels | Turns current TLS-HLA-II spatial support into direct HT-overlap mechanism rather than general thyroid tumor immune context. |
| Both | Orthogonal HLA typing QC: arcasHLA vs NGS/SBT on a subset | Removes reviewer attack on RNA-seq HLA imputation and DQB1 subtype artifacts. |

## Claims To Avoid

- Do not say HLA alleles predict thyroid cancer risk or outcome from the current cancer cohorts.
- Do not join HLA allele tables to TCGA-THCA survival, RAI response, BRAF/RAS/TERT, DM1/DM2, stage, lymph node, or recurrence.
- Do not call DPB1*05:01 an HT-PTC allele finding; current HT-overlap data are null/weak for DPB1*05:01.
- Do not call DRB1*04:05 validated; it is a power-calibrated candidate.
- Do not sell HLA-G as selective immune escape in THCA; HLA-G rises proportionally with classical HLA-I.

## Outputs

- `hypothesis_validation_matrix.tsv` contains {len(hypotheses)} hypotheses with source-backed key results.
- `paper_priority_table.tsv` contains the two-paper claim spine and required validation datasets.

## Source Index

Primary source directory: `project/results/hla_deepdive_2026_05_08/`

Key files used: Track 1, 3, 4, 5, 6, 10, 11, 12, 18, 19, 21, 25, 26, 27, 28, 29, 31, 33, 34 summaries/tables.
"""

    report_path = OUT / "HLA_TWO_PAPER_NCOMM_SYNTHESIS.md"
    report_path.write_text(md)
    (REPORTS / "2026_05_09_HLA_TWO_PAPER_NCOMM_SYNTHESIS.md").write_text(md)

    print(report_path)
    print(OUT / "hypothesis_validation_matrix.tsv")
    print(REPORTS / "2026_05_09_HLA_TWO_PAPER_NCOMM_SYNTHESIS.md")


if __name__ == "__main__":
    main()
