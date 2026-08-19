import React from "react";

/* ================================================================
   Verification log — every key claim with verification status,
   source file, audit date. For PI review before external sharing.
   ================================================================ */

type Row = {
  claim: string;
  value: string;
  source: string;
  status: "verified" | "audit-locked" | "needs-PI-check" | "illustrative";
  date: string;
  note?: string;
};

const ROWS: Row[] = [
  {
    claim: "TCGA-THCA primary PTC sample size",
    value: "n = 504",
    source: "GDC TCGA-THCA harmonized RNA-seq · STAR-aligned RSEM TPM",
    status: "verified",
    date: "2026-05-08",
    note: "post-QC n = 504; d4p2 master table"
  },
  {
    claim: "DM1 prevalence in TCGA",
    value: "28.4 % (140 / 504)",
    source: "8-gene KMeans k=2 · d4p2_tcga_hashimoto_signature/tcga_with_clinical_mutations.tsv",
    status: "audit-locked",
    date: "2026-05-08"
  },
  {
    claim: "Pan-genome ARI",
    value: "0.92 (top-5000 MAD vs 8-gene reference)",
    source: "p4_pangenome_vs_tiera67/ari_comparison.tsv",
    status: "audit-locked",
    date: "2026-05-08",
    note: "Hypergeometric TIERA67 ∩ top-100 MAD p = 3 × 10⁻⁴"
  },
  {
    claim: "BRAF V600E single-feature AUC vs DM1",
    value: "AUC = 0.602  ·  d = −0.044  ·  MW p = 0.57",
    source: "p1_driver_mrna_audit/driver_neutrality.tsv",
    status: "audit-locked",
    date: "2026-04-29"
  },
  {
    claim: "DM1 prevalence in BRAF V600E+",
    value: "0.7 %  (2 / 273)",
    source: "d4p2 master + cBioPortal MAF",
    status: "audit-locked",
    date: "2026-05-20"
  },
  {
    claim: "DM1 prevalence in RAS-mutant",
    value: "96.4 %  (53 / 55)",
    source: "d4p2 master + cBioPortal MAF",
    status: "audit-locked",
    date: "2026-05-20"
  },
  {
    claim: "DM1 prevalence in BRAF/RAS-negative",
    value: "49.1 %  (85 / 173)",
    source: "d4p2 master + cBioPortal MAF",
    status: "audit-locked",
    date: "2026-05-20"
  },
  {
    claim: "TPO HM450 promoter β Cohen's d (DM1 vs DM2)",
    value: "d = 2.30  ·  p = 1.9 × 10⁻¹⁸",
    source: "audit_2026_04_30/round5/r5_2_per_gene_methylation_DM.tsv",
    status: "audit-locked",
    date: "2026-04-30"
  },
  {
    claim: "Mean 8-gene β: DM1 vs DM2",
    value: "0.385 vs 0.253  (+52 %)",
    source: "audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv",
    status: "audit-locked",
    date: "2026-04-30"
  },
  {
    claim: "DM1 fusion+ enrichment (Fisher OR)",
    value: "OR = 7.41  [4.38, 12.55]  ·  p = 1.9 × 10⁻¹³",
    source: "audit_2026_04_30/round3/cbio_sv_thca.tsv",
    status: "audit-locked",
    date: "2026-04-30",
    note: "DM1 fusion+ 76.8 % (63/82) vs DM2 30.9 % · chi² MAR p = 0.56"
  },
  {
    claim: "FVPTC enrichment (Fisher OR)",
    value: "OR = 17.9  ·  p = 3 × 10⁻³¹",
    source: "audit_2026_04_30/round3/round3_FVPTC_OR.tsv",
    status: "needs-PI-check",
    date: "2026-05-08",
    note: "TCGA pathology call 기준; 95 % CI 재계산 필요"
  },
  {
    claim: "Pooled OS HR (TCGA + MSK)",
    value: "HR = 2.53  [1.31, 4.89]  ·  I² = 0 %",
    source: "audit_2026_04_30/p6_per_study_hr.tsv + DerSimonian-Laird random-effects",
    status: "audit-locked",
    date: "2026-04-30",
    note: "TCGA HR 2.30 [0.77, 6.88] 단독 NS · MSK 가 견인"
  },
  {
    claim: "FFPE vs FF Kolmogorov-Smirnov",
    value: "KS p = 0.44  (no detectable shift)",
    source: "GSE213647 FFPE n=632 vs TCGA FF n=504",
    status: "audit-locked",
    date: "2026-04-30"
  },
  {
    claim: "Master cross-cohort forest mean Cohen's d",
    value: "mean 2.81  ·  median 2.37  ·  all 14 entries ≥ 1.56",
    source: "manuscript_biorxiv_2026_05_20/assets/fig_cross_cohort_master_forest.png + underlying TSV",
    status: "audit-locked",
    date: "2026-05-21"
  },
  {
    claim: "Per-gene × cohort × contrast matrix direction-consistency",
    value: "80 / 80 cell direction-consistent",
    source: "manuscript_biorxiv_2026_05_20/assets/fig_pergene_xcohort_matrix.png + underlying TSV",
    status: "audit-locked",
    date: "2026-05-21",
    note: "8 genes × 10 contrast × 4 cohort"
  },
  {
    claim: "Lee 2024 (GSE213647) within-cohort KMeans d",
    value: "all-cohort d = 5.93  ·  tumour-only n=370 d = 5.48",
    source: "v17_korean/GSE213647_panel_score.tsv + within-cohort KMeans",
    status: "audit-locked",
    date: "2026-05-21"
  },
  {
    claim: "K2 (PRJEB11591) within-cohort KMeans d",
    value: "d = 1.94  ·  DM1 prevalence 60.8 %",
    source: "v17_korean K2 within-cohort unsupervised analysis",
    status: "audit-locked",
    date: "2026-05-21",
    note: "transfer label AUC 0.275 (calibration mismatch) vs continuous AUC 0.883"
  },
  {
    claim: "Mun 2025 proteogenomic thyroid_diff d",
    value: "d = −1.91  ·  7 / 7 panel proteins sign-consistent (NKX2-1 absent)",
    source: "proteogenomic_v1/paper3_mun2025_dediff_layer/module_dediff_trend.tsv",
    status: "audit-locked",
    date: "2026-05-08",
    note: "n = 336 · within-cohort k=2 KMeans Cohen's d = 2.86"
  },
  {
    claim: "Landa 2016 silenced-list overlap with 8-gene panel",
    value: "5 / 8 overlap  (TG · TSHR · TPO · PAX8 · DIO1)",
    source: "Landa I et al. 2016 J Clin Invest 126(3):1052-1066 · supplementary list",
    status: "verified",
    date: "2026-05-04",
    note: "Independent design — panel 정의 후 cross-reference (reverse-causality lock)"
  },
  {
    claim: "GPL570 4-cohort Spearman ρ (DM1 × nonoverlap)",
    value: "4 / 4 cohorts ρ ≤ −0.84",
    source: "GSE33630 ρ=−0.93 · GSE29265 ρ=−0.85 · GSE65144 ρ=−0.94 · GSE53157 ρ=−0.84",
    status: "audit-locked",
    date: "2026-05-20"
  },
  {
    claim: "GSE151179 pre vs post-RAI thyroid_diff d",
    value: "d ≈ −1.01  ·  MW p ≈ 1 × 10⁻⁴",
    source: "aggressive_sprint_2026_05_06/gse151179_rai_scores.tsv (pre n=35 / post n=17)",
    status: "audit-locked",
    date: "2026-05-06",
    note: "sample size 작음 · CI 넓음"
  },
  {
    claim: "Lu 2023 single-cell thyrocyte cell count",
    value: "n = 14,624 cells (KRT8 ∩ KRT19 ∩ EPCAM filter)",
    source: "v17_lu2023/GSE193581_cell_panel_score.tsv",
    status: "verified",
    date: "2026-05-08"
  },
  {
    claim: "Pu 2021 per-patient Spearman r range",
    value: "r = 0.798–0.886  ·  all p < 10⁻¹⁰ (Bonferroni)",
    source: "audit_2026_04_30/p5_lu2023_per_patient_r.tsv (also used for Pu 2021 GSE184362)",
    status: "audit-locked",
    date: "2026-04-30"
  },
  {
    claim: "Selpercatinib eligible population estimate",
    value: "≈ 48 / 1,000 PTC patients",
    source: "TCGA RET-fusion 6.5 % × DM1 capture 81.8 % × FDA approval criteria",
    status: "illustrative",
    date: "2026-05-04",
    note: "population estimate · prospective trial 결과 아님"
  },
  {
    claim: "Permutation AUC null",
    value: "empirical p = 0.0002 (1 / 5,000 shuffles)",
    source: "audit_2026_04_30/permutation_null/perm_AUC_null.tsv",
    status: "audit-locked",
    date: "2026-04-30"
  },
  {
    claim: "n−1 knock-out floor (LOO sensitivity)",
    value: "min Cohen's d = 1.54 (DIO1 제외 시)  ·  8 / 8 d > 1.5",
    source: "audit_2026_04_30/knock_out_n_minus_1.tsv",
    status: "audit-locked",
    date: "2026-04-30"
  },
  {
    claim: "Cluster stability bootstrap (TCGA k=2)",
    value: "97.4 % tumors ≥ 90 % bootstrap iteration agreement",
    source: "500 × bootstrap re-clustering",
    status: "audit-locked",
    date: "2026-04-30"
  }
];

const STATUS_LABEL: Record<Row["status"], string> = {
  "verified":          "✓ Verified (external source)",
  "audit-locked":      "✓ Audit-locked (manuscript v8)",
  "needs-PI-check":    "⚠ Needs PI re-check",
  "illustrative":      "○ Illustrative estimate"
};

export const VerificationLog: React.FC = () => (
  <div>
    <table className="verif-table">
      <thead>
        <tr>
          <th>Claim</th>
          <th>Value</th>
          <th>Source</th>
          <th>Status</th>
          <th>Audit date</th>
        </tr>
      </thead>
      <tbody>
        {ROWS.map((r, i) => (
          <tr key={i} className={`verif-${r.status}`}>
            <td className="verif-claim">{r.claim}</td>
            <td className="verif-value">{r.value}</td>
            <td className="verif-source"><code>{r.source}</code></td>
            <td className="verif-status">{STATUS_LABEL[r.status]}</td>
            <td className="verif-date">{r.date}</td>
          </tr>
        ))}
      </tbody>
    </table>
    <div className="verif-summary">
      <span className="verif-pill verified">✓ Verified  {ROWS.filter(r => r.status === "verified").length}</span>
      <span className="verif-pill audit-locked">✓ Audit-locked  {ROWS.filter(r => r.status === "audit-locked").length}</span>
      <span className="verif-pill needs-PI-check">⚠ Needs PI check  {ROWS.filter(r => r.status === "needs-PI-check").length}</span>
      <span className="verif-pill illustrative">○ Illustrative  {ROWS.filter(r => r.status === "illustrative").length}</span>
      <span className="verif-total">Total <b>{ROWS.length}</b> tracked claims</span>
    </div>
    <div className="verif-foot">
      모든 수치는 manuscript audit log 와 cross-checked.  ⚠ 표시는 PI 검토 전 외부 공유 보류.
      Audit-locked = manuscript_v8/14_self_verification_report.md 참조.
    </div>
  </div>
);
