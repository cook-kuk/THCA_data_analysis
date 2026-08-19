import React from "react";

/**
 * PresubInquiry — Nature Medicine (또는 NC) pre-submission inquiry 자료.
 * Cover letter + 200-word 요약 + editorial 전략 + 추천 reviewer.
 */

export const PresubInquiry: React.FC = () => (
  <div className="presub">

    {/* Cover */}
    <div className="presub-cover">
      <div className="presub-cover-badge">★ Editorial pre-submission inquiry template</div>
      <div className="presub-cover-title">Nature Medicine (primary) / Nature Communications (fallback) Presubmission Package</div>
      <div className="presub-cover-sub">
        <b>Editorial strategy:</b> Nature Medicine 정식 inquiry → 반응 강하면 direct submission · 약하면 NC 로 전환.
        분당 미팅 (2026-06-25) 후 <b>2 주 내</b> inquiry 발송 target.
      </div>
    </div>

    {/* Editorial strategy */}
    <div className="presub-strategy">
      <div className="presub-strategy-title">📊 Editorial Strategy · 3-tier journal ladder</div>
      <div className="presub-ladder">
        <div className="presub-tier best">
          <div className="presub-tier-name">Tier 1  ·  Nature Medicine</div>
          <div className="presub-tier-cond">
            <b>Condition:</b> A2 interaction (BRAF+ predictive) 결과 highlighted · IHC 3-plex clinical readiness 강조<br/>
            <b>Timing:</b> Pre-submission inquiry 우선 (200-word summary + cover letter)<br/>
            <b>Success signal:</b> Editor 가 "please submit" 응답 시
          </div>
        </div>
        <div className="presub-tier mid">
          <div className="presub-tier-name">Tier 2  ·  Nature Communications</div>
          <div className="presub-tier-cond">
            <b>Condition:</b> Nature Med editor 반응 약하거나 무응답 시 (3 주 후)<br/>
            <b>Timing:</b> Direct 정식 submission<br/>
            <b>Frame:</b> Discovery + multi-modality convergence + IHC translation
          </div>
        </div>
        <div className="presub-tier alt">
          <div className="presub-tier-name">Tier 3  ·  JCI Insight (fallback)</div>
          <div className="presub-tier-cond">
            <b>Condition:</b> NC reject 시 · rebuttal 어려운 methodology critique 발생 시<br/>
            <b>Timing:</b> 6 개월 이내<br/>
            <b>Frame:</b> Translational focus 강조 · IHC readiness 중심
          </div>
        </div>
      </div>
    </div>

    {/* 200-word summary */}
    <div className="presub-block">
      <div className="presub-block-title">1 · 200-word Pre-submission Summary (Nature Medicine)</div>
      <div className="presub-body">
        <p>
          Radioiodine is central to differentiated thyroid cancer management, yet a substantial subset of BRAF V600E-mutant papillary thyroid cancer patients receive repeated high-dose radioiodine despite molecularly compromised iodine-handling biology, incurring cumulative toxicity without therapeutic benefit. We define a compact eight-gene thyroid differentiation and iodine-handling axis (TG, TPO, TSHR, SLC5A5, DIO1, PAX8, NKX2-1, FOXE1) and identify a significant driver-specific interaction between this axis and BRAF V600E status on progression-free interval in TCGA-THCA (Cox interaction p = 0.022; BRAF+ subset log-rank p = 7.6 × 10⁻³, n = 287). A three-marker immunohistochemistry combination (TG + PAX8 + TTF-1) — already in routine clinical use across thyroid pathology laboratories — reproduces this stratification with log-rank p = 1.7 × 10⁻⁴, positioning the finding for immediate translational deployment without new assay development. The axis replicates across five external cohorts including Lee 2024 (Korean n = 370, Mann-Whitney p = 2.7 × 10⁻⁷). Monte Carlo simulation using TCGA-derived priors predicts &gt;90% power for prospective replication at n = 200. We propose this observation as a candidate treatment-selection biomarker for the BRAF V600E-mutant subset, with a defined path from retrospective institutional validation through a planned prospective cohort at Seoul National University Bundang Hospital.
        </p>
      </div>
      <div className="presub-wc">Word count: ~ 200 · Fits typical pre-submission inquiry length</div>
    </div>

    {/* Cover letter */}
    <div className="presub-block">
      <div className="presub-block-title">2 · Cover Letter Draft</div>
      <div className="presub-body">
        <p><em>[Editor name],</em></p>
        <p><em>[Nature Medicine editorial address]</em></p>
        <p>Dear Dr. [Editor name],</p>

        <p className="presub-para-hook">
          <b>[AUTHOR HOOK SLOT — protected voice section.]</b>{" "}
          Insert opening paragraph in author's voice. Recommended: start from the clinical paradox that most differentiated thyroid cancers have excellent survival, yet a subset receives repeated radioiodine despite molecularly poor iodine-handling biology — with cumulative toxicity but no therapeutic gain. Frame the paper as identifying, for the first time in a defined driver subset (BRAF V600E+), a routine-immunohistochemistry-deployable classifier that can pre-emptively identify this subset.
        </p>

        <p>
          We submit for pre-submission consideration a manuscript that defines a compact eight-gene thyroid differentiation and iodine-handling axis and identifies a driver-specific interaction with BRAF V600E on progression-free interval that, if confirmed prospectively, would elevate the axis from a prognostic risk-stratification tool to a candidate treatment-selection biomarker.
        </p>

        <p>
          The manuscript's key advances relative to prior work — including the Cancer Genome Atlas Research Network (2014 Cell), Landa et al. (2016 JCI), and Chakravarty et al. (2011 JCI) / Ho et al. (2013 NEJM) / Rothenberg et al. (2015 CCR) — are three-fold.  First, we quantify a driver-specific predictive interaction (DM1 × BRAF V600E, Cox interaction p = 0.022, n = 472), whereas prior work established the biology of MAPK-driven lineage silencing and the clinical possibility of radioiodine restoration but did not identify a pre-treatment classifier for the responder subset.  Second, we demonstrate that a three-marker immunohistochemistry combination corresponding to routine thyroid diagnostic pathology (TG + PAX8 + TTF-1) reproduces the BRAF+ progression-free interval stratification with log-rank p = 1.7 × 10⁻⁴, providing an immediate translational deployment path that does not require new assay development or reagent qualification.  Third, we replicate the underlying axis across five external cohorts, most strongly in a Korean papillary thyroid cancer cohort (Lee 2024, n = 370 tumours, Mann-Whitney p = 2.7 × 10⁻⁷), establishing generalizability beyond the discovery cohort.
        </p>

        <p>
          We believe the work fits Nature Medicine because it satisfies four criteria that distinguish translational discovery from prognostic biomarker development: (i) driver-specific predictive interaction (BRAF+-restricted DM1 effect), (ii) immediate clinical deployment path (routine IHC 3-plex), (iii) external multi-population replication (Korean + microarray + proteomic + single-cell), and (iv) a defined prospective validation plan at n = 200 with pre-computed &gt;90% power under TCGA-derived priors.  The clinical actionability is not aspirational: the reflex algorithm we describe can be tested tomorrow at any institution with routine thyroid diagnostic pathology and BRAF genotyping.
        </p>

        <p>
          The full analysis package includes eight main figures, extended data figures, complete methods, source data, and a public reproducibility repository.  All source datasets are publicly available.  We would welcome your assessment of whether the manuscript falls within Nature Medicine's scope, and we would be grateful for guidance on any additional analyses that would strengthen the pre-submission position.
        </p>

        <p>Sincerely,</p>

        <p>
          Seungho Cook, Ph.D.<br/>
          Hyeong-Won Yu, M.D., Ph.D.<br/>
          Corresponding authors<br/>
          Seoul National University Bundang Hospital<br/>
          [institutional email · phone]
        </p>
      </div>
    </div>

    {/* Suggested reviewers */}
    <div className="presub-block">
      <div className="presub-block-title">3 · Suggested Reviewers  (author-suggested — editor may or may not use)</div>
      <div className="presub-reviewer-grid">
        <div className="presub-reviewer">
          <div className="presub-rev-role">Molecular biology · Thyroid TF / lineage</div>
          <div className="presub-rev-hint">Author of foundational thyroid lineage TF work (PAX8 / NKX2-1) — thyroid organogenesis + differentiation programme expertise</div>
          <div className="presub-rev-conflict"><b>Non-conflict check:</b> No shared institution · No co-authorship in last 5 years</div>
        </div>
        <div className="presub-reviewer">
          <div className="presub-rev-role">Genomics / Multi-omics thyroid</div>
          <div className="presub-rev-hint">Corresponding author of TCGA-THCA 2014 or Landa 2016 JCI — molecular thyroid oncology genomics leadership</div>
          <div className="presub-rev-conflict"><b>Non-conflict check:</b> Applies</div>
        </div>
        <div className="presub-reviewer">
          <div className="presub-rev-role">Radioiodine biology / Redifferentiation</div>
          <div className="presub-rev-hint">Author of MAPK inhibitor radioiodine restoration work (Chakravarty · Ho · Rothenberg lineage) — mechanistic biology + clinical translation</div>
          <div className="presub-rev-conflict"><b>Non-conflict check:</b> Applies</div>
        </div>
        <div className="presub-reviewer">
          <div className="presub-rev-role">Clinical thyroid oncology</div>
          <div className="presub-rev-hint">Thyroid cancer guideline committee member (ATA / European Thyroid Association) — clinical decision-making expertise</div>
          <div className="presub-rev-conflict"><b>Non-conflict check:</b> Applies</div>
        </div>
      </div>
      <div className="presub-note">
        <b>Note:</b> 실제 name 은 authors 가 최종 review 후 submission 시 fill.  Editorial guidelines 준수 (author-suggested reviewer 는 editor 재량).
      </div>
    </div>

    {/* Novelty highlights */}
    <div className="presub-block">
      <div className="presub-block-title">4 · Novelty Highlights  (editor 눈에 즉시 들어와야 하는 3 개)</div>
      <div className="presub-highlights">
        <div className="presub-highlight best">
          <div className="presub-h-num">01</div>
          <div className="presub-h-body">
            <b>Predictive interaction with BRAF V600E</b> — DM1 × BRAF Cox interaction p = 0.022.
            Prior work established prognostic axes but did not identify a driver-context-specific predictive classifier.
            This is the paper's <em>core novelty</em>.
          </div>
        </div>
        <div className="presub-highlight best">
          <div className="presub-h-num">02</div>
          <div className="presub-h-body">
            <b>Immediate clinical translation via routine IHC 3-plex</b> — TG + PAX8 + TTF-1 reproduce BRAF+ PFI stratification (log-rank p = 1.7 × 10⁻⁴).
            No new assay development required; deployment timeline is measured in months, not years.
          </div>
        </div>
        <div className="presub-highlight">
          <div className="presub-h-num">03</div>
          <div className="presub-h-body">
            <b>Cross-population multi-omics generalizability</b> — Lee 2024 (Korean n = 370, p = 2.7 × 10⁻⁷) + GPL570 microarray + Mun 2025 proteomics (7/7 concordant) + single-cell thyrocyte-restricted concordance.
            5-modality Spearman ρ = 0.71-0.94.
          </div>
        </div>
      </div>
    </div>

    {/* Data + reproducibility commitment */}
    <div className="presub-block">
      <div className="presub-block-title">5 · Data · Code · Reproducibility Commitment</div>
      <ul className="presub-ul">
        <li><b>Data availability:</b> All source cohorts publicly available (TCGA, cBioPortal, GEO, ENA) — accession list attached as Supplementary Table 1.</li>
        <li><b>Code availability:</b> Complete analysis and figure-generation scripts to be archived at Zenodo with DOI at submission.</li>
        <li><b>Processed data:</b> Per-sample scores, methylation summaries, interaction Cox outputs, Monte Carlo seeds to be deposited.</li>
        <li><b>Reproducibility statement:</b> All figures regenerable from raw public data + our deposited processing scripts.</li>
        <li><b>Author contribution transparency:</b> S.C. computational work; H.-W.Y. clinical framing; independent code review by [collaborator TBD].</li>
      </ul>
    </div>

    {/* Timeline */}
    <div className="presub-block">
      <div className="presub-block-title">6 · Timeline · 미팅 후 3 개월</div>
      <table className="presub-tbl">
        <thead>
          <tr><th>Week</th><th>Action</th><th>Owner</th><th>Deliverable</th></tr>
        </thead>
        <tbody>
          <tr><td>W 1</td><td>미팅 후 rebuttal 검토 · 200-word 최종 확정</td><td>Cook · Yu · Kang</td><td>Final summary + cover letter</td></tr>
          <tr><td>W 2</td><td>Nature Medicine editor 이메일 (pre-submission inquiry)</td><td>Yu (corresponding)</td><td>Editor 응답 대기</td></tr>
          <tr><td>W 3-4</td><td>Editor 응답 대기 · 병렬로 NC direct submission 준비</td><td>Cook · Yu</td><td>NC submission ready</td></tr>
          <tr><td>W 5 (case A · Nature Med 긍정)</td><td>Nature Medicine 정식 submission</td><td>Cook · Yu</td><td>Submission ID</td></tr>
          <tr><td>W 5 (case B · Nature Med 무응답 / 부정)</td><td>NC direct submission</td><td>Cook · Yu</td><td>Submission ID</td></tr>
          <tr><td>W 6-10</td><td>Review 대기 (typical Nature Medicine ~ 8 주 · NC ~ 6 주)</td><td>—</td><td>Reviewer comments</td></tr>
          <tr><td>W 11-12</td><td>Rebuttal 작성 · 이 웹의 Reviewer Q&A 시뮬레이션 참조</td><td>Cook · Yu</td><td>Rebuttal letter + revised MS</td></tr>
        </tbody>
      </table>
    </div>

    {/* Final risk */}
    <div className="presub-final">
      <div className="presub-final-title">🎯 Journal Decision Risk Matrix</div>
      <div className="presub-risk-grid">
        <div className="presub-risk-card best">
          <div className="presub-risk-name">Nature Medicine Accept</div>
          <div className="presub-risk-prob">~ 15-25%</div>
          <div className="presub-risk-cond">A2 interaction + IHC translation + Lee replication 3 축이 강하면.  Editor 의 clinical 판단 결정적.</div>
        </div>
        <div className="presub-risk-card mid">
          <div className="presub-risk-name">Nature Comm Accept</div>
          <div className="presub-risk-prob">~ 40-55%</div>
          <div className="presub-risk-cond">Multi-modality convergence + external replication + IHC pathway 로 realistic target.  주력 target.</div>
        </div>
        <div className="presub-risk-card low">
          <div className="presub-risk-name">JCI Insight Accept</div>
          <div className="presub-risk-prob">~ 65-80%</div>
          <div className="presub-risk-cond">NC reject 시 fallback.  translational thyroid biomarker 로 fit 매우 좋음.</div>
        </div>
      </div>
      <div className="presub-final-note">
        <b>결정적 upside:</b> 만약 분당 pilot 30 (IHC handoff §1-8) 이 6 개월 내 성공하면 revised MS 에 <em>"prospective institutional pilot confirms retrospective finding"</em> 문장 추가 → Nature Medicine 상방 재도전.
      </div>
    </div>

  </div>
);
