# Synthetic lethality 정리 — 2026-05-09

## Bottom line

현재 로컬 근거로 가능한 가장 강한 표현은 **validated synthetic lethality**가 아니라 **DM1/lineage-silenced state에서의 perturbation vulnerability 후보 우선순위화**입니다. Paper 1에는 reviewer-reserve/actionability block으로만 두고, 본문 핵심 claim으로 올리면 과합니다.

## Claim boundary

- **Allowed:** pharmacologic vulnerability, CRISPR dependency prioritization, wet-lab perturbation roadmap, reviewer-reserve actionability.
- **Forbidden:** validated synthetic lethality, clinical treatment recommendation, patient-selection biomarker, treatment-induced thyroid-gene restoration.
- **Best one-line wording:** DM1-high or lineage-silenced models show public-data vulnerability signals, especially MAPK-axis PRISM sensitivity and MYC/NAMPT dependency, motivating focused perturbation experiments rather than establishing synthetic lethality.

## Highest-signal facts

- PRISM: 1518 drugs, FDR<0.05 hits 11; canonical MAPK hits 7/11; top 7 by FDR are 7/7 canonical MAPK.
- Top PRISM hit: AZD-0364, Cohen's d=-0.594, delta LFC=-0.937, FDR=5.88e-07.
- Continuous PRISM MAPK score: DM1 score vs canonical MAPK-inhibitor mean LFC rho=-0.236, p=3.31e-10, n=690; excluding thyroid rho=-0.240, p=2.27e-10.
- DepMap CRISPR: MYC d=-0.499, p=1.00e-11; NAMPT d=-0.445, p=1.53e-09; n_high=372, n_low=394.
- Paper9 wet-lab queue: SLC1A5 > GLUD1 > GLS. SLC1A5 rho=0.102, p=5.33e-04; GLUD1 rho=0.075, p=0.011; GLS rho=0.088, p=0.003.

## Candidate matrix

| tier | candidate_program | perturbation | headline_metrics | disposition |
| --- | --- | --- | --- | --- |
| T1 | DM1-high / MAPK-axis pharmacologic vulnerability | MEK/RAF/ERK inhibitors | PRISM 1518 drugs; FDR<0.05 hits 11, canonical MAPK 7/11; top7 MAPK 7/7; AZD-0364 d=-0.594, FDR=5.88e-07; DM1 score x MAPK LFC rho=-0.236, p=3.31e-10 | KEEP AS REVIEWER-RESERVE, NOT MAIN SYNTHETIC-LETHAL CLAIM |
| T2 | DM1-high CRISPR dependency | MYC and NAMPT genetic perturbation; inhibitor follow-up only as separate validation | n_high=372, n_low=394; MYC d=-0.499, p=1.00e-11; NAMPT d=-0.445, p=1.53e-09 | KEEP AS FUNCTIONAL PRIORITIZATION; MYC HAS PAN-ESSENTIAL RISK |
| T3A | Glutamine transporter vulnerability | SLC1A5 KD/KO/CRISPRi; glutamine withdrawal; V-9302/GPNA only if metadata supports | SLC1A5 dependency rho=0.102, p=5.33e-04, n=1140; pan-cancer corrected rho=0.078, p=0.008; strict thyroid n=11 NS | TOP PAPER9 WET-LAB PRIORITY |
| T3B | Glutamate/TCA anaplerosis vulnerability | GLUD1 KD/KO/CRISPRi; inhibitor only if specificity acceptable | GLUD1 dependency rho=0.075, p=0.011; pan-cancer corrected rho=0.062, p=0.037; strict thyroid n=11 NS | SECOND PAPER9 WET-LAB PRIORITY |
| T3C | Glutaminase vulnerability | GLS KD/KO/CRISPRi; BPTES/CB-839 style pharmacology only with caveat | GLS dependency rho=0.088, p=0.003; pan-cancer corrected rho=0.081, p=0.006; strict thyroid n=11 NS | CANDIDATE, NOT HEADLINE |
| BOUNDARY | Spatial MAPK-panel test | None; observational spatial stress test | v16/v17 did not rescue spatial anti-correlation; 0/45 subset/adjustment tests negative in v16; full+detection residual rho ~0 in v17. | USE AS HONEST NEGATIVE CAVEAT ONLY |
| BOUNDARY | ICI vulnerability / immune phenotype | Immune checkpoint therapy context | Project boundary map already separates immune vulnerability from synthetic-lethal and synthetic-rescue hypotheses. | KEEP OUT OF SYNTHETIC-LETHALITY CLAIM SET |

## Deployment rule

Use this as a **factual data block / dossier**, not as voice-protected Discussion, Limitations, Cover Letter, Hook, Aim, or Q9 prose.
