# v13 Analysis Summary

## Scope

v13 extends the thyroid cancer pipeline from expression/subtype analysis into an in-silico therapeutic landscape layer. The objective was not to claim validated drug efficacy, but to rank practical next steps across repurposing, structure-informed triage, and modality fit for the eight BRAF-like targets carried forward from the earlier pipeline.

Targets:

- `CYP1B1`
- `TACSTD2`
- `TMPRSS4`
- `PLEKHA6`
- `LDLR`
- `GABRB2`
- `B3GNT3`
- `PTPRE`

## What Was Actually Run

1. Structure acquisition
   Experimental PDB/cryo-EM when available, AlphaFold DB fallback otherwise.

2. Pocket detection
   In-repo alpha-sphere-lite heuristic plus HETATM ligand rescue.

3. Virtual screening
   Similarity-based FDA-library screen against curated literature binders, not full GNINA/Vina docking.

4. Multi-modality matrix
   Biology-informed modality scoring across 10 therapeutic formats, augmented with ClinicalTrials.gov evidence retrieval.

5. TROP2 ADC review
   Dedicated thyroid/TROP2 literature and trial pull.

6. ADMET triage
   RDKit proxy model with rule-based toxicity flags.

7. Composite prioritization
   Binding + ADMET + maturity + thyroid evidence + modality fit + repurposing bias.

## Key Output Counts

- Structures: `8`
  - PDB/cryo-EM: `5`
  - AF2 fallback: `3`
- Virtual screening library: `2443` FDA-phase small molecules
- Screened target space: `2443 × 8`
- Modality matrix: `8 × 10 = 80` cells
- TROP2 literature records: `60`
- TROP2 trial records indexed: `194`
- Thyroid-relevant TROP2 trials: `2`
- Ranked candidates: `400`
- Tier A: `59`
- Tier B: `283`
- Tier C: `56`
- Tier D: `2`

## Structure Layer

Best-pocket PLB and preferred modality summary:

| Target | Best PLB | Tier | Preferred modalities |
|---|---:|---|---|
| `CYP1B1` | 0.850 | High | `SM` |
| `LDLR` | 0.700 | High | `mAb, ADC, siRNA, mRNA_replacement` |
| `TACSTD2` | 0.833 | High | `ADC, mAb, bispecific, CAR_T, radioligand` |
| `TMPRSS4` | 0.507 | Medium | `SM, mAb` |
| `GABRB2` | 0.887 | High | `SM` |
| `PLEKHA6` | 0.807 | High | `PROTAC, siRNA, ASO` |
| `PTPRE` | 0.392 | Low | `SM, PROTAC, siRNA` |
| `B3GNT3` | 0.575 | Medium | `SM, siRNA, ASO` |

Interpretation:

- `TACSTD2`, `LDLR`, `CYP1B1`, `GABRB2`, and `PLEKHA6` look actionable in principle.
- `PTPRE` remains the weakest small-molecule target in this run because the top pocket score stayed below 0.4.
- `PLEKHA6` stands out as a better degraders/oligo target than a classic repurposing target.

## Virtual Screening

Top-ranked examples by target:

- `GABRB2`: `propofol`, `etomidate`, benzodiazepine-class chemistry
- `LDLR`: `simvastatin`, `lovastatin`, `pravastatin`
- `CYP1B1`: `dronabinol`, `hydroquinone`, flavonoid-like scaffolds
- `TMPRSS4`: `lidocaine`, `acetaminophen`, `clofibrate`
- `PLEKHA6`: `thalidomide`
- `PTPRE`: `indomethacin`
- `B3GNT3`: `fludeoxyglucose F18`, `dextrose`, `ribavirin`
- `TACSTD2`: `topotecan`, `topotecan hydrochloride`, `estrone`

Important limitation:

- These are similarity-screen outputs, not docking scores and not biochemical validation.
- `TACSTD2` is biologically strongest as an ADC/mAb-style target, so its small-molecule list should not be over-read.

## Modality Landscape

Strongest immediate categories:

- `TACSTD2`: ADC / mAb / bispecific / CAR-T / radioligand
- `LDLR`: mAb / siRNA / metabolic intervention / indirect repurposing
- `CYP1B1`: small molecule
- `TMPRSS4`: small molecule + exploratory mAb
- `PLEKHA6`: PROTAC / siRNA / ASO
- `B3GNT3`: siRNA / ASO

Counts by modality evidence class show:

- approved evidence exists for `SM`, `mAb`, `ADC`, `siRNA`, `ASO`, `radioligand`
- hypothesis-only cells are concentrated in `PROTAC`, `ADC`, `bispecific`
- `siRNA` is the broadest non-small-molecule exploratory class across these targets

## TROP2 / TACSTD2 Summary

This remains the highest-value translational finding in v13.

What is supported:

- `TACSTD2` is structurally and biologically compatible with surface-targeted modalities.
- The modality matrix supports ADC as the most mature thyroid-relevant actionable path.
- Literature retrieval found `60` TROP2/thyroid-related PubMed records.
- ClinicalTrials.gov retrieval found `194` total TROP2-ADC-related studies and `2` thyroid-relevant trials.

What is still inferential:

- The statement that BRAF-like THCA is the best enrichment cohort is a biomarker-driven translational hypothesis, not prospective clinical proof.
- The actual response prediction requires archived-tissue correlation or new prospective biomarker gating.

## Prioritization Readout

Tier A count by target:

- `LDLR`: `16`
- `CYP1B1`: `15`
- `GABRB2`: `13`
- `TMPRSS4`: `8`
- `B3GNT3`: `5`
- `PLEKHA6`: `1`
- `PTPRE`: `1`

Notably:

- `TACSTD2` did not produce Tier A small-molecule hits, which is biologically consistent with the claim that it is primarily an ADC-class opportunity.
- The strongest global Tier A scores were dominated by CNS/lipid/small-molecule repurposing chemistry rather than ADC payload logic.

## What Is Reliable vs What Is Not

More reliable:

- Structure source provenance
- Pocket tier ranking at the coarse triage level
- Modality fit by target biology
- Existence of TROP2 thyroid trials
- Relative translational attractiveness of `TACSTD2` as an ADC target

Less reliable:

- Exact affinity or binding claims from the screen
- Any numeric efficacy interpretation from the similarity scores
- ADMET favorability as a substitute for experimental PK/Tox
- Small-molecule relevance for extracellular targets such as `TACSTD2`

## Practical Takeaways

1. `TACSTD2/TROP2` should be treated as the flagship immediate repurposing axis.
2. `CYP1B1` remains the cleanest classic small-molecule target in this set.
3. `PLEKHA6`, `B3GNT3`, and `PTPRE` are better framed as next-wave modality programs than immediate drug repurposing wins.
4. `LDLR` and `GABRB2` produce strong repurposing signals numerically, but their thyroid-oncology specificity needs stricter contextual review before prioritizing them experimentally.

## Files

- Results root: `/opt/thyroid-dash/project/results/v13_drug_discovery`
- Paper section: `/opt/thyroid-dash/project/reports/v13/v13_paper_section.md`
- TROP2 draft: `/opt/thyroid-dash/project/reports/v13/v13_trop2_standalone_draft.md`
- Executive brief: `/opt/thyroid-dash/project/reports/v13/v13_executive_brief.md`
- Dashboard: `/opt/thyroid-dash/project/reports/html/pages/v13_drug_discovery.html`
