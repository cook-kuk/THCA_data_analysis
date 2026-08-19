# CROSS-Neo endpoint energy gate

## Verdict

This layer defines how to score endpoint interaction energies from existing trajectories. It is a cheap ranking gate, not a proof of immunogenicity.

## What to compute

- mutant ensemble mean interaction energy
- WT ensemble mean interaction energy
- decoy ensemble mean interaction energy
- mutant-minus-WT bootstrap confidence interval
- mutant-minus-decoy bootstrap confidence interval

## Rows to prioritize

|   endpoint_order | peptide    | hla_4digit   | physics_tier        | md_label   |   md_score | launch_sequence                                                                                          | endpoint_step                                                        | compare_against                                     | output_fields                                                     | decision_use                                                                     | hold_rule                                                                  |
|-----------------:|:-----------|:-------------|:--------------------|:-----------|-----------:|:---------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------|:----------------------------------------------------|:------------------------------------------------------------------|:---------------------------------------------------------------------------------|:---------------------------------------------------------------------------|
|                1 | HMTEVVRHC  | HLA-A*02:01  | P0_FLAGSHIP_PHYSICS | NA         |          0 | endpoint energy now -> PMF if endpoint supports mutant specificity -> FEP only after WT mapping is clean | compare mutant vs WT vs decoy endpoint energies over ensemble frames | mutant_minus_wt and mutant_minus_decoy bootstrap CI | mean_energy, delta_energy, bootstrap_ci, frame_count, control_gap | If mutant is consistently more favorable than WT/decoy, keep physics gate alive. | No PMF/FEP unless endpoint energies are favorable and WT mapping is clean. |
|                2 | GADGVGKSAL | HLA-C*08:02  | P0_FLAGSHIP_PHYSICS | NA         |          0 | endpoint energy now -> PMF only if TCR-pMHC model stays stable                                           | compare mutant vs WT vs decoy endpoint energies over ensemble frames | mutant_minus_wt and mutant_minus_decoy bootstrap CI | mean_energy, delta_energy, bootstrap_ci, frame_count, control_gap | If mutant is consistently more favorable than WT/decoy, keep physics gate alive. | No PMF/FEP unless endpoint energies are favorable and WT mapping is clean. |
|                3 | ILDKVLVHL  | HLA-A*02:01  | P3_NO_PHYSICS_NOW   | NA         |          0 | cheap geometry/endpoint check only                                                                       | geometry sanity only                                                 | no expensive endpoint energy yet                    | mean_energy, delta_energy, bootstrap_ci, frame_count, control_gap | Hold until assay or MD uncertainty justifies more physics budget.                | No PMF/FEP unless endpoint energies are favorable and WT mapping is clean. |
|                4 | SYLDSGIHF  | HLA-A*24:02  | P3_NO_PHYSICS_NOW   | NA         |          0 | cheap geometry/endpoint check only                                                                       | geometry sanity only                                                 | no expensive endpoint energy yet                    | mean_energy, delta_energy, bootstrap_ci, frame_count, control_gap | Hold until assay or MD uncertainty justifies more physics budget.                | No PMF/FEP unless endpoint energies are favorable and WT mapping is clean. |
|                5 | KLILWRGLK  | HLA-A*03:01  | P3_NO_PHYSICS_NOW   | NA         |          0 | cheap geometry/endpoint check only                                                                       | geometry sanity only                                                 | no expensive endpoint energy yet                    | mean_energy, delta_energy, bootstrap_ci, frame_count, control_gap | Hold until assay or MD uncertainty justifies more physics budget.                | No PMF/FEP unless endpoint energies are favorable and WT mapping is clean. |
|                6 | YVDFREYEYY | HLA-A*01:01  | P3_NO_PHYSICS_NOW   | NA         |          0 | cheap geometry/endpoint check only                                                                       | geometry sanity only                                                 | no expensive endpoint energy yet                    | mean_energy, delta_energy, bootstrap_ci, frame_count, control_gap | Hold until assay or MD uncertainty justifies more physics budget.                | No PMF/FEP unless endpoint energies are favorable and WT mapping is clean. |
|                7 | VVGAVGVGK  | HLA-A*11:01  | P3_NO_PHYSICS_NOW   | NA         |          0 | cheap geometry/endpoint check only                                                                       | geometry sanity only                                                 | no expensive endpoint energy yet                    | mean_energy, delta_energy, bootstrap_ci, frame_count, control_gap | Hold until assay or MD uncertainty justifies more physics budget.                | No PMF/FEP unless endpoint energies are favorable and WT mapping is clean. |
|                8 | KLYGLDWAEL | HLA-A*02:01  | P3_NO_PHYSICS_NOW   | NA         |          0 | cheap geometry/endpoint check only                                                                       | geometry sanity only                                                 | no expensive endpoint energy yet                    | mean_energy, delta_energy, bootstrap_ci, frame_count, control_gap | Hold until assay or MD uncertainty justifies more physics budget.                | No PMF/FEP unless endpoint energies are favorable and WT mapping is clean. |

## Claim boundary

Endpoint energy can support a specificity ranking. It cannot prove activation, killing, or clinical utility.

