# NeoImmune-Stack / CLEAN-Neo++ repo audit

## Decision
Build an integration layer, not a new foundation model. Existing local algorithms are reusable and should be separated from frozen public predictors.

## Existing local algorithms to integrate
- Structure_LR: found or referenced
- Wave8_TCR_SelfSim_full: registry-required; direct artifact not yet found
- ESM2_Bayesian: registry-required; direct artifact not yet found
- GP_quantum: registry-required; direct artifact not yet found
- VQC: registry-required; direct artifact not yet found
- W7A_QK_only: registry-required; direct artifact not yet found
- W7A_full: registry-required; direct artifact not yet found
- W7B_stacked: registry-required; direct artifact not yet found
- quantum_kernel_no_anchor_gamma1.0: registry-required; direct artifact not yet found
- quantum_kernel_compact_gamma0.5: registry-required; direct artifact not yet found
- Stack_mean: registry-required; direct artifact not yet found
- Stack_median: registry-required; direct artifact not yet found
- Stack_LR: registry-required; direct artifact not yet found
- BAR-Neo: registry-required; direct artifact not yet found
- BAR-Neo-BMA: registry-required; direct artifact not yet found
- BAR-Neo-X: registry-required; direct artifact not yet found
- stress_guarded_clean_stack: registry-required; direct artifact not yet found

## Existing public/external comparator artifacts
- NetMHCpan_4.1_EL: existing output or raw source found
- NetMHCpan_4.1_BA: existing output or raw source found
- MHCflurry_2.0_presentation: existing output or raw source found
- MHCflurry_2.0_affinity: existing output or raw source found
- BigMHC_EL: existing output or raw source found
- BigMHC_IM: existing output or raw source found
- PRIME: existing output or raw source found
- PRIME2.1: pending adapter/status stub
- HLApollo: pending adapter/status stub
- pVACtools_parser: pending adapter/status stub
- NeoDisc_parser: pending adapter/status stub
- Topiary_pending: pending adapter/status stub
- MixMHCpred_pending: existing output or raw source found
- NetMHCIIpan_pending: pending adapter/status stub

## High-value reusable artifacts
- `project/results/cross_neo_immunogenicity_algorithm_compare_2026_05_10/bigmhc_el_predictions.csv`
- `project/results/cross_neo_immunogenicity_algorithm_compare_2026_05_10/bigmhc_im_predictions.csv`
- `project/results/cross_neo_immunogenicity_algorithm_compare_2026_05_10/bigmhc_input_candidates.csv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_bigmhc_work/dbpepneo2_mhci_im.csv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_bigmhc_work/dbpepneo2_mhci_input.csv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_bigmhc_work/mcpas_im.csv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_bigmhc_work/mcpas_input.csv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_bigmhc_work/neodb_im.csv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_bigmhc_work/neodb_input.csv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_bigmhc_work/nepdb_im.csv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_bigmhc_work/nepdb_input.csv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/scripts/05_score_bigmhc.py`
- `project/results/p_neo_bayesian_2026_05_09/wave11/wave11_predictions/BigMHC_IM__dbpepneo2_mhci.tsv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/wave11_predictions/BigMHC_IM__itsndb.tsv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/wave11_predictions/BigMHC_IM__mcpas.tsv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/wave11_predictions/BigMHC_IM__neodb.tsv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/wave11_predictions/BigMHC_IM__nepdb.tsv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_bigmhc_work/cedar_partial_im.csv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_bigmhc_work/cedar_partial_input.csv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/wave11_predictions/BigMHC_IM__cedar_partial.tsv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_bigmhc_work/tesla_mmc4_im.csv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_bigmhc_work/tesla_mmc4_input.csv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_bigmhc_work/tesla_mmc7_im.csv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_bigmhc_work/tesla_mmc7_input.csv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/wave11_predictions/BigMHC_IM__tesla_mmc4.tsv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/wave11_predictions/BigMHC_IM__tesla_mmc7.tsv`
- `project/results/cross_neo_product_demo_2026_05_10/selected_trainsets/02_CEDAR_include_positive_prior_downweighted.tsv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/test_bundles/cedar_partial.tsv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/wave11_predictions/RF_biophys__cedar_partial.tsv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/wave11_predictions/TransPHLA__cedar_partial.tsv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/test_bundles/improve_cedar.tsv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/scripts/03_score_mhcflurry.py`
- `project/results/p_neo_bayesian_2026_05_09/wave11/wave11_predictions/MHCflurry__dbpepneo2_mhci.tsv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/wave11_predictions/MHCflurry__itsndb.tsv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/wave11_predictions/MHCflurry__mcpas.tsv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/wave11_predictions/MHCflurry__neodb.tsv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/wave11_predictions/MHCflurry__nepdb.tsv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/wave11_predictions/MHCflurry__cedar_partial.tsv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/wave11_predictions/MHCflurry__tesla_mmc4.tsv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/wave11_predictions/MHCflurry__tesla_mmc7.tsv`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A0101.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A0201.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A0202.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A0203.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A0205.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A0206.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A0211.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A0224.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A0301.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A1101.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A2301.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A2402.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A2403.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A2404.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A2501.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A2601.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A2902.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A3001.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A3002.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A3101.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A3201.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A3303.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A6801.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A6802.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_A6901.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_B0702.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_B0801.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_B1402.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_B1501.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_B1507.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_B1517.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_B1801.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_B2705.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_B3501.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_B3701.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_B3801.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_B4001.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_B4002.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_B4102.txt`
- `project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work/dbpepneo2_mhci/in_B4402.txt`

## Missing dependencies / execution risk
- netMHCpan: not on PATH
- netmhcpan: not on PATH
- mhcflurry-predict: available
- python: available
- Rscript: available
- pvacseq: not on PATH
- pvacbind: not on PATH

## High-risk assumptions
- Public predictor scores may overlap with public training data; they are barred from the clean science track.
- Presentation labels and immunogenicity labels must remain separated.
- Patient-level top-N metrics are only meaningful when patient IDs and positive labels exist.
- Existing strong random-split signals are not accepted as mechanistic or deployable evidence without source/patient/peptide leakage controls.
- LLM usage is restricted to narrative audit, rationale summarization, and report drafting; it is not a vaccine efficacy predictor.

## Integration map
- Canonical candidates: `clean_neobench_master.tsv`, `cross_neo_v0/master_table.tsv`, `/data/neoantigen_vaccine_hub/data_processed/neoantigen_master.csv`.
- Local scores: `clean_neobench_method_scores.tsv`, BAR-Neo score files, cross-neo score files, wave11 local model predictions.
- External scores: wave11 predictions for BigMHC/MHCflurry/PRIME/NetMHCpan plus adapter stubs for unavailable tools.
- Leakage evidence: clean_neobench overlap flags, cross_neo retrieval audit, neoantigen hub leakage JSON.
