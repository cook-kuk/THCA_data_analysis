#!/usr/bin/env bash
set -euo pipefail
mkdir -p baker_rosetta_interface_out
# InterfaceAnalyzer chain grouping may need manual adjustment for multi-chain TCR-pMHC complexes.
InterfaceAnalyzer.default.linuxgccrelease -s /home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v2_sota_sprint_2026_05_09/tcr_extension/md_escalation/p0_structures/pilot_complexes/6VRN_HMTEVVRHC_HLA-A0201_chainP.pdb -out:file:scorefile baker_rosetta_interface_out/BAKERFILTER_0000.sc || true
# InterfaceAnalyzer chain grouping may need manual adjustment for multi-chain TCR-pMHC complexes.
InterfaceAnalyzer.default.linuxgccrelease -s /home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v2_sota_sprint_2026_05_09/tcr_extension/md_escalation/p0_structures/pilot_complexes/6VRN_HMTEVVRHC_HLA-A0201_chainP.pdb -out:file:scorefile baker_rosetta_interface_out/BAKERFILTER_0001.sc || true
# InterfaceAnalyzer chain grouping may need manual adjustment for multi-chain TCR-pMHC complexes.
InterfaceAnalyzer.default.linuxgccrelease -s /home/seungho/personal/THCA_data_analysis/project/results/cross_neo_md_audit_2026_05_10/counterfactual_design/positive_control_pdbs/extracted/5NQK_ELAGIGILTV_HLA-A0201_chainP.pdb -out:file:scorefile baker_rosetta_interface_out/BAKERFILTER_0002.sc || true
# InterfaceAnalyzer chain grouping may need manual adjustment for multi-chain TCR-pMHC complexes.
InterfaceAnalyzer.default.linuxgccrelease -s /home/seungho/personal/THCA_data_analysis/project/results/cross_neo_md_audit_2026_05_10/counterfactual_design/positive_control_pdbs/extracted/5NQK_ELAGIGILTV_HLA-A0201_chainP.pdb -out:file:scorefile baker_rosetta_interface_out/BAKERFILTER_0003.sc || true
# InterfaceAnalyzer chain grouping may need manual adjustment for multi-chain TCR-pMHC complexes.
InterfaceAnalyzer.default.linuxgccrelease -s /home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v2_sota_sprint_2026_05_09/tcr_extension/md_escalation/p0_structures/pilot_complexes/6UON_GADGVGKSAL_HLA-C0802_chainF.pdb -out:file:scorefile baker_rosetta_interface_out/BAKERFILTER_0004.sc || true
# InterfaceAnalyzer chain grouping may need manual adjustment for multi-chain TCR-pMHC complexes.
InterfaceAnalyzer.default.linuxgccrelease -s /home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v2_sota_sprint_2026_05_09/tcr_extension/md_escalation/p0_structures/pilot_complexes/6UON_GADGVGKSAL_HLA-C0802_chainF.pdb -out:file:scorefile baker_rosetta_interface_out/BAKERFILTER_0005.sc || true
# InterfaceAnalyzer chain grouping may need manual adjustment for multi-chain TCR-pMHC complexes.
InterfaceAnalyzer.default.linuxgccrelease -s /home/seungho/personal/THCA_data_analysis/project/results/cross_neo_md_audit_2026_05_10/counterfactual_design/positive_control_pdbs/extracted/6ULR_GADGVGKSA_HLA-C0802_chainC.pdb -out:file:scorefile baker_rosetta_interface_out/BAKERFILTER_0006.sc || true
