# Top Quantum Permutation Check

Candidate: `quantum_kernel_no_anchor_gamma1.0`

Features: `GP_quantum,VQC,W7A_QK_only,W7A_full,mean_pLDDT_peptide,min_pLDDT_peptide,mean_pLDDT_HLA,anchor_pLDDT,interface_contacts_8A,interface_contacts_10A,n_buried_residues_8A,mean_min_pep_to_hla_CA_dist,max_min_pep_to_hla_CA_dist,radius_of_gyration_peptide,end_to_end_CA_dist,peptide_helicity_proxy,tcr_motif_score`

Observed AUROC: 0.785
Observed AUPRC: 0.603

Permutation n=100

- AUROC null mean: 0.495; 95% range: 0.310-0.696; empirical p=0.010
- AUPRC null mean: 0.276; 95% range: 0.176-0.477; empirical p=0.010

Interpretation: permutation checks whether this small strict-set internal CV signal is obviously label-random. It does not replace an external leakage-controlled benchmark.
