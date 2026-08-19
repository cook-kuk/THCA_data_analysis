# Counterfactual Embedding Report

No PLM was fine-tuned.
WT peptide is missing in the current master table, so WT vectors are zero and `wt_missing_flag=1`.
The representation is a deterministic AA-property + k-mer fallback with HLA pseudo-sequence encoding.
Matrix shape: (2715, 353)
