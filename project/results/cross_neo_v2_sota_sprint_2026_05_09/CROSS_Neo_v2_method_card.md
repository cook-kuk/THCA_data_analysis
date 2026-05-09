# CROSS-Neo v2 Method Card

Inputs: mutant peptide, WT peptide where available, HLA, source/study metadata, source window, structure/missingness flags, split-safe derived features.

Excluded: public predictor scores from MHCflurry, NetMHCpan, NetMHCstabpan, BigMHC, PRIME, MixMHCpred.

Model families: clean counterfactual anchors, frozen hashed PLM pilot, multimodal feature models, source-balanced proxies, rank normalization, prespecified MoE, and imported v1 QK fallback diagnostics.

Selection: outer test labels are not used for fitting, thresholding, calibration, fusion, or model promotion.
