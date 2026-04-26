# v4 Track A - ComBat batch-correction rescue

- Cohorts pooled: 5
- Excluded: [['GSE213647', 'ENSG-format index; no symbol mapping']]
- Shared genes: 11731
- Samples: 679
- Pre identifiability macro-AUC: 1.000
- Post identifiability macro-AUC: 0.487
- Pre LODO mean AUC: 0.997
- Post LODO mean AUC: 0.011
- **Verdict: UNRECOVERABLE**

Decision rule: RESCUED if post_ident<0.70 AND mean_post_auc>0.85; UNRECOVERABLE if post_ident<0.70 AND post_auc<0.75; INSUFFICIENT otherwise.
