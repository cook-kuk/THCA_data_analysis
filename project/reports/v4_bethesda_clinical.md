# v4 Track B - Bethesda clinical decision simulation

- 6 prevalences x 20,000 synthetic patients
- Costs: lobectomy 3.5M KRW, FU 0.25M KRW, panel 0.3M KRW, missed 30M KRW **(Korean single-payer estimates, NOT generalizable)**

## Operating points
- prev=5% | Se>=0.95 (rule-out): Se=0.952 Sp=0.809 PPV=0.208 NPV=0.997 thr=0.43
- prev=5% | balanced (max Se+Sp): Se=0.894 Sp=0.897 PPV=0.314 NPV=0.994 thr=0.51
- prev=5% | Sp>=0.95 (rule-in): Se=0.812 Sp=0.954 PPV=0.482 NPV=0.990 thr=0.59
- prev=10% | Se>=0.95 (rule-out): Se=0.956 Sp=0.794 PPV=0.340 NPV=0.994 thr=0.42
- prev=10% | balanced (max Se+Sp): Se=0.893 Sp=0.890 PPV=0.475 NPV=0.987 thr=0.50
- prev=10% | Sp>=0.95 (rule-in): Se=0.785 Sp=0.953 PPV=0.648 NPV=0.976 thr=0.59
- prev=15% | Se>=0.95 (rule-out): Se=0.952 Sp=0.780 PPV=0.433 NPV=0.989 thr=0.41
- prev=15% | balanced (max Se+Sp): Se=0.899 Sp=0.888 PPV=0.586 NPV=0.980 thr=0.50
- prev=15% | Sp>=0.95 (rule-in): Se=0.790 Sp=0.952 PPV=0.745 NPV=0.963 thr=0.59
- prev=20% | Se>=0.95 (rule-out): Se=0.952 Sp=0.782 PPV=0.522 NPV=0.985 thr=0.41
- prev=20% | balanced (max Se+Sp): Se=0.875 Sp=0.908 PPV=0.704 NPV=0.967 thr=0.52
- prev=20% | Sp>=0.95 (rule-in): Se=0.803 Sp=0.951 PPV=0.803 NPV=0.951 thr=0.58
- prev=30% | Se>=0.95 (rule-out): Se=0.952 Sp=0.782 PPV=0.652 NPV=0.974 thr=0.41
- prev=30% | balanced (max Se+Sp): Se=0.895 Sp=0.880 PPV=0.761 NPV=0.951 thr=0.49
- prev=30% | Sp>=0.95 (rule-in): Se=0.780 Sp=0.953 PPV=0.876 NPV=0.910 thr=0.59
- prev=40% | Se>=0.95 (rule-out): Se=0.953 Sp=0.780 PPV=0.743 NPV=0.962 thr=0.41
- prev=40% | balanced (max Se+Sp): Se=0.895 Sp=0.878 PPV=0.831 NPV=0.926 thr=0.49
- prev=40% | Sp>=0.95 (rule-in): Se=0.770 Sp=0.951 PPV=0.913 NPV=0.861 thr=0.59

## Caveats
- Synthetic cohort; not prospective.
- TCGA is surgical tissue, not FNA.
- Cost estimates single-payer-specific.
