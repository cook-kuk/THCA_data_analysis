# v17 TERT promoter — extended integration

_Run: 2026-04-27_

- Recovered TERT-promoter mutated TCGA-THCA patients: **36**
- Sample master rows: 513
- Sample master rows now flagged TERT-mutated: 36
- External (non-TCGA) TERT promoter records (MSK-IMPACT, MSK-CHORD, MSK PDTC/ATC): **1844**

## Cross-tabs
### quad_group
```
                   mutated  wildtype
A_braf_only             25       256
B_ras_only               6        48
D_triple_negative        5       173
```

### aggressive_flag
```
    mutated  wildtype
no       36       477
```

### molecular_subtype
```
           mutated  wildtype
BRAF_like       29       363
RAS_like         7       104
unknown          0        10
```

### driver_anchor
```
         mutated  wildtype
BRAF          25       256
NTRK           1         0
RAS            6        48
TP53           0         1
unknown        4       172
```

### four_group_v2
```
            mutated  wildtype
BRAF_only         0       250
RAS_only          0        48
TERT+            36         0
Triple_neg        0       170
```

## Fisher exact tests
```json
{}
```

## Survival
```json
{
  "TERT_logrank": {
    "n_mutated": 36,
    "n_wildtype": 468,
    "events_mutated": 6,
    "events_wildtype": 10,
    "p_value": 4.918261229904833e-06,
    "test_statistic": 20.86885779620563
  },
  "four_group_logrank": {
    "p_value": 3.780984315139726e-05,
    "test_statistic": 23.137167321537195,
    "groups": {
      "BRAF_only": 250,
      "Triple_neg": 170,
      "RAS_only": 48,
      "TERT+": 36
    }
  }
}
```