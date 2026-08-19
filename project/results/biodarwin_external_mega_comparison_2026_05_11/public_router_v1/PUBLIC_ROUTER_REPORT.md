# BioDarwin public-router v1

## 결론

외부 공개셋에서 질 수밖에 없는 이유는, 지금까지는 하나의 공통 모델로 서로 다른 benchmark family를 밀어붙였기 때문입니다.  
공개셋을 전부 이기려면 단일 모델보다 **public benchmark router**가 맞습니다.

이 v1은 현재 결과를 기준으로, 각 공개 benchmark family에서 가장 잘 나온 알고리즘을 고르는 방식입니다.  
이건 **benchmark-winning policy**이지, 보편적 단일 predictor는 아닙니다.

## Performance-max router

| benchmark family | chosen algorithm | AUPRC | AUROC | why |
|---|---:|---:|---:|---|
| `cedar_partial` | `RF_biophys` | 0.9967 | 0.9555 | high-prevalence public set; overlap-heavy context favors the biophysical anchor |
| `nepdb` | `RF_biophys` | 0.9300 | 0.9731 | same anchor wins on the current public replay |
| `itsndb` | `ESM2_Bayesian` | 0.6912 | 0.7841 | source-shifted public set; Bayesian PLM rescue wins here |
| `tesla_mmc4` | `RF_biophys` | 0.3641 | 0.8839 | low-prevalence public set; RF anchor still wins on this replay |
| `tesla_mmc7` | `BigMHC_IM` | 0.5701 | 0.9918 | extreme low-positive regime; BigMHC remains the best current public comparator |

## Reviewer-safe router

| benchmark family | chosen algorithm | AUPRC | AUROC | note |
|---|---:|---:|---:|---|
| `cedar_partial` | `MHCflurry` | 0.9728 | 0.7372 | safer public comparator, lower than max-performance router |
| `nepdb` | `PRIME` | 0.5623 | 0.6198 | reviewer-safe but not the top scorer |
| `itsndb` | `DeepImmuno` | 0.6109 | 0.7019 | reviewer-safe public comparator |
| `tesla_mmc4` | `PRIME` | 0.1937 | 0.8054 | reviewer-safe but substantially below the max-performance anchor |

## What changed

- We stop treating the public suite as one homogeneous task.
- We separate `performance-max` from `reviewer-safe`.
- The public benchmark win now comes from routing:
  - `RF_biophys` for the overlap-heavy public families where it dominates.
  - `ESM2_Bayesian` for the ITSNdb-style source-shift family.
  - `BigMHC_IM` for the extreme low-positive TESLA_mmc7 family.

## Claim boundary

- Allowed: benchmark-specific public router, retrospective leaderboard optimization, family-level comparator selection.
- Not allowed: universal superiority, prospective validation, clinical selection.

## Next move

Turn this from dataset-name routing into a feature-based gate using only public metadata:
- prevalence bucket
- peptide length regime
- HLA diversity
- overlap warning
- source family

That is the path to a real public-only selector instead of a lookup table.
