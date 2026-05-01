---
title: "DEPRECATED (partial) — d4p1_panasian_meta"
date: 2026-05-01
status: partially superseded
supersededBy: project/results/p2_pillar1_forest/
note: "이 디렉토리의 일부 파일만 deprecated. 정확한 분류는 § 0 표 참조."
---

# ⚠️ PARTIAL DEPRECATION — d4p1_panasian_meta

**이 디렉토리의 일부 파일만 deprecated.** Per-sample arcasHLA genotype 파일들 (`korean_PTC_pool_n908.tsv`, `GSE286332_arcasHLA_genotypes.tsv` 등) 은 canonical 분석 (`v17_paper2_pillar1_forest.py`) 과 manuscript supplementary (S9) 에서 계속 사용 중 — 그대로 유지.

## 0. 파일별 상태 분류

| 파일 | 상태 | 이유 |
|---|---|---|
| `D4P1_summary.json` | ⚠️ DEPRECATED | `Chen2018_*` field name + 잘못된 OR/freq 숫자 |
| `panasian_forest_meta.tsv` | ⚠️ DEPRECATED | 같은 이유 (column name 에 `Chen2018_*`) |
| `chen2018_han_chinese_GD_published.tsv` | ⚠️ DEPRECATED | 잘못된 published values (실제 Chu 2018 값과 다름) |
| `korean_PTC_pool_n908.tsv` | ✅ KEEP | Per-sample 4-digit alleles, `v17_paper2_pillar1_forest.py:53` 에서 read, manuscript S9 |
| `korean_PTC_pool_focus_freq.tsv` | ✅ KEEP | Korean PTC pool focus allele freq (Chen 의존 없음) |
| `GSE286332_arcasHLA_genotypes.tsv` | ✅ KEEP | GSE286332 per-sample genotypes (raw arcasHLA output) |
| `GSE286332_focus_freq.tsv` | ✅ KEEP | GSE286332 focus allele freq |
| `GSE286332_PTC_vs_PTCHT_fisher.tsv` | ✅ KEEP | PTC vs PTC+HT Fisher test |

## 1. 왜 일부 deprecated 인가

위 표의 ⚠️ 표시 파일들은 다음 두 가지 stale issue 가 있음:

### 1. Citation 오류
`Chen2018_*` field name + `Chen et al. 2018 PMC6161647 / Front Endocrinol 9:467` 인용은 **잘못된 cite**.

- **잘못 표기**: "Chen 2018", "Front Endocrinol 9:467"
- **실제 paper**: Chu X et al. 2018 *J Med Genet* 55(10):685–692, doi:10.1136/jmedgenet-2017-105146 (PMC 6161647 은 동일하나 author/journal 잘못 attribute)

### 2. 숫자 오류
DPB1*05:01 의 reported Chen 2018 OR / freq 가 실제 Chu 2018 published values 와 다름:

| Field | D4P1 (잘못) | Chu 2018 actual |
|---|---|---|
| GD freq | 0.61 | **0.44** |
| ctrl freq | 0.39 | **0.313** |
| OR (GD vs ctrl) | 2.45 | **1.90** |
| p | 1e−30 | **1.7e−26** |

이 차이 때문에 D4P1 의 forest meta 결과는 wrong baseline 으로 계산됨.

## 정정된 분석 위치

✅ **Use instead**: `project/results/p2_pillar1_forest/`

| 파일 | 용도 |
|---|---|
| `PILLAR1_FOREST_SUMMARY.md` | 2026-05-03 PARTIAL → STRONG decision summary |
| `P2_PILLAR1_summary.json` | Full structured data (Chu 2018 정확 cite) |
| `chu2018_allele_summary.tsv` | Chu 2018 published 8 alleles values |
| `forest_meta_results.tsv` | 6 alleles × 3 arms (Chu ctrl / Korean PTC / Chu GD) |
| `random_effects_pooled.tsv` | DerSimonian-Laird pooled OR + Cochran Q + I² |
| `korean_subcohort_heterogeneity.tsv` | K2 / Lee / GSE286332 sub-cohort heterogeneity |
| `sensitivity_4scenarios.tsv` | 4-scenario robustness check |
| `methods_paragraph.md` | Cell Press paste-ready Methods |
| `discussion_paragraph.md` | Cell Press paste-ready Discussion |

## 코드 chain

| 단계 | 파일 | 상태 |
|---|---|---|
| OLD analysis | `project/notebooks_or_scripts/v17_D4P1_forest_meta.py` | ⚠️ deprecated header 추가됨 |
| NEW analysis | `project/notebooks_or_scripts/v17_paper2_pillar1_forest.py` | ✅ canonical |
| OLD output | `project/results/d4p1_panasian_meta/` | ⚠️ 이 directory (deprecated) |
| NEW output | `project/results/p2_pillar1_forest/` | ✅ canonical |

## Audit trail

- 2026-05-03: Citation correction discovered via PMC 6161647 web fetch — actual author/journal/numbers verified.
- 2026-05-03: Pillar 1 STRONG sprint regenerated with correct Chu 2018 values (`v17_paper2_pillar1_forest.py`).
- 2026-05-01: 본 DEPRECATED note 추가, p2_advisor_discussion.html 통합 작업 중 발견.

## 이 디렉토리의 파일을 reference 하는 다른 코드

```bash
$ grep -rln "d4p1_panasian_meta" project/ --include="*.py" --include="*.md" --include="*.html"
```

→ 이 디렉토리를 reference 하는 모든 후속 분석 / figure / manuscript 는 `p2_pillar1_forest/` 로 update 필요.

---

*DEPRECATION notice 작성 2026-05-01. p2_pillar1_forest/ canonical reference 사용 권장.*
