# Dark Matter — Phase 1+2 Web Dashboard

자체 완결 HTML 페이지. CDN (Bootstrap 5 + DataTables + jQuery + **Plotly 2.32**)만 외부 의존.

## 두 가지 페이지

- **v2.html** (권장, 1090 lines, 92KB) — 12 sections · 18 static figures + **4 Plotly interactive figures** · 27 tables · 81 glossary terms
- **index.html** (legacy, 43KB) — 초기 dashboard

## 로컬 서빙

```bash
cd project/results/dark_matter_phase2/web
python3 -m http.server 9999 --bind 127.0.0.1
# 브라우저: http://localhost:9999/v2.html
```

또는 `python3 -m http.server 9999 --bind 0.0.0.0` 으로 LAN에 노출.

## v2.html 인터랙티브 figure (4개)

1. **I-1 sc 8-gene ↔ FVPTC scatter** — Phase 1 GSE241184 thyrocytes 2,500 cell subsample, hover로 sample/score
2. **I-2 Cross-cohort r forest** — 5 sc context r 비교, 95% CI 에러바
3. **I-3 Multi-site trajectory** — P→T→LN per-patient line, hover unified
4. **I-4 Driver class sunburst** — TCGA n=482 hierarchy

## 포함 내용

- **단일 index.html (~30KB)** — Bootstrap UI, sticky TOC, 18개 figure embed, sortable tables
- **figures/** — 18개 PNG (Phase 1 + Phase 2 모든 deliverable)
- **data/** — 6개 TSV/JSON (interactive table용 + machine-readable summaries)

## 주요 figure 인덱스

| 파일 | 출처 | 설명 |
|---|---|---|
| `fig5D_headline_8gene_vs_fvptc.png` | Phase 1 | 단일환자 sc r=0.905 |
| `fig5D_v2_external_validation.png` | Phase 2 P2-A1 | GSE193581 13-sample r=0.893 |
| `p2a2_gse184362_summary.png` | Phase 2 P2-A2 | 6 adult PTC r=0.889 |
| `fig5F_multisite_trajectory.png` | Phase 2 multi-site | P→T→LN dedifferentiation 4 patients |
| `fig8_KM_curves.png` | Phase 2 extras | KM 5-class p=0.035 + DM1/DM2 |
| `figS3_dialu_stability.png` | Phase 2 P2-G | Bootstrap ARI=0.994 cluster stability |
| `fig_S1_driver_class_x_cluster.png` | Phase 2 P2-D | Driver class × cluster stacked bar |

## 약어 / 용어 사전

페이지 내 #glossary 섹션에 50+ 항목 정리.

## References

PMID 직접 링크된 14개 references in #refs section.
