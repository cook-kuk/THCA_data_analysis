# K-Thyro Public Pilot PPT Evidence Package

이 패키지는 삼성미래기술육성사업 PPT 제작용 public-data pilot 근거 묶음입니다.

핵심 framing:

> TCGA 505명 + GSE250521 spatial 16 slide에서 치료취약성 niche가 분리되고, 그 niche가 공간적으로 응집된다.

주의:

> Public-data pilot supports hypothesis generation and validation design, not clinical deployment.

## Package Contents

- `docs/KTHYRO_SAMSUNG_PPT_MASTER_CONTEXT.md`  
  전체 상황, 핵심 숫자, 10-slide storyline, claim boundary, figure path.

- `docs/ANALYSIS_FLOW_COMPREHENSIVE.md`  
  TCGA, spatial, scRNA, external validation, drug resource, mutation-only insufficiency 분석흐름.

- `docs/FIGURE_TABLE_GUIDE.md`  
  PPT에 넣을 figure/table별 설명, caption, allowed/forbidden wording.

- `docs/SLIDE_COPY_KR.md`  
  10장 PPT에 바로 붙일 한국어 제목, 메시지, bullets.

- `docs/CLAIM_BOUNDARY_BOXES.md`  
  TCGA/spatial/drug/clinical validation의 allowed claim / forbidden claim.

- `figures/`  
  PPT용 핵심 PNG figure 사본.

- `tables/`  
  PPT 숫자 검산용 핵심 TSV table 사본.

- `reports/`  
  기존 report 및 executive summary 사본.

## Go/No-Go

최종 판단: **GO**

근거:

- TCGA-THCA 505명에서 therapeutic vulnerability label 분리.
- GSE250521 spatial 16개 slide, 57,144 spots에서 16/16 slide non-random niche coherence.
- External public validation matrix v2: 301 tests, strong 144, moderate 30.
- Exact K-Thyro external GPL570 validation: 206 samples, 168 tests, strong 65, moderate 12.
- BRAF group 내부에서도 7개 vulnerability label로 분산되어 mutation-only stratification 한계 확인.

