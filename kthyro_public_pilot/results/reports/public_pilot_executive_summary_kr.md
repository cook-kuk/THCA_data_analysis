# K-Thyro Public Pilot Executive Summary

## 한 줄 결론

공개데이터 pilot 결과, 갑상선암은 평균 유전체 변이만으로 치료반응을 설명하기 어렵고, RAI 분화상태, HLA/APM 면역가시성, CD8 exclusion, myeloid/CAF barrier, drug-delivery failure가 서로 다른 치료취약성 niche로 분리될 가능성이 있다.

## 삼성육성과제용 메시지

본 과제는 갑상선암을 평균적인 driver mutation 질환이 아니라, 수술 전후 perturbation과 공간오믹스·병리·기능실험으로 읽는 치료취약성 생태계로 재정의한다.

## Pilot 근거

TCGA-THCA 환자 수준에서는 RAI 분화, HLA/APM, cytotoxic T cell, myeloid/CAF, hypoxia/vascular proxy, proliferation 축을 분리해 정량화했다. GSE250521 공간전사체에서는 spot 수준 치료취약성 score와 niche label을 만들고, slide 수준 요약과 공간 coherence를 계산했다. DepMap/PRISM 기반 후보 표는 redifferentiation, HLA/APM rescue, epigenetic/IFN 축 perturbation을 기능실험 후보로 제안한다.

## Claim boundaries

- Public pilot supports hypothesis generation, not clinical deployment.
- Spatial transcriptomics supports spatial expression states, not direct peptide presentation.
- Drug-delivery failure score is proxy unless validated by fluorescent drug/nanoparticle imaging.
- RAI-restorable score is hypothesis unless iodide uptake assay validates it.
- HLA/APM-low state requires protein validation by mIHC/IF or GeoMx.

## Experimental validation plan

FFPE cohort, mIHC/IF panel, GeoMx/ROI validation, fresh tissue organoid/slice culture, iodide uptake assay, HLA/APM rescue assay, PBMC/TIL co-culture, fluorescent liposome/nanoparticle distribution imaging을 연결한다.

## Go/no-go

GO. 공개데이터는 제안서의 pilot evidence와 figure를 만들기에 충분하지만, drug-delivery와 RAI-restorable claim은 반드시 fresh tissue 기능실험으로 검증해야 한다.
