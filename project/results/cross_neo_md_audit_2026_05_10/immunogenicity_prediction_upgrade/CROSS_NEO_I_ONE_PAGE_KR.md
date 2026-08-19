# CROSS-Neo-I 면역원성 예측 업그레이드 요약

임팩트를 더 올리려면 지금의 후보선별/MD 패키지를 **면역원성 예측 모델 계약서**로 바꿔야 합니다.

- endpoint는 `mutant-specific T-cell activation > WT/decoy controls` 입니다.
- BigMHC-like presentation, PRIME/IMPROVE-like foreignness, PISTE/EPACT-like TCR expert, OpenMM MD, tumor context, calibration을 분리합니다.
- unknown을 negative로 학습하면 안 됩니다. assay-confirmed negative만 negative입니다.
- HMTEVVRHC와 GADGVGKSAL은 CROSS-Neo-I의 flagship label-generation cases입니다.

한 줄 결론: **이제 대박 포인트는 예측 점수 하나가 아니라, 면역원성 label을 만들고 calibrate할 수 있는 모델/실험 contract를 갖췄다는 점입니다.**
