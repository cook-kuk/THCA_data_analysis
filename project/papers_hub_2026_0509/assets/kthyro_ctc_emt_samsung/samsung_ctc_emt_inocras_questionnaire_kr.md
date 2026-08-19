# Inocras Vendor Questionnaire

## Scope Fit

1. PTC thyroidectomy research cohort에서 matched tumor-normal sequencing 지원 가능 여부.
2. 적합한 서비스: CancerVision, MRDVision, custom research WGS, target-enhanced WGS, bioinformatics-only 중 무엇인지.
3. FASTQ, BAM/CRAM, VCF, CNV, SV, fusion, TMB/MSI/HRD/signature output 제공 가능 여부.
4. BRAF, TERT promoter, RET/NTRK/ALK fusions, RAS, EIF1AX, PPM1D, CHEK2, CNV/SV를 research report에 포함 가능한지.

## Sample Requirements

5. FFPE input: curls/slides 수, tumor purity threshold, H&E requirement, macrodissection 지원.
6. matched normal input: blood, buccal, saliva 중 허용 범위.
7. cfDNA input: tube, blood volume, plasma volume, minimum cfDNA mass, processing window.
8. CTC-enriched pellet 또는 low-input/single-cell CTC DNA/RNA 처리 가능 여부.
9. 가능하다면 WGA/library method와 SNV/CNV/SV/fusion 한계.
10. 불가능하다면 research-only exploratory input으로 받을 수 있는 대안.

## Performance / Validation

11. cfDNA/MRDVision LOD와 low tumor fraction에서의 성능.
12. PTC/thyroid cancer validation data 포함 여부.
13. early-stage/low-shedding tumor에서 sensitivity 저하 예상.
14. TERT promoter, fusion, structural variant 처리 방식.
15. FFPE tumor-informed variant로 serial cfDNA interpretation 가능한지.
16. germline, clonal hematopoiesis, tumor-derived variant 구분 방식.

## Operations / Governance

17. sample type별 TAT.
18. tumor-normal WGS, serial cfDNA, CTC-enriched sample, bioinformatics-only 단가.
19. batch requirement와 Year 2 scale-up 가능성.
20. failed sample policy와 rerun cost.
21. 데이터 저장 위치, 국외 이전 여부, encryption, access control, IRB 문구.
22. CAP/CLIA status와 research-use output boundary.
23. 한국 병원 데이터가 국내에 머물 수 있는지.
24. IRB 전 assay design 논의할 scientific contact.

## Deliverables

25. de-identified report 예시.
26. raw data manifest 예시.
27. genome build와 annotation database version.
28. pipeline versioning/reproducibility statement.
29. molecular tumor board/expert interpretation 가능 여부.
30. publication policy와 acknowledgment requirement.
