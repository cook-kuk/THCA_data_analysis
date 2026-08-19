# Paper 2 submission overview pack

This is scaffolding only: figure/table/data-use roadmap and caption skeletons. It does not write protected manuscript sections.

## Representative main figure flow
- Fig 1: Whole-study overview and clinical translation path -- H&E -> foundation model -> DM1/RAI spatial state -> reflex molecular testing
- Fig 2: Image-DM1 classifier performance -- UNI LOTO AUC 0.852; RAS-like AUC 0.718
- Fig 3: GSE250521 spatial RNA recovery -- slide-centered rho 0.440; PTC rho 0.487; LPTC rho 0.471
- Fig 4: External thyroid validation in GSE230424 -- H&E rho 0.644; residual-target rho 0.232
- Fig 5: Robustness, hotspot, dose-response, and caveat moat -- hotspot lift 3.5x/3.3x; decile rho 0.407/0.607
- Fig 6: Impact and validation unlock map -- Cell Reports Medicine primary; K2/Bundang H&E unlock for Nat Commun attempt

## Data-use inventory
- TCGA-THCA H&E WSI: Primary image-DM1 classifier and attention heatmaps; used in Main Fig 2; Main Fig 5; Table 1; Table S1-S2.
- TCGA-THCA RNA labels: Ground-truth molecular label for H&E model; used in Main Fig 2; Table 1; Table S1.
- GSE250521: Spatial transcriptomic recovery and stage generalization; used in Main Fig 3; Main Fig 4; Table 1; Table S3-S5.
- GSE230424: External thyroid H&E-to-DM1/RAI support plus residual/QC caveat; used in Main Fig 4; Main Fig 5; Table 1; Table S4-S6.
- GSE248205: Negative-control/caveat lock for local H&E immune-axis overclaim; used in Main Fig 5; Table S6; reviewer response table.
- Korean K2 / Bundang H&E: Highest-impact revision unlock and Nat Commun ceiling raiser; used in Main Fig 6 roadmap; Table 4; presubmission plan.

## Files
- `F00_paper2_whole_story_overview.png`
- `F01_paper2_data_use_map.png`
- `F02_paper2_figure_table_roadmap.png`
- `F03_paper2_representative_figure_bundle.png`
