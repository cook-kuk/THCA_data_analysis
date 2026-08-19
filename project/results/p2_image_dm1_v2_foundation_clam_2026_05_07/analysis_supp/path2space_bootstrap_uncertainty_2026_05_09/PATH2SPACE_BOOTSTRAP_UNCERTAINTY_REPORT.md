# Path2Space-inspired bootstrap uncertainty

## Verdict

- Bootstrap replicates per evidence layer: 2000.
- GSE250521 smoothed spot-level DM1/RAI rho: **0.440**; 95% bootstrap CI **[0.408, 0.475]**.
- GSE250521 coordinate-domain DM1/RAI rho: **0.581**; 95% bootstrap CI **[0.381, 0.729]**.
- GSE230424 smoothed spot-level DM1/low-RAI rho: **0.644**; 95% bootstrap CI **[0.634, 0.655]**.
- GSE230424 coord+QC residual-target rho: **0.232**; 95% bootstrap CI **[0.217, 0.248]**.

## Interpretation

Bootstrap intervals remain above zero for the primary spot-level and domain-level DM1/RAI evidence layers. These intervals are descriptive clustered/stratified uncertainty checks, not a replacement for independent-cohort validation.

## Summary Table

| evidence                                           | level   |     n |   n_groups |   observed_centered_rho |   boot_mean |   boot_ci_low |   boot_ci_high |   boot_p05 |   boot_p95 |   n_boot |
|:---------------------------------------------------|:--------|------:|-----------:|------------------------:|------------:|--------------:|---------------:|-----------:|-----------:|---------:|
| GSE250521 UNI DM1/RAI smoothed spots               | spot    |  3200 |         16 |                  0.4404 |      0.4408 |        0.4078 |         0.4747 |     0.4141 |     0.4682 |     2000 |
| GSE250521 UNI DM1/RAI coordinate domains           | domain  |    96 |         16 |                  0.5811 |      0.5649 |        0.3811 |         0.7292 |     0.409  |     0.708  |     2000 |
| GSE230424 H&E DM1/low-RAI smoothed spots           | spot    | 15489 |          4 |                  0.6444 |      0.6443 |        0.6342 |         0.6546 |     0.636  |     0.653  |     2000 |
| GSE230424 H&E DM1/low-RAI coord+QC residual target | spot    | 15489 |          4 |                  0.2323 |      0.2326 |        0.2175 |         0.2477 |     0.22   |     0.2455 |     2000 |
| GSE230424 H&E DM1/low-RAI coordinate domains       | domain  |    39 |          4 |                  0.9099 |      0.8949 |        0.7884 |         0.9562 |     0.8115 |     0.9509 |     2000 |
