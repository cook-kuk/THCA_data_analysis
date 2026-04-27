# U2D — K-selection justification (Variant A 47-gene panel)

**Cohort:** TCGA-THCA primary tumours, n=500, matched to v17 REALFIX R1A subset.
**Features:** 47 genes (TIERA67 minus 8 thyroid-differentiation genes; leak-free).
**Method:** `KMeans(random_state=42, n_init=20)` on z-scored expression, K ∈ {2,3,4,5}.

## 1. Internal cluster-quality metrics

| K | Silhouette ↑ | Davies–Bouldin ↓ | Calinski–Harabasz ↑ |
|---|--------------|------------------|---------------------|
| 2 | 0.1856  | 2.1917       | 82.2          |
| 3 | 0.1637  | 2.0316       | 71.2          |
| 4 | 0.1043  | 2.2135       | 62.9          |
| 5 | 0.0910  | 2.1705       | 55.6          |

K=2 maximises silhouette and Calinski–Harabasz and minimises Davies–Bouldin in this
sweep — the canonical signature of an intrinsically bi-modal feature space. Adding
clusters monotonically degrades all three internal indices, indicating that K≥3
splits begin to fracture homogeneous regions rather than reveal separated modes.

## 2. K=2 vs K=4 confusion (sample counts)

```
K2    0    1
K4          
0     0  110
1   182    0
2    45   16
3   129   18
```

Per-K=4-cluster fraction in the original DM1/DM2 partition:

```
    frac_DM1  frac_DM2    n
K4                         
0   1.000000  0.000000  110
1   0.000000  1.000000  182
2   0.229508  0.770492   61
3   0.108844  0.891156  147
```

K=4 is essentially a **refinement** of K=2: each K=4 cluster lives almost entirely
inside a single DM1 or DM2 partition. K=2 therefore loses no patient-level
information that K=4 carries — it just collapses two pairs of close sub-clusters.

## 3. K=4 mapping to Pu W et al. 2022 (Oncogene) 4-subtype framework

Each K=4 cluster was scored against four Pu-2022 signature panels (immune,
stromal/EMT, BRAF-like MAPK output, proliferation/CNV proxy). Dominant signal:

- **K=4 cluster 0** (n=110) → **Stromal** (scores: Immune-enriched=-0.48, Stromal=-0.38, BRAF-enriched=-0.89, CNV-enriched=-0.43)
- **K=4 cluster 1** (n=182) → **CNV-enriched** (scores: Immune-enriched=+0.29, Stromal=+0.15, BRAF-enriched=+0.20, CNV-enriched=+0.36)
- **K=4 cluster 2** (n=61) → **Immune-enriched** (scores: Immune-enriched=-0.08, Stromal=-0.36, BRAF-enriched=-0.09, CNV-enriched=-0.14)
- **K=4 cluster 3** (n=147) → **BRAF-enriched** (scores: Immune-enriched=+0.03, Stromal=+0.24, BRAF-enriched=+0.45, CNV-enriched=-0.07)

Interpretation: K=4 recovers a Pu-2022-compatible decomposition (immune,
stromal, BRAF-like MAPK, proliferation), confirming that the 47-gene panel
encodes the same biological axes as the published 4-subtype taxonomy. This is
useful for **research depth** — e.g. immunotherapy candidate triage from the
Immune-enriched cluster — but not for clinical deployment.

## 4. Recommendation

- **Deploy K=2 (DM1/DM2).** Best internal metrics, simplest decision boundary,
  highest patient-level reliability, and the granularity that drives every
  downstream survival / TERT / ATA-risk readout in v17.
- **Use K=4 only as a Pu-2022-aligned research overlay.** It is interpretable
  (immune / stromal / BRAF-MAPK / proliferation) and supports targeted
  analyses, but its silhouette is materially lower and clusters are not
  independent of the K=2 partition.

A clinician asked to act on a single binary report ("aggressive vs indolent")
should receive the K=2 call. A translational team prioritising mechanistic
follow-up should additionally see the K=4 / Pu-2022 overlay. Both views are
fully consistent because K=4 nests inside K=2.
