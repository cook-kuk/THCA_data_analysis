# next_steps

1. Lock the exact BRS71 and TierA67 gene lists before interpreting panel coverage or panel-based classifier performance.
2. Add controlled-access cohorts only after data use approval, then keep them in separate ingestion and provenance layers.
3. For TCGA, add fusion and copy-number annotations if the downstream question requires molecular subtype refinement beyond histology.
4. Re-run microarray from raw CEL files with RMA if cross-study reproducibility is a hard requirement.
5. Evaluate whether GSE213647 tumor-only PTC vs follicular-neoplasm labels are clinically comparable enough to serve as true external validation.
