# Literature review — RAI response genomics in thyroid cancer

## 1. Why RAI response matters clinically

Radioactive iodine (I-131) therapy is the cornerstone adjuvant and salvage treatment for differentiated thyroid cancer (DTC), particularly after total thyroidectomy. Successful RAI eradicates residual thyroid tissue and iodine-avid metastases. Failure of RAI — **RAI-refractory thyroid cancer (RAIR-DTC)** — defines a clinical phenotype with markedly worse prognosis: 10-year survival drops from > 90% in RAI-responsive disease to ~10–14% in distant-metastatic RAIR-DTC.

The definition of "RAI-refractory" remains operational and heterogeneous (ATA 2015 criteria, Schlumberger 2014). Subtypes commonly include:

- No uptake on diagnostic or post-therapy scan (no-avidity refractoriness).
- Mixed avidity (some lesions avid, others not).
- Avidity present but no structural response (functional refractoriness).
- Initial response then progression (acquired refractoriness).

Identifying patients at risk for RAIR before exposure to repeated I-131 doses and cumulative radiation burden is a major unmet need.

## 2. Why RAI failure is not fully explained by NIS alone

The historical paradigm focused on **NIS (SLC5A5)** as the rate-limiting step. Lower NIS expression in tumor vs. normal thyroid, partial NIS mis-trafficking to subcellular compartments, and post-translational silencing have all been invoked. But:

- Many BRAF-mutant tumors retain measurable NIS mRNA yet fail to concentrate iodine functionally.
- TG, TPO, DUOX1/2, IYD, PAX8, NKX2-1, FOXE1, TSHR, and DIO1/DIO2 are coordinately downregulated in a substantial fraction of dedifferentiated tumors.
- Re-induction of iodine uptake by MAPK inhibition restores **multiple** genes in this program, not only NIS.

This motivates a **multi-gene differentiation / iodide-handling program** view, captured in the field's "Thyroid Differentiation Score (TDS)" and in our 8-gene panel.

## 3. Thyroid differentiation score and iodide-handling genes

The TDS was introduced in the TCGA-THCA paper (Cell 2014) and refined by Landa et al. (JCI 2016, Cancer 2019). The 16-gene TDS captures the broader iodine-handling program. Our 8-gene panel is a parsimonious subset:

**Panel**: SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1.

Coverage logic:
- **SLC5A5**: iodide uptake at the basolateral membrane.
- **TG**: iodine storage scaffold and apical secretion.
- **TPO**: iodide organification onto TG tyrosines.
- **TSHR**: upstream signal that maintains differentiation.
- **PAX8, NKX2-1, FOXE1**: thyroid lineage TFs that bind TG/TPO/SLC5A5 promoters.
- **DIO1**: thyroid hormone metabolism, marker of mature thyrocyte program.

The score is calibrated to detect **silencing of the program as a whole**, which empirically tracks RAI failure across multiple cohorts.

## 4. MAPK / BRAF / RAS / TERT / TP53 / PI3K associations with RAI refractoriness

Drivers consistently associated with RAIR or aggressive thyroid cancer include:

- **BRAF V600E** — strongest single driver of TDS loss and PTC dedifferentiation; ~60% of classic PTC, near-universal in tall-cell variants.
- **RAS family** (NRAS Q61R/K, HRAS Q61R, KRAS G12) — characterizes follicular variant PTC and follicular thyroid cancer; preserves some differentiation but is over-represented in FTC and ATC.
- **TERT promoter** C228T / C250T — late event; strongly co-occurs with BRAF or RAS in aggressive disease (Liu 2014, Melo 2014). **TERT × BRAF** double-positive has the highest mortality.
- **TP53** — common in PDTC and ATC.
- **PIK3CA / PTEN** — ATC and PDTC; defines a sub-population with poor RAI response.
- **RET / NTRK / ALK / PPARG fusions** — actionable but most associated with younger PTC; RAI response variable.

Each driver maps onto our two-axis frame differently (Paper 1 R17 finding):
- BRAF: silences RAI machinery without HT-overlap.
- RAS: preserves RAI machinery with HT-overlap.
- TERT+: enriches in the silenced + immune-active "dark-matter" zone.

## 5. Evidence for molecular gray zones and non-binary RAI avidity

Mu et al. JCEM 2024 explicitly identified **four uptake patterns** in 220 distant-metastatic DTC patients:

1. Initially RAI-refractory.
2. Continually RAI-avid.
3. Gradually RAI-refractory (avidity decays over RAI cycles).
4. Partly RAI-refractory (avid in some metastases, not others).

This refutes a simple binary view. The 8-gene panel and TDS-like signatures should therefore be benchmarked against **ordinal or multi-class** outcomes, not just RAI-avid vs RAI-refractory.

Boucai 2023 CCR similarly framed "exceptional responder" as one tail of a distribution, not a categorical label.

## 6. Public datasets with RAI labels

Comprehensive list in `docs/dataset_inventory.md` and `config/datasets.yaml`. Top entries:

- **GSE151179 / GSE151181** (PTC, mRNA, RAI-avid vs refractory). Public, downloadable.
- **GSE151180** (PTC, miRNA companion). Public.
- **GSE299988** (small PTC, RAI-avid vs non-avid). Public.
- **Boucai 2023 CCR** (exceptional responders). Supplement + request.
- **Mu 2024 JCEM HRA004166** (4-class uptake patterns, n=220). Controlled access.
- **GSE184362** (Pu 2021 scRNA, includes RAI-refractory distant metastasis sample). Public, heavy.
- **GSE112202** (digoxin redifferentiation). Public.
- **Redifferentiation trials** (selumetinib, dabrafenib/trametinib, MERAIODE). Trial supplement only.

## 7. Limitations of existing datasets

- TCGA-THCA has **no RAI response field**. Only ~12% of TCGA-THCA patients develop recurrence; OS is dominated by non-thyroid causes.
- GSE151179-class datasets are small (n ≤ 50) and unbalanced.
- Mu 2024 is targeted NGS only — no expression, can only validate the **driver-mutation half** of the manuscript.
- Redifferentiation trials have N ~ 10–20.
- No public dataset combines (a) large N, (b) RAI response endpoint, (c) RNA-seq, (d) full driver panel.

This limitation **is the point** — the manuscript argues that integrating an expression-based panel with the small RAI-labeled cohorts adds value beyond any single dataset.

## 8. How an eight-gene panel can be positioned safely

Position the panel as:

> "A parsimonious 8-gene readout of the thyroid differentiation / iodide-handling program, validated across multiple cohorts as **associated with** RAI refractoriness phenotypes, and proposed as a **risk-stratification readout** complementary to driver mutation status."

Avoid:

- Overclaiming as a stand-alone predictor.
- Comparing AUC against larger published signatures without like-for-like CV.
- Implying clinical utility without prospective data.

Affirm:

- Cross-cohort consistency.
- Driver-stratified value (BRAF / RAS / driver-negative).
- Methylation + protein corroboration (Paper 1 4-pillar argument).
- Single-cell substrate (Lu 2023 dark-matter zone, 38.3% of ATC malignant cells).
- Explicit Tier mapping for every claim.

## 9. Reviewer risks and defensive wording

See `docs/reviewer_risk_register.md` for the full register. Key defensive lines:

- "TCGA-THCA is used as discovery; RAI-labeled cohorts (GSE151179, Boucai 2023, Mu 2024) are used for label-anchored validation."
- "The panel is associated with RAI avidity (Tier 2) and is consistent with RAI refractoriness, not validated as a clinical RAI response predictor (which would require Tier 1 prospective data)."
- "RAI refractoriness is modelled as a tiered/gray-zone phenotype, not a binary, leveraging Mu 2024's 4-class uptake pattern data."

## 10. Recommended figure order (manuscript)

1. **Figure 1**: Biological model — RAI as a multi-gate program (uptake → lineage TFs → organification → storage → killing).
2. **Figure 2**: Discovery — 8-gene differentiation-silence axis in TCGA-THCA + integrated bulk cohorts, with driver stratification.
3. **Figure 3**: Validation against RAI-labeled public data (GSE151179 + GSE299988 + Mu 2024 4-class).
4. **Figure 4**: Genomic context — driver overlap with panel silencing; molecular gray zone framing.
5. **Figure 5**: Redifferentiation — panel score increases with MAPK inhibition; restored iodide handling.
6. **Supp**: scRNA cellular substrate (Lu 2023), proteomic corroboration (Mun 2025), methylation (HM450 β), Paper 1 R17 cross-reference.
