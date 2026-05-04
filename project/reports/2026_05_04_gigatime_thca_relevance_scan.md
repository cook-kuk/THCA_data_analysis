# GigaTIME (Valanarasu 2026 Cell) — THCA relevance scan

**Date:** 2026-05-04 · **Owner:** Seungho Cook
**Reference:** Valanarasu JMJ et al. *Multimodal AI generates virtual population for tumor microenvironment modeling.* **Cell** 189(2):386–400.e19, Jan 2026. **DOI:** 10.1016/j.cell.2025.11.016. **PMID:** 41371214. **License:** CC-BY (open access).
**Purpose:** Bound what we can claim about GigaTIME's THCA relevance from publicly accessible metadata. Source for citations in Paper 2A Methods supplement, Paper 1 Reviewer Q scaffold, and future-work design doc.

---

## 1. Confirmed metadata

| field | value | source |
|---|---|---|
| Citation | Cell 189(2):386–400.e19 | OpenAlex, CrossRef |
| First author | Jeya Maria Jose Valanarasu (Microsoft Research) | OpenAlex |
| Senior author | Hoifung Poon (Microsoft Research) | OpenAlex |
| Total authors | 21 (Microsoft Research + Providence Cancer Institute + Providence Research Network) | OpenAlex |
| Open access | CC-BY | OpenAlex |
| Citation impact (FWCI) | 10.47 (top 10%) | OpenAlex |
| Cohort | 14,256 patients across 51 hospitals + ~1,000 clinics in 7 US states (Providence Health) | abstract |
| Training data scale | 40 × 10⁶ paired cells (H&E ↔ mIF) | abstract |
| Protein panel | 21 proteins (composition NOT enumerated in public metadata) | abstract |
| Virtual mIF generated | 299,376 datasets | abstract |
| Cancer types | 24 (composition NOT enumerated in public metadata) | abstract |
| Cancer subtypes | 306 | abstract |
| Independent validation | 10,200 TCGA patients | abstract |
| Statistical associations | 1,234 protein-biomarker-staging-survival | abstract |
| Code/model release URLs | NOT located in publicly accessible search (GitHub microsoft/GigaTIME → 404; HuggingFace search → none) | search 2026-05-04 |
| Supplementary file URLs | NOT located in publicly accessible search; full PDF gated behind cell.com Cloudflare challenge | scan 2026-05-04 |

## 2. THCA inclusion — UNKNOWN (publicly)

The abstract states "24 cancer types and 306 subtypes" but does NOT enumerate the cancer type list. Search across OpenAlex metadata, EuropePMC abstract, public Microsoft Research project pages, and Google search with "GigaTIME thyroid" returned **no statement confirming or denying thyroid cancer inclusion**.

**Implications**:
- We cannot claim a direct comparison "GigaTIME does X for THCA, we do Y for THCA"
- We CAN claim "GigaTIME provides a general-cancer multimodal AI baseline; THCA-specific transcriptional axis prediction is not addressed in their reported associations"
- For Methods supplement and Reviewer Q: frame the comparison at **task class** (image-to-protein-translation vs continuous transcriptional regression) rather than cancer overlap

## 3. Protein panel — UNKNOWN (publicly)

21 proteins not enumerated. From Microsoft Research's prior publications and the Providence cohort context, **likely candidates** (based on standard ImmunoSEQ / Akoya / Lunaphore TIME panels) include:
- **Immune compartment**: CD3, CD4, CD8, FoxP3 (T cells); CD20, CD138 (B/plasma); CD68, CD163 (macrophages); CD56 (NK); HLA-DR; PD-1, PD-L1, CTLA-4 (checkpoint)
- **Tumor structural**: pan-CK, vimentin, αSMA, Ki67
- **TLS markers** (possibly): CD20 + CD3 patches

**Likely NOT in panel** (transcriptional axis / thyroid lineage markers): FOXE1, NKX2-1, PAX8, SLC5A5/NIS, TG, TPO, DIO1, TSHR, STAT3, FOSL1, DNMT1, AICDA — these are RNA / TF activity / methylation / B-cell-clonality readouts that are not standard mIF-imaged proteins.

**Implication**: Even if our 8-gene RAI / DM1 axis genes were in their 21-protein panel (highly unlikely), the readout would be at protein-translation-from-H&E level, NOT at depth-residualized transcriptional regression level — different task class as stated above.

## 4. Cited-text claims that are confirmed (safe to use)

✅ "GigaTIME translates H&E into virtual multiplex-IF images at population scale (40M paired cells, 14,256 patients, 24 cancer types, 306 subtypes)" — abstract verbatim
✅ "Independent validation on 10,200 TCGA patients corroborated findings" — abstract verbatim
✅ "1,234 statistically significant associations linking proteins, biomarkers, staging, and survival" — abstract verbatim
✅ "Microsoft Research + Providence Health collaboration" — author affiliations
✅ Cell 189(2), Jan 2026, CC-BY — bibliographic record

## 5. Claims that are NOT yet supported (do NOT make without further verification)

❌ "GigaTIME includes thyroid cancer" — not confirmed
❌ "GigaTIME panel includes TROP2 / FOXE1 / NKX2-1 / etc." — not confirmed
❌ "GigaTIME outperforms ResNet50 ImageNet on transcriptional regression" — not confirmed (different task)
❌ "GigaTIME model weights are publicly available" — not located
❌ Specific AUC / accuracy benchmarks vs ResNet50 — not located in available metadata

## 6. How this scan is used by other strategies

- **S1 Methods supp comparison table** — uses §1 confirmed numbers only; explicitly notes §2/§3 unknowns as "scope distinction" rather than gap
- **S2 Reviewer Q scaffold** — pre-emptively addresses "Why didn't you compare directly to GigaTIME on THCA?" with §2 honesty
- **S4 Future-work design doc** — uses §5 unknowns to define what THCA-paired training would need (RNA + H&E + 14k+ patients to match GigaTIME power)
- **S5 Cover letter scaffold** — uses §1 + §4 to position our axis vs their translator (orthogonal complementarity)
- **S6 Venue justification memo** — uses §1 publication context (Cell 189(2) Jan 2026) for venue ladder argument

## 7. Future re-verification (not blocking)

Per `2026_05_04_image_dm1_final_nogo_decision.md` §5 condition #4 (foundation model access acquired), if GigaTIME model code/weights are released later (Microsoft Research often releases ~6-12 months after publication):
- check `https://github.com/microsoft` for `GigaTIME` repo
- check `https://huggingface.co/microsoft` for GigaTIME-* models
- if/when available, re-evaluate §3 panel inclusion + §2 THCA inclusion

This memo can be revised at that time. For current marathon (5/4–6/13), the public metadata is sufficient for citation context.

---

*Scan based on publicly accessible metadata as of 2026-05-04. Cell paper full PDF behind Cloudflare; Europe PMC ingestion pending; no GitHub/HuggingFace public release found. All citations and comparisons in downstream strategies stay within §4 confirmed claims; §5 unsupported claims are explicitly excluded.*
