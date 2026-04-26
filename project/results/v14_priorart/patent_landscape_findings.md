# Patent landscape — TROP2 ADC × thyroid × BRAF (DONE 2026-04-25)

_Author: Seungho Cook · Source: Google Patents via WebSearch + WebFetch on JP7525633B2, US20250152731A1._
_Replaces `patent_scan_status.md`._

## TL;DR

**No granted or pending patent claims a BRAF-stratified TROP2 ADC use in thyroid cancer.** The TROP2-ADC IP space is dense with composition-of-matter and broad method-of-use claims, but BRAF / BRAF V600E stratification specifically for thyroid TROP2 ADC therapy is **unclaimed in publicly searchable filings as of 2026-04-25**. Our v14 stratification proposal does not infringe nor is it scooped by IP.

## Key relevant patents

### Composition-of-matter / antibody backbone (Daiichi Sankyo)

| Publication | Subject | Comment |
|---|---|---|
| `EP3088419` (and `PT3088419T`, `ES2703903T3`, `RU2743077C2`, `CA2933666C`) | Anti-TROP2 ADC with camptothecin-class payload | Datopotamab deruxtecan backbone IP. Broad. |
| `US9850312B2`, `US9770517B2`, `US20140377287` | Anti-TROP2 ADC + SN-38 / paclitaxel / cytotoxic drug | Older, broad. Cover sacituzumab govitecan space. |
| `WO2018086239A1` | Anti-TACSTD2 ADC (China) | Linker-optimised; broad indication coverage. |
| `EP2594589A1` | Anti-Trop-2 antibody (Japanese groups) | Antibody-only, no ADC. |

### Method-of-use / combination (closer to our space)

| Publication | Title | Assignee | Priority | Thyroid? | BRAF? | Verdict |
|---|---|---|---|---|---|---|
| **US20250152731A1** / **WO2023060283A2** | Anti-TROP2 antibody combination therapies and methods of use | **UT System** | **2021-10-08** | **YES** (in 16-cancer list) | **No** | **Closest neighbour.** Claims TROP2-ADC + DNMT-i (decitabine) or Zeb1-i for **TROP2-low** tumours (opposite end from our BRAF-high → TROP2-high gate). Different stratification angle; **does not block our hypothesis**. Worth one Discussion citation. |
| **JP7525633B2** | Biomarkers for sacituzumab govitecan therapy | **Immunomedics** | **2020-03-20**; granted **2024-07-30** | **No** (thyroid not enumerated) | **No** (BRAF not in biomarker set) | Highest-name-similarity hit. Actually claims **DNA-damage-repair (DDR) gene biomarkers**: BRCA1/2, CHEK2, MSH2/6, TP53, CDKN1A, etc. Indications: breast, urothelial, lung, GI, GU, gyn, head/neck. **No thyroid, no BRAF.** **Does not block.** |
| `US20170224837A1` / `US10954305B2` / `CA3011372A1` | ABCG2 inhibitor + sacituzumab govitecan for SN-38 resistance in TROP2 cancers | Immunomedics | 2016 | Indirect | No | Resistance-rescue angle, not stratification. **No conflict.** |
| `US10413539B2` / `US20180110772A1` / `US20190381032A1` | Sacituzumab govitecan in metastatic urothelial cancer | Immunomedics | 2017 | No | No | Urothelial-specific. **No conflict.** |

### BRAF detection (general, not thyroid TROP2 specific)

| Publication | Title | Comment |
|---|---|---|
| `WO2013062976A1` | Methods of detecting BRAF mutations in cancer | AS-PCR assay for BRAF V600E / Y600E; analytical only, no therapy claim. **No conflict.** |

### NSCLC TROP2 patent (orthogonal)

| Publication | Comment |
|---|---|
| `WO2022011197A1` | TROP2 detection in lung cancer — diagnostic, not thyroid. |

## Where the v14 proposal sits in the IP map

```
   Composition-of-matter (anti-TROP2 ADC)         ← held by Daiichi Sankyo, Immunomedics/Gilead
       │
       ├── Method-of-use, TROP2-low + DNMT-i      ← UT System (US20250152731A1, includes thyroid, BRAF-agnostic)
       ├── Method-of-use, urothelial              ← Immunomedics
       ├── Combination, ABCG2-i resistance        ← Immunomedics
       └── Biomarker, DDR-gene panel              ← Immunomedics (JP7525633B2, no thyroid, no BRAF)

   ╔═══════════════════════════════════════════╗
   ║  EMPTY SPACE: BRAF-stratified TROP2 ADC    ║
   ║  use in thyroid cancer                     ║   ← v14 proposal sits here
   ║  (no claims found 2026-04-25)              ║
   ╚═══════════════════════════════════════════╝
```

## Implication for the manuscript

1. **No IP barrier** to publishing the v14 BRAF-like × TACSTD2-high two-gate stratification hypothesis.
2. **Discussion citation worth adding:** `US20250152731A1` (UT System, 2021) covers thyroid as an indication for TROP2 ADC + DNMT-i combination — a complementary stratification angle (TROP2-low rescue) distinct from ours (BRAF-high → TROP2-high enrichment). One sentence acknowledging this in §6.X positions our work as orthogonal.
3. **Honest disclosure:** the broad TROP2-ADC composition patents (Daiichi Sankyo, Immunomedics) cover the agents in the live trials; this is normal commercial-development context and not a novelty issue for our paper.
4. **No need to file IP** on the v14 proposal — academic publication is sufficient and intended.

## Caveats

- WebSearch sees Google Patents indexed listings; very recent (last 30 days) WIPO/PCT filings may have ~30-90 day lag.
- Korean (KIPO) and Chinese (CNIPA) regional filings not separately swept beyond what Google Patents surfaces.
- Defensive publications (IP.com, defensive disclosures) not searched.

## Time spent
~10 min WebSearch + ~5 min WebFetch on top two candidate patents. This closes the §E.1 patent action item.
