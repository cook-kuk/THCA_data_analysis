# Korean NGS HLA Control Bundle

## Executive verdict

Korean NGS healthy/control HLA sources do exist. The best immediately usable open source is Baek et al. 2021 PLOS ONE, which typed 173 healthy South Koreans across 11 HLA loci by NGS and provides supporting DOCX tables for allele frequencies. Jung et al. 2023 HLA is larger (n=339 unrelated healthy Koreans, 11 loci NGS), but in this pass only PubMed abstract values were accessible without a full table.

This bundle is useful for Paper 2 source inspection and Korean NGS baseline comparison. It does not solve the carrier-frequency problem by itself, because the accessible Baek/Jung/Choe values are allele frequencies over 2n, not individual carrier counts.

## Files downloaded

- `project/external_refs/korean_ngs_hla_controls/baek2021_plos/article.html`
- `project/external_refs/korean_ngs_hla_controls/baek2021_plos/S9_HLA_A_B_C_freq.docx`
- `project/external_refs/korean_ngs_hla_controls/baek2021_plos/S10_HLA_classII_freq.docx`
- `project/external_refs/korean_ngs_hla_controls/baek2021_plos/S1_haplotype_A_B_DRB1.docx`
- `project/external_refs/korean_ngs_hla_controls/baek2021_plos/S5_haplotype_DRB1_DQB1_DPB1.docx`
- parsed plain text copies under `project/external_refs/korean_ngs_hla_controls/baek2021_plos/parsed/`

## Extracted target allele values

| Source | n | Metric | A*02:07 | B*46:01 | C*01:02 | DRB1*07:01 | DQB1*02:01 | DPB1*05:01 |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| Baek 2021 PLOS ONE NGS | 173 | allele frequency over 2n=346 | 18/346, 5.20% | 21/346, 6.07% | 72/346, 20.81% | 27/346, 7.80% | 8/346, 2.31% | 118/346, 34.10% |
| Jung 2023 HLA NGS abstract | 339 | allele frequency over 2n=678 | not in abstract | not in abstract | 18.44% | not in abstract | not in abstract | 35.10% |
| Baek 2023 HLA amplicon NGS | not visible in abstract | expected allele frequency | full table needed | full table needed | full table needed | full table needed | full table needed | full table needed |
| Choe 2021 Ann Lab Med NGS | 128 | allele frequency over 2n=256, A/B/C/DRB1 only | full table needed | full table needed | 19.9% in abstract | full table needed | not covered | not covered |

## Additional Korean NGS sources found

| Source | What it adds | Limitation |
|---|---|---|
| Baek 2023 HLA, doi `10.1111/tan.14981` | Amplicon-based NGS, 11 loci, healthy donors from South Korea; technical match to the 11-locus NGS baseline idea | PubMed abstract confirms design but not exact target allele table; full journal table or institutional access needed |
| Jung 2023 HLA, doi `10.1111/tan.14980` | Larger n=339 unrelated healthy Koreans, 11-locus NGS; abstract confirms C*01:02 and DPB1*05:01 frequencies | Full table not retrieved locally in this pass |
| Choe 2021 Ann Lab Med, doi `10.3343/alm.2021.41.3.310` | n=128 healthy unrelated Korean adults, 8-digit NGS for A/B/C/DRB1; C*01:02=19.9% in abstract | No DQB1/DPB1; full table needed for A*02:07/B*46:01/DRB1*07:01 |

## Interpretation for Paper 2

- Yes, Korean NGS controls exist and can be used as source evidence.
- Baek 2021 is open-access and locally bundled for direct inspection.
- Jung 2023 is a strong larger NGS cross-check source, but full table access still needs journal/supplement retrieval.
- Baek 2023 is likely the next best 11-locus NGS control, but exact target values were not open from PubMed alone.
- Choe 2021 is useful as an A/B/C/DRB1 high-resolution NGS cross-check only.
- These are allele-frequency controls, not individual-level carrier-frequency controls.
- Therefore, they can replace or support In/AFND as a cleaner modern NGS baseline, but they should not be mixed with PTC carrier frequency without the same caveat.

## Best next use

1. Use Baek 2021 as the open Korean NGS 11-locus healthy baseline source table.
2. Use Jung 2023 abstract values as a DPB1*05:01/C*01:02 cross-check until the full table is retrieved.
3. Add Baek 2023 and Choe 2021 as source targets for manual table retrieval.
4. Keep Kim 2014 separate as the only currently parsed Korean individual-level carrier control-like reference.
5. If the user wants to manually inspect source tables, start with:
   - `project/external_refs/korean_ngs_hla_controls/baek2021_plos/parsed/S9.txt`
   - `project/external_refs/korean_ngs_hla_controls/baek2021_plos/parsed/S10.txt`

## Sources

- Baek IC et al. PLOS ONE 2021. doi: `10.1371/journal.pone.0253619`.
- Jung K et al. HLA 2023. doi: `10.1111/tan.14980`; PubMed PMID `36719349`.
- Baek IC et al. HLA 2023. doi: `10.1111/tan.14981`; PubMed PMID `36720674`.
- Choe W et al. Ann Lab Med 2021. doi: `10.3343/alm.2021.41.3.310`.
