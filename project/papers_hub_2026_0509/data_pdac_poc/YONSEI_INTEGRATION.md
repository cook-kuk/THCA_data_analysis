# Yonsei (or any institutional) PDAC cohort · drop-in integration template

This document is the contract between Lumenix and the Yonsei team.
When the Yonsei data arrives in the format below, the entire downstream
analysis runs in **one command**, in **&lt; 60 seconds**, and produces
manuscript-ready text + an updated meta-forest figure.

---

## 1 · Required input format

A single TSV or CSV file. **Five required columns** (case-insensitive, exact names):

| Column | Type | Description |
|---|---|---|
| `sample_id` | str | patient identifier (any string) |
| `os_months` | float | overall survival in months (0 &lt; OS &lt; 200) |
| `os_event` | int | 1 = death/event, 0 = alive/censored |
| `kras_allele` | str | one of: `G12D`, `G12V`, `G12R`, `G12C`, `G13D`, `Q61H`, `Q61R`, `KRAS_other`, `WT` |
| `age` | float | age at diagnosis (years) |

**Optional columns** (enrich the analysis):

| Column | Type | What it unlocks |
|---|---|---|
| `moffitt_call` | str | basal-like / classical → cross-cohort paradox replication |
| `ancestry` | str | Korean / Asian / European → HLA-bias context |
| `hla_a`, `hla_b`, `hla_c`, `hla_drb1`, `hla_dpb1`, `hla_dqb1` | str | per-patient joint mut × HLA actionability |
| `rna_path` or `rna_csv` | str | path to mRNA z-score table (one row per patient) → TME paradox per-patient score |

## 2 · Example input row

```tsv
sample_id    os_months    os_event    kras_allele    age    moffitt_call    ancestry
YS_PDAC_001    14.2    1    G12D    62    basal-like    Korean
YS_PDAC_002    23.7    1    G12V    71    classical    Korean
YS_PDAC_003    8.4     0    G12R    58    basal-like    Korean
...
```

## 3 · Run command

```bash
python /data/pdac_poc/scripts/21_yonsei_dropin.py /path/to/yonsei.tsv \
    --label "Yonsei PDAC 2026"
```

**This produces in one run:**
1. `/data/pdac_poc/results/yonsei/SUMMARY.json` — Cox + meta JSON
2. `/data/pdac_poc/results/yonsei/yonsei_updated_forest.png` — meta-forest with Yonsei row added
3. `/data/pdac_poc/results/yonsei/YONSEI_RESULTS.md` — manuscript-ready paragraph

## 4 · What the output answers

- **Q1: Does Yonsei replicate the G12D HR > 1 finding?**
  → Yonsei per-cohort Cox table with HR, 95% CI, p for each KRAS allele
- **Q2: Does Yonsei strengthen the meta-analysis?**
  → updated pooled HR (95% CI, p) across 6 cohorts (TCGA + MSK Cancer Cell + MSK Nat Med + paad_tcga + PRINCE + Yonsei)
- **Q3: Where does Yonsei fall relative to the other cohorts?**
  → updated forest plot — Yonsei row appears in violet
- **Q4: How big does Yonsei need to be to flip the meta?**
  → drop-one sensitivity replays automatically

## 5 · Privacy + IRB

- **No raw RNA matrix is uploaded outside Yonsei.** Optionally, Yonsei pre-computes z-scores locally and shares the per-sample summary; no individual-patient genomic file leaves Yonsei.
- Anonymized `sample_id` is required only as a key for joining; no patient-identifying information needed.
- The script runs entirely on Yonsei or Lumenix infrastructure depending on agreement; data sharing limited to the joint summary statistics for the manuscript.

## 6 · Expected timeline

| Step | Lumenix | Yonsei |
|---|---|---|
| Format check | — | 1 day · sample 5-row TSV |
| First run | 1 hour | — |
| Manuscript paragraph drafted | 1 day | — |
| Joint methods section | 3 days | review |
| Manuscript-ready submission package | 1 week | sign-off |

## 7 · NComm strengthening this enables

Adding Yonsei (n ≥ 100) to the meta-analysis would:
1. **Push k from 5 to 6 cohorts**, total n with OS from 3,113 → 3,213+
2. **Add a Korean-population validation** — directly relevant to our Korean off-the-shelf cassette claim (off-shelf coverage 14.93%)
3. **Test the paradox in a Korean cohort** — does the inflamed-AND-suppressed quadrant replicate?
4. **Optionally**: with HLA columns, validate the joint mut × HLA actionability calculation that currently uses AFND priors

## 8 · If the data is not yet available

The script tolerates incomplete columns gracefully — for example, an early Yonsei dump with only `sample_id`, `os_months`, `os_event`, `kras_allele` (no age) still produces a Kaplan-Meier per-allele plot + 3-way logrank, with a clear "Cox not run, awaiting age covariate" flag.

---

**Contact**: Seungho Cook, Lumenix · `kukshomr@gmail.com` · 2026-05-07.
