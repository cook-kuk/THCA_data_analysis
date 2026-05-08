"""Inject per-page customized data-definitions block.

Each page only shows the cohorts / scores / groups / paper-numbering rows
that are actually relevant to its content. Falls back to the full dictionary
for hub-level pages (index, master_view, etc.).

Idempotent: removes any existing DATA_DEFINITIONS_BLOCK first, then re-inserts.
"""
from __future__ import annotations
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# Master rows (id, *cells)
# ---------------------------------------------------------------------------
COHORTS = {
    "TCGA":   ("TCGA-THCA",       "TCGA / GDC",                 "RNA-seq + clinical + mut + meth", "≈496",     "Paper 1 main cohort; 8-gene + DM1/DM2 derivation"),
    "HM450":  ("TCGA-THCA HM450", "TCGA Illumina HM450",        "methylation array",               "503",      "Paper 1 promoter hypermethylation arm"),
    "CBIO":   ("cBioPortal `thca_tcga_pub`", "cBioPortal",      "mutation + clinical",             "TERT⁺ n=36 recovered", "Paper 1 TERT⁺ Cox HR=4.33 (p=4.9e-6)"),
    "MSK":    ("MSK-IMPACT thyroid", "Landa 2016 Cell, MSKCC", "targeted DNA + clinical",          "117 (84 PDTC + 33 ATC)", "Paper 1 advanced-disease anchor; Cox HR 2.67"),
    "K2":     ("K2 / KOREAN_K2",  "PRJEB11591 (Yoo 2016, SNU-GMI)", "Korean PTC RNA-seq",          "≈260",     "Korean replication; ≠ Bundang SNUH"),
    "LEE24":  ("Lee 2024 / GSE213647", "GEO, Korean PTC",       "RNA-seq",                         "632",      "Paper 1 Korean external main; 22–28% Hashimoto-like ≈ TCGA 18–20%"),
    "286332": ("GSE286332",       "GEO Korean PTC ± Hashimoto", "RNA-seq",                         "18",       "Paper 2 PTC+HT signature; 8-gene d=−1.6, HLA-II d=+3.65"),
    "76039":  ("GSE76039",        "Landa 2016 JCI",             "PDTC + ATC microarray",           "≈100 (37 PDTC+ATC RNA-seq subset)", "Discussion 3.1 reverse-causality citation; advanced lineage benchmark"),
    "33630":  ("GSE33630",        "GEO, GPL570 Affy U133+2",    "microarray",                      "105 (11 ATC / 49 PTC / 45 N)",       "External RAI/DM1 main replication panel"),
    "65144":  ("GSE65144",        "GEO, GPL570",                "microarray",                      "25 (12 ATC / 13 N)",                 "External ATC vs normal supporting"),
    "29265":  ("GSE29265",        "GEO, GPL570",                "microarray",                      "49 (9 ATC / 20 PTC / 20 paired-N)",  "External paired-normal supp"),
    "53157":  ("GSE53157",        "GEO, GPL570",                "microarray",                      "26 (5 PDTC / 7 PTC / 8 FVPTC / 4 FTC / 2 N)", "PDTC sensitivity; underpowered"),
    "60542":  ("GSE60542",        "GEO, GPL570",                "PTC primary + nodal mets array",  "92 (17 PTC N+ / 11 PTC N0 / 17 LN-met / 24 paired-N / 4 LN)", "Paper 1 supp lymph-node confounding"),
    "126698": ("GSE126698",       "GEO RNA-seq",                "SRA only",                        "28",       "DROP — Series Matrix metadata-only; alignment forbidden"),
    "250521": ("GSE250521",       "GEO Visium ST",              "spatial transcriptomics + H&E",   "PT→ATC trajectory", "Spatial validation; Lumenix demo cohort"),
    "146003": ("GSE146003",       "GEO methylation",            "EPIC array",                      "—",        "Nature Cancer feasibility; delta_beta probe panel"),
    "151179": ("GSE151179",       "GEO post-RAI",               "RNA-seq",                         "52",       "Paper 3 RAI-dediff axis; thyroid_diff d=−1.01 (p=1e-4)"),
    "184362": ("GSE184362",       "Pu 2021 Nat Commun, Fudan",  "scRNA",                           "6 PTC patients / 158K cells", "Paper 1 thyrocyte-intrinsic; per-patient r=0.798–0.886"),
    "193581": ("GSE193581",       "Lu 2023",                    "scRNA + bulk cell line",          "23 samples / 32 GSM", "Paper 1 single-cell external"),
    "241184": ("GSE241184",       "Phase 1 single-cell",        "scRNA",                           "1 patient", "Paper 1 author-independence cross-cohort"),
    "78220":  ("GSE78220",        "Hugo 2016 melanoma ICI",     "bulk RNA-seq + clinical",         "28",       "Paper 3 Track B-lite anchor (component of Hugo+Riaz pooled)"),
    "91061":  ("GSE91061",        "Riaz 2017 nivolumab",        "bulk RNA-seq + WES + clinical",   "109",      "Paper 3 Track B-lite anchor (component of Hugo+Riaz pooled)"),
    "HUGO":   ("Hugo+Riaz pooled","8-cohort melanoma ICI",      "RNA-seq + response",              "n=415 pooled", "Paper 3 Track B-lite; DIAL-lite 4/7 module FLIP"),
    "POZ":    ("Pozdeyev 2018 CCR", "publication supp",         "mutation + clinical",             "779 advanced DTC+ATC", "Disc 3.1 보조 (supp table only; raw → corresponding author)"),
    "BUND":   ("Bundang SNUH",    "Korean retrospective FFPE",  "outreach-stage",                  "no data yet","Future Korean cohort; ≠ K2"),
    "AFND":   ("AFND South Korea","Allele Freq Net DB",         "HLA pop allele freq",             "pooled",   "Paper 2 / Paper 4 Korean baseline (replaces Chu 2018)"),
}

SCORES = {
    "8GENE":   ("8-gene mini-index", "RandomForest-ranked from curated 55-gene pool with drivers excluded by design (NOT unsupervised genome-wide)", "8 genes (Paper 1)",                                         "DM1-like / dedifferentiated"),
    "RAI8":    ("RAI_8",              "RAI-uptake / lineage panel mean of 8 thyroid-differentiation genes",                                            "SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1",          "RAI-avid / differentiated"),
    "DM1S":    ("DM1_like_score",     "−RAI_8 (definition; not an independent measurement)",                                                          "−mean(RAI_8 z)",                                            "DM1-like / dedifferentiated"),
    "NONO":    ("THYROID_NONOVERLAP", "Orthogonal lineage panel sharing NO genes with RAI_8 (panel-overlap artifact control)",                        "SLC26A4, IYD, DUOX1, DUOX2, TFF3, HHEX, GLIS3, DIO2",       "differentiated"),
    "TDS":     ("TDS / TDS_TF / TF_collapse", "Thyroid Differentiation Score TF subset",                                                              "FOXE1, NKX2-1, PAX8, HHEX",                                 "differentiated"),
    "MECH":    ("STAT3_AP1_DNMT",     "Mechanism-arm panel",                                                                                          "STAT3, FOSL1, JUNB, DNMT1, DNMT3B",                         "silencing-active / dedifferentiated"),
    "TIERA":   ("TIERA67",            "67-gene driver-orthogonal index, ARI=0.90 ≈ pan-genome top-5000 ARI=0.92",                                     "67 genes",                                                  "DM1-like"),
    "TDIFF":   ("thyroid_diff",       "Module score for thyroid differentiation (Paper 3)",                                                            "module",                                                    "differentiated"),
    "MYEL":    ("myeloid",            "Paper 3 module; ATC d=+2.53 (4/4 sign across 4 thyroid cohorts)",                                              "module",                                                    "myeloid-infiltrated"),
    "HLA":     ("HLA-II / IFN-γ",     "Paper 2 PTC+HT signature",                                                                                     "module / FDR=2e-4",                                         "antigen-presenting / inflamed"),
    "COH":     ("Cohen's d",          "standardized mean difference",                                                                                 "(μ1−μ2)/sp",                                                "positive ⇒ higher in group1"),
    "DIAL":    ("DIAL",               "Domain Invariance via Adversarial Learning leak metric (LODO ComBat); 0=true biology, 1=leak",                  "5 THCA classifiers under LODO ComBat = 0.000",              "leak (avoid)"),
}

GROUPS = {
    "DM1":     ("DM1",        "Molecular Dark Matter cluster 1 — driver-orthogonal dedifferentiated subgroup (main Paper 1 finding)", "TCGA-THCA derived"),
    "DM1B":    ("DM1 sub-B",  "DM1 subset n=56, 96% mutation-negative; TCGA equivalent of K2 NBNR + 4× Hashimoto-like rate",          "TCGA → K2 bridge"),
    "DM2":     ("DM2",        "Hashimoto-like cluster, OR up to 5× (p=6e-10); resolved 18/18 paradox via GSE286332 → TCGA transfer",   "TCGA-THCA"),
    "DTC":     ("DTC",        "Differentiated Thyroid Cancer = PTC ∪ FVPTC ∪ FTC",                                                     "WHO histology"),
    "PTC":     ("PTC / FVPTC / FTC", "Papillary / Follicular-variant PTC / Follicular Thyroid Carcinoma",                              "WHO histology"),
    "PDTC":    ("PDTC",       "Poorly Differentiated Thyroid Carcinoma (Turin criteria)",                                              "WHO histology"),
    "ATC":     ("ATC",        "Anaplastic Thyroid Carcinoma",                                                                          "WHO histology"),
    "ADV":     ("\"advanced\"", "ATC ∪ PDTC (sweep convention)",                                                                       "analysis-defined"),
    "NBNR":    ("NBNR",       "No-BRAF No-RAS Korean K2 cohort cluster",                                                               "K2 / PRJEB11591"),
    "PTHT":    ("PTC+HT",     "Papillary thyroid carcinoma with concurrent Hashimoto's thyroiditis (Paper 2 main scope)",              "GSE286332 + TCGA transfer"),
    "TERT":    ("TERT⁺",      "TERT promoter mutation positive (cBioPortal `thca_tcga_pub`, n=36 recovered); Cox HR=4.33, p=4.9e-6",   "TCGA-THCA"),
}

PAPER_NUM = [
    ("1",     "Paper 1", "DM1 molecular dark matter (8-gene + DM1)",                                  "Active main"),
    ("2",     "Paper 2", "H&E → DM1 / pathology projection / TCGA validation; HT-isolated",           "Active"),
    ("3",     "Paper 3", "ICI vulnerability dark thyroid cancer",                                     "Track A FROZEN; Track B-lite findings 2026-05-06"),
    ("4",     "Paper 4", "Korean GD HLA / Pan-Asian (backlog)",                                       "Backlog (4/4 gating)"),
    ("5-12",  "Paper 5–12", "Atlas / pan-cancer / network / phaseB drug platform / etc.",             "Per-paper hub pages"),
]

# Section F — per-figure data sourcing (Paper 1 manuscript v8)
# (figure_id, caption_short, dataset, n)
FIGURES = [
    # Main figures
    ("Fig 1A", "TIERA67 → 8-gene curation Sankey",                 "(schematic; analyst-curated, no data)",                "—"),
    ("Fig 1B", "UMAP DM1/DM2 KMeans on 8-gene space",              "TCGA-THCA",                                            "504"),
    ("Fig 1C", "BRAF V600E vs WT mRNA neutrality + driver AUC",    "TCGA-THCA RNA-seq + cBioPortal mut",                   "273 V600E + 182 WT"),
    ("Fig 1D", "Pan-genome ARI ladder (8 / TIERA67 / top-5000)",   "TCGA-THCA",                                            "504"),
    ("Fig 2A", "Xing 2014 dark-matter rescue Sankey 131/180",      "TCGA-THCA + Xing 2014 list",                           "180 → 131 rescued"),
    ("Fig 2B", "DM1 prevalence ancestry (28.4% vs 37.8%)",         "TCGA + K2 + Lee/GSE213647 + GSE286332",                "504 + 874"),
    ("Fig 2C", "KM OS DM1 vs DM2 within Xing dark matter",         "TCGA-THCA dark-matter",                                "180"),
    ("Fig 3A", "Per-patient Spearman r tumor↔normal thyrocyte",    "GSE184362 (Pu 2021)",                                  "6 patients / 158K cells"),
    ("Fig 3B", "UMAP thyrocyte cluster DM1 score",                 "GSE193581 (Lu 2023)",                                  "23 samples"),
    ("Fig 3C", "Author-independence cross-cohort sc",              "GSE241184 (Phase 1) vs GSE184362",                     "1 vs 6"),
    ("Fig 4A", "8-cell stack BRAF × TERT × DM",                    "TCGA-THCA + cBioPortal `thca_tcga_pub` (TERT⁺ n=36)",  "504"),
    ("Fig 4B", "KM OS BRAF_TERT⁺ vs BRAF/TERT-neg dark matter × DM","TCGA-THCA",                                           "504"),
    ("Fig 4C", "OTHER_TERT⁺ N=4 caveat box",                       "TCGA-THCA",                                            "4"),
    ("Fig 5A", "Korean 8-gene DM-score distribution",              "GSE213647 (Lee 2024)",                                 "632"),
    ("Fig 5B", "K2 NBNR DM cluster + vascular invasion",           "K2 / PRJEB11591",                                      "260"),
    ("Fig 5C", "FFPE vs FF concordance KS p=0.44",                 "TCGA-THCA + ⚠ FFPE source verify",                     "—"),
    ("Fig 6A", "Forest Cox HR DM1 vs DM2",                         "TCGA-THCA + MSK-IMPACT (Landa 2016)",                  "504 + 117"),
    ("Fig 6B", "I² heterogeneity panel",                           "(derived from 6A)",                                    "—"),
    ("Fig 6C", "DM1 sub-A vs sub-B Cox",                           "TCGA-THCA DM1",                                        "91 (4 events)"),
    ("Fig 7A", "DM1 sub-A vs sub-B silhouette",                    "TCGA-THCA DM1",                                        "72 + 19"),
    ("Fig 7B", "Fusion-positivity DM1 vs DM2 (OR 7.41)",           "TCGA-THCA + cBioPortal SV",                            "489 SV-tested"),
    ("Fig 7C", "Fusion partner stack within DM1",                  "TCGA-THCA fusion+ DM1",                                "63"),
    ("Fig 7D", "DM1 sub-A vs sub-B age/stage/CD8/IFNγ/checkpoint", "TCGA-THCA DM1",                                        "72 + 19"),
    ("Fig 7E", "DM1 captures 81.8% TCGA RET-fusion+",              "TCGA-THCA RET subset",                                 "27/33"),
    ("Fig 8A", "HM450 promoter β heatmap (8-gene + DIO2 + SLC26A4)","TCGA-THCA HM450",                                     "503"),
    ("Fig 8B", "Mean panel β by DM cluster",                       "TCGA-THCA HM450",                                      "503"),
    ("Fig 8C", "Methylation fusion-independence within DM1",       "TCGA-THCA HM450 ∩ DM1",                                "63 vs 19"),
    # Supplementary core
    ("S1",  "8-gene heatmap sorted by P_DM1; TCGA vs K2",            "TCGA + K2",                                          "504 + 260"),
    ("S2",  "8-gene panel similarity matrix sc",                     "GSE184362 + GSE193581 + GSE241184",                  "6 + 23 + 1"),
    ("S3",  "Pan-genome ARI ladder (8/16/67/200/1000/5000)",         "TCGA-THCA",                                          "504"),
    ("S4",  "Immune-residualization Δd (1.78 → 0.87)",               "TCGA-THCA + immune signature",                       "504"),
    ("S5",  "SV missingness sensitivity",                            "TCGA-THCA + cBioPortal SV",                          "489 + 15"),
    ("S6",  "DM1 sub-B age/stage/immune detail",                    "TCGA-THCA DM1",                                      "72 + 19"),
    ("S7",  "Cross-cohort DM-score portability",                    "K2 + Lee/GSE213647 + ⚠ small Korean ref",            "260 + 632 + ?"),
    ("S8",  "Decitabine + I-131 schematic + lit meta",              "NCT00085293 + NCT01065090 (citation only)",          "n/a"),
    ("S9",  "K2 mini-index calibration diagnostic R4-4",            "K2 vs Yoo 2016 reference",                           "260"),
    # Additional supp — Part A v17 audit
    ("SA1", "4-way BRAF × TERT revalidation",                       "TCGA-THCA + cBioPortal TERT⁺",                       "504 (BRAF−/TERT⁺ n=4)"),
    ("SA2", "FFPE vs FF QC (KS p=0.44)",                            "TCGA-THCA + ⚠ FFPE source",                          "—"),
    ("SA3", "Panel-size sensitivity (8/10/12/16)",                  "TCGA-THCA",                                          "504"),
    ("SA4", "sc wrap-up (3 cohorts)",                                "GSE184362 + GSE193581 + GSE241184",                  "6 + 23 + 1"),
    ("SA5", "Pseudotime trajectory along 8-gene",                   "TCGA-THCA",                                          "504"),
    ("SA6", "MSK-IMPACT bias panel (84 PDTC + 33 ATC)",             "MSK-IMPACT (Landa 2016)",                            "117"),
    # Part B / C / D
    ("SB1", "sc UMAP overview by patient",                          "GSE184362 (Pu 2021)",                                "6"),
    ("SB2", "sc 8-gene + HLA-II + thyroid TFs + proliferation",     "GSE184362",                                          "6"),
    ("SB3", "Thyrocyte-restricted UMAP",                            "GSE184362",                                          "6"),
    ("SC1", "Per-patient r forest tumor↔normal thyrocyte",          "GSE184362",                                          "6"),
    ("SC2", "Multisite trajectory DM1/DM2 score",                   "TCGA + MSK-IMPACT + K2 + Lee/GSE213647",             "504 + 117 + 260 + 632"),
    ("SC3", "Joint sc UMAP (⚠ integration method verify)",           "GSE184362 + GSE193581 + GSE241184",                  "6 + 23 + 1"),
    ("SC4", "GSE184362 per-patient summary table",                  "GSE184362",                                          "6"),
    ("SC5", "Pooled scatter tumor vs normal",                       "GSE184362 + GSE193581",                              "6 + 23"),
    ("SC6", "Pooled scatter DM1 prob vs 8-gene",                    "GSE184362 + GSE193581",                              "6 + 23"),
    ("SC7", "K2 vs Yoo 2016 reference R4-4 audit",                   "K2 / PRJEB11591 vs Yoo 2016 panel",                  "260"),
    ("SD1", "GSE76039 differentiation transcript heatmap",          "GSE76039 (Landa 2016)",                              "~37 PDTC+ATC"),
]

# Per-paper figure subsets — which figures to display per paper
FIGS_PAPER1 = [f[0] for f in FIGURES]  # all
FIGS_NONE   = []  # paper2/3/4 pages don't show Paper-1 figures

# Section E — request-required / restricted-access data
REQUEST = {
    "DBGAP_TCGA_BAM": ("TCGA-THCA paired tumor-normal BAM", "dbGaP",                  "WES BAM",         "≈496",  "LOHHLA / HLA LOH / neoantigen (Paper 3 Module C)", "Application 4–8 wk; user 결정 대기"),
    "DBGAP_LIU":      ("phs000452 (Liu melanoma ICI)",      "dbGaP",                  "WES + RNA-seq",   "—",     "DIAL audit Tier 2 anchor (Paper 3)",               "Application 4–8 wk"),
    "DBGAP_CDR":      ("phs000178 (TCGA-CDR)",              "dbGaP",                  "clinical",        "—",     "Survival downstream",                              "Standard application"),
    "EGA_3540":       ("EGAS00001003540",                   "EGA",                    "TBD verify",      "—",     "Verify usage in manuscript",                       "DAC approval"),
    "EGA_2556":       ("EGAS00001002556",                   "EGA",                    "TBD verify",      "—",     "Verify usage",                                     "DAC approval"),
    "EGA_4845":       ("EGAD00001004845",                   "EGA",                    "TBD verify",      "—",     "Verify usage",                                     "Dataset-level access"),
    "EGA_YOO":        ("Yoo SK 2019 Korean ATC",            "EGA likely",             "RNA-seq",         "Korean ATC", "Asian dedifferentiation generalization (Paper 1)", "Accession verify + 4–8 wk"),
    "KOGES":          ("KARE / KoGES Korean SNP",           "KCDC Korean Genome",     "~10K SNP / ImmunoChip-like", "—", "Korean GD HLA imputation reference (Paper 4)", "Korean Genome 신청 + IRB"),
    "KEY158":         ("KEYNOTE-158 thyroid",               "Merck (NCT02628067)",    "raw RNA-seq 비공개", "—",   "K1 해제 candidate (response prediction claim)",     "Merck IIS — months; possibly never"),
    "KEY028":         ("KEYNOTE-028 thyroid",               "Merck (NCT02054806)",    "raw RNA-seq 비공개", "—",   "동일",                                              "Merck IIS"),
    "SPART":          ("Spartalizumab ATC",                 "Novartis (NCT02404441)", "raw 비공개",      "—",     "K1 해제 candidate",                                 "Novartis ETOP"),
    "NIVO_IPI":       ("nivo+ipi aggressive thyroid",       "(NCT03246958)",          "raw 비공개",      "—",     "Citation only",                                    "Trial PI request"),
    "ATEZO_TT":       ("Atezolizumab matched-TT ATC",       "Genentech (NCT03181100)","raw 비공개",      "—",     "Citation only",                                    "Genentech IIS"),
    "DIERKS":         ("Dierks lenva+pembro / Sehgal / Cabanillas thyroid ICI", "publications", "clinical only", "—", "Citation only", "Corresponding author request"),
    "BUNDOUT":        ("Bundang SNUH Graves' BTC + retrospective FFPE", "Korean institutional", "outreach 단계 0%", "—", "Paper 4 backlog activation; Paper 2 FFPE Korean", "Yu 교수 outreach (6+ wk wait)"),
    "YONS":           ("연세 / 서울대 / 삼성 thyroid retrospective FFPE", "Korean institutional", "미접촉",      "—",  "Paper 2 Korean FFPE pack",                          "Yu 교수 협력"),
    "CHEN18":         ("Chen 2018 Han Chinese GD HLA fine-mapping", "PMC6161647 supp 또는 dbGaP", "summary stats", "—", "Paper 4 Pan-Asian", "Supp 추출 또는 dbGaP"),
    "POZRAW":         ("Pozdeyev 2018 raw (n=779 advanced)", "corresponding author",  "mutation + clinical", "779", "Disc 3.1 raw 보강", "PI request"),
}

# ---------------------------------------------------------------------------
# Per-page mapping — what each page should display
# Keys are filenames; values are dicts with cohorts/scores/groups/papers (lists of IDs)
# Special: empty dict = full dictionary (hub-level pages).
# ---------------------------------------------------------------------------
ALL = {"cohorts": list(COHORTS), "scores": list(SCORES), "groups": list(GROUPS), "papers": ["1","2","3","4","5-12"], "requests": list(REQUEST), "figures": FIGS_PAPER1}

# default request subset per paper (added on top of cohort mapping below)
REQ_PAPER1 = ["DBGAP_TCGA_BAM","HM450" if False else "EGA_YOO","CHEN18","POZRAW"]  # we'll also append request rows specific to a paper inside build_block
REQ_PAPER1 = ["DBGAP_TCGA_BAM","EGA_YOO","POZRAW","EGA_3540","EGA_2556","EGA_4845"]
REQ_PAPER2 = ["BUNDOUT","YONS"]
REQ_PAPER3 = ["DBGAP_TCGA_BAM","DBGAP_LIU","KEY158","KEY028","SPART","NIVO_IPI","ATEZO_TT","DIERKS"]
REQ_PAPER4 = ["BUNDOUT","KOGES","CHEN18"]
REQ_PAPERS_ALL = list(REQUEST)

PAGE_MAP: dict[str, dict] = {
    # Hub-level — full dictionary
    "index.html": ALL,
    "master_view.html": ALL,
    "master_situation_2026_05_04.html": ALL,
    "trajectory.html": ALL,

    # Paper 1
    "paper1.html":                       {"cohorts":["TCGA","HM450","CBIO","MSK","K2","LEE24","76039","33630","65144","29265","53157","60542","184362","193581","241184"], "scores":["8GENE","RAI8","DM1S","NONO","TDS","MECH","TIERA","COH","DIAL"], "groups":["DM1","DM1B","DM2","DTC","PTC","PDTC","ATC","ADV","TERT"], "papers":["1"], "requests": REQ_PAPER1},
    "portfolio_paper1.html":             {"cohorts":["TCGA","HM450","CBIO","MSK","K2","LEE24","76039","33630","65144","29265","53157","184362","193581"], "scores":["8GENE","RAI8","DM1S","NONO","TDS","MECH","COH"], "groups":["DM1","DM1B","DM2","DTC","PTC","ATC","ADV","TERT"], "papers":["1"], "requests": REQ_PAPER1},
    "paper1_nature_cancer_board.html":   {"cohorts":["TCGA","HM450","CBIO","146003","MSK","76039","33630","65144","29265","53157","LEE24","184362","193581"], "scores":["8GENE","RAI8","DM1S","NONO","TDS","COH"], "groups":["DM1","DM1B","DM2","ATC","ADV","TERT"], "papers":["1"], "requests": REQ_PAPER1},
    "gpl570_validation_pack.html":       {"cohorts":["33630","65144","29265","53157","60542","126698"],     "scores":["RAI8","DM1S","NONO","TDS","MECH","COH"],                       "groups":["DTC","PTC","PDTC","ATC","ADV"], "papers":["1"]},
    "intro_strategy_audit.html":         {"cohorts":["TCGA","MSK","76039","LEE24"],                         "scores":["8GENE","RAI8","DM1S","TIERA"],                                  "groups":["DM1","DM2","ATC","ADV"], "papers":["1","2"]},
    "portfolio_paper0.html":             {"cohorts":["TCGA","HM450","CBIO"],                                "scores":["8GENE","DM1S","DIAL"],                                          "groups":["DM1","DM2"], "papers":["1"]},

    # Paper 2 (HT-isolated)
    "paper2a.html":                      {"cohorts":["286332","TCGA","K2","LEE24"],   "scores":["8GENE","DM1S","HLA","COH"], "groups":["DM1","DM1B","DM2","PTHT","NBNR"], "papers":["2"], "requests": REQ_PAPER2},
    "paper2b.html":                      {"cohorts":["286332","TCGA","K2","LEE24"],   "scores":["8GENE","HLA","COH"],        "groups":["DM2","PTHT","NBNR"], "papers":["2"], "requests": REQ_PAPER2},
    "paper2_hla.html":                   {"cohorts":["AFND","K2","BUND"],     "scores":["HLA","COH"],                "groups":["PTHT"], "papers":["2","4"], "requests": REQ_PAPER2 + REQ_PAPER4},
    "paper2_master_view.html":           {"cohorts":["286332","TCGA","K2","LEE24","AFND"], "scores":["8GENE","DM1S","HLA","COH"], "groups":["DM1","DM2","PTHT","NBNR"], "papers":["2"], "requests": REQ_PAPER2},
    "portfolio_paper2.html":             {"cohorts":["286332","TCGA","K2","LEE24","AFND"], "scores":["8GENE","HLA","COH"],   "groups":["DM2","PTHT","NBNR"], "papers":["2"], "requests": REQ_PAPER2},
    "pillar1_v2_spec.html":              {"cohorts":["AFND","K2"],            "scores":["HLA","COH"],                "groups":["PTHT"], "papers":["2","4"], "requests": REQ_PAPER2 + REQ_PAPER4},
    "ht_vs_gd_audit.html":               {"cohorts":["AFND","K2","286332"],   "scores":["HLA","COH"],                "groups":["PTHT"], "papers":["2","4"]},
    "yu_send_package.html":              {"cohorts":["AFND","K2","BUND","286332"], "scores":["HLA","COH"],           "groups":["PTHT"], "papers":["2","4"], "requests": REQ_PAPER2 + REQ_PAPER4},

    # Paper 3 (ICI dark)
    "paper3.html":                       {"cohorts":["TCGA","MSK","HUGO","78220","91061","151179","76039"], "scores":["TDIFF","MYEL","HLA","DIAL","COH"], "groups":["ATC","ADV","DM1"], "papers":["3"], "requests": REQ_PAPER3},
    "portfolio_paper3.html":             {"cohorts":["TCGA","HUGO","78220","91061","151179"], "scores":["TDIFF","MYEL","DIAL","COH"],       "groups":["ATC","ADV"], "papers":["3"], "requests": REQ_PAPER3},
    "cancer_vaccine_agent.html":         {"cohorts":["TCGA","HUGO","78220","91061"],          "scores":["TDIFF","MYEL","HLA","COH"],        "groups":["ATC","DM1","PTHT"], "papers":["3"], "requests": REQ_PAPER3},

    # Paper 4 (Korean GD HLA backlog)
    "paper4.html":                       {"cohorts":["AFND","K2","BUND"],     "scores":["HLA","COH"], "groups":["PTHT"], "papers":["4"], "requests": REQ_PAPER4},
    "paper4_gd_hla.html":                {"cohorts":["AFND","K2","BUND"],     "scores":["HLA","COH"], "groups":["PTHT"], "papers":["4"], "requests": REQ_PAPER4},
    "portfolio_paper4.html":             {"cohorts":["AFND","K2","BUND"],     "scores":["HLA","COH"], "groups":["PTHT"], "papers":["4"], "requests": REQ_PAPER4},
    "paperG.html":                       {"cohorts":["AFND","K2","BUND"],     "scores":["HLA","COH"], "groups":["PTHT"], "papers":["4"], "requests": REQ_PAPER4},
    "portfolio_paperG.html":             {"cohorts":["AFND","K2","BUND"],     "scores":["HLA","COH"], "groups":["PTHT"], "papers":["4"], "requests": REQ_PAPER4},
    "hla_advisor_packet.html":           {"cohorts":["AFND","K2","BUND"],     "scores":["HLA","COH"], "groups":["PTHT"], "papers":["4"], "requests": REQ_PAPER4},
    "hla_beta_plan.html":                {"cohorts":["AFND","K2","BUND"],     "scores":["HLA","COH"], "groups":["PTHT"], "papers":["4"], "requests": REQ_PAPER4},
    "hla_situation.html":                {"cohorts":["AFND","K2","BUND"],     "scores":["HLA","COH"], "groups":["PTHT"], "papers":["4"], "requests": REQ_PAPER4},
    "k2_bd_report.html":                 {"cohorts":["K2","AFND","BUND"],     "scores":["HLA","COH","DIAL"], "groups":["NBNR","PTHT"], "papers":["4"], "requests": REQ_PAPER4},
    "k2_t1_cards.html":                  {"cohorts":["K2","AFND"],            "scores":["HLA","COH"],        "groups":["NBNR","PTHT"], "papers":["4"], "requests": REQ_PAPER4},

    # Paper 5–7 (atlas / etc.)
    "paper5.html":                       {"cohorts":["TCGA","MSK","76039"],   "scores":["8GENE","DM1S","TDS","COH"], "groups":["DM1","DM2","ATC","ADV"], "papers":["1","5-12"]},
    "paper6.html":                       {"cohorts":["TCGA","MSK","76039"],   "scores":["8GENE","DM1S","TDS","COH"], "groups":["DM1","DM2","ATC","ADV"], "papers":["5-12"]},
    "paper7.html":                       {"cohorts":["TCGA","MSK","76039"],   "scores":["8GENE","DM1S","TDS","COH"], "groups":["DM1","DM2","ATC","ADV"], "papers":["5-12"]},

    # Paper 9 (CDK7/THZ1)
    "paper9.html":                       {"cohorts":["TCGA","MSK"],           "scores":["8GENE","DM1S","COH"], "groups":["ATC","DM1","ADV"], "papers":["5-12"]},
    "paper9_first_pass.html":            {"cohorts":["TCGA","MSK"],           "scores":["8GENE","DM1S","COH"], "groups":["ATC","DM1","ADV"], "papers":["5-12"]},
    "paper9_plan.html":                  {"cohorts":["TCGA","MSK"],           "scores":["8GENE","DM1S","COH"], "groups":["ATC","DM1","ADV"], "papers":["5-12"]},
    "portfolio_paper9.html":             {"cohorts":["TCGA","MSK"],           "scores":["8GENE","DM1S","COH"], "groups":["ATC","DM1","ADV"], "papers":["5-12"]},

    # Paper 10–12 (atlas / pancan / network)
    "paper10_atlas.html":                {"cohorts":["TCGA","MSK","76039"],   "scores":["8GENE","DM1S","TDS","COH"], "groups":["DM1","DM2","ATC","ADV"], "papers":["5-12"]},
    "paper11_pancancer.html":            {"cohorts":["TCGA"],                 "scores":["8GENE","DM1S","COH"],       "groups":["DM1","DM2"],             "papers":["5-12"]},
    "paper12_network.html":              {"cohorts":["TCGA"],                 "scores":["8GENE","DM1S","MECH","COH"], "groups":["DM1","DM2"],            "papers":["5-12"]},

    # PaperB / N / drug platform
    "paperB.html":                       {"cohorts":["TCGA","HUGO","78220","91061"], "scores":["8GENE","TDIFF","COH"], "groups":["DM1","ATC","ADV"], "papers":["3","5-12"], "requests": REQ_PAPER3},
    "phaseB_drug_platform.html":         {"cohorts":["TCGA","HUGO","78220","91061"], "scores":["8GENE","TDIFF","COH"], "groups":["DM1","ATC","ADV"], "papers":["3","5-12"], "requests": REQ_PAPER3},
    "portfolio_paperB.html":             {"cohorts":["TCGA","HUGO","78220","91061"], "scores":["8GENE","TDIFF","COH"], "groups":["DM1","ATC","ADV"], "papers":["3","5-12"], "requests": REQ_PAPER3},
    "paperN.html":                       {"cohorts":["TCGA","HM450","146003"],"scores":["8GENE","DM1S","COH"],  "groups":["DM1","DM2"],       "papers":["5-12"]},
    "portfolio_paperN.html":             {"cohorts":["TCGA","HM450","146003"],"scores":["8GENE","DM1S","COH"],  "groups":["DM1","DM2"],       "papers":["5-12"]},

    # Lumenix demo (spatial)
    "lumenix_demo_report.html":          {"cohorts":["250521","TCGA"],        "scores":["8GENE","DM1S","DIAL","COH"], "groups":["DM1","ATC","ADV"], "papers":["1","5-12"]},
}


# ---------------------------------------------------------------------------
# HTML generation
# ---------------------------------------------------------------------------
TH = ("border:1px solid #c9c9c9;padding:5px 8px;background:#efe7e0")
TD = ("border:1px solid #c9c9c9;padding:5px 8px")

def _row(*cells: str) -> str:
    return "<tr>" + "".join(f'<td style="{TD}">{c}</td>' for c in cells) + "</tr>"

def _table(headers: list[str], rows: list[tuple[str, ...]]) -> str:
    if not rows:
        return ""
    head = "<tr>" + "".join(f'<th style="{TH};text-align:left">{h}</th>' for h in headers) + "</tr>"
    body = "\n".join(_row(*r) for r in rows)
    return f'<table style="border-collapse:collapse;width:100%;font-size:12px"><thead>{head}</thead><tbody>{body}</tbody></table>'

def _section(title: str, table_html: str) -> str:
    if not table_html:
        return ""
    return (
        f'<h4 style="margin:14px 0 6px;font-family:\'JetBrains Mono\',monospace;font-size:11px;letter-spacing:.08em;color:#7B1F2A">{title}</h4>'
        + table_html
    )

def build_block(filename: str) -> str:
    cfg = PAGE_MAP.get(filename, ALL)
    cohort_rows  = [COHORTS[k] for k in cfg.get("cohorts", []) if k in COHORTS]
    score_rows   = [SCORES[k]  for k in cfg.get("scores",  []) if k in SCORES]
    group_rows   = [GROUPS[k]  for k in cfg.get("groups",  []) if k in GROUPS]
    paper_rows   = [(p[1], p[2], p[3]) for p in PAPER_NUM if p[0] in cfg.get("papers", [])]
    request_rows = [REQUEST[k] for k in cfg.get("requests", []) if k in REQUEST]
    fig_ids      = cfg.get("figures", [])
    figure_rows  = [f for f in FIGURES if f[0] in fig_ids]

    is_full = (cfg is ALL)
    scope = "전체 portfolio 공통 사전" if is_full else f"이 페이지에 사용된 정의만 표시 ({filename})"

    a = _section("A. COHORTS / DATASETS",                     _table(["ID","Source","Modality","n","Use"], cohort_rows))
    b = _section("B. SCORES / GENE PANELS",                   _table(["Score","Definition","Genes / formula","Higher ="], score_rows))
    c = _section("C. PATIENT GROUPS / LABELS",                _table(["Label","Meaning","Source"], group_rows))
    d = _section("D. PAPER NUMBERING (CANONICAL 2026-05-04)", _table(["#","Title / scope","Status"], paper_rows))
    e = _section("E. REQUEST-REQUIRED / RESTRICTED-ACCESS DATA (★ 결정 대기)",
                 _table(["Dataset","Authority","Modality","n","Use","Path / lag"], request_rows))

    note_e = '<p style="margin:8px 0 0;font-size:11px;color:#7B1F2A">★ 위 항목은 <b>아직 다운로드 불가</b> — dbGaP / EGA / KoGES application, PI/pharma request, 또는 Korean institutional outreach 필요. 분석 명령 아니라 registry. 자세한 결정 행렬은 <code>project/reports/2026_05_07_data_audit_and_request_registry.md</code>.</p>' if request_rows else ''

    # Section F — per-figure data sourcing (collapsed by default to reduce visual load)
    if figure_rows:
        f_table = _table(["Figure","Caption (short)","Dataset","n"], figure_rows)
        f_section = (
            '<h4 style="margin:14px 0 6px;font-family:\'JetBrains Mono\',monospace;font-size:11px;letter-spacing:.08em;color:#7B1F2A">F. PER-FIGURE DATA SOURCING (Paper 1 manuscript v8)</h4>'
            + '<details style="margin:0 0 8px;"><summary style="cursor:pointer;font-size:11px;color:#7B1F2A;font-family:\'JetBrains Mono\',monospace">▶ click to expand figure-by-figure dataset table (Fig 1–8 + S1–S9 + SA/SB/SC/SD)</summary>'
            + f'<div style="margin-top:10px">{f_table}</div></details>'
            + '<p style="margin:8px 0 0;font-size:11px;color:#777">자세한 figure 단위 audit + verify-필요 ⚠ 항목은 <code>project/reports/2026_05_07_per_figure_data_audit.md</code>.</p>'
        )
    else:
        f_section = ''

    return f"""<!-- DATA_DEFINITIONS_BLOCK_START v4 2026-05-07 {filename} -->
<section id="data-definitions" style="margin:0 0 24px;padding:18px 22px;border:1px solid #c9c9c9;border-left:4px solid #7B1F2A;background:#fbf7f4;font-family:'Newsreader','Noto Sans KR',serif;">
  <details open style="margin:0;">
    <summary style="cursor:pointer;font-family:'JetBrains Mono',monospace;font-size:12px;letter-spacing:.06em;color:#7B1F2A;font-weight:700;">데이터 정의서 / DATA DEFINITIONS v4 — {scope} (open by default — click to fold)</summary>
    <div style="margin-top:14px;font-size:13px;line-height:1.55;">
      <p style="margin:0 0 10px;color:#444">아래 표는 <b>이 페이지가 실제로 참조하는</b> cohort / score / group / paper 번호 / 다운로드 불가 데이터 / <b>per-figure data sourcing</b> 정의입니다. v4 (2026-05-07) — Section F per-figure audit 추가 (Paper 1 manuscript v8 Fig 1–8 + S1–S9 + SA/SB/SC/SD).</p>
      {a}{b}{c}{d}{e}{note_e}{f_section}
      <p style="margin:12px 0 0;font-size:11px;color:#777;font-family:'JetBrains Mono',monospace">v4 · 2026-05-07 · per-page filtered · source: <code>_inject_dict.py</code> · audits: <code>2026_05_07_data_audit_and_request_registry.md</code> · <code>2026_05_07_per_figure_data_audit.md</code></p>
    </div>
  </details>
</section>
<!-- DATA_DEFINITIONS_BLOCK_END v4 -->
"""


# ----------------------------------------------------------------------------
# Attach Paper 1 figure list to EVERY page (only Paper 1 has finalized
# figure captions in manuscript v8; Paper 2/3/4 figures are TBD).
# ----------------------------------------------------------------------------
for _name, _cfg in PAGE_MAP.items():
    if _cfg is ALL:
        continue
    _cfg.setdefault("figures", FIGS_PAPER1)


START_RE = re.compile(
    r"<!-- DATA_DEFINITIONS_BLOCK_START.*?<!-- DATA_DEFINITIONS_BLOCK_END[^>]*-->\s*",
    re.DOTALL,
)
BODY_RE = re.compile(r"(<body[^>]*>)", re.IGNORECASE)

SKIP_NAMES = {"_data_dictionary.html"}

def process(p: Path) -> str:
    src = p.read_text(encoding="utf-8", errors="replace")
    cleaned = START_RE.sub("", src)
    m = BODY_RE.search(cleaned)
    if not m:
        return "no-body-tag"
    block = build_block(p.name)
    end = m.end()
    out = cleaned[:end] + "\n" + block + cleaned[end:]
    if out == src:
        return "unchanged"
    p.write_text(out, encoding="utf-8")
    return "updated"

def main():
    files = sorted(ROOT.glob("*.html"))
    counts = {}
    for f in files:
        if f.name in SKIP_NAMES or f.name.endswith(".bak-immune-readiness-2026-05-06"):
            counts["skipped"] = counts.get("skipped", 0) + 1
            continue
        scope = "FULL" if f.name in {"index.html","master_view.html","master_situation_2026_05_04.html","trajectory.html"} else (
            "MAPPED" if f.name in PAGE_MAP else "FULL-fallback"
        )
        status = process(f)
        counts[status] = counts.get(status, 0) + 1
        print(f"  {status:12s} [{scope:14s}]  {f.name}")
    print("\nsummary:", counts)

if __name__ == "__main__":
    main()
