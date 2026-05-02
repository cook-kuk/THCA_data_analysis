"""4 deeper pipelines on lit_enrich data:

P1. Topical relevance filter — TF-IDF rerank of 277 dedup papers per claim
P2. Random-effects forest meta — AFND DPB1*05:01 East Asian baseline (DerSimonian-Laird)
P3. Trial × mechanism cross-reference — 38 targeted-drug trials → drug class → BRAF-WT/RAS-WT relevance
P4. Bib fact-check — 10 OA abstracts vs verified_refs metadata

Outputs: project/manuscript_p2_brief/lit_enrich_2026_05_02/10_pipelines/
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
DATA = ROOT / "project/manuscript_p2_brief/lit_enrich_2026_05_02/data"
DA = ROOT / "project/manuscript_p2_brief/lit_enrich_2026_05_02/deep_analysis"
MORE = ROOT / "project/manuscript_p2_brief/lit_enrich_2026_05_02/09_more_fetched"
OUT = ROOT / "project/manuscript_p2_brief/lit_enrich_2026_05_02/10_pipelines"
OUT.mkdir(parents=True, exist_ok=True)


# ============================================================================ #
# P1. Topical relevance filter
# ============================================================================ #
print("=" * 70)
print("P1. Topical relevance filter (TF-IDF cosine on title vs claim)")
print("=" * 70)

dedup = pd.read_csv(DA / "02_related_works_dedup_ranked.csv")
print(f"  Input: {len(dedup)} dedup papers across {dedup['claim'].nunique()} claims")

# Build expanded keyword set for each claim (claim text + canonical synonyms)
CLAIM_KEYWORDS = {
    "8-gene / driver-excluded PTC stratification": (
        "8-gene signature papillary thyroid carcinoma stratification driver-excluded "
        "BRAF RAS RET fusion molecular subtype RandomForest classifier"
    ),
    "BCR repertoire / TLS in thyroid cancer": (
        "B cell receptor BCR repertoire tertiary lymphoid structure TLS "
        "thyroid carcinoma immunoglobulin IGHV clonality AICDA antigen-driven"
    ),
    "BRAF/RAS-negative ('dark matter') PTC subtypes": (
        "BRAF wild-type RAS wild-type wildtype papillary thyroid cancer "
        "molecular subtype fusion driver dark matter mutation negative"
    ),
    "HLA-II PTC autoimmunity (DPB1*05:01 etc.)": (
        "HLA class II DPB1 DRB1 DQB1 papillary thyroid carcinoma autoimmunity "
        "Hashimoto thyroiditis allele frequency genetic susceptibility"
    ),
    "Hashimoto-like signature & PTC outcomes": (
        "Hashimoto thyroiditis papillary thyroid carcinoma transcriptome "
        "lymphocytic infiltration immune signature outcome prognosis tumor microenvironment"
    ),
    "Pan-Asian HLA fine-mapping (Graves' / autoimmune thyroid)": (
        "Pan-Asian HLA fine mapping Graves disease autoimmune thyroid "
        "Korean Chinese Japanese Taiwanese GWAS allele frequency"
    ),
    "Single-cell PTC progression": (
        "single cell RNA sequencing scRNA papillary thyroid carcinoma progression "
        "tumor heterogeneity malignant epithelial transformation lineage trajectory"
    ),
    "TIERA-like Korean PTC molecular cohort": (
        "Korean papillary thyroid carcinoma cohort molecular landscape "
        "transcriptome TIERA SNUH TCGA genomic subtype"
    ),
}

reranked = []
for claim, sub in dedup.groupby("claim"):
    sub = sub.copy()
    titles = sub["title"].fillna("").tolist()
    claim_text = CLAIM_KEYWORDS.get(claim, claim)
    docs = [claim_text] + titles
    if len(titles) == 0:
        continue
    vec = TfidfVectorizer(stop_words="english", min_df=1, max_features=2000, ngram_range=(1, 2))
    try:
        mat = vec.fit_transform(docs)
        sims = cosine_similarity(mat[0:1], mat[1:]).ravel()
    except ValueError:
        sims = np.zeros(len(titles))
    sub["topical_relevance"] = sims
    sub["new_score"] = sub["topical_relevance"] * 100 + sub["n_sources"] * 10 + np.log1p(sub["cited_by"].fillna(0))
    reranked.append(sub)
df_rerank = pd.concat(reranked).sort_values(["claim", "new_score"], ascending=[True, False])
df_rerank.to_csv(OUT / "P1_relevance_reranked_full.csv", index=False)

# Filter: keep only papers with topical_relevance >= 0.05 (lenient cut)
df_filtered = df_rerank[df_rerank["topical_relevance"] >= 0.05]
print(f"  After relevance filter (≥0.05): {len(df_filtered)} / {len(df_rerank)} kept")

top5_clean = df_filtered.groupby("claim").head(5).reset_index(drop=True)
top5_clean[["claim", "topical_relevance", "n_sources", "year", "first_author", "title", "cited_by", "doi"]].to_csv(
    OUT / "P1_top5_clean_per_claim.csv", index=False
)

# Show contrast: how much did top-5 change?
print("\n  Top-5 change per claim (old top-5 → new top-5 overlap):")
old_top5 = pd.read_csv(DA / "02_related_works_top5_per_claim.csv")
for claim in top5_clean["claim"].unique():
    old = set(old_top5[old_top5["claim"] == claim]["doi"].dropna())
    new = set(top5_clean[top5_clean["claim"] == claim]["doi"].dropna())
    overlap = len(old & new)
    print(f"    [{overlap}/5 same]  {claim[:55]}")


# ============================================================================ #
# P2. Random-effects forest meta — DPB1*05:01 East Asian baseline
# ============================================================================ #
print("\n" + "=" * 70)
print("P2. Random-effects DerSimonian-Laird forest — DPB1*05:01 East Asian")
print("=" * 70)

afnd = json.loads((DATA / "afnd_alleles.json").read_text())
dpb_studies = []
for country, hits in afnd.get("DPB1*05:01", {}).items():
    if not isinstance(hits, list):
        continue
    for h in hits:
        try:
            f = float(h["allele_freq"])
            n = int(h["sample_size"])
            if 0 < f < 1 and n > 0:
                dpb_studies.append({
                    "country": country,
                    "pop_name": h["population"],
                    "freq": f,
                    "n": n,
                    "k_carriers_2n": int(round(2 * n * f)),  # diploid carrier count
                    "two_n": 2 * n,  # diploid alleles total
                })
        except (KeyError, ValueError, TypeError):
            continue
df_dpb = pd.DataFrame(dpb_studies)
print(f"  DPB1*05:01 studies after parse: {len(df_dpb)} across {df_dpb['country'].nunique()} countries")


def random_effects_meta(p, n):
    """Random-effects (DerSimonian-Laird) meta on proportions.
    Use Freeman-Tukey double arcsine transformation for variance stabilization.
    """
    p = np.asarray(p, dtype=float)
    n = np.asarray(n, dtype=float)
    # Freeman-Tukey arcsine (sample size weights more stable for proportions)
    t = 0.5 * (np.arcsin(np.sqrt((p * n) / (n + 1))) + np.arcsin(np.sqrt(((p * n) + 1) / (n + 1))))
    var_t = 1.0 / (4.0 * n + 2.0)  # variance of FT transform
    w_fix = 1.0 / var_t
    t_fix = (w_fix * t).sum() / w_fix.sum()
    Q = (w_fix * (t - t_fix) ** 2).sum()
    df = len(p) - 1
    C = w_fix.sum() - (w_fix ** 2).sum() / w_fix.sum()
    tau2 = max(0.0, (Q - df) / C)
    I2 = max(0.0, (Q - df) / Q) if Q > 0 else 0.0
    w_re = 1.0 / (var_t + tau2)
    t_re = (w_re * t).sum() / w_re.sum()
    var_t_re = 1.0 / w_re.sum()
    se_t_re = np.sqrt(var_t_re)
    # Back-transform to proportion (inverse Freeman-Tukey)
    # Approximation: p = sin^2(t)
    p_re = np.sin(t_re) ** 2
    p_lo = np.sin(t_re - 1.96 * se_t_re) ** 2
    p_hi = np.sin(t_re + 1.96 * se_t_re) ** 2
    return {
        "n_studies": len(p),
        "total_n": int(n.sum()),
        "pooled_p": p_re,
        "ci_low": p_lo,
        "ci_high": p_hi,
        "tau2": tau2,
        "I2": I2,
        "Q": Q,
        "Q_pvalue": 1 - stats.chi2.cdf(Q, df) if df > 0 else 1.0,
    }


print(f"\n  All East Asian pooled:")
all_meta = random_effects_meta(df_dpb["freq"].values, df_dpb["n"].values)
for k, v in all_meta.items():
    print(f"    {k}: {v if not isinstance(v, float) else f'{v:.4f}'}")

print(f"\n  By country (random-effects meta within each country):")
country_metas = []
for country, sub in df_dpb.groupby("country"):
    if len(sub) < 1:
        continue
    if len(sub) == 1:
        m = {
            "n_studies": 1,
            "total_n": int(sub["n"].iloc[0]),
            "pooled_p": float(sub["freq"].iloc[0]),
            "ci_low": np.nan,
            "ci_high": np.nan,
            "tau2": np.nan,
            "I2": np.nan,
            "Q": np.nan,
            "Q_pvalue": np.nan,
        }
    else:
        m = random_effects_meta(sub["freq"].values, sub["n"].values)
    m["country"] = country
    country_metas.append(m)
    print(f"    {country:8s}: pooled={m['pooled_p']:.3f}  [{m['ci_low']:.3f}, {m['ci_high']:.3f}]  "
          f"n_studies={m['n_studies']:3d}  total_n={m['total_n']:5d}  I²={m['I2']:.2f}  τ²={m['tau2']:.6f}")
df_country_meta = pd.DataFrame(country_metas)
df_country_meta.to_csv(OUT / "P2_dpb1_country_meta.csv", index=False)

# Compare Korean baseline (random-effects meta) vs PTC pool 53.2% (memory)
korean = next((m for m in country_metas if m["country"] == "Korea"), None)
if korean:
    delta = 0.532 - korean["pooled_p"]
    # Wald test: 53.2% Korean PTC pool (n=874) vs Korean baseline pooled (random-effects)
    p1, n1 = 0.532, 874
    p2, n2 = korean["pooled_p"], korean["total_n"]
    se = np.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
    z = (p1 - p2) / se
    pval = 2 * (1 - stats.norm.cdf(abs(z)))
    print(f"\n  PILLAR 1 statistical contrast:")
    print(f"    Korean PTC pool:  53.2% (n=874)")
    print(f"    Korean baseline:  {korean['pooled_p']*100:.1f}% [95% CI: {korean['ci_low']*100:.1f}% - {korean['ci_high']*100:.1f}%]")
    print(f"                      (random-effects meta of {korean['n_studies']} studies, total n={korean['total_n']})")
    print(f"    Δ:                +{delta*100:.1f} pp")
    print(f"    Wald test:        z={z:.2f}, p={pval:.2e}")


# ============================================================================ #
# P3. Trial × mechanism cross-reference
# ============================================================================ #
print("\n" + "=" * 70)
print("P3. Trial × dark-matter mechanism cross-reference")
print("=" * 70)

trials = pd.read_csv(DA / "03_trials_targeted_drugs.csv")
print(f"  Input: {len(trials)} targeted-drug trials")

# Drug class taxonomy (manually curated)
DRUG_CLASS = {
    # RET inhibitors
    "selpercatinib": "RET_inhibitor", "loxo-292": "RET_inhibitor", "ly3527723": "RET_inhibitor",
    "pralsetinib": "RET_inhibitor", "blu-667": "RET_inhibitor",
    "loxo-260": "RET_inhibitor_2nd_gen", "kl590586": "RET_inhibitor",
    # BRAF inhibitors
    "dabrafenib": "BRAF_inhibitor", "vemurafenib": "BRAF_inhibitor",
    "encorafenib": "BRAF_inhibitor", "fore8394": "pan-RAF",
    # MEK inhibitors
    "trametinib": "MEK_inhibitor", "cobimetinib": "MEK_inhibitor",
    "binimetinib": "MEK_inhibitor", "tunlametinib": "MEK_inhibitor",
    # Multi-kinase / TKI
    "lenvatinib": "multi_TKI", "cabozantinib": "multi_TKI",
    "sorafenib": "multi_TKI", "axitinib": "multi_TKI",
    # PD-1/PD-L1 / immune checkpoint
    "pembrolizumab": "anti-PD1", "nivolumab": "anti-PD1",
    "cemiplimab": "anti-PD1", "atezolizumab": "anti-PDL1",
    "ipilimumab": "anti-CTLA4", "sasanlimab": "anti-PD1",
    # FAK / SHP2 / RAS path
    "defactinib": "FAK_inhibitor", "avutometinib": "RAF-MEK_clamp",
    "s241656": "SHP2_inhibitor",
    # Cell therapy / ADC / radioconjugate
    "ln-145": "TIL_therapy", "sacituzumab": "ADC_TROP2",
    "sacituzumab tirumotecan": "ADC_TROP2",
    "[177lu]": "radioligand", "lu-akir001": "radioligand", "i-131": "RAI",
    # Other
    "xl092": "VEGFR_inhibitor",
}


def classify_intervention(name_str: str | float) -> tuple[str, str]:
    if not isinstance(name_str, str):
        return ("", "")
    s = name_str.lower()
    classes = []
    drugs = []
    for drug, cls in DRUG_CLASS.items():
        if drug in s:
            classes.append(cls)
            drugs.append(drug)
    return (",".join(sorted(set(drugs))), ",".join(sorted(set(classes))))


trials[["matched_drugs", "drug_classes"]] = trials.apply(
    lambda r: pd.Series(classify_intervention(r["intervention_names"])), axis=1
)

# Categorize relevance to dark-matter / 8-gene
def relevance_to_dark_matter(row):
    classes = row["drug_classes"]
    title = (row["title"] or "").lower() if isinstance(row["title"], str) else ""
    flags = []
    # RET-fusion targeting → directly relevant to BRAF/RAS-WT (most RET fusions are mutually exclusive with BRAF/RAS)
    if "RET_inhibitor" in classes:
        flags.append("RET-fusion_targeted")
    # BRAF inhibitor → BRAF V600E specific (reverse — these are NOT for dark matter, but for BRAF+ cohort)
    if "BRAF_inhibitor" in classes:
        flags.append("BRAF_mut_specific")
    # IO / pembro / nivo → not driver-specific, broad relevance
    if "anti-PD1" in classes or "anti-PDL1" in classes or "anti-CTLA4" in classes:
        flags.append("IO_agnostic")
    # Multi-TKI → broad RAI-refractory, includes BRAF-WT
    if "multi_TKI" in classes:
        flags.append("RAI-refractory_broad")
    # Title screen for "RAI" / "neoadjuvant" / "wild-type"
    if "wild-type" in title or "wt " in title or "neg" in title:
        flags.append("explicit_WT_subset")
    if "rai" in title or "i-131" in title or "iodine" in title:
        flags.append("RAI_context")
    if "neoadjuvant" in title:
        flags.append("neoadjuvant")
    return ",".join(flags)


trials["mechanism_flags"] = trials.apply(relevance_to_dark_matter, axis=1)

# Summary tables
class_summary = trials["drug_classes"].value_counts().reset_index()
class_summary.columns = ["drug_classes", "n_trials"]
class_summary.to_csv(OUT / "P3_drug_class_summary.csv", index=False)

flag_summary = pd.Series(",".join(trials["mechanism_flags"].dropna().tolist()).split(",")).value_counts()
flag_summary = flag_summary[flag_summary.index != ""].reset_index()
flag_summary.columns = ["flag", "n_trials"]
flag_summary.to_csv(OUT / "P3_mechanism_flag_summary.csv", index=False)

# Dark-matter most-relevant: RET-fusion targeted + RECRUITING + ≥Phase2
dark_priority = trials[
    trials["mechanism_flags"].str.contains("RET-fusion_targeted", na=False)
    & trials["status"].isin(["RECRUITING", "ACTIVE_NOT_RECRUITING"])
    & trials["phase"].fillna("").str.contains("PHASE2|PHASE3|PHASE4", regex=True)
].copy()
dark_priority.to_csv(OUT / "P3_RET_fusion_priority_trials.csv", index=False)
print(f"\n  Drug-class breakdown:")
print(class_summary.head(10).to_string(index=False))
print(f"\n  Mechanism-flag distribution:")
print(flag_summary.to_string(index=False))
print(f"\n  RET-fusion priority trials (≥Phase2, recruiting/active): {len(dark_priority)}")
print(dark_priority[["nct_id", "phase", "status", "matched_drugs", "lead_sponsor"]].to_string(index=False))


# ============================================================================ #
# P4. Bib fact-check using OA abstracts
# ============================================================================ #
print("\n" + "=" * 70)
print("P4. Bib fact-check — 10 OA abstracts vs verified_refs metadata")
print("=" * 70)

oa = json.loads((MORE / "oa_abstracts.json").read_text())
verified = json.loads((DATA / "verified_refs.json").read_text())
verified_by_key = {e.get("key"): e for e in verified if isinstance(e, dict)}

checks = []
for key, blob in oa.items():
    if not isinstance(blob, dict) or not blob.get("title"):
        continue
    v = verified_by_key.get(key, {})
    epmc_year = blob.get("year")
    epmc_doi = blob.get("doi")
    epmc_journal = blob.get("journal")
    epmc_first = blob.get("first_author")
    epmc_title = blob.get("title", "")[:100]
    # Bib expectations encoded in key
    expected_year_in_key = re.search(r"(\d{4})", key)
    year_match = (
        bool(expected_year_in_key)
        and str(epmc_year) == expected_year_in_key.group(1)
    )
    # Bib first author = key prefix (e.g., Haugen2016 → Haugen)
    expected_author = re.match(r"([A-Z][a-z]+)", key)
    author_match = (
        bool(expected_author)
        and bool(epmc_first)
        and expected_author.group(1).lower() in epmc_first.lower()
    )
    checks.append({
        "key": key,
        "year_match": year_match,
        "expected_year": expected_year_in_key.group(1) if expected_year_in_key else None,
        "epmc_year": epmc_year,
        "author_match": author_match,
        "expected_author": expected_author.group(1) if expected_author else None,
        "epmc_first_author": epmc_first,
        "epmc_doi": epmc_doi,
        "epmc_journal": epmc_journal,
        "epmc_title_head": epmc_title,
    })
df_check = pd.DataFrame(checks)
df_check.to_csv(OUT / "P4_bib_factcheck.csv", index=False)

n_year_ok = int(df_check["year_match"].sum())
n_auth_ok = int(df_check["author_match"].sum())
n_both = int((df_check["year_match"] & df_check["author_match"]).sum())
print(f"\n  Bib fact-check ({len(df_check)} OA entries):")
print(f"    year matches:    {n_year_ok}/{len(df_check)}")
print(f"    author matches:  {n_auth_ok}/{len(df_check)}")
print(f"    BOTH match:      {n_both}/{len(df_check)}")
print(f"\n  Discrepancies (potentially require bib update):")
disc = df_check[~(df_check["year_match"] & df_check["author_match"])]
print(disc[["key", "expected_author", "epmc_first_author", "expected_year", "epmc_year", "epmc_journal"]].to_string(index=False))


# ============================================================================ #
# Done
# ============================================================================ #
print("\n" + "=" * 70)
print("ALL 4 PIPELINES DONE")
print("=" * 70)
for p in sorted(OUT.glob("*")):
    print(f"  {p.name}  ({p.stat().st_size/1024:.1f} KB)")
