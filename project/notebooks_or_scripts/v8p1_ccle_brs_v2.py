#!/usr/bin/env python3
"""v8.1 — CCLE thyroid BRS52 + mutation-truth validation (v2).

Improvements over v1:
  - Pull CCLE mutation calls from cBioPortal for BRAF/NRAS/KRAS/HRAS
  - Use mutation truth (7 BRAF-mut + 4 RAS-mut = 11 lines) as the
    primary direction-test labels; BRS labels are reported as a
    sanity-check of the signature on cell lines.
  - Compute BRS classification accuracy on CCLE vs mutation truth
    (analogous to the 95.2 % TCGA validation).

Outputs (results/v8p1_rigor/f_ccle_validation/, overwriting v1):
  ccle_thyroid_metadata.tsv
  ccle_thyroid_mutation_truth.tsv  — 23 lines × {has_BRAF, has_RAS, label}
  ccle_brs_labels.tsv              — 13 lines with expression × BRS label
  ccle_brs_vs_mutation.tsv         — BRS vs mutation confusion (subset with both)
  ccle_concordance.tsv             — direction concordance (mutation-truth labels)
  ccle_validation_analysis.md      — narrative
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np
import pandas as pd
import requests

ROOT = Path("/opt/thyroid-dash/project")
RES  = ROOT / "results" / "v8p1_rigor"
OUT  = RES / "f_ccle_validation"
OUT.mkdir(parents=True, exist_ok=True)
LOG  = ROOT / "logs" / "v8p1_ccle_brs_v2.log"

API = "https://www.cbioportal.org/api"
STUDY = "ccle_broad_2019"
EXPR_PROFILE = "ccle_broad_2019_rna_seq_mrna"   # RPKM
MUT_PROFILE  = "ccle_broad_2019_mutations"

_log_f = open(LOG, "w")
def log(msg: str) -> None:
    stamp = time.strftime("[%Y-%m-%d %H:%M:%S]")
    print(f"{stamp} {msg}", flush=True)
    _log_f.write(f"{stamp} {msg}\n"); _log_f.flush()

log("v8.1 CCLE thyroid BRS52 + mutation-truth validation (v2)")
log("=" * 60)

# --------------------------------------------------------------------
# 1. Reuse v1 metadata + expression
# --------------------------------------------------------------------
meta_df = pd.read_csv(OUT / "ccle_thyroid_metadata.tsv", sep="\t")
expr = pd.read_csv(OUT / "ccle_thyroid_expression.tsv", sep="\t", index_col=0)
thyroid_lines = meta_df["sample_id"].tolist()
log(f"[reuse] metadata: {len(meta_df)} lines  expression: {expr.shape}")

# --------------------------------------------------------------------
# 2. Pull mutation truth for BRAF/NRAS/KRAS/HRAS
# --------------------------------------------------------------------
r = requests.post(f"{API}/genes/fetch", params={"geneIdType":"HUGO_GENE_SYMBOL"},
                  json=["BRAF","NRAS","KRAS","HRAS"], timeout=30)
gene2eid = {x["hugoGeneSymbol"]: x["entrezGeneId"] for x in r.json()}
log(f"[mut] entrez: {gene2eid}")

r = requests.post(f"{API}/molecular-profiles/{MUT_PROFILE}/mutations/fetch",
                  json={"sampleIds": thyroid_lines,
                        "entrezGeneIds": list(gene2eid.values())},
                  timeout=60)
muts = pd.DataFrame(r.json())
if len(muts):
    eid2gene = {v:k for k,v in gene2eid.items()}
    muts["gene"] = muts["entrezGeneId"].map(eid2gene)
log(f"[mut] hotspot mutations found: {len(muts)} rows across {muts['sampleId'].nunique() if len(muts) else 0} lines")

mut_truth_rows = []
for s in thyroid_lines:
    sub = muts[muts["sampleId"]==s] if len(muts) else pd.DataFrame()
    has_braf = bool(("BRAF" in sub["gene"].values) if len(sub) else False)
    has_ras  = bool(len(set(["NRAS","KRAS","HRAS"]) & set(sub["gene"].values)) if len(sub) else False)
    if has_braf and not has_ras:    label = "BRAF"
    elif has_ras and not has_braf:  label = "RAS"
    elif has_braf and has_ras:      label = "double_mutant"
    else:                           label = "other"
    proteins = ";".join(sub["proteinChange"].astype(str).tolist()) if len(sub) else ""
    mut_truth_rows.append(dict(sample_id=s, has_BRAF=has_braf,
                               has_RAS=has_ras, mutation_label=label,
                               protein_changes=proteins))
mut_truth = pd.DataFrame(mut_truth_rows)
mut_truth.to_csv(OUT / "ccle_thyroid_mutation_truth.tsv", sep="\t", index=False)
log(f"[mut] label counts: {mut_truth['mutation_label'].value_counts().to_dict()}")

# --------------------------------------------------------------------
# 3. BRS centroids on TCGA-THCA + apply to CCLE
# --------------------------------------------------------------------
brs = pd.read_csv(ROOT/"results"/"v12_literature"/"crosscheck"/
                  "reference_gene_lists"/"Chakravarty_2011_BRS52.tsv",
                  sep="\t")["gene"].tolist()
TCGA_TPM = Path("/data/thca/v5_cross_cancer/raw/THCA/tcga_log2tpm.tsv.gz")
TCGA_LBL = Path("/data/thca/v5_cross_cancer/raw/THCA/tcga_samples_labels.tsv")
tpm = pd.read_csv(TCGA_TPM, sep="\t", index_col=0)
labels = pd.read_csv(TCGA_LBL, sep="\t").set_index("sample_id")["label"]

brs_in_tcga = [g for g in brs if g in tpm.index]
samples_in = [s for s in tpm.columns if s in labels.index]
tpm_brs = tpm.loc[brs_in_tcga, samples_in]
y = labels.loc[samples_in]
c_braf = tpm_brs.loc[:, (y=="BRAF").values].mean(axis=1)
c_ras  = tpm_brs.loc[:, (y=="RAS").values].mean(axis=1)

expr_log = np.log2(expr.fillna(0) + 1.0)
common = expr_log.index.intersection(c_braf.index)
log(f"[brs-ccle] BRS genes overlapping CCLE: {len(common)}")

Xc = expr_log.loc[common].values - expr_log.loc[common].values.mean(axis=0, keepdims=True)
nx = np.sqrt((Xc**2).sum(axis=0))
cb_c = c_braf.loc[common].values - c_braf.loc[common].mean()
cr_c = c_ras.loc[common].values - c_ras.loc[common].mean()
nb = np.sqrt((cb_c**2).sum())
nr = np.sqrt((cr_c**2).sum())
corr_b = (Xc * cb_c[:, None]).sum(axis=0) / (nx*nb + 1e-12)
corr_r = (Xc * cr_c[:, None]).sum(axis=0) / (nx*nr + 1e-12)
score = corr_b - corr_r

brs_df = pd.DataFrame({
    "sample": expr_log.columns,
    "corr_braf": corr_b, "corr_ras": corr_r,
    "brs_score": score,
    "brs_label": np.where(score>0, "BRAF_like", "RAS_like"),
}).set_index("sample")
brs_df = brs_df.join(meta_df.set_index("sample_id")[
    ["name","cancer_type_detailed","oncotree","hist_subtype"]
], how="left")
brs_df = brs_df.join(mut_truth.set_index("sample_id")[
    ["mutation_label","protein_changes"]
], how="left")
brs_df.to_csv(OUT / "ccle_brs_labels.tsv", sep="\t")
log(f"[brs] BRS labels (n={len(brs_df)}): "
    f"{brs_df['brs_label'].value_counts().to_dict()}")

# --------------------------------------------------------------------
# 4. BRS vs mutation truth
# --------------------------------------------------------------------
both = brs_df.dropna(subset=["mutation_label"]).copy()
both = both[both["mutation_label"].isin(["BRAF","RAS"])]
log(f"[brs-vs-mut] lines with both BRS score and BRAF/RAS mutation: {len(both)}")
if len(both):
    confusion = pd.crosstab(both["mutation_label"], both["brs_label"])
    log(f"[brs-vs-mut] confusion:\n{confusion.to_string()}")
    agree = ((both["mutation_label"]=="BRAF") & (both["brs_label"]=="BRAF_like")).sum() + \
            ((both["mutation_label"]=="RAS")  & (both["brs_label"]=="RAS_like")).sum()
    log(f"[brs-vs-mut] BRS accuracy on CCLE: {agree}/{len(both)} "
        f"({100*agree/len(both):.1f}%)  (TCGA reference: 95.2%)")
    confusion.to_csv(OUT / "ccle_brs_vs_mutation.tsv", sep="\t")

# --------------------------------------------------------------------
# 5. Mutation-truth-based direction concordance
# --------------------------------------------------------------------
de = pd.read_csv(RES/"a_deseq2_raw"/"deseq2_tcga_raw.tsv", sep="\t")
de = de.dropna(subset=["gene_name"]).drop_duplicates("gene_name")
de_sig = de[de["significant"]].set_index("gene_name")

braf_mut = mut_truth[mut_truth["mutation_label"]=="BRAF"]["sample_id"].tolist()
ras_mut  = mut_truth[mut_truth["mutation_label"]=="RAS"]["sample_id"].tolist()
braf_mut = [s for s in braf_mut if s in expr_log.columns]
ras_mut  = [s for s in ras_mut  if s in expr_log.columns]
log(f"[concord-mut] BRAF-mut lines with expression: {len(braf_mut)} {braf_mut}")
log(f"[concord-mut] RAS-mut  lines with expression: {len(ras_mut)} {ras_mut}")

ccle_concord_pct = None
overlap = pd.Index([])
agree = 0
delta = pd.Series(dtype=float)
if len(braf_mut) >= 2 and len(ras_mut) >= 2:
    delta = (expr_log[braf_mut].mean(axis=1) - expr_log[ras_mut].mean(axis=1))
    overlap = de_sig.index.intersection(delta.index)
    log(f"[concord-mut] sig × CCLE overlap: {len(overlap)}")
    if len(overlap):
        tcga_sgn = (de_sig.loc[overlap, "log2FoldChange"] > 0).values
        ccle_sgn = (delta.loc[overlap] > 0).values
        agree = int((tcga_sgn == ccle_sgn).sum())
        ccle_concord_pct = 100 * agree / len(overlap)
        log(f"[concord-mut] BRAF-up direction concordance "
            f"(mutation-truth labels): {agree}/{len(overlap)} "
            f"({ccle_concord_pct:.1f}%)")
else:
    log(f"[concord-mut] insufficient mutation-truth lines per class")

# --------------------------------------------------------------------
# 6. 4-cohort match-count update with CCLE
# --------------------------------------------------------------------
con3 = pd.read_csv(RES/"e_brs_validation"/"three_cohort_per_gene_consensus.tsv",
                   sep="\t").set_index("gene")
log(f"[4cohort] 3-cohort consensus rows: {len(con3)}")

four_match_pct = {}
if ccle_concord_pct is not None and len(overlap):
    overlap_4 = con3.index.intersection(overlap)
    log(f"[4cohort] genes in all 4 cohorts: {len(overlap_4)}")
    if len(overlap_4):
        ccle_match = (((de_sig.loc[overlap_4, "log2FoldChange"] > 0).values ==
                       (delta.loc[overlap_4] > 0).values).astype(int))
        match4 = con3.loc[overlap_4, "match_count_3cohort"].values + ccle_match
        from collections import Counter
        cnt = Counter(match4.tolist())
        log(f"[4cohort] 4-cohort match-count distribution: {dict(sorted(cnt.items()))}")
        for k in [4, 3, 2, 1, 0]:
            n = (match4 == k).sum()
            four_match_pct[k] = 100 * n / len(match4)
            log(f"[4cohort] match {k}/4: {n} ({four_match_pct[k]:.1f}%)")
        log(f"[4cohort] genes matching ≥3 of 4: "
            f"{(match4>=3).sum()} ({100*(match4>=3).sum()/len(match4):.1f}%)")
        out4 = pd.DataFrame({
            "gene": list(overlap_4),
            "match_count_4cohort": match4,
        })
        out4.to_csv(OUT / "four_cohort_per_gene_consensus.tsv", sep="\t", index=False)

# --------------------------------------------------------------------
# 7. Final summary
# --------------------------------------------------------------------
summary = pd.DataFrame({
    "metric": [
        "ccle_thyroid_lines_total",
        "ccle_lines_with_expression",
        "ccle_BRAF_mutation",
        "ccle_RAS_mutation",
        "ccle_double_or_other",
        "BRS_accuracy_on_CCLE_vs_mutation",
        "deseq2_sig_in_ccle_universe",
        "BRAF_up_direction_concordance_TCGA_vs_CCLE_mutation_truth",
        "n_genes_in_all_4_cohorts",
        "consensus_4_of_4_pct",
        "consensus_>=3_of_4_pct",
    ],
    "value": [
        len(thyroid_lines),
        expr.shape[1],
        int((mut_truth["mutation_label"]=="BRAF").sum()),
        int((mut_truth["mutation_label"]=="RAS").sum()),
        int(mut_truth["mutation_label"].isin(["double_mutant","other"]).sum()),
        f"{100*agree/len(both):.1f}%" if len(both) else "n/a",
        len(overlap) if len(overlap) else 0,
        f"{ccle_concord_pct:.1f}%" if ccle_concord_pct is not None else "n/a",
        len(con3.index.intersection(overlap)) if len(overlap) else 0,
        f"{four_match_pct.get(4, 0):.1f}%" if four_match_pct else "n/a",
        f"{four_match_pct.get(4, 0) + four_match_pct.get(3, 0):.1f}%" if four_match_pct else "n/a",
    ],
})
summary.to_csv(OUT / "ccle_concordance.tsv", sep="\t", index=False)
log(f"[write] {OUT/'ccle_concordance.tsv'}")

log("v8.1 CCLE BRS+mutation-truth validation complete")
_log_f.close()
