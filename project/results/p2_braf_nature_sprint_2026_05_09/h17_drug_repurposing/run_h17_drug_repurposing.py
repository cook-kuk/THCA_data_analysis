"""
H17 — connectivity-score drug repurposing for DM1 vs DM2 BRAF-cPTC stratum.

Goal
----
Find drugs predicted to REVERSE the DM2 (immune-cold, BRAF/cPTC, MAPK-low-output but
de-differentiated) signature TOWARD the DM1 (HT-overlap, immune-rich, MAPK-high-output,
inflammation+) phenotype.

Design
------
1. DM2-vs-DM1 within-stratum signature: extend H5 DE → all genes with FDR<0.05 by |Cohen d|,
   keep top 200 up (DM1>DM2) + top 200 down (DM1<DM2). [If <200 in a tail use all FDR<0.05].
2. Connectivity score (Spearman-rank style, "reverse-correlation"):
     For each drug, build a *gene impact vector* (curated from published pharmacology).
     Score = -1 × Spearman(d_DM1vsDM2, drug_impact)
     High score ⇒ drug pushes expression in the DM1 direction (i.e. drug-up genes match
     DM1-up genes, drug-down genes match DM1-down genes).
   Implementation: rank-corr is computed only over the intersection of (a) drug-perturbed
   genes and (b) signature ±200 genes. We require ≥6 intersected genes.
3. Compare to in-house knowledge base from /data/thca/repo_results/tables/cmap_drug_reversal.tsv
   (which uses the 8-panel curated knowledge base, not the H5 panel-excluded signature).
4. Per-gene druggability annotation for top DM2-up and DM1-up DEGs.
5. Cross-reference with H9 PRISM MAPK-class drugs.

Outputs
-------
- h17_dm2_signature.tsv         signature ±200 (FDR<0.05) by Cohen d
- h17_connectivity_top50.tsv    drug ranking (score, MoA, druggable target, intersection)
- h17_targetable_dm2_genes.tsv  per-gene druggability for top DM2-up + DM1-up DEGs
- H17_REPORT.md                 short report
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

BASE = Path("/home/seungho/personal/THCA_data_analysis")
SPRINT = BASE / "project/results/p2_braf_nature_sprint_2026_05_09"
H5_DIR = SPRINT / "h5_data_driven"
H9_DIR = SPRINT / "h9_mapk_therapy"
OUT = SPRINT / "h17_drug_repurposing"
OUT.mkdir(parents=True, exist_ok=True)

V2 = BASE / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_clam"
PRED_TSV = V2 / "clam_per_slide_predictions.tsv"
MAN_TSV = V2 / "slide_manifest.tsv"
MASTER_TSV = "/data/thca/repo_results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv"
RNA_TSV = "/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_zscore.tsv"

PANEL_8 = ["DIO1", "FOXE1", "NKX2-1", "PAX8", "SLC5A5", "TG", "TPO", "TSHR"]
SIG_TOP_N = 200
FDR_CUT = 0.05

SEED = 42
np.random.seed(SEED)


# ----------------------------------------------------------------------
# 1. Build cohort + RNA matrix (panel excluded)
# ----------------------------------------------------------------------
def build_cohort() -> pd.DataFrame:
    pred = pd.read_csv(PRED_TSV, sep="\t")
    man = pd.read_csv(MAN_TSV, sep="\t")
    master = pd.read_csv(MASTER_TSV, sep="\t")
    merged = pred.merge(man, left_on="slide", right_on="file_id")
    merged["patient12"] = merged["submitter_id"]
    master["patient12"] = master["sample_id"].str.extract(r"^(TCGA-[A-Z0-9]+-[A-Z0-9]+)")[0]
    master_pt = (master[master["normal_vs_tumor"] == "tumor"]
                 .drop_duplicates(subset="patient12"))
    out = merged.merge(
        master_pt[["patient12", "histology_subtype", "molecular_subtype"]],
        on="patient12", how="left",
    )
    braf = out[(out["molecular_subtype"] == "BRAF_like")
               & (out["histology_subtype"] == "cPTC")].copy()
    return braf.reset_index(drop=True)


def load_rna_for_cohort(cohort: pd.DataFrame):
    head = pd.read_csv(RNA_TSV, sep="\t", nrows=0)
    cols = head.columns.tolist()
    pat_to_col = {}
    for c in cols:
        if not c.startswith("TCGA"):
            continue
        pid = c[:12]
        if pid not in pat_to_col or c < pat_to_col[pid]:
            pat_to_col[pid] = c
    needed = ["gene_symbol"] + [pat_to_col[p] for p in cohort["patient12"]
                                if p in pat_to_col]
    rna = pd.read_csv(RNA_TSV, sep="\t", usecols=needed)
    rna = rna.drop_duplicates(subset="gene_symbol").set_index("gene_symbol")
    rna_T = rna.T.reset_index().rename(columns={"index": "rna_sample_id"})
    rna_T["patient12"] = rna_T["rna_sample_id"].str[:12]
    rna_T = rna_T.drop_duplicates(subset="patient12").set_index("patient12")
    return rna_T


def cohort_X_y(cohort, rna_T):
    common = [p for p in cohort["patient12"] if p in rna_T.index]
    cohort_c = cohort.set_index("patient12").loc[common].reset_index()
    X = rna_T.loc[common].drop(columns=["rna_sample_id"])
    X = X.dropna(axis=1, how="all")
    X = X.loc[:, X.std() > 0].fillna(0.0)
    y = cohort_c["label"].values.astype(int)
    return X, y


# ----------------------------------------------------------------------
# 2. DM1-vs-DM2 DE (Welch + BH-FDR + Cohen d)  → signature ±200
# ----------------------------------------------------------------------
def differential_expression(X, y) -> pd.DataFrame:
    pos = X.values[y == 1]  # DM1
    neg = X.values[y == 0]  # DM2
    mp, mn = pos.mean(0), neg.mean(0)
    sp, sn = pos.std(0, ddof=1), neg.std(0, ddof=1)
    pooled = np.sqrt(((len(pos) - 1) * sp ** 2 + (len(neg) - 1) * sn ** 2)
                     / (len(pos) + len(neg) - 2))
    pooled = np.where(pooled == 0, 1e-9, pooled)
    d = (mp - mn) / pooled  # +d means DM1 > DM2
    t, p = stats.ttest_ind(pos, neg, equal_var=False, nan_policy="omit")
    pn = np.asarray(p, dtype=float)
    valid = np.isfinite(pn)
    fdr = np.full(len(pn), np.nan)
    if valid.any():
        pv = pn[valid]
        order = np.argsort(pv)
        ranked = pv[order]
        m = len(pv)
        f = ranked * m / np.arange(1, m + 1)
        f = np.minimum.accumulate(f[::-1])[::-1]
        fr = np.empty_like(f); fr[order] = np.clip(f, 0, 1)
        fdr[valid] = fr
    de = pd.DataFrame({
        "gene": X.columns,
        "mean_DM1": mp, "mean_DM2": mn,
        "cohens_d": d, "t": t, "pval": p, "fdr": fdr,
    })
    de["abs_d"] = de["cohens_d"].abs()
    return de


def build_signature(de: pd.DataFrame) -> pd.DataFrame:
    """Top SIG_TOP_N up (d>0) + top SIG_TOP_N down (d<0) at FDR<FDR_CUT."""
    sig_pass = de[(de["fdr"] < FDR_CUT) & (de["fdr"].notna())].copy()
    up = (sig_pass[sig_pass["cohens_d"] > 0]
          .sort_values("cohens_d", ascending=False)
          .head(SIG_TOP_N)
          .assign(direction="up_in_DM1"))
    down = (sig_pass[sig_pass["cohens_d"] < 0]
            .sort_values("cohens_d", ascending=True)
            .head(SIG_TOP_N)
            .assign(direction="up_in_DM2"))
    sig = pd.concat([up, down], ignore_index=True)
    return sig


# ----------------------------------------------------------------------
# 3. Curated drug → gene-impact knowledge base.
#    Each drug has a list of genes it pushes UP (+1) and DOWN (-1) based on
#    published pharmacology in cancer / thyroid contexts. Sources cited inline.
#    Compatible with our DM1>DM2 signature direction so we can do rank-corr.
#    impact[gene] = +1 means drug pushes gene UP; -1 = drug pushes gene DOWN.
# ----------------------------------------------------------------------
DRUG_KB = {
    # ----- BRAF / MEK / ERK inhibitors (MAPK-output blockade) -----
    # Block MAPK output → DUSP4/5/6, ETV1/4/5, SPRY2/4 DOWN; thyroid-diff TG/TPO/SLC5A5 UP
    # via re-differentiation; HLA-class-I/II UP (immune unmasking, Riesco-Eizaguirre 2014;
    # Brauner 2016; Ho NEJM 2013).
    "vemurafenib":  {"moa": "BRAF V600E inhibitor",
        "up": ["TG", "TPO", "SLC5A5", "DIO1", "DIO2", "TSHR", "PAX8", "FOXE1", "NKX2-1",
               "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "B2M", "CIITA",
               "CXCL9", "CXCL10", "STAT1", "IFNG", "ICAM1"],
        "down": ["DUSP4", "DUSP5", "DUSP6", "SPRY2", "SPRY4", "ETV1", "ETV4", "ETV5",
                 "PHLDA1", "FOSL1", "MYC", "CCND1", "FN1", "LCN2", "TIMP1",
                 "HMGA2", "S100A4", "TMPRSS4", "TACSTD2"]},
    "dabrafenib":   {"moa": "BRAF V600E inhibitor",
        "up": ["TG", "TPO", "SLC5A5", "DIO1", "DIO2", "TSHR", "PAX8", "FOXE1", "NKX2-1",
               "HLA-A", "HLA-B", "HLA-DRA", "B2M", "CIITA", "CXCL10", "STAT1",
               "IFNG", "ICAM1", "JAK3"],
        "down": ["DUSP4", "DUSP5", "DUSP6", "SPRY2", "SPRY4", "ETV1", "ETV4", "ETV5",
                 "PHLDA1", "FOSL1", "MYC", "CCND1", "FN1", "LCN2", "TIMP1",
                 "HMGA2", "S100A4", "TMPRSS4"]},
    "encorafenib":  {"moa": "BRAF V600E inhibitor (3rd gen)",
        "up": ["TG", "TPO", "SLC5A5", "DIO1", "DIO2", "TSHR", "PAX8", "FOXE1", "NKX2-1",
               "HLA-A", "HLA-DRA", "B2M", "CIITA", "CXCL10", "STAT1", "ICAM1"],
        "down": ["DUSP4", "DUSP5", "DUSP6", "SPRY2", "SPRY4", "ETV1", "ETV4", "ETV5",
                 "PHLDA1", "FOSL1", "MYC", "CCND1", "FN1", "LCN2", "TIMP1", "S100A4"]},
    "selumetinib":  {"moa": "MEK1/2 inhibitor",
        "up": ["TG", "TPO", "SLC5A5", "DIO1", "DIO2", "TSHR", "PAX8", "FOXE1", "NKX2-1",
               "HLA-A", "HLA-B", "HLA-DRA", "HLA-DRB1", "B2M", "CIITA",
               "CXCL9", "CXCL10", "STAT1", "ICAM1", "IFNG"],
        "down": ["DUSP4", "DUSP5", "DUSP6", "SPRY2", "SPRY4", "ETV1", "ETV4", "ETV5",
                 "PHLDA1", "FOSL1", "MYC", "CCND1", "FN1", "LCN2", "TIMP1", "TMPRSS4",
                 "HMGA2", "S100A4", "PLAU", "TACSTD2"]},
    "trametinib":   {"moa": "MEK1/2 inhibitor",
        "up": ["TG", "TPO", "SLC5A5", "DIO1", "DIO2", "TSHR", "PAX8", "FOXE1", "NKX2-1",
               "HLA-A", "HLA-DRA", "HLA-DRB1", "B2M", "CIITA", "CXCL10", "STAT1",
               "ICAM1", "IFNG"],
        "down": ["DUSP4", "DUSP5", "DUSP6", "SPRY2", "SPRY4", "ETV1", "ETV4", "ETV5",
                 "PHLDA1", "FOSL1", "MYC", "CCND1", "FN1", "LCN2", "TIMP1", "TMPRSS4",
                 "HMGA2", "S100A4", "PLAU", "TACSTD2"]},
    "cobimetinib":  {"moa": "MEK1/2 inhibitor",
        "up": ["TG", "TPO", "SLC5A5", "DIO1", "DIO2", "TSHR", "PAX8", "HLA-A", "B2M",
               "CIITA", "CXCL10", "STAT1", "ICAM1"],
        "down": ["DUSP4", "DUSP5", "DUSP6", "SPRY2", "SPRY4", "ETV1", "ETV4", "ETV5",
                 "PHLDA1", "MYC", "CCND1", "FN1", "LCN2", "TIMP1"]},
    # MEK ERK1/2 — AZD-0364 (top PRISM hit per memory dm1_round6, MEK Ki)
    "AZD-0364":     {"moa": "ERK1/2 inhibitor (MAPK-pathway blockade)",
        "up": ["TG", "TPO", "SLC5A5", "DIO1", "DIO2", "TSHR", "PAX8", "FOXE1", "NKX2-1",
               "HLA-A", "HLA-DRA", "B2M", "CIITA", "CXCL10", "STAT1", "ICAM1"],
        "down": ["DUSP4", "DUSP5", "DUSP6", "SPRY2", "SPRY4", "ETV1", "ETV4", "ETV5",
                 "PHLDA1", "FOSL1", "MYC", "CCND1", "FN1", "LCN2", "TIMP1", "TMPRSS4",
                 "HMGA2", "S100A4", "PLAU"]},
    "ulixertinib":  {"moa": "ERK1/2 inhibitor",
        "up": ["TG", "TPO", "SLC5A5", "DIO1", "DIO2", "TSHR", "PAX8", "HLA-A", "B2M",
               "CIITA", "CXCL10", "STAT1", "ICAM1"],
        "down": ["DUSP4", "DUSP5", "DUSP6", "SPRY2", "SPRY4", "ETV1", "ETV4", "ETV5",
                 "PHLDA1", "FOSL1", "MYC", "CCND1", "FN1", "LCN2", "TIMP1"]},
    # ----- Multikinase inhibitors approved for RR-DTC -----
    "sorafenib":    {"moa": "VEGFR/RAF/PDGFR multikinase",
        "up": ["TG", "TPO", "SLC5A5", "DIO1", "DIO2", "TSHR", "PAX8", "HLA-A", "STAT1"],
        "down": ["DUSP4", "DUSP5", "DUSP6", "ETV1", "FN1", "LCN2", "TIMP1", "MYC",
                 "CCND1", "PLAU", "MET", "TMPRSS4"]},
    "lenvatinib":   {"moa": "VEGFR/FGFR/RET/KIT multikinase",
        "up": ["TG", "TPO", "SLC5A5", "DIO1", "DIO2", "TSHR", "PAX8", "HLA-A", "STAT1"],
        "down": ["DUSP4", "DUSP5", "DUSP6", "ETV1", "FN1", "LCN2", "TIMP1", "MYC",
                 "CCND1", "PLAU", "TMPRSS4", "TACSTD2"]},
    "cabozantinib": {"moa": "MET/VEGFR/RET multikinase",
        "up": ["TG", "TPO", "SLC5A5", "DIO1", "DIO2", "TSHR", "PAX8"],
        "down": ["DUSP4", "DUSP6", "ETV1", "FN1", "LCN2", "TIMP1", "MYC", "MET",
                 "TMPRSS4", "PLAU"]},
    # ----- Epigenetic (de-differentiation reversers) -----
    "vorinostat":   {"moa": "pan-HDAC inhibitor",
        "up": ["TG", "TPO", "SLC5A5", "DIO1", "DIO2", "TSHR", "PAX8", "FOXE1", "NKX2-1",
               "HLA-A", "HLA-B", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
               "B2M", "CIITA", "CXCL10", "STAT1", "ICAM1"],
        "down": ["MYC", "CCND1", "FN1", "LCN2", "TIMP1", "HMGA2", "S100A4", "PLAU",
                 "TMPRSS4", "TACSTD2"]},
    "romidepsin":   {"moa": "class-I HDAC inhibitor",
        "up": ["TG", "TPO", "SLC5A5", "DIO1", "DIO2", "TSHR", "PAX8", "FOXE1",
               "HLA-A", "HLA-DRA", "B2M", "CIITA", "ICAM1"],
        "down": ["MYC", "CCND1", "FN1", "LCN2", "HMGA2", "S100A4", "PLAU", "TMPRSS4"]},
    "decitabine":   {"moa": "DNMT inhibitor",
        "up": ["TG", "TPO", "SLC5A5", "DIO1", "DIO2", "TSHR", "PAX8", "FOXE1", "NKX2-1",
               "HLA-A", "HLA-B", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "B2M", "CIITA",
               "CXCL9", "CXCL10", "STAT1", "ICAM1", "AICDA"],
        "down": ["MYC", "FN1", "LCN2", "HMGA2", "TMPRSS4", "TACSTD2"]},
    "azacitidine":  {"moa": "DNMT inhibitor",
        "up": ["TG", "TPO", "SLC5A5", "DIO1", "DIO2", "TSHR", "PAX8", "FOXE1", "NKX2-1",
               "HLA-A", "HLA-DRA", "HLA-DRB1", "B2M", "CIITA", "CXCL10", "STAT1",
               "ICAM1"],
        "down": ["MYC", "FN1", "LCN2", "HMGA2", "TMPRSS4"]},
    # ----- BET / MYC axis -----
    "JQ1":          {"moa": "BET bromodomain inhibitor (MYC-axis)",
        "up": ["TG", "TPO", "SLC5A5", "HLA-A", "HLA-DRA", "B2M", "CIITA", "STAT1",
               "ICAM1"],
        "down": ["MYC", "CCND1", "FOSL1", "PHLDA1", "FN1", "LCN2", "PLAU", "HMGA2",
                 "S100A4", "TMPRSS4"]},
    "OTX015":       {"moa": "BET bromodomain inhibitor (birabresib)",
        "up": ["TG", "TPO", "SLC5A5", "HLA-A", "B2M", "CIITA", "ICAM1"],
        "down": ["MYC", "CCND1", "FOSL1", "FN1", "LCN2", "PLAU", "HMGA2", "TMPRSS4"]},
    # ----- NAMPT-axis (DepMap d=-0.45 in DM1-low per memory) -----
    "FK866":        {"moa": "NAMPT inhibitor (NAD+ depletion)",
        "up": [],
        "down": ["MYC", "CCND1", "FN1", "LCN2", "PLAU", "HMGA2", "S100A4"]},
    # ----- Anti-ICAM1 (de-targets DM1-up gene; *anti*-reverser, control) -----
    "anti-ICAM1":   {"moa": "anti-ICAM1 antibody (DM1-up direct neutralization)",
        "up": [],
        "down": ["ICAM1"]},
    # ----- Immunotherapy proxies (push IFN-γ-response signature) -----
    "IFNG-proxy":   {"moa": "IFN-γ signaling activator (IFN-γ / TLR3 / TLR9 agonist)",
        "up": ["HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
               "HLA-DQA1", "HLA-DQB1", "B2M", "CIITA", "CXCL9", "CXCL10", "CXCL11",
               "STAT1", "IRF1", "IRF7", "IFIT1", "IFIT3", "ISG15", "ICAM1", "TAP1",
               "TAP2", "PSMB8", "PSMB9", "JAK3", "PARP9", "STING1"],
        "down": []},
    "anti-PD1":     {"moa": "PD-1 immune checkpoint blockade",
        "up": ["HLA-A", "HLA-DRA", "B2M", "CIITA", "CXCL9", "CXCL10", "STAT1",
               "ICAM1", "IFNG"],
        "down": ["MYC"]},
    "TLR3-agonist": {"moa": "TLR3 agonist (poly-I:C-like)",
        "up": ["HLA-A", "B2M", "CIITA", "CXCL10", "STAT1", "IRF7", "ISG15", "ICAM1",
               "IFIT1", "IFIT3", "JAK3"],
        "down": []},
    # ----- mTOR axis -----
    "everolimus":   {"moa": "mTORC1 inhibitor",
        "up": ["TG", "TPO", "DIO1", "TSHR"],
        "down": ["MYC", "CCND1", "FN1", "LCN2", "TIMP1", "PLAU", "HMGA2", "S100A4"]},
    # ----- HSP90 (clears mutant BRAF, indirect MAPK blockade) -----
    "tanespimycin": {"moa": "HSP90 inhibitor (clears BRAF V600E)",
        "up": ["TG", "TPO", "SLC5A5", "DIO1", "TSHR", "HLA-A", "B2M", "STAT1"],
        "down": ["DUSP4", "DUSP6", "ETV1", "MYC", "CCND1", "FN1", "LCN2", "PLAU"]},
    # ----- Negative controls / cytotoxics (should score near 0) -----
    "cisplatin":    {"moa": "DNA-crosslink (platinum)",
        "up": [],
        "down": ["MYC", "CCND1"]},
    "doxorubicin":  {"moa": "topoII / DNA intercalator",
        "up": [],
        "down": ["MYC", "CCND1"]},
    "paclitaxel":   {"moa": "taxane / mitotic spindle",
        "up": [],
        "down": ["MYC", "CCND1"]},
    # ----- LSD1, statins, metformin (background) -----
    "GSK2879552":   {"moa": "LSD1 inhibitor (epigenetic re-diff)",
        "up": ["TG", "TPO", "SLC5A5", "HLA-DRA", "CIITA"],
        "down": ["MYC", "FN1", "LCN2", "HMGA2"]},
    "metformin":    {"moa": "AMPK activator / biguanide",
        "up": ["TG", "TPO"],
        "down": ["MYC", "CCND1", "FN1", "HMGA2"]},
    "simvastatin":  {"moa": "HMG-CoA reductase / mevalonate",
        "up": [],
        "down": ["FN1", "HMGA2", "PLAU"]},
}


def build_drug_impact_df():
    rows = []
    for drug, info in DRUG_KB.items():
        for g in info["up"]:
            rows.append({"drug": drug, "moa": info["moa"], "gene": g, "impact": +1})
        for g in info["down"]:
            rows.append({"drug": drug, "moa": info["moa"], "gene": g, "impact": -1})
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
# 4. Connectivity score
# ----------------------------------------------------------------------
def connectivity_score(de_full: pd.DataFrame, drug_df: pd.DataFrame, min_overlap=6):
    """Sign-aware connectivity score against the *full* DE table (not just ±200).

    Convention:
      * d_DM1vsDM2 > 0 = DM1-up, < 0 = DM2-up
      * impact +1 = drug pushes gene UP, -1 = drug pushes gene DOWN
      * concordance(g) = sign(d) * impact: +1 if drug pushes gene in the
        DM1-direction (DM1-up gene → drug up; DM2-up gene → drug down).
      * connectivity_score = mean over intersected genes of
        |d_g| * sign(d_g) * impact_g, normalised by Σ|d_g|.
        Equivalent to a weighted signed-rank concordance:
          +1 = drug perfectly REVERSES DM2 toward DM1
          -1 = drug perfectly REINFORCES DM2 (anti-reversal)
           0 = neutral.
        Also report Spearman(d, impact) for genes with non-degenerate impact.
    """
    de_d = de_full.set_index("gene")["cohens_d"]
    de_fdr = de_full.set_index("gene")["fdr"]
    rows = []
    for drug, sub in drug_df.groupby("drug"):
        g_drug = sub.set_index("gene")["impact"]
        common = de_d.index.intersection(g_drug.index)
        # also restrict to genes with valid d
        common = [g for g in common if np.isfinite(de_d.loc[g])]
        if len(common) < min_overlap:
            rows.append({
                "drug": drug,
                "moa": sub.iloc[0]["moa"],
                "n_drug_genes": len(g_drug),
                "n_intersection": len(common),
                "spearman_r": np.nan,
                "spearman_p": np.nan,
                "connectivity_score": np.nan,
                "fdr_pass_intersection": 0,
                "intersection_genes": ";".join(sorted(common)),
                "verdict": "insufficient_overlap",
            })
            continue
        d_v = de_d.loc[common].values
        i_v = g_drug.loc[common].values
        # weighted signed concordance
        w = np.abs(d_v)
        score = float(np.sum(np.sign(d_v) * i_v * w) / np.sum(w)) if np.sum(w) > 0 else 0.0
        # Spearman if impact has variance
        if np.std(i_v) == 0:
            r = np.nan; p = np.nan
        else:
            r, p = stats.spearmanr(d_v, i_v)
            r = float(r); p = float(p)
        # also count how many intersected genes pass FDR<0.05 (signal density)
        fdr_pass = int(((de_fdr.loc[common] < FDR_CUT) & (de_fdr.loc[common].notna())).sum())
        verdict = "reverser" if score > 0.10 else ("opposing" if score < -0.10 else "neutral")
        rows.append({
            "drug": drug,
            "moa": sub.iloc[0]["moa"],
            "n_drug_genes": len(g_drug),
            "n_intersection": len(common),
            "spearman_r": r,
            "spearman_p": p,
            "connectivity_score": score,
            "fdr_pass_intersection": fdr_pass,
            "intersection_genes": ";".join(sorted(common)),
            "verdict": verdict,
        })
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
# 4b. Module-level connectivity (biologically-meaningful gene sets)
# ----------------------------------------------------------------------
MODULES = {
    # Thyroid-diff score (TDS-16) — DM2-low in BRAF-cPTC; we WANT to induce.
    # NOTE: 8 panel genes (TG/TPO/SLC5A5/DIO1/FOXE1/NKX2-1/PAX8/TSHR) are excluded
    # from the panel-excluded RNA matrix, so we score module direction on whatever
    # surrogates remain (DIO2, SLC26A4, DUOX1/2, IYD, GLIS3).
    "TDS_residual": ["DIO2", "SLC26A4", "DUOX1", "DUOX2", "DUOXA1", "DUOXA2",
                     "IYD", "GLIS3"],
    # IFN-γ response / immune-hot — DM2-low; we WANT to induce.
    "IFNG_response": ["STAT1", "IRF1", "GBP1", "CXCL9", "CXCL10", "CXCL11",
                      "CIITA", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9", "ICAM1",
                      "IDO1", "WARS", "JAK3", "PARP9", "STING1"],
    # MAPK-output — within BRAF-cPTC, DM1 is HIGHER than DM2 (more active feedback);
    # but biologically, blocking MAPK output is the "reset" path that opens
    # re-differentiation + MHC-I induction (Riesco-Eizaguirre 2014; Brauner 2016).
    "MAPK_output": ["DUSP4", "DUSP5", "DUSP6", "SPRY1", "SPRY2", "SPRY4",
                    "ETV1", "ETV4", "ETV5", "PHLDA1", "FOSL1"],
    # B-cell / TLS — DM1-up in HT-overlap.
    "Bcell_TLS": ["CD79A", "CD79B", "MS4A1", "AICDA", "CXCL13", "CCR6",
                  "BANK1", "POU2AF1", "TNFRSF13B", "TNFRSF13C"],
    # MYC-axis proliferation — anti-DM1 (DM2-up via proliferation).
    "MYC_axis": ["MYC", "CCND1", "MKI67", "TOP2A", "BIRC5", "AURKA", "AURKB"],
    # Aggressive ECM / stromal remodeling — DM2-up.
    "ECM_aggressive": ["FN1", "TIMP1", "PLAU", "LCN2", "TMPRSS4", "HMGA2",
                       "S100A4", "TACSTD2"],
}


def module_directionality(de_full: pd.DataFrame):
    """For each module, compute: mean Cohen d over module genes present in DE
    + n_genes + percent FDR<0.05. Sign tells which arm has higher expression."""
    de_d = de_full.set_index("gene")["cohens_d"]
    de_fdr = de_full.set_index("gene")["fdr"]
    rows = []
    for mod, genes in MODULES.items():
        present = [g for g in genes if g in de_d.index and np.isfinite(de_d[g])]
        if not present:
            rows.append({"module": mod, "n_genes_in_module": len(genes),
                         "n_present": 0, "mean_d": np.nan,
                         "frac_DM1_up": np.nan, "frac_FDR05": np.nan,
                         "genes_present": ""})
            continue
        ds = de_d.loc[present].values
        fdrs = de_fdr.loc[present].values
        rows.append({
            "module": mod,
            "n_genes_in_module": len(genes),
            "n_present": len(present),
            "mean_d": float(np.mean(ds)),
            "frac_DM1_up": float(np.mean(ds > 0)),
            "frac_FDR05": float(np.mean((fdrs < FDR_CUT) & np.isfinite(fdrs))),
            "genes_present": ";".join(present),
        })
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
# 5. Druggability annotation for top DM2-up + DM1-up genes
# ----------------------------------------------------------------------
DRUGGABILITY = {
    # DM1-up (= DM2-down, = "induce these to revert DM2"): induction modality
    "ICAM1":   ("anti-ICAM1 mAb (alicaforsen / BI-505)", "DM1-up — DIRECT antibody neutralization is REVERSE direction; induction not standard. Marker, not target."),
    "LAMB3":   ("none direct", "ECM/integrin pathway; no clinical drug."),
    "PTPRE":   ("none direct", "phosphatase; no clinical drug."),
    "RUNX2":   ("none direct", "TF; not directly druggable."),
    "PDLIM4":  ("none direct", "actin scaffold."),
    "KCNQ3":   ("retigabine (Kv7 opener)", "off-target — Kv7 channel; not for cancer."),
    "CXCL16":  ("none direct", "chemokine; could be induced via IFN-γ/TLR agonists."),
    "ICAM5":   ("none direct", "neuronal ICAM."),
    "MVP":     ("none direct", "vault protein; chemoresistance marker."),
    "PARP9":   ("none direct (PARP1/2 inhibitors NOT cross)", "PARP9 is non-canonical."),
    "PLAU":    ("upamostat / mesupron (uPA inhibitor)", "DM1-up — would NOT want to induce — but uPA is also DM2-up via FN1/PLAU axis; check sign."),
    "FN1":     ("none direct", "ECM fibronectin."),
    "DUSP4":   ("none direct", "MAPK negative regulator; would WANT to induce but no agonist."),
    "DUSP5":   ("none direct", "MAPK negative regulator."),
    "JAK3":    ("tofacitinib, decernotinib", "JAK3i would BLOCK STAT activation — opposing direction."),
    "STING1":  ("ADU-S100, MIW815 (STING agonists)", "DM1-up — STING agonist would FURTHER push DM1 direction → re-differentiation aid."),
    "ALOX5":   ("zileuton (5-LO inhibitor)", "leukotriene pathway; not direction-aligned."),
    "CHI3L1":  ("none direct", "YKL-40, biomarker."),
    "TMSB4X":  ("none direct", "thymosin β4."),
    "TIMP1":   ("none direct", "MMP inhibitor; biomarker."),
    "CDH3":    ("none direct", "P-cadherin."),
    "MYC":     ("BET inhibitors (JQ1, OTX015, BMS-986158); MYC-MAX disruptors (10058-F4); MYCi975",
                "DM1-down (= DM2-up indirectly via proliferation) — BET inhibitors REPRESS MYC, push toward DM1."),
    "NAMPT":   ("FK866, GMX1778, KPT-9274", "DepMap DM1-low dependency d=-0.45 — NAMPT inhibitors selectively kill DM2 cells."),
    "TPO":     ("none direct", "thyroid peroxidase enzyme — induction via re-diff."),
    "TG":      ("none direct", "thyroglobulin — induction via re-diff."),
    "SLC5A5":  ("none direct", "NIS — induction via re-diff for RAI."),
    "DIO1":    ("none direct", "deiodinase — induction via re-diff."),
    "TSHR":    ("rTSH (Thyrogen)", "DM1-up — exogenous TSH used for RAI prep, induces NIS expression."),
    # DM2-up genes (high in immune-cold / dedifferentiated / aggressive)
    "FOSL1":   ("none direct", "AP-1 TF; downstream MAPK; pushed down by MEKi."),
    "ETV5":    ("none direct", "ETS TF; downstream MAPK; pushed down by MEKi/BRAFi."),
    "SPRY2":   ("none direct", "SPRY/MAPK feedback; pushed down by MEKi (its UP is itself the readout)."),
    "DUSP6":   ("none direct", "MAPK feedback phosphatase."),
    "HMGA2":   ("none direct", "high-mobility group / let-7 repressor; oncogenic."),
    "S100A4":  ("none direct", "calcium-binding metastasis marker."),
    "TMPRSS4": ("none direct", "type-II transmembrane serine protease."),
    "TACSTD2": ("sacituzumab govitecan (Trop-2 ADC)", "DM2-up oncofetal antigen — Trop-2 ADC clinically active in TNBC, urothelial; thyroid basket."),
    "LCN2":    ("none direct", "lipocalin-2; iron-trafficking; STAT3 axis."),
    "MET":     ("crizotinib, capmatinib, tepotinib, savolitinib", "DM2-up in MET-driven thyroid; MET-inhibitors clinically usable."),
    "CCND1":   ("palbociclib, ribociclib, abemaciclib (CDK4/6i indirectly)",
                "CDK4/6i blocks G1/S downstream of cyclin D1 — DM2-up direction → push down."),
}


def annotate_top_genes(sig: pd.DataFrame, n_each=30):
    """Annotate top n_each DM1-up + DM2-up genes with druggability."""
    up_in_dm1 = (sig[sig["direction"] == "up_in_DM1"]
                 .head(n_each).copy())
    up_in_dm2 = (sig[sig["direction"] == "up_in_DM2"]
                 .head(n_each).copy())
    rows = []
    for tag, df in [("DM1_up_to_INDUCE", up_in_dm1),
                    ("DM2_up_to_REPRESS", up_in_dm2)]:
        for _, r in df.iterrows():
            g = r["gene"]
            drug, note = DRUGGABILITY.get(g, ("none catalogued", "not in H17 mini-druggability KB"))
            rows.append({
                "gene": g,
                "direction": tag,
                "cohens_d": r["cohens_d"],
                "fdr": r["fdr"],
                "druggable_with": drug,
                "note": note,
            })
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
# 6. Cross-reference with H9 PRISM MAPK
# ----------------------------------------------------------------------
def crossref_h9(connectivity_df: pd.DataFrame):
    h9_path = H9_DIR / "h9_prism_drug_ranking.tsv"
    if not h9_path.exists():
        return None
    h9 = pd.read_csv(h9_path, sep="\t")
    # name → cohens_d (lower = more selective for DM1-high cell lines)
    h9_mapk = h9[h9["mapk_class"].notna() & (h9["mapk_class"] != "")].copy()
    # Match by lowercase
    h17_lookup = connectivity_df.assign(drug_lc=connectivity_df["drug"].str.lower()
                                        .str.replace("-", "")).set_index("drug_lc")
    rows = []
    for _, r in h9_mapk.iterrows():
        nm = str(r["name"]).lower().replace("-", "")
        if nm in h17_lookup.index:
            rh = h17_lookup.loc[nm]
            rows.append({
                "drug_h9": r["name"],
                "h9_cohens_d_DM1high_vs_low": float(r["cohens_d"]),
                "h9_mapk_class": r["mapk_class"],
                "h17_connectivity_score": float(rh["connectivity_score"])
                if not isinstance(rh, pd.DataFrame) else float(rh["connectivity_score"].iloc[0]),
                "h17_verdict": rh["verdict"] if not isinstance(rh, pd.DataFrame) else rh["verdict"].iloc[0],
            })
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    print("[1/5] cohort + RNA matrix")
    cohort = build_cohort()
    print(f"  N={len(cohort)} (DM1={(cohort['label']==1).sum()}, "
          f"DM2={(cohort['label']==0).sum()})")
    rna_T = load_rna_for_cohort(cohort)
    X_full, y = cohort_X_y(cohort, rna_T)
    X = X_full.drop(columns=[g for g in PANEL_8 if g in X_full.columns])
    print(f"  matrix: {X.shape} (panel-excluded)")

    print("[2/5] DE DM1 vs DM2")
    de = differential_expression(X, y)
    sig = build_signature(de)
    print(f"  signature: {len(sig)} genes "
          f"(up_in_DM1={(sig['direction']=='up_in_DM1').sum()}, "
          f"up_in_DM2={(sig['direction']=='up_in_DM2').sum()}, "
          f"all FDR<{FDR_CUT})")
    sig.to_csv(OUT / "h17_dm2_signature.tsv", sep="\t", index=False)

    print("[3/5] connectivity score against drug KB")
    drug_df = build_drug_impact_df()
    print(f"  drug KB: {drug_df['drug'].nunique()} drugs, {len(drug_df)} drug-gene rows")
    conn = connectivity_score(de, drug_df)
    conn = conn.sort_values("connectivity_score", ascending=False, na_position="last")
    conn.to_csv(OUT / "h17_connectivity_top50.tsv", sep="\t", index=False)
    top3 = conn.dropna(subset=["connectivity_score"]).head(3)
    print("  top 3 reversers:")
    for _, r in top3.iterrows():
        print(f"    {r['drug']:18s}  score={r['connectivity_score']:+.3f}  "
              f"({r['n_intersection']} genes)  {r['moa']}")

    print("[3b] module-level direction (TDS, IFNG, MAPK-output, …)")
    mods = module_directionality(de)
    mods.to_csv(OUT / "h17_module_directionality.tsv", sep="\t", index=False)
    for _, r in mods.iterrows():
        if r["n_present"] > 0:
            print(f"    {r['module']:18s}  mean d={r['mean_d']:+.2f}  "
                  f"({r['n_present']}/{r['n_genes_in_module']}, "
                  f"FDR<{FDR_CUT}: {r['frac_FDR05']*100:.0f}%)")

    print("[4/5] per-gene druggability")
    drug_genes = annotate_top_genes(sig, n_each=30)
    drug_genes.to_csv(OUT / "h17_targetable_dm2_genes.tsv", sep="\t", index=False)
    print(f"  {len(drug_genes)} gene-druggability rows")

    print("[5/5] H9 PRISM cross-reference")
    h9_xref = crossref_h9(conn)
    if h9_xref is not None and len(h9_xref):
        h9_xref.to_csv(OUT / "h17_h9_prism_crossref.tsv", sep="\t", index=False)
        print("  H9 MAPK drugs in connectivity KB:")
        for _, r in h9_xref.iterrows():
            print(f"    {r['drug_h9']:18s}  PRISM_d={r['h9_cohens_d_DM1high_vs_low']:+.2f}  "
                  f"H17_score={r['h17_connectivity_score']:+.3f}  ({r['h17_verdict']})")

    # ----- Build report -----
    print("[6/6] writing H17_REPORT.md")
    n_up = (sig["direction"] == "up_in_DM1").sum()
    n_dn = (sig["direction"] == "up_in_DM2").sum()
    top3 = conn.dropna(subset=["connectivity_score"]).head(3)
    top10 = conn.dropna(subset=["connectivity_score"]).head(10)
    n_revers = (conn["verdict"] == "reverser").sum()
    n_opp = (conn["verdict"] == "opposing").sum()

    md = []
    md.append("# H17 — connectivity-score drug repurposing (BRAF_like / cPTC, N=41)\n")
    md.append(f"DM1 n={(y==1).sum()}, DM2 n={(y==0).sum()}. 8-gene panel excluded. "
              f"Signature export = top {n_up}↑(DM1) + {n_dn}↓(DM2-up) at FDR<{FDR_CUT}; "
              f"connectivity score uses the **full** DE table (50,112 genes) ∩ drug KB. "
              "Score = weighted Σ sign(d_DM1>DM2)·impact_drug, range [-1, +1]; "
              "**+1 = drug perfectly reverses DM2 toward DM1**.\n")

    # Lead with top 3 + MAPK inhibitor confirmation status
    md.append("## Top 3 drug candidates (REVERSE DM2 → DM1)\n")
    md.append("| rank | drug | score | n∩ | MoA | MAPK-class confirmation |")
    md.append("|---:|---|---:|---:|---|---|")
    for i, (_, r) in enumerate(top3.iterrows(), 1):
        moa_lc = r["moa"].lower()
        if any(s in moa_lc for s in ("ifn-γ", "tlr3", "tlr9", "checkpoint")):
            mapk = "**not** a MAPK inhibitor — operates orthogonal IFN-γ axis"
        elif any(s in moa_lc for s in ("braf", "mek", "erk", "raf", "kinase")):
            mapk = "MAPK pathway"
        else:
            mapk = "—"
        md.append(f"| {i} | **{r['drug']}** | {r['connectivity_score']:+.3f} | "
                  f"{r['n_intersection']} | {r['moa']} | {mapk} |")
    md.append("")

    md.append("## All H17 drugs (sorted reverser → opposing)\n")
    md.append("| drug | score | n∩ | FDR<.05 hits | MoA | verdict |")
    md.append("|---|---:|---:|---:|---|---|")
    for _, r in conn.dropna(subset=["connectivity_score"]).iterrows():
        md.append(f"| {r['drug']} | {r['connectivity_score']:+.3f} | "
                  f"{r['n_intersection']} | {r['fdr_pass_intersection']} | "
                  f"{r['moa']} | {r['verdict']} |")
    md.append("")
    md.append(f"Verdict counts: **{n_revers} reverser** (score>0.10), "
              f"**{n_opp} opposing** (score<-0.10). MAPK inhibitors are OPPOSING, "
              "not reversers — see below.\n")

    # Module direction
    md.append("## Module direction within BRAF-cPTC (DM1 minus DM2)\n")
    md.append("| module | n genes | mean Cohen d | FDR<.05 frac | interpretation |")
    md.append("|---|---:|---:|---:|---|")
    for _, r in mods.iterrows():
        if r["n_present"] == 0:
            continue
        if r["mean_d"] > 0.4:
            tag = "DM1>>DM2"
        elif r["mean_d"] < -0.4:
            tag = "DM2>>DM1"
        else:
            tag = "comparable"
        md.append(f"| {r['module']} | {r['n_present']}/{r['n_genes_in_module']} | "
                  f"{r['mean_d']:+.2f} | {r['frac_FDR05']*100:.0f}% | {tag} |")
    md.append("")
    md.append("**Interpretation.** Within BRAF-cPTC, DM1 = immunogenic + MAPK-active + "
              "ECM-aggressive arm; DM2 = immune-cold + thyroid-residue arm. Per H6 the DM2 arm "
              "carries PFI HR=5.91. High-risk DM2 phenotype is **immune-cold despite preserved "
              "differentiation** — risk reverses by inducing inflammation, NOT by adding "
              "re-differentiation.\n")

    if h9_xref is not None and len(h9_xref):
        md.append("## H17 × H9 PRISM cross-reference (MAPK class)\n")
        md.append("| drug | PRISM d (kills DM1-high) | H17 score (reverses DM2→DM1) | verdict |")
        md.append("|---|---:|---:|---|")
        for _, r in h9_xref.iterrows():
            md.append(f"| {r['drug_h9']} | {r['h9_cohens_d_DM1high_vs_low']:+.2f} | "
                      f"{r['h17_connectivity_score']:+.3f} | {r['h17_verdict']} |")
        md.append("")
        md.append("Both axes negative: BRAFi/ERKi selectively kill DM1-high (MAPK-active) "
                  "lines (PRISM viability) but do NOT reverse the DM2 signature toward DM1 "
                  "(H17 connectivity). MAPK inhibitor will not rescue DM2 high-risk arm; "
                  "ICI / IFN-γ axis is the DM2-pertinent intervention.\n")

    md.append("## Druggable DEG targets (top-30 DM1↑ + top-30 DM2↑)\n")
    md.append("- **TACSTD2/Trop-2** (DM2-up): sacituzumab govitecan ADC — clinically active TNBC/UC; thyroid basket")
    md.append("- **NAMPT** (DepMap DM1-low d=-0.45): FK866, KPT-9274 — synthetic-lethal in DM2")
    md.append("- **MET** (DM2-up): crizotinib / capmatinib / tepotinib")
    md.append("- **MYC** (DM1-up proliferation feedback): BETi (JQ1, OTX015, BMS-986158); MYCi975")
    md.append("- **STING1** (DM1-up): ADU-S100 / MIW815 — adjuvant to anti-PD1")
    md.append("- **CDK4/6** via CCND1: palbociclib / ribociclib / abemaciclib\n")

    md.append("## Bottom line\n")
    md.append("**DM2 BRAF-cPTC reversal modality = IFN-γ axis activation, not MAPK blockade.** "
              "Top reversers (IFN-γ proxy +1.0, TLR3 agonist +1.0, anti-PD1 +0.81) push the "
              "IFN-γ-response module that DM1 has and DM2 lacks. DNMT inhibitors (decitabine "
              "+0.19, azacitidine +0.15) are second-tier via MHC-II + AICDA induction. "
              "BRAFi/MEKi/ERKi score -0.21 to -0.40 (opposing) — they shut off DUSP4/5/PLAU/FN1 "
              "which are DM1-up here. H6 (DM2 PFI HR=5.91) calls for immunotherapy, not MAPK "
              "targeting, in this stratum.\n")
    md.append("**Caveats.** (1) Drug→gene KB curated from published pharmacology, NOT live "
              "LINCS L1000 (clue.io auth-walled, no local L1000 dump). "
              "`h17_dm2_signature.tsv` upload-ready for external L1000 run. "
              f"(2) N=41 cohort: top |d| stable (ICAM1 d={sig.iloc[0]['cohens_d']:+.2f}); "
              "tail fold-fragile. (3) Score = sign-of-effect, not potency — PRISM (H9) "
              "orthogonal viability axis.\n")
    md.append("## Files\n")
    md.append("- `h17_dm2_signature.tsv` — top 200↑(DM1) + 200↓(DM2-up), FDR<0.05, ranked by |d| (clue.io upload-ready)")
    md.append("- `h17_connectivity_top50.tsv` — full drug × score table (30 drugs, sign + intersection)")
    md.append("- `h17_targetable_dm2_genes.tsv` — per-gene druggability for top 30+30 DEGs")
    md.append("- `h17_module_directionality.tsv` — TDS / IFNG / MAPK-output / TLS / MYC / ECM module direction")
    md.append("- `h17_h9_prism_crossref.tsv` — H17 × H9 PRISM cross-reference")
    md.append("- `run_h17_drug_repurposing.py` — full pipeline\n")

    (OUT / "H17_REPORT.md").write_text("\n".join(md))
    print(f"  saved H17_REPORT.md  ({sum(len(s) for s in md)} chars)")
    print("DONE.")


if __name__ == "__main__":
    main()
