#!/usr/bin/env python3
"""H20 — Paper 1+2 BRAF Nature sprint: HT-13 axis pan-cancer prognostic generalization.

Reads:
  - /data/thca/repo_data/raw/TCGA_pancan/pancan_geneExp.gz  (HGNC-symbol RSEM log2)
  - /data/thca/repo_data/raw/TCGA_pancan/survival.tsv       (Liu 2018 PanCancerAtlas)
  - project/results/paper11_pancancer/pancan_dm1_scored.tsv (sample → lineage map + DM1_like)
  - project/results/paper11_pancancer/phase_B_survival/cox_per_lineage.tsv  (DM1_like Cox; reuse for scatter)

For each lineage:
  - z-score HT-13 expression within-cohort, mean → HT13_score
  - Cox PH on OS, DSS, PFI with covariates {age, stage} where available
  - BH-FDR across 33 lineages per endpoint

Outputs to the same H20 dir.
"""
from __future__ import annotations
import gzip, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project")
P11 = ROOT / "results" / "paper11_pancancer"
SCORES = P11 / "pancan_dm1_scored.tsv"
PANCAN_EXPR = Path("/data/thca/repo_data/raw/TCGA_pancan/pancan_geneExp.gz")
SURV = Path("/data/thca/repo_data/raw/TCGA_pancan/survival.tsv")
DM1_COX = P11 / "phase_B_survival" / "cox_per_lineage.tsv"
OUT = ROOT / "results" / "p2_braf_nature_sprint_2026_05_09" / "h20_pancancer"
OUT.mkdir(parents=True, exist_ok=True)

HT13 = [
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
    "CD79A", "CD79B", "MS4A1", "AICDA", "CXCL13", "CCR6", "IFNG",
]

ICI_LINEAGES = {
    "skin cutaneous melanoma", "uveal melanoma",
    "lung adenocarcinoma", "lung squamous cell carcinoma",
    "liver hepatocellular carcinoma", "head & neck squamous cell carcinoma",
    "bladder urothelial carcinoma", "esophageal carcinoma",
    "stomach adenocarcinoma", "kidney clear cell carcinoma",
    "kidney papillary cell carcinoma",
}


def stream_pancan_expr(path: Path, target_genes):
    target = set(target_genes)
    print(f"[H20] streaming {path} for {len(target)} genes …")
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        sample_ids = header[1:]
        out = {}
        n = 0
        for line in f:
            n += 1
            if n % 5000 == 0:
                print(f"  …{n}")
            tab = line.find("\t")
            if tab < 0:
                continue
            gene = line[:tab].strip()
            if gene in target:
                vals = line[tab + 1:].rstrip("\n").split("\t")
                out[gene] = pd.to_numeric(pd.Series(vals), errors="coerce").values
                if len(out) == len(target):
                    print(f"  matched all {len(target)} genes by row {n}")
                    break
    df = pd.DataFrame(out, index=sample_ids)
    df.index.name = "sample"
    print(f"[H20] expr matrix: {df.shape}, missing genes: {sorted(target - set(out))}")
    return df


def parse_stage(s):
    if pd.isna(s) or not isinstance(s, str):
        return np.nan
    s = s.lower().strip()
    if "stage iv" in s:
        return 4
    if "stage iii" in s:
        return 3
    if "stage ii" in s:
        return 2
    if "stage i" in s:
        return 1
    return np.nan


def bh_fdr(pvals: pd.Series) -> pd.Series:
    p = pvals.astype(float)
    out = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.notna()
    if valid.sum() == 0:
        return out
    pv = p[valid].values
    order = np.argsort(pv)
    ranked = pv[order]
    m = len(ranked)
    bh = ranked * m / np.arange(1, m + 1)
    bh = np.minimum.accumulate(bh[::-1])[::-1]
    bh = np.clip(bh, 0, 1)
    fdr = np.empty_like(bh)
    fdr[order] = bh
    out.loc[valid] = fdr
    return out


def main():
    if not SCORES.exists() or not PANCAN_EXPR.exists() or not SURV.exists():
        sys.exit("[H20] missing required inputs")

    expr = stream_pancan_expr(PANCAN_EXPR, HT13)
    found = [g for g in HT13 if g in expr.columns]
    print(f"[H20] found {len(found)}/13: {found}")

    scores = pd.read_csv(SCORES, sep="\t")
    surv = pd.read_csv(SURV, sep="\t", low_memory=False)
    print(f"[H20] scores={len(scores)} surv={len(surv)}")

    # join expr + lineage + DM1_like
    expr = expr.reset_index()
    df = scores.merge(expr, on="sample", how="left")
    print(f"[H20] joined scores+expr: {len(df)}")
    n_with_expr = df[found].notna().any(axis=1).sum()
    print(f"[H20] samples with HT-13 expression: {n_with_expr}")

    # within-lineage z then mean → HT13_score
    z_blocks = []
    for lin, sub in df.groupby("lineage"):
        sub = sub.copy()
        for g in found:
            col = sub[g]
            mu, sd = col.mean(), col.std()
            sub[g + "_z"] = (col - mu) / (sd if (sd and sd > 0) else 1.0)
        z_cols = [g + "_z" for g in found]
        sub["HT13_score"] = sub[z_cols].mean(axis=1)
        z_blocks.append(sub[["sample", "lineage", "DM1_like", "HT13_score"]])
    big = pd.concat(z_blocks, ignore_index=True)
    big.to_csv(OUT / "h20_ht13_scored_per_sample.tsv", sep="\t", index=False)

    # Merge survival
    big = big.merge(surv, on="sample", how="left")
    big["age"] = pd.to_numeric(big.get("age_at_initial_pathologic_diagnosis"), errors="coerce")
    big["stage_int"] = big.get("ajcc_pathologic_tumor_stage").apply(parse_stage)

    try:
        from lifelines import CoxPHFitter
    except ImportError:
        sys.exit("[H20] lifelines not installed")

    rows = []
    for lin, sub in big.groupby("lineage"):
        for outcome, t_col, e_col in [("OS", "OS.time", "OS"),
                                       ("DSS", "DSS.time", "DSS"),
                                       ("PFI", "PFI.time", "PFI")]:
            cox = sub[[outcome, t_col, "HT13_score", "DM1_like", "age", "stage_int"]].copy()
            cox = cox.dropna(subset=[outcome, t_col, "HT13_score"])
            cox = cox[cox[t_col] > 0]
            if len(cox) < 20 or cox[outcome].sum() < 5:
                continue
            covars = ["HT13_score"]
            if cox["age"].notna().sum() / len(cox) > 0.5:
                covars.append("age")
            if cox["stage_int"].notna().sum() / len(cox) > 0.3:
                covars.append("stage_int")
            cox_uv = cox[[outcome, t_col] + covars].dropna()
            if len(cox_uv) < 20:
                continue
            cph = CoxPHFitter(penalizer=0.001)
            try:
                cph.fit(cox_uv, duration_col=t_col, event_col=outcome)
                s = cph.summary.loc["HT13_score"]
                row = {
                    "lineage": lin, "outcome": outcome,
                    "n": len(cox_uv), "events": int(cox_uv[outcome].sum()),
                    "covariates": "+".join(covars),
                    "HR": float(s["exp(coef)"]),
                    "HR_lower95": float(s["exp(coef) lower 95%"]),
                    "HR_upper95": float(s["exp(coef) upper 95%"]),
                    "p": float(s["p"]),
                    "concordance_HT13only": np.nan,
                    "concordance_full": float(cph.concordance_index_),
                    "ICI_lineage": lin in ICI_LINEAGES,
                }
                # HT13 alone C-index
                try:
                    cph_only = CoxPHFitter(penalizer=0.001)
                    cph_only.fit(cox[[outcome, t_col, "HT13_score"]].dropna(),
                                 duration_col=t_col, event_col=outcome)
                    row["concordance_HT13only"] = float(cph_only.concordance_index_)
                except Exception:
                    pass
                # DM1+HT13 combined C-index (no other covars, like-for-like)
                try:
                    combined = cox[[outcome, t_col, "HT13_score", "DM1_like"]].dropna()
                    if len(combined) >= 20 and combined[outcome].sum() >= 5:
                        cph_c = CoxPHFitter(penalizer=0.001)
                        cph_c.fit(combined, duration_col=t_col, event_col=outcome)
                        row["concordance_HT13_plus_DM1"] = float(cph_c.concordance_index_)
                        row["DM1_HR_in_combined"] = float(cph_c.summary.loc["DM1_like", "exp(coef)"])
                        row["DM1_p_in_combined"] = float(cph_c.summary.loc["DM1_like", "p"])
                        row["HT13_HR_in_combined"] = float(cph_c.summary.loc["HT13_score", "exp(coef)"])
                        row["HT13_p_in_combined"] = float(cph_c.summary.loc["HT13_score", "p"])
                    # DM1-alone C-index
                    dm1_only = cox[[outcome, t_col, "DM1_like"]].dropna()
                    if len(dm1_only) >= 20 and dm1_only[outcome].sum() >= 5:
                        cph_d = CoxPHFitter(penalizer=0.001)
                        cph_d.fit(dm1_only, duration_col=t_col, event_col=outcome)
                        row["concordance_DM1only"] = float(cph_d.concordance_index_)
                except Exception as e:
                    row["combined_error"] = str(e)[:80]
                rows.append(row)
            except Exception as e:
                rows.append({
                    "lineage": lin, "outcome": outcome,
                    "n": len(cox_uv), "events": int(cox_uv[outcome].sum()),
                    "covariates": "+".join(covars),
                    "HR": np.nan, "p": np.nan, "error": str(e)[:80],
                    "ICI_lineage": lin in ICI_LINEAGES,
                })

    res = pd.DataFrame(rows)

    # FDR per outcome
    pieces = []
    for outc, g in res.groupby("outcome"):
        g = g.copy()
        g["fdr"] = bh_fdr(g["p"])
        pieces.append(g)
    res = pd.concat(pieces, ignore_index=True).sort_values(["outcome", "fdr"])
    res.to_csv(OUT / "h20_pancancer_ht_cox.tsv", sep="\t", index=False)
    print(f"[H20] cox table → h20_pancancer_ht_cox.tsv ({len(res)} rows)")

    # Build DM1 vs HT13 scatter (per lineage, OS preferred)
    dm1_cox = pd.read_csv(DM1_COX, sep="\t")
    scatter_rows = []
    for outc in ("OS", "DSS", "PFI"):
        ht = res[res.outcome == outc].set_index("lineage")
        dm = dm1_cox[dm1_cox.outcome == outc].set_index("lineage")
        common = sorted(set(ht.index) & set(dm.index))
        for lin in common:
            scatter_rows.append({
                "outcome": outc, "lineage": lin,
                "HT13_HR": ht.loc[lin, "HR"], "HT13_p": ht.loc[lin, "p"],
                "HT13_fdr": ht.loc[lin, "fdr"],
                "DM1_HR": dm.loc[lin, "HR"], "DM1_p": dm.loc[lin, "p"],
                "DM1_fdr": dm.loc[lin, "fdr"] if "fdr" in dm.columns else np.nan,
                "n": ht.loc[lin, "n"], "events": ht.loc[lin, "events"],
                "ICI_lineage": lin in ICI_LINEAGES,
            })
    scat = pd.DataFrame(scatter_rows)
    scat.to_csv(OUT / "h20_dm1_vs_ht_scatter.tsv", sep="\t", index=False)
    print(f"[H20] scatter table → h20_dm1_vs_ht_scatter.tsv ({len(scat)} rows)")

    # Spearman log(HR) DM1 vs HT13 per outcome
    summary = {"genes_found": found, "genes_missing": sorted(set(HT13) - set(found)),
               "n_lineages_per_outcome": {}, "n_significant_fdr10": {},
               "ici_significant_fdr10": {},
               "dm1_vs_ht13_logHR_correlation": {}, "headline": {}}
    for outc in ("OS", "DSS", "PFI"):
        g = res[res.outcome == outc].dropna(subset=["HR", "fdr"])
        summary["n_lineages_per_outcome"][outc] = int(len(g))
        sig = g[g.fdr < 0.1]
        summary["n_significant_fdr10"][outc] = int(len(sig))
        ici_sig = sig[sig.ICI_lineage]
        summary["ici_significant_fdr10"][outc] = {
            "n": int(len(ici_sig)),
            "lineages": ici_sig["lineage"].tolist(),
        }
        s2 = scat[scat.outcome == outc].dropna(subset=["HT13_HR", "DM1_HR"])
        if len(s2) >= 5:
            r, p = stats.spearmanr(np.log(s2.HT13_HR), np.log(s2.DM1_HR))
            summary["dm1_vs_ht13_logHR_correlation"][outc] = {
                "n_lineages": int(len(s2)), "spearman_r": float(r), "p": float(p),
            }
        if len(sig):
            summary["headline"][outc] = {
                "n_sig": int(len(sig)),
                "top": sig.sort_values("fdr").head(8)[
                    ["lineage", "HR", "HR_lower95", "HR_upper95", "p", "fdr",
                     "n", "events", "ICI_lineage"]
                ].to_dict(orient="records"),
            }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str)[:3500])

    # Forest + scatter plots
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        for outc in ("OS", "DSS", "PFI"):
            g = res[res.outcome == outc].dropna(subset=["HR"]).copy()
            if not len(g):
                continue
            g = g.sort_values("HR")
            fig, ax = plt.subplots(figsize=(7, max(5, 0.22 * len(g))))
            y = np.arange(len(g))
            ax.errorbar(g.HR, y,
                        xerr=[g.HR - g.HR_lower95, g.HR_upper95 - g.HR],
                        fmt="o", capsize=3, color="black", ecolor="grey")
            colors = ["red" if i else "black" for i in g.ICI_lineage]
            for i, (yi, hr, lin, fdr, ici) in enumerate(zip(y, g.HR, g.lineage, g.fdr, g.ICI_lineage)):
                mk = "**" if fdr < 0.05 else ("*" if fdr < 0.1 else "")
                ax.text(max(g.HR_upper95.iloc[i], hr) * 1.05, yi, mk, va="center")
            ax.axvline(1, color="r", ls=":")
            ax.set_yticks(y)
            ax.set_yticklabels([f"{l} {'[ICI]' if i else ''}" for l, i in zip(g.lineage, g.ICI_lineage)],
                              fontsize=7)
            ax.set_xscale("log")
            ax.set_xlabel(f"HR (HT13_score, {outc} Cox; * FDR<0.1, ** FDR<0.05)")
            ax.set_title(f"H20 — Pan-cancer HT-13 axis prognostic effect ({outc})")
            plt.tight_layout()
            plt.savefig(OUT / f"forest_HT13_{outc}.png", dpi=140)
            plt.close()
        # scatter
        for outc in ("OS", "DSS", "PFI"):
            s2 = scat[scat.outcome == outc].dropna(subset=["HT13_HR", "DM1_HR"])
            if len(s2) < 4:
                continue
            fig, ax = plt.subplots(figsize=(6, 6))
            for ici, sub in s2.groupby("ICI_lineage"):
                ax.scatter(np.log2(sub.DM1_HR), np.log2(sub.HT13_HR),
                           label=("ICI lineage" if ici else "other"),
                           c=("red" if ici else "steelblue"), s=40, alpha=0.85)
                for _, r in sub.iterrows():
                    ax.annotate(r.lineage[:14], (np.log2(r.DM1_HR), np.log2(r.HT13_HR)),
                                fontsize=6)
            ax.axhline(0, color="grey", ls=":")
            ax.axvline(0, color="grey", ls=":")
            ax.set_xlabel("log2 HR  (DM1_like Cox)")
            ax.set_ylabel("log2 HR  (HT13_score Cox)")
            corr = summary["dm1_vs_ht13_logHR_correlation"].get(outc, {})
            ax.set_title(f"H20 DM1 vs HT-13 effect-size scatter ({outc})  ρ={corr.get('spearman_r', np.nan):.2f}")
            ax.legend(loc="best", fontsize=8)
            plt.tight_layout()
            plt.savefig(OUT / f"scatter_dm1_vs_ht13_{outc}.png", dpi=140)
            plt.close()
        print("[H20] plots written")
    except Exception as e:
        print(f"[H20] plot failed: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
