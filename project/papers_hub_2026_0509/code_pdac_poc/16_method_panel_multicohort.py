"""
Apply our full method (Moffitt + KRAS allele × Moffitt × paradox + Cox)
to EVERY fetchable PDAC cBioPortal cohort with sufficient RNA + survival.

Per cohort:
  1. clinical fetch
  2. KRAS mutation fetch + allele assignment
  3. mRNA z-score fetch on Moffitt 47-gene panel + immune modules
  4. Moffitt classifier
  5. Cox PH (allele-binary covariates + age)
  6. Paradox quantification (basal × inflamed × suppress)

Outputs:
  /data/pdac_poc/results/method_panel/<study>.json
  /data/pdac_poc/results/method_panel/META_ANALYSIS.json — fixed-effect HR pooling
"""

import json
import time
import math
from pathlib import Path
import pandas as pd
import numpy as np
import requests
from lifelines import CoxPHFitter

API = "https://www.cbioportal.org/api"
OUT = Path("/data/pdac_poc/results/method_panel")
OUT.mkdir(parents=True, exist_ok=True)
S = requests.Session()
S.headers.update({"Accept": "application/json"})

COHORTS = [
    ("paad_tcga_pan_can_atlas_2018", "TCGA-PAAD anchor"),
    ("paad_qcmg_uq_2016",            "Bailey 2016 Nature QCMG"),
    ("paad_icgc",                    "ICGC PACA"),
    ("paad_utsw_2015",               "Witkiewicz 2015 NatComm"),
    ("paad_cptac_2021",              "Cao 2021 CPTAC Cell"),
    ("paad_iatlas_prince_2022",      "PRINCE Trial Nat Med"),
    ("pancreas_msk_2024",            "MSK 2024 Cancer Cell"),
    ("pdac_msk_2024",                "MSK 2024 Nat Med"),
    ("paad_msk_2025",                "MSK 2025"),
    ("paad_tcga",                    "TCGA Firehose Legacy"),
]

MOFFITT_BASAL = ["VGLL1","UCA1","S100A2","LY6D","SPRR3","SPRR1B","LEMD1","LYPD3",
                 "KRT15","DHRS9","AREG","CST6","SERPINB3","KRT6A","SERPINB4",
                 "FAM83A","SCEL","FGFBP1","KRT7","KRT17","GPR87","TNS4","SLC2A1"]
MOFFITT_CLASS = ["BTNL8","FAM3D","AGR3","CTSE","LYZ","TFF2","TFF1","ANXA10",
                 "LGALS4","ECT2","CLRN3","MYO1A","CLDN18","LRRC31","TFF3","CDX2",
                 "SERPINA10","VSIG2","TSPAN8","ST6GALNAC1","AGR2","TOX3"]
PARA_INFL = ["CXCL13","IFNG","HLA-DRA","HLA-DRB1","CD8A"]
PARA_SUPP = ["CD68","CD163","CSF1R","CD274","PDCD1","TIGIT","LAG3"]


def safe_get(url, **kw):
    try:
        r = S.get(url, params=kw, timeout=60); r.raise_for_status(); return r.json()
    except Exception:
        return None


def safe_post(url, body, **kw):
    try:
        r = S.post(url, json=body, params=kw, timeout=120,
                    headers={"Content-Type":"application/json"}); r.raise_for_status(); return r.json()
    except Exception:
        return None


def list_profiles(study):
    return safe_get(f"{API}/studies/{study}/molecular-profiles") or []


def find_z_profile(profiles):
    for p in profiles:
        pid = p.get("molecularProfileId","")
        if pid.endswith("_mrna_median_Zscores") or pid.endswith("_rna_seq_v2_mrna_median_Zscores") \
           or pid.endswith("_rna_seq_mrna_median_Zscores") \
           or "Zscores" in pid:
            return pid
    return None


def fetch_clinical(study):
    samples = safe_get(f"{API}/studies/{study}/samples") or []
    if not samples: return None
    sdf = pd.DataFrame(samples)
    pcl = safe_get(f"{API}/studies/{study}/clinical-data", clinicalDataType="PATIENT") or []
    scl = safe_get(f"{API}/studies/{study}/clinical-data", clinicalDataType="SAMPLE") or []
    pdf = (pd.DataFrame(pcl).pivot_table(index="patientId", columns="clinicalAttributeId",
            values="value", aggfunc="first").reset_index() if pcl else pd.DataFrame())
    sdf2 = (pd.DataFrame(scl).pivot_table(index="sampleId", columns="clinicalAttributeId",
            values="value", aggfunc="first").reset_index() if scl else pd.DataFrame())
    if not sdf2.empty: sdf = sdf.merge(sdf2, on="sampleId", how="left")
    if not pdf.empty:  sdf = sdf.merge(pdf, on="patientId", how="left")
    return sdf


def fetch_kras(study):
    g = safe_post(f"{API}/genes/fetch", body=["KRAS"], geneIdType="HUGO_GENE_SYMBOL")
    if not g: return None
    eid = g[0]["entrezGeneId"]
    body = {"entrezGeneIds":[eid], "sampleListId": f"{study}_sequenced"}
    muts = safe_post(f"{API}/molecular-profiles/{study}_mutations/mutations/fetch",
                     body=body, projection="DETAILED")
    return muts or []


def fetch_z(study, profile_id, genes):
    if not profile_id: return pd.DataFrame()
    g = safe_post(f"{API}/genes/fetch", body=genes, geneIdType="HUGO_GENE_SYMBOL") or []
    if not g: return pd.DataFrame()
    eid_to_sym = {x["entrezGeneId"]: x["hugoGeneSymbol"] for x in g}
    sample_list = f"{study}_rna_seq_v2_mrna"
    # fall-back to all_complete or _all if needed
    body = {"entrezGeneIds": list(eid_to_sym.keys()), "sampleListId": sample_list}
    data = safe_post(f"{API}/molecular-profiles/{profile_id}/molecular-data/fetch", body=body)
    if not data:
        # try alternative sample list
        for alt in [f"{study}_rna_seq_mrna", f"{study}_complete", f"{study}_all"]:
            body["sampleListId"] = alt
            data = safe_post(f"{API}/molecular-profiles/{profile_id}/molecular-data/fetch", body=body)
            if data: break
    if not data: return pd.DataFrame()
    df = pd.DataFrame(data)
    df["gene"] = df["entrezGeneId"].map(eid_to_sym)
    return df.pivot_table(index="sampleId", columns="gene", values="value", aggfunc="first")


def assign_allele(muts, sample_ids):
    out = pd.Series("WT", index=list(sample_ids))
    if not muts: return out
    df = pd.DataFrame(muts)
    if "gene" in df.columns:
        df["hugo"] = df["gene"].apply(lambda g: g.get("hugoGeneSymbol") if isinstance(g, dict) else None)
    df = df[(df.get("mutationType","").fillna("") == "Missense_Mutation") &
            (df.get("hugo","KRAS") == "KRAS")]
    for s, g in df.groupby("sampleId"):
        if s not in out.index: continue
        pcs = g["proteinChange"].dropna().tolist()
        for p in ["G12D","G12V","G12R","G12C","Q61H","Q61R"]:
            if any(pc == p for pc in pcs): out.loc[s] = p; break
        else:
            if any(pc.startswith("G12") for pc in pcs): out.loc[s] = "G12_other"
            elif any(pc.startswith("Q61") for pc in pcs): out.loc[s] = "Q61"
            elif pcs: out.loc[s] = "KRAS_other"
    return out


def best_os_cols(clin):
    pairs = [("OS_MONTHS","OS_STATUS"),("OS_MONTHS","STATUS"),("OS_MONTHS","VITAL_STATUS"),
             ("DAYS_TO_LAST_FOLLOWUP","STATUS"),("DAYS_TO_LAST_FOLLOWUP","VITAL_STATUS"),
             ("OVERALL_SURVIVAL_MONTHS","OVERALL_SURVIVAL_STATUS")]
    cols = set(clin.columns)
    for t,s in pairs:
        if t in cols and s in cols: return t, s
    return None, None


def parse_event(v):
    if isinstance(v,str):
        v=v.lower()
        if v.startswith("1") or v in ("dead","deceased","yes"): return 1
        if v.startswith("0") or v in ("alive","living","no"): return 0
    try: return int(float(v))
    except: return None


def cohens_d(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if len(a)<2 or len(b)<2: return float("nan")
    v = ((len(a)-1)*a.var(ddof=1) + (len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2)
    sp = np.sqrt(v) if v>0 else np.nan
    return (a.mean() - b.mean())/sp if sp else float("nan")


def analyze(study, label):
    out = {"study": study, "label": label}
    clin = fetch_clinical(study)
    if clin is None or clin.empty:
        out["status"] = "no clinical"; return out
    muts = fetch_kras(study)
    allele = assign_allele(muts, clin["sampleId"].tolist())
    out["n_samples"] = int(len(clin))
    out["allele_counts"] = allele.value_counts().to_dict()

    # Cox
    t_col, s_col = best_os_cols(clin)
    if t_col and s_col:
        df = clin.set_index("sampleId").copy()
        df["allele"] = df.index.map(allele)
        df["t"] = pd.to_numeric(df[t_col], errors="coerce")
        df["e"] = df[s_col].apply(parse_event)
        if t_col.endswith("DAYS") or t_col.endswith("FOLLOWUP"):
            df["t"] = df["t"]/30.4375
        df = df[(df["t"]>0) & df["e"].notna()].dropna(subset=["t","e"])
        df["G12D"] = (df["allele"]=="G12D").astype(int)
        df["G12V"] = (df["allele"]=="G12V").astype(int)
        df["G12R"] = (df["allele"]=="G12R").astype(int)
        if len(df) >= 30:
            try:
                cph = CoxPHFitter(penalizer=0.01).fit(df[["t","e","G12D","G12V","G12R"]],
                                                       duration_col="t", event_col="e")
                s = cph.summary[["coef","exp(coef)","p","exp(coef) lower 95%","exp(coef) upper 95%"]].round(4)
                out["cox"] = s.reset_index().to_dict(orient="records")
                out["concordance"] = float(cph.concordance_index_)
                g = next((r for r in out["cox"] if r["covariate"]=="G12D"), None)
                if g:
                    out["G12D_HR"] = g["exp(coef)"]
                    out["G12D_p"]  = g["p"]
                    out["G12D_log_HR"] = math.log(g["exp(coef)"]) if g["exp(coef)"]>0 else None
                    # se(log HR) ≈ (log(upper) - log(lower)) / (2*1.96)
                    if g["exp(coef) upper 95%"]>0 and g["exp(coef) lower 95%"]>0:
                        out["G12D_se_log_HR"] = (math.log(g["exp(coef) upper 95%"]) -
                                                  math.log(g["exp(coef) lower 95%"])) / 3.92
                out["n_with_OS"] = int(len(df))
            except Exception as e:
                out["cox_error"] = str(e)

    # mRNA
    profiles = list_profiles(study)
    z_profile = find_z_profile(profiles)
    out["z_profile"] = z_profile
    if z_profile:
        genes = sorted(set(MOFFITT_BASAL + MOFFITT_CLASS + PARA_INFL + PARA_SUPP))
        z = fetch_z(study, z_profile, genes)
        if not z.empty:
            avail_b = [g for g in MOFFITT_BASAL if g in z.columns]
            avail_c = [g for g in MOFFITT_CLASS if g in z.columns]
            if avail_b and avail_c:
                z["score_basal"] = z[avail_b].mean(axis=1)
                z["score_class"] = z[avail_c].mean(axis=1)
                z["moffitt"] = np.where(z["score_basal"] > z["score_class"], "basal", "classical")
                out["n_with_RNA"] = int(len(z))
                out["fraction_basal"] = float((z["moffitt"]=="basal").mean())
                # paradox
                avail_inf = [g for g in PARA_INFL if g in z.columns]
                avail_sup = [g for g in PARA_SUPP if g in z.columns]
                if avail_inf and avail_sup:
                    z["inflamed"] = z[avail_inf].mean(axis=1)
                    z["suppress"] = z[avail_sup].mean(axis=1)
                    b = z.query("moffitt=='basal'"); c = z.query("moffitt=='classical'")
                    out["paradox_d_inflamed"]  = round(cohens_d(b["inflamed"], c["inflamed"]),3)
                    out["paradox_d_suppress"]  = round(cohens_d(b["suppress"], c["suppress"]),3)
                    out["fraction_basal_BOTH_high"] = round(
                        float(((b["inflamed"]>0) & (b["suppress"]>0)).mean()), 3)
                    out["fraction_classical_BOTH_high"] = round(
                        float(((c["inflamed"]>0) & (c["suppress"]>0)).mean()), 3)
                # KRAS allele × Moffitt
                z["allele"] = z.index.map(allele)
                ct = pd.crosstab(z["allele"], z["moffitt"])
                if "basal" in ct.columns:
                    bf = {}
                    for al in ["G12D","G12V","G12R","G12C"]:
                        if al in ct.index:
                            tot = ct.loc[al].sum()
                            bf[al] = round(ct.loc[al,"basal"]/tot, 3) if tot else None
                    out["basal_fraction_by_allele"] = bf
    return out


def fixed_effect_meta(rows):
    """Inverse-variance fixed-effect meta-analysis on log HR."""
    obs = [(r["G12D_log_HR"], r["G12D_se_log_HR"], r["label"], r["n_with_OS"])
            for r in rows
            if r.get("G12D_log_HR") is not None and r.get("G12D_se_log_HR") and r["G12D_se_log_HR"] > 0]
    if not obs:
        return None
    weights = [1/(s*s) for _, s, *_ in obs]
    wsum = sum(weights)
    pooled_log = sum(w*l for w,(l,*_) in zip(weights,obs))/wsum
    pooled_se = math.sqrt(1/wsum)
    pooled_HR = math.exp(pooled_log)
    z = pooled_log/pooled_se
    # two-sided p
    from math import erf
    p = 2*(1 - 0.5*(1+erf(abs(z)/math.sqrt(2))))
    ci_lo = math.exp(pooled_log - 1.96*pooled_se)
    ci_hi = math.exp(pooled_log + 1.96*pooled_se)
    return {"pooled_HR": round(pooled_HR,3), "CI_lo": round(ci_lo,3), "CI_hi": round(ci_hi,3),
            "p_value": float(f"{p:.4g}"), "z": round(z,3),
            "n_studies": len(obs), "n_total_OS": sum(o[3] for o in obs),
            "studies_used": [o[2] for o in obs]}


def main():
    rows = []
    for sid, label in COHORTS:
        print(f"\n[16] {sid}  ({label})")
        try:
            r = analyze(sid, label)
            rows.append(r)
            with open(OUT/f"{sid}.json","w") as f:
                json.dump(r,f,indent=2,default=str)
            tags = []
            if r.get("G12D_HR") is not None:
                tags.append(f"HR={r['G12D_HR']:.2f} p={r['G12D_p']:.3f}")
            if r.get("fraction_basal") is not None:
                tags.append(f"basal={r['fraction_basal']*100:.0f}%")
            if r.get("paradox_d_suppress") is not None:
                tags.append(f"paradox_d_suppress={r['paradox_d_suppress']:+.2f}")
            print(f"   n={r.get('n_samples')}  OS_n={r.get('n_with_OS','—')}  RNA_n={r.get('n_with_RNA','—')}  {' · '.join(tags) if tags else 'no_signal'}")
        except Exception as e:
            print(f"   ERR {sid}: {e}")
            rows.append({"study": sid, "label": label, "error": str(e)})
        time.sleep(0.4)

    meta = fixed_effect_meta(rows)
    out = {"generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
           "n_cohorts_attempted": len(COHORTS),
           "n_cohorts_with_HR": sum(1 for r in rows if r.get("G12D_HR") is not None),
           "n_cohorts_with_RNA": sum(1 for r in rows if r.get("n_with_RNA")),
           "meta_analysis": meta,
           "rows": rows}
    with open(OUT/"META_ANALYSIS.json","w") as f:
        json.dump(out, f, indent=2, default=str)

    print("\n" + "="*70)
    print("[16] CROSS-COHORT METHOD-PANEL SUMMARY")
    print("="*70)
    print(f"{'cohort':45s} {'n_samp':>7s} {'n_OS':>6s} {'n_RNA':>6s} {'G12D_HR':>9s} {'paradox':>9s}")
    for r in rows:
        hr = f"{r['G12D_HR']:.2f}*" if (r.get('G12D_HR') and r.get('G12D_p',1)<0.05) else (f"{r['G12D_HR']:.2f}" if r.get('G12D_HR') else "—")
        par = f"{r['paradox_d_suppress']:+.2f}" if r.get('paradox_d_suppress') is not None else "—"
        print(f"{(r.get('label') or r.get('study'))[:45]:45s} {(r.get('n_samples') or 0):>7d} "
              f"{(r.get('n_with_OS') or 0):>6d} {(r.get('n_with_RNA') or 0):>6d} {hr:>9s} {par:>9s}")
    print()
    if meta:
        print(f"  META-ANALYSIS  pooled HR = {meta['pooled_HR']}  "
              f"95% CI {meta['CI_lo']}–{meta['CI_hi']}  "
              f"p = {meta['p_value']}  ({meta['n_studies']} studies, n={meta['n_total_OS']})")


if __name__ == "__main__":
    main()
