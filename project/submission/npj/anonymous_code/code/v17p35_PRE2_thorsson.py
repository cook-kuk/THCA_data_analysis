"""v17p35 PRE-2 — Thorsson immune subtype via cBioPortal."""
from __future__ import annotations
import json, requests
from pathlib import Path
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact

PROJECT = Path("/opt/thyroid-dash/project")
TAB = PROJECT / "results/v17p35/tables"
FIG = PROJECT / "results/v17p35/figs"
UA = {"User-Agent": "thca-v17p35", "Accept": "application/json"}


def main():
    url = ("https://www.cbioportal.org/api/studies/thca_tcga_pan_can_atlas_2018/"
           "clinical-data?clinicalDataType=PATIENT&projection=DETAILED")
    r = requests.get(url, headers=UA, timeout=60)
    print(f"cBioPortal status: {r.status_code}, bytes: {len(r.content)}")
    if r.status_code != 200:
        print("FAIL")
        return
    data = r.json()
    df = pd.DataFrame(data)
    print(f"  rows: {len(df)}, attrs: {df['clinicalAttributeId'].nunique() if 'clinicalAttributeId' in df.columns else 'NA'}")
    print(f"  unique attrs (first 30): {sorted(df['clinicalAttributeId'].unique())[:30] if 'clinicalAttributeId' in df.columns else 'NA'}")

    # filter for Thorsson immune subtype attribute
    target_attrs = [a for a in df["clinicalAttributeId"].unique()
                    if "IMMUNE" in a.upper() or "SUBTYPE" in a.upper() or "THORSSON" in a.upper()]
    print(f"  candidate attrs: {target_attrs}")

    if not target_attrs:
        print("  no immune attribute found, dumping all distinct attr values for SUBTYPE")
        return

    immune_df = df[df["clinicalAttributeId"].isin(target_attrs)].copy()
    print(f"  immune attr rows: {len(immune_df)}")
    print(f"  values: {immune_df['value'].value_counts().head(15).to_dict()}")

    # Pivot to per-patient
    wide = immune_df.pivot_table(index="patientId", columns="clinicalAttributeId", values="value", aggfunc="first")
    wide.to_csv(TAB / "FIX2_thorsson_subtype_v2.tsv", sep="\t")
    print(f"  wrote FIX2_thorsson_subtype_v2.tsv ({len(wide)} patients)")

    # Merge with DM cluster labels
    dm_path = PROJECT / "results/v17p3/tables/A2_dm_score_full_cohort.tsv"
    if not dm_path.exists():
        print("  A2 dm_score not found, skip Fisher")
        return
    dm = pd.read_csv(dm_path, sep="\t")
    # patient = first 12 chars of TCGA barcode
    dm["patient_id"] = dm["sample"].str[:12].str.upper()
    wide.index = wide.index.str.upper()
    merged = dm.merge(wide, left_on="patient_id", right_index=True, how="inner")
    print(f"  DM × Thorsson merged: {len(merged)} samples")

    # Cross-tab + Fisher
    immune_col = [c for c in wide.columns if "SUBTYPE_IMMUNE" in c.upper() or "IMMUNE_MODEL" in c.upper()]
    if not immune_col:
        immune_col = wide.columns[:1].tolist()
    if not immune_col or "dm_like" not in merged.columns:
        return
    col = immune_col[0]
    print(f"  Using attribute: {col}")
    ct = pd.crosstab(merged["dm_like"], merged[col])
    print(f"\nCross-tab DM × Thorsson:\n{ct}")
    ct.to_csv(TAB / "FIX2_thorsson_dm_distribution.tsv", sep="\t")

    if ct.shape[0] >= 2 and ct.shape[1] >= 2:
        chi2, p, dof, expected = chi2_contingency(ct)
        print(f"\nChi2 p = {p:.4g}, dof = {dof}")

    # Per-subtype fisher
    fisher_rows = []
    for subt in ct.columns:
        ct2 = pd.DataFrame({
            subt: ct[subt],
            "other": ct.sum(axis=1) - ct[subt],
        })
        if ct2.shape[0] == 2 and ct2.shape[1] == 2:
            try:
                odds, pf = fisher_exact(ct2.values)
                fisher_rows.append({"subtype": subt, "odds_ratio": odds, "pvalue": pf})
            except Exception:
                pass
    if fisher_rows:
        fisher_df = pd.DataFrame(fisher_rows).sort_values("pvalue")
        fisher_df.to_csv(TAB / "FIX2_thorsson_fisher.tsv", sep="\t", index=False)
        print(f"\nFisher per subtype:\n{fisher_df.to_string(index=False)}")

    summary = {
        "method": "cBioPortal_thca_tcga_pan_can_atlas_2018",
        "attribute_used": col,
        "n_patients_with_subtype": int(len(wide)),
        "n_dm_thorsson_merged": int(len(merged)),
        "subtype_distribution": ct.to_dict(),
    }
    (TAB / "FIX2_summary_v2.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
