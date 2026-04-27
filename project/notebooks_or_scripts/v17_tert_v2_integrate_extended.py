"""Extended integration — apply 36 recovered TERT promoter calls to full
sample_master_v17_tert.tsv and produce all relevant cross-tabs.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from v17_tert_v2_common import V2_DIR, configure_logger

LOG = configure_logger("integrate_extended")

REPO = Path(__file__).resolve().parents[1].parent
SAMPLE_MASTER = REPO / "project" / "results" / "v17" / "tables" / "sample_master_v17_tert.tsv"

RE_SHORT = re.compile(r"(TCGA-[0-9A-Z]{2}-[0-9A-Z]{4})", re.I)


def short_barcode(s) -> str:
    if pd.isna(s):
        return ""
    m = RE_SHORT.search(str(s))
    return m.group(1).upper() if m else ""


def fisher_safe(table: list[list[int]]) -> dict:
    try:
        from scipy.stats import fisher_exact

        odds, p = fisher_exact(table)
        return {"table": table, "odds": float(odds), "p": float(p)}
    except Exception as e:  # noqa: BLE001
        return {"table": table, "error": str(e)}


def main():
    final = pd.read_csv(V2_DIR / "FINAL_tert_status_integrated.tsv", sep="\t")
    LOG.info("Recovered TERT-mut patients (TCGA): %d", len(final))
    mutated_set = set(final["tcga_short"].astype(str).str.upper())

    # Also load the raw promoter records to grab non-TCGA samples (MSK cohorts) for external validation
    raw = pd.read_csv(V2_DIR / "FINAL_promoter_records_all_sources.tsv", sep="\t")
    n_msk = (~raw["sample_barcode"].astype(str).str.startswith("TCGA-")).sum() if "sample_barcode" in raw.columns else 0
    LOG.info("Non-TCGA promoter records (external cohorts): %d", n_msk)

    sm = pd.read_csv(SAMPLE_MASTER, sep="\t")
    LOG.info("sample_master rows: %d", len(sm))
    sm["tcga_short"] = sm["sample_id"].apply(short_barcode)

    # Build new tert columns from recovered data
    sm["tert_promoter_v2"] = sm["tcga_short"].apply(
        lambda b: "mutated" if (b in mutated_set and b) else ""
    )
    # Preserve old column for comparison
    if "tert_promoter" in sm.columns:
        sm["tert_promoter_old"] = sm["tert_promoter"].astype(str)
    # Final integrated column: mutated if either old="mutation_verified" or new=mutated
    def integrated(row):
        old = str(row.get("tert_promoter_old", "")).lower()
        new = str(row.get("tert_promoter_v2", "")).lower()
        if new == "mutated" or "mutation_verified" in old:
            return "mutated"
        return "wildtype"

    sm["tert_promoter_integrated"] = sm.apply(integrated, axis=1)

    n_mut = (sm["tert_promoter_integrated"] == "mutated").sum()
    n_total = len(sm)
    LOG.info("Final mutated: %d / %d (%.1f%%)", n_mut, n_total, 100 * n_mut / n_total)

    # Save updated master
    out_master = V2_DIR / "sample_master_v17_tert_v2.tsv"
    sm.to_csv(out_master, sep="\t", index=False)
    LOG.info("Wrote %s", out_master)

    crosstabs = {}
    fishers = {}

    # 1. quad_group cross-tab
    if "quad_group" in sm.columns:
        ct = pd.crosstab(sm["quad_group"], sm["tert_promoter_integrated"])
        crosstabs["quad_group"] = ct.to_dict()
        out = V2_DIR / "FINAL_crosstab_quad_group.tsv"
        ct.to_csv(out, sep="\t")
        LOG.info("quad_group crosstab:\n%s", ct)

    # 2. v17_dark_cluster cross-tab (DM1 vs DM2)
    if "v17_dark_cluster" in sm.columns:
        ct = pd.crosstab(sm["v17_dark_cluster"], sm["tert_promoter_integrated"])
        crosstabs["v17_dark_cluster"] = ct.to_dict()
        out = V2_DIR / "FINAL_crosstab_dm_cluster.tsv"
        ct.to_csv(out, sep="\t")
        LOG.info("DM cluster crosstab:\n%s", ct)
        # Fisher between DM1 and DM2 if both exist
        if "DM1" in ct.index and "DM2" in ct.index and "mutated" in ct.columns and "wildtype" in ct.columns:
            a = ct.loc["DM1", "mutated"]
            b = ct.loc["DM1", "wildtype"]
            c = ct.loc["DM2", "mutated"]
            d = ct.loc["DM2", "wildtype"]
            fishers["DM1_vs_DM2"] = fisher_safe([[int(a), int(b)], [int(c), int(d)]])

    # 3. aggressive_flag cross-tab
    if "aggressive_flag" in sm.columns:
        ct = pd.crosstab(sm["aggressive_flag"].astype(str), sm["tert_promoter_integrated"])
        crosstabs["aggressive_flag"] = ct.to_dict()
        out = V2_DIR / "FINAL_crosstab_aggressive.tsv"
        ct.to_csv(out, sep="\t")

    # 4. molecular_subtype cross-tab
    if "molecular_subtype" in sm.columns:
        ct = pd.crosstab(sm["molecular_subtype"].astype(str), sm["tert_promoter_integrated"])
        crosstabs["molecular_subtype"] = ct.to_dict()
        out = V2_DIR / "FINAL_crosstab_molecular_subtype.tsv"
        ct.to_csv(out, sep="\t")

    # 5. driver_anchor (BRAF/RAS/etc) cross-tab
    if "driver_anchor" in sm.columns:
        ct = pd.crosstab(sm["driver_anchor"].astype(str), sm["tert_promoter_integrated"])
        crosstabs["driver_anchor"] = ct.to_dict()
        out = V2_DIR / "FINAL_crosstab_driver_anchor.tsv"
        ct.to_csv(out, sep="\t")

    # 6. Survival comparison (if os_event/os_days exist)
    surv = {}
    if {"os_event", "os_days"}.issubset(sm.columns):
        s2 = sm[(sm["os_event"].notna()) & (sm["os_days"].notna())].copy()
        if len(s2) > 30:
            try:
                from lifelines import KaplanMeierFitter
                from lifelines.statistics import logrank_test

                mut = s2[s2["tert_promoter_integrated"] == "mutated"]
                wt = s2[s2["tert_promoter_integrated"] == "wildtype"]
                if len(mut) >= 3 and len(wt) >= 3:
                    lr = logrank_test(mut["os_days"], wt["os_days"], mut["os_event"], wt["os_event"])
                    surv["TERT_logrank"] = {
                        "n_mutated": int(len(mut)),
                        "n_wildtype": int(len(wt)),
                        "events_mutated": int(mut["os_event"].sum()),
                        "events_wildtype": int(wt["os_event"].sum()),
                        "p_value": float(lr.p_value),
                        "test_statistic": float(lr.test_statistic),
                    }
                    LOG.info("Survival logrank TERT mut vs wt: p=%.4g", lr.p_value)
            except ImportError:
                surv["TERT_logrank"] = "lifelines not installed; manual KM analysis needed"
            except Exception as e:  # noqa: BLE001
                surv["TERT_logrank"] = f"error: {e}"

    # 4-group survival: BRAF/RAS/+TERT/triple-neg
    if "quad_group" in sm.columns and {"os_event", "os_days"}.issubset(sm.columns):
        s2 = sm[(sm["os_event"].notna()) & (sm["os_days"].notna())].copy()
        # Define new 4-group: incorporate TERT
        def four_group(r):
            base = str(r.get("quad_group", ""))
            tert = r.get("tert_promoter_integrated", "wildtype") == "mutated"
            if tert:
                return "TERT+"
            if "braf" in base.lower():
                return "BRAF_only"
            if "ras" in base.lower():
                return "RAS_only"
            return "Triple_neg"

        s2["four_group_v2"] = s2.apply(four_group, axis=1)
        ct4 = pd.crosstab(s2["four_group_v2"], s2["tert_promoter_integrated"])
        ct4.to_csv(V2_DIR / "FINAL_crosstab_4group.tsv", sep="\t")
        crosstabs["four_group_v2"] = ct4.to_dict()
        try:
            from lifelines.statistics import multivariate_logrank_test

            lr4 = multivariate_logrank_test(s2["os_days"], s2["four_group_v2"], s2["os_event"])
            surv["four_group_logrank"] = {
                "p_value": float(lr4.p_value),
                "test_statistic": float(lr4.test_statistic),
                "groups": s2["four_group_v2"].value_counts().to_dict(),
            }
            LOG.info("4-group logrank p=%.4g", lr4.p_value)
        except Exception as e:  # noqa: BLE001
            surv["four_group_logrank"] = f"error: {e}"

    summary = {
        "n_recovered_tcga_thca_tert_promoter": len(final),
        "n_total_sample_master": int(len(sm)),
        "n_mutated_in_sample_master": int(n_mut),
        "external_promoter_records_msk_etc": int(n_msk),
        "crosstabs": {k: {str(kk): {str(kkk): int(vv) for kkk, vv in v.items()} for kk, v in val.items()} for k, val in crosstabs.items()},
        "fishers": fishers,
        "survival": surv,
    }
    (V2_DIR / "FINAL_extended_summary.json").write_text(
        json.dumps(summary, indent=2, default=lambda o: int(o) if hasattr(o, "__int__") else str(o))
    )

    # Update audit
    audit = ["# v17 TERT promoter — extended integration", "", f"_Run: 2026-04-27_", ""]
    audit.append(f"- Recovered TERT-promoter mutated TCGA-THCA patients: **{len(final)}**")
    audit.append(f"- Sample master rows: {len(sm)}")
    audit.append(f"- Sample master rows now flagged TERT-mutated: {n_mut}")
    audit.append(f"- External (non-TCGA) TERT promoter records (MSK-IMPACT, MSK-CHORD, MSK PDTC/ATC): **{n_msk}**")
    audit.append("")
    audit.append("## Cross-tabs")
    for name, ct in crosstabs.items():
        audit.append(f"### {name}")
        audit.append("```")
        df_ct = pd.DataFrame(ct)
        audit.append(df_ct.to_string())
        audit.append("```")
        audit.append("")

    audit.append("## Fisher exact tests")
    audit.append("```json")
    audit.append(json.dumps(fishers, indent=2, default=str))
    audit.append("```")
    audit.append("")
    audit.append("## Survival")
    audit.append("```json")
    audit.append(json.dumps(surv, indent=2, default=str))
    audit.append("```")

    out = V2_DIR / "FINAL_extended_audit.md"
    out.write_text("\n".join(audit), encoding="utf-8")
    print("\n".join(audit))


if __name__ == "__main__":
    main()
