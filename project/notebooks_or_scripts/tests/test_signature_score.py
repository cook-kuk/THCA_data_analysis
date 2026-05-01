#!/usr/bin/env python3
"""Smoke test — signature score reproducibility."""
import json
from pathlib import Path
import numpy as np
import pandas as pd

PROJ = Path("/opt/thyroid-dash/project")


def test_pillar1_korean_pool_size():
    """Korean PTC pool n=874 (235 + 630 + 9)."""
    pool = pd.read_csv(PROJ / "results/d4p1_panasian_meta/korean_PTC_pool_n908.tsv", sep="\t")
    assert len(pool) == 874, f"Expected n=874, got {len(pool)}"
    cohort_counts = pool["cohort"].value_counts().to_dict()
    assert cohort_counts.get("K2") == 235, f"K2 expected 235, got {cohort_counts.get('K2')}"
    assert cohort_counts.get("Lee2024") == 630, f"Lee expected 630, got {cohort_counts.get('Lee2024')}"
    assert cohort_counts.get("GSE286332_PTC") == 9, f"GSE286332-PTC expected 9, got {cohort_counts.get('GSE286332_PTC')}"


def test_pillar1_dpb1_freq():
    """DPB1*05:01 carrier freq 53.2% (465/874) Korean PTC pool."""
    pool = pd.read_csv(PROJ / "results/d4p1_panasian_meta/korean_PTC_pool_n908.tsv", sep="\t")
    a1 = pool["DPB1_a1_4d"].fillna("")
    a2 = pool["DPB1_a2_4d"].fillna("")
    carriers = ((a1 == "DPB1*05:01") | (a2 == "DPB1*05:01")).sum()
    assert carriers == 465, f"Expected 465 DPB1*05:01 carriers, got {carriers}"
    freq = carriers / len(pool)
    assert abs(freq - 0.532) < 0.001, f"Expected freq 0.532, got {freq:.3f}"


def test_pillar2_DEG_count():
    """GSE286332 PyDESeq2 DEGs at padj<0.05 should be 10,380."""
    deg = pd.read_csv(PROJ / "results/p3_gse286332/deg_ptcht_vs_ptc.tsv", sep="\t")
    n_sig = ((deg["padj"] < 0.05)).sum()
    assert n_sig == 10380, f"Expected 10,380 DEGs, got {n_sig}"


def test_pillar2_8gene_per_gene():
    """8-gene per-gene Cohen d direction (PAX8 < TPO < SLC5A5)."""
    df = pd.read_csv(PROJ / "results/p3_gse286332/8gene_per_gene_compare.tsv", sep="\t")
    df_idx = df.set_index("gene")
    # Direction: most strongly down (PAX8 d=-2.32) < TPO (-1.11) < SLC5A5 (+0.25)
    assert df_idx.loc["PAX8", "cohen_d"] < -2.0
    assert df_idx.loc["NKX2-1", "cohen_d"] < -1.5
    assert df_idx.loc["FOXE1", "cohen_d"] < -1.5
    assert df_idx.loc["SLC5A5", "cohen_d"] > 0  # preserved iodide transporter


def test_pillar3_braf_neutrality():
    """BRAF mRNA × V600E status Cohen d ~ -0.04 (negligible)."""
    df = pd.read_csv(PROJ / "results/p1_driver_mrna_audit/driver_mrna_mutation_audit.tsv", sep="\t")
    braf_row = df[df["gene"] == "BRAF"]
    assert len(braf_row) == 1
    assert abs(braf_row.iloc[0]["cohen_d"] - (-0.044)) < 0.01


def test_pillar4_pangenome_top5000_ARI():
    """Pan-genome top-5000 MAD ARI should be 0.918."""
    ari = pd.read_csv(PROJ / "results/p4_pangenome_vs_tiera67/ari_comparison.tsv", sep="\t")
    top5000 = ari[ari["label"] == "Pan-genome top 5000 MAD"]
    assert len(top5000) == 1
    assert abs(top5000.iloc[0]["ARI"] - 0.918) < 0.01


def test_pillar4_driver_only_ARI_zero():
    """Driver_anchor 12-only cluster ARI ~ -0.007 (random)."""
    ari = pd.read_csv(PROJ / "results/p4_pangenome_vs_tiera67/ari_comparison.tsv", sep="\t")
    drivers = ari[ari["label"] == "Driver_anchor 12"]
    assert len(drivers) == 1
    assert abs(drivers.iloc[0]["ARI"]) < 0.05


def test_pillar5_mediation_HLA2_140pct():
    """HLA-II mediation %=140% (over-mediation/full-pathway)."""
    med = json.loads((PROJ / "results/d3p5_pdm1_gradient/mediation_results.json").read_text())
    assert "HLA_II" in med
    assert abs(med["HLA_II"]["pct_mediated"] - 140.0) < 5.0  # 140% ± 5
    assert med["HLA_II"]["boot_p_emp"] < 0.05


def test_pillar5_TCGA_hashimoto_top30_OR():
    """TCGA Hashimoto-like top-30% threshold OR=0.20, p=6e-10."""
    summary = json.loads((PROJ / "results/d4p2_tcga_hashimoto_signature/D4P2_summary.json").read_text())
    top30 = summary["crosstabs"]["hashi_top30"]
    assert abs(float(top30["OR"]) - 0.20) < 0.01
    assert float(top30["fisher_p"]) < 1e-9


def test_pillar5_korean_replication():
    """Korean GSE213647 Hashimoto-like: GMM 22.8%, Otsu 28.2%."""
    summary = json.loads((PROJ / "results/d8b_korean_replication/D8B_summary.json").read_text())
    assert abs(summary["hashi_pct"]["GMM"] - 22.8) < 1.0
    assert abs(summary["hashi_pct"]["Otsu"] - 28.2) < 1.0


def test_chu2018_published_OR():
    """Chu 2018 published values: DPB1*05:01 OR=1.90, B*46:01 OR=2.38."""
    chu = pd.read_csv(PROJ / "results/p2_pillar1_forest/chu2018_allele_summary.tsv", sep="\t")
    chu_idx = chu.set_index("allele")
    assert abs(chu_idx.loc["DPB1*05:01", "OR"] - 1.90) < 0.01
    assert abs(chu_idx.loc["B*46:01", "OR"] - 2.38) < 0.01
    # 4 alleles direction-consistent up + 2 protective down
    assert (chu["OR"] > 1.0).sum() >= 5
    assert (chu["OR"] < 1.0).sum() >= 3


if __name__ == "__main__":
    print("=== Running smoke tests ===")
    test_pillar1_korean_pool_size()
    print("✓ Pillar 1 Korean pool size n=874")
    test_pillar1_dpb1_freq()
    print("✓ Pillar 1 DPB1*05:01 freq 53.2%")
    test_pillar2_DEG_count()
    print("✓ Pillar 2 DEG count 10,380")
    test_pillar2_8gene_per_gene()
    print("✓ Pillar 2 8-gene direction")
    test_pillar3_braf_neutrality()
    print("✓ Pillar 3 BRAF neutrality d=-0.04")
    test_pillar4_pangenome_top5000_ARI()
    print("✓ Pillar 4 pan-genome top-5000 ARI=0.92")
    test_pillar4_driver_only_ARI_zero()
    print("✓ Pillar 4 driver-only ARI~0")
    test_pillar5_mediation_HLA2_140pct()
    print("✓ Pillar 5 HLA-II 140% mediation")
    test_pillar5_TCGA_hashimoto_top30_OR()
    print("✓ Pillar 5 TCGA Hashimoto top30 OR=0.20")
    test_pillar5_korean_replication()
    print("✓ Pillar 5 Korean replication GMM=22.8%")
    test_chu2018_published_OR()
    print("✓ Chu 2018 published OR values")
    print("\n✅ ALL SMOKE TESTS PASSED")
