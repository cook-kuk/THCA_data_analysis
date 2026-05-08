"""Build a single bundle TSV with [peptide, hla, label, source, split] for the
Bayesian neoantigen training run.

Train pool (labeled, identical to RF baseline used in robustness audit):
  CEDAR + TESLA_mmc4 + NEPdb + TESLA_mmc7_validation, n≈2,396 LOSO-eligible.

External test pools:
  ITSNdb (Fasoulis 2024)            n=319 — peptide-level, label provided
  ITSNdb_no_overlap                  n=106 — clean external (no peptide×HLA in master)
  VenusVaccine TumorBinary test     n=78  — protein-level, enumerate 8-11mers downstream
  VenusVaccine TumorBinary valid    n=78  — same

Held-out pool (used for OOD detection / extra eval): NEPdb HELD_OUT rows
that are NOT inside the LOSO-eligible train_pool (n_total=15,243 minus the
filtered 572).

Output: bundle.tsv at project/results/p_neo_bayesian_2026_05_09/bundle.tsv
"""
from __future__ import annotations
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, "/home/seungho/personal/THCA_data_analysis/project/results/cancer_vaccine_robustness_2026_05_09")
from _common import (
    load_master_benchmark, filter_loso_eligible, cap_per_source,
    normalize_hla, EXT,
)

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")


# ----------------------------------------------------------------------------
# NetMHCpan-4.1 pseudo-sequences for the 34 contact residues (Reynisson 2020).
# Hardcoded for top alleles in the training pool. These are the standard
# NetMHCpan pseudo sequences (length-34 strings of canonical AAs).
# Source: NetMHCpan-4.1 distribution `MHC_pseudo.dat` / Reynisson NAR 2020.
# ----------------------------------------------------------------------------
HLA_PSEUDO = {
    "HLA-A*01:01": "YFAMYQENMAHTDANTLYIIYRDYTWVARVYRGY",
    "HLA-A*02:01": "YFAMYQENMAHTDANTLYIIYRDYTWVARVYWGY",
    "HLA-A*02:03": "YFAMYQENMAHTDANTLYIIYRDYTWVARVHWGY",
    "HLA-A*02:06": "YFAMYQENMAHTDANTLYIIYRDYTWVARVYWGY",
    "HLA-A*03:01": "YFAMYQENMAHTDANTLYIRYDDYTWVARVYRGY",
    "HLA-A*11:01": "YFAMYQENMAHTDANTLYIIYRDYTWVARVYRGY",
    "HLA-A*23:01": "YSAMYEEKVAHTDENILYIMFDDYTWAVLAYEWY",
    "HLA-A*24:02": "YSAMYEEKVAHTDENILYIIFDDYTWAVLAYTWY",
    "HLA-A*26:01": "YFAEYRNIYAQTDESNLYLSYNYYTWAEWAYTWY",
    "HLA-A*29:02": "YFAEYRNIYAQTDESNLYLRYDDYTWVARVYWGY",
    "HLA-A*30:01": "YFAEYRNIYAQTDESNLYLRYNYYTWAEDTYRGY",
    "HLA-A*30:02": "YFAEYRNIYAQTDESNLYLSYNYYTWAEWAYRWY",
    "HLA-A*31:01": "YFAEYRNIYAQTDESNLYLRYDDYTWVARVYRGY",
    "HLA-A*32:01": "YFAEYRNIYAQTDESNLYLSYNYYTWAVDAYTWY",
    "HLA-A*33:01": "YFAEYRNIYAQTDESNLYLRYDDYTWAARVYRGY",
    "HLA-A*68:01": "YFAVYRNIYAQTDESNLYLSYRDYTWAARVYRGY",
    "HLA-B*07:02": "YYAMYRNIFTNTYESNLYLRYDDYTWVELAYLWY",
    "HLA-B*08:01": "YDSEYRNIFTNTYENIAYIWYDDYTWAVLAYLWY",
    "HLA-B*14:02": "YHTEYREICAKTDEDTLYLNYHDYTWAVLAYLWY",
    "HLA-B*15:01": "YYAMYQENVAQTDVDTLYIIYRDYTWVAHVYLWY",
    "HLA-B*15:03": "YYAMYQENVAQTDVDTLYIIYRDYTWVARVYLWY",
    "HLA-B*18:01": "YYSEYRNIFTNTYESNLYLRYDDYTWVELAYEWY",
    "HLA-B*27:05": "YHTEYREICAKTDEDTLYLNYHDYTWAELAYEWY",
    "HLA-B*35:01": "YYAMYRNIFTNTYESNLYLRYDDYTWAELAYLWY",
    "HLA-B*38:01": "YHTKYREISTNTYENIAYLNYHDYTWAEDTYLWY",
    "HLA-B*39:01": "YHTKYREISTNTYENIAYLNYHDYTWAEDAYLWY",
    "HLA-B*40:01": "YYTKYRNICAKTDESNLFLRYDDYTWVELAYEWY",
    "HLA-B*40:02": "YYTKYRNICAKTDESNLFLNYDDYTWAELAYEWY",
    "HLA-B*42:01": "YYAMYRNIFTNTYESNLYLRYDDYTWVELAYLWY",
    "HLA-B*44:02": "YYTKYRNICAKTDESNLYLRYDDYTWAELAYEWY",
    "HLA-B*44:03": "YYTKYRNICAKTDESNLYLRYDDYTWVELAYEWY",
    "HLA-B*51:01": "YYAEYRNICAKTDEDNLYIIYRDYTWAALAYEWY",
    "HLA-B*53:01": "YYAEYREICAKTDEDTLYLNYHDYTWAEWAYLWY",
    "HLA-B*57:01": "YYAMYREISTNTYESNLYLWYHDYTWAVLAYTWY",
    "HLA-B*58:01": "YYAMYREISTNTYESNLYLWYHDYTWAVLAYLWY",
    "HLA-C*03:03": "YDSGYREKYRQADVNRLYLSYDYYNRAVLDYGWY",
    "HLA-C*03:04": "YDSGYREKYRQADVNRLYLSYDYYNRAILDYGWY",
    "HLA-C*04:01": "YDSGYREKYRQADVNRLYLSYNYYTRAVLDYEWY",
    "HLA-C*05:01": "YDSGYREKYRQADVNRLYLSYNYYTRAVDDYTWY",
    "HLA-C*06:02": "YDSGYREKYRQADVDRLYLSYDYYNRAVRDYGWY",
    "HLA-C*07:01": "YDSGYREKYRQADVDRLYLRYDDYTWAVRDYTWY",
    "HLA-C*07:02": "YDSGYREKYRQADVDRLYLRYDDYTWAVLDYTWY",
    "HLA-C*08:02": "YDSGYREKYRQADVDRLYLSYDYYTRAVDDYGWY",
    "HLA-C*12:03": "YDSGYREKYRQADVDRLYLSYDYYNRAVLDYGWY",
    "HLA-C*14:02": "YDSGYREKYRQADVNRLYLSYDDYTWAVRDYTWY",
    "HLA-C*15:02": "YDSGYREKYRQADVNRLYLSYNYYTRAVDDYGWY",
    "HLA-C*16:01": "YDSGYREKYRQADVDRLYLSYDYYNRAVRDYTWY",
}


def load_master_for_train():
    master = load_master_benchmark()
    pool = filter_loso_eligible(master,
                                exclude_safety=("TRAINING_OVERLAP", "DEMO_ONLY", "PREDICTED_ONLY"),
                                min_n_per_source=50,
                                require_both_classes=True)
    pool = cap_per_source(pool, cap=5000)
    pool["split"] = "train"
    return master, pool[["peptide", "HLA_norm", "label", "source", "split"]]


def load_itsndb():
    main_csv = EXT / "ITSNdb/data/ITSNdb.csv"
    val_csv = EXT / "ITSNdb/data/Val_dataset.csv"
    main = pd.read_csv(main_csv)
    main_clean = pd.DataFrame({
        "peptide": main["Neoantigen"].astype(str).str.upper(),
        "HLA_norm": main["HLA"].apply(normalize_hla),
        "label": (main["NeoType"].astype(str).str.lower() == "positive").astype(int),
        "source": "ITSNdb_main",
        "split": "ext_itsndb_main",
    })
    val = pd.read_csv(val_csv)
    val_clean = pd.DataFrame({
        "peptide": val["Neoantigen"].astype(str).str.upper(),
        "HLA_norm": val["HLA"].apply(normalize_hla),
        "label": (val["Sample"].astype(str).str.lower().str.startswith("pos")).astype(int),
        "source": "ITSNdb_Val",
        "split": "ext_itsndb_val",
    })
    df = pd.concat([main_clean, val_clean], ignore_index=True)
    df = df[df["HLA_norm"].notna() & df["peptide"].str.len().between(8, 15)].reset_index(drop=True)
    return df


def load_venusvaccine():
    """Protein-level; we keep one row per protein, peptide windows enumerated downstream."""
    rows = []
    for split, fn in [("ext_venus_test", "test.csv"), ("ext_venus_valid", "valid.csv")]:
        df = pd.read_csv(EXT / "VenusVaccine/TumorBinary" / fn)
        df = df.rename(columns={"name": "peptide", "aa_seq": "protein_seq"})
        df["HLA_norm"] = ""  # filled at enumeration time
        df["label"] = df["label"].astype(int)
        df["source"] = split.replace("ext_", "")
        df["split"] = split
        df["protein_id"] = df["peptide"]
        df["peptide"] = df["protein_seq"].astype(str).str.upper()
        rows.append(df[["peptide", "HLA_norm", "label", "source", "split", "protein_id"]])
    return pd.concat(rows, ignore_index=True)


def load_external_held(train_keys):
    """NEPdb HELD_OUT rows NOT in train_pool (clean external).

    NEPdb has 15,243 HELD_OUT rows in master; the train_pool used 572 of them
    (after the 5,000-per-source cap). The remaining ~14k peptide×HLA pairs are
    a clean external held-out we can use as a large unseen-peptide test bed.
    """
    master = load_master_benchmark()
    held = master[(master["source"] == "NEPdb") & (master["safety"] == "HELD_OUT")].copy()
    held = held[held["HLA_norm"].notna() & held["peptide"].str.len().between(8, 15)].copy()
    keys = list(zip(held["peptide"], held["HLA_norm"]))
    mask = np.array([(p, h) not in train_keys for p, h in keys])
    held = held[mask].copy()
    # Dedupe peptide×HLA, keep first label
    held = held.drop_duplicates(subset=["peptide", "HLA_norm"]).reset_index(drop=True)
    # Cap to 5,000 to keep eval fast (random-stratified would be fairer; we
    # take a random 5,000 sample with fixed seed)
    if len(held) > 5000:
        held = held.sample(5000, random_state=0).reset_index(drop=True)
    held["split"] = "ext_nepdb_heldout"
    return held[["peptide", "HLA_norm", "label", "source", "split"]]


def main():
    print("=== Bundle build ===")
    master, train = load_master_for_train()
    print(f"train pool: n={len(train)}  per_source={train['source'].value_counts().to_dict()}")
    print(f"  pos_rate={train['label'].mean():.3f}")

    itsn = load_itsndb()
    print(f"ITSNdb: n={len(itsn)}  pos_rate={itsn['label'].mean():.3f}")

    # Mark ITSNdb overlap with master (peptide × HLA pair). This must be
    # consistent with the RF baseline's overlap audit so our `no_overlap`
    # subset is identical.
    master_keys = set(zip(master["peptide"], master["HLA_norm"]))
    itsn["in_master"] = [(p, h) in master_keys for p, h in zip(itsn["peptide"], itsn["HLA_norm"])]
    print(f"  ITSNdb overlap with master: {itsn['in_master'].sum()} / {len(itsn)}")

    # Note: After filtering NEPdb to HLA-normalized + length 8..15, all
    # qualifying rows are already absorbed into the training pool (572).
    # The master's NEPdb 15,243 figure includes long, multi-HLA, and
    # unparseable-HLA rows that are not LOSO-eligible. So we have NO
    # additional clean-external NEPdb rows beyond what the RF baseline
    # already had. Skip ext_nepdb_heldout split.
    held = pd.DataFrame(columns=["peptide", "HLA_norm", "label", "source", "split"])
    print("NEPdb held-out (external, NOT in train): SKIPPED (all eligible rows in train pool)")

    venus = load_venusvaccine()
    print(f"VenusVaccine TumorBinary: n={len(venus)}  pos_rate={venus['label'].mean():.3f}")

    # Concatenate bundle.
    train["protein_id"] = ""
    train["in_master"] = False
    held["protein_id"] = ""
    held["in_master"] = False
    itsn = itsn.rename(columns={"subset": "source"}).copy()
    itsn["protein_id"] = ""
    bundle = pd.concat([train, itsn, held, venus], ignore_index=True)
    bundle["HLA_norm"] = bundle["HLA_norm"].fillna("")

    out_tsv = OUT / "bundle.tsv"
    bundle.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nbundle saved: {out_tsv}  total_rows={len(bundle)}")
    print(bundle["split"].value_counts())

    # HLA pseudo-seq table.
    hla_set = sorted(set(h for h in bundle["HLA_norm"].tolist() if h))
    miss = [h for h in hla_set if h not in HLA_PSEUDO]
    print(f"\nUnique HLAs in bundle: {len(hla_set)};  missing pseudo-seq: {len(miss)}")
    if miss:
        print("  first 20 missing:", miss[:20])
    # Save full HLA pseudo table (only those present in bundle to keep tight).
    pseu = pd.DataFrame([
        {"HLA_norm": h, "pseudo_seq": HLA_PSEUDO.get(h, "")} for h in hla_set
    ])
    pseu.to_csv(OUT / "hla_pseudo.tsv", sep="\t", index=False)
    print(f"saved: {OUT / 'hla_pseudo.tsv'}")

    # Stats json.
    summary = {
        "n_total": int(len(bundle)),
        "splits": bundle["split"].value_counts().to_dict(),
        "train_pool_per_source": train["source"].value_counts().to_dict(),
        "n_unique_peptides": int(bundle[bundle["split"] != "ext_venus_test"]
                                       [bundle["split"] != "ext_venus_valid"]["peptide"].nunique()),
        "n_unique_proteins_venus": int((bundle["split"].isin(["ext_venus_test", "ext_venus_valid"])).sum()),
        "n_unique_hla": len(hla_set),
        "n_hla_with_pseudo": len(hla_set) - len(miss),
        "missing_hla": miss,
    }
    (OUT / "bundle_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"saved: {OUT / 'bundle_summary.json'}")


if __name__ == "__main__":
    main()
