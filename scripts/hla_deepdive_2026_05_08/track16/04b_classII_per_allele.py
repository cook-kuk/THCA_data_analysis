#!/usr/bin/env python3
"""Single-allele HLA-II prediction subprocess.
Usage: 04b_classII_per_allele.py <allele> <pep_file> <out_csv>
"""
import sys, os
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
import warnings; warnings.filterwarnings("ignore")
from mhcnuggets.src.predict import predict as mhcn_predict


def to_token(allele: str) -> str:
    return "HLA-" + allele.replace("*", "")


def main():
    allele, pep_file, out_csv = sys.argv[1], sys.argv[2], sys.argv[3]
    mhcn_predict(class_="II",
                 peptides_path=pep_file,
                 mhc=to_token(allele),
                 output=out_csv,
                 rank_output=True)


if __name__ == "__main__":
    main()
