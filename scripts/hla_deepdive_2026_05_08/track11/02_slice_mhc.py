#!/usr/bin/env python3
"""Track 11 step 2 — slice MHC (chr6:28-34Mb) per endpoint via remote tabix streaming.

Inputs:  T01_endpoint_manifest.tsv
Outputs: cache/<source>__<phenocode>__<ancestry>.tsv.gz   (per-endpoint MHC slice in canonical 10 cols)
         tables/T02_mhc_slice_summary.tsv                 (per-endpoint count + status)
"""
from __future__ import annotations
import csv
import gzip
import os
import re
import sys
from pathlib import Path

os.environ['CURL_CA_BUNDLE'] = '/etc/ssl/certs/ca-certificates.crt'
import pysam  # noqa: E402

OUT_DIR = Path('/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track11_finngen_panukbb')
CACHE = OUT_DIR / 'cache'
TABLES = OUT_DIR / 'tables'
CACHE.mkdir(parents=True, exist_ok=True)
TABLES.mkdir(parents=True, exist_ok=True)

MHC_GRCh38 = ('6', 28510120, 33480577)
MHC_GRCh37 = ('6', 28477797, 33448354)

# Canonical output columns (per-variant slice).
CANON_COLS = ['chrom', 'pos', 'ref', 'alt', 'rsid', 'gene', 'beta', 'se', 'pval', 'af', 'n']


def read_bgz_first_line(url: str) -> str:
    """Read first line of a remote bgz file — gives header for sumstats with single-line headers."""
    bgz = pysam.BGZFile(url, 'r')
    return bgz.readline().decode()


def fetch_finngen(url: str, region) -> list[dict]:
    """FinnGen R12 columns: #chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt af_alt_cases af_alt_controls"""
    chrom, start, stop = region
    tbx = pysam.TabixFile(url)
    out = []
    # FinnGen tabix uses '6' (not 'chr6')
    try:
        rows = tbx.fetch(chrom, start, stop)
    except Exception as e:
        # try 'chr6'
        rows = tbx.fetch('chr' + chrom, start, stop)
    for r in rows:
        f = r.split('\t')
        if len(f) < 11:
            continue
        try:
            beta = float(f[8]) if f[8] not in ('NA', '') else None
            se   = float(f[9]) if f[9] not in ('NA', '') else None
            p    = float(f[6]) if f[6] not in ('NA', '') else None
            af   = float(f[10]) if f[10] not in ('NA', '') else None
        except ValueError:
            continue
        if p is None:
            continue
        out.append({
            'chrom': f[0], 'pos': int(f[1]),
            'ref': f[2], 'alt': f[3],
            'rsid': f[4] or '.',
            'gene': f[5] or '.',
            'beta': beta, 'se': se, 'pval': p, 'af': af,
            'n': '',  # FinnGen embeds case/control counts separately at endpoint level
        })
    return out


def fetch_panukbb(url: str, region, ancestry: str) -> list[dict]:
    """Pan-UKBB columns: chr pos ref alt af_cases_<A> af_controls_<A> beta_<A> se_<A> neglog10_pval_<A> low_confidence_<A> [...repeats per ancestry]
    Header line tells us where the ancestry block lives."""
    chrom, start, stop = region
    header = read_bgz_first_line(url).rstrip('\n').split('\t')
    # Column indices for this ancestry
    def ci(col):
        try:
            return header.index(col)
        except ValueError:
            return None
    ix = {
        'af_case': ci(f'af_cases_{ancestry}'),
        'af_ctrl': ci(f'af_controls_{ancestry}'),
        'beta':    ci(f'beta_{ancestry}'),
        'se':      ci(f'se_{ancestry}'),
        'mlog':    ci(f'neglog10_pval_{ancestry}'),
        'lowconf': ci(f'low_confidence_{ancestry}'),
    }
    if ix['beta'] is None or ix['mlog'] is None:
        return []
    tbx = pysam.TabixFile(url)
    out = []
    rows = tbx.fetch(chrom, start, stop)
    for r in rows:
        f = r.split('\t')
        try:
            mlog = float(f[ix['mlog']]) if f[ix['mlog']] not in ('NA','') else None
            beta = float(f[ix['beta']]) if f[ix['beta']] not in ('NA','') else None
            se   = float(f[ix['se']])   if f[ix['se']]   not in ('NA','') else None
        except (ValueError, IndexError):
            continue
        if mlog is None or beta is None:
            continue
        # convert -log10(p) -> p
        try:
            p = 10 ** (-mlog)
        except OverflowError:
            p = 0.0
        af = None
        if ix['af_case'] is not None and ix['af_ctrl'] is not None:
            try:
                af_case = float(f[ix['af_case']]) if f[ix['af_case']] not in ('NA','') else None
                af_ctrl = float(f[ix['af_ctrl']]) if f[ix['af_ctrl']] not in ('NA','') else None
                if af_case is not None and af_ctrl is not None:
                    af = af_ctrl  # control AF as proxy
            except ValueError:
                pass
        # rsid not in Pan-UKBB flat — leave '.' (we'll annotate downstream from Pan-UKBB variant manifest if needed)
        out.append({
            'chrom': f[0], 'pos': int(f[1]),
            'ref': f[2], 'alt': f[3],
            'rsid': '.',
            'gene': '.',
            'beta': beta, 'se': se, 'pval': p, 'af': af,
            'n': '',
        })
    return out


def write_slice(rows, out_path):
    with gzip.open(out_path, 'wt') as f:
        w = csv.DictWriter(f, fieldnames=CANON_COLS, delimiter='\t')
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    manifest = TABLES / 'T01_endpoint_manifest.tsv'
    summary_rows = []
    with open(manifest) as f:
        endpoints = list(csv.DictReader(f, delimiter='\t'))
    print(f'Slicing {len(endpoints)} endpoints', flush=True)
    for i, ep in enumerate(endpoints):
        tag = ep['tag']
        url = ep['url']
        if ep['source'] == 'PanUKBB':
            url = url.replace('s3://pan-ukb-us-east-1/', 'https://pan-ukb-us-east-1.s3.amazonaws.com/')
        out_path = CACHE / f"{ep['source']}__{ep['phenocode'].replace(':','_').replace('/','_')}__{ep['ancestry']}.tsv.gz"
        print(f'[{i+1}/{len(endpoints)}] {ep["source"]} {tag} ({ep["ancestry"]}) -> {out_path.name}', flush=True)
        if out_path.exists() and out_path.stat().st_size > 200:
            # already done
            n_rows = sum(1 for _ in gzip.open(out_path, 'rt')) - 1
            summary_rows.append({**ep, 'n_variants': n_rows, 'status': 'cached', 'slice_path': str(out_path)})
            print(f'    cached, {n_rows} variants', flush=True)
            continue
        try:
            if ep['source'] == 'FinnGen_R12':
                region = MHC_GRCh38
                rows = fetch_finngen(url, region)
            else:
                region = MHC_GRCh37  # Pan-UKBB sumstats are GRCh37
                rows = fetch_panukbb(url, region, ep['ancestry'])
            write_slice(rows, out_path)
            summary_rows.append({**ep, 'n_variants': len(rows), 'status': 'ok', 'slice_path': str(out_path)})
            print(f'    OK, {len(rows)} variants', flush=True)
        except Exception as e:
            err = f'{type(e).__name__}: {e}'
            summary_rows.append({**ep, 'n_variants': 0, 'status': f'fail: {err[:120]}', 'slice_path': ''})
            print(f'    FAIL {err}', flush=True)

    out_tsv = TABLES / 'T02_mhc_slice_summary.tsv'
    fields = list(endpoints[0].keys()) + ['n_variants', 'status', 'slice_path']
    with open(out_tsv, 'w') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        for r in summary_rows:
            w.writerow(r)
    n_ok = sum(1 for r in summary_rows if r['status'] in ('ok', 'cached'))
    print(f'Done: {n_ok}/{len(summary_rows)} sliced -> {out_tsv}')


if __name__ == '__main__':
    main()
