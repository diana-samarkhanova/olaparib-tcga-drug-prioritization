#!/usr/bin/env python3
"""Recover analysis-cohort purity from the archived matched clinical output.

This creates a derived input, not a replacement for the original GDC snapshot.
Run after restoring the archive and before rerunning analyze_tcga.py. It is only
intended to reproduce the archived cohort; obtain original data for new cohorts.
"""
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path,
                        default=ROOT / 'matched_clinical_cohort.csv')
    args = parser.parse_args()
    source = args.source.resolve()
    clinical = pd.read_csv(source)
    required = {'SAMPLE_ID', 'PURITY'}
    if not required.issubset(clinical.columns):
        raise ValueError('Archived clinical source must contain SAMPLE_ID and PURITY.')
    if clinical.SAMPLE_ID.isna().any() or not clinical.SAMPLE_ID.is_unique:
        raise ValueError('Archived clinical source must contain unique, nonmissing sample IDs.')
    if clinical.PURITY.notna().sum() == 0:
        raise ValueError('Archived clinical source contains no purity values to recover.')
    derived = clinical[['SAMPLE_ID', 'PURITY']].rename(
        columns={'SAMPLE_ID': 'Sample', 'PURITY': 'Purity'})
    raw = ROOT / 'raw'
    raw.mkdir(parents=True, exist_ok=True)
    destination = raw / 'purity_archive_derived.tsv'
    derived.to_csv(destination, sep='\t', index=False)
    provenance = {
        'purpose': 'Archive-cohort reproduction from a saved intermediate; not original GDC data',
        'source_file': source.name,
        'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'derived_file': destination.name,
        'derived_sha256': hashlib.sha256(destination.read_bytes()).hexdigest(),
        'sample_rows': len(derived),
        'nonmissing_purity': int(derived.Purity.notna().sum()),
        'transformation': 'Select SAMPLE_ID and PURITY; rename to Sample and Purity; retain missing values',
    }
    (raw / 'purity_archive_derived.provenance.json').write_text(
        json.dumps(provenance, indent=2) + '\n')
    print(json.dumps(provenance, indent=2))


if __name__ == '__main__':
    main()
