"""Export the supplied Chemical Data sheet without changing its rows or values."""

import argparse
from pathlib import Path

import pandas as pd


def main():
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--input', type=Path,
        default=root / 'upload/L1021-DiscoveryProbe-FDA-approved-Drug-Library (1).xlsx',
        help='Original library workbook containing the Chemical Data sheet.',
    )
    parser.add_argument(
        '--output', type=Path,
        default=Path(__file__).resolve().parent / 'library_raw.csv',
        help='CSV export used by the curated-alias stage.',
    )
    args = parser.parse_args()
    if not args.input.is_file():
        parser.error(f'Library workbook not found: {args.input}')
    catalog = pd.read_excel(args.input, sheet_name='Chemical Data')
    required = {'CatalogNumber', 'Item Name', 'CAS Number'}
    missing = required.difference(catalog.columns)
    if missing:
        parser.error(f'Missing Chemical Data columns: {sorted(missing)}')
    # Preserve all source rows, including the three rows without an Item Name
    # in the archived L1021 workbook. Downstream stages apply their own filters.
    args.output.parent.mkdir(parents=True, exist_ok=True)
    catalog.to_csv(args.output, index=False)
    print(f'Exported {len(catalog)} rows '
          f'({catalog["Item Name"].notna().sum()} named) to {args.output}')


if __name__ == '__main__':
    main()
