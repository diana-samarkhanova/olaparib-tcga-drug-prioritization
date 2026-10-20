"""Download original official PRISM release metadata and extract breast rows.

Uses only Python's standard library and pandas. No portal scraping.
The original full curve/matrix files are streamed and not stored in full.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import Request, urlopen
import csv
import io
import json
import pandas as pd

P = Path(__file__).resolve().parent
P.mkdir(parents=True, exist_ok=True)

FILES = {
    'primary_metadata.csv': 20237718,
    'secondary_metadata.csv': 20237769,
    'primary_readme.txt': 20237700,
    'secondary_readme.txt': 20238123,
    'primary_treatment_info.csv': 20237715,
    'secondary_treatment_info.csv': 20237763,
}

def open_url(url):
    return urlopen(Request(url, headers={'User-Agent':'PRISM-secondary-data-analysis/1.0'}),timeout=90)

def small_file(item):
    name, file_id = item
    with open_url(f'https://ndownloader.figshare.com/files/{file_id}') as response:
        (P / name).write_bytes(response.read())
    print('Downloaded',name,flush=True)

with ThreadPoolExecutor(max_workers=3) as pool:
    list(pool.map(small_file, FILES.items()))
with open_url('https://api.figshare.com/v2/articles/9393293') as response:
    (P / 'article_9393293.json').write_bytes(response.read())

secondary_metadata = pd.read_csv(P / 'secondary_metadata.csv')
breast = secondary_metadata[secondary_metadata.primary_tissue.eq('breast')]
breast.to_csv(P / 'secondary_breast_metadata.csv',index=False)
breast_ids = set(breast.depmap_id.dropna().astype(str))
curve_url = 'https://ndownloader.figshare.com/files/20237739'
selected = []
total = 0
with open_url(curve_url) as response:
    reader = csv.DictReader(io.TextIOWrapper(response,encoding='utf-8',newline=''))
    for record in reader:
        total += 1
        if record.get('depmap_id') in breast_ids:
            selected.append(record)
curves = pd.DataFrame(selected)
curves.to_csv(P / 'secondary_curve_parameters_breast.csv',index=False)
(P / 'curve_download_provenance.json').write_text(json.dumps({
    'source_url': curve_url,
    'landing_url': 'https://depmap.org/repurposing/',
    'article_id': 9393293,
    'records_total': total,
    'records_breast': len(curves),
    'breast_cell_line_ids': sorted(breast_ids),
},indent=2))

primary_metadata = pd.read_csv(P / 'primary_metadata.csv')
breast_barcodes = set(primary_metadata.loc[primary_metadata.primary_tissue.eq('breast'),'row_name'].dropna().astype(str))
primary_url = 'https://ndownloader.figshare.com/files/20237709'
rows = []
total_primary = 0
with open_url(primary_url) as response:
    reader = csv.reader(io.TextIOWrapper(response,encoding='utf-8',newline=''))
    header = next(reader)
    for record in reader:
        total_primary += 1
        if record[0] in breast_barcodes:
            rows.append(record)
wide = pd.DataFrame(rows,columns=header).rename(columns={header[0]:'row_name'})
wide.to_csv(P / 'primary_logfold_change_breast_wide.csv',index=False)
long = wide.melt(id_vars='row_name',var_name='column_name',value_name='log2_fold_change')
long['log2_fold_change'] = pd.to_numeric(long.log2_fold_change,errors='coerce')
long = long[long.log2_fold_change.notna()].copy()
long['viability_relative_DMSO'] = 2**long.log2_fold_change
treatment = pd.read_csv(P / 'primary_treatment_info.csv')
long = long.merge(treatment,on='column_name',how='left',validate='many_to_one')
long = long.merge(primary_metadata[['row_name','depmap_id','ccle_name','secondary_tissue']],on='row_name',how='left',validate='many_to_one')
long.to_csv(P / 'primary_logfold_change_breast_long.csv',index=False)
long[long.depmap_id.isin(['ACH-000849','ACH-000768','ACH-000573'])].to_csv(P / 'primary_logfold_change_468_231_436.csv',index=False)
(P / 'primary_download_provenance.json').write_text(json.dumps({
    'source_url': primary_url,
    'total_matrix_rows': total_primary,
    'selected_breast_rows': len(rows),
    'numeric_measurements_breast': len(long),
},indent=2))
print('Extracted secondary breast curves:',len(curves),'primary breast measurements:',len(long))
