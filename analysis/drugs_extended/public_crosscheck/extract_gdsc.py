"""Extract actual official GDSC release 8.5 fitted data without pooling datasets.

Reads the official source workbooks and the uploaded local library; writes only
compact CSVs. Name-only mappings are exact-normalized or provisional; no structure
identity is implied. Duplicate drug IDs remain separate for auditability.
"""
from pathlib import Path
import json
import hashlib
import re
import datetime
import numpy as np
import pandas as pd

P = Path(__file__).resolve().parent
ROOT = P.parents[2]

# The fitted workbooks are external inputs, not downloaded by this extractor.
required = [
    ROOT / 'analysis/prism/secondary_breast_metadata.csv',
    ROOT / 'analysis/prism/library_prism_name_coverage.csv',
    P / 'GDSC1_fitted_dose_response_27Oct23.xlsx',
]
gdsc2_file = P / 'GDSC2_fitted_dose_response_27Oct23.xlsx'
if not gdsc2_file.exists() and not (P / 'gdsc2_probe.xlsx').exists():
    required.append(gdsc2_file)
missing = [str(path) for path in required if not path.is_file()]
if missing:
    raise FileNotFoundError(
        'Missing GDSC extraction inputs:\n  ' + '\n  '.join(missing)
        + '\nThis script does not download the workbooks. See '
        'docs/DRUG_REPRODUCIBILITY.md for the release 8.5 URLs, input '
        'checksums, and prerequisite PRISM steps.'
    )

def name_key(s):
    return re.sub(r'[^a-z0-9]', '', str(s).lower())

def alias_key(s):
    s = str(s).strip().lower()
    s = re.sub(r'(?<=\w)\s+\([^)]*\)', '', s)
    salts = r'(?:hydrochloride|dihydrochloride|trihydrochloride|hcl|2hcl|hydrobromide|maleate|dimaleate|mesylate|methanesulfonate|citrate|fumarate|tartrate|ditartrate|sodium salt|potassium salt|monohydrate|dihydrate|trihydrate|hemihydrate|hydrate|free base)'
    while True:
        clean = re.sub(r'\s+' + salts + r'\s*$', '', s).strip()
        if clean == s:
            return name_key(s)
        s = clean

prism_meta = pd.read_csv(ROOT / 'analysis/prism/secondary_breast_metadata.csv')
prism_meta = prism_meta[prism_meta.ccle_name.notna()].copy()
prism_meta['cell_key'] = prism_meta.ccle_name.str.replace('_BREAST', '', regex=False).map(name_key)
prism_meta = prism_meta[['cell_key', 'depmap_id', 'ccle_name', 'secondary_tissue']].rename(columns={'secondary_tissue':'prism_subtype'})

library = pd.read_csv(ROOT / 'analysis/prism/library_prism_name_coverage.csv')
library = library[['library_excel_row','CatalogNumber','Item Name','CAS Number','Plate Location','Rack Number']].rename(columns={
    'Item Name':'library_name','CAS Number':'library_cas','Plate Location':'library_plate','Rack Number':'library_rack'})

def map_library(compounds):
    source = compounds.copy()
    source['exact_key'] = source.DRUG_NAME.map(name_key)
    source['alias_key'] = source.DRUG_NAME.map(alias_key)
    records = []
    for _, row in library.iterrows():
        m = source[source.exact_key.eq(name_key(row.library_name))]
        method = 'exact_normalized_name'
        if m.empty:
            m = source[source.alias_key.eq(alias_key(row.library_name))]
            method = 'provisional_salt_or_parenthetical_alias'
        if m.DRUG_NAME.nunique() > 1:
            method = 'ambiguous_name_candidates'
        for _, drug in m.iterrows():
            records.append({**row.to_dict(), 'match_method':method,
                            **drug.drop(['exact_key','alias_key']).to_dict()})
    return pd.DataFrame(records)

counts = {}
provenance = []
for dataset in ['GDSC1','GDSC2']:
    file = P / f'{dataset}_fitted_dose_response_27Oct23.xlsx'
    if not file.exists() and dataset == 'GDSC2':
        (P/'gdsc2_probe.xlsx').rename(file)
    print('Reading', file.name, flush=True)
    data = pd.read_excel(file)
    assert data.DATASET.eq(dataset).all()
    assert data.NLME_CURVE_ID.is_unique
    data['cell_key'] = data.CELL_LINE_NAME.map(name_key)
    breast = data[data.TCGA_DESC.eq('BRCA')].copy()
    breast = breast.merge(prism_meta, on='cell_key', how='left', validate='many_to_one')
    breast['ic50_fitted_uM'] = np.exp(breast.LN_IC50)
    # The source fits can extrapolate outside tested doses. We flag extrapolation
    # and retain the original fit; no zero/ceiling substitution is applied.
    breast['ic50_range_flag'] = np.select([
        ~np.isfinite(breast.ic50_fitted_uM),
        breast.ic50_fitted_uM.lt(breast.MIN_CONC),
        breast.ic50_fitted_uM.gt(breast.MAX_CONC)],
        ['missing_or_nonfinite','below_tested_range','above_tested_range'], default='within_tested_range')
    breast['ic50_in_tested_range_uM'] = breast.ic50_fitted_uM.where(breast.ic50_range_flag.eq('within_tested_range'))
    breast['source_release'] = '8.5; fitted workbook 27Oct2023'
    breast['assay_duration_hours'] = 72
    breast['AUC_orientation'] = 'lower AUC = greater loss of viability; dataset-native fraction'
    breast['IC50_unit'] = 'micromolar; exp(LN_IC50)'
    breast.to_csv(P/f'{dataset.lower()}_breast_curve_records.csv', index=False)
    breast[breast.cell_key.isin(['mdamb468','mdamb231','mdamb436'])].to_csv(P/f'{dataset.lower()}_468_231_436_curve_records.csv', index=False)
    drug_cols = ['DATASET','DRUG_ID','DRUG_NAME','PUTATIVE_TARGET','PATHWAY_NAME']
    compounds = data[drug_cols].drop_duplicates()
    compounds.to_csv(P/f'{dataset.lower()}_compound_metadata.csv', index=False)
    mapping = map_library(compounds)
    mapping.to_csv(P/f'library_{dataset.lower()}_name_mapping.csv', index=False)
    matched = mapping.merge(breast, on=['DATASET','DRUG_ID'], suffixes=('_mapping',''), how='inner')
    matched = matched.drop(columns=[x for x in matched if x.endswith('_mapping')])
    matched.to_csv(P/f'library_{dataset.lower()}_breast_records.csv', index=False)
    m468 = matched[matched.cell_key.eq('mdamb468')]
    m468.to_csv(P/f'library_{dataset.lower()}_mda_mb_468_records.csv', index=False)
    counts[dataset] = {
        'whole_release_curve_rows':len(data),
        'whole_release_models':int(data.SANGER_MODEL_ID.nunique()),
        'whole_release_drug_ids':int(data.DRUG_ID.nunique()),
        'whole_release_drug_names':int(data.DRUG_NAME.nunique()),
        'breast_curve_rows':len(breast),
        'breast_models':int(breast.SANGER_MODEL_ID.nunique()),
        'breast_models_with_prism_TNBC_annotation':int(breast.loc[breast.prism_subtype.eq('TNBC'),'SANGER_MODEL_ID'].nunique()),
        'mda_mb_468_curve_rows':int(breast.cell_key.eq('mdamb468').sum()),
        'mda_mb_468_drug_names':int(breast.loc[breast.cell_key.eq('mdamb468'),'DRUG_NAME'].nunique()),
        'library_compounds_metadata_matched':int(mapping.library_excel_row.nunique()),
        'library_compounds_metadata_exact_matched':int(mapping.loc[mapping.match_method.eq('exact_normalized_name'),'library_excel_row'].nunique()),
        'library_compounds_468_matched':int(m468.library_excel_row.nunique()),
        'library_compounds_468_exact_matched':int(m468.loc[m468.match_method.eq('exact_normalized_name'),'library_excel_row'].nunique()),
        '468_ic50_in_range_records':int(breast.loc[breast.cell_key.eq('mdamb468'),'ic50_range_flag'].eq('within_tested_range').sum()),
        '468_rmse_above_0p3_records':int(breast.loc[breast.cell_key.eq('mdamb468'),'RMSE'].gt(0.3).sum())}
    url=f'https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/{file.name}'
    provenance.append({'dataset':dataset,'source_url':url,'source_release':'8.5','source_filename':file.name,
                       'download_date_utc':datetime.datetime.now(datetime.timezone.utc).date().isoformat(),
                       'bytes':file.stat().st_size,'sha256':hashlib.sha256(file.read_bytes()).hexdigest()})
    print(json.dumps(counts[dataset]),flush=True)

(P/'gdsc_coverage_summary.json').write_text(json.dumps(counts,indent=2)+'\n')
(P/'gdsc_download_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
