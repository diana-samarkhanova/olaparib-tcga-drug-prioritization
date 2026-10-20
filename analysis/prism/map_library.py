from pathlib import Path
import json
import re
import numpy as np
import pandas as pd

P = Path(__file__).resolve().parent
ROOT = P.parents[1]
library = pd.read_excel(ROOT / 'upload/L1021-DiscoveryProbe-FDA-approved-Drug-Library (1).xlsx', sheet_name='Chemical Data')
library = library[library['Item Name'].notna()].copy()
library['library_excel_row'] = library.index + 2

def key(s):
    return re.sub(r'[^a-z0-9]', '', str(s).lower())

def alias_key(s):
    s = str(s).strip().lower()
    # Retain leading stereochemical descriptors. Remove trailing parenthetical
    # catalog/brand aliases only as a provisional name-based match.
    s = re.sub(r'(?<=\w)\s+\([^)]*\)', '', s)
    salts = r'(?:hydrochloride|dihydrochloride|trihydrochloride|hcl|2hcl|hydrobromide|maleate|dimaleate|mesylate|methanesulfonate|citrate|fumarate|tartrate|ditartrate|sodium salt|potassium salt|monohydrate|dihydrate|trihydrate|hemihydrate|hydrate|free base)'
    while True:
        cleaned = re.sub(r'\s+' + salts + r'\s*$', '', s).strip()
        if cleaned == s:
            break
        s = cleaned
    return key(s)

def metadata(path):
    d = pd.read_csv(path)
    cols = [x for x in ['broad_id', 'name', 'moa', 'target', 'disease.area', 'indication', 'smiles', 'phase'] if x in d]
    return d[cols].drop_duplicates().copy()

secondary = metadata(P / 'secondary_treatment_info.csv')
primary = metadata(P / 'primary_treatment_info.csv')

def match_sources(source):
    source = source.copy()
    source['exact_name_key'] = source['name'].map(key)
    source['alias_name_key'] = source['name'].map(alias_key)
    records = []
    for _, row in library.iterrows():
        matches = source[source.exact_name_key.eq(key(row['Item Name']))]
        method = 'exact_normalized_name'
        if matches.empty:
            matches = source[source.alias_name_key.eq(alias_key(row['Item Name']))]
            method = 'provisional_salt_or_parenthetical_alias'
        distinct_names = sorted(matches['name'].dropna().unique())
        if len(distinct_names) > 1:
            method = 'ambiguous_name_candidates'
        for _, match in matches.iterrows():
            records.append({
                'library_excel_row': row['library_excel_row'],
                'CatalogNumber': row['CatalogNumber'],
                'library_name': row['Item Name'],
                'library_cas': row['CAS Number'],
                'library_plate': row['Plate Location'],
                'library_rack': row['Rack Number'],
                'match_method': method,
                **{str(k): v for k, v in match.items() if k not in ['exact_name_key', 'alias_name_key']},
            })
    return pd.DataFrame(records)

sm = match_sources(secondary)
pm = match_sources(primary)
sm.to_csv(P / 'library_secondary_name_mapping.csv', index=False)
pm.to_csv(P / 'library_primary_name_mapping.csv', index=False)

coverage = library[['library_excel_row','CatalogNumber','Item Name','CAS Number','Plate Location','Rack Number']].copy()
for prefix, mapping in [('secondary', sm), ('primary', pm)]:
    g = mapping.groupby('library_excel_row').agg(
        matched_names=('name', lambda x: '; '.join(sorted(set(x.dropna())))),
        match_methods=('match_method', lambda x: '; '.join(sorted(set(x)))),
        broad_id_count=('broad_id', 'nunique'))
    coverage = coverage.join(g.rename(columns=lambda x: prefix + '_' + x), on='library_excel_row')
    coverage[prefix + '_broad_id_count'] = coverage[prefix + '_broad_id_count'].fillna(0).astype(int)
coverage.to_csv(P / 'library_prism_name_coverage.csv', index=False)

curves = pd.read_csv(P / 'secondary_curve_parameters_breast.csv')
treatment = pd.read_csv(P / 'secondary_treatment_info.csv')
doses = treatment.groupby(['broad_id', 'screen_id'], as_index=False).agg(
    dose_min_uM=('dose', 'min'), dose_max_uM=('dose', 'max'), number_of_doses=('dose','nunique'))
curves = curves.merge(doses, on=['broad_id','screen_id'], how='left', validate='many_to_one')
curves['auc_outside_0_1'] = (curves.auc < 0) | (curves.auc > 1)
curves['r2_nonpositive'] = curves.r2.le(0) | curves.r2.isna()
curves['ic50_nonfinite_or_missing'] = ~np.isfinite(curves.ic50)
curves['ic50_outside_tested_range'] = np.isfinite(curves.ic50) & ((curves.ic50 < curves.dose_min_uM) | (curves.ic50 > curves.dose_max_uM))
curves['ec50_outside_tested_range'] = np.isfinite(curves.ec50) & ((curves.ec50 < curves.dose_min_uM) | (curves.ec50 > curves.dose_max_uM))
curves['has_MTS010_same_name_cell'] = curves[['name','depmap_id']].apply(tuple, axis=1).isin(set(curves.loc[curves.screen_id.eq('MTS010'), ['name','depmap_id']].itertuples(index=False, name=None)))
curves['recommended_screen_available'] = ~curves.has_MTS010_same_name_cell | curves.screen_id.eq('MTS010')
curves.to_csv(P / 'secondary_curve_parameters_breast_annotated.csv', index=False)

merged = sm.merge(curves.drop(columns=['name','moa','target','disease.area','indication','smiles','phase']), on='broad_id', how='inner')
merged.to_csv(P / 'library_secondary_curves_breast_all_matches.csv', index=False)
merged[merged.depmap_id.eq('ACH-000849')].to_csv(P / 'library_secondary_curves_MDA_MB_468_all_matches.csv', index=False)
curves[curves.depmap_id.isin(['ACH-000849','ACH-000768','ACH-000573'])].to_csv(P / 'secondary_curve_parameters_468_231_436.csv', index=False)

summary = {
    'library_compounds_with_names': len(library),
    'source_secondary_unique_names': secondary.name.nunique(),
    'source_primary_unique_names': primary.name.nunique(),
    'secondary_library_names_matched_any': int(sm.library_excel_row.nunique()),
    'secondary_library_names_exact_normalized_name': int(sm.loc[sm.match_method.eq('exact_normalized_name'), 'library_excel_row'].nunique()),
    'secondary_library_names_provisional_salt_or_parenthetical_alias': int(sm.loc[sm.match_method.eq('provisional_salt_or_parenthetical_alias'), 'library_excel_row'].nunique()),
    'secondary_library_names_ambiguous': int(sm.loc[sm.match_method.eq('ambiguous_name_candidates'), 'library_excel_row'].nunique()),
    'primary_library_names_matched_any': int(pm.library_excel_row.nunique()),
    'primary_library_names_exact_normalized_name': int(pm.loc[pm.match_method.eq('exact_normalized_name'), 'library_excel_row'].nunique()),
    'primary_library_names_provisional_salt_or_parenthetical_alias': int(pm.loc[pm.match_method.eq('provisional_salt_or_parenthetical_alias'), 'library_excel_row'].nunique()),
    'primary_library_names_ambiguous': int(pm.loc[pm.match_method.eq('ambiguous_name_candidates'), 'library_excel_row'].nunique()),
    'breast_cells_with_curves': int(curves.depmap_id.nunique()),
    'MDA_MB_468_library_names_matched_any': int(merged.loc[merged.depmap_id.eq('ACH-000849'), 'library_excel_row'].nunique()),
    'MDA_MB_468_library_names_exact_normalized_name': int(merged.loc[merged.depmap_id.eq('ACH-000849') & merged.match_method.eq('exact_normalized_name'), 'library_excel_row'].nunique()),
    'note': 'Name-based mappings only. Provisional salt/alias mappings require manual chemical identity review. FDA approval is not inferred from vendor library title or Broad phase=Launched.',
}
(P / 'coverage_summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False))
print(json.dumps(summary, indent=2, ensure_ascii=False))
