from pathlib import Path
import json
import numpy as np
import pandas as pd

P = Path(__file__).resolve().parent
curves = pd.read_csv(P / 'secondary_curve_parameters_breast_annotated.csv')
selected = curves[curves['recommended_screen_available']].copy()
selected.to_csv(P / 'secondary_selected_curve_records_breast.csv', index=False)
metadata = pd.read_csv(P / 'secondary_breast_metadata.csv').dropna(subset=['depmap_id']).drop_duplicates('depmap_id')

def summarize(g):
    valid_ic = g.loc[~g.ic50_nonfinite_or_missing & ~g.ic50_outside_tested_range, 'ic50']
    return pd.Series({
        'ccle_name': g.ccle_name.iloc[0],
        'n_selected_curve_records': len(g),
        'selected_screen_ids': '; '.join(sorted(g.screen_id.unique())),
        'selected_broad_ids': '; '.join(sorted(g.broad_id.unique())),
        'selection_rule': 'MTS010 technical redo when available; otherwise median across available non-MTS010 batches/screens',
        'auc_median': g.auc.median(),
        'auc_min': g.auc.min(),
        'auc_max': g.auc.max(),
        'r2_median': g.r2.median(),
        'n_r2_nonpositive_or_missing': int(g.r2_nonpositive.sum()),
        'n_auc_outside_0_1': int(g.auc_outside_0_1.sum()),
        'n_ic50_missing_or_nonfinite': int(g.ic50_nonfinite_or_missing.sum()),
        'n_ic50_outside_tested_range': int(g.ic50_outside_tested_range.sum()),
        'n_ic50_finite_within_tested_range': len(valid_ic),
        'ic50_uM_median_within_tested_range': valid_ic.median() if len(valid_ic) else np.nan,
        'dose_min_uM_across_selected_records': g.dose_min_uM.min(),
        'dose_max_uM_across_selected_records': g.dose_max_uM.max(),
        'has_at_least_one_positive_r2_curve': bool((~g.r2_nonpositive).any()),
        'all_selected_records_positive_r2_and_auc_0_1': bool((~g.r2_nonpositive & ~g.auc_outside_0_1).all()),
    })

summary = selected.groupby(['name', 'depmap_id'], sort=True).apply(summarize, include_groups=False).reset_index()
summary = summary.merge(metadata[['depmap_id','secondary_tissue','passed_str_profiling']], on='depmap_id', how='left', validate='many_to_one')
summary.to_csv(P / 'secondary_compound_cell_summary.csv', index=False)

mapping = pd.read_csv(P / 'library_secondary_name_mapping.csv')
# A compound may have multiple screened batches; retain one library-to-name
# association for per-compound aggregation and keep batch details in raw files.
mapcols = ['library_excel_row','CatalogNumber','library_name','library_cas','library_plate','library_rack','match_method','name','moa','target','phase']
mapping = mapping[mapcols].drop_duplicates(['library_excel_row','name'])
joined = mapping.merge(summary, on='name', how='inner')
fda_path = P.parent / 'fda/library_fda_annotation.csv'
if fda_path.exists():
    fda = pd.read_csv(fda_path)
    fcols = ['Item Name','FDA_record_match','FDA_current_product_record','FDA_ingredient_matches','FDA_mapping_note']
    joined = joined.merge(fda[fcols].rename(columns={'Item Name':'library_name'}), on='library_name', how='left', validate='many_to_one')
joined.to_csv(P / 'library_secondary_compound_cell_summary.csv', index=False)
joined[joined.depmap_id.eq('ACH-000849')].sort_values('auc_median').to_csv(P / 'library_secondary_MDA_MB_468_summary.csv', index=False)

candidate_names = ['lapatinib', 'verapamil', 'temozolomide', 'fludarabine', 'exemestane', 'belinostat', 'vorinostat', 'gemcitabine', 'panobinostat', 'omacetaxine', 'chlorambucil']
candidate = summary[summary.name.str.lower().isin(candidate_names) & summary.depmap_id.isin(['ACH-000849','ACH-000768','ACH-000573'])].copy()
candidate.to_csv(P / 'requested_candidates_468_231_436_summary.csv', index=False)

if 'FDA_record_match' in joined:
    top = joined[joined.depmap_id.eq('ACH-000849') & joined.FDA_record_match.fillna(False) & joined.all_selected_records_positive_r2_and_auc_0_1].copy()
    top = top.sort_values(['auc_median','library_name']).drop_duplicates('name')
    top.to_csv(P / 'FDA_record_matched_468_activity_summary.csv', index=False)
    print('TOP FDA record matched drugs, exploratory monotherapy activity:')
    print(top[['library_name','name','auc_median','r2_median','selected_screen_ids','ic50_uM_median_within_tested_range']].head(20).to_string(index=False))
print('Requested candidates MDA-MB-468:')
print(candidate[candidate.depmap_id.eq('ACH-000849')][['name','auc_median','r2_median','selected_screen_ids','ic50_uM_median_within_tested_range','n_ic50_outside_tested_range']].to_string(index=False))
print('Requested candidate names missing all breast curves:', sorted(set(candidate_names) - set(summary.name.str.lower())))
