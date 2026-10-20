from pathlib import Path
import pandas as pd

P = Path(__file__).resolve().parent
d = pd.read_csv(P / 'primary_logfold_change_breast_long.csv')
metadata = pd.read_csv(P / 'primary_metadata.csv')
d = d.merge(metadata[['row_name','passed_str_profiling']], on='row_name', how='left', validate='many_to_one')
# Exclude unassigned or misidentified barcodes from summaries. Raw files retain
# the original extracted rows, including the documented identity metadata.
d = d[d.depmap_id.notna() & d.passed_str_profiling.eq(True)].copy()

summary = d.groupby(['name','depmap_id','ccle_name'], as_index=False).agg(
    primary_log2_fold_change_median=('log2_fold_change','median'),
    primary_log2_fold_change_min=('log2_fold_change','min'),
    primary_log2_fold_change_max=('log2_fold_change','max'),
    primary_n_measurements=('log2_fold_change','size'),
    primary_screen_ids=('screen_id', lambda x: '; '.join(sorted(x.unique()))),
    primary_dose_min_uM=('dose','min'),
    primary_dose_max_uM=('dose','max'),
    primary_n_distinct_doses=('dose','nunique'),
)
summary['primary_viability_relative_DMSO_from_median_LFC'] = 2**summary.primary_log2_fold_change_median
summary.to_csv(P / 'primary_compound_cell_summary.csv', index=False)

mapping = pd.read_csv(P / 'library_primary_name_mapping.csv')
cols = ['library_excel_row','CatalogNumber','library_name','library_cas','library_plate','library_rack','match_method','name','moa','target','phase']
mapping = mapping[cols].drop_duplicates(['library_excel_row','name'])
joined = mapping.merge(summary,on='name',how='inner')
fda_file = P.parent / 'fda/library_fda_annotation.csv'
if fda_file.exists():
    fda = pd.read_csv(fda_file)
    cols = ['Item Name','FDA_record_match','FDA_current_product_record','FDA_ingredient_matches','FDA_mapping_note']
    joined = joined.merge(fda[cols].rename(columns={'Item Name':'library_name'}), on='library_name', how='left', validate='many_to_one')
joined.to_csv(P / 'library_primary_compound_cell_summary.csv',index=False)
joined[joined.depmap_id.eq('ACH-000849')].to_csv(P / 'library_primary_MDA_MB_468_summary.csv',index=False)
for cell in ['ACH-000849','ACH-000768','ACH-000573']:
    z=joined[joined.depmap_id.eq(cell)]
    print(cell,'primary names:',z.library_excel_row.nunique(),'exact:',z.loc[z.match_method.eq('exact_normalized_name'),'library_excel_row'].nunique())
