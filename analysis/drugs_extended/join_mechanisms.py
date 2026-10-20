from pathlib import Path
import json,re
import pandas as pd

P=Path(__file__).resolve().parent
ROOT=P.parents[1]
lib=pd.read_csv(P/'whole_catalog_audit.csv')
curated=pd.read_csv(P/'curated_development_code_aliases.csv')
catalog_by_id=lib.set_index('CatalogNumber')

def key(s):return re.sub(r'[^a-z0-9]','',str(s).lower())
def generic(s):
    s=re.sub(r'\([^)]*\)',' ',str(s).lower())
    s=re.sub(r'\b(?:\d*hcl|hydrochloride|hydrobromide|mesylate|tosylate|ditosylate|maleate|dimaleate|tartrate|ditartrate|fumarate|sodium|disodium|potassium|dihydrate|trihydrate|monohydrate|hydrate|anhydrous|free base|free acid|salt)\b',' ',s)
    return key(s)
aliases={}
for _,r in lib.iterrows():
    names=[re.sub(r'\([^)]*\)',' ',str(r['Item Name']))]
    names += [x for text in re.findall(r'\(([^)]*)\)',str(r['Item Name'])) for x in re.split(r'[,;/]',text)]
    if pd.notna(r.curated_generic_name) and r.curated_generic_name:names.append(r.curated_generic_name)
    for n in names:
        if generic(n):aliases.setdefault(generic(n),set()).add(r.CatalogNumber)

evidence=[]
filenames=['evidence.csv'] if (ROOT/'analysis/mechanisms_extended/evidence.csv').exists() else ['mechanism_review_rows.csv','derl2_nme6_evidence.csv','activation_core_evidence.csv']
for filename in filenames:
    d=pd.read_csv(ROOT/'analysis/mechanisms_extended'/filename)
    d['mechanism_source_file']=filename
    d['evidence_id']=[filename.replace('.csv','')+'_'+str(i+1).zfill(3) for i in range(len(d))]
    evidence.append(d)
ev=pd.concat(evidence,ignore_index=True)
ev.to_csv(P/'gene_mechanism_evidence_used.csv',index=False)
edges=[];missing=[];unassigned=[]
name_only_compounds=['bortezomib','carfilzomib','ixazomib','chloroquine','hydroxychloroquine','gemcitabine','olaparib','metformin','docetaxel','paclitaxel']
for _,r in ev.iterrows():
    examples='' if pd.isna(r.compound_examples) else str(r.compound_examples)
    ids=set(re.findall(r'\b[A-Z]\d{4}\b',examples))
    ids.discard('L1021')
    valid=ids & set(catalog_by_id.index)
    for invalid in ids-set(catalog_by_id.index):missing.append({'evidence_id':r.evidence_id,'gene':r.gene,'invalid_catalog_id':invalid,'compound_examples':examples})
    if not ids:
        for name in name_only_compounds:
            if re.search(r'\b'+re.escape(name)+r'\b',examples,re.I):
                valid.update(aliases.get(generic(name),set()))
                if name=='ixazomib':valid.update(set(['A4007','A4008']) & set(catalog_by_id.index))
        # Explicit lyso-drug phosphate salts need named catalog-form retention;
        # they are not globally equated to all phosphate-containing molecules.
        if re.search(r'\bchloroquine\b',examples,re.I):
            valid.update(lib.loc[lib['Item Name'].str.match(r'^Chloroquine\b',case=False,na=False),'CatalogNumber'])
        if re.search(r'\bhydroxychloroquine\b',examples,re.I):
            valid.update(lib.loc[lib['Item Name'].str.match(r'^Hydroxychloroquine\b',case=False,na=False),'CatalogNumber'])
    if not valid:
        unassigned.append(r.to_dict());continue
    for cat in sorted(valid):
        lr=catalog_by_id.loc[cat]
        if isinstance(lr,pd.DataFrame):lr=lr.iloc[0]
        row=r.to_dict()
        row.update({'CatalogNumber':cat,'library_name':lr['Item Name'],'library_CAS':lr['CAS Number'],'catalog_row_id':lr.catalog_row_id,'gene_compound_link_method':'explicit catalog ID in reviewed evidence' if cat in ids else 'named compound matched to catalog active-moiety name; form retained','chemical_structure_validated':False})
        gene=str(r.gene)
        name=generic(lr.curated_generic_name) if pd.notna(lr.curated_generic_name) and lr.curated_generic_name else generic(lr['Item Name'])
        band='directional pathway/proteostasis/metabolic probe or contextual comparator'
        if gene=='FOXQ1' and 'imatinib' in name:band='gene-linked drug phenotype in external breast models; no olaparib proof'
        if gene=='CERK' and ('rapamycin' in name or 'sirolimus' in name):band='gene-linked downstream-drug phenotype in TNBC; migration endpoint'
        if gene=='DERL2' and 'gemcitabine' in name:band='gene-linked collateral drug phenotype in another tumour lineage'
        if gene=='PARG' and 'temozolomide' in name:band='gene-linked collateral drug phenotype in BRCA2-deficient external models'
        if gene=='ABCB1':band='transporter-function bypass hypothesis; assay substrate context required'
        if gene=='PHGDH' and 'disulfiram' in name:band='nonselective biochemical target-associated probe; context counterevidence'
        if 'counterevidence' in str(r.drug_role).lower():band='counterevidence/context exclusion comparator'
        row['evidence_priority_category']=band
        row['olaparib_resensitization_in_user_model_demonstrated']=False
        row['FDA_record_match_with_curated_alias']=lr.FDA_record_match_with_curated_alias
        row['regulatory_review_category']=lr.regulatory_review_category
        row['regulatory_evidence_note']=lr.regulatory_evidence_note
        edges.append(row)

edges=pd.DataFrame(edges)
edges.to_csv(P/'gene_drug_evidence_edges.csv',index=False)
pd.DataFrame(unassigned).to_csv(P/'mechanism_evidence_without_catalog_drug.csv',index=False)
pd.DataFrame(missing).to_csv(P/'unresolved_catalog_ids_in_mechanism_review.csv',index=False)

# One catalog compound profile. Keep ambiguous source names and dose windows
# separate rather than averaging native AUC across drug IDs or platforms.
primary=pd.read_csv(P/'prism_primary_catalog_breast_records.csv')
secondary=pd.read_csv(P/'prism_secondary_catalog_breast_records.csv')
gdsc={name:pd.read_csv(P/'public_crosscheck'/f'{name.lower()}_breast_curve_records.csv') for name in ['GDSC1','GDSC2']}
gdsc_catalog={name:pd.read_csv(P/(name.lower()+'_catalog_breast_records.csv')) for name in ['GDSC1','GDSC2']}
profiles=[]
profile_catalogs=list(edges.CatalogNumber.drop_duplicates())
if 'B1464' in catalog_by_id.index and 'B1464' not in profile_catalogs: profile_catalogs.append('B1464')
for cat in profile_catalogs:
    lr=catalog_by_id.loc[cat]
    if isinstance(lr,pd.DataFrame):lr=lr.iloc[0]
    er=edges[edges.CatalogNumber.eq(cat)]
    out={'CatalogNumber':cat,'library_name':lr['Item Name'],'library_CAS':lr['CAS Number'],'curated_generic_name':lr.curated_generic_name,'gene_axes':'; '.join(sorted(set(er.gene))), 'evidence_priority_categories':'; '.join(sorted(set(er.evidence_priority_category))), 'drug_roles':'; '.join(sorted(set(er.drug_role))), 'FDA_record_match_with_curated_alias':lr.FDA_record_match_with_curated_alias,'FDA_marketing_record_current_with_curated_alias':lr.FDA_marketing_record_current_with_curated_alias,'regulatory_review_category':lr.regulatory_review_category,'regulatory_evidence_note':lr.regulatory_evidence_note,'chemical_structure_validated':False}
    if cat=='B1464':
        out['gene_axes']='Empirical olaparib combination comparator'
        out['evidence_priority_categories']='exact parental468 combination evidence; no target-specific attribution'
        out['drug_roles']='Selinexor/XPO1 empirical olaparib-combination comparator; not assigned to an individual screen target'
    out['comparison_role']='reference olaparib treatment; not a new sensitizer' if 'olaparib' in generic(lr['Item Name']) else 'gene-linked nominated compound or mechanistic/context comparator; see drug_roles'
    for name,data in [('primary',primary),('secondary',secondary)]:
        z=data[data.CatalogNumber.eq(cat)&data.depmap_id.eq('ACH-000849')].drop_duplicates('source_compound_name')
        out[name+'_468_source_names']='; '.join(sorted(set(z.source_compound_name)))
        out[name+'_468_identity_mapping_methods']='; '.join(sorted(set(z.identity_mapping_method)))
        out[name+'_468_n_source_names']=len(z)
        if name=='primary':
            out['primary_468_log2FC']=z.primary_log2_fold_change_median.iloc[0] if len(z)==1 else None
            out['primary_468_relative_viability']=z.primary_viability_relative_DMSO_from_median_LFC.iloc[0] if len(z)==1 else None
            out['primary_468_dose_min_uM']=z.primary_dose_min_uM.iloc[0] if len(z)==1 else None
            out['primary_468_dose_max_uM']=z.primary_dose_max_uM.iloc[0] if len(z)==1 else None
        else:
            fields=['auc_median','r2_median','selected_screen_ids','dose_min_uM_across_selected_records','dose_max_uM_across_selected_records','all_selected_records_positive_r2_and_auc_0_1','ic50_uM_median_within_tested_range']
            for f in fields:out['secondary_468_'+f]=z[f].iloc[0] if len(z)==1 else None
    for name,data in gdsc_catalog.items():
        z=data[data.CatalogNumber.eq(cat)&data.depmap_id.eq('ACH-000849')].drop_duplicates(['DRUG_ID','NLME_CURVE_ID'])
        out[name+'_468_n_distinct_drug_IDs']=z.DRUG_ID.nunique()
        out[name+'_468_n_curves']=len(z)
        out[name+'_468_n_in_range_IC50']=int(z.ic50_in_tested_range_uM.notna().sum())
        detailcols=['DRUG_ID','DRUG_NAME','MIN_CONC','MAX_CONC','AUC','ic50_fitted_uM','ic50_range_flag','ic50_in_tested_range_uM']
        out[name+'_468_curve_details_JSON']=z[detailcols].to_json(orient='records') if len(z) else '[]'
    profiles.append(out)
pd.DataFrame(profiles).to_csv(P/'gene_linked_catalog_compound_profiles.csv',index=False)

# Panel summaries are stratified by assay concentration context, not pooled
# across GDSC IDs or unrelated PRISM dose windows.
pm=pd.read_csv(ROOT/'analysis/prism/primary_metadata.csv').dropna(subset=['depmap_id'])
tnbc_ids=set(pm.loc[pm.secondary_tissue.eq('TNBC')&pm.passed_str_profiling.eq(True),'depmap_id'])
tnbc_primary=primary[primary.depmap_id.isin(tnbc_ids)].drop_duplicates(['catalog_row_id','source_compound_name','depmap_id'])
tp=tnbc_primary.groupby(['catalog_row_id','CatalogNumber','library_name','source_compound_name','primary_dose_min_uM','primary_dose_max_uM'],as_index=False).agg(TNBC_models=('depmap_id','nunique'),primary_log2FC_TNBC_median=('primary_log2_fold_change_median','median'),primary_log2FC_TNBC_min=('primary_log2_fold_change_median','min'),primary_log2FC_TNBC_max=('primary_log2_fold_change_median','max'))
tp.to_csv(P/'prism_primary_TNBC_concentration_context_summary.csv',index=False)
ts=secondary[secondary.secondary_tissue.eq('TNBC')].drop_duplicates(['catalog_row_id','source_compound_name','depmap_id'])
ts.groupby(['catalog_row_id','CatalogNumber','library_name','source_compound_name','selected_screen_ids','dose_min_uM_across_selected_records','dose_max_uM_across_selected_records'],as_index=False).agg(TNBC_models=('depmap_id','nunique'),secondary_auc_TNBC_median=('auc_median','median'),secondary_auc_TNBC_min=('auc_median','min'),secondary_auc_TNBC_max=('auc_median','max'),n_models_all_curve_diagnostics_pass=('all_selected_records_positive_r2_and_auc_0_1','sum')).to_csv(P/'prism_secondary_TNBC_concentration_context_summary.csv',index=False)

print('Gene evidence rows',len(ev),'edges',len(edges),'catalog compounds',edges.CatalogNumber.nunique(),'genes',edges.gene.nunique(),'unresolved IDs',len(missing))
print('Profile genes:',sorted(set(edges.gene)))
print('No catalog edge genes:',sorted(set(ev.gene)-set(edges.gene)))
