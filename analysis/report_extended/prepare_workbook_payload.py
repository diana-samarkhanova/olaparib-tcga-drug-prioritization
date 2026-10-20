from pathlib import Path
import json
import pandas as pd
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
T=ROOT/'analysis/tcga_extended'; D=ROOT/'analysis/drugs_extended'; A=ROOT/'analysis/astra_extended'
O=ROOT/'outputs/TCGA_Drug_Analysis_EN_20261013'
genes=['LAMTOR2','WDR24','DERL2','NME6','TP53BP1','PARG','TP53I11','FOXQ1','CERK','PHF19','SLC35D1','PHLDA2','WT1','HEPHL1','LINC03007','SLC47A1','ELL2','INHBA','PHGDH','ABCB1']

def clean(df):
    return df.astype(object).where(pd.notna(df),None)
def spec(name,df,title=None,widths=None,formats=None,preview_columns=8,wrap_columns=None):
    numeric=[i for i,c in enumerate(df.columns) if pd.api.types.is_numeric_dtype(df[c]) and not pd.api.types.is_bool_dtype(df[c])]
    fm={str(i):'0.000E+00' if any(v in df.columns[i].lower() for v in ['fdr','p_value','pvalue','ic50','conc','dose_']) else ('0' if pd.api.types.is_integer_dtype(df.iloc[:,i]) or df.columns[i].lower().startswith('n_') or any(v in df.columns[i].lower() for v in ['_n','events','sgrna','sample_n','year','source_excel_row','drug_id']) else '0.000') for i in numeric}
    if formats:fm.update({str(k):v for k,v in formats.items()})
    return {'name':name,'title':title or name,'output':title is not None,'headers':list(df.columns),'rows':clean(df).values.tolist(),'widths':widths,'numeric_columns':numeric,'formats':fm,'preview_columns':preview_columns,'wrap_columns':wrap_columns or []}

scr=pd.read_csv(O/'tables/Selected_CRISPR_panel_20_genes.csv')
expression=pd.read_csv(T/'expression_contrasts.csv')
coverage=pd.read_csv(T/'gene_coverage.csv')
assert coverage.query("gene=='LAMTOR2'").source_symbol.iloc[0]=='ROBLD3'
assert coverage.rna_coverage.eq('available').sum()==19
cox=pd.read_csv(T/'survival_continuous_Cox.csv')
cn=pd.read_csv(T/'CNA_prevalence_and_expression.csv')
gse=pd.read_csv(A/'gse165914_candidate_context.csv')
priorities=pd.read_csv(D/'final_model_priorities.csv')
scientific=pd.read_csv(A/'final_gene_priority_matrix.csv')
summary=scientific[['gene','perturbation','user_validation_group','scientific_priority','mechanistic_program']].copy()
for time in ['early','late']:
    rows=scr.groupby('gene').nth(0 if time=='early' else 1).reset_index()
    cols=rows[['gene','log2FC','positive_FDR','supporting_sgRNAs']].rename(columns={c:f'{time}_{c}' for c in ['log2FC','positive_FDR','supporting_sgRNAs']})
    summary=summary.merge(cols,on='gene',validate='one_to_one')
e=expression.query("comparison=='clinical_TNBC_vs_nonTNBC'")[['gene','median_delta_log2_RSEM_plus1','BH_FDR']].rename(columns={'median_delta_log2_RSEM_plus1':'TCGA_TNBC_median_delta_log2_RSEM_plus1','BH_FDR':'TCGA_TNBC_expression_BH_FDR'})
summary=summary.merge(e,on='gene',how='left',validate='one_to_one')
summary=summary.merge(gse[['gene','resistant_15uM_vs_parental_untreated_author_log2FC','resistant_15uM_vs_parental_untreated_author_FDR','contextual_direction_vs_screen']],on='gene',how='left',validate='one_to_one')
summary=summary.merge(priorities[['gene','nomination_category','primary_catalog_number','action','interpretation_gate']],on='gene',how='left',validate='one_to_one')
summary=summary.set_index('gene').loc[genes].reset_index()
summary.to_csv(O/'tables/Integrated_gene_summary.csv',index=False)

sheets=[]
sheets.append(spec('Gene summary',summary,'CRISPR candidate and drug evidence',widths=[16,14,24,28,44]+[18]*(len(summary.columns)-5),preview_columns=5))
sheets.append(spec('Screen results',scr,widths=[16,14,27,16,14,16,18,18,16,50,20],preview_columns=9))
react=pd.read_csv(ROOT/'analysis/crispr/enrichment/KO_reactome_ora_all_tested.csv')
first=['pathway_id','pathway_name','overlap_genes','foreground_pathway_overlap','fold_enrichment','bh_fdr','p_value']
react=react[first+[c for c in react.columns if c not in first]]
sheets.append(spec('KO Reactome enrichment',react,widths=[22,65,35,24,22,22,22]+[24]*(len(react.columns)-7),preview_columns=7,wrap_columns=[1]))
sheets.append(spec('TCGA expression',expression,widths=[38,16,16,18,18,20,20,24,20,20,20],preview_columns=8))
sheets.append(spec('TCGA CNA',cn,widths=[16,16]+[20]*(len(cn.columns)-2),preview_columns=9))
sheets.append(spec('TCGA mutations',pd.read_csv(T/'mutation_prevalence.csv'),widths=[16,16,35]+[24]*6,preview_columns=9))
sheets.append(spec('TCGA somatic variants',pd.read_csv(T/'selected_somatic_variants.csv'),widths=None,preview_columns=7))
sheets.append(spec('TCGA pathways',pd.read_csv(T/'pathway_associations.csv'),widths=[16,16,38]+[22]*11,preview_columns=8))
sheets.append(spec('TCGA RPPA',pd.read_csv(T/'RPPA_associations.csv'),widths=[16,16,43,16,20,20,20],preview_columns=7))
sheets.append(spec('TCGA survival',cox,widths=[16,16,16,27]+[22]*(len(cox.columns)-4),preview_columns=9))
sheets.append(spec('PFI subtype sensitivity',pd.read_csv(T/'survival_PAM50_sensitivity.csv'),widths=[32,18,18,36,20,20,100]+[23]*9,preview_columns=6))
sheets.append(spec('Resistance RNA',gse,widths=[16,16,24]+[28]*(len(gse.columns)-3),preview_columns=7))
selected=priorities[['gene','perturbation','nomination_category','primary_catalog_number','library_name','action','interpretation_gate','primary_evidence_title','primary_evidence_url','regulatory_review_category','regulatory_evidence_note','edited_MDA_MB_468_resensitization_demonstrated']]
sheets.append(spec('Drug priorities',selected,widths=[16,16,32,20,32,70,75,75,60,70,80,36],preview_columns=5,wrap_columns=[5,6,7,9,10]))
profiles=pd.read_csv(D/'final_focus_compound_profiles.csv')
profile_first=['CatalogNumber','library_name','gene_axes','primary_468_relative_viability','primary_468_dose_min_uM','secondary_468_auc_median','secondary_468_r2_median','secondary_468_ic50_uM_median_within_tested_range']
profiles=profiles[profile_first+[c for c in profiles.columns if c not in profile_first]]
sheets.append(spec('Drug profiles',profiles,widths=[20,40,40]+[24]*(len(profiles.columns)-3),preview_columns=8))
sheets.append(spec('GDSC dose curves',pd.read_csv(D/'final_focus_GDSC_curves.csv'),widths=None,preview_columns=9))
sheets.append(spec('Olaparib combinations',pd.read_csv(ROOT/'analysis/combinations_extended/combo_evidence.csv'),widths=[20,35,20,40,85,65,22,100,55,60,110,75,30,95,30],preview_columns=5,wrap_columns=[3,4,5,7,9,10,13]))
sheets.append(spec('Published screen hits',pd.read_csv(ROOT/'analysis/combinations_extended/Yusoh_nine_library_hits.csv'),widths=None,preview_columns=8))
lib=pd.read_csv(D/'whole_catalog_audit.csv')
keep=['CatalogNumber','Item Name','CAS Number','Plate Location','Rack Number','curated_generic_name','regulatory_review_category','PRISM_primary_MDA_MB_468_numeric_record','PRISM_secondary_MDA_MB_468_numeric_record','GDSC1_MDA_MB_468_numeric_record','GDSC2_MDA_MB_468_numeric_record','FDA_source','URL']
sheets.append(spec('Library coverage',lib[keep],widths=[20,45,20,18,18,35,65,35,35,35,35,60,60],preview_columns=6))
evidence=pd.read_csv(ROOT/'analysis/mechanisms_extended/evidence.csv')
evidence['pmid']=evidence['pmid'].map(lambda x:'' if pd.isna(x) else str(int(x)))
keep=['gene','branch','evidence_level','model','perturbation','observed_direction','drug_hypothesis','compound_examples','drug_role','parp_specific_evidence','counterevidence','readout','primary_title','year','pmid','doi','url']
sheets.append(spec('Mechanism evidence',evidence[keep],widths=[16,45,40,55,65]+[70]*(len(keep)-5),preview_columns=4,wrap_columns=list(range(1,13))))
sheets.append(spec('Gene interpretation',scientific,widths=[16,16,28,48,32]+[70]*(len(scientific.columns)-5),preview_columns=5))
counts=json.loads((T/'cohort_counts.json').read_text())
methods=pd.DataFrame([
 ['TCGA expression','RSEM batch-normalized log2(RSEM+1), one primary tumour per patient. LAMTOR2 mapped to ROBLD3 Entrez28956; LINC03007 unavailable.','https://github.com/cBioPortal/datahub/tree/master/public/brca_tcga_pan_can_atlas_2018'],
 ['Clinical TNBC','ER-negative, PR-negative and resolved HER2-negative; discordant HER2 excluded from negative classification.','https://github.com/cBioPortal/datahub/tree/master/public/brca_tcga'],
 ['Pathways','Hallmark mean gene z scores; candidate excluded. Partial rank correlation adjusts age, E2F, target CNA and purity where covered; BRCA additionally PAM50.','https://www.gsea-msigdb.org/gsea/msigdb/human/collection_details.jsp#H'],
 ['Survival','Continuous-expression Cox HR per SD; BRCA age/stage adjustment, subgroup univariable due low event counts. PFI primary, DFI supporting, OS exploratory.','https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018'],
 ['GSE165914','Targeted extraction of author genome-wide DE, with exposure/resistance confounding. No new DE model.','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165914'],
 ['PRISM','Primary fixed-dose log2 viability, secondary source AUC/fit parameters. MTS010 redo preferred. Salt/development-code mappings provisional.','https://depmap.org/repurposing/'],
 ['GDSC','Release 8.5 fitted 72-hour curves. exp(LN_IC50) in µM; within-range fits separate from extrapolations. AUCs not pooled.','https://www.cancerrxgene.org/downloads/bulk_download'],
 ['Drug evidence','No matched cell-line expression-response model, docking or new combination synergy fit. Public parental sensitivity is not engineered-model selectivity.','https://depmap.org/repurposing/'],
 ['FDA status','Ingredient/product records are partial automated mapping. Tazemetostat/panobinostat historical withdrawn; elacridar/zosuquidar investigational.','https://www.fda.gov/drugs/resources-information-approved-drugs/withdrawn-cancer-accelerated-approvals'],
 ['Numerical meaning','Missing numerical results are blank, not zero; Δmedian log2(RSEM+1) is not log2 fold change. Study-native FDR families retained.',''],
],columns=['Resource','Method and limitation','Source URL'])
sheets.append(spec('Methods',methods,widths=[25,115,80],preview_columns=3,wrap_columns=[1]))
payload={'actual_results_ready':True,'sheets':sheets}
(ROOT/'analysis/report_extended/workbook_payload.json').write_text(json.dumps(payload,allow_nan=False,default=lambda x: x.item() if isinstance(x,np.generic) else str(x)))
print(json.dumps({'sheets':len(sheets),'gene_rows':len(summary),'TCGA_genes':int(coverage.rna_coverage.eq('available').sum()),'library_rows':len(lib)}))
