"""Audit the entire uploaded catalog against actual public monotherapy data.

Identity matches remain name/synonym based. No combined score, expression-based
response predictor or predicted synergy is produced.
"""
from pathlib import Path
import json,re
import pandas as pd

P=Path(__file__).resolve().parent
ROOT=P.parents[1]
PRISM=ROOT/'analysis/prism'
lib=pd.read_csv(ROOT/'analysis/fda/library_fda_annotation.csv')
lib['catalog_row_id']=range(1,len(lib)+1)
curated=pd.read_csv(P/'curated_development_code_aliases.csv')

def exact_key(s):return re.sub(r'[^a-z0-9]','',str(s).lower())
def clean(s):
    s=str(s).lower().replace('α','alpha').replace('β','beta')
    s=re.sub(r'\([^)]*\)',' ',s)
    # This is a provisional active-moiety name heuristic, never chemical proof.
    s=re.sub(r'\b(?:\d*hcl|hydrochloride|hydrobromide|mesylate|tosylate|ditosylate|besylate|maleate|dimaleate|tartrate|ditartrate|fumarate|sodium|disodium|potassium|dihydrate|trihydrate|monohydrate|hydrate|anhydrous|free base|free acid|salt)\b',' ',s)
    return exact_key(s)

def aliases(name):
    result={clean(name)}
    for inner in re.findall(r'\(([^)]*)\)',str(name)):
        result.update(clean(x) for x in re.split(r'[,;/]',inner))
    return result-{''}

alias_by_catalog={}
for _,r in curated.iterrows():
    keys={clean(r.generic_alias),clean(r.development_code)}
    f=P/'identity_sources'/f'{r.CatalogNumber}_pubchem_synonyms.json'
    if f.exists() and r.CatalogNumber not in ['A4007','A4008']:
        j=json.loads(f.read_text())
        for entry in j.get('InformationList',{}).get('Information',[]):
            keys.update(clean(x) for x in entry.get('Synonym',[]))
    alias_by_catalog[r.CatalogNumber]=keys-{''}

sources={
 'PRISM_primary':pd.read_csv(PRISM/'primary_compound_cell_summary.csv'),
 'PRISM_secondary':pd.read_csv(PRISM/'secondary_compound_cell_summary.csv'),
 'GDSC1':pd.read_csv(P/'public_crosscheck/gdsc1_breast_curve_records.csv'),
 'GDSC2':pd.read_csv(P/'public_crosscheck/gdsc2_breast_curve_records.csv'),
}

mappings=[];joined_sources={}
for source,df in sources.items():
    namecol='DRUG_NAME' if source.startswith('GDSC') else 'name'
    names=df[namecol].dropna().unique()
    source_exact={};source_clean={}
    for name in names:
        source_exact.setdefault(exact_key(name),set()).add(name)
        for k in aliases(name):source_clean.setdefault(k,set()).add(name)
    entries=[]
    for _,r in lib.iterrows():
        direct=source_exact.get(exact_key(r['Item Name']),set())
        method='direct_normalized_name'
        matched=set(direct)
        if not matched:
            for k in aliases(r['Item Name']):matched.update(source_clean.get(k,set()))
            method='provisional_declared_alias_or_salt_name'
        if not matched and r.CatalogNumber in alias_by_catalog:
            for k in alias_by_catalog[r.CatalogNumber]:matched.update(source_clean.get(k,set()))
            method='curated_development_code_or_PubChem_synonym; form unverified'
        for name in sorted(matched):
            entry={'catalog_row_id':r.catalog_row_id,'CatalogNumber':r.CatalogNumber,'library_name':r['Item Name'],'library_CAS':r['CAS Number'],'source_dataset':source,'source_compound_name':name,'identity_mapping_method':method,'multiple_source_names_for_catalog_compound':len(matched)>1,'chemical_structure_validated':False}
            entries.append(entry);mappings.append(entry)
    mapping=pd.DataFrame(entries)
    merged=mapping.merge(df,left_on='source_compound_name',right_on=namecol,how='inner')
    merged.to_csv(P/(source.lower()+'_catalog_breast_records.csv'),index=False)
    merged[merged.depmap_id.eq('ACH-000849')].to_csv(P/(source.lower()+'_catalog_MDA_MB_468.csv'),index=False)
    joined_sources[source]=merged

mapping=pd.DataFrame(mappings)
mapping.to_csv(P/'catalog_source_identity_mapping.csv',index=False)

# Supplement FDA ingredient-record matching for development-code aliases only.
# The published withdrawal table supersedes a stale current marketing record.
products=pd.read_csv(ROOT/'analysis/fda/Products.txt',sep='\t',encoding='cp1252',low_memory=False)
marketing=pd.read_csv(ROOT/'analysis/fda/MarketingStatus.txt',sep='\t',encoding='cp1252')
lookup=pd.read_csv(ROOT/'analysis/fda/MarketingStatus_Lookup.txt',sep='\t',encoding='cp1252')
sub=pd.read_csv(ROOT/'analysis/fda/Submissions.txt',sep='\t',encoding='cp1252',low_memory=False)
approved_apps=set(sub.loc[sub.SubmissionStatus.eq('AP'),'ApplNo'])
products=products[products.ApplNo.isin(approved_apps)].merge(marketing,on=['ApplNo','ProductNo'],how='left').merge(lookup,on='MarketingStatusID',how='left')
ingredients={}
for _,r in products.iterrows():
    for ingredient in str(r.ActiveIngredient).split(';'):
        ingredients.setdefault(clean(ingredient),[]).append((ingredient.strip(),r.ApplNo,r.MarketingStatusDescription))
alias_rows=curated.set_index('CatalogNumber').to_dict('index')
lib['curated_generic_name']=''
lib['FDA_record_match_with_curated_alias']=lib.FDA_record_match
lib['FDA_marketing_record_current_with_curated_alias']=lib.FDA_current_product_record
lib['FDA_curated_ingredient_matches']=lib.FDA_ingredient_matches.fillna('')
lib['regulatory_evidence_note']='Drugs@FDA ingredient/product records only; current indication not independently adjudicated'
lib['regulatory_review_category']=lib.FDA_record_match.map({True:'FDA ingredient record; indication-specific review required',False:'FDA nonmatch unresolved; no approval claim'})
for i,r in lib.iterrows():
    alias=alias_rows.get(r.CatalogNumber)
    if alias:
        lib.at[i,'curated_generic_name']=alias['generic_alias']
        hits=ingredients.get(clean(alias['generic_alias']),[])
        # MLN2238 is an active form distinct from the marketed prodrug.
        if not hits and alias['generic_alias']=='ixazomib':hits=ingredients.get(clean('ixazomib citrate'),[])
        if not hits and alias['generic_alias']=='sunitinib':hits=ingredients.get(clean('sunitinib malate'),[])
        if hits:
            lib.at[i,'FDA_record_match_with_curated_alias']=True
            lib.at[i,'FDA_marketing_record_current_with_curated_alias']=any(x[2] in ['Prescription','Over-the-counter'] for x in hits)
            lib.at[i,'FDA_curated_ingredient_matches']='; '.join(sorted({x[0] for x in hits}))
            lib.at[i,'regulatory_review_category']='FDA active-moiety record via curated alias; catalog form not proven marketed form'
            lib.at[i,'regulatory_evidence_note']=alias['chemical_form_relation']
    k=clean(r['Item Name']);generic=clean(alias['generic_alias']) if alias else k
    if generic=='tazemetostat' or k=='epz6438':
        lib.at[i,'regulatory_review_category']='historical FDA approval; all Tazverik approvals withdrawn 2026-06-22; research probe only'
        lib.at[i,'regulatory_evidence_note']='Official FDA withdrawn cancer accelerated approvals table; overrides product-record matching. https://www.fda.gov/drugs/resources-information-approved-drugs/withdrawn-cancer-accelerated-approvals'
    if generic=='panobinostat':
        lib.at[i,'regulatory_review_category']='historical FDA approval; Farydak withdrawn 2022-03-24; separate historical/probe branch'
        lib.at[i,'regulatory_evidence_note']='Official FDA withdrawn cancer accelerated approvals table. https://www.fda.gov/drugs/resources-information-approved-drugs/withdrawn-cancer-accelerated-approvals'
    if generic in ['romidepsin','idelalisib','duvelisib','selinexor']:
        lib.at[i,'regulatory_evidence_note']='Some oncology indications withdrawn; do not infer all approvals withdrawn from one indication. Review current label by indication.'

# Annotation-driven whole-catalog inventory: these are retrieval classes only.
# Individual gene/drug claims are made separately in the evidence-edge table.
patterns={
 'PDGFR_receptor_kinases':r'PDGFR|platelet.derived growth',
 'MEK_ERK_MAPK':r'\bMEK\b|\bERK\b|MAPK inhibitor',
 'PI3K_AKT_mTOR':r'PI3K|\bAKT\b|mTOR',
 'proteasome_VCP_ERAD':r'Proteasome|\bVCP\b|\bp97\b',
 'HSP90_ER_stress_UPR':r'Hsp90|heat shock|HSP antagonist|ER stress|endoplasmic|unfolded protein',
 'lysosome_autophagy':r'Autophagy|lysosom',
 'EZH2_PRC2_chromatin':r'EZH2|Histone Methyltransferase|PRC2',
 'HDAC_DNMT_epigenetic':r'HDAC|histone deacetyl|DNA Methyltransferase',
 'DNA_replication_antimetabolite':r'DNA Synthesis|Thymidylate|Ribonucleotide reductase|purine antagonist|pyrimidine antagonist|folate',
 'DDR_checkpoint_PARPi':r'\bATR\b|\bATM\b|WEE1|CHK1|CHK2|DNA.PK|PARP',
 'ceramide_sphingolipid':r'ceramide|CERK|sphingosine|sphingolipid',
 'PHGDH_serine_onecarbon':r'PHGDH|serine biosynth|SHMT|one.carbon',
 'glucose_glycolysis':r'Hexokinase|glucose transport|GLUT1|GAPDH|glycolysis',
 'mitochondria_OXPHOS':r'mitochond|electron transport|OXPHOS',
 'activin_TGFbeta_SMAD':r'TGF.?[βb]|TGF beta|Activin|SMAD|ALK4|ALK5|ALK7',
 'ABCB1_MATE_transport':r'ABCB1|P.glycoprotein|MATE1|SLC47A1|multidrug resistance',
 'metal_iron_copper_homeostasis':r'iron chelat|copper chelat|metal chelat|ferroxidase',
 'transcription_CDK7_CDK9':r'CDK7|CDK9|transcription inhibitor',
}

primarymeta=pd.read_csv(PRISM/'primary_treatment_info.csv')
meta=primarymeta[['name','moa','target']].drop_duplicates('name').set_index('name')
bycat=mapping[mapping.source_dataset.eq('PRISM_primary')].groupby('catalog_row_id').source_compound_name.apply(list).to_dict()
classrows=[]
for i,r in lib.iterrows():
    names=bycat.get(r.catalog_row_id,[])
    broad_moa='; '.join(str(meta.at[x,'moa']) for x in names if x in meta.index and pd.notna(meta.at[x,'moa']))
    broad_target='; '.join(str(meta.at[x,'target']) for x in names if x in meta.index and pd.notna(meta.at[x,'target']))
    text=' '.join(str(r.get(c,'')) for c in ['Target','Information','Pathway'])+' '+broad_moa
    found=[class_name for class_name,regex in patterns.items() if re.search(regex,text,re.I)]
    lib.at[i,'annotation_retrieval_classes']='; '.join(found)
    lib.at[i,'Broad_MOA_annotations']=broad_moa
    lib.at[i,'Broad_target_annotations']=broad_target
    for class_name in found:
        classrows.append({'catalog_row_id':r.catalog_row_id,'CatalogNumber':r.CatalogNumber,'library_name':r['Item Name'],'retrieval_class':class_name,'classification_method':'keyword retrieval in vendor Target/Information/Pathway and Broad MOA; not gene-specific drug evidence'})

# Coverage is counted by actual matched numerical records, not mere metadata.
coverage={}
for source,df in joined_sources.items():
    subset=df[df.depmap_id.eq('ACH-000849')]
    valid_ids=set(subset.catalog_row_id)
    lib[source+'_MDA_MB_468_numeric_record']=lib.catalog_row_id.isin(valid_ids)
    exact_ids=set(subset.loc[subset.identity_mapping_method.eq('direct_normalized_name'),'catalog_row_id'])
    lib[source+'_MDA_MB_468_direct_name']=lib.catalog_row_id.isin(exact_ids)
    coverage[source]={'catalog_compounds_with_actual_468_records':len(valid_ids),'direct_normalized_name_matches':len(exact_ids),'structure_validated_matches':0}

lib.to_csv(P/'whole_catalog_audit.csv',index=False)
pd.DataFrame(classrows).to_csv(P/'catalog_mechanism_class_inventory.csv',index=False)
(P/'expanded_coverage_summary.json').write_text(json.dumps({'catalog_rows':len(lib),'source_coverage':coverage,'chemical_identity_note':'All joins are name/synonym-based. No canonical chemical structures or catalog salt forms validated. Provisional matches are not counted as chemically confirmed.','no_expression_association':'No matched cell-line expression matrix was analyzed.','no_combined_metric':'PRISM single-dose LFC, PRISM AUC and GDSC AUC are retained separately; no synergy prediction or combined score.'},indent=2))
print(json.dumps(coverage,indent=2))
