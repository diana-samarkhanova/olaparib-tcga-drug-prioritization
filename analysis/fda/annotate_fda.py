from pathlib import Path
import pandas as pd
import re, json

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'analysis/fda'
lib=pd.read_excel(ROOT/'upload/L1021-DiscoveryProbe-FDA-approved-Drug-Library (1).xlsx',sheet_name='Chemical Data')
lib=lib.loc[lib['Item Name'].notna()].copy()
lib['CatalogNumber']=lib['CatalogNumber'].str.strip()
def norm(s):
    s=str(s).lower().strip().replace('α','alpha').replace('β','beta')
    s=re.sub(r'\([^)]*\)',' ',s)
    s=re.sub(r'\b(?:\d*hcl|hydrochloride|hydrobromide|hbr|mesylate|tosylate|ditosylate|besylate|maleate|tartrate|citrate|succinate|fumarate|phosphate|diphosphate|sulfate|sulphate|lactate|sodium|disodium|trisodium|potassium|calcium|magnesium|dihydrate|trihydrate|monohydrate|hydrate|anhydrous|free base|free acid|salt)\b',' ',s)
    return re.sub(r'[^a-z0-9]','',s)
def aliases(s):
    al=[re.sub(r'\([^)]*\)',' ',str(s))]
    for inner in re.findall(r'\(([^)]*)\)',str(s)):
        al+=re.split(r'[,;/]',inner)
    return {norm(a) for a in al if norm(a)}
products=pd.read_csv(OUT/'Products.txt',sep='\t',encoding='cp1252',low_memory=False)
marketing=pd.read_csv(OUT/'MarketingStatus.txt',sep='\t',encoding='cp1252')
lookup=pd.read_csv(OUT/'MarketingStatus_Lookup.txt',sep='\t',encoding='cp1252')
sub=pd.read_csv(OUT/'Submissions.txt',sep='\t',encoding='cp1252',low_memory=False)
approved_apps=set(sub.loc[sub['SubmissionStatus'].eq('AP'),'ApplNo'])
products=products.loc[products['ApplNo'].isin(approved_apps)].merge(marketing,on=['ApplNo','ProductNo'],how='left').merge(lookup,on='MarketingStatusID',how='left')
ingredient_map={}
for _,row in products.iterrows():
    for ingredient in str(row['ActiveIngredient']).split(';'):
        k=norm(ingredient)
        if k: ingredient_map.setdefault(k,[]).append({'ingredient':ingredient.strip(),'application':int(row['ApplNo']),'status':row.get('MarketingStatusDescription')})
# Manually resolved active-moiety synonym (same chemical identity; not a spelling approximation).
manual={'B6185':'acalabrutinib','N1504':'omacetaxine mepesuccinate'}
rows=[]
for _,r in lib.iterrows():
    keys=aliases(r['Item Name'])
    if r['CatalogNumber'] in manual: keys.add(norm(manual[r['CatalogNumber']]))
    matches=[]
    for k in keys: matches.extend(ingredient_map.get(k,[]))
    current=any(m['status'] in ['Prescription','Over-the-counter'] for m in matches)
    result=dict(r)
    result.update({'FDA_record_match':bool(matches),'FDA_current_product_record':current,'FDA_ingredient_matches':'; '.join(sorted({m['ingredient'] for m in matches})),'FDA_application_numbers':'; '.join(f'{n:06d}' for n in sorted({m['application'] for m in matches})),'FDA_marketing_status':'; '.join(sorted({str(m['status']) for m in matches})),'FDA_mapping_note':'Exact normalized ingredient or declared synonym; salts collapsed. No match means unresolved, not proof of nonapproval.','FDA_source':'https://www.fda.gov/drugs/drug-approvals-and-databases/drugsfda-data-files'})
    rows.append(result)
data=pd.DataFrame(rows)
data.to_csv(OUT/'library_fda_annotation.csv',index=False)
summary={'library_compounds':len(data),'FDA_ingredient_record_matches':int(data['FDA_record_match'].sum()),'FDA_current_product_record_matches':int(data['FDA_current_product_record'].sum()),'FDA_no_match_unresolved':int((~data['FDA_record_match']).sum()),'FDA_records_downloaded':'2026-10-12','FDA_page_last_updated':'2026-10-09','interpretation':'Automated active-ingredient matching; not comprehensive regulatory adjudication. Nonmatches retained as unresolved; current marketing records do not establish approval of a combination with olaparib.'}
(OUT/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
print(json.dumps(summary,ensure_ascii=False))
patterns='Lapatinib|Exemestane|Elacridar|Zosuquidar|Panobinostat|Fludarabine|Chlorambucil|Plerixafor|Tamoxifen|Abiraterone|ACP-196|Homoharringtonine'
print(data.loc[data['Item Name'].str.contains(patterns,case=False),['CatalogNumber','Item Name','FDA_record_match','FDA_current_product_record','FDA_ingredient_matches','FDA_marketing_status']].to_string(index=False))
