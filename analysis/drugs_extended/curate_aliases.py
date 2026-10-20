from pathlib import Path
from urllib.request import Request,urlopen
from urllib.parse import quote
from concurrent.futures import ThreadPoolExecutor
import json,re
import pandas as pd

P=Path(__file__).resolve().parent
ROOT=P.parents[1]
P.mkdir(parents=True,exist_ok=True)
RAW=P/'identity_sources';RAW.mkdir(exist_ok=True)
lib=pd.read_csv(ROOT/'analysis/library/library_raw.csv')
defs=[
 ('A8221','tazemetostat','EPZ-6438','same development compound; FDA marketed hydrobromide form','All Tazverik approvals withdrawn 2026-06-22'),
 ('A8346','alpelisib','BYL-719','same development compound',''),
 ('A4007','ixazomib citrate','MLN9708','clinically marketed citrate ester/prodrug form; distinct from MLN2238',''),
 ('A4008','ixazomib','MLN2238','active boronic-acid form; FDA product contains ixazomib citrate, not the exact catalog form',''),
 ('A3617','niraparib','MK-4827','same active moiety; marketed tosylate form may differ',''),
 ('A1794','abemaciclib','LY2835219','catalog CAS differs from free base; retain formulation distinction',''),
 ('A3575','abemaciclib','LY2835219','explicit free-base catalog entry',''),
 ('A8194','venetoclax','ABT-199','same development compound',''),
 ('B4882','lorlatinib','PF-06463922','same development compound',''),
 ('B7805','ivosidenib','AG-120','same development compound',''),
 ('B2266','sonidegib','LDE225','same active moiety; marketed phosphate form may differ',''),
 ('A2664','tezacaftor','VX-661','same development compound',''),
 ('A8351','lumacaftor','VX-809','same development compound',''),
 ('A3618','grazoprevir','MK-5172','same development compound',''),
 (None,'azacitidine','5-azacytidine','common chemical synonym; catalog presence checked by exact name',''),
 ('A8255','sunitinib','sunitinib','catalog malate counterion form; active-moiety mapping to public generic compound remains provisional',''),
 ('B1045','sunitinib','sunitinib','catalog free base versus FDA marketed malate form; not an exact marketed-form match',''),
 ('A3887','trametinib','trametinib','catalog DMSO-solvate form; public free-compound response mapping is provisional',''),
 ('A8628','chloroquine','chloroquine','catalog diphosphate counterion form; public generic-compound mapping is provisional',''),
 ('B4874','hydroxychloroquine','hydroxychloroquine','catalog sulfate counterion form; public generic-compound mapping is provisional',''),
]

def query(t):
 cat,generic,code,relation,override=t
 row=lib[lib.CatalogNumber.eq(cat)]
 if row.empty and generic=='azacitidine': row=lib[lib['Item Name'].eq('5-Azacytidine')]
 if row.empty:return None
 r=row.iloc[0];cat=r.CatalogNumber
 cas=str(r['CAS Number']).strip()
 url='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+quote(code,safe='')+'/synonyms/JSON'
 result={'CatalogNumber':cat,'library_name':r['Item Name'],'library_cas':cas,'generic_alias':generic,'development_code':code,'chemical_form_relation':relation,'regulatory_override':override,'identity_source_url':url,'identity_method':'curated literature/development-code alias; structure unverified','pubchem_CIDs':'','catalog_CAS_in_synonyms':False,'pubchem_lookup_status':'not retrieved'}
 try:
  cached=RAW/(str(cat)+'_pubchem_synonyms.json')
  if cached.exists():b=cached.read_bytes()
  else:
   with urlopen(Request(url,headers={'User-Agent':'compound-identity-research/1.0'}),timeout=20) as response:b=response.read()
   cached.write_bytes(b)
  data=json.loads(b);entries=data.get('InformationList',{}).get('Information',[])
  result['pubchem_CIDs']='; '.join(str(x['CID']) for x in entries)
  result['catalog_CAS_in_synonyms']=any(cas in x.get('Synonym',[]) for x in entries)
  result['pubchem_lookup_status']='retrieved'
  if result['catalog_CAS_in_synonyms']:result['identity_method']='curated development-code alias with catalog CAS present in PubChem synonyms; structure unverified'
 except Exception as e:result['pubchem_lookup_status']=type(e).__name__+': '+str(e)[:160]
 return result

with ThreadPoolExecutor(max_workers=3) as ex:rows=[x for x in ex.map(query,defs) if x]
d=pd.DataFrame(rows);d.to_csv(P/'curated_development_code_aliases.csv',index=False)
print(d[['CatalogNumber','library_name','generic_alias','catalog_CAS_in_synonyms','pubchem_lookup_status']].to_string(index=False))
