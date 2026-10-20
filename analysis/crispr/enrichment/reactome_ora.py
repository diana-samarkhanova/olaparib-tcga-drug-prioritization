from pathlib import Path
from zipfile import ZipFile
from datetime import datetime, timezone
import hashlib
import json
import pandas as pd
import numpy as np
import scipy
from scipy.stats import hypergeom

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / 'enrichment'
MIN_SIZE = 10
MAX_SIZE = 500
ARCHIVE = OUT / 'ReactomePathways.gmt.zip'
EXPECTED_GMT_SHA256 = '8c1dbc8578431da5d2d5118262718c60b553a9be3398e93658daa069e4a9afd4'
if not ARCHIVE.is_file():
    raise FileNotFoundError(
        f'Required historical Reactome input is missing: {ARCHIVE}. '
        'The original archive omitted this ZIP. The source was '
        'https://reactome.org/download/current/ReactomePathways.gmt.zip '
        '(retrieved 2026-10-12); the current download may have changed. '
        f'Exact reproduction requires SHA256 {EXPECTED_GMT_SHA256}. '
        'See docs/CRISPR_REPRODUCIBILITY.md.'
    )
if hashlib.sha256(ARCHIVE.read_bytes()).hexdigest() != EXPECTED_GMT_SHA256:
    raise ValueError(
        'Reactome GMT differs from the historical snapshot. '
        f'Expected SHA256 {EXPECTED_GMT_SHA256}. '
        'A new database release must be analysed and documented separately.'
    )
with ZipFile(ARCHIVE) as z:
    text = z.read('ReactomePathways.gmt').decode('utf-8')
    zip_timestamp = z.getinfo('ReactomePathways.gmt').date_time

pathways = {}
for line in text.splitlines():
    parts = line.strip().split('\t')
    if len(parts)<3 or not parts[1].startswith('R-HSA-'):
        continue
    name,pid,*genes = parts
    if pid in pathways: raise ValueError('Duplicate Reactome stable ID '+pid)
    pathways[pid]={'name':name,'genes':set(g for g in genes if g)}
annotation_genes = set.union(*(p['genes'] for p in pathways.values()))

def bh(p):
    p=np.asarray(p,dtype=float)
    order=np.argsort(p)
    adjusted=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    output=np.empty_like(p)
    output[order]=np.minimum(adjusted,1)
    return output

metadata=json.loads((OUT/'download_metadata.json').read_text())
metadata.update({
    'analysis_utc':datetime.now(timezone.utc).isoformat(),
    'database_version_api_url':'https://reactome.org/ContentService/data/database/version',
    'database_version_api_response':(OUT/'reactome_version_response.txt').read_text().strip(),
    'database_release_date_verified_on_homepage':'2026-06-30',
    'database_homepage_url':'https://reactome.org/',
    'gmt_member_timestamp':list(zip_timestamp),
    'human_gene_sets_total':len(pathways),'human_gene_symbols_or_labels_union':len(annotation_genes),
    'software':{'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__},
    'method':'One-sided hypergeometric over-representation analysis, Benjamini-Hochberg adjustment independently across every eligible pathway for each library, including zero-overlap pathways. Foreground is recurrent positive-selection FDR<0.05 plus LFC>0 at both timepoints; background is all genes assayed in that library intersected by exact symbol with the human Reactome GMT annotation universe.',
    'min_pathway_size_in_assayed_annotated_background':MIN_SIZE,
    'max_pathway_size_in_assayed_annotated_background':MAX_SIZE,
    'mapping':'Exact, case-sensitive symbol matching only; no alias conversion or mapping from previous identifiers. Unmapped genes are exported and excluded from the annotated primary analysis.',
    'libraries':{},
    'limitations':['ORA identifies over-represented annotations; it does not establish pathway activation, mechanism, drug sensitivity, synergy, or direction of a drug effect.',
                   'Related Reactome pathways share genes; significant child and parent pathways are not independent mechanisms.',
                   'Annotated gene background is smaller than the entire assayed library. Annotation coverage and excluded genes are reported.',
                   'Repeated screen timepoints are not independent wet-lab validations. The small SAM recurrent foreground limits power.',
                   'Current GMT includes human pathways with some viral gene labels; these labels do not match the assayed human symbol background.']
})
results=[]
for library in ('KO','SAM'):
    table=pd.read_csv(BASE/f'{library}_all_genes_two_timepoints.csv')
    assayed=set(table['gene'])
    foreground=set(table.loc[table['recurrent_positive_direction_consistent'],'gene'])
    background=assayed & annotation_genes
    annotated_foreground=foreground & background
    n=len(annotated_foreground);N=len(background)
    mapping=pd.DataFrame({'gene':sorted(assayed)})
    mapping['is_recurrent_foreground']=mapping.gene.isin(foreground)
    mapping['is_reactome_annotated']=mapping.gene.isin(background)
    mapping.to_csv(OUT/f'{library}_gene_mapping.csv',index=False)
    rows=[]
    for pid,p in pathways.items():
        pathway_genes=p['genes'] & background
        K=len(pathway_genes)
        if not MIN_SIZE<=K<=MAX_SIZE: continue
        hit_genes=annotated_foreground & pathway_genes
        k=len(hit_genes)
        expected=n*K/N
        pval=float(hypergeom.sf(k-1,N,K,n)) if n else 1.0
        rows.append({'library':library,'pathway_id':pid,'pathway_name':p['name'],
                     'pathway_url':'https://reactome.org/content/detail/'+pid,
                     'background_assayed_genes':len(assayed),'background_annotated_genes':N,
                     'foreground_input_genes':len(foreground),'foreground_annotated_genes':n,
                     'pathway_total_gmt_genes':len(p['genes']),
                     'pathway_genes_in_annotated_background':K,'foreground_pathway_overlap':k,
                     'expected_overlap':expected,'fold_enrichment':k/expected if expected else np.nan,
                     'p_value':pval,'overlap_genes':';'.join(sorted(hit_genes))})
    df=pd.DataFrame(rows)
    df['bh_fdr']=bh(df['p_value'])
    df['fdr_below_0_05']=df.bh_fdr<0.05
    df=df.sort_values(['bh_fdr','p_value','foreground_pathway_overlap'],ascending=[True,True,False])
    df.to_csv(OUT/f'{library}_reactome_ora_all_tested.csv',index=False)
    df.loc[df.fdr_below_0_05].to_csv(OUT/f'{library}_reactome_ora_significant.csv',index=False)
    df.head(30).to_csv(OUT/f'{library}_reactome_ora_top30.csv',index=False)
    metadata['libraries'][library]={
        'assayed_genes':len(assayed),'annotated_background_genes':N,'background_annotation_percent':100*N/len(assayed),
        'foreground_input_genes':len(foreground),'foreground_annotated_genes':n,
        'foreground_annotated_gene_symbols':sorted(annotated_foreground),
        'foreground_unmapped_gene_symbols':sorted(foreground-background),
        'pathways_tested':len(df),'significant_pathways_bh_fdr_0_05':int(df.fdr_below_0_05.sum()),
        'top_results':df.head(10).to_dict(orient='records')}
    results.append(df)
    print('\n',library,metadata['libraries'][library])
    print(df.loc[df.fdr_below_0_05,['pathway_id','pathway_name','foreground_pathway_overlap','overlap_genes','p_value','bh_fdr']].to_string(index=False))

pd.concat(results,ignore_index=True).to_csv(OUT/'reactome_ora_both_libraries_all_tested.csv',index=False)
(OUT/'ora_metadata.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2))
