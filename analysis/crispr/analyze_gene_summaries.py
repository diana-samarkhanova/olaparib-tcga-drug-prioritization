from pathlib import Path
import json
import re
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'analysis/crispr'
OUT.mkdir(parents=True, exist_ok=True)
FDR_CUTOFF = 0.05
SPECS = {
    'KO_T16': ('KO_RRA_T16_OL_vs_DMSO_NEWcnv.gene_summary (1).xlsx', 'KO', 'T16'),
    'KO_T24': ('KO_RRA_T24_OL_vs_DMSO_NEWcnv.gene_summary (2).xlsx', 'KO', 'T24'),
    'SAM_T14': ('SAM_RRA_T14_OL_vs_DMSO_NEWcnv_outlier_removed.gene_summary (3).xlsx', 'SAM', 'T14'),
    'SAM_T23': ('SAM_RRA_T23_OL_vs_DMSO_NEWcnv_outlier.gene_summary (2).xlsx', 'SAM', 'T23'),
}
EXPECTED = ['id', 'num', 'neg|score', 'neg|p-value', 'neg|fdr', 'neg|rank', 'neg|goodsgrna', 'neg|lfc',
            'pos|score', 'pos|p-value', 'pos|fdr', 'pos|rank', 'pos|goodsgrna', 'pos|lfc']
TARGETS = ['PARG','TP53BP1','LAMTOR2','WDR24','ABCB1','PARP1','PARP2','ADPRHL2','PHLDA2','NANOS2',
           'TP53I11','TUBB2A','TDO2','WT1','DUSP22','SLC2A1','ACACA','EIF4A1','MLST8','GNB1','GNB2','GNG5']

summary = {
    'fdr_cutoff': FDR_CUTOFF,
    'method': 'Read one complete 14-column MAGeCK RRA gene-summary worksheet per comparison; classify FDR < 0.05 separately for positive and negative selection. Direction-consistent hits additionally require pos LFC > 0 or neg LFC < 0. No reanalysis of counts, recomputation of FDR, or pooling of time-point P-values.',
    'limitations': ['MAGeCK gene summaries do not retain sample identifiers, replicate structure, count matrices, or outlier-removal logs.',
                    'Recurrent calls indicate reproducibility across supplied time points, not independent biological validation.',
                    'Vendor target/pathway labels are annotations; exact labels do not establish binding specificity, activity on these cells, FDA status, synergy, or resensitization.',
                    'Inhibition of a gene whose KO is enriched under olaparib cannot be treated as a strategy to reverse resistance solely from that KO result.'],
    'comparisons': {}, 'libraries': {}, 'selected_targets': {}, 'excluded_worksheets': []
}
frames = {}
for key, (fname, library, timepoint) in SPECS.items():
    path = ROOT / 'upload' / fname
    xls = pd.ExcelFile(path)
    candidates = []
    for sheet in xls.sheet_names:
        header = pd.read_excel(path, sheet_name=sheet, nrows=0)
        if header.columns.tolist() == EXPECTED:
            candidates.append(sheet)
        else:
            summary['excluded_worksheets'].append({'comparison':key, 'sheet':sheet, 'column_count':len(header.columns), 'columns':header.columns.tolist()})
    if len(candidates) != 1:
        raise ValueError(f'{key}: expected one full sheet, found {candidates}')
    df = pd.read_excel(path, sheet_name=candidates[0])
    if df['id'].isna().any() or df['id'].duplicated().any():
        raise ValueError(f'{key}: missing/duplicated gene ID')
    for c in EXPECTED[1:]:
        df[c] = pd.to_numeric(df[c], errors='raise')
    rename = {c:c.replace('|','_').replace('-','_') for c in EXPECTED}
    df = df.rename(columns=rename)
    df['gene'] = df.pop('id').astype(str).str.strip()
    df['comparison'] = key
    df['library'] = library
    df['timepoint'] = timepoint
    df['source_file'] = fname
    df['source_sheet'] = candidates[0]
    df['source_excel_row'] = np.arange(2,len(df)+2)
    df['pos_significant'] = df['pos_fdr'] < FDR_CUTOFF
    df['neg_significant'] = df['neg_fdr'] < FDR_CUTOFF
    df['pos_direction_consistent'] = df['pos_significant'] & (df['pos_lfc'] > 0)
    df['neg_direction_consistent'] = df['neg_significant'] & (df['neg_lfc'] < 0)
    df['pos_at_least_two_goodsgrna'] = df['pos_direction_consistent'] & (df['pos_goodsgrna'] >= 2)
    df['neg_at_least_two_goodsgrna'] = df['neg_direction_consistent'] & (df['neg_goodsgrna'] >= 2)
    df['selection_class'] = np.select([
        df['pos_significant'] & df['neg_significant'],
        df['pos_direction_consistent'], df['neg_direction_consistent'],
        df['pos_significant'] | df['neg_significant']],
        ['both_directions_significant', 'positive_selection','negative_selection','significant_lfc_discordant'],default='not_significant')
    for c in ('pos_fdr','neg_fdr','pos_p_value','neg_p_value'):
        if not df[c].between(0,1).all(): raise ValueError(f'{key}: invalid {c}')
    audit = {
        'source_file':fname,'source_sheet':candidates[0], 'genes':len(df),
        'positive_fdr_hits':int(df['pos_significant'].sum()),
        'negative_fdr_hits':int(df['neg_significant'].sum()),
        'positive_direction_consistent_hits':int(df['pos_direction_consistent'].sum()),
        'negative_direction_consistent_hits':int(df['neg_direction_consistent'].sum()),
        'positive_hits_at_least_two_goodsgrna':int(df['pos_at_least_two_goodsgrna'].sum()),
        'negative_hits_at_least_two_goodsgrna':int(df['neg_at_least_two_goodsgrna'].sum()),
        'both_directions_significant':int((df['pos_significant'] & df['neg_significant']).sum()),
        'positive_fdr_hits_nonpositive_lfc':df.loc[df['pos_significant'] & ~(df['pos_lfc']>0),'gene'].tolist(),
        'negative_fdr_hits_nonnegative_lfc':df.loc[df['neg_significant'] & ~(df['neg_lfc']<0),'gene'].tolist(),
        'pos_neg_lfc_equal_all_genes':bool(np.isclose(df['pos_lfc'],df['neg_lfc']).all()),
        'num_sgrna_distribution':{str(k):int(v) for k,v in df['num'].value_counts().sort_index().items()},
        'missing_numeric_values':int(df[list(rename.values())[1:]].isna().sum().sum()),
    }
    frames[key] = df
    summary['comparisons'][key] = audit

all_long = pd.concat(frames.values(),ignore_index=True)
all_long.to_csv(OUT/'all_comparisons_genes_long.csv',index=False)
all_long.loc[all_long.pos_significant | all_long.neg_significant].to_csv(OUT/'all_fdr_hits_long.csv',index=False)

wide_frames = {}
for library,keys in {'KO':['KO_T16','KO_T24'],'SAM':['SAM_T14','SAM_T23']}.items():
    parts=[]
    for key in keys:
        d=frames[key].drop(columns=['comparison','library','timepoint']).set_index('gene')
        d.columns=[key+'_'+c for c in d.columns]
        parts.append(d)
    wide = pd.concat(parts,axis=1).reset_index()
    a,b=keys
    flags={
        'recurrent_positive_fdr':wide[a+'_pos_significant'].fillna(False) & wide[b+'_pos_significant'].fillna(False),
        'recurrent_negative_fdr':wide[a+'_neg_significant'].fillna(False) & wide[b+'_neg_significant'].fillna(False),
        'recurrent_positive_direction_consistent':wide[a+'_pos_direction_consistent'].fillna(False) & wide[b+'_pos_direction_consistent'].fillna(False),
        'recurrent_negative_direction_consistent':wide[a+'_neg_direction_consistent'].fillna(False) & wide[b+'_neg_direction_consistent'].fillna(False),
        'recurrent_positive_two_goodsgrna':wide[a+'_pos_at_least_two_goodsgrna'].fillna(False) & wide[b+'_pos_at_least_two_goodsgrna'].fillna(False),
        'recurrent_negative_two_goodsgrna':wide[a+'_neg_at_least_two_goodsgrna'].fillna(False) & wide[b+'_neg_at_least_two_goodsgrna'].fillna(False),
        'direction_switch':(wide[a+'_pos_direction_consistent'].fillna(False) & wide[b+'_neg_direction_consistent'].fillna(False)) | (wide[a+'_neg_direction_consistent'].fillna(False) & wide[b+'_pos_direction_consistent'].fillna(False)),
    }
    for name,flag in flags.items():wide[name]=flag
    wide['recurrent_class']=np.select([
        wide['direction_switch'],wide['recurrent_positive_direction_consistent'],wide['recurrent_negative_direction_consistent']],
        ['direction_switch','recurrent_positive_selection','recurrent_negative_selection'],default='not_recurrent')
    wide['best_pos_fdr']=wide[[a+'_pos_fdr',b+'_pos_fdr']].min(axis=1)
    wide['worst_pos_fdr']=wide[[a+'_pos_fdr',b+'_pos_fdr']].max(axis=1)
    wide['best_neg_fdr']=wide[[a+'_neg_fdr',b+'_neg_fdr']].min(axis=1)
    wide['worst_neg_fdr']=wide[[a+'_neg_fdr',b+'_neg_fdr']].max(axis=1)
    wide=wide.sort_values(['recurrent_positive_direction_consistent','recurrent_negative_direction_consistent','best_pos_fdr','best_neg_fdr'],ascending=[False,False,True,True])
    wide.to_csv(OUT/f'{library}_all_genes_two_timepoints.csv',index=False)
    hit_mask=wide[[key+'_'+s+'_significant' for key in keys for s in ('pos','neg')]].fillna(False).any(axis=1)
    wide.loc[hit_mask].to_csv(OUT/f'{library}_hits_any_timepoint.csv',index=False)
    recurrent_mask=wide.recurrent_positive_direction_consistent | wide.recurrent_negative_direction_consistent | wide.direction_switch
    wide.loc[recurrent_mask].to_csv(OUT/f'{library}_recurrent_hits.csv',index=False)
    wide_frames[library]=wide
    lib={'genes_union':len(wide),'all_hit_genes_union':int(hit_mask.sum())}
    for name in flags:lib[name+'_count']=int(wide[name].sum())
    lib['recurrent_positive_genes']=wide.loc[wide.recurrent_positive_direction_consistent,'gene'].tolist()
    lib['recurrent_negative_genes']=wide.loc[wide.recurrent_negative_direction_consistent,'gene'].tolist()
    lib['direction_switch_genes']=wide.loc[wide.direction_switch,'gene'].tolist()
    summary['libraries'][library]=lib

for target in TARGETS:
    rows=all_long.loc[all_long.gene==target,['comparison','num','pos_fdr','pos_lfc','pos_goodsgrna','neg_fdr','neg_lfc','neg_goodsgrna','selection_class']]
    summary['selected_targets'][target]=rows.to_dict(orient='records')
all_long.loc[all_long.gene.isin(TARGETS)].to_csv(OUT/'selected_targets_all_comparisons.csv',index=False)

# Exact vendor gene/known target-family text searches are evidence annotations only.
chem_path=ROOT/'upload/L1021-DiscoveryProbe-FDA-approved-Drug-Library (1).xlsx'
chem=pd.read_excel(chem_path,sheet_name='Chemical Data')
groups={
    'ABCB1':(r'\bABCB1\b|P[- ]?glycoprotein|\bP[- ]?gp\b|MDR1', 'ABCB1 or established synonymous target label'),
    'PARP1/PARP2':(r'\bPARP\b|PARP1|PARP2', 'PARP-family label, not isoform-specific'),
    'LAMTOR2/WDR24/MLST8':(r'\bmTOR\b|PI3K|AKT', 'Pathway-level PI3K AKT mTOR association, not binding to LAMTOR2/WDR24/MLST8'),
    'SLC2A1':(r'\bSLC2A1\b|\bGLUT1\b', 'SLC2A1 or GLUT1 label'),
    'ACACA':(r'\bACACA\b|acetyl[- ]?CoA carboxylase|\bACC1\b', 'ACACA or ACC target-family label'),
    'EIF4A1':(r'\bEIF4A1\b|\beIF4A\b', 'EIF4A target-family label'),
    'TDO2':(r'\bTDO2\b|tryptophan 2,3[- ]?dioxygenase|\bTDO\b', 'TDO target-family label'),
    'TUBB2A':(r'tubulin|microtubule', 'Tubulin target family, not TUBB2A specific'),
    'PARG':(r'\bPARG\b|poly\(ADP[- ]ribose\) glycohydrolase', 'PARG label'),
    'TP53BP1':(r'\bTP53BP1\b|\b53BP1\b', 'TP53BP1 or 53BP1 label'),
}
maprows=[]
for gene,(pattern,relationship) in groups.items():
    # Search vendor Target + Pathway fields only; Information may be misleading free text.
    targetmatch=chem['Target'].fillna('').astype(str).str.contains(pattern,case=False,regex=True)
    pathwaymatch=chem['Pathway'].fillna('').astype(str).str.contains(pattern,case=False,regex=True)
    for _,row in chem.loc[targetmatch|pathwaymatch].iterrows():
        out=row.to_dict();out['screen_gene_group']=gene;out['relationship']=relationship
        out['matched_field']='Target' if bool(re.search(pattern,str(row['Target']),re.I)) else 'Pathway'
        maprows.append(out)
pd.DataFrame(maprows).to_csv(OUT/'vendor_target_pathway_text_matches.csv',index=False)
summary['vendor_annotation_text_matches']={gene:sum(r['screen_gene_group']==gene for r in maprows) for gene in groups}
summary['vendor_annotation_rows']=len(chem)

axis_specs = [
    ('PAR metabolism and DNA damage response','KO',['PARP2','ADPRHL2','PARG'],
     'PARP2 and ADPRHL2 recurrent positive selection; PARG positive selection only T24',
     'Loss is associated with enrichment; a compound inhibiting the same protein does not by itself imply reversal of resistance.'),
    ('Lysosomal amino-acid sensing and mTOR signaling','KO',['LAMTOR2','WDR24','MLST8'],
     'Three recurrent positive-selection genes in the nutrient-sensing/mTOR axis',
     'Vendor mTOR/PI3K/AKT labels are pathway-level matches, not LAMTOR2 or WDR24 binding evidence. mTOR inhibition may phenocopy a resistance-associated change; response must be tested.'),
    ('Drug efflux','SAM',['ABCB1'],
     'Recurrent positive selection with all 3 guides supporting enrichment at both timepoints',
     'An efflux inhibitor or an alternative non-substrate can be investigated. Vendor P-gp hit is zosuquidar; approval status must be verified separately.'),
    ('Translation initiation','KO',['EIF4A1'],
     'Recurrent positive selection with all 4 guides at both timepoints',
     'Inhibitor activity on resistant cells is not demonstrated by this loss-of-function resistance hit.'),
    ('Glucose uptake and lipid synthesis','KO',['SLC2A1','ACACA'],
     'Recurrent positive selection for both genes',
     'Loss-associated resistance prevents automatic prioritization of target inhibitors as resensitizers.'),
    ('ER-associated degradation and glycosylation','KO',['DERL2','DPM1','MPDU1'],
     'Recurrent positive selection for all three genes',
     'Biological grouping only; no statistical pathway enrichment has been performed in this analysis.'),
    ('Cell-cycle transcription and nucleotide synthesis','KO',['LIN54','LIN52','E2F7','TYMS'],
     'Recurrent positive selection for all four genes',
     'Cell-cycle slowing can affect drug selection; a corresponding cytotoxic drug needs independent pharmacological evidence.'),
    ('Transcription chromatin and RNA processing','KO',['SMCHD1','SETD1B','LUC7L2','PAPOLA'],
     'Recurrent positive selection for all four genes',
     'Biological grouping only; no statistical pathway enrichment has been performed.'),
    ('Additional activation-screen candidates','SAM',['TP53I11','CMTR2'],
     'Recurrent positive selection at both timepoints',
     'No direct vendor target match identified; functional mechanisms and druggable dependencies need investigation.'),
]
axisrows=[]
for axis,library,genes,evidence,caveat in axis_specs:
    for gene in genes:
        for _,row in all_long.loc[(all_long.gene==gene)&(all_long.library==library)].iterrows():
            axisrows.append({'axis':axis,'gene':gene,'library':library,'comparison':row['comparison'],
                             'num':row['num'],'pos_fdr':row['pos_fdr'],'pos_lfc':row['pos_lfc'],
                             'pos_goodsgrna':row['pos_goodsgrna'],'selection_class':row['selection_class'],
                             'evidence_summary':evidence,'interpretation_limit':caveat,
                             'grouping_method':'Manually assigned biological categories; not a pathway-enrichment result'})
pd.DataFrame(axisrows).to_csv(OUT/'candidate_gene_axes.csv',index=False)
summary['candidate_axes']=[{'axis':a,'library':lib,'genes':genes,'evidence':ev,'interpretation_limit':cav} for a,lib,genes,ev,cav in axis_specs]

with (OUT/'summary.json').open('w') as f:json.dump(summary,f,indent=2,ensure_ascii=False)
print(json.dumps({k:summary[k] for k in ('comparisons','libraries','vendor_annotation_text_matches')},indent=2,ensure_ascii=False))
print('\nSelected targets:')
print(all_long.loc[all_long.gene.isin(TARGETS),['gene','comparison','num','pos_fdr','pos_lfc','pos_goodsgrna','neg_fdr','neg_lfc','neg_goodsgrna','selection_class']].sort_values(['gene','comparison']).to_string(index=False))
