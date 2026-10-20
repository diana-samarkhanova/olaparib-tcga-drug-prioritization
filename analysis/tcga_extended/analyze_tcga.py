#!/usr/bin/env python3
"""Reproducible targeted TCGA-BRCA analysis; data snapshots hashed separately."""
from pathlib import Path
import json, warnings, hashlib, re, platform
import numpy as np
import pandas as pd
from scipy import stats
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

ROOT=Path(__file__).resolve().parent; RAW=ROOT/'raw'; FIG=ROOT/'figures'; FIG.mkdir(exist_ok=True)
KO=['LAMTOR2','WDR24','DERL2','NME6','TP53BP1','PARG']
ACT=['TP53I11','FOXQ1','CERK','PHF19','SLC35D1','PHLDA2','WT1','HEPHL1','LINC03007','SLC47A1','ELL2','INHBA','PHGDH','ABCB1']
GENES=KO+ACT
PATHS=['MTORC1_SIGNALING','UNFOLDED_PROTEIN_RESPONSE','OXIDATIVE_PHOSPHORYLATION','GLYCOLYSIS','DNA_REPAIR','G2M_CHECKPOINT','TGF_BETA_SIGNALING','EPITHELIAL_MESENCHYMAL_TRANSITION','APOPTOSIS','HYPOXIA']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
def bh(p):
 p=np.asarray(p,float); out=np.full(len(p),np.nan);ok=np.isfinite(p);v=p[ok];order=np.argsort(v); n=len(v)
 if n:out[np.where(ok)[0][order]]=np.minimum(1,np.minimum.accumulate((v[order]*n/np.arange(1,n+1))[::-1])[::-1])
 return out
def readclinical(f):return pd.read_csv(RAW/f,sep='\t',comment='#',low_memory=False)
def save(d,f):d.to_csv(ROOT/f,index=False)
sample=readclinical('data_clinical_sample.txt');patient=readclinical('data_clinical_patient.txt');rec=readclinical('data_clinical_firehose_patient.txt')
sample=sample.loc[sample.SAMPLE_TYPE.eq('Primary')].sort_values('SAMPLE_ID').drop_duplicates('PATIENT_ID')
clinical=sample.merge(patient,on='PATIENT_ID',validate='one_to_one').merge(rec[['PATIENT_ID','ER_STATUS_BY_IHC','PR_STATUS_BY_IHC','IHC_HER2','HER2_FISH_STATUS']],on='PATIENT_ID',how='left',validate='one_to_one')
def finalher2(r):
 a,b=r.IHC_HER2,r.HER2_FISH_STATUS
 if a in ['Negative','Positive'] and b in ['Negative','Positive'] and a!=b:return 'Discordant'
 if b in ['Negative','Positive']:return b
 if a in ['Negative','Positive']:return a
 return 'Unresolved'
clinical['HER2_FINAL']=clinical.apply(finalher2,axis=1)
clinical['TNBC_STATUS']='Unclassified'
negative=clinical.ER_STATUS_BY_IHC.eq('Negative')&clinical.PR_STATUS_BY_IHC.eq('Negative')&clinical.HER2_FINAL.eq('Negative')
positive=clinical.ER_STATUS_BY_IHC.eq('Positive')|clinical.PR_STATUS_BY_IHC.eq('Positive')|clinical.HER2_FINAL.eq('Positive')
clinical.loc[negative,'TNBC_STATUS']='TNBC';clinical.loc[positive,'TNBC_STATUS']='non-TNBC'
clinical['BASAL_STATUS']=clinical.SUBTYPE.map(lambda x:'Basal' if x=='BRCA_Basal' else ('Other PAM50' if isinstance(x,str) else 'Unclassified'))
cdr=pd.read_excel(RAW/'TCGA-CDR-SupplementalTableS1.xlsx',sheet_name='TCGA-CDR')
cdr=cdr.loc[cdr.type.eq('BRCA')].rename(columns={'bcr_patient_barcode':'PATIENT_ID','age_at_initial_pathologic_diagnosis':'CDR_AGE','ajcc_pathologic_tumor_stage':'CDR_STAGE'})
clinical=clinical.merge(cdr[['PATIENT_ID','CDR_AGE','CDR_STAGE','PFI','PFI.time','DFI','DFI.time','OS','OS.time','Redaction']],on='PATIENT_ID',how='left',validate='one_to_one')
def stage(s):
 s=str(s).upper()
 if re.match(r'^STAGE IV',s):return 1.
 if re.match(r'^STAGE III',s):return 1.
 if re.match(r'^STAGE II',s) or re.match(r'^STAGE I($|[AB])',s):return 0.
 return np.nan
clinical['STAGE_ADVANCED']=clinical.CDR_STAGE.map(stage)
purity_input=RAW/'purity.tsv'
if not purity_input.exists():
 purity_input=RAW/'purity_archive_derived.tsv'
 if not purity_input.exists():
  raise FileNotFoundError('Tumour purity is required to reproduce the reported adjusted analysis. Retrieve raw/purity.tsv with download_sources.py, or restore archived matched_clinical_cohort.csv and run reconstruct_archived_purity.py before rerunning analyze_tcga.py.')
 warnings.warn('Using explicitly derived archive purity input: '+str(purity_input),RuntimeWarning)
pu=pd.read_csv(purity_input,sep='\t');pc=[x for x in pu.columns if x.lower()=='purity'];sc=[x for x in pu.columns if x.lower()=='sample']
if len(pc)!=1 or len(sc)!=1:
 raise ValueError('Purity input must contain exactly one Sample and one Purity column (case-insensitive).')
pu['SAMPLE_ID']=pu[sc[0]].astype(str).str[:15];pu=pu.sort_values(sc[0]).drop_duplicates('SAMPLE_ID')
clinical=clinical.merge(pu[['SAMPLE_ID',pc[0]]].rename(columns={pc[0]:'PURITY'}),on='SAMPLE_ID',how='left',validate='one_to_one')
if clinical.PURITY.notna().sum()==0:
 raise ValueError('No purity values matched the clinical cohort; refusing to silently omit the adjustment.')
rna_input=RAW/'data_mrna_seq_v2_rsem.txt'
if not rna_input.exists():rna_input=RAW/'RNA_required_gene_subset.tsv.gz'
rna=pd.read_csv(rna_input,sep='\t',low_memory=False)
rna['Hugo_Symbol']=rna.Hugo_Symbol.fillna('')
alias={g:(g if g in set(rna.Hugo_Symbol) else ('LOC646588' if g=='LINC03007' and 'LOC646588' in set(rna.Hugo_Symbol) else ('ROBLD3' if g=='LAMTOR2' and 'ROBLD3' in set(rna.Hugo_Symbol) else None))) for g in GENES}
maps=[]
for g,s in alias.items():
 rows=rna.loc[rna.Hugo_Symbol.eq(s)] if s else rna.iloc[0:0]
 maps.append({'gene':g,'source_symbol':s,'rna_rows':len(rows),'entrez_gene_id':str(rows.iloc[0].Entrez_Gene_Id) if len(rows)==1 else '', 'rna_coverage':'available' if len(rows)==1 else 'unavailable_or_ambiguous'})
save(pd.DataFrame(maps),'gene_coverage.csv')
rna=rna.drop_duplicates('Hugo_Symbol').set_index('Hugo_Symbol').drop(columns='Entrez_Gene_Id').apply(pd.to_numeric,errors='coerce')
rna=rna.loc[:,[s for s in clinical.SAMPLE_ID if s in rna.columns]]
clinical=clinical.loc[clinical.SAMPLE_ID.isin(rna.columns)].set_index('SAMPLE_ID').loc[rna.columns]
logrna=np.log2(rna+1)
expr=pd.DataFrame({g:logrna.loc[s] for g,s in alias.items() if s in logrna.index})
expr.index.name='SAMPLE_ID';expr.to_csv(ROOT/'selected_expression_log2_RSEM_plus1.csv')
sets={}
for line in (RAW/'h.all.v2025.1.Hs.symbols.gmt').read_text().splitlines():
 a=line.split('\t');sets[a[0].removeprefix('HALLMARK_')]=a[2:]
zrna=logrna.sub(logrna.mean(axis=1),axis=0).div(logrna.std(axis=1,ddof=1).replace(0,np.nan),axis=0)
scores=pd.DataFrame({p:zrna.loc[[g for g in sets[p] if g in zrna.index]].mean(axis=0) for p in PATHS+['E2F_TARGETS']})
scores.to_csv(ROOT/'hallmark_mean_z_scores.csv')
pathcov=pd.DataFrame([{'pathway':p,'source_genes':len(sets[p]),'rna_matched_genes':sum(g in zrna.index for g in sets[p]),'source_url':'https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_'+p} for p in PATHS+['E2F_TARGETS']]);save(pathcov,'hallmark_coverage.csv')
clinical.to_csv(ROOT/'matched_clinical_cohort.csv',index_label='SAMPLE_ID')
counts={'clinical_primary_patients':len(sample),'rna_matched_primary_patients':len(clinical),'TNBC_status':clinical.TNBC_STATUS.value_counts().to_dict(),'PAM50':clinical.SUBTYPE.fillna('Unclassified').value_counts().to_dict(),'HER2_final':clinical.HER2_FINAL.value_counts().to_dict(),'purity_available':int(clinical.PURITY.notna().sum()),'gene_coverage':maps,'endpoint_counts':[]}
for cohort,mask in [('BRCA',clinical.index==clinical.index),('TNBC',clinical.TNBC_STATUS.eq('TNBC')),('Basal',clinical.BASAL_STATUS.eq('Basal'))]:
 for ep in ['PFI','DFI','OS']:
  d=clinical.loc[mask,[ep,ep+'.time']].dropna();d=d.loc[d[ep+'.time']>0]
  counts['endpoint_counts'].append({'cohort':cohort,'endpoint':ep,'n':len(d),'events':int(d[ep].sum())})
(ROOT/'cohort_counts.json').write_text(json.dumps(counts,indent=2))
print('COHORT_COUNTS',json.dumps(counts),flush=True)

# Two prespecified expression-contrast families, each corrected over available panel genes.
contrasts=[]
for comparison,col,a,b in [('clinical_TNBC_vs_nonTNBC','TNBC_STATUS','TNBC','non-TNBC'),('PAM50_Basal_vs_other','BASAL_STATUS','Basal','Other PAM50')]:
 for g in GENES:
  if g not in expr:continue
  x=expr.loc[clinical[col].eq(a),g].dropna();y=expr.loc[clinical[col].eq(b),g].dropna();u,p=stats.mannwhitneyu(x,y,alternative='two-sided',method='asymptotic')
  contrasts.append({'comparison':comparison,'gene':g,'panel':'KO' if g in KO else 'CRISPRa','n_group1':len(x),'n_group2':len(y),'median_group1':x.median(),'median_group2':y.median(),'median_delta_log2_RSEM_plus1':x.median()-y.median(),'rank_biserial':2*u/(len(x)*len(y))-1,'p_value':p})
contrasts=pd.DataFrame(contrasts);contrasts['BH_FDR']=contrasts.groupby('comparison').p_value.transform(bh);save(contrasts,'expression_contrasts.csv')

# GISTIC discrete CNA prevalence and dosage-expression association.
cna_input=RAW/'data_cna.txt'
if not cna_input.exists():cna_input=RAW/'CNA_selected_gene_subset.tsv.gz'
cna_raw=pd.read_csv(cna_input,sep='\t',low_memory=False)
cna_alias={g:(g if g in set(cna_raw.Hugo_Symbol) else (str(cna_raw.loc[pd.to_numeric(cna_raw.Entrez_Gene_Id,errors='coerce').eq(int(next(x['entrez_gene_id'] for x in maps if x['gene']==g and x['entrez_gene_id']))),'Hugo_Symbol'].iloc[0]) if any(x['gene']==g and x['entrez_gene_id'] for x in maps) and pd.to_numeric(cna_raw.Entrez_Gene_Id,errors='coerce').eq(int(next(x['entrez_gene_id'] for x in maps if x['gene']==g and x['entrez_gene_id']))).sum()==1 else None)) for g in GENES}
cna=cna_raw.drop_duplicates('Hugo_Symbol').set_index('Hugo_Symbol').drop(columns='Entrez_Gene_Id').apply(pd.to_numeric,errors='coerce')
cnarows=[]; cnamat=pd.DataFrame(index=clinical.index)
for g,s in cna_alias.items():
 if s not in cna.index:continue
 cnamat[g]=cna.loc[s].reindex(clinical.index)
 for cohort,mask in [('BRCA',clinical.index==clinical.index),('TNBC',clinical.TNBC_STATUS.eq('TNBC')),('Basal',clinical.BASAL_STATUS.eq('Basal'))]:
  v=cnamat.loc[mask,g].dropna();m=pd.concat([v,expr[g]],axis=1,keys=['cna','expr']).dropna();r,p=stats.spearmanr(m.cna,m.expr)
  cnarows.append({'cohort':cohort,'gene':g,'n_profiled':len(v),'deep_deletion_n':int((v==-2).sum()),'shallow_loss_n':int((v==-1).sum()),'diploid_n':int((v==0).sum()),'gain_n':int((v==1).sum()),'amplification_n':int((v==2).sum()),'dosage_expression_n':len(m),'dosage_expression_spearman_rho':r,'p_value':p})
cnamat.to_csv(ROOT/'selected_CNA_GISTIC.csv',index_label='SAMPLE_ID');cnarows=pd.DataFrame(cnarows);cnarows['BH_FDR']=cnarows.groupby('cohort').p_value.transform(bh);save(cnarows,'CNA_prevalence_and_expression.csv')
coverage=pd.read_csv(ROOT/'gene_coverage.csv');coverage['CNA_source_symbol']=coverage.gene.map(cna_alias);coverage['CNA_coverage']=coverage.gene.map(lambda g:'available' if g in cnamat else 'unavailable');save(coverage,'gene_coverage.csv')

# Raw and covariate-adjusted rank associations. These are co-expression, not KO dependencies.
def partialrank(x,y,cov):
 d=pd.concat([x.rename('x'),y.rename('y'),cov],axis=1).dropna(); n=len(d)
 if n<30:return n,np.nan,np.nan
 a=np.column_stack([np.ones(n),d.iloc[:,2:].to_numpy(float)]); k=np.linalg.matrix_rank(a)-1
 xr=stats.rankdata(d.x);yr=stats.rankdata(d.y)
 xr=xr-a@np.linalg.lstsq(a,xr,rcond=None)[0];yr=yr-a@np.linalg.lstsq(a,yr,rcond=None)[0]
 r=np.corrcoef(xr,yr)[0,1];df=n-k-2;p=2*stats.t.sf(abs(r*np.sqrt(df/max(1-r*r,1e-15))),df)
 return n,r,p
assoc=[]
for cohort,mask in [('BRCA',clinical.index==clinical.index),('TNBC',clinical.TNBC_STATUS.eq('TNBC')),('Basal',clinical.BASAL_STATUS.eq('Basal'))]:
 idx=clinical.index[mask]
 for g in expr:
  cov=pd.DataFrame({'age':stats.rankdata(clinical.loc[idx,'CDR_AGE'],nan_policy='omit'),'E2F':stats.rankdata(scores.loc[idx,'E2F_TARGETS'],nan_policy='omit')},index=idx)
  if clinical.loc[idx,'PURITY'].notna().sum()>=30:cov['purity']=clinical.loc[idx,'PURITY'].rank()
  if g in cnamat:cov['target_CNA']=cnamat.loc[idx,g].rank()
  if cohort=='BRCA':cov=cov.join(pd.get_dummies(clinical.loc[idx,'SUBTYPE'].fillna('Unclassified'),drop_first=True,dtype=float))
  for p in PATHS:
   members=[x for x in sets[p] if x in zrna.index and x!=alias[g]];y=zrna.loc[members,idx].mean(axis=0);x=expr.loc[idx,g]
   raw=pd.concat([x,y],axis=1).dropna();r,pv=stats.spearmanr(raw.iloc[:,0],raw.iloc[:,1]);n,pr,pp=partialrank(x,y,cov)
   assoc.append({'cohort':cohort,'gene':g,'pathway':p,'signature_gene_n':len(members),'raw_n':len(raw),'raw_spearman_rho':r,'raw_p_value':pv,'adjusted_n':n,'partial_rank_rho':pr,'adjusted_p_value':pp,'covariates':';'.join(cov.columns),'target_excluded_from_score':True})
assoc=pd.DataFrame(assoc);assoc['raw_BH_FDR']=assoc.groupby('cohort').raw_p_value.transform(bh);assoc['adjusted_BH_FDR']=assoc.groupby('cohort').adjusted_p_value.transform(bh);save(assoc,'pathway_associations.csv')

rppa=pd.read_csv(RAW/'data_rppa.txt',sep='\t').set_index('Composite.Element.REF').apply(pd.to_numeric,errors='coerce')
rpnames=['RPS6KB1|p70S6K_pT389','RPS6|S6_pS235_S236','RPS6|S6_pS240_S244','EIF4EBP1|4E-BP1_pT37T46','MAPK1 MAPK3|MAPK_pT202_Y204','AKT1 AKT2 AKT3|Akt_pS473','TP53BP1|53BP1','RAD51|Rad51','XBP1|XBP1']
rprows=[]
for cohort,mask in [('BRCA',clinical.index==clinical.index),('TNBC',clinical.TNBC_STATUS.eq('TNBC')),('Basal',clinical.BASAL_STATUS.eq('Basal'))]:
 idx=clinical.index[mask]
 for g in expr:
  for rp in rpnames:
   if rp not in rppa.index:continue
   y=rppa.loc[rp].reindex(idx);d=pd.concat([expr.loc[idx,g],y],axis=1).dropna()
   if len(d)<15:continue
   r,p=stats.spearmanr(d.iloc[:,0],d.iloc[:,1]);rprows.append({'cohort':cohort,'gene':g,'RPPA_readout':rp,'n':len(d),'spearman_rho':r,'p_value':p})
rprows=pd.DataFrame(rprows);rprows['BH_FDR']=rprows.groupby('cohort').p_value.transform(bh);save(rprows,'RPPA_associations.csv')

# Continuous-expression Cox; no optimized cut points. PFI primary, DFI supporting, OS exploratory.
surv=[]
for cohort,mask in [('BRCA',clinical.index==clinical.index),('TNBC',clinical.TNBC_STATUS.eq('TNBC')),('Basal',clinical.BASAL_STATUS.eq('Basal'))]:
 idx=clinical.index[mask]
 for ep in ['PFI','DFI','OS']:
  for g in expr:
   base=pd.DataFrame({'time_days':clinical.loc[idx,ep+'.time'],'event':clinical.loc[idx,ep],'gene_z':(expr.loc[idx,g]-expr.loc[idx,g].mean())/expr.loc[idx,g].std(ddof=1),'age_per10':clinical.loc[idx,'CDR_AGE']/10,'stage_III_IV':clinical.loc[idx,'STAGE_ADVANCED']})
   for model,cols,min_events in [('univariable',['gene_z'],15),('age_stage_adjusted',['gene_z','age_per10','stage_III_IV'],30)]:
    d=base[['time_days','event']+cols].dropna();d=d.loc[d.time_days>0];ev=int(d.event.sum());row={'cohort':cohort,'endpoint':ep,'gene':g,'model':model,'n':len(d),'events':ev,'HR_per_expression_SD':np.nan,'CI95_lower':np.nan,'CI95_upper':np.nan,'p_value':np.nan,'PH_rank_p_value':np.nan,'status':'insufficient_events' if ev<min_events else 'pending'}
    if ev>=min_events:
     try:
      with warnings.catch_warnings():
       warnings.simplefilter('ignore');fit=CoxPHFitter().fit(d,duration_col='time_days',event_col='event');s=fit.summary.loc['gene_z'];ph=proportional_hazard_test(fit,d,time_transform='rank').summary.loc['gene_z','p']
      row.update(HR_per_expression_SD=float(s['exp(coef)']),CI95_lower=float(s['exp(coef) lower 95%']),CI95_upper=float(s['exp(coef) upper 95%']),p_value=float(s.p),PH_rank_p_value=float(ph),status='fitted')
     except Exception as e:row['status']='fit_failed: '+str(e)[:130]
    surv.append(row)
surv=pd.DataFrame(surv);surv['BH_FDR']=surv.groupby(['cohort','endpoint','model']).p_value.transform(bh);surv['PH_BH_FDR']=surv.groupby(['cohort','endpoint','model']).PH_rank_p_value.transform(bh);save(surv,'survival_continuous_Cox.csv')

def figsave(f,n):f.savefig(FIG/(n+'.png'),dpi=300,bbox_inches='tight');plt.close(f)
# F1 contrasts, independent clinical and intrinsic-subtype definitions.
f,ax=plt.subplots(1,2,figsize=(11,8.5),sharey=True)
for a,comp,title in zip(ax,['clinical_TNBC_vs_nonTNBC','PAM50_Basal_vs_other'],['Clinical TNBC versus receptor-positive tumours','PAM50 basal versus other assigned subtypes']):
 d=contrasts.loc[contrasts.comparison.eq(comp)].set_index('gene').reindex(GENES);a.axvline(0,color='.5',lw=.8)
 for i,g in enumerate(GENES):
  if g not in expr:
   a.text(0,i,'NA',color='.5',ha='center',va='center');continue
  v=d.loc[g,'median_delta_log2_RSEM_plus1'];q=d.loc[g,'BH_FDR'];a.barh(i,v,color='#A44C45' if g in KO else '#527C91',height=.7);a.text(v+(.04 if v>=0 else -.04),i,'*' if q<.05 else '',ha='left' if v>=0 else 'right',va='center')
 a.set_title(title,fontsize=11);a.set_xlabel('Difference in median log₂(RSEM + 1)');a.grid(axis='x',alpha=.15)
ax[0].set_yticks(range(len(GENES)),GENES);ax[0].invert_yaxis();f.suptitle('TCGA-BRCA: candidate expression context',fontsize=15);f.text(.02,.015,'* BH FDR < 0.05 across the available candidate panel within each contrast. Missing LINC03007 is not plotted.\nDifferences are transformed-expression contrasts, not log₂ fold changes; tumour associations do not establish drug response.',fontsize=9);f.tight_layout(rect=[0,.07,1,.95]);figsave(f,'TCGA_F01_expression_contrasts')
# F2 CNA prevalence.
f,axs=plt.subplots(1,2,figsize=(12,8.5),sharey=True)
for a,cohort in zip(axs,['BRCA','TNBC']):
 d=cnarows.loc[cnarows.cohort.eq(cohort)].set_index('gene').reindex(GENES);left=np.zeros(len(GENES))
 for field,label,color in [('deep_deletion_n','Deep deletion','#35628A'),('shallow_loss_n','Shallow loss','#A3BED5'),('diploid_n','Diploid','#E7E7E7'),('gain_n','Gain','#E8B4A5'),('amplification_n','Amplification','#A7473B')]:
  v=(100*d[field]/d.n_profiled).fillna(0).to_numpy();a.barh(range(len(GENES)),v,left=left,color=color,label=label,height=.7);left+=v
 a.set_xlim(0,100);a.set_title(cohort);a.set_xlabel('Percentage of CNA-profiled tumours')
 for i,g in enumerate(GENES):
  if g not in cnamat:a.text(50,i,'NA',color='.5',ha='center',va='center')
axs[0].set_yticks(range(len(GENES)),GENES);axs[0].invert_yaxis();axs[1].legend(loc='upper left',bbox_to_anchor=(1,1));f.suptitle('TCGA-BRCA: copy-number context',fontsize=15);f.text(.03,.015,'GISTIC −2/−1/0/+1/+2. Shallow loss does not establish biallelic knockout; amplification does not establish CRISPRa-equivalent activity.',fontsize=9);f.tight_layout(rect=[0,.06,1,.95]);figsave(f,'TCGA_F02_CNA_prevalence')
# F3 prespecified adjusted within-TNBC co-expression.
f,a=plt.subplots(figsize=(13,8.5));d=assoc.loc[assoc.cohort.eq('TNBC')];mat=d.pivot(index='gene',columns='pathway',values='partial_rank_rho').reindex(index=GENES,columns=PATHS);q=d.pivot(index='gene',columns='pathway',values='adjusted_BH_FDR').reindex(index=GENES,columns=PATHS)
im=a.imshow(mat,aspect='auto',cmap='RdBu_r',norm=TwoSlopeNorm(vmin=-.75,vcenter=0,vmax=.75));a.set_yticks(range(len(GENES)),GENES);a.set_xticks(range(len(PATHS)),[p.replace('_',' ').title() for p in PATHS],rotation=45,ha='right')
for i in range(len(GENES)):
 for j in range(len(PATHS)):
  if pd.notna(mat.iloc[i,j]):a.text(j,i,('*' if q.iloc[i,j]<.05 else ''),ha='center',va='center',fontsize=13)
f.colorbar(im,ax=a,label='Partial rank correlation');a.set_title('Clinical TNBC: candidate–pathway co-expression after covariate adjustment',pad=16);f.text(.02,.005,'Adjusted for age, E2F proliferation score, purity where covered, and target CNA. Candidate excluded from its pathway score.\n* BH FDR < 0.05 across all gene–pathway tests in TNBC. These observational correlations are not synthetic-lethal dependencies.',fontsize=9);f.tight_layout(rect=[0,.07,1,1]);figsave(f,'TCGA_F03_TNBC_pathway_associations')
# F4 PFI forest plots, adjusted where sufficient event count.
f,axs=plt.subplots(1,2,figsize=(12,9),sharey=True)
for a,cohort in zip(axs,['BRCA','TNBC']):
 d=surv.loc[surv.cohort.eq(cohort)&surv.endpoint.eq('PFI')];chosen=[]
 for i,g in enumerate(GENES):
  sub=d.loc[d.gene.eq(g)&d.status.eq('fitted')];s=sub.loc[sub.model.eq('age_stage_adjusted')]
  if s.empty:s=sub.loc[sub.model.eq('univariable')]
  if s.empty:continue
  s=s.iloc[0];col='#A44C45' if g in KO else '#527C91';a.errorbar(s.HR_per_expression_SD,i,xerr=[[s.HR_per_expression_SD-s.CI95_lower],[s.CI95_upper-s.HR_per_expression_SD]],fmt='o',color=col,capsize=2);a.text(3.65,i,'*' if s.BH_FDR<.05 else '',va='center')
  chosen.append(s)
 a.axvline(1,color='.5',lw=.8);a.set_xscale('log');a.set_xlim(.35,4);a.set_xlabel('PFI hazard ratio per expression SD (95% CI)');a.set_title(cohort+(' — age/stage adjusted' if any(x.model=='age_stage_adjusted' for x in chosen) else ' — univariable'))
axs[0].set_yticks(range(len(GENES)),GENES);axs[0].invert_yaxis();f.suptitle('TCGA-CDR: continuous-expression prognosis',fontsize=15);f.text(.02,.015,'* BH FDR < 0.05 within cohort/endpoint/model. No optimized cut points. PFI is not an olaparib-specific outcome.\nModel sample/event counts and proportional-hazards checks are provided in the supplementary table.',fontsize=9);f.tight_layout(rect=[0,.07,1,.95]);figsave(f,'TCGA_F04_PFI_forest')
# F5 selected RPPA correlations.
f,a=plt.subplots(figsize=(12,8.5));d=rprows.loc[rprows.cohort.eq('TNBC')];mat=d.pivot(index='gene',columns='RPPA_readout',values='spearman_rho').reindex(index=GENES,columns=rpnames);q=d.pivot(index='gene',columns='RPPA_readout',values='BH_FDR').reindex(index=GENES,columns=rpnames)
im=a.imshow(mat,aspect='auto',cmap='RdBu_r',norm=TwoSlopeNorm(vmin=-.6,vcenter=0,vmax=.6));a.set_yticks(range(len(GENES)),GENES);a.set_xticks(range(len(rpnames)),[x.split('|')[1] for x in rpnames],rotation=45,ha='right')
for i in range(len(GENES)):
 for j in range(len(rpnames)):
  if pd.notna(mat.iloc[i,j]) and q.iloc[i,j]<.05:a.text(j,i,'*',ha='center',va='center',fontsize=13)
f.colorbar(im,ax=a,label='Spearman correlation');a.set_title('Clinical TNBC: RNA associations with RPPA protein/phosphoprotein readouts',pad=16);f.text(.02,.015,'* BH FDR < 0.05 across available gene–readout pairs within TNBC. Correlations are unadjusted and exploratory.\nPhosphoprotein abundance is not normalized to total protein here; no intervention or drug response is measured.',fontsize=9);f.tight_layout(rect=[0,.08,1,1]);figsave(f,'TCGA_F05_TNBC_RPPA_associations')

summary={'expression_significant':contrasts.loc[contrasts.BH_FDR<.05].to_dict('records'),'PFI_significant':surv.loc[surv.endpoint.eq('PFI')&surv.BH_FDR.lt(.05)].to_dict('records'),'TNBC_pathway_significant':assoc.loc[assoc.cohort.eq('TNBC')&assoc.adjusted_BH_FDR.lt(.05)].to_dict('records'),'TNBC_RPPA_significant':rprows.loc[rprows.cohort.eq('TNBC')&rprows.BH_FDR.lt(.05)].to_dict('records')}
(ROOT/'analysis_summary.json').write_text(json.dumps(summary,indent=2,default=str))
print('FINISHED',json.dumps({'expression_sig':len(summary['expression_significant']),'PFI_sig':len(summary['PFI_significant']),'TNBC_pathway_sig':len(summary['TNBC_pathway_significant']),'TNBC_RPPA_sig':len(summary['TNBC_RPPA_significant'])}),flush=True)
