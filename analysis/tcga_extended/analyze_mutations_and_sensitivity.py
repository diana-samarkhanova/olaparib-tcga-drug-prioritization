from pathlib import Path
import json,re,warnings,hashlib,datetime
import pandas as pd,numpy as np
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent
(R/'figures').mkdir(exist_ok=True)
C=pd.read_csv(R/'matched_clinical_cohort.csv',index_col='SAMPLE_ID');E=pd.read_csv(R/'selected_expression_log2_RSEM_plus1.csv',index_col='SAMPLE_ID')
genes=['LAMTOR2','WDR24','DERL2','NME6','TP53BP1','PARG','TP53I11','FOXQ1','CERK','PHF19','SLC35D1','PHLDA2','WT1','HEPHL1','LINC03007','SLC47A1','ELL2','INHBA','PHGDH','ABCB1']
case=(R/'raw/mutation_cases.txt').read_text();sequenced=set(re.search(r'case_list_ids:\s*(.*)',case).group(1).split())
rows=[];variants=[]
lof={'Frame_Shift_Del','Frame_Shift_Ins','Nonsense_Mutation','Splice_Site'}
for g in genes+['BRCA1','BRCA2','TP53']:
 f=R/'raw/mutations_api'/(g+'.json');available=f.exists();j=json.loads(f.read_text()) if available else []
 for x in j:
  variants.append({'gene':g,'sample_id':x.get('sampleId'),'patient_id':x.get('patientId'),'mutation_type':x.get('mutationType'),'protein_change':x.get('proteinChange'),'putative_truncating_or_splice':x.get('mutationType') in lof,'chr':x.get('chr'),'start_position':x.get('startPosition'),'reference_allele':x.get('referenceAllele'),'variant_allele':x.get('variantAllele')})
 for cohort,mask in [('BRCA',C.index==C.index),('TNBC',C.TNBC_STATUS.eq('TNBC')),('Basal',C.BASAL_STATUS.eq('Basal'))]:
  ss=set(C.index[mask])&sequenced;mutated={x.get('sampleId') for x in j}&ss;trunc={x.get('sampleId') for x in j if x.get('mutationType') in lof}&ss
  rows.append({'cohort':cohort,'gene':g,'mutation_data_status':'available_coding_WXS' if available else 'not_retrieved_or_not_coding_WXS_gene','n_sequenced':len(ss) if available else np.nan,'n_mutated_tumours':len(mutated) if available else np.nan,'mutated_percent':100*len(mutated)/len(ss) if available and ss else np.nan,'n_putative_truncating_splice_tumours':len(trunc) if available else np.nan,'putative_truncating_splice_percent':100*len(trunc)/len(ss) if available and ss else np.nan})
V=pd.DataFrame(variants);V.to_csv(R/'selected_somatic_variants.csv',index=False);M=pd.DataFrame(rows);M.to_csv(R/'mutation_prevalence.csv',index=False)

# Supporting BRCA prognostic analysis controlling for intrinsic subtype.
def bh(x):
 x=np.asarray(x,float);out=np.full(len(x),np.nan);ok=np.isfinite(x);v=x[ok];o=np.argsort(v)
 if len(v):out[np.where(ok)[0][o]]=np.minimum(1,np.minimum.accumulate((v[o]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1])
 return out
out=[]
for ep in ['PFI','DFI','OS']:
 for g in E:
  d=pd.DataFrame({'time':C[ep+'.time'],'event':C[ep],'gene_z':(E[g]-E[g].mean())/E[g].std(ddof=1),'age_per10':C.CDR_AGE/10,'stage_III_IV':C.STAGE_ADVANCED,'subtype':C.SUBTYPE})
  d=d.dropna();d=d.loc[d.time>0];d=d.join(pd.get_dummies(d.subtype,drop_first=True,dtype=float)).drop(columns='subtype');n=len(d);ev=int(d.event.sum());min_events=10*(len(d.columns)-2);s={'cohort':'BRCA_PAM50_assigned','endpoint':ep,'gene':g,'model':'age_stage_PAM50_adjusted','n':n,'events':ev,'model_covariates':';'.join(d.columns[2:]),'minimum_events_required':min_events,'HR_per_expression_SD':np.nan,'CI95_lower':np.nan,'CI95_upper':np.nan,'p_value':np.nan,'PH_rank_p_value':np.nan,'status':'insufficient_events' if ev<min_events else 'pending'}
  if ev>=min_events:
   try:
    with warnings.catch_warnings():
     warnings.simplefilter('ignore');f=CoxPHFitter().fit(d,'time','event');a=f.summary.loc['gene_z'];ph=proportional_hazard_test(f,d,time_transform='rank').summary.loc['gene_z','p']
    s.update(HR_per_expression_SD=a['exp(coef)'],CI95_lower=a['exp(coef) lower 95%'],CI95_upper=a['exp(coef) upper 95%'],p_value=a.p,PH_rank_p_value=ph,status='fitted')
   except Exception as e:s['status']='fit_failed: '+str(e)[:100]
  out.append(s)
S=pd.DataFrame(out);S['BH_FDR']=S.groupby('endpoint').p_value.transform(bh);S['PH_BH_FDR']=S.groupby('endpoint').PH_rank_p_value.transform(bh);S.to_csv(R/'survival_PAM50_sensitivity.csv',index=False)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
f,axs=plt.subplots(1,2,figsize=(11,8),sharey=True)
for a,co in zip(axs,['BRCA','TNBC']):
 d=M.loc[M.cohort.eq(co)&M.gene.isin(genes)].set_index('gene').reindex(genes);a.barh(range(len(genes)),d.mutated_percent,color='#8496A2',label='Any source-reported coding mutation');a.barh(range(len(genes)),d.putative_truncating_splice_percent,color='#A34F45',label='Putative truncating/splice subset');a.set_xlabel('Percentage of sequenced tumours');a.set_title(co+' (n='+str(int(d.n_sequenced.dropna().iloc[0]))+' sequenced)');a.set_xlim(0,3)
 for i,g in enumerate(genes):
  if pd.isna(d.loc[g,'n_sequenced']):a.text(1.5,i,'NA',color='.5',ha='center',va='center')
axs[0].set_yticks(range(len(genes)),genes);axs[0].invert_yaxis();axs[1].legend(loc='upper left',bbox_to_anchor=(1,1));f.suptitle('TCGA-BRCA: somatic mutation context',fontsize=15);f.text(.02,.01,'Coding mutations are generally rare in this panel. LINC03007 is unavailable in this coding WXS query.\nPutative truncating/splice calls do not establish biallelic functional knockout; germline variants and epigenetic loss are not assessed.',fontsize=9);f.tight_layout(rect=[0,.07,1,.95]);f.savefig(R/'figures/TCGA_F06_mutation_prevalence.png',dpi=300,bbox_inches='tight');plt.close(f)
print('MUTATION_COVERAGE',M.loc[M.cohort.eq('BRCA'),['gene','n_sequenced','n_mutated_tumours','n_putative_truncating_splice_tumours']].to_dict('records'),flush=True)
print('SUBTYPE_SENSITIVITY_PFI',S.loc[S.endpoint.eq('PFI')&S.BH_FDR.lt(.05)].to_dict('records'),flush=True)
