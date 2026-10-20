"""Re-extract published edgeR results, audit sample-level processed counts, plot context.
No de novo DE fitting: supplied counts are fractional processed values and are not
asserted to be raw integer read counts. Author FDR preserved; no candidate-only FDR.
"""
from pathlib import Path
import pandas as pd, numpy as np, hashlib, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent
matrix=pd.read_csv(P/'final_gene_priority_matrix.csv')
de=pd.read_csv(P/'raw/GSE165914_M468P_OR_Diff_Expression_mshekhar.csv.gz')
counts=pd.read_csv(P/'raw/GSE165914_SUM_M468_Parental_OR_counts_mshekhar.csv.gz',index_col=0)
rows=[];long=[]
for _,gene in matrix.iterrows():
 g=gene.gene; r={'gene':g,'perturbation':gene.perturbation,'present_in_processed_table':g in counts.index}
 for label,col in [('acute_parental_3uM_vs_untreated',0),('resistant_15uM_vs_parental_untreated',20)]:
  a=de.loc[de.iloc[:,col].eq(g)]
  for field,offset in [('author_log2FC',1),('author_PValue',2),('author_FDR',3)]:r[label+'_'+field]=float(a.iloc[0,col+offset]) if len(a) else np.nan
  r[label+'_author_FDR_lt_0_05']=bool(len(a) and a.iloc[0,col+3]<.05)
 for label,prefix in [('parental_untreated','M468PCtrl'),('parental_acute_3uM','M468P3Ola'),('resistant_continuous_15uM','M468OR15Ola')]:
  cols=[prefix+str(i) for i in (1,2,3)]
  vals=counts.loc[g,cols] if g in counts.index else pd.Series([np.nan]*3)
  r[label+'_mean_supplied_processed_count']=vals.mean()
  for i,v in enumerate(vals):long.append({'gene':g,'condition':label,'sample':cols[i],'supplied_processed_count':v})
 fc=r['resistant_15uM_vs_parental_untreated_author_log2FC'];sig=r['resistant_15uM_vs_parental_untreated_author_FDR_lt_0_05']
 r['contextual_direction_vs_screen']=('not_available' if pd.isna(fc) else 'not_FDR_significant' if not sig else 'same_direction' if ((fc<0)==(gene.perturbation=='KO')) else 'opposite_direction')
 rows.append(r)
x=pd.DataFrame(rows);x.to_csv(P/'gse165914_candidate_context.csv',index=False);pd.DataFrame(long).to_csv(P/'gse165914_candidate_sample_values.csv',index=False)
fig,ax=plt.subplots(figsize=(10,10));y=np.arange(len(x))
for label,off,color,name in [('acute_parental_3uM_vs_untreated',-.15,'#84909c','Parental: acute 3 µM olaparib'),('resistant_15uM_vs_parental_untreated',.15,'#146b83','Acquired resistant: continuous 15 µM')]:
 vals=x[label+'_author_log2FC']; sig=x[label+'_author_FDR_lt_0_05']
 ax.scatter(vals,y+off,c=[color if s else 'white' for s in sig],edgecolors=color,s=55,label=name,zorder=3)
ax.axvline(0,c='#555',lw=.8);ax.set_yticks(y,x.gene);ax.invert_yaxis();ax.set_xlabel('Published log2 fold change versus untreated parental MDA-MB-468')
ax.set_title('GSE165914: context for the 20 screen candidates',loc='left',fontsize=15,pad=18);ax.grid(axis='x',alpha=.18);ax.spines[['top','right']].set_visible(False);ax.legend(loc='upper right',fontsize=9)
fig.text(.13,.022,'Filled = author genome-wide FDR < 0.05; open = not significant; LINC03007 unavailable.\nResistance and exposure are confounded. These are author DE estimates, not isogenic gene perturbations.',fontsize=9)
fig.subplots_adjust(left=.16,bottom=.11,top=.92);fig.savefig(P/'GSE165914_candidate_context.png',dpi=220);plt.close(fig)
manifest={'source':'GSE165914','retrieved_utc':'2026-10-13','analysis':'New extraction of published edgeR glmQLFTest/BH results; descriptive means only from supplied processed counts','count_matrix_shape':list(counts.shape),'n_candidate_genes':20,'n_mapped':int(x.present_in_processed_table.sum()),'condition_n':3,'no_new_DE_pvalues':True,'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (P/'raw').glob('*.gz')}}
(P/'gse165914_provenance.json').write_text(json.dumps(manifest,indent=2))
print(x[['gene','contextual_direction_vs_screen']].to_string(index=False))
