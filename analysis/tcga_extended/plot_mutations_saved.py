from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent
(R/'figures').mkdir(exist_ok=True)
M=pd.read_csv(R/"mutation_prevalence.csv")
genes=["LAMTOR2","WDR24","DERL2","NME6","TP53BP1","PARG","TP53I11","FOXQ1","CERK","PHF19","SLC35D1","PHLDA2","WT1","HEPHL1","LINC03007","SLC47A1","ELL2","INHBA","PHGDH","ABCB1"]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
f,axs=plt.subplots(1,2,figsize=(11,8),sharey=True)
for a,co in zip(axs,['BRCA','TNBC']):
 d=M.loc[M.cohort.eq(co)&M.gene.isin(genes)].set_index('gene').reindex(genes);a.barh(range(len(genes)),d.mutated_percent,color='#8496A2',label='Any source-reported coding mutation');a.barh(range(len(genes)),d.putative_truncating_splice_percent,color='#A34F45',label='Putative truncating/splice subset');a.set_xlabel('Percentage of sequenced tumours');a.set_title(co+' (n='+str(int(d.n_sequenced.dropna().iloc[0]))+' sequenced)');a.set_xlim(0,3)
 for i,g in enumerate(genes):
  if pd.isna(d.loc[g,'n_sequenced']):a.text(1.5,i,'NA',color='.5',ha='center',va='center')
axs[0].set_yticks(range(len(genes)),genes);axs[0].invert_yaxis();axs[1].legend(loc='upper left',bbox_to_anchor=(1,1));f.suptitle('TCGA-BRCA: somatic mutation context',fontsize=15);f.text(.02,.01,'Coding mutations are generally rare in this panel. LINC03007 is unavailable in this coding WXS query.\nPutative truncating/splice calls do not establish biallelic functional knockout; germline variants and epigenetic loss are not assessed.',fontsize=9);f.tight_layout(rect=[0,.07,1,.95]);f.savefig(R/'figures/TCGA_F06_mutation_prevalence.png',dpi=300,bbox_inches='tight');plt.close(f)
