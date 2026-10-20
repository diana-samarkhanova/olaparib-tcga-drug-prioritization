from pathlib import Path
import io,textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/TCGA_Drug_Analysis_EN_20261013'
(OUT/'figures').mkdir(parents=True,exist_ok=True)
(OUT/'tables').mkdir(parents=True,exist_ok=True)
d=pd.read_csv(ROOT/'analysis/crispr/enrichment/KO_reactome_ora_significant.csv').sort_values('bh_fdr',ascending=False)
d.to_csv(OUT/'tables/Reactome_KO_recurrent_hits_significant.csv',index=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(figsize=(12,6.7))
y=np.arange(len(d));v=-np.log10(d.bh_fdr.to_numpy(float));colors=['#d2a25a' if 'mTORC1' in x else '#386c73' for x in d.pathway_name]
ax.barh(y,v,color=colors,height=.65)
for i,(x,row) in enumerate(zip(v,d.itertuples())):
 ax.text(x+.025,i,f'q={row.bh_fdr:.3g}; {row.foreground_pathway_overlap} genes',va='center',fontsize=9)
ax.set_yticks(y,['\n'.join(textwrap.wrap(x,47)) for x in d.pathway_name],fontsize=9)
ax.invert_yaxis();ax.set_xlabel('−log10 BH FDR');ax.set_xlim(0,2.7)
ax.axvline(-np.log10(.05),ls='--',color='#666666',lw=.8)
ax.set_title('Reactome enrichment of recurrent positive-selection KO hits',fontsize=14,pad=15)
fig.text(.02,.012,'Genome-wide foreground: 24 recurrent positive hits, 20 mapped; assayed annotated background: 10,831 genes.\nOne-sided hypergeometric tests; BH correction across all 1,713 eligible pathways. Reactome release 97.\nThe mTORC1 pathway contains LAMTOR2, MLST8 and WDR24. Other overlapping annotations share genes and do not imply endothelial cell identity.\nThis is enrichment of the full recurrent KO hit set, not an over-representation test of the six selected candidates.',fontsize=8)
fig.subplots_adjust(left=.40,right=.98,bottom=.20,top=.89)
b=io.BytesIO();fig.savefig(b,format='png',dpi=300,bbox_inches='tight');(OUT/'figures/F08_Reactome_KO_enrichment.png').write_bytes(b.getvalue())
fig.savefig(OUT/'figures/F08_Reactome_KO_enrichment.pdf',bbox_inches='tight');plt.close(fig)
print('Created Reactome enrichment figure')
