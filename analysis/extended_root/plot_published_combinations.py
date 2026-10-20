from pathlib import Path
import io
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/TCGA_Drug_Analysis_EN_20261013'
(OUT/'figures').mkdir(parents=True,exist_ok=True)
(OUT/'tables').mkdir(parents=True,exist_ok=True)
data=pd.read_csv(ROOT/'analysis/olaparib_screen/L1021_screen_overlap.csv')
hits=pd.read_csv(ROOT/'analysis/combinations_extended/Yusoh_nine_library_hits.csv')
codes=hits.catalog_number.tolist()+['B1399','A8167','B2171','A8217','B1464']
# Vorinostat catalog number is resolved from the source overlap rather than assumed.
v=data[data.library_name.str.contains('Vorinostat',case=False)].catalog_number.iloc[0]
codes=[v if c=='A8217' else c for c in codes]
sub=data.set_index('catalog_number').reindex(codes)
if sub.ci_0_1uM.isna().any(): raise ValueError('Requested compound lacks published numeric record')
sub.to_csv(OUT/'tables/Published_combination_figure_values.csv',index=True)
labels=['Chlorambucil','Tamoxifen citrate','Fludarabine phosphate','Exemestane','Abiraterone',
        'Omacetaxine / homoharringtonine','Panobinostat [withdrawn US]','Plerixafor','Acalabrutinib',
        'Temozolomide','Sirolimus','Imatinib','Vorinostat','Selinexor']
cols=['ci_0_1uM','ci_1uM','ci_10uM']
a=sub[cols].to_numpy(float)
fig,ax=plt.subplots(figsize=(10.8,8.5))
im=ax.imshow(np.log10(a),cmap='BrBG_r',norm=TwoSlopeNorm(0,-1.3,1.3),aspect='auto')
for i in range(len(a)):
    for j in range(3):
        ax.text(j,i,f'{a[i,j]:.2f}',ha='center',va='center',color='white' if abs(np.log10(a[i,j]))>.9 else 'black',fontsize=10)
ax.set(yticks=range(len(labels)),yticklabels=labels,xticks=range(3),xticklabels=['0.1 µM partner','1 µM partner','10 µM partner'])
ax.tick_params(length=0);ax.axhline(8.5,color='black',lw=1.2)
for sp in ax.spines.values(): sp.set_visible(False)
fig.colorbar(im,ax=ax,label='log10 of published combination index',fraction=.045,pad=.045)
fig.suptitle('Published olaparib combinations in parental MDA-MB-231',fontsize=14,y=.98)
fig.text(.02,.012,'Yusoh et al. 2025 supplementary Table S2; fixed olaparib 10 µM, 24-hour assay. Values shown are published CI, not newly calculated synergy.\nNine author-selected library hits above the dividing line; relevant comparators below. CI is dose and model dependent.\nRapamycin has positive evidence in another 72-hour resistant-468 study despite these 231 results; vorinostat has negative 468 evidence.\nThese measurements do not establish responses of the engineered KO or activation models.',fontsize=8)
fig.subplots_adjust(left=.38,right=.90,bottom=.17,top=.91)
b=io.BytesIO();fig.savefig(b,format='png',dpi=300,bbox_inches='tight');(OUT/'figures/F07_Published_olaparib_combinations.png').write_bytes(b.getvalue())
fig.savefig(OUT/'figures/F07_Published_olaparib_combinations.pdf',bbox_inches='tight');plt.close(fig)
print('Created published combination PNG and exact numeric table')
