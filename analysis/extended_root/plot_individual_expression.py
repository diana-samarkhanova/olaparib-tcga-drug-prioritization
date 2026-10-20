from pathlib import Path
import io
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/TCGA_Drug_Analysis_EN_20261013/figures/genes'
OUT.mkdir(parents=True,exist_ok=True)
expr=pd.read_csv(ROOT/'analysis/tcga_extended/selected_expression_log2_RSEM_plus1.csv',index_col=0)
clinical=pd.read_csv(ROOT/'analysis/tcga_extended/matched_clinical_cohort.csv',index_col=0).reindex(expr.index)
contrast=pd.read_csv(ROOT/'analysis/tcga_extended/expression_contrasts.csv')
genes=['LAMTOR2','WDR24','DERL2','NME6','TP53BP1','PARG','TP53I11','FOXQ1','CERK','PHF19','SLC35D1','PHLDA2','WT1','HEPHL1','SLC47A1','ELL2','INHBA','PHGDH','ABCB1']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
rng=np.random.default_rng(4381)
for i,gene in enumerate(genes,1):
    fig,axes=plt.subplots(1,2,figsize=(11,5),sharey=True,gridspec_kw={'width_ratios':[1,1.6]})
    groups=[('TNBC_STATUS',['non-TNBC','TNBC'],['Receptor-positive','Clinical TNBC']),
            ('SUBTYPE',['BRCA_LumA','BRCA_LumB','BRCA_Her2','BRCA_Basal','BRCA_Normal'],['Luminal A','Luminal B','HER2-enriched','Basal-like','Normal-like'])]
    for ax,(column,categories,labels) in zip(axes,groups):
        values=[expr.loc[clinical[column].eq(c),gene].dropna().values for c in categories]
        boxes=ax.boxplot(values,positions=range(len(values)),widths=.48,showfliers=False,patch_artist=True,
                         medianprops={'color':'#111111','linewidth':1.5},boxprops={'linewidth':1},
                         whiskerprops={'linewidth':1},capprops={'linewidth':1})
        for j,(box,v) in enumerate(zip(boxes['boxes'],values)):
            box.set_facecolor('#d2a25a' if categories[j] in ['TNBC','BRCA_Basal'] else '#c8d7d9')
            ax.scatter(j+rng.uniform(-.18,.18,len(v)),v,s=7,alpha=.19,color='#465556',rasterized=True,zorder=1)
        ax.set_xticks(range(len(values)),[f'{label}\nn={len(v)}' for label,v in zip(labels,values)],rotation=0,fontsize=9)
        ax.set_xlim(-.55,len(values)-.45)
    d=contrast.loc[(contrast.gene==gene)&(contrast.comparison=='clinical_TNBC_vs_nonTNBC')].iloc[0]
    axes[0].set_title(f"TNBC median difference {d.median_delta_log2_RSEM_plus1:+.2f}\nBH FDR {d.BH_FDR:.3g}",fontsize=11)
    axes[1].set_title('PAM50 tumour subtypes',fontsize=11)
    axes[0].set_ylabel('log2(RSEM + 1)')
    fig.suptitle(f'{gene} expression in TCGA-BRCA primary tumours',fontsize=14)
    fig.text(.02,.018,'One primary tumour per patient. Boxes show median and IQR; whiskers extend to 1.5 × IQR; points show all observations.\nClinical contrast excludes unknown receptor status. Normal-like denotes a tumour subtype, not normal tissue.\nBH correction covers the 19 measurable panel genes; the median difference is not a log2 fold change.',fontsize=8)
    fig.subplots_adjust(left=.085,right=.985,bottom=.25,top=.78,wspace=.22)
    b=io.BytesIO();fig.savefig(b,format='png',dpi=300,bbox_inches='tight')
    p=OUT/f'G{i:02d}_{gene}_TCGA_expression.png';p.write_bytes(b.getvalue());plt.close(fig)
print('Created',len(genes),'individual expression PNG figures.')
