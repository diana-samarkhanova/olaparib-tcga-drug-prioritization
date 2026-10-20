from pathlib import Path
import json
import io
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'outputs/TCGA_Drug_Analysis_EN_20261013'
FIG = OUT / 'figures'
TAB = OUT / 'tables'
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.facecolor': 'white'})
KO = ['LAMTOR2','WDR24','DERL2','NME6','TP53BP1','PARG']
ACT = ['TP53I11','FOXQ1','CERK','PHF19','SLC35D1','PHLDA2','WT1',
       'HEPHL1','LINC03007','SLC47A1','ELL2','INHBA','PHGDH','ABCB1']
ROLES = {g:'Core activation' for g in ACT[:6]}
ROLES.update({'WT1':'Alternative','PHGDH':'Pharmacology branch','ABCB1':'Positive control'})
for g in ACT[7:12]: ROLES[g]='Exploratory'
for g in KO: ROLES[g]='KO candidate' if g != 'TP53BP1' else 'Canonical resistance control'

def save(fig, name):
    buffer=io.BytesIO()
    fig.savefig(buffer, format='png', dpi=300, bbox_inches='tight')
    temporary=FIG / (name + '.png.tmp')
    temporary.write_bytes(buffer.getvalue())
    temporary.replace(FIG / (name + '.png'))
    fig.savefig(FIG / (name + '.pdf'), bbox_inches='tight')
    plt.close(fig)

panel = []
for library, genes, times in [('KO',KO,[16,24]),('SAM',ACT,[14,23])]:
    d = pd.read_csv(ROOT/f'analysis/crispr/{library}_all_genes_two_timepoints.csv').set_index('gene')
    for gene in genes:
        source_gene=gene if gene in d.index else ('LOC646588' if gene=='LINC03007' else gene)
        if source_gene not in d.index:
            panel.append({'gene':gene,'library':library,'role':ROLES[gene],'source_gene':source_gene,'missing':True})
            continue
        row=d.loc[source_gene]
        for t in times:
            prefix=f'{library}_T{t}'
            panel.append({'gene':gene,'library':library,'role':ROLES[gene],
                          'source_gene':source_gene,'timepoint':f'T{t}',
                          'log2FC':row[f'{prefix}_pos_lfc'],
                          'positive_FDR':row[f'{prefix}_pos_fdr'],
                          'supporting_sgRNAs':row[f'{prefix}_pos_goodsgrna'],
                          'total_sgRNAs':row[f'{prefix}_num'],
                          'source_file':row[f'{prefix}_source_file'],
                          'source_excel_row':row[f'{prefix}_source_excel_row']})
panel=pd.DataFrame(panel)
panel.to_csv(TAB/'Selected_CRISPR_panel_20_genes.csv',index=False)
for library,genes,times,title,name in [
    ('KO',KO,['T16','T24'],'Knockout candidates under olaparib selection','F01_KO_screen_evidence'),
    ('SAM',ACT,['T14','T23'],'Activation candidates under olaparib selection','F02_Activation_screen_evidence')]:
    sub=panel[panel.library==library]
    vals=sub.pivot(index='gene',columns='timepoint',values='log2FC').reindex(index=genes,columns=times)
    fdr=sub.pivot(index='gene',columns='timepoint',values='positive_FDR').reindex(index=genes,columns=times)
    h=5.2 if library=='KO' else 8.4
    fig,axes=plt.subplots(1,2,figsize=(9.4,h),gridspec_kw={'width_ratios':[1.25,1]})
    vmax=max(1,np.nanmax(np.abs(vals.values)))
    im=axes[0].imshow(vals.values,cmap='RdBu_r',norm=TwoSlopeNorm(0,-vmax,vmax),aspect='auto')
    for i,g in enumerate(genes):
        for j,t in enumerate(times):
            v=vals.loc[g,t]; q=fdr.loc[g,t]
            mark='*' if q<.05 and v>0 else ''
            axes[0].text(j,i,f'{v:+.2f}{mark}',ha='center',va='center',color='white' if abs(v)>.7*vmax else 'black')
    axes[0].set(yticks=range(len(genes)),yticklabels=genes,xticks=range(2),xticklabels=times,title='Gene log2FC')
    axes[0].tick_params(length=0)
    fig.colorbar(im,ax=axes[0],fraction=.055,pad=.035,label='log2 fold change')
    axes[1].set_xlim(-.5,1.5); axes[1].set_ylim(len(genes)-.5,-.5)
    axes[1].set(xticks=[0,1],xticklabels=times,yticks=[],title='Positive-selection FDR')
    axes[1].tick_params(length=0)
    for i,g in enumerate(genes):
        for j,t in enumerate(times):
            q=fdr.loc[g,t]
            axes[1].text(j,i,f'{q:.3g}',ha='center',va='center',weight='bold' if q<.05 else 'normal')
    for ax in axes:
        if library=='SAM':
            for y in [5.5,6.5,11.5,12.5]: ax.axhline(y,color='#b8b8b8',lw=.8)
        for sp in ax.spines.values(): sp.set_visible(False)
    fig.suptitle(title,fontsize=14,y=.985)
    note='* Positive FDR < 0.05 and log2FC > 0 in the supplied comparison. FDR is from the genome-wide MAGeCK analysis.'
    if library=='SAM': note+='\nCore six first; WT1 alternative; five exploratory; PHGDH branch; ABCB1 control. Outlier-removed source tables.'
    else: note+='\nTP53BP1 retained as a canonical control; pooled enrichment is not individual KO validation.'
    fig.text(.03,.015,note,fontsize=8,va='bottom')
    fig.subplots_adjust(bottom=.14,top=.90,wspace=.5,left=.17,right=.97)
    save(fig,name)

coverage=json.loads((ROOT/'analysis/drugs_extended/expanded_coverage_summary.json').read_text())
source=coverage['source_coverage']
names=list(source)
direct=np.array([source[n]['direct_normalized_name_matches'] for n in names])
total=np.array([source[n]['catalog_compounds_with_actual_468_records'] for n in names])
fig,ax=plt.subplots(figsize=(9,4.8))
y=np.arange(len(names))
ax.barh(y,direct,color='#386c73',label='Direct normalized-name mapping')
ax.barh(y,total-direct,left=direct,color='#d2a25a',label='Provisional alias or chemical-form mapping')
for i,n in enumerate(total): ax.text(n+15,i,f'{n:,} / 2,321',va='center',fontsize=10)
ax.set(yticks=y,yticklabels=['PRISM primary','PRISM secondary','GDSC1','GDSC2'],xlabel='Catalog entries with parental MDA-MB-468 records',xlim=(0,2000))
ax.invert_yaxis(); ax.legend(loc='lower right',frameon=False,fontsize=9)
ax.set_title('Actual public drug-response coverage',fontsize=14)
fig.text(.02,.015,'All mappings are name or synonym based; chemical structures and physical library reagents were not validated.\nPRISM and GDSC assays differ in duration, dose windows and response metrics.',fontsize=8)
fig.subplots_adjust(bottom=.20,left=.18)
save(fig,'F03_Drug_resource_coverage')

profiles=pd.read_csv(ROOT/'analysis/drugs_extended/gene_linked_catalog_compound_profiles.csv')
codes=['B1464','A1402','A2614','A1933','A1805','A8167','A8169','A8314','A8218','B1867','B1399','B1970']
selected=profiles.set_index('CatalogNumber').reindex(codes)
selected.to_csv(TAB/'Focused_parental_468_drug_profiles.csv',index=True)
labels=['Selinexor (KPT-330)','Gemcitabine HCl','Bortezomib','Carfilzomib','Imatinib mesylate','Sirolimus','Everolimus','Temsirolimus','Lapatinib','Verapamil HCl','Temozolomide','Metformin HCl']
fig,ax=plt.subplots(figsize=(9.5,6.5))
v=selected.primary_468_log2FC.astype(float).values
ax.barh(np.arange(len(v)),v,color=np.where(v<0,'#386c73','#d2a25a'))
for i,x in enumerate(v): ax.text(x+(.06 if x>=0 else -.06),i,f'{x:.2f}',va='center',ha='left' if x>=0 else 'right',fontsize=9)
ax.axvline(0,color='#333333',lw=.8)
ax.set(yticks=np.arange(len(v)),yticklabels=labels,xlabel='PRISM primary log2 relative viability versus DMSO',xlim=(-5.4,1.6))
ax.invert_yaxis();ax.set_title('Baseline drug activity in parental MDA-MB-468',fontsize=14)
fig.text(.02,.015,'Selected compound-level median barcode readouts. Source dose usually 2.5 µM and exposure 5 days; actual dose is retained in the table.\nActivity in parental cells does not establish altered-genotype selectivity, combination benefit or olaparib resensitization.\nSalt and parent-form matches remain provisional where applicable.',fontsize=8)
fig.subplots_adjust(left=.25,bottom=.19,top=.92)
save(fig,'F04_Parental_468_PRISM_primary')

fig,ax=plt.subplots(figsize=(9.5,6.3))
auc=pd.to_numeric(selected.secondary_468_auc_median,errors='coerce').values
r2=pd.to_numeric(selected.secondary_468_r2_median,errors='coerce').values
for i,(a,r) in enumerate(zip(auc,r2)):
    if np.isfinite(a): ax.scatter(a,i,s=65,color='#386c73' if r>0 else '#d2a25a',marker='o' if r>0 else 'x',zorder=3)
    else: ax.text(.03,i,'Not covered',va='center',color='#888888',fontsize=9)
    if np.isfinite(a): ax.text(a+.025,i,f'R² {r:.2f}',va='center',fontsize=8)
ax.set(yticks=np.arange(len(labels)),yticklabels=labels,xlabel='Native PRISM secondary AUC (lower indicates stronger integrated inhibition)',xlim=(0,1.35),ylim=(len(labels)-.5,-.5))
ax.set_title('PRISM dose-response evidence and fit diagnostics',fontsize=14)
ax.scatter([],[],color='#386c73',label='Positive source R²');ax.scatter([],[],color='#d2a25a',marker='x',label='Nonpositive source R²')
ax.legend(loc='lower right',frameon=False,fontsize=9)
fig.text(.02,.015,'Native source AUCs are not clamped and are not pooled across datasets. Different dose windows affect AUC.\nPositive R² alone does not certify a strong curve fit. Lapatinib has very low R²; several rapalog fits are poor.\nMissing or extrapolated IC50 values are not replaced by zero.',fontsize=8)
fig.subplots_adjust(left=.25,bottom=.19,top=.92)
save(fig,'F05_Parental_468_PRISM_secondary')

gdsc=pd.read_csv(ROOT/'analysis/drugs_extended/gdsc2_catalog_MDA_MB_468.csv')
gdsc=gdsc[gdsc.CatalogNumber.isin(codes)].drop_duplicates(['DRUG_ID','NLME_CURVE_ID']).copy()
gdsc=gdsc.sort_values(['DRUG_NAME','DRUG_ID'])
gdsc.to_csv(TAB/'Focused_GDSC2_468_curve_records.csv',index=False)
if len(gdsc):
    fig,ax=plt.subplots(figsize=(10,max(4.8,.47*len(gdsc)+2)))
    for i,row in enumerate(gdsc.itertuples()):
        ax.plot([row.MIN_CONC,row.MAX_CONC],[i,i],color='#a9a9a9',lw=3,zorder=1)
        inrange=row.ic50_range_flag=='within_tested_range'
        # Source labels are audited below rather than assumed to use one spelling.
        inrange=np.isfinite(pd.to_numeric(row.ic50_in_tested_range_uM,errors='coerce'))
        ax.scatter(row.ic50_fitted_uM,i,color='#386c73' if inrange else '#b9762f',marker='o' if inrange else 'x',s=60,zorder=3)
    ax.set_xscale('log');ax.set(yticks=np.arange(len(gdsc)),yticklabels=[f'{r.DRUG_NAME}  (ID {r.DRUG_ID})' for r in gdsc.itertuples()],xlabel='Concentration (µM)',ylim=(len(gdsc)-.5,-.5))
    ax.set_title('GDSC2 fitted IC50 relative to the tested dose window',fontsize=14)
    ax.plot([],[],color='#a9a9a9',lw=3,label='Tested concentration interval');ax.scatter([],[],color='#386c73',label='Fit within tested interval');ax.scatter([],[],color='#b9762f',marker='x',label='Extrapolated fit')
    ax.legend(loc='upper center',bbox_to_anchor=(.5,-.25),ncol=3,frameon=False,fontsize=8)
    fig.text(.02,.015,'Parental MDA-MB-468, GDSC2 release 27 October 2023, 72-hour assay. exp(LN_IC50) is the source fitted estimate.\nAn estimate outside the source dose interval is extrapolated and cannot be described as a measured IC50.\nName/synonym mapping does not validate catalog salt or chemical structure.',fontsize=8)
    fig.subplots_adjust(left=.31,bottom=.39,top=.90)
    save(fig,'F06_GDSC2_IC50_dose_windows')
print('Created project evidence figures and exact plot tables in',OUT)
