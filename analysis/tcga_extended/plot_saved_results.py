from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
ROOT=Path(__file__).resolve().parent;FIG=ROOT/"figures";FIG.mkdir(exist_ok=True)
KO=["LAMTOR2","WDR24","DERL2","NME6","TP53BP1","PARG"]
ACT=["TP53I11","FOXQ1","CERK","PHF19","SLC35D1","PHLDA2","WT1","HEPHL1","LINC03007","SLC47A1","ELL2","INHBA","PHGDH","ABCB1"]
GENES=KO+ACT
PATHS=["MTORC1_SIGNALING","UNFOLDED_PROTEIN_RESPONSE","OXIDATIVE_PHOSPHORYLATION","GLYCOLYSIS","DNA_REPAIR","G2M_CHECKPOINT","TGF_BETA_SIGNALING","EPITHELIAL_MESENCHYMAL_TRANSITION","APOPTOSIS","HYPOXIA"]
contrasts=pd.read_csv(ROOT/"expression_contrasts.csv");cnarows=pd.read_csv(ROOT/"CNA_prevalence_and_expression.csv");assoc=pd.read_csv(ROOT/"pathway_associations.csv");surv=pd.read_csv(ROOT/"survival_continuous_Cox.csv");rprows=pd.read_csv(ROOT/"RPPA_associations.csv");expr=pd.read_csv(ROOT/"selected_expression_log2_RSEM_plus1.csv",index_col=0);cnamat=pd.read_csv(ROOT/"selected_CNA_GISTIC.csv",index_col=0)
rpnames=["RPS6KB1|p70S6K_pT389","RPS6|S6_pS235_S236","RPS6|S6_pS240_S244","EIF4EBP1|4E-BP1_pT37T46","MAPK1 MAPK3|MAPK_pT202_Y204","AKT1 AKT2 AKT3|Akt_pS473","TP53BP1|53BP1","RAD51|Rad51","XBP1|XBP1"]
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,"axes.spines.right":False,"savefig.facecolor":"white"})
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
