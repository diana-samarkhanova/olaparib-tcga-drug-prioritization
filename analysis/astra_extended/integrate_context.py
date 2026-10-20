"""Join actual source screen and TCGA context to mechanism matrix; preserve missing values."""
from pathlib import Path
import pandas as pd
P=Path(__file__).resolve().parent;root=P.parent.parent
m=pd.read_csv(P/'final_gene_priority_matrix.csv')
m=m[[c for c in m if not c.startswith(('screen_','tcga_','gse_'))]]
s=pd.read_csv(root/'tables/Selected_CRISPR_panel_20_genes.csv')
s['strict_positive']=(s.positive_FDR<.05)&(s.log2FC>0)
a=[]
for g,x in s.groupby('gene'):
 a.append(dict(gene=g,screen_timepoints=';'.join(x.timepoint),screen_strict_positive_timepoints=';'.join(x.loc[x.strict_positive,'timepoint']),screen_all_displayed_timepoints_strict_positive=bool(x.strict_positive.all()),screen_direction_caveat='Positive-test FDR alone is insufficient; T14 gene log2FC is negative' if g=='SLC47A1' else ''))
m=m.merge(pd.DataFrame(a),on='gene',how='left')
t=P.parent/'tcga_extended'
if (t/'expression_contrasts.csv').exists():
 x=pd.read_csv(t/'expression_contrasts.csv');x=x[x.comparison.eq('clinical_TNBC_vs_nonTNBC')][['gene','median_delta_log2_RSEM_plus1','BH_FDR','n_group1','n_group2']];x=x.rename(columns={c:'tcga_TNBC_vs_nonTNBC_'+c for c in x if c!='gene'});m=m.merge(x,on='gene',how='left')
if (t/'CNA_prevalence_and_expression.csv').exists():
 x=pd.read_csv(t/'CNA_prevalence_and_expression.csv');x=x[x.cohort.eq('TNBC')][['gene','n_profiled','deep_deletion_n','amplification_n','dosage_expression_spearman_rho','BH_FDR']];x=x.rename(columns={c:'tcga_TNBC_CNA_'+c for c in x if c!='gene'});m=m.merge(x,on='gene',how='left')
if (P/'gse165914_candidate_context.csv').exists():
 x=pd.read_csv(P/'gse165914_candidate_context.csv')[['gene','contextual_direction_vs_screen','resistant_15uM_vs_parental_untreated_author_log2FC','resistant_15uM_vs_parental_untreated_author_FDR']];x=x.rename(columns={c:'gse_'+c for c in x if c!='gene'});m=m.merge(x,on='gene',how='left')

if (t/'pathway_associations.csv').exists():
 x=pd.read_csv(t/'pathway_associations.csv');x=x[x.cohort.eq('TNBC')]
 programmes={'LAMTOR2':'MTORC1_SIGNALING','WDR24':'MTORC1_SIGNALING','DERL2':'UNFOLDED_PROTEIN_RESPONSE','NME6':'OXIDATIVE_PHOSPHORYLATION','TP53BP1':'DNA_REPAIR','PARG':'DNA_REPAIR','TP53I11':'HYPOXIA','FOXQ1':'EPITHELIAL_MESENCHYMAL_TRANSITION','CERK':'MTORC1_SIGNALING','PHF19':'G2M_CHECKPOINT','SLC35D1':'UNFOLDED_PROTEIN_RESPONSE','PHLDA2':'MTORC1_SIGNALING','WT1':'EPITHELIAL_MESENCHYMAL_TRANSITION','PHGDH':'GLYCOLYSIS','INHBA':'TGF_BETA_SIGNALING','ELL2':'DNA_REPAIR'}
 rows=[]
 for g,program in programmes.items():
  z=x[x.gene.eq(g)&x.pathway.eq(program)]
  if len(z):
   r=z.iloc[0];rows.append({'gene':g,'tcga_TNBC_context_program':program,'tcga_TNBC_context_partial_rank_rho':r.partial_rank_rho,'tcga_TNBC_context_adjusted_BH_FDR':r.adjusted_BH_FDR,'tcga_TNBC_context_n':r.adjusted_n})
 m=m.merge(pd.DataFrame(rows),on='gene',how='left')

m.to_csv(P/'final_gene_priority_matrix.csv',index=False)
print(m.shape)
