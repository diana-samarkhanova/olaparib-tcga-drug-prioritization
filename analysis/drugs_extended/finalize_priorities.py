"""Produce auditable model-specific nominations from reviewed evidence and source records.
No TCGA result, synergy prediction or edited-cell drug response is manufactured.
"""
from pathlib import Path
import pandas as pd,json
P=Path(__file__).resolve().parent
profiles=pd.read_csv(P/'gene_linked_catalog_compound_profiles.csv').set_index('CatalogNumber')
ev=pd.read_csv(P/'gene_mechanism_evidence_used.csv')
fdasource='https://www.fda.gov/drugs/drug-approvals-and-databases/drugsfda-data-files'
withdrawals='https://www.fda.gov/drugs/resources-information-approved-drugs/withdrawn-cancer-accelerated-approvals'
# Gene, library perturbation, role, compound ID, action, biological gate, evidence title term.
rows=[
('LAMTOR2','KO','Conditional pathway probe','A3018','Trametinib only if ERK compensation is observed; sirolimus and hydroxychloroquine are directional comparators','Measure mTORC1, ERK and basal proliferation; further inhibition can reproduce resistance-associated pathway loss','LAMTOR2 regulates'),
('WDR24','KO','Conditional pathway probe','A8167','Sirolimus is a nutrient-signalling comparator, not a nominated reversal drug','Establish actual mTORC1/lysosomal response and compensatory dependency after KO','GATOR2 Component'),
('DERL2','KO','Named collateral-drug hypothesis','A1402','Prioritize gemcitabine alone first; evaluate olaparib combination separately','External KO evidence is in QBC939 cholangiocarcinoma; establish selective effect relative to parental468','stabilizes BAG6'),
('NME6','KO','No specific approved sensitizer assigned','','Investigate pyrimidine salvage, mtRNA/OXPHOS and glycolytic compensation before assigning a drug','Further OXPHOS inhibition need not reverse an OXPHOS-defective state; metformin is a contextual comparator only','Ribonucleotide synthesis'),
('TP53BP1','KO','Research bypass outside approved shortlist','','ATR inhibitors provide a mechanistic comparator; no selective FDA-library ATR agent nominated','BRCA1-deficient HR restoration literature is not automatically applicable to BRCA-proficient468','ATR inhibition'),
('PARG','KO','Named collateral-drug hypothesis','B1399','Prioritize temozolomide as an alternative cytotoxicity branch','Published evidence is primarily BRCA2-deficient; TMZ is neither PARG replacement nor established olaparib resensitizer','Selective Loss of PARG'),
('TP53I11','CRISPRa','No specific approved sensitizer assigned','','Prioritize growth/EMT/ER-calcium mechanism validation; doxorubicin is an alternative-response comparator','Activation may slow proliferation rather than induce a generic pro-survival pathway','suppress'),
('FOXQ1','CRISPRa','Named downstream-drug hypothesis','A1805','Prioritize imatinib conditionally after verifying PDGFRA/PDGFRB pathway induction','Breast-model evidence concerns doxorubicin/paclitaxel chemoresistance; no engineered468 olaparib proof','PDGFRalpha'),
('CERK','CRISPRa','Named downstream-drug hypothesis','A8167','Prioritize sirolimus as the strongest reviewed approved downstream probe','Published rescue endpoint is migration; establish whether mTOR dependence also explains olaparib response','Ceramide kinase mediates'),
('PHF19','CRISPRa','Historical mechanistic probe only','A8221','EPZ-6438/tazemetostat may probe EZH2 biology but is excluded from current-approved shortlist','Catalytic EZH2 inhibition differs from PHF19/PRC2 clustering; fork-protection direction can oppose desired sensitization','PHF19 drives'),
('SLC35D1','CRISPRa','No specific approved sensitizer assigned','','Retain discovery branch; quantify nucleotide-sugar/glycosylation consequences','Shared substrates do not establish a UXS1/UGDH drug dependency','Nucleotide-sugar transporter'),
('PHLDA2','CRISPRa','No unconditional approved candidate assigned','','PI3K/AKT/mTOR drugs remain conditional after defining direction','Breast EGFR/ErbB2 negative-feedback evidence conflicts with positive AKT direction in other tumour lineages','negative feedback'),
('WT1','CRISPRa','Exploratory adjacent-pathway comparator','A1805','Imatinib only as an exploratory comparator if PDGFR changes are demonstrated','Proteomic association after WT1 silencing is weaker than causal FOXQ1-PDGFR-drug evidence','Proteomics analysis'),
('HEPHL1','CRISPRa','Exploratory functional probe','B8692','Penicillamine is a copper-homeostasis comparator, not a demonstrated HEPHL1 inhibitor','Establish iron/copper biology and phenotype before attributing drug response','Biallelic HEPHL1'),
('LINC03007','CRISPRa','No specific approved sensitizer assigned','','Confirm locus identity and neighbouring promoter effects; no defensible drug nominated','LOC646588 is the former identifier; transcript attribution and functional mechanism remain unresolved','NCBI Gene'),
('SLC47A1','CRISPRa','Transport-function comparator','A4353','Pyrimethamine can probe MATE1 function; no ABCB1-style olaparib reversal claim','Olaparib has inhibitor interaction with MATE1, not a validated MATE1-mediated efflux mechanism in this analysis','Discovery of potent'),
('ELL2','CRISPRa','No specific approved sensitizer assigned','','Investigate transcription elongation and NHEJ; SEC inhibitors are research comparators','No specific approved in-library ELL2/NHEJ bypass was established','ELL2 regulates'),
('INHBA','CRISPRa','Research mechanism first','','Validate activin-SMAD and transporter induction; ALK4/5/7 research inhibition is outside approved arm','No direct approved candidate; activin can reduce ixazomib sensitivity in another breast context','Activin-A signaling'),
('PHGDH','CRISPRa','Nonselective biochemical repurposing probe','A4015','Disulfiram is a conditional repurposing comparator; selective PHGDH inhibition is the mechanistic reference','Multiple drug actions and copper conditions confound attribution; PHGDH-low ovarian PARPi resistance provides counterdirection','disulfiram inhibits'),
('ABCB1','CRISPRa positive control','Transport-bypass validation priority','B1867','Prioritize verapamil plus olaparib; lapatinib is an orthogonal approved-active-moiety comparator, elacridar an investigational control','Confirm functional efflux and reversal at minimally cytotoxic inhibitor exposure','ABCB1 (MDR1) induction'),
]
out=[]
for gene,pert,role,cat,action,gate,term in rows:
 e=ev[ev.gene.eq(gene)]
 chosen=e[e.primary_title.str.contains(term,case=False,na=False,regex=False)]
 if chosen.empty: chosen=e
 r={'gene':gene,'perturbation':pert,'nomination_category':role,'primary_catalog_number':cat,'action':action,'interpretation_gate':gate,'primary_evidence_title':chosen.primary_title.iloc[0] if len(chosen) else '', 'primary_evidence_url':chosen.url.iloc[0] if len(chosen) else '', 'all_reviewed_gene_urls':'; '.join(sorted(set(e.url.dropna()))), 'edited_MDA_MB_468_resensitization_demonstrated':False,'TCGA_integration':'See separately computed TCGA result tables; drug nomination is not a TCGA response prediction'}
 if cat:
  x=profiles.loc[cat]
  for k in ['library_name','library_CAS','regulatory_review_category','regulatory_evidence_note','chemical_structure_validated','primary_468_source_names','primary_468_identity_mapping_methods','primary_468_relative_viability','primary_468_dose_min_uM','primary_468_dose_max_uM','secondary_468_source_names','secondary_468_identity_mapping_methods','secondary_468_auc_median','secondary_468_r2_median','secondary_468_ic50_uM_median_within_tested_range','secondary_468_dose_min_uM_across_selected_records','secondary_468_dose_max_uM_across_selected_records','GDSC1_468_curve_details_JSON','GDSC2_468_curve_details_JSON']:r[k]=x[k]
  r['FDA_regulatory_source_url']=withdrawals if cat=='A8221' else fdasource
 else:
  r.update({'library_name':'No nominated library compound','regulatory_review_category':'Not applicable','FDA_regulatory_source_url':''})
 out.append(r)
out=pd.DataFrame(out)
assert len(out)==20 and out.gene.is_unique
out.to_csv(P/'final_model_priorities.csv',index=False)
# Root may plot each endpoint in a separate labelled panel; native AUCs stay separate.
focus=['A1402','A1805','A8167','B1399','B1867','A8218','A3384','A4015','A8221','A2614','A1933','A3018','B4874','A4353','B1464']
profiles.loc[focus].reset_index().to_csv(P/'final_focus_compound_profiles.csv',index=False)
long=[]
for cat in focus:
 x=profiles.loc[cat]
 for ds in ['GDSC1','GDSC2']:
  for c in json.loads(x[ds+'_468_curve_details_JSON']):
   long.append({'CatalogNumber':cat,'library_name':x.library_name,'source':ds,**c})
pd.DataFrame(long).to_csv(P/'final_focus_GDSC_curves.csv',index=False)
print('Written 20-gene model table, 15-compound profile table and source-specific GDSC curve table')
