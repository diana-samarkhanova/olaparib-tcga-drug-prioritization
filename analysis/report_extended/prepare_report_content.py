from pathlib import Path
import json,re
import pandas as pd

R=Path(__file__).resolve().parents[2]; T=R/'analysis/tcga_extended'; D=R/'analysis/drugs_extended'; A=R/'analysis/astra_extended'
O=R/'outputs/TCGA_Drug_Analysis_EN_20261013'; blocks=[]
def p(text):blocks.append({'type':'paragraph','text':text})
def h(text,level=1):blocks.append({'type':'heading','text':text,'level':level})
def tab(headers,rows,widths,fs=9.5):blocks.append({'type':'table','headers':headers,'rows':rows,'widths':widths,'font_size':fs})
def fig(path,caption,width=6.65):blocks.append({'type':'figure','path':str(Path(path).relative_to(R)),'caption':caption,'width':width})
def md(text):
    paragraphs=re.split(r'\n\s*\n',text.strip())
    for x in paragraphs:
        x=x.strip()
        if x.startswith('#'):
            lines=x.splitlines(); line=lines.pop(0);h(line.lstrip('# ').strip(),2 if line.startswith('###') else 1)
            if lines:p(' '.join(lines))
        elif x.startswith('|'):
            rows=[[c.strip() for c in line.strip('|').split('|')] for line in x.splitlines() if line.startswith('|')]
            if len(rows)>2:
                cols=len(rows[0]);w=[7/cols]*cols;tab(rows[0],rows[2:],w,9)
        else:p(' '.join(x.splitlines()).replace('`','').replace("the user's",'the project').replace("user's",'project'))

counts=json.loads((T/'cohort_counts.json').read_text())
cov=pd.read_csv(T/'gene_coverage.csv');assert cov.rna_coverage.eq('available').sum()==19
cna=pd.read_csv(T/'CNA_prevalence_and_expression.csv'); assert 'LAMTOR2' in cna.gene.values
exp=pd.read_csv(T/'expression_contrasts.csv');paths=pd.read_csv(T/'pathway_associations.csv');cox=pd.read_csv(T/'survival_continuous_Cox.csv')
mut=pd.read_csv(T/'mutation_prevalence.csv');sens=pd.read_csv(T/'survival_PAM50_sensitivity.csv');rppa=pd.read_csv(T/'RPPA_associations.csv')
drugcov=json.loads((D/'expanded_coverage_summary.json').read_text())
scr=pd.read_csv(O/'tables/Selected_CRISPR_panel_20_genes.csv')
ko=['LAMTOR2','WDR24','DERL2','NME6','TP53BP1','PARG'];act=['TP53I11','FOXQ1','CERK','PHF19','SLC35D1','PHLDA2','WT1','HEPHL1','LINC03007','SLC47A1','ELL2','INHBA','PHGDH','ABCB1']

p('The most defensible pharmacological branches are ABCB1 transport blockade as a reversal benchmark, FOXQ1–PDGFR inhibition as a conditional adjacent-target hypothesis, DERL2–gemcitabine and PARG–temozolomide as collateral-sensitivity hypotheses, and sirolimus or selinexor as empirical olaparib-combination comparators with published MDA-MB-468 evidence. These branches address different outcomes and should not be collapsed into one sensitizer ranking.')
p('The analysis combines the six KO candidates and fourteen activation candidates with newly executed TCGA-BRCA profiling, published acquired-resistance RNA data, all 2,321 named L1021 catalog entries, PRISM, GDSC, primary mechanism studies and published olaparib-combination experiments. No public result demonstrates restored olaparib sensitivity in the individual engineered derivatives. The value of the analysis is to identify which hypotheses have a causal anchor, which have tumour context, and which still require mechanism before pharmacology.')
h('Scope and evidence definitions')
tab(['Panel','Genes and roles'],[['Knockout','LAMTOR2, WDR24, DERL2, NME6, TP53BP1 and PARG'],['Core activation','TP53I11, FOXQ1, CERK, PHF19, SLC35D1 and PHLDA2'],['Additional activation','WT1 alternative; PHGDH pharmacology branch; HEPHL1, LINC03007, SLC47A1, ELL2 and INHBA exploratory'],['Positive control','ABCB1 activation, kept distinct from discovery candidates']],[1.5,5.5])
p('Three endpoints are retained throughout. Selective monotherapy means preferential activity of the second drug in the altered model. Combination benefit means a response beyond a declared single-agent reference. Olaparib resensitization means recovery of an olaparib response beyond partner-alone toxicity. A strong parental-cell IC50, tumour co-expression, or a favorable published CI from a different model is insufficient to establish the third endpoint.')

h('Resources and analytical methods')
tab(['Resource','Actual analysis','Coverage'],[
['Project CRISPR','Supplied CNV-corrected MAGeCK RRA gene summaries; positive-selection FDR and guide support retained','40 candidate–timepoint records'],
['TCGA-BRCA','RNA, GISTIC CNA, public somatic mutations, RPPA, purity and clinical endpoints','1,082 RNA-matched primary patients; 19/20 genes'],
['GSE165914','Extraction of published genome-wide differential-expression statistics and processed sample values','19/20 candidates; parental, acute and resistant-exposure contexts'],
['PRISM','Original primary fixed-dose and secondary dose-response resources','1,502 and 531 catalog entries linked to actual parental 468 records'],
['GDSC1 and GDSC2','Official release 8.5 fitted dose-response workbooks','76 and 89 catalog entries linked to parental 468 records'],
['Combination studies','Exact-line and external-line published interaction or response evidence','18 source-specific evidence rows; nine in-library author hits from the 166-drug screen'],
['Regulatory and primary literature','FDA indication-aware audit and model-specific mechanism evidence','Current, withdrawn, investigational and unresolved statuses separated'],
],[1.4,3.75,1.85],9.2)
p('TCGA data were obtained from the public PanCancer Atlas breast-cancer study, supplemented by original Firehose receptor annotations. One primary tumour per patient was retained. The frozen batch-normalized RSEM matrix was transformed with log2(RSEM+1); it was not mixed with TPM or another expression pipeline. LAMTOR2 maps to historical ROBLD3 in RNA and current LAMTOR2 in CNA, both Entrez 28956. LINC03007/LOC646588 is unavailable in this matrix. Clinical TNBC requires ER-negative, PR-negative and resolved HER2-negative status; discordant HER2 results cannot qualify as HER2-negative. Clinical TNBC and PAM50 basal-like cohorts are analysed separately. [TCGA public study](https://github.com/cBioPortal/datahub/tree/master/public/brca_tcga_pan_can_atlas_2018) [Receptor supplement](https://github.com/cBioPortal/datahub/tree/master/public/brca_tcga) [GDC PanCancer resources](https://gdc.cancer.gov/about-data/publications/pancanatlas)')
p('Expression contrasts use two-sided Mann–Whitney tests, median differences and rank-biserial effect sizes. BH correction is applied across the 19 available candidates within each contrast. GISTIC −2/−1/0/+1/+2 defines deep deletion, shallow loss, diploid, gain and amplification. Dosage–expression associations use Spearman correlation with cohort-specific BH correction. Somatic coding variants and putative truncating/splice variants are reported descriptively against sequenced denominators; variant annotation alone does not establish biallelic functional loss.')
p('Ten prespecified Hallmark programmes are scored as mean gene-level z scores, using MSigDB v2025.1.Hs. They cover mTORC1, UPR, OXPHOS, glycolysis, DNA repair, G2M checkpoint, TGF-beta, EMT, apoptosis and hypoxia. The candidate is excluded from its own score. Partial rank associations adjust age, E2F proliferation score, purity where covered and target CNA; pooled BRCA models also include PAM50. BH correction covers all candidate–pathway tests within a cohort. These are co-expression analyses, not intervention or synthetic-lethality experiments. RPPA associations are exploratory unadjusted Spearman tests, corrected over candidate–readout pairs; phosphoprotein levels are not normalized to total protein. [Hallmark collection](https://www.gsea-msigdb.org/gsea/msigdb/human/collection_details.jsp#H)')
p('TCGA-CDR PFI is the primary prognostic endpoint, DFI is supporting and OS exploratory. Cox models use continuous expression scaled per cohort SD. Univariable models require at least 15 events and age/stage-adjusted models at least 30. Pooled-BRCA sensitivity models additionally include PAM50 when event counts permit. Model-wise BH correction and proportional-hazards diagnostics are retained. No optimized expression cut point is selected. PFI is not an olaparib-specific outcome. [Clinical resource and endpoint guidance](https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018)')
p('PRISM primary is a barcode-derived log2 response relative to DMSO at the retained source dose; 2^LFC is a derived relative readout. Secondary AUC, R-squared, concentration windows and source batch remain separate, with the documented MTS010 technical-redo preference. GDSC uses exp(LN_IC50) in micromolar and retains each drug ID, native AUC, dose range and range flag. Within-range values are fitted estimates from a tested curve; extrapolated values are not measured potency. No AUC is pooled across sources or dose windows. All identity joins are name/synonym based, with zero canonical structures or physical reagent forms validated. [PRISM](https://depmap.org/repurposing/) [GDSC downloads](https://www.cancerrxgene.org/downloads/bulk_download)')

h('Project screen evidence')
p('LAMTOR2, WDR24, DERL2 and NME6 satisfy positive FDR <0.05 at both KO timepoints with four supporting guides at each. NME6 T24 FDR is 0.048809, close to the threshold. PARG is significant at T24 only; TP53BP1 is retained as a repair-context control despite not passing the positive FDR threshold in either supplied KO comparison. Repeated timepoints are supporting screen evidence, not independent validation of isolated clones.')
fig(O/'figures/F01_KO_screen_evidence.png','Figure 1. Six selected KO candidates in the supplied MAGeCK RRA tables. LFC describes guide enrichment under olaparib, not gene expression. Statistical thresholds, guide support and timepoint-specific uncertainty remain separate.')
p('The activation panel preserves the predefined biological prioritization rather than restricting selection to hits recurrent at both timepoints. Early and late candidates can have different trajectories. In particular, SLC35D1 and LINC03007 have negative early LFC and positive late selection, while INHBA changes from positive early to negative late LFC. ABCB1 is the positive-control branch. These patterns require timepoint-aware individual validation; external RNA concordance cannot erase a direction switch.')
fig(O/'figures/F02_Activation_screen_evidence.png','Figure 2. Fourteen activation candidates with supplied positive-selection FDR and supporting-guide counts. Early-only, late-only, recurrent and direction-changing patterns are shown without creating a combined efficacy score.')
p('The previously executed full-screen Reactome v97 enrichment provides an independent annotation layer. Its foreground is the 24 recurrent positive KO hits, of which 20 map; the background is 10,831 annotated genes from 19,112 assayed genes. One-sided hypergeometric tests cover 1,713 eligible pathways, with six BH-significant results. Amino acids regulate mTORC1 includes LAMTOR2, MLST8 and WDR24 (q=0.0303); G1/S transcription and integration of energy metabolism are also enriched. Shared endothelial/shear-stress annotations reflect overlapping signalling genes rather than establishing an endothelial identity for 468. This enrichment is not an ORA of the six selected genes and does not determine a drug’s direction. The full table and separate figure are supplied. [Reactome](https://reactome.org/)')

h('TCGA cohort and expression context')
p('The primary RNA cohort contains 1,082 patients: 157 clinical TNBC, 848 non-TNBC and 77 unclassified. PAM50 assigns 171 basal-like tumours; 101 have unassigned subtype. Purity is available for 1,021 patients. The classification rules and missing categories remain visible in the cohort table; basal-like and TNBC are overlapping biological definitions rather than interchangeable sample labels.')
ce=exp.query("comparison=='clinical_TNBC_vs_nonTNBC'").set_index('gene')
vals=[]
for g in ko:
    row=ce.loc[g];vals.append(f'{g} {row.median_delta_log2_RSEM_plus1:+.3f} (q={row.BH_FDR:.3g})')
p('For the KO panel, TNBC-versus-non-TNBC median transformed-expression differences are '+', '.join(vals)+'. LAMTOR2 and DERL2 do not meet q<0.05 in this contrast. Lower tumour expression of a KO candidate gives loss-related context but is not evidence of complete gene deletion or of resistance to olaparib.')
p('Among core activation candidates, FOXQ1, CERK, PHF19 and PHLDA2 are higher in clinical TNBC, while TP53I11 and SLC35D1 are lower, all at candidate-panel BH FDR <0.05. PHGDH has the largest positive median difference in the panel. A low TNBC average does not exclude an activation-induced phenotype in 468; tumour subtype enrichment and intervention direction answer different questions.')
fig(T/'figures/TCGA_F01_expression_contrasts.png','Figure 3. Newly calculated TCGA candidate expression contrasts under separate clinical and intrinsic-subtype definitions. Δ is a difference in median log2(RSEM+1), not a log2 fold change. Asterisks use candidate-panel BH FDR <0.05; LINC03007 is unavailable.')

h('Copy number and somatic loss context')
rows=[]
for g in ko:
    x=cna.query("cohort=='TNBC' and gene==@g").iloc[0]
    y=mut.query("cohort=='TNBC' and gene==@g").iloc[0]
    rows.append([g,int(x.n_profiled),int(x.deep_deletion_n),int(x.shallow_loss_n),int(y.n_mutated_tumours),int(y.n_putative_truncating_splice_tumours)])
tab(['KO gene','CNA n','Deep loss','Shallow loss','Mutated n','Truncating or splice n'],rows,[1.05,.7,1,1.05,1,2.2],9.3)
p('Deep deletions and coding loss-type variants are uncommon in this TNBC subset. Somatic-variant denominators are 157 sequenced TNBC tumours and 1,064 pooled-BRCA tumours. LAMTOR2 instead shows frequent gains/amplifications, despite a small nonsignificant lower median RNA contrast. This is a direct warning against interpreting a tumour expression comparison as an engineered loss model. Shallow loss cannot establish a biallelic KO. The public somatic calls also do not provide a complete assessment of germline status, methylation, repair function or allele-specific inactivation.')
fig(T/'figures/TCGA_F02_CNA_prevalence.png','Figure 4. GISTIC copy-number distributions in pooled BRCA and clinical TNBC. Percentages use gene-specific profiled denominators. Copy-number categories provide dosage context, not proof of CRISPR-equivalent KO or activation.')

h('Within TNBC pathway and protein associations')
ap=paths.query("cohort=='TNBC'")
focus=[('LAMTOR2','MTORC1_SIGNALING'),('WDR24','MTORC1_SIGNALING'),('DERL2','UNFOLDED_PROTEIN_RESPONSE'),('NME6','OXIDATIVE_PHOSPHORYLATION'),('FOXQ1','EPITHELIAL_MESENCHYMAL_TRANSITION'),('CERK','MTORC1_SIGNALING'),('INHBA','TGF_BETA_SIGNALING'),('PHLDA2','GLYCOLYSIS')]
rows=[]
for g,path in focus:
    x=ap.query('gene==@g and pathway==@path').iloc[0];rows.append([g,path.replace('_',' ').title(),int(x.adjusted_n),f'{x.partial_rank_rho:+.3f}',f'{x.adjusted_BH_FDR:.3g}'])
tab(['Gene','Programme','Adjusted n','Partial rho','BH FDR'],rows,[1,3.15,.85,.85,1.15],9.3)
p('The patient data support contextual ER-stress and mitochondrial programmes for DERL2/NME6, and strong activin/TGF-beta context for INHBA. They do not support a simple positive mTORC1 RNA relationship for LAMTOR2/WDR24, and FOXQ1–EMT is not significant after the predefined adjustment. CERK associates negatively with the mTORC1 transcript programme, opposite to the signalling direction in the overexpression experiments. Transcript programme scores and acute phosphosignalling are different measurements; this contradiction should constrain extrapolation rather than erase the perturbation evidence.')
fig(T/'figures/TCGA_F03_TNBC_pathway_associations.png','Figure 5. Adjusted candidate–programme co-expression within clinical TNBC. Candidate genes are excluded from their corresponding scores; age, proliferation, purity and target CNA are included where covered. Stars use BH correction across the full within-TNBC test family. Associations do not establish compensatory dependence.')
rp=rppa.query("cohort=='TNBC' and BH_FDR<0.05")
p('Only two exploratory RNA–RPPA associations survive the TNBC correction: TP53BP1 RNA with 53BP1 protein (rho +0.368, q=0.00149) and PHF19 RNA with AKT pS473 (rho −0.391, q=0.000802). The full matrix is supplied as a separate PNG. These unadjusted correlations neither demonstrate pathway activation by a candidate nor identify a drug-response mechanism.')

h('Prognosis and its limits')
p('PFI is available for 1,067 pooled-BRCA patients with 145 events, 157 TNBC patients with 28 events and 169 basal-like patients with 25 events. The subgroup event counts permit predefined univariable models but not the age/stage-adjusted three-covariate models. No TNBC or basal-like candidate meets PFI BH FDR <0.05. Absence of significance in these small event sets is not proof of no clinical relevance.')
p('In pooled BRCA, age/stage-adjusted HEPHL1 associates with higher PFI hazard (HR 1.236 per SD, 95% CI 1.071–1.427, q=0.0359) and ABCB1 with lower hazard (HR 0.706, 95% CI 0.604–0.825, q=0.000226). In the PAM50-adjusted sensitivity model (944 patients, 120 events), the HEPHL1 signal loses support; ABCB1 retains an inverse association (HR 0.675, 95% CI 0.563–0.810, q=0.000420). Its rank-time proportional-hazards diagnostic has raw p=0.0369 and PH-FDR=0.702, so a constant HR warrants caution. These are prognostic observations in heterogeneous patient treatment contexts. The ABCB1 result does not negate its experimentally established olaparib efflux mechanism and does not support calling high ABCB1 a clinical olaparib-response biomarker.')
fig(T/'figures/TCGA_F04_PFI_forest.png','Figure 6. Continuous-expression PFI models, using age/stage adjustment in pooled BRCA and univariable subgroup models where adjustment event requirements are unmet. Full sample counts, event counts, BH FDR, proportional-hazards diagnostics and subtype sensitivity models are retained in the tables.')

h('Acquired resistance RNA context')
text=(A/'scientific_findings_en.md').read_text()
gsetext=text.split('## Additional acquired-resistance transcript context: GSE165914',1)[1].split('## Patient context',1)[0]
md(gsetext)
fig(A/'GSE165914_candidate_context.png','Figure 7. Targeted extraction of the published GSE165914 differential-expression statistics and processed sample context. Resistant cultures are under continuous 15 µM olaparib, confounding resistance and ongoing exposure. Values preserve author genome-wide FDR rather than a newly fitted candidate-only DE model.')

md((R/'analysis/report_extended/mechanisms_compact.md').read_text())

h('Drug priorities and public parental activity')
p('The entire L1021 catalog was mined before model-specific selection. Actual parental 468 numerical coverage is 1,502 PRISM primary entries (933 direct normalized-name links), 531 PRISM secondary entries (315 direct links), 76 GDSC1 entries (19 direct links) and 89 GDSC2 entries (28 direct links). Aliases, salts and development codes account for the remaining links. Different catalog forms can map to one active moiety; these are not chemically verified unique-compound counts. CTRP archival access was unsuccessful and supplies no numerical result in this report.')
tab(['Branch','Primary candidate','Evidence role'],[['ABCB1 activation','Verapamil B1867; lapatinib comparator','External causal olaparib transport reversal'],['FOXQ1 activation','Imatinib A1805/B2171','Conditional PDGFR blockade; other chemotherapy endpoints'],['DERL2 knockout','Gemcitabine A1402','Gene-loss collateral sensitivity in another lineage'],['PARG knockout','Temozolomide B1399','Alternative-drug sensitivity in a BRCA2-deficient setting'],['CERK activation','Sirolimus A8167','Downstream migration rescue and separate empirical 468 combination'],['Any model as comparator','Selinexor KPT-330 B1464','Published parental 468 olaparib-combination evidence']],[1.55,2.4,3.05],9.3)
p('Gemcitabine has strong parental 468 activity in multiple resources. Its PRISM secondary fitted IC50 is 0.00184 µM, while separate within-range GDSC fits are 0.0156/0.0221 µM in GDSC1 and 0.165 µM in GDSC2. These estimates are not pooled. Strong baseline killing supports feasibility but can obscure genotype selectivity. Conversely, the temozolomide GDSC fits of approximately 218/498 µM exceed the 30 µM test maximum; they are extrapolations and cannot be presented as measured potency or used to dismiss a KO-created sensitivity.')
p('Imatinib parental 468 primary relative readout is approximately 0.477 at 2.5 µM, whereas secondary AUC is 0.887 with R-squared 0.174 and no retained in-range IC50. Sirolimus has within-range fitted estimates of 0.0516 µM in PRISM and 0.0186 µM in GDSC2, the latter from a very narrow tested interval. These source-specific results determine assay feasibility, not which engineered model will respond.')
fig(O/'figures/F05_Parental_468_PRISM_secondary.png','Figure 8. Source-native PRISM secondary AUC and fit diagnostics for selected catalog compounds in parental MDA-MB-468. Dose windows, quality flags and unavailable in-range IC50 remain explicit. AUC is not combined with CRISPR enrichment or TCGA expression.')
fig(O/'figures/F06_GDSC2_IC50_dose_windows.png','Figure 9. Parental MDA-MB-468 GDSC2 fitted IC50 estimates relative to each tested concentration interval. Out-of-range values are extrapolations. Separate drug IDs and assay contexts remain visible; the plot does not estimate genotype-specific sensitivity.')

md((R/'analysis/combinations_extended/combo_results_en.md').read_text())

h('Regulatory status and chemical interpretation')
p('Automated Drugs@FDA ingredient and product matches are an audit layer, not complete adjudication of every catalog form or current indication. Tazemetostat/Tazverik full-NDA approval was withdrawn on 22 June 2026, and panobinostat/Farydak on 24 March 2022. They belong in the historical mechanistic-probe arm. Elacridar and zosuquidar remain investigational controls. Selinexor retains current multiple-myeloma indications despite withdrawal of its DLBCL indication. Missing automated FDA matches remain unresolved rather than being classified as nonapproved. [FDA withdrawn indications](https://www.fda.gov/drugs/resources-information-approved-drugs/withdrawn-cancer-accelerated-approvals) [Tazverik full-NDA notice](https://www.govinfo.gov/content/pkg/FR-2026-06-22/pdf/2026-12367.pdf) [Drugs@FDA data](https://www.fda.gov/drugs/drug-approvals-and-databases/drugsfda-data-files)')
p('Free base, salt, solvate, prodrug and active form are not automatically interchangeable. MLN9708/ixazomib citrate and MLN2238/active ixazomib remain distinct. Fludarabine phosphate is kept separate from fludarabine. Current drug identity is checked against actual catalog IDs and retrieved NIH synonym records; no docking-based affinity or predicted cellular IC50 is introduced.')

h('Recommended order of experimental follow up')
p('First, establish the resistance phenotype and untreated division/death response in each individual derivative. ABCB1 provides a useful pharmacological reversal benchmark. Use FOXQ1–imatinib only when PDGFR activation is demonstrated. For DERL2 and PARG, examine partner-alone genotype selectivity before adding olaparib. Sirolimus and selinexor are the most useful exact-468 empirical combination comparators, with the caveat that sirolimus could have a different direction in nutrient-signalling KOs. LAMTOR2, WDR24 and NME6 remain high-interest mechanistic models without a validated approved-drug rescue from the available resources.')
p('Second, align the endpoint with the hypothesis. A minimally toxic transporter-blocking condition is informative for ABCB1; metabolic or lysosomal branches require basal-growth and pathway measurements; repair branches require functional HR/replication-stress context. BRCA sequence status does not establish HR competence. Genotype specificity requires comparison of the interaction or differential response across edited and control cells, rather than significance in one model and nonsignificance in another.')
p('Third, retain contradictory data. Vorinostat antagonism in parental 468, TMZ antagonism in the parental 231 screen, CERK’s negative patient mTORC1 association, PHGDH counterdirection and activation candidates with opposite resistant-RNA expression all limit universal sensitizer claims. Untreated growth and metabolic assay effects are especially important for NME6, PHGDH and nutrient-signalling KOs. Growth-rate metrics require baseline and endpoint counts and cannot be reconstructed from the public summaries used here. [Hafner et al.](https://lincs.hms.harvard.edu/hafner-natmethods-2016/)')
p('The resulting priorities are evidence roles rather than a computed efficacy score. TCGA adds tumour context and exposes confounding or opposing directions. Primary perturbation and combination studies supply the causal anchors. Actual public drug-response curves supply baseline pharmacology. Together they support a focused follow-up screen without claiming that any compound has already reversed resistance in the engineered models.')

h('References and resource links')
# Normalize source-note shorthand in the reader-facing report.
replacements={
    'parental468':'parental MDA-MB-468','resistant468':'resistant MDA-MB-468',
    'parental231':'parental MDA-MB-231','positive231':'positive MDA-MB-231',
    'universal468':'universal MDA-MB-468','validated468':'validated MDA-MB-468',
    'not468':'not MDA-MB-468','current2026 label':'current 2026 label',
    'the project individual KO':'the individual KO',
}
for b in blocks:
    if 'text' in b:
        for old,new in replacements.items():b['text']=b['text'].replace(old,new)
        b['text']=re.sub(r'(?<=\d)(µM|µm|hours|hours\b)',r' \1',b['text'])
        b['text']=re.sub(r'PMID(?=\d)','PMID ',b['text'])
alltext='\n'.join(b.get('text','') for b in blocks)
pat=re.compile(r'\[([^\]]+)\]\((https?://(?:[^\s()]|\([^)]*\))+)\)')
links=[]
for label,url in pat.findall(alltext):
    if url not in [x[1] for x in links]:links.append((label,url))
sources=pd.read_csv(A/'sources.csv').fillna('')
sourceby={row.url:row for _,row in sources.iterrows()}
for i,(label,url) in enumerate(links,1):
    if url in sourceby:
        x=sourceby[url]; year=str(x.year).removesuffix('.0');doi=f' DOI {x.doi}.' if x.doi else ''
        text=f'{i}. {x.title}. {year}.{doi} [Source]({url})'
    else:text=f'{i}. {label}. [Source]({url})'
    blocks.append({'type':'reference','text':text})
content={'actual_results_ready':True,'title':'TCGA and drug resource analysis of olaparib resistance candidates','subtitle':'Six knockout candidates and fourteen activation candidates in MDA MB 468\n13 October 2026','blocks':blocks}
(R/'analysis/report_extended/report_content.json').write_text(json.dumps(content,indent=2))
print(json.dumps({'blocks':len(blocks),'references':len(links),'words':sum(len(b.get('text','').split()) for b in blocks)}))
