"""Small data-selection helpers shared by the report exports."""

import pandas as pd


def format_cox_result(row):
    """Format the fitted result actually present in the input table."""
    if row.status != "fitted":
        return f"{row.gene}: model not fitted ({row.status})"
    return (f"{row.gene}: HR {row.HR_per_expression_SD:.3f} per expression SD "
            f"(95% CI {row.CI95_lower:.3f}–{row.CI95_upper:.3f}; "
            f"BH FDR {row.BH_FDR:.3g}; n={int(row.n)}, events={int(row.events)})")


def tcga_result_paragraphs(counts, cna, exp, cox, mut, sens, rppa):
    """Generate numerical TCGA prose from results instead of archived constants."""
    n = counts['rna_matched_primary_patients']
    status, subtypes = counts['TNBC_status'], counts['PAM50']
    cohort = (f"The eligible RNA cohort contains {n:,} patients: "
              f"{status.get('TNBC',0)} clinical TNBC, {status.get('non-TNBC',0)} non-TNBC "
              f"and {status.get('Unclassified',0)} unclassified. PAM50 assigns "
              f"{subtypes.get('BRCA_Basal',0)} basal-like tumours; "
              f"{subtypes.get('Unclassified',0)} have unassigned subtype. "
              f"Purity is available for {counts['purity_available']:,} patients. "
              "Clinical TNBC and PAM50 basal-like are analysed separately.")
    excluded = counts.get('redacted_RNA_profiled_patients_excluded',0)
    if excluded:
        cohort += (f" A conservative quality-control rule excludes {excluded} RNA-profiled "
                   "cases flagged for redaction in TCGA-CDR. These flags have heterogeneous "
                   "reasons, including sample-identity or eligibility concerns; the rule does "
                   "not imply that every flagged case is invalid or reflects consent withdrawal. "
                   "The unfiltered archived cohort is retained separately for audit comparison.")
    cohort += (" This is an observational primary-tumour cohort, not a cohort selected for "
               "absence of treatment or stratified by olaparib response.")
    endpoints = [r for r in counts['endpoint_counts'] if r['endpoint']=='PFI']
    pfi = "PFI availability: " + '; '.join(
        f"{r['cohort']} n={r['n']:,}, events={r['events']}" for r in endpoints) + '. '
    for group in ['TNBC','Basal']:
        sub=cox.loc[cox.cohort.eq(group)&cox.endpoint.eq('PFI')&cox.model.eq('univariable')]
        significant=sub.loc[sub.BH_FDR.lt(.05),'gene'].tolist()
        pfi += (f"{group}: " + (', '.join(significant) if significant else 'no candidate')
                + ' meets the candidate-panel PFI BH FDR <0.05 threshold. ')
    pfi += "Low subgroup event counts limit precision and preclude the prespecified age/stage-adjusted models."
    fits=cox.loc[cox.cohort.eq('BRCA')&cox.endpoint.eq('PFI')&cox.model.eq('age_stage_adjusted')].set_index('gene',drop=False)
    sensitivity=sens.loc[sens.endpoint.eq('PFI')].set_index('gene',drop=False)
    survival = 'Pooled-BRCA age/stage-adjusted estimates: ' + '; '.join(
        format_cox_result(fits.loc[g]) for g in ['HEPHL1','ABCB1']) + '. '
    survival += 'With PAM50 adjustment: ' + '; '.join(
        format_cox_result(sensitivity.loc[g]) for g in ['HEPHL1','ABCB1']) + '. '
    a=sensitivity.loc['ABCB1']
    survival += (f"ABCB1 rank-time proportional-hazards diagnostic: raw p={a.PH_rank_p_value:.3g}, "
                 f"PH-FDR={a.PH_BH_FDR:.3g}. These are prognostic associations in heterogeneous "
                 "treatment contexts. They do not measure or predict clinical olaparib response, "
                 "and do not negate a transporter mechanism established in cell experiments.")
    rp=rppa.loc[rppa.cohort.eq('TNBC')&rppa.BH_FDR.lt(.05)]
    protein=f"{len(rp)} exploratory RNA–RPPA associations meet the TNBC BH FDR threshold. "
    protein+='; '.join(f"{r.gene} RNA with {r.RPPA_readout.split('|')[-1]} "
                      f"(rho {r.spearman_rho:+.3f}, q={r.BH_FDR:.3g})" for r in rp.itertuples())+'. '
    protein+='These unadjusted correlations neither demonstrate candidate-induced pathway activation nor establish a drug-response mechanism.'
    denominator='Somatic-variant denominators: '+', '.join(
        f"{int(mut.loc[mut.cohort.eq(g),'n_sequenced'].dropna().iloc[0]):,} {g} tumours" for g in ['TNBC','BRCA'])+'. '
    denominator+='Shallow copy loss or a coding variant does not establish biallelic functional knockout.'
    ko=['LAMTOR2','WDR24','DERL2','NME6','TP53BP1','PARG']
    contrast=exp.loc[exp.comparison.eq('clinical_TNBC_vs_nonTNBC')].set_index('gene')
    nonsig=[g for g in ko if contrast.loc[g,'BH_FDR']>=.05]
    expression='For the KO panel, TNBC-versus-non-TNBC median transformed-expression differences are '+', '.join(
        f"{g} {contrast.loc[g,'median_delta_log2_RSEM_plus1']:+.3f} (q={contrast.loc[g,'BH_FDR']:.3g})" for g in ko)+'. '
    expression+=(', '.join(nonsig)+' do not meet q<0.05 in this contrast. ' if nonsig else '')
    expression+='Lower tumour expression is not evidence of complete gene loss or resistance to olaparib.'
    return dict(cohort=cohort,pfi=pfi,survival=survival,protein=protein,denominator=denominator,expression=expression)


def select_screen_timepoint(screen: pd.DataFrame, phase: str) -> pd.DataFrame:
    """Select the declared early/late comparison, independent of row order.

    Missing or duplicated gene comparisons are errors: silently dropping a gene
    from the integrated report would hide an incomplete input panel.
    """
    times = {"early": {"KO": "T16", "SAM": "T14"},
             "late": {"KO": "T24", "SAM": "T23"}}
    if phase not in times:
        raise ValueError(f"Unknown screen phase: {phase}")
    expected = screen["library"].map(times[phase])
    if expected.isna().any():
        raise ValueError("Unknown screen library in selected candidate panel")
    selected = screen.loc[screen["timepoint"].eq(expected)].copy()
    if selected["gene"].duplicated().any():
        raise ValueError(f"Duplicate {phase} comparison for a candidate gene")
    missing = set(screen["gene"]) - set(selected["gene"])
    if missing:
        raise ValueError(f"Missing {phase} comparison for: {', '.join(sorted(missing))}")
    return selected
