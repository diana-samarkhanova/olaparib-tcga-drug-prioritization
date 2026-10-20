"""End-to-end tests of analyze_tcga.run on the synthetic snapshot.

They check the output-file contract relied on by downstream scripts, the
alias/coverage logic, consistency of the summary counts, and recovery of
planted signals. They are synthetic integration checks, not validation of
archived patient-level results or clinical drug-response predictions.
"""
import json

import numpy as np
import pandas as pd
import pytest

from analysis.tcga_extended import analyze_tcga as m
import synthetic_tcga

OUTPUT_FILES = ["selected_expression_log2_RSEM_plus1.csv", "hallmark_mean_z_scores.csv",
                "hallmark_coverage.csv", "matched_clinical_cohort.csv", "cohort_counts.json",
                "expression_contrasts.csv", "selected_CNA_GISTIC.csv",
                "CNA_prevalence_and_expression.csv", "gene_coverage.csv",
                "pathway_associations.csv", "RPPA_associations.csv",
                "survival_continuous_Cox.csv", "analysis_summary.json", "RNA_symbol_resolution.csv"]

# Columns read by analyze_mutations_and_sensitivity.py, integrate_context.py,
# plot_saved_results.py and report_extended/*.
DOWNSTREAM_COLUMNS = {
    "matched_clinical_cohort.csv": ["SAMPLE_ID", "PATIENT_ID", "TNBC_STATUS", "BASAL_STATUS", "SUBTYPE",
                                    "CDR_AGE", "STAGE_ADVANCED", "PURITY", "PFI", "PFI.time", "DFI",
                                    "DFI.time", "OS", "OS.time", "HER2_FINAL"],
    "expression_contrasts.csv": ["comparison", "gene", "panel", "n_group1", "n_group2",
                                 "median_delta_log2_RSEM_plus1", "p_value", "BH_FDR"],
    "CNA_prevalence_and_expression.csv": ["cohort", "gene", "n_profiled", "deep_deletion_n",
                                          "shallow_loss_n", "diploid_n", "gain_n", "amplification_n",
                                          "dosage_expression_spearman_rho", "BH_FDR"],
    "pathway_associations.csv": ["cohort", "gene", "pathway", "partial_rank_rho", "adjusted_BH_FDR",
                                 "adjusted_n", "raw_spearman_rho", "raw_BH_FDR", "covariates"],
    "RPPA_associations.csv": ["cohort", "gene", "RPPA_readout", "n", "spearman_rho", "BH_FDR"],
    "survival_continuous_Cox.csv": ["cohort", "endpoint", "gene", "model", "n", "events",
                                    "HR_per_expression_SD", "CI95_lower", "CI95_upper", "p_value",
                                    "PH_rank_p_value", "status", "BH_FDR", "PH_BH_FDR"],
    "gene_coverage.csv": ["gene", "source_symbol", "rna_rows", "entrez_gene_id", "rna_coverage",
                          "CNA_source_symbol", "CNA_coverage"],
}


def test_all_contract_files_and_columns_are_written(pipeline_outputs):
    root, _, _ = pipeline_outputs
    for name in OUTPUT_FILES:
        assert (root / name).exists(), name
    for name, columns in DOWNSTREAM_COLUMNS.items():
        header = pd.read_csv(root / name, nrows=0).columns
        missing = [c for c in columns if c not in header]
        assert not missing, f"{name} lacks {missing}"


def test_gene_coverage_uses_aliases_and_reports_absent_gene(pipeline_outputs):
    root, tables, _ = pipeline_outputs
    cov = pd.read_csv(root / "gene_coverage.csv").set_index("gene")
    assert list(cov.index) == m.GENES
    assert cov.loc["LAMTOR2", "source_symbol"] == "ROBLD3"
    assert cov.loc["LAMTOR2", "CNA_source_symbol"] == "ROBLD3"
    assert cov.loc["LINC03007", "rna_coverage"] == "unavailable_or_ambiguous"
    assert cov.loc["LINC03007", "CNA_coverage"] == "unavailable"
    assert cov.rna_coverage.eq("available").sum() == 19  # asserted by prepare_report_content.py
    assert list(tables["expr"].columns) == [g for g in m.GENES if g != "LINC03007"]


def test_expression_table_is_log2_rsem_plus_one(pipeline_outputs, synthetic_raw):
    root, _, _ = pipeline_outputs
    raw, _ = synthetic_raw
    expr = pd.read_csv(root / "selected_expression_log2_RSEM_plus1.csv", index_col="SAMPLE_ID")
    rna = pd.read_csv(raw / "data_mrna_seq_v2_rsem.txt", sep="\t").set_index("Hugo_Symbol")
    sample = expr.index[5]
    assert expr.loc[sample, "ABCB1"] == pytest.approx(np.log2(rna.loc["ABCB1", sample] + 1))
    assert expr.loc[sample, "LAMTOR2"] == pytest.approx(np.log2(rna.loc["ROBLD3", sample] + 1))


def test_clinical_classifications_are_internally_consistent(pipeline_outputs):
    root, _, _ = pipeline_outputs
    c = pd.read_csv(root / "matched_clinical_cohort.csv")
    assert c.SAMPLE_ID.is_unique and c.PATIENT_ID.is_unique
    assert c.SAMPLE_ID.str.len().eq(15).all()
    recomputed = [m.classify_tnbc(e, p, h) for e, p, h in
                  zip(c.ER_STATUS_BY_IHC, c.PR_STATUS_BY_IHC, c.HER2_FINAL, strict=True)]
    assert recomputed == c.TNBC_STATUS.tolist()
    assert c.loc[c.TNBC_STATUS.eq("TNBC"), ["ER_STATUS_BY_IHC", "PR_STATUS_BY_IHC", "HER2_FINAL"]].eq("Negative").all().all()
    assert c.loc[c.BASAL_STATUS.eq("Basal"), "SUBTYPE"].eq("BRCA_Basal").all()
    assert c.PURITY.notna().sum() > 0 and c.PURITY.notna().sum() < len(c)
    counts = json.loads((root / "cohort_counts.json").read_text())
    assert counts["rna_matched_primary_patients"] == len(c)
    assert counts["clinical_primary_patients"] == len(c)  # counted after one-primary-per-patient dedup
    assert not c.SAMPLE_ID.str.endswith(("-06", "-01B")).any()  # metastatic and duplicate rows dropped
    assert counts["TNBC_status"] == c.TNBC_STATUS.value_counts().to_dict()
    assert counts["purity_available"] == int(c.PURITY.notna().sum())
    cdr_other = counts["PAM50"]
    assert sum(cdr_other.values()) == len(c)


def test_contrasts_recover_planted_tnbc_signal(pipeline_outputs):
    root, _, _ = pipeline_outputs
    x = pd.read_csv(root / "expression_contrasts.csv")
    assert len(x) == 2 * 19
    assert set(x.comparison) == {"clinical_TNBC_vs_nonTNBC", "PAM50_Basal_vs_other"}
    row = x.set_index(["comparison", "gene"]).loc[("clinical_TNBC_vs_nonTNBC", "PHGDH")]
    assert row.median_delta_log2_RSEM_plus1 > 0.8 and row.BH_FDR < 0.05 and row.rank_biserial > 0.3
    assert row.panel == "CRISPRa"
    assert x.loc[x.gene.eq("LAMTOR2"), "panel"].eq("KO").all()
    assert (x.n_group1 + x.n_group2 <= 320).all()
    assert x.groupby("comparison").BH_FDR.apply(lambda q: (q >= 0).all() and (q <= 1).all()).all()


def test_cna_prevalence_counts_sum_to_profiled(pipeline_outputs):
    root, _, _ = pipeline_outputs
    cna = pd.read_csv(root / "CNA_prevalence_and_expression.csv")
    assert len(cna) == 3 * 19 and "LAMTOR2" in cna.gene.values
    bands = cna[["deep_deletion_n", "shallow_loss_n", "diploid_n", "gain_n", "amplification_n"]].sum(axis=1)
    assert (bands == cna.n_profiled).all()
    assert (cna.dosage_expression_n <= cna.n_profiled).all()
    sizes = cna.groupby("cohort").n_profiled.max()
    assert sizes["BRCA"] > sizes["TNBC"] and sizes["BRCA"] > sizes["Basal"]
    gistic = pd.read_csv(root / "selected_CNA_GISTIC.csv", index_col="SAMPLE_ID")
    assert set(gistic.stack().unique()) <= {-2.0, -1.0, 0.0, 1.0, 2.0}


def test_pathway_associations_structure_and_covariates(pipeline_outputs):
    root, _, truth = pipeline_outputs
    a = pd.read_csv(root / "pathway_associations.csv")
    assert len(a) == 3 * 19 * 10
    assert a.target_excluded_from_score.all()
    n_dna = len(truth["sets"]["DNA_REPAIR"])
    assert a.set_index(["cohort", "gene", "pathway"]).loc[("BRCA", "TP53BP1", "DNA_REPAIR"), "signature_gene_n"] == n_dna - 1
    assert a.set_index(["cohort", "gene", "pathway"]).loc[("BRCA", "ABCB1", "DNA_REPAIR"), "signature_gene_n"] == n_dna
    brca_cov = a.loc[a.cohort.eq("BRCA"), "covariates"].iloc[0].split(";")
    tnbc_cov = a.loc[a.cohort.eq("TNBC"), "covariates"].iloc[0].split(";")
    assert brca_cov[:4] == ["age", "E2F", "purity", "target_CNA"] and any(c.startswith("BRCA_") for c in brca_cov)
    assert tnbc_cov == ["age", "E2F", "purity", "target_CNA"]
    assert (a.adjusted_n <= a.raw_n).all()
    assert a.partial_rank_rho.abs().max() <= 1


def test_rppa_rows_use_available_samples_only(pipeline_outputs):
    root, _, truth = pipeline_outputs
    r = pd.read_csv(root / "RPPA_associations.csv")
    assert len(r) == 3 * 19 * 9 and set(r.RPPA_readout) == set(m.RPPA_READOUTS)
    assert r.loc[r.cohort.eq("BRCA"), "n"].max() == len(truth["rppa_samples"])
    assert (r.n >= m.MIN_RPPA_N).all()


def test_survival_models_recover_planted_wt1_effect(pipeline_outputs):
    root, _, _ = pipeline_outputs
    s = pd.read_csv(root / "survival_continuous_Cox.csv")
    assert len(s) == 3 * 3 * 19 * 2
    assert set(s.status) <= {"fitted", "insufficient_events"}
    assert list(s.columns[-3:]) == ["BH_FDR", "PH_BH_FDR", "fit_warnings"]
    minimum = dict((name, k) for name, _, k in m.COX_MODELS)
    for _, row in s.iterrows():
        assert (row.status == "insufficient_events") == (row.events < minimum[row.model])
    fitted = s.loc[s.status.eq("fitted")]
    assert (fitted.CI95_lower <= fitted.HR_per_expression_SD).all()
    assert (fitted.HR_per_expression_SD <= fitted.CI95_upper).all()
    key = s.set_index(["cohort", "endpoint", "gene", "model"])
    for model in ("univariable", "age_stage_adjusted"):
        wt1 = key.loc[("BRCA", "PFI", "WT1", model)]
        assert wt1.status == "fitted" and 1.3 < wt1.HR_per_expression_SD < 2.0 and wt1.BH_FDR < 0.05
    # PFI rows with time 0 are excluded from n
    assert key.loc[("BRCA", "PFI", "WT1", "univariable"), "n"] < 320


def test_summary_json_lists_significant_rows(pipeline_outputs):
    root, tables, _ = pipeline_outputs
    summary = json.loads((root / "analysis_summary.json").read_text())
    assert set(summary) == {"expression_significant", "PFI_significant",
                            "TNBC_pathway_significant", "TNBC_RPPA_significant"}
    assert {r["gene"] for r in summary["expression_significant"]} >= {"PHGDH"}
    assert "WT1" in {r["gene"] for r in summary["PFI_significant"]}
    assert tables["summary"] == summary


def test_legacy_mode_cannot_silently_disable_current_qc(synthetic_raw, tmp_path):
    raw, _ = synthetic_raw
    with pytest.raises(ValueError, match="[Ll]egacy"):
        m.run(tmp_path, raw, figures=False, legacy=True)
    assert not (tmp_path / "analysis_summary.json").exists()


def test_unresolved_duplicated_candidate_symbol_stops_analysis(tmp_path):
    raw = tmp_path / "raw"
    synthetic_tcga.make_dataset(raw, n=120, seed=1, duplicate_candidate="WT1")
    rna = m.load_rna(raw)
    alias, maps = m.gene_coverage(rna)
    assert {x["rna_coverage"] for x in maps if x["gene"] == "WT1"} == {"unavailable_or_ambiguous"}
    clinical, _ = m.build_clinical_cohort(raw)
    clinical = m.attach_purity(clinical, m.load_purity(raw))
    with pytest.raises(ValueError, match="Ambiguous required RNA symbol WT1"):
        m.prepare_expression(rna, clinical, alias, maps)


def test_missing_purity_is_an_error(tmp_path):
    raw = tmp_path / "raw"
    synthetic_tcga.make_dataset(raw, n=60, seed=2)
    (raw / "purity.tsv").unlink()
    with pytest.raises(FileNotFoundError, match="purity"):
        m.load_purity(raw)
    (raw / "purity.tsv").write_text("Sample\tPurity\nTCGA-ZZ-9999-01A\t0.5\n")
    clinical, _ = m.build_clinical_cohort(raw)
    with pytest.raises(ValueError, match="No purity values matched"):
        m.attach_purity(clinical, m.load_purity(raw))


def test_derived_purity_fallback_warns(tmp_path):
    raw = tmp_path / "raw"
    synthetic_tcga.make_dataset(raw, n=60, seed=3)
    (raw / "purity.tsv").rename(raw / "purity_archive_derived.tsv")
    with pytest.warns(RuntimeWarning, match="derived archive purity"):
        purity = m.load_purity(raw)
    assert list(purity.columns) == ["SAMPLE_ID", "PURITY"]


def test_cli_runs_and_writes_figures(synthetic_raw, tmp_path):
    raw, _ = synthetic_raw
    assert m.main(["--root", str(tmp_path), "--raw", str(raw), "--log-level", "WARNING"]) == 0
    figures = sorted(p.name for p in (tmp_path / "figures").glob("*.png"))
    assert figures == ["TCGA_F01_expression_contrasts.png", "TCGA_F02_CNA_prevalence.png",
                       "TCGA_F03_TNBC_pathway_associations.png", "TCGA_F04_PFI_forest.png",
                       "TCGA_F05_TNBC_RPPA_associations.png"]
