# Olaparib resistance: TCGA context and drug prioritization

Research code for integrating CRISPR knockout/activation screen summaries,
TCGA-BRCA tumour data, published resistant-cell expression results, and public
drug-response resources to prioritize follow-up experiments with an FDA-approved
drug library. This repository preserves the analysis accompanying
`Complete_Analysis_EN.zip`, with portability fixes and explicit input checks.

**This is a computational prioritization project.** It does not report a completed
wet-lab drug screen on the project's olaparib-resistant or edited cell models.
Public parental-cell responses, patient-tumour associations, published combination
data and mechanistic hypotheses remain separate evidence types.

## Scope

- Knockout panel: **LAMTOR2, WDR24, DERL2, NME6, TP53BP1, PARG**.
- Activation panel: **TP53I11, FOXQ1, CERK, PHF19, SLC35D1, PHLDA2, WT1,
  HEPHL1, LINC03007, SLC47A1, ELL2, INHBA, PHGDH, ABCB1**.
- Catalog: **L1021 DiscoveryProbe**, 2,321 named entries. The analysis is specific
  to the supplied catalog, not every product sold as an FDA-approved library.
- TCGA after conservative CDR-flag quality control: 1,077 RNA-matched primary
  patients; 156 clinical TNBC, 844 non-TNBC, 77 unclassified. PAM50 basal (170)
  is a separate classification. The archived 1,082-patient analysis is retained
  as the pre-audit comparison, not overwritten. This is not an untreated cohort
  or a cohort stratified by olaparib response.
- Public resources: TCGA/cBioPortal, TCGA-CDR, GSE165914, Reactome, PRISM,
  GDSC release 8.5, Drugs@FDA, and source-linked literature curation.

## Installation

Use Python 3.12. The requirements pin the direct packages used by the analysis;
they are not a full transitive environment lock.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Optional DOCX rendering requires
`python -m pip install -r analysis/report_extended/requirements.txt`.

For the synthetic-data and regression test suite:

```bash
python -m pip install -r requirements-dev.txt
python -m pytest
```

The tests do not require patient records or the proprietary catalog. They check
statistical calculations, planted synthetic signals, output contracts and known
regressions; passing them does not validate olaparib-response prediction.

## Restore the analysis workspace

The repository contains source code, provenance manifests and documentation.
User workbooks, patient-level matrices, vendor catalogs, generated figures and
result tables are excluded from version control. Keep the original analysis
archive separately; it contains both source inputs and saved results.

```bash
python scripts/restore_inputs.py --archive /path/to/Complete_Analysis_EN.zip
python analysis/tcga_extended/reconstruct_archived_purity.py
```

The restore helper retains existing files and never replaces Python/R/JavaScript
source. It also restores saved output tables; restoring them is **not** a rerun.
The purity helper recovers the archived cohort's purity covariate from its saved
clinical table because the original purity download was omitted from the archive.
It writes a separately identified derived input and provenance.

## Run the analyses

These are independent stages with explicit prerequisites, not a single automated
end-to-end workflow. Run from the repository root after restoring inputs.

| Stage | Entry point | Instructions |
| --- | --- | --- |
| Supplied CRISPR gene summaries | `analysis/crispr/analyze_gene_summaries.py` | [CRISPR, enrichment and figures](docs/CRISPR_REPRODUCIBILITY.md) |
| TCGA expression, CNA, pathways and survival | `analysis/tcga_extended/analyze_tcga.py` | [TCGA and GSE reproduction](docs/TCGA_REPRODUCIBILITY.md) |
| Mutation and survival sensitivity analysis | `analysis/tcga_extended/analyze_mutations_and_sensitivity.py` | [TCGA and GSE reproduction](docs/TCGA_REPRODUCIBILITY.md) |
| Published resistance expression context | `analysis/astra_extended/analyze_gse165914.py` | [TCGA and GSE reproduction](docs/TCGA_REPRODUCIBILITY.md) |
| Catalog mapping, PRISM, GDSC and prioritization | `analysis/library/prepare_catalog.py` and drug modules | [Drug resources and run order](docs/DRUG_REPRODUCIBILITY.md) |
| Project figures | `analysis/extended_root/plot_project_evidence.py` and companion scripts | [Figure commands](docs/CRISPR_REPRODUCIBILITY.md#project-figures) |

For a short verified rerun of the screen-summary stage:

```bash
python analysis/crispr/analyze_gene_summaries.py
python analysis/library/prepare_catalog.py
```

For TCGA and context integration:

```bash
python analysis/tcga_extended/analyze_tcga.py
python analysis/tcga_extended/analyze_mutations_and_sensitivity.py
python analysis/astra_extended/build_synthesis.py
python analysis/astra_extended/analyze_gse165914.py
python analysis/extended_root/plot_project_evidence.py
python analysis/astra_extended/integrate_context.py
```

Historical directory names and the dated figure-output directory are retained
to preserve dependencies between scripts and archived tables.

## Reproducibility limits

- Exact Reactome ORA recomputation requires a historical GMT ZIP absent from
  the archive. The script checks its recorded hash; archived results can be plotted.
- Full GDSC fitting workbooks and some published-supplement/OCR intermediates
  are absent. Their extraction stages require additional files documented in
  the drug guide; restored tables support downstream analysis.
- The consolidated mechanism evidence table is a manually curated input.
  Component generators do not recreate its final consolidation automatically.
- Gene priorities and drug nominations are qualitative literature curation.
  TCGA results are appended as tumour context; they do not train a predictor,
  calculate olaparib resistance probabilities or determine the compound list.
- GSE165914 uses published differential-expression estimates. Continued olaparib
  exposure is confounded with acquired resistance; this code does not refit DE
  models from raw sequencing reads.
- Public drug measurements are parental-cell single-agent observations. They
  do not establish edited-model selectivity, combination synergy or restoration
  of olaparib sensitivity. Endpoint scales and source fit/range flags are retained.
- Name/synonym matching is not chemical-structure identity validation. Catalog
  labeling and FDA annotations are historical snapshots, not current approval guarantees.
- The original environment-specific Excel renderer is excluded; ordinary CSV
  tables remain usable without it. Report prose is specific to this snapshot.

## Validation performed during packaging

- Full portable CRISPR rerun reproduced four archived tables exactly: 19,112 KO
  rows, 18,886 activation rows, 75,996 long-format rows and 46 candidate-axis rows.
- Independent recurrent-hit checks matched 24 KO and 3 activation genes.
- Catalog export matched all 2,324 source rows (2,321 named entries).
- All 114 available published GSE estimates matched the supplied table exactly.
- Independent retrieval of original TCGA RNA, CNA and purity files matched
  their frozen hashes. All compact matrix cells and all archived purity values
  matched the original sources exactly.
- A full original TCGA rerun, including Cox models, reproduced the archived
  numerical outputs. The audited rerun applies conservative CDR-redaction QC,
  canonical QSOX1 mapping and complete-case rank adjustment; see the audit note
  for changed estimates and interpretation limits.
- Drug replay preserves source/dose context and restores the deterministic
  selinexor alias. Incompatible PRISM dose windows are not pooled.
- Regression checks cover stale integration inputs, incompatible dose windows,
  gene/timepoint selection and TCGA quality-control rules.

See [provenance and rights](docs/PROVENANCE.md) for source locations and the
changes made while preparing this repository.

See the [20 October code and scientific audit](docs/AUDIT_20261020.md) for
confirmed corrections, before/after results and the boundary between tumour
context and olaparib-response evidence.

The [refactor comparison](docs/REFACTOR_REVIEW_20261020.md) explains which changes
from the separately supplied TCGA refactor were adopted and why the scientific
audit corrections remain part of the default analysis.
