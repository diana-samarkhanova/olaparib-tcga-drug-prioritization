# CRISPR summaries, curated mechanisms and figures

These scripts reproduce the project-specific analysis of supplied MAGeCK RRA
gene summaries. They do not reprocess sequencing reads, fit MAGeCK models, or
recover sample/replicate metadata absent from the supplied workbooks. Restore
the archived inputs as described in the root README before running commands.
Run commands from the repository root.

## Screen summary analysis

```bash
python analysis/crispr/analyze_gene_summaries.py
```

The script expects the following files in `upload/`, with their original names:

- `KO_RRA_T16_OL_vs_DMSO_NEWcnv.gene_summary (1).xlsx`
- `KO_RRA_T24_OL_vs_DMSO_NEWcnv.gene_summary (2).xlsx`
- `SAM_RRA_T14_OL_vs_DMSO_NEWcnv_outlier_removed.gene_summary (3).xlsx`
- `SAM_RRA_T23_OL_vs_DMSO_NEWcnv_outlier.gene_summary (2).xlsx`
- `L1021-DiscoveryProbe-FDA-approved-Drug-Library (1).xlsx`

Each screen workbook must contain exactly one worksheet with the complete
14-column MAGeCK gene-summary schema. The catalog's `Chemical Data` sheet is
required by the vendor annotation stage at the end of the same script. Catalog
text matches are annotations, not evidence of direct target binding or drug
response. Raw workbooks and generated result tables are not versioned here.

Outputs are written into `analysis/crispr/`. Positive recurrent hits require
the supplied positive-selection FDR < 0.05 and positive LFC at both timepoints.
The original genome-wide FDR values are retained; timepoint P values are not
pooled. Repeated timepoints are not independent biological validation.

## Reactome enrichment: missing historical input

```bash
python analysis/crispr/enrichment/reactome_ora.py
```

This command additionally requires
`analysis/crispr/enrichment/ReactomePathways.gmt.zip`. The deliverable archive
omitted the ZIP, although it retained the ORA result tables and provenance.
The script checks the snapshot hash before running:

| Item | Archived provenance |
| --- | --- |
| Download URL | <https://reactome.org/download/current/ReactomePathways.gmt.zip> |
| Retrieved UTC | 2026-10-12T07:54:07.085104+00:00 |
| SHA-256 | `8c1dbc8578431da5d2d5118262718c60b553a9be3398e93658daa069e4a9afd4` |
| Bytes | 298479 |

The `current` URL can change. A new release cannot silently substitute for
the historical input. Exact recomputation is blocked until a matching file is
available; plotting the archived ORA table is available after input restoration.
The script also reads restored `download_metadata.json` and
`reactome_version_response.txt` from the same directory.

ORA uses the full recurrent foreground and the assayed, Reactome-annotated
background, with BH adjustment across every eligible pathway, including
zero-overlap pathways. The archived KO foreground contains 24 input genes,
20 mapped genes and 10,831 annotated background genes. The candidate panel is
not substituted for the full foreground.

## Literature curation

`analysis/mechanisms_extended/evidence.csv` is a restored, manually curated
56-row input. `derl2_nme6_evidence.csv` is also a curated input. No generator for
their final consolidation was present in the source archive. Do not describe
these files as an automated literature search or new quantitative experiment.

```bash
python analysis/mechanisms_extended/build_mechanism_evidence.py
python analysis/mechanisms_extended/build_activation_core.py
```

These commands regenerate the component `mechanism_review_rows.csv` and
`activation_core_evidence.csv` tables from the literature annotations embedded
in the scripts. They do not update the consolidated `evidence.csv` automatically.
Publication URLs, DOIs and PMIDs are the reusable references. Historical
`verification_ref` values are internal audit labels, not public citations.
Regulatory statements retain the original analysis date and require a fresh
status check before a later prospective drug shortlist is used.

## Project figures

After restoring archived results, these commands regenerate the figures:

```bash
python analysis/extended_root/plot_project_evidence.py
python analysis/extended_root/plot_individual_expression.py
python analysis/extended_root/plot_published_combinations.py
python analysis/extended_root/plot_reactome.py
```

Outputs are written into `outputs/TCGA_Drug_Analysis_EN_20261013/`. The first
command combines archived CRISPR and drug-resource tables and creates
`tables/Selected_CRISPR_panel_20_genes.csv`, required by the report builders.
The individual expression plots use TCGA expression, clinical and contrast
tables. The combination plot uses the published-screen extraction and curated
combination table. The Reactome figure uses the saved significant KO ORA table;
it does not rerun enrichment.

## Optional report export

```bash
python -m pip install -r analysis/report_extended/requirements.txt
python analysis/report_extended/prepare_report_content.py
python analysis/report_extended/build_report.py
```

These scripts format the historical English report after all result tables,
curated text and figure assets have been restored/generated. They contain
snapshot-specific prose and assertions, so they are not a generic report engine
for arbitrary new cohorts. The DOCX builder uses the public `python-docx` package.

`prepare_workbook_payload.py` optionally writes the historical workbook's JSON
payload and integrated summary table. The original workbook renderer used the
environment-specific `@oai/artifact-tool` Node package and is excluded from this
portable Python repository. The archived XLSX and ordinary CSV result tables
remain the records of that export. Workbook rendering is not required for any
scientific calculation. The former `package_analysis.py` deliverable bundler is
also excluded because it only copied files and assembled ZIP archives.

## Validation during repository preparation

The four supplied screen workbooks were independently reread and recurrent
positive calls recomputed from the original FDR/LFC columns. The results matched
the archived flags exactly: 24 KO genes and 3 SAM genes. A complete rerun of the
portable screen-summary script also reproduced the archived KO table (19,112
rows), SAM table (18,886 rows), long table (75,996 rows), and candidate-axis table
(46 rows) exactly. The Reactome figure was regenerated from the archived table,
and the missing-snapshot diagnostic was checked.

Portability changes alter filesystem roots, input diagnostics and
output-directory creation. The activation curation generator also incorporates
the archived PHLDA2 reference's distinction between its 2015 online publication
and 2018 issue date; this aligns the generator with the saved curation table.
Scientific calculations are unchanged.
