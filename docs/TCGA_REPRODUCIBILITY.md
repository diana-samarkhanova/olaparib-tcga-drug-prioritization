# TCGA and published resistance-context reproduction

These scripts reproduce the targeted TCGA-BRCA analysis and extract published
GSE165914 differential-expression estimates. The repository stores code and
provenance; restore the analysis inputs before running the commands below.

## Environment

Use Python 3.12 and install the archived direct-dependency versions:

```bash
python -m pip install -r analysis/tcga_extended/requirements.txt
```

The list pins pandas 2.2.3, numpy 2.3.5, scipy 1.17.0, lifelines 0.30.3,
matplotlib 3.10.8, and openpyxl 3.1.5. It is not a complete transitive lockfile.

## TCGA inputs

Restore `analysis/tcga_extended/raw/` and its source manifests from the analysis
archive. Required inputs are the clinical patient, primary-sample, and Firehose
receptor-status tables; TCGA-CDR supplemental table; RNA and CNA matrices;
MSigDB Hallmark 2025.1.Hs definitions; RPPA matrix; mutation case list; and the
per-gene mutation JSON snapshots. The compact RNA and CNA inputs in the archive
are accepted when their full matrices are absent. Their hashes are recorded in
`computation_subset_manifest.json`; the RNA subset retains the required genes
for all eleven signatures and the candidate panel.

The archive omits the original `raw/purity.tsv`, although the reported adjusted
pathway analysis used purity for 1,021 of 1,082 matched patients. **Do not omit
purity when reproducing that analysis.** If the original source is unavailable,
recover the archived cohort's values from its saved clinical intermediate:

```bash
python analysis/tcga_extended/reconstruct_archived_purity.py
```

This creates `raw/purity_archive_derived.tsv` and an explicit provenance JSON
with source/output hashes. It does not recreate the full original GDC file.
Run it immediately after restoring `matched_clinical_cohort.csv`, before
overwriting that intermediate with a new run. Original `raw/purity.tsv` takes
precedence; the derived fallback prints a warning. Missing or unusable purity
now stops analysis instead of silently changing adjustment covariates.

For source retrieval rather than archive restoration:

```bash
python analysis/tcga_extended/download_sources.py
```

The downloader verifies SHA-256 before replacing local files and stops on a
changed snapshot. It retrieves all used source entries, including the full RNA
and CNA matrices (more than 200 MB combined); it does not need to run for an
offline compact-input rerun. `--include-unused` also retrieves provenance-only
sources. Live URLs are not guaranteed to keep serving the recorded snapshots.

## Run order

Run from the repository root after restoring inputs:

```bash
python analysis/tcga_extended/reconstruct_archived_purity.py
python analysis/tcga_extended/analyze_tcga.py
python analysis/tcga_extended/analyze_mutations_and_sensitivity.py
python analysis/astra_extended/build_synthesis.py
python analysis/astra_extended/analyze_gse165914.py
python analysis/astra_extended/integrate_context.py
```

Skip the first command if using the original purity input. The main TCGA script
creates its tables and figures; mutation analysis adds mutation and PAM50
survival-sensitivity results. To redraw existing TCGA tables without refitting:

```bash
python analysis/tcga_extended/plot_saved_results.py
python analysis/tcga_extended/plot_mutations_saved.py
```

GSE/context integration additionally requires:

- `analysis/mechanisms_extended/evidence.csv` for the curated synthesis;
- `analysis/astra_extended/raw/GSE165914_M468P_OR_Diff_Expression_mshekhar.csv.gz`;
- `analysis/astra_extended/raw/GSE165914_SUM_M468_Parental_OR_counts_mshekhar.csv.gz`;
- `tables/Selected_CRISPR_panel_20_genes.csv` restored from the archive's top-level tables.

GEO source URLs are recorded in the archived `GSE165914_metadata.txt`.
`build_synthesis.py` writes curated mechanistic interpretation, not a learned
ranking. It overwrites the priority matrix, so rerun extraction and integration
after it to restore the appended contextual columns.

## Expected checks and interpretation

The archived analysis has 1,082 unique RNA-matched primary patients: 157 clinical
TNBC, 848 non-TNBC, and 77 unclassified. PAM50 basal is a separate classification
(171 patients). Nineteen coding candidates are represented; LINC03007 is absent
from the selected RNA/CNA platform. LAMTOR2 uses the source alias ROBLD3.

Expected TCGA table sizes are 570 pathway-association rows, 342 main Cox rows,
57 PAM50-adjusted sensitivity rows, and 69 mutation-prevalence rows including
BRCA1, BRCA2, and TP53 context. Preserve missing values and model status flags.

The packaging audit verified archived source hashes and the compact-input
hashes. All 114 available GSE author log2FC/P/FDR values matched the supplied
source table exactly. These checks are not a claim that a new full Cox analysis
was run during packaging.

TCGA associations describe patient-tumour context, not olaparib response or
gene-specific drug sensitivity. GSE165914 extraction preserves the authors'
genome-wide differential-expression results; it does not refit DE models from
the supplied fractional processed counts. Continuous 15-micromolar olaparib
exposure and acquired resistance are confounded in that comparison. These
observations do not establish treatment efficacy in the user's edited models.
