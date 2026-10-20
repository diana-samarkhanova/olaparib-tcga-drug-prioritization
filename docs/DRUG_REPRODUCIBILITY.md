# Drug-resource analysis: inputs and run order

Run the commands below from the repository root after installing the Python
requirements. These are source-specific analyses of public pharmacology and a
supplied catalog. They do not measure drug responses in the user's edited cells.
Scripts write derived tables beside their source scripts, so keep a separate
copy of the delivered results when comparing a rerun.

## Two reproduction routes

The private `Complete_Analysis_EN.zip` delivery contains most small inputs,
public-data extracts, curated tables, and historical results under
`reproducibility/`. Restore its non-code contents using the repository's restore
instructions. That route supports downstream recomputation without downloading
the large public matrices again. The archive and uploaded workbooks are not
tracked in this code repository.

Starting from public sources instead requires the downloads and prerequisites
below. There is no fully automatic fresh-download workflow for every stage.
In particular, the GDSC extractor reads already downloaded workbooks, FDA data
must be supplied separately, and the published-supplement extraction has missing
OCR intermediates in the delivered archive. The restored result tables remain
usable for downstream integration despite that last limitation.

## 1. Catalog and FDA snapshot

Place the supplied workbook at:

```text
upload/L1021-DiscoveryProbe-FDA-approved-Drug-Library (1).xlsx
```

The exact sheet name is `Chemical Data`. This is the L1021 catalog used in the
analysis, not the separate HYL022P catalog. Its original 2,324 rows contain 2,321
nonempty `Item Name` entries. Export the complete sheet for the alias stage:

```bash
python analysis/library/prepare_catalog.py
```

This export preserves the original rows and columns. Downstream scripts apply
their own nonempty-name filters. `--input` and `--output` allow alternative paths
for this export; the existing analysis stages expect the default paths.

The FDA annotation stage requires `Products.txt`, `MarketingStatus.txt`,
`MarketingStatus_Lookup.txt`, and `Submissions.txt` in `analysis/fda/`.
The delivery retains the snapshot used on 12 October 2026. The official source
is [Drugs@FDA data files](https://www.fda.gov/drugs/drug-approvals-and-databases/drugsfda-data-files).
Then run:

```bash
python analysis/fda/annotate_fda.py
```

The script's date metadata describes that historical snapshot. If using newly
downloaded FDA files, update the provenance and independently review the curated
regulatory statements before interpreting the output as current. A matched
ingredient/product record does not establish a present indication, exact catalog
form, or approval of a combination with olaparib.

## 2. PRISM original release

Numerical results use the [PRISM 19Q4 dataset, version 4](https://doi.org/10.6084/m9.figshare.9393293.v4).
They do not use 20Q2 or 24Q2 numerical measurements. The source dataset metadata
records CC BY 4.0 and credits Broad DepMap, Steven Corsello, Mustafa Kocak, and
Todd Golub. Source methods: [Corsello et al.](https://doi.org/10.1038/s43018-019-0018-6).

| Input | Figshare file ID |
| --- | --- |
| Primary replicate-collapsed log2 fold-change matrix | 20237709 |
| Primary treatment metadata | 20237715 |
| Primary cell metadata | 20237718 |
| Secondary fitted curve parameters | 20237739 |
| Secondary treatment metadata | 20237763 |
| Secondary cell metadata | 20237769 |

Each download URL is `https://ndownloader.figshare.com/files/FILE_ID`.
The download script streams the large matrices and retains breast-model rows.
Skip it if restoring the matching inputs from the delivery archive.

```bash
python analysis/prism/download_sources.py
python analysis/prism/map_library.py
python analysis/prism/summarize_curves.py
python analysis/prism/summarize_primary.py
```

Run FDA annotation first if its columns are needed in the PRISM summaries.
The summary implementation uses pandas `GroupBy.apply(include_groups=False)`;
use the repository's pinned environment rather than an older pandas release.

Primary log2 fold change and its derived `2**LFC` readout are not IC50 values.
Secondary IC50 is already in micromolar, whereas AUC retains the source dose
window. MTS010 technical redos are preferred when available; otherwise the
script retains the available records. A median AUC or IC50 is reported only
when selected dose windows agree within 0.1% (to accommodate source rounding).
When windows differ or are missing, pooled AUC, R² and IC50 are left missing and
the individual records remain in `selected_curve_details_JSON` and
`secondary_compound_cell_dose_context_summary.csv`. Extrapolated or missing
IC50 values remain flagged. Within-range fitted IC50 does not by itself imply
that a fit is reliable; inspect R² and the retained curve diagnostics.

The October 2026 audit corrected an earlier pooling error. For selinexor in
MDA-MB-468, HTS002 tested up to 10 µM and fitted IC50 0.509477 µM, whereas MTS006
tested up to 2.5 µM and fitted IC50 0.020076 µM. The old pooled IC50 0.264777 µM
and AUC 0.612029 should not be used. The corrected outputs retain the two source
contexts without assigning a combined potency. Three compounds (selinexor,
doxycycline and U-0126) have materially different selected windows across 21
breast cell lines each; the most pronounced window difference is selinexor.

## 3. GDSC release 8.5

First complete the PRISM download and catalog-mapping stages. They supply
`secondary_breast_metadata.csv` and `library_prism_name_coverage.csv`.
Download these two workbooks into `analysis/drugs_extended/public_crosscheck/`:

- [GDSC1_fitted_dose_response_27Oct23.xlsx](https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC1_fitted_dose_response_27Oct23.xlsx)
- [GDSC2_fitted_dose_response_27Oct23.xlsx](https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC2_fitted_dose_response_27Oct23.xlsx)

The delivery omitted these workbooks and retained the breast-model extracts.
Historical SHA-256 values are:

```text
GDSC1: 837b0686500fde75179e490de08f034abd9f882d8b0253d637bafe83e156dafd
GDSC2: f950a7027be265f8a7a74220a27fd18cbd368485349bd8c2048e88bb1cd07560
```

```bash
python analysis/drugs_extended/public_crosscheck/extract_gdsc.py
```

This command does not download the files. It computes IC50 in micromolar as
`exp(LN_IC50)`, retains source drug IDs and dose windows, and flags extrapolation.
Its regenerated provenance records the execution date under the legacy field
`download_date_utc`; retain the original manifest for the historical download
date. GDSC1 and GDSC2 stay separate. Source methods and field definitions:
[Sanger documentation](https://depmap.sanger.ac.uk/documentation/datasets/drug-sensitivity/)
and [release field description](https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC_Fitted_Data_Description.pdf).

## 4. Catalog mining and mechanism integration

After stages 1–3, run:

```bash
python analysis/drugs_extended/curate_aliases.py
python analysis/drugs_extended/mine_catalog.py
python analysis/drugs_extended/join_mechanisms.py
python analysis/drugs_extended/finalize_priorities.py
```

`curate_aliases.py` uses `analysis/library/library_raw.csv` and retrieves selected
PubChem synonyms unless the corresponding JSON responses are already cached in
`analysis/drugs_extended/identity_sources/`. Source links and lookup status are
retained. Synonym expansion can differ if a live resource changes; use the
archived cache for the historical mapping.
The builder now explicitly includes KPT-330/selinexor (B1464). This reviewed
alias was present in the original result table but missing from its builder,
so an earlier fresh run would silently discard the selinexor links.

`join_mechanisms.py` uses the reviewed
`analysis/mechanisms_extended/evidence.csv` restored from the delivery. Its
fallback expects all three of `mechanism_review_rows.csv`,
`derl2_nme6_evidence.csv`, and `activation_core_evidence.csv`. Two builder scripts
exist for parts of that curation, but they are not a substitute for the final
reviewed evidence table; the DERL2/NME6 curated table must also be supplied.

Names, developmental codes, and salt heuristics provide provisional chemical
links. No canonical structure validation was performed. Coverage counts are
catalog-entry counts rather than counts of unique chemical entities. Drug
nomination remains evidence-based curation, not a fitted TCGA response predictor.
The drug profile tables retain PRISM source-curve JSON and dose-compatibility
flags. Full catalog-linked contexts are also exported to
`prism_secondary_catalog_dose_context_records.csv` and its MDA-MB-468 subset;
TNBC context summaries preserve screen, batch, and dose window.

## 5. Published combination evidence

The delivered `analysis/olaparib_screen/L1021_screen_overlap.csv` can be used
directly after restoring the archive. With that table and
`analysis/drugs_extended/whole_catalog_audit.csv` present, run:

```bash
python analysis/combinations_extended/build_combo_evidence.py
```

The source extraction script is retained for transparency:

```bash
python analysis/olaparib_screen/extract_yusoh.py
```

However, it requires files not included in the delivery:

- `MOL2-19-3387-s001.docx` from the [published supplement](https://pmc-oa-opendata.s3.amazonaws.com/PMC12591312.1/MOL2-19-3387-s001.docx);
- the original `PMC12591312.1.json` metadata containing the supplement MD5;
- `S2_block1.txt` through `S2_block4.txt`, the original OCR text intermediates.

The supplement's historical MD5 is `c2de36e3f764e6362ecf10c35dbd1719`.
Recreating the OCR files requires matching the original block segmentation and
row order and checking values against the raster image. The script documents
three manual corrections and six transcribed final-block rows. The missing text
must not be synthesized from the already processed table and presented as an
independent extraction. The original article is
[Yusoh et al.](https://doi.org/10.1002/1878-0261.70093).

These are author-reported combination-index and IC50 summaries in parental
MDA-MB-231, not measurements in the user's MDA-MB-468 derivatives. The code does
not reconstruct raw viability, recompute Bliss/ZIP synergy, or replace censored
IC50 values with point estimates.

## Upstream R references

The delivered archive included two upstream scripts used to understand PRISM
processing. They were not executed to generate this project's results and are
linked rather than redistributed here:

- [secondary_processing_script.R](https://github.com/broadinstitute/repurposing/blob/master/secondary_processing_script.R), credited to Mustafa Anil Kocak, 2019;
- [depmap_reprocessing_script_for_secondary.R](https://github.com/broadinstitute/repurposing/blob/master/depmap_reprocessing_script_for_secondary.R), for the separate 20Q2 reprocessing workflow.

The latter contains a user-specific download path and batch variable. Neither
is a required stage of the Python workflow above. The CC BY metadata for the
Figshare datasets should not be assumed to license separately hosted code.

## Interpretation boundaries

PRISM and GDSC provide public parental-cell monotherapy data. These endpoints
do not demonstrate olaparib resensitization, synergy, or genotype-selective
killing in the user's KO or activation models. Do not average native AUCs across
sources or across incompatible concentration windows. Mechanistic hypotheses,
published combination evidence, TCGA associations, and public monotherapy
activity remain separate evidence layers. No numerical CTRP analysis was
completed in the archived project.
