# Provenance and rights

Source: the project's `Complete_Analysis_EN.zip`, analysis snapshot from
12–13 October 2026. Repository preparation: 20 October 2026.

Source archive SHA-256:
`0ec9e3aa21e483a9db811281f5a39fba8cd031b14573ae27b413b2945ad0bd11`.

The scientific scripts retain their original analysis structure. Packaging
changes resolve filesystem roots relative to scripts, create missing output
directories, add explicit input/hash checks, and supply archive-restoration and
catalog-export helpers. One activation-reference date annotation was aligned
with the archived table. Missing purity is no longer silently omitted.
Statistical formulas, screen thresholds and drug-response endpoints were retained.

## Source manifests

- `analysis/tcga_extended/source_snapshot_manifest.json`: original URLs,
  SHA-256 values and snapshot dates; compact inputs have a separate manifest.
- `analysis/astra_extended/raw/GSE165914_metadata.txt`: GEO supplementary-file
  URLs and study metadata.
- `analysis/crispr/enrichment/download_metadata.json`: historical Reactome GMT
  URL and checksum.
- `analysis/prism/*provenance.json`: source release and measurement coverage.
- `analysis/drugs_extended/public_crosscheck/gdsc_download_provenance.json`:
  GDSC source workbooks and checksums.
- [Drug guide](DRUG_REPRODUCIBILITY.md): resource URLs, prerequisites and
  source-specific limitations.

## Source data and third-party code

Raw and processed biological measurements, vendor workbooks and publication
supplements are not redistributed in this code repository. Their original terms
continue to apply when retrieved or restored for local use. Restored source and
result files are ignored by Git.

The upstream PRISM R scripts in the original archive are reference material,
not executed pipeline stages. They are not redistributed here; see
<https://github.com/broadinstitute/repurposing>. The source dataset's license
must not be assumed to license separate scripts.

This private project does not apply a new open-source license. A suitable
code license and source-data redistribution terms can be chosen before a public
release. Literature identifiers and URLs embedded in curation scripts identify
evidence sources; internal verification labels are not public citations.
