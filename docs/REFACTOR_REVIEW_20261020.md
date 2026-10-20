# Review of the supplied TCGA refactor

The supplied patch, Python file and ZIP describe the same TCGA refactoring
proposal. The patch records edits, the Python file is the replacement main
module, and the ZIP adds tests, configuration and documentation. They do not
contain a new olaparib-response dataset or a new drug-prioritization model.

The standalone Python file is byte-identical to its ZIP counterpart (SHA-256
`16564f1c010bd15555d1be96f2a9fe1654d669ef85fd446ee47f85da84672355`).
The supplied ZIP SHA-256 is
`9900ac21241ece221dec367496322eaf6c7dd6ec2749d9b51860519c7907a9c9`.

## Decision

Adopt the engineering improvements and combine them with the independently
audited scientific corrections. Relative to the original compact script, the
refactor is easier to inspect, test and reuse. Relative to the audited pipeline,
the supplied bundle alone is incomplete.

| Area | Supplied refactor | Integrated decision |
| --- | --- | --- |
| Code organization | Named functions, importable module and explicit entry point | Adopt; imports must not start an analysis |
| Input and output locations | CLI options for raw inputs, table outputs and figures; figures can be skipped | Adopt and document |
| Cox diagnostics | Preserve convergence and collinearity warnings in an output column | Adopt; warnings are diagnostics, not proof of a reliable fit |
| Statistical checks | Synthetic examples, BH reference checks, planted effects and output contracts | Adopt, retain existing regressions and add scientific-QC coverage |
| Missing-data ranking | Rank continuous covariates after complete-case selection | Retain the audited rule |
| Ambiguous identifiers | Exclude duplicated candidates, but pathway rows still use their first occurrence | Retain explicit QSOX1 Entrez 5768 resolution and fail on unresolved required-symbol ambiguity |
| CDR flags | Merged but not excluded | Retain declared conservative QC and report exclusions |
| Drug and integration fixes | Not included; ancillary changes remove unused imports | Preserve deterministic selinexor mapping, dose-window separation and current-input checks |
| Historical mode | Partial reversions controlled by a legacy flag | Use the archived commit for historical reproduction; do not expose a flag that can silently bypass current QC |

## What the tests establish

The supplied synthetic tests exercise calculations and recover effects deliberately
inserted into artificial data. Their original fixture does not contain the five
real CDR flags or the QSOX1 historical duplicate. Therefore it cannot by itself
detect those defects. Test counts must not be presented as clinical validation,
and agreement between two implementations on synthetic data does not establish
correctness of the real-data cohort or identifiers.

The integrated tests additionally cover the audited QC rules, deterministic
drug aliases and dose-context handling, explicit screen timepoints, and use of
the regenerated integration panel. Real-data reproduction is a separate audit
step, recorded in [AUDIT_20261020.md](AUDIT_20261020.md).

The supplied suite passed all 66 tests when run independently. The integrated
suite passed 87 tests; the original 14 audit regressions also passed through
their unittest entry point. These are software verification results, not
additional biological observations.

Three additional cases protect undefined partial correlations when inputs are
constant, fully explained by adjustment variables, or have no residual degrees
of freedom. This safeguard changes no real TCGA result. Scoped Ruff checks pass
for the integrated TCGA modules and tests; the wider legacy analysis still has
12 lint findings, so repository-wide lint cleanliness is not claimed.

## Biological interpretation

Neither version tests olaparib resistance in TCGA. The patient analysis provides
expression, copy-number, pathway, protein and prognosis context for the chosen
genes. Drug nominations remain qualitative literature-based hypotheses,
supplemented by public parental-cell drug measurements. A cleaner implementation
and stronger software tests improve reproducibility; they do not turn these
observations into evidence of edited-model drug efficacy.
