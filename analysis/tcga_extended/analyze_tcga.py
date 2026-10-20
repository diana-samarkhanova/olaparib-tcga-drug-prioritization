#!/usr/bin/env python3
"""Targeted TCGA-BRCA context analysis for the olaparib-resistance candidate panel.

The script joins cBioPortal PanCanAtlas clinical, RNA, CNA and RPPA snapshots
with TCGA-CDR survival endpoints and produces, for the 20-gene KO/CRISPRa panel:

* expression contrasts (clinical TNBC vs non-TNBC; PAM50 basal vs other),
* GISTIC copy-number prevalence and dosage-expression correlation,
* covariate-adjusted partial rank correlations with Hallmark programme scores,
* unadjusted RNA-RPPA correlations within TNBC,
* continuous-expression Cox models (PFI primary, DFI supporting, OS exploratory),
* five figures and two JSON summaries.

Output file names and column names are a contract with downstream scripts
(``analyze_mutations_and_sensitivity.py``, ``astra_extended/integrate_context.py``,
``report_extended/*``) and must not change without updating them.

Run from any directory::

    python analysis/tcga_extended/analyze_tcga.py [--raw DIR] [--root DIR]

The audited default excludes TCGA-CDR-redacted cases under a conservative
flag-based QC policy, resolves required duplicate RNA symbols explicitly, and
ranks variables on the same complete cases. Historical output is reproduced
with its archived code and inputs, not by partially reverting this pipeline.
"""
from __future__ import annotations

import argparse
import json
import logging
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

if __package__:
    from .statistics_qc import (REDACTION_POLICY, exclude_cdr_redactions,
                               resolve_required_rna_symbols, partialrank)
else:
    from statistics_qc import (REDACTION_POLICY, exclude_cdr_redactions,
                               resolve_required_rna_symbols, partialrank)

log = logging.getLogger("analyze_tcga")

# --------------------------------------------------------------------------- #
# Panel and analysis constants
# --------------------------------------------------------------------------- #
KO = ["LAMTOR2", "WDR24", "DERL2", "NME6", "TP53BP1", "PARG"]
ACT = ["TP53I11", "FOXQ1", "CERK", "PHF19", "SLC35D1", "PHLDA2", "WT1", "HEPHL1",
       "LINC03007", "SLC47A1", "ELL2", "INHBA", "PHGDH", "ABCB1"]
GENES = KO + ACT
PATHS = ["MTORC1_SIGNALING", "UNFOLDED_PROTEIN_RESPONSE", "OXIDATIVE_PHOSPHORYLATION",
         "GLYCOLYSIS", "DNA_REPAIR", "G2M_CHECKPOINT", "TGF_BETA_SIGNALING",
         "EPITHELIAL_MESENCHYMAL_TRANSITION", "APOPTOSIS", "HYPOXIA"]
E2F = "E2F_TARGETS"
SCORED_PATHWAYS = PATHS + [E2F]
ENDPOINTS = ["PFI", "DFI", "OS"]
# Previous-identifier fallbacks for symbols absent from the RNA/CNA platform.
RNA_ALIASES = {"LINC03007": ["LOC646588"], "LAMTOR2": ["ROBLD3"]}
RPPA_READOUTS = ["RPS6KB1|p70S6K_pT389", "RPS6|S6_pS235_S236", "RPS6|S6_pS240_S244",
                 "EIF4EBP1|4E-BP1_pT37T46", "MAPK1 MAPK3|MAPK_pT202_Y204",
                 "AKT1 AKT2 AKT3|Akt_pS473", "TP53BP1|53BP1", "RAD51|Rad51", "XBP1|XBP1"]
COX_MODELS = [("univariable", ["gene_z"], 15),
              ("age_stage_adjusted", ["gene_z", "age_per10", "stage_III_IV"], 30)]
CONTRASTS = [("clinical_TNBC_vs_nonTNBC", "TNBC_STATUS", "TNBC", "non-TNBC"),
             ("PAM50_Basal_vs_other", "BASAL_STATUS", "Basal", "Other PAM50")]
MIN_PARTIAL_N = 30
MIN_PURITY_N = 30
MIN_RPPA_N = 15
HALLMARK_GMT = "h.all.v2025.1.Hs.symbols.gmt"
MSIGDB_URL = "https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_"


def cohort_masks(clinical: pd.DataFrame) -> list[tuple[str, pd.Series]]:
    """Ordered analysis cohorts; BRCA is the pooled cohort."""
    return [("BRCA", pd.Series(True, index=clinical.index)),
            ("TNBC", clinical.TNBC_STATUS.eq("TNBC")),
            ("Basal", clinical.BASAL_STATUS.eq("Basal"))]


# --------------------------------------------------------------------------- #
# Statistics helpers
# --------------------------------------------------------------------------- #
def benjamini_hochberg(p) -> np.ndarray:
    """BH-adjusted p-values; non-finite inputs stay NaN and are not counted."""
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    v = p[ok]
    n = len(v)
    if n:
        order = np.argsort(v)
        ranked = v[order] * n / np.arange(1, n + 1)
        adjusted = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
        out[np.where(ok)[0][order]] = adjusted
    return out


def partial_rank_correlation(x: pd.Series, y: pd.Series, cov: pd.DataFrame,
                             rank_columns=(), min_n: int = MIN_PARTIAL_N):
    """Partial Spearman correlation of ``x`` and ``y`` given ``cov``.

    Rows with any missing value are dropped first. Columns in ``rank_columns``
    are then rank-transformed within the retained rows; the remaining columns
    (for example, indicator dummies) enter unchanged. Returns ``(n, rho, p)``;
    undefined residual correlations or insufficient observations return NaN.
    """
    return partialrank(x, y, cov, rank_columns=rank_columns, min_n=min_n)


# --------------------------------------------------------------------------- #
# Clinical cohort
# --------------------------------------------------------------------------- #
def read_clinical_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#", low_memory=False)


def classify_her2(ihc, fish) -> str:
    """Resolve IHC and FISH HER2 calls; FISH wins when only one is definitive."""
    definitive = ("Negative", "Positive")
    if ihc in definitive and fish in definitive and ihc != fish:
        return "Discordant"
    if fish in definitive:
        return fish
    if ihc in definitive:
        return ihc
    return "Unresolved"


def classify_tnbc(er, pr, her2_final) -> str:
    """Clinical TNBC requires three definitive negatives; any positive excludes."""
    calls = (er, pr, her2_final)
    if any(c == "Positive" for c in calls):
        return "non-TNBC"
    if all(c == "Negative" for c in calls):
        return "TNBC"
    return "Unclassified"


def classify_basal(subtype) -> str:
    if subtype == "BRCA_Basal":
        return "Basal"
    return "Other PAM50" if isinstance(subtype, str) else "Unclassified"


def stage_advanced(label) -> float:
    """1.0 for AJCC stage III/IV, 0.0 for stage I/II, NaN otherwise."""
    s = str(label).upper()
    if re.match(r"^STAGE IV", s) or re.match(r"^STAGE III", s):
        return 1.0
    if re.match(r"^STAGE II", s) or re.match(r"^STAGE I($|[AB])", s):
        return 0.0
    return np.nan


def build_clinical_cohort(raw: Path, rna_sample_ids=None) -> tuple[pd.DataFrame, int]:
    """One primary sample per patient with receptor status, PAM50 and CDR endpoints.

    Applies conservative CDR flag filtering before molecular analyses. Returns
    the eligible cohort and the original primary-patient count; aggregate QC
    metadata are stored in ``clinical.attrs['qc_metadata']``. Supplying RNA
    column identifiers additionally permits the pre-exclusion RNA count.
    """
    sample = read_clinical_table(raw / "data_clinical_sample.txt")
    patient = read_clinical_table(raw / "data_clinical_patient.txt")
    firehose = read_clinical_table(raw / "data_clinical_firehose_patient.txt")
    sample = (sample.loc[sample.SAMPLE_TYPE.eq("Primary")]
              .sort_values("SAMPLE_ID").drop_duplicates("PATIENT_ID"))
    receptor_cols = ["PATIENT_ID", "ER_STATUS_BY_IHC", "PR_STATUS_BY_IHC",
                     "IHC_HER2", "HER2_FISH_STATUS"]
    clinical = (sample.merge(patient, on="PATIENT_ID", validate="one_to_one")
                .merge(firehose[receptor_cols], on="PATIENT_ID", how="left",
                       validate="one_to_one"))
    clinical["HER2_FINAL"] = [classify_her2(a, b) for a, b in
                              zip(clinical.IHC_HER2, clinical.HER2_FISH_STATUS, strict=True)]
    clinical["TNBC_STATUS"] = [classify_tnbc(e, p, h) for e, p, h in
                               zip(clinical.ER_STATUS_BY_IHC, clinical.PR_STATUS_BY_IHC,
                                   clinical.HER2_FINAL, strict=True)]
    clinical["BASAL_STATUS"] = clinical.SUBTYPE.map(classify_basal)

    cdr = pd.read_excel(raw / "TCGA-CDR-SupplementalTableS1.xlsx", sheet_name="TCGA-CDR")
    cdr = cdr.loc[cdr.type.eq("BRCA")].rename(columns={
        "bcr_patient_barcode": "PATIENT_ID",
        "age_at_initial_pathologic_diagnosis": "CDR_AGE",
        "ajcc_pathologic_tumor_stage": "CDR_STAGE"})
    cdr_cols = ["PATIENT_ID", "CDR_AGE", "CDR_STAGE", "PFI", "PFI.time", "DFI", "DFI.time",
                "OS", "OS.time", "Redaction"]
    clinical = clinical.merge(cdr[cdr_cols], on="PATIENT_ID", how="left", validate="one_to_one")
    eligible, excluded = exclude_cdr_redactions(clinical)
    rna_excluded = (int(clinical.loc[excluded, "SAMPLE_ID"].isin(rna_sample_ids).sum())
                    if rna_sample_ids is not None else None)
    clinical = eligible
    clinical["STAGE_ADVANCED"] = clinical.CDR_STAGE.map(stage_advanced)
    clinical.attrs["qc_metadata"] = {
        "clinical_primary_patients_before_redaction_exclusion": len(sample),
        "CDR_redacted_primary_patients_excluded": int(excluded.sum()),
        "clinical_primary_patients": len(clinical),
        "redacted_RNA_profiled_patients_excluded": rna_excluded,
        "redaction_exclusion_policy": REDACTION_POLICY,
    }
    return clinical, len(sample)


def load_purity(raw: Path) -> pd.DataFrame:
    """Tumour purity keyed by 15-character sample ID; the derived file is a fallback."""
    path = raw / "purity.tsv"
    if not path.exists():
        path = raw / "purity_archive_derived.tsv"
        if not path.exists():
            raise FileNotFoundError(
                "Tumour purity is required to reproduce the reported adjusted analysis. "
                "Retrieve raw/purity.tsv with download_sources.py, or restore archived "
                "matched_clinical_cohort.csv and run reconstruct_archived_purity.py "
                "before rerunning analyze_tcga.py.")
        warnings.warn("Using explicitly derived archive purity input: " + str(path),
                      RuntimeWarning, stacklevel=2)
    pu = pd.read_csv(path, sep="\t")
    purity_col = [c for c in pu.columns if c.lower() == "purity"]
    sample_col = [c for c in pu.columns if c.lower() == "sample"]
    if len(purity_col) != 1 or len(sample_col) != 1:
        raise ValueError("Purity input must contain exactly one Sample and one Purity "
                         "column (case-insensitive).")
    pu["SAMPLE_ID"] = pu[sample_col[0]].astype(str).str[:15]
    pu = pu.sort_values(sample_col[0]).drop_duplicates("SAMPLE_ID")
    return pu[["SAMPLE_ID", purity_col[0]]].rename(columns={purity_col[0]: "PURITY"})


def attach_purity(clinical: pd.DataFrame, purity: pd.DataFrame) -> pd.DataFrame:
    clinical = clinical.merge(purity, on="SAMPLE_ID", how="left", validate="one_to_one")
    if clinical.PURITY.notna().sum() == 0:
        raise ValueError("No purity values matched the clinical cohort; refusing to "
                         "silently omit the adjustment.")
    return clinical


# --------------------------------------------------------------------------- #
# Expression, aliases and Hallmark scores
# --------------------------------------------------------------------------- #
def resolve_symbol(gene: str, symbols: set, aliases: dict = RNA_ALIASES):
    """Platform symbol for ``gene``: itself, a listed alias, or None."""
    if gene in symbols:
        return gene
    for alt in aliases.get(gene, []):
        if alt in symbols:
            return alt
    return None


def _entrez_int(value):
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    return int(f) if np.isfinite(f) and f == int(f) else None


def gene_coverage(rna_raw: pd.DataFrame, genes=GENES) -> tuple[dict, list[dict]]:
    """Map each panel gene to an RNA symbol and record whether it is usable.

    ``unavailable_or_ambiguous`` covers both absent and duplicated symbols.
    """
    symbols = set(rna_raw.Hugo_Symbol)
    alias = {g: resolve_symbol(g, symbols) for g in genes}
    maps = []
    for g, s in alias.items():
        rows = rna_raw.loc[rna_raw.Hugo_Symbol.eq(s)] if s else rna_raw.iloc[0:0]
        unique = len(rows) == 1
        maps.append({"gene": g, "source_symbol": s, "rna_rows": len(rows),
                     "entrez_gene_id": str(rows.iloc[0].Entrez_Gene_Id) if unique else "",
                     "rna_coverage": "available" if unique else "unavailable_or_ambiguous"})
    return alias, maps


def load_rna(raw: Path) -> pd.DataFrame:
    path = raw / "data_mrna_seq_v2_rsem.txt"
    if not path.exists():
        path = raw / "RNA_required_gene_subset.tsv.gz"
    rna = pd.read_csv(path, sep="\t", low_memory=False)
    rna["Hugo_Symbol"] = rna.Hugo_Symbol.fillna("")
    return rna


def numeric_matrix(table: pd.DataFrame) -> pd.DataFrame:
    """Symbol-indexed numeric matrix; the first row of a duplicated symbol is kept."""
    return (table.drop_duplicates("Hugo_Symbol").set_index("Hugo_Symbol")
            .drop(columns="Entrez_Gene_Id").apply(pd.to_numeric, errors="coerce"))


def prepare_expression(rna_raw: pd.DataFrame, clinical: pd.DataFrame, alias: dict,
                       maps: list[dict], legacy: bool = False, required_symbols=None):
    """Resolve required symbols and match samples; return log2 RNA and panel table.

    The pipeline supplies all candidate and scored-signature symbols. Standalone
    use retains unique background rows, required candidate rows, and the reviewed
    QSOX1 duplicate; unrelated ambiguous background symbols are not selected.
    """
    if legacy:
        raise ValueError("legacy=True is not supported; use the archived code and inputs for historical reproduction.")
    if required_symbols is None:
        unique = set(rna_raw.loc[~rna_raw.Hugo_Symbol.duplicated(keep=False), "Hugo_Symbol"])
        required_symbols = unique | {s for s in alias.values() if s is not None} | {"QSOX1"}
    selected = rna_raw.loc[rna_raw.Hugo_Symbol.isin(required_symbols)].copy()
    selected, resolution = resolve_required_rna_symbols(selected)
    rna = numeric_matrix(selected)
    rna = rna.loc[:, [s for s in clinical.SAMPLE_ID if s in rna.columns]]
    clinical = (clinical.loc[clinical.SAMPLE_ID.isin(rna.columns)]
                .set_index("SAMPLE_ID").loc[rna.columns])
    logrna = np.log2(rna + 1)
    usable = {m["gene"] for m in maps if m["rna_coverage"] == "available"}
    expr = {}
    for g, s in alias.items():
        if s is None or s not in logrna.index:
            continue
        if g not in usable:
            raise ValueError(f"Candidate {g} was not uniquely mapped in gene_coverage; resolve its stable identifier before analysis.")
        expr[g] = logrna.loc[s]
    expr = pd.DataFrame(expr)
    expr.index.name = "SAMPLE_ID"
    logrna.attrs["symbol_resolution"] = resolution.to_dict("records")
    return clinical, logrna, expr


def load_hallmark(path: Path) -> dict[str, list[str]]:
    sets = {}
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        sets[fields[0].removeprefix("HALLMARK_")] = fields[2:]
    return sets


def hallmark_scores(logrna: pd.DataFrame, sets: dict, pathways=SCORED_PATHWAYS):
    """Gene-wise z-scores across samples and mean-z programme scores per sample."""
    sd = logrna.std(axis=1, ddof=1).replace(0, np.nan)
    zrna = logrna.sub(logrna.mean(axis=1), axis=0).div(sd, axis=0)
    scores = pd.DataFrame({p: zrna.loc[[g for g in sets[p] if g in zrna.index]].mean(axis=0)
                           for p in pathways})
    coverage = pd.DataFrame([{"pathway": p, "source_genes": len(sets[p]),
                              "rna_matched_genes": sum(g in zrna.index for g in sets[p]),
                              "source_url": MSIGDB_URL + p} for p in pathways])
    return zrna, scores, coverage


def cohort_counts(clinical: pd.DataFrame, n_primary: int, maps: list[dict],
                  qc_metadata: dict | None = None) -> dict:
    qc = qc_metadata or clinical.attrs.get("qc_metadata", {})
    excluded_rna = qc.get("redacted_RNA_profiled_patients_excluded", 0)
    counts = {"clinical_primary_patients_before_redaction_exclusion": n_primary,
              "CDR_redacted_primary_patients_excluded": qc.get("CDR_redacted_primary_patients_excluded", 0),
              "clinical_primary_patients": qc.get("clinical_primary_patients", n_primary),
              "rna_matched_primary_patients_before_redaction_exclusion":
                  len(clinical) + excluded_rna if excluded_rna is not None else None,
              "redacted_RNA_profiled_patients_excluded": excluded_rna,
              "rna_matched_primary_patients": len(clinical),
              "redaction_exclusion_policy": REDACTION_POLICY,
              "TNBC_status": clinical.TNBC_STATUS.value_counts().to_dict(),
              "PAM50": clinical.SUBTYPE.fillna("Unclassified").value_counts().to_dict(),
              "HER2_final": clinical.HER2_FINAL.value_counts().to_dict(),
              "purity_available": int(clinical.PURITY.notna().sum()),
              "gene_coverage": maps, "endpoint_counts": []}
    for cohort, mask in cohort_masks(clinical):
        for ep in ENDPOINTS:
            d = clinical.loc[mask, [ep, ep + ".time"]].dropna()
            d = d.loc[d[ep + ".time"] > 0]
            counts["endpoint_counts"].append({"cohort": cohort, "endpoint": ep,
                                              "n": len(d), "events": int(d[ep].sum())})
    return counts


# --------------------------------------------------------------------------- #
# Analyses
# --------------------------------------------------------------------------- #
def expression_contrasts(expr: pd.DataFrame, clinical: pd.DataFrame) -> pd.DataFrame:
    """Mann-Whitney contrasts; BH over the available panel within each comparison."""
    rows = []
    for comparison, col, a, b in CONTRASTS:
        for g in GENES:
            if g not in expr:
                continue
            x = expr.loc[clinical[col].eq(a), g].dropna()
            y = expr.loc[clinical[col].eq(b), g].dropna()
            u, p = stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic")
            rows.append({"comparison": comparison, "gene": g,
                         "panel": "KO" if g in KO else "CRISPRa",
                         "n_group1": len(x), "n_group2": len(y),
                         "median_group1": x.median(), "median_group2": y.median(),
                         "median_delta_log2_RSEM_plus1": x.median() - y.median(),
                         "rank_biserial": 2 * u / (len(x) * len(y)) - 1, "p_value": p})
    out = pd.DataFrame(rows)
    out["BH_FDR"] = out.groupby("comparison").p_value.transform(benjamini_hochberg)
    return out


def load_cna(raw: Path) -> pd.DataFrame:
    path = raw / "data_cna.txt"
    if not path.exists():
        path = raw / "CNA_selected_gene_subset.tsv.gz"
    return pd.read_csv(path, sep="\t", low_memory=False)


def resolve_cna_aliases(cna_raw: pd.DataFrame, maps: list[dict]) -> dict:
    """CNA symbol per gene: direct symbol, else a unique Entrez ID match from RNA."""
    symbols = set(cna_raw.Hugo_Symbol)
    entrez = pd.to_numeric(cna_raw.Entrez_Gene_Id, errors="coerce")
    alias = {}
    for m in maps:
        g = m["gene"]
        if g in symbols:
            alias[g] = g
            continue
        eid = _entrez_int(m["entrez_gene_id"])
        hit = cna_raw.loc[entrez.eq(eid), "Hugo_Symbol"] if eid is not None else []
        alias[g] = str(hit.iloc[0]) if len(hit) == 1 else None
    return alias


def cna_analysis(cna_raw: pd.DataFrame, cna_alias: dict, clinical: pd.DataFrame,
                 expr: pd.DataFrame):
    """GISTIC prevalence per cohort and Spearman dosage-expression association."""
    cna = numeric_matrix(cna_raw)
    cnamat = pd.DataFrame(index=clinical.index)
    rows = []
    for g, s in cna_alias.items():
        if s is None or s not in cna.index:
            continue
        cnamat[g] = cna.loc[s].reindex(clinical.index)
        for cohort, mask in cohort_masks(clinical):
            v = cnamat.loc[mask, g].dropna()
            if g in expr:
                m = pd.concat([v, expr[g]], axis=1, keys=["cna", "expr"]).dropna()
                r, p = stats.spearmanr(m.cna, m.expr)
            else:  # copy-number prevalence is still reported without an RNA row
                m, r, p = v.iloc[0:0], np.nan, np.nan
            rows.append({"cohort": cohort, "gene": g, "n_profiled": len(v),
                         "deep_deletion_n": int((v == -2).sum()),
                         "shallow_loss_n": int((v == -1).sum()),
                         "diploid_n": int((v == 0).sum()), "gain_n": int((v == 1).sum()),
                         "amplification_n": int((v == 2).sum()),
                         "dosage_expression_n": len(m),
                         "dosage_expression_spearman_rho": r, "p_value": p})
    out = pd.DataFrame(rows)
    out["BH_FDR"] = out.groupby("cohort").p_value.transform(benjamini_hochberg)
    return cnamat, out


def covariate_frame(clinical: pd.DataFrame, scores: pd.DataFrame, cnamat: pd.DataFrame,
                    idx: pd.Index, cohort: str, gene: str, legacy: bool = False):
    """Adjustment covariates for one cohort/gene; returns the frame and rank columns."""
    age = clinical.loc[idx, "CDR_AGE"]
    e2f = scores.loc[idx, E2F]
    if legacy:
        raise ValueError("legacy=True is not supported; use the archived code and inputs for historical reproduction.")
    cov = pd.DataFrame({"age": age.astype(float), "E2F": e2f.astype(float)}, index=idx)
    if clinical.loc[idx, "PURITY"].notna().sum() >= MIN_PURITY_N:
        purity = clinical.loc[idx, "PURITY"]
        cov["purity"] = purity.astype(float)
    if gene in cnamat:
        target = cnamat.loc[idx, gene]
        cov["target_CNA"] = target.astype(float)
    if cohort == "BRCA":
        dummies = pd.get_dummies(clinical.loc[idx, "SUBTYPE"].fillna("Unclassified"),
                                 drop_first=True, dtype=float)
        cov = cov.join(dummies)
    rank_columns = tuple(c for c in ("age", "E2F", "purity", "target_CNA") if c in cov)
    return cov, rank_columns


def pathway_associations(expr: pd.DataFrame, zrna: pd.DataFrame, scores: pd.DataFrame,
                         clinical: pd.DataFrame, cnamat: pd.DataFrame, sets: dict,
                         alias: dict, legacy: bool = False) -> pd.DataFrame:
    """Raw and adjusted rank associations between each gene and each programme score.

    The candidate gene is removed from its programme score to avoid circularity.
    """
    rows = []
    for cohort, mask in cohort_masks(clinical):
        idx = clinical.index[mask]
        for g in expr:
            cov, rank_columns = covariate_frame(clinical, scores, cnamat, idx, cohort, g, legacy)
            for p in PATHS:
                members = [m for m in sets[p] if m in zrna.index and m != alias[g]]
                y = zrna.loc[members, idx].mean(axis=0)
                x = expr.loc[idx, g]
                raw = pd.concat([x, y], axis=1).dropna()
                r, pv = stats.spearmanr(raw.iloc[:, 0], raw.iloc[:, 1])
                n, pr, pp = partial_rank_correlation(x, y, cov, rank_columns)
                rows.append({"cohort": cohort, "gene": g, "pathway": p,
                             "signature_gene_n": len(members), "raw_n": len(raw),
                             "raw_spearman_rho": r, "raw_p_value": pv, "adjusted_n": n,
                             "partial_rank_rho": pr, "adjusted_p_value": pp,
                             "covariates": ";".join(cov.columns),
                             "target_excluded_from_score": True})
    out = pd.DataFrame(rows)
    out["raw_BH_FDR"] = out.groupby("cohort").raw_p_value.transform(benjamini_hochberg)
    out["adjusted_BH_FDR"] = out.groupby("cohort").adjusted_p_value.transform(benjamini_hochberg)
    return out


def rppa_associations(raw: Path, expr: pd.DataFrame, clinical: pd.DataFrame) -> pd.DataFrame:
    """Unadjusted Spearman RNA-RPPA correlations for the selected readouts."""
    rppa = (pd.read_csv(raw / "data_rppa.txt", sep="\t").set_index("Composite.Element.REF")
            .apply(pd.to_numeric, errors="coerce"))
    rows = []
    for cohort, mask in cohort_masks(clinical):
        idx = clinical.index[mask]
        for g in expr:
            for rp in RPPA_READOUTS:
                if rp not in rppa.index:
                    continue
                d = pd.concat([expr.loc[idx, g], rppa.loc[rp].reindex(idx)], axis=1).dropna()
                if len(d) < MIN_RPPA_N:
                    continue
                r, p = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])
                rows.append({"cohort": cohort, "gene": g, "RPPA_readout": rp, "n": len(d),
                             "spearman_rho": r, "p_value": p})
    out = pd.DataFrame(rows)
    out["BH_FDR"] = out.groupby("cohort").p_value.transform(benjamini_hochberg)
    return out


def fit_cox(d: pd.DataFrame) -> dict:
    """Cox model of ``gene_z`` (plus covariates present) with a PH check.

    Warnings raised by lifelines (convergence, collinearity) are kept in
    ``fit_warnings`` instead of being silenced.
    """
    from lifelines import CoxPHFitter
    from lifelines.statistics import proportional_hazard_test
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fit = CoxPHFitter().fit(d, duration_col="time_days", event_col="event")
        s = fit.summary.loc["gene_z"]
        ph = proportional_hazard_test(fit, d, time_transform="rank").summary.loc["gene_z", "p"]
    messages = sorted({str(w.message).split("\n")[0][:120] for w in caught})
    return {"HR_per_expression_SD": float(s["exp(coef)"]),
            "CI95_lower": float(s["exp(coef) lower 95%"]),
            "CI95_upper": float(s["exp(coef) upper 95%"]),
            "p_value": float(s.p), "PH_rank_p_value": float(ph), "status": "fitted",
            "fit_warnings": " | ".join(messages)}


def survival_models(expr: pd.DataFrame, clinical: pd.DataFrame) -> pd.DataFrame:
    """Continuous-expression Cox models per cohort, endpoint, gene and model."""
    rows = []
    for cohort, mask in cohort_masks(clinical):
        idx = clinical.index[mask]
        for ep in ENDPOINTS:
            for g in expr:
                e = expr.loc[idx, g]
                base = pd.DataFrame({"time_days": clinical.loc[idx, ep + ".time"],
                                     "event": clinical.loc[idx, ep],
                                     "gene_z": (e - e.mean()) / e.std(ddof=1),
                                     "age_per10": clinical.loc[idx, "CDR_AGE"] / 10,
                                     "stage_III_IV": clinical.loc[idx, "STAGE_ADVANCED"]})
                for model, cols, min_events in COX_MODELS:
                    d = base[["time_days", "event"] + cols].dropna()
                    d = d.loc[d.time_days > 0]
                    events = int(d.event.sum())
                    row = {"cohort": cohort, "endpoint": ep, "gene": g, "model": model,
                           "n": len(d), "events": events, "HR_per_expression_SD": np.nan,
                           "CI95_lower": np.nan, "CI95_upper": np.nan, "p_value": np.nan,
                           "PH_rank_p_value": np.nan,
                           "status": "insufficient_events" if events < min_events else "pending",
                           "fit_warnings": ""}
                    if events >= min_events:
                        try:
                            row.update(fit_cox(d))
                        except Exception as error:  # noqa: BLE001 - recorded per row
                            row["status"] = "fit_failed: " + str(error)[:130]
                    rows.append(row)
    out = pd.DataFrame(rows)
    fit_warnings = out.pop("fit_warnings")
    group = out.groupby(["cohort", "endpoint", "model"])
    out["BH_FDR"] = group.p_value.transform(benjamini_hochberg)
    out["PH_BH_FDR"] = group.PH_rank_p_value.transform(benjamini_hochberg)
    out["fit_warnings"] = fit_warnings  # appended last to keep the archived column order
    return out


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
def make_figures(fig_dir: Path, expr, contrasts, cnarows, cnamat, assoc, surv, rprows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    fig_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.facecolor": "white"})
    ko_colour, act_colour = "#A44C45", "#527C91"

    def save(f, name):
        f.savefig(fig_dir / (name + ".png"), dpi=300, bbox_inches="tight")
        plt.close(f)

    # F1 expression contrasts under independent clinical and intrinsic definitions.
    f, ax = plt.subplots(1, 2, figsize=(11, 8.5), sharey=True)
    titles = ["Clinical TNBC versus receptor-positive tumours",
              "PAM50 basal versus other assigned subtypes"]
    for a, (comp, *_), title in zip(ax, CONTRASTS, titles, strict=True):
        d = contrasts.loc[contrasts.comparison.eq(comp)].set_index("gene").reindex(GENES)
        a.axvline(0, color=".5", lw=.8)
        for i, g in enumerate(GENES):
            if g not in expr:
                a.text(0, i, "NA", color=".5", ha="center", va="center")
                continue
            v, q = d.loc[g, "median_delta_log2_RSEM_plus1"], d.loc[g, "BH_FDR"]
            a.barh(i, v, color=ko_colour if g in KO else act_colour, height=.7)
            a.text(v + (.04 if v >= 0 else -.04), i, "*" if q < .05 else "",
                   ha="left" if v >= 0 else "right", va="center")
        a.set_title(title, fontsize=11)
        a.set_xlabel("Difference in median log₂(RSEM + 1)")
        a.grid(axis="x", alpha=.15)
    ax[0].set_yticks(range(len(GENES)), GENES)
    ax[0].invert_yaxis()
    f.suptitle("TCGA-BRCA: candidate expression context", fontsize=15)
    f.text(.02, .015, "* BH FDR < 0.05 across the available candidate panel within each "
           "contrast. Missing LINC03007 is not plotted.\nDifferences are transformed-expression "
           "contrasts, not log₂ fold changes; tumour associations do not establish drug response.",
           fontsize=9)
    f.tight_layout(rect=[0, .07, 1, .95])
    save(f, "TCGA_F01_expression_contrasts")

    # F2 CNA prevalence.
    f, axs = plt.subplots(1, 2, figsize=(12, 8.5), sharey=True)
    bands = [("deep_deletion_n", "Deep deletion", "#35628A"),
             ("shallow_loss_n", "Shallow loss", "#A3BED5"), ("diploid_n", "Diploid", "#E7E7E7"),
             ("gain_n", "Gain", "#E8B4A5"), ("amplification_n", "Amplification", "#A7473B")]
    for a, cohort in zip(axs, ["BRCA", "TNBC"], strict=True):
        d = cnarows.loc[cnarows.cohort.eq(cohort)].set_index("gene").reindex(GENES)
        left = np.zeros(len(GENES))
        for field, label, colour in bands:
            v = (100 * d[field] / d.n_profiled).fillna(0).to_numpy()
            a.barh(range(len(GENES)), v, left=left, color=colour, label=label, height=.7)
            left += v
        a.set_xlim(0, 100)
        a.set_title(cohort)
        a.set_xlabel("Percentage of CNA-profiled tumours")
        for i, g in enumerate(GENES):
            if g not in cnamat:
                a.text(50, i, "NA", color=".5", ha="center", va="center")
    axs[0].set_yticks(range(len(GENES)), GENES)
    axs[0].invert_yaxis()
    axs[1].legend(loc="upper left", bbox_to_anchor=(1, 1))
    f.suptitle("TCGA-BRCA: copy-number context", fontsize=15)
    f.text(.03, .015, "GISTIC −2/−1/0/+1/+2. Shallow loss does not establish biallelic "
           "knockout; amplification does not establish CRISPRa-equivalent activity.", fontsize=9)
    f.tight_layout(rect=[0, .06, 1, .95])
    save(f, "TCGA_F02_CNA_prevalence")

    # F3 adjusted within-TNBC co-expression heatmap.
    f, a = plt.subplots(figsize=(13, 8.5))
    d = assoc.loc[assoc.cohort.eq("TNBC")]
    mat = (d.pivot(index="gene", columns="pathway", values="partial_rank_rho")
           .reindex(index=GENES, columns=PATHS))
    q = (d.pivot(index="gene", columns="pathway", values="adjusted_BH_FDR")
         .reindex(index=GENES, columns=PATHS))
    im = a.imshow(mat, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-.75, vcenter=0, vmax=.75))
    a.set_yticks(range(len(GENES)), GENES)
    a.set_xticks(range(len(PATHS)), [p.replace("_", " ").title() for p in PATHS], rotation=45, ha="right")
    for i in range(len(GENES)):
        for j in range(len(PATHS)):
            if pd.notna(mat.iloc[i, j]):
                a.text(j, i, "*" if q.iloc[i, j] < .05 else "", ha="center", va="center", fontsize=13)
    f.colorbar(im, ax=a, label="Partial rank correlation")
    a.set_title("Clinical TNBC: candidate–pathway co-expression after covariate adjustment", pad=16)
    f.text(.02, .005, "Adjusted for age, E2F proliferation score, purity where covered, and "
           "target CNA. Candidate excluded from its pathway score.\n* BH FDR < 0.05 across all "
           "gene–pathway tests in TNBC. These observational correlations are not synthetic-lethal "
           "dependencies.", fontsize=9)
    f.tight_layout(rect=[0, .07, 1, 1])
    save(f, "TCGA_F03_TNBC_pathway_associations")

    # F4 PFI forest plots, adjusted where the event count allows.
    f, axs = plt.subplots(1, 2, figsize=(12, 9), sharey=True)
    for a, cohort in zip(axs, ["BRCA", "TNBC"], strict=True):
        d = surv.loc[surv.cohort.eq(cohort) & surv.endpoint.eq("PFI")]
        chosen = []
        for i, g in enumerate(GENES):
            sub = d.loc[d.gene.eq(g) & d.status.eq("fitted")]
            s = sub.loc[sub.model.eq("age_stage_adjusted")]
            if s.empty:
                s = sub.loc[sub.model.eq("univariable")]
            if s.empty:
                continue
            s = s.iloc[0]
            a.errorbar(s.HR_per_expression_SD, i,
                       xerr=[[s.HR_per_expression_SD - s.CI95_lower],
                             [s.CI95_upper - s.HR_per_expression_SD]],
                       fmt="o", color=ko_colour if g in KO else act_colour, capsize=2)
            a.text(3.65, i, "*" if s.BH_FDR < .05 else "", va="center")
            chosen.append(s)
        a.axvline(1, color=".5", lw=.8)
        a.set_xscale("log")
        a.set_xlim(.35, 4)
        a.set_xlabel("PFI hazard ratio per expression SD (95% CI)")
        adjusted = any(x.model == "age_stage_adjusted" for x in chosen)
        a.set_title(cohort + (" — age/stage adjusted" if adjusted else " — univariable"))
    axs[0].set_yticks(range(len(GENES)), GENES)
    axs[0].invert_yaxis()
    f.suptitle("TCGA-CDR: continuous-expression prognosis", fontsize=15)
    f.text(.02, .015, "* BH FDR < 0.05 within cohort/endpoint/model. No optimized cut points. "
           "PFI is not an olaparib-specific outcome.\nModel sample/event counts and "
           "proportional-hazards checks are provided in the supplementary table.", fontsize=9)
    f.tight_layout(rect=[0, .07, 1, .95])
    save(f, "TCGA_F04_PFI_forest")

    # F5 selected RPPA correlations within TNBC.
    f, a = plt.subplots(figsize=(12, 8.5))
    d = rprows.loc[rprows.cohort.eq("TNBC")]
    mat = (d.pivot(index="gene", columns="RPPA_readout", values="spearman_rho")
           .reindex(index=GENES, columns=RPPA_READOUTS))
    q = (d.pivot(index="gene", columns="RPPA_readout", values="BH_FDR")
         .reindex(index=GENES, columns=RPPA_READOUTS))
    im = a.imshow(mat, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-.6, vcenter=0, vmax=.6))
    a.set_yticks(range(len(GENES)), GENES)
    a.set_xticks(range(len(RPPA_READOUTS)), [x.split("|")[1] for x in RPPA_READOUTS], rotation=45, ha="right")
    for i in range(len(GENES)):
        for j in range(len(RPPA_READOUTS)):
            if pd.notna(mat.iloc[i, j]) and q.iloc[i, j] < .05:
                a.text(j, i, "*", ha="center", va="center", fontsize=13)
    f.colorbar(im, ax=a, label="Spearman correlation")
    a.set_title("Clinical TNBC: RNA associations with RPPA protein/phosphoprotein readouts", pad=16)
    f.text(.02, .015, "* BH FDR < 0.05 across available gene–readout pairs within TNBC. "
           "Correlations are unadjusted and exploratory.\nPhosphoprotein abundance is not "
           "normalized to total protein here; no intervention or drug response is measured.",
           fontsize=9)
    f.tight_layout(rect=[0, .08, 1, 1])
    save(f, "TCGA_F05_TNBC_RPPA_associations")


# --------------------------------------------------------------------------- #
# Pipeline
# --------------------------------------------------------------------------- #
def run(root: Path, raw: Path | None = None, fig_dir: Path | None = None,
        legacy: bool = False, figures: bool = True) -> dict:
    """Execute the full analysis; write tables into ``root`` and return them."""
    if legacy:
        raise ValueError("legacy=True is not supported; use the archived code and inputs for historical reproduction.")
    root = Path(root)
    raw = Path(raw) if raw else root / "raw"
    fig_dir = Path(fig_dir) if fig_dir else root / "figures"
    root.mkdir(parents=True, exist_ok=True)

    rna_raw = load_rna(raw)
    clinical, n_primary = build_clinical_cohort(raw, rna_sample_ids=rna_raw.columns)
    qc_metadata = clinical.attrs["qc_metadata"]
    clinical = attach_purity(clinical, load_purity(raw))
    alias, maps = gene_coverage(rna_raw)
    sets = load_hallmark(raw / HALLMARK_GMT)
    required_symbols = {g for p in SCORED_PATHWAYS for g in sets[p]} | {s for s in alias.values() if s is not None}
    clinical, logrna, expr = prepare_expression(rna_raw, clinical, alias, maps,
                                               required_symbols=required_symbols)
    resolution = pd.DataFrame(logrna.attrs["symbol_resolution"], columns=[
        "gene", "source_entrez_ids", "selected_entrez_id", "resolution"])
    resolution.to_csv(root / "RNA_symbol_resolution.csv", index=False)
    log.info("Matched %d primary patients; %d of %d panel genes with expression",
             len(clinical), expr.shape[1], len(GENES))
    expr.to_csv(root / "selected_expression_log2_RSEM_plus1.csv")

    zrna, scores, path_coverage = hallmark_scores(logrna, sets)
    scores.to_csv(root / "hallmark_mean_z_scores.csv")
    path_coverage.to_csv(root / "hallmark_coverage.csv", index=False)
    clinical.to_csv(root / "matched_clinical_cohort.csv", index_label="SAMPLE_ID")

    counts = cohort_counts(clinical, n_primary, maps, qc_metadata)
    (root / "cohort_counts.json").write_text(json.dumps(counts, indent=2))
    print("COHORT_COUNTS", json.dumps(counts), flush=True)

    contrasts = expression_contrasts(expr, clinical)
    contrasts.to_csv(root / "expression_contrasts.csv", index=False)

    cna_raw = load_cna(raw)
    cna_alias = resolve_cna_aliases(cna_raw, maps)
    cnamat, cnarows = cna_analysis(cna_raw, cna_alias, clinical, expr)
    cnamat.to_csv(root / "selected_CNA_GISTIC.csv", index_label="SAMPLE_ID")
    cnarows.to_csv(root / "CNA_prevalence_and_expression.csv", index=False)
    coverage = pd.DataFrame(maps)
    coverage["CNA_source_symbol"] = coverage.gene.map(cna_alias)
    coverage["CNA_coverage"] = coverage.gene.map(lambda g: "available" if g in cnamat else "unavailable")
    coverage.to_csv(root / "gene_coverage.csv", index=False)

    assoc = pathway_associations(expr, zrna, scores, clinical, cnamat, sets, alias, legacy)
    assoc.to_csv(root / "pathway_associations.csv", index=False)
    rprows = rppa_associations(raw, expr, clinical)
    rprows.to_csv(root / "RPPA_associations.csv", index=False)
    surv = survival_models(expr, clinical)
    surv.to_csv(root / "survival_continuous_Cox.csv", index=False)

    if figures:
        make_figures(fig_dir, expr, contrasts, cnarows, cnamat, assoc, surv, rprows)

    tnbc_assoc = assoc.cohort.eq("TNBC") & assoc.adjusted_BH_FDR.lt(.05)
    tnbc_rppa = rprows.cohort.eq("TNBC") & rprows.BH_FDR.lt(.05)
    summary = {
        "expression_significant": contrasts.loc[contrasts.BH_FDR < .05].to_dict("records"),
        "PFI_significant": surv.loc[surv.endpoint.eq("PFI") & surv.BH_FDR.lt(.05)].to_dict("records"),
        "TNBC_pathway_significant": assoc.loc[tnbc_assoc].to_dict("records"),
        "TNBC_RPPA_significant": rprows.loc[tnbc_rppa].to_dict("records"),
    }
    (root / "analysis_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print("FINISHED", json.dumps({k.replace("_significant", "_sig"): len(v)
                                  for k, v in summary.items()}), flush=True)
    return {"clinical": clinical, "expr": expr, "scores": scores, "contrasts": contrasts,
            "cnamat": cnamat, "cnarows": cnarows, "coverage": coverage, "assoc": assoc,
            "rprows": rprows, "surv": surv, "counts": counts, "summary": summary}


def main(argv=None) -> int:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--root", type=Path, default=here,
                        help="Output directory for tables (default: script directory).")
    parser.add_argument("--raw", type=Path, default=None,
                        help="Input snapshot directory (default: ROOT/raw).")
    parser.add_argument("--figures-dir", type=Path, default=None,
                        help="Figure output directory (default: ROOT/figures).")
    parser.add_argument("--no-figures", action="store_true", help="Skip figure rendering.")
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args(argv)
    logging.basicConfig(level=args.log_level, format="%(levelname)s %(message)s")
    run(args.root, args.raw, args.figures_dir, figures=not args.no_figures)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
