"""Unit tests for the statistical and classification helpers in analyze_tcga."""
import numpy as np
import pandas as pd
import pytest
from scipy import stats

from analysis.tcga_extended import analyze_tcga as m


# --------------------------------------------------------------------------- #
# Benjamini-Hochberg
# --------------------------------------------------------------------------- #
def test_bh_matches_scipy_reference():
    rng = np.random.default_rng(1)
    p = rng.uniform(size=200)
    expected = stats.false_discovery_control(p, method="bh")
    assert np.allclose(m.benjamini_hochberg(p), expected)


def test_bh_ignores_nan_and_keeps_positions():
    p = np.array([0.01, np.nan, 0.04, 0.5, np.nan])
    out = m.benjamini_hochberg(p)
    assert np.isnan(out[1]) and np.isnan(out[4])
    finite = stats.false_discovery_control([0.01, 0.04, 0.5], method="bh")
    assert np.allclose(out[[0, 2, 3]], finite)


def test_bh_handles_empty_and_all_nan():
    assert m.benjamini_hochberg([]).shape == (0,)
    assert np.isnan(m.benjamini_hochberg([np.nan, np.nan])).all()


def test_bh_is_monotone_and_bounded():
    rng = np.random.default_rng(2)
    p = np.sort(rng.uniform(size=50))
    q = m.benjamini_hochberg(p)
    assert (np.diff(q) >= -1e-15).all() and q.max() <= 1 and (q >= p - 1e-12).all()


# --------------------------------------------------------------------------- #
# Partial rank correlation
# --------------------------------------------------------------------------- #
def _residual(v, w):
    a = np.column_stack([np.ones(len(v)), w])
    return v - a @ np.linalg.lstsq(a, v, rcond=None)[0]


def test_partial_rank_equals_rank_residual_correlation():
    rng = np.random.default_rng(3)
    n = 150
    c = rng.normal(size=n)
    x = 0.7 * c + rng.normal(size=n)
    y = 0.5 * x + 0.9 * c + rng.normal(size=n)
    cov = pd.DataFrame({"c": c})
    n_out, rho, p = m.partial_rank_correlation(pd.Series(x), pd.Series(y), cov, rank_columns=("c",))
    rx, ry, rc = (stats.rankdata(v) for v in (x, y, c))
    expected = np.corrcoef(_residual(rx, rc), _residual(ry, rc))[0, 1]
    assert n_out == n and rho == pytest.approx(expected) and 0 < p < 1
    # the partial correlation removes the shared confounder, so it is well below the raw one
    assert rho < stats.spearmanr(x, y)[0]


def test_partial_rank_without_covariates_is_spearman():
    rng = np.random.default_rng(4)
    x, y = rng.normal(size=80), rng.normal(size=80)
    _, rho, p = m.partial_rank_correlation(pd.Series(x), pd.Series(y), pd.DataFrame(index=range(80)))
    ref_rho, ref_p = stats.spearmanr(x, y)
    assert rho == pytest.approx(ref_rho) and p == pytest.approx(ref_p, rel=1e-6)


def test_partial_rank_returns_nan_below_minimum_n():
    x = pd.Series(np.arange(10.0))
    n, rho, p = m.partial_rank_correlation(x, x, pd.DataFrame(index=x.index))
    assert n == 10 and np.isnan(rho) and np.isnan(p)


@pytest.mark.parametrize("degeneracy", ["constant", "explained", "exhausted_df"])
def test_partial_rank_returns_nan_when_no_residual_variation_or_degrees_of_freedom(degeneracy):
    n = 40
    trend = np.arange(n, dtype=float)
    other = np.sin(trend)
    if degeneracy == "constant":
        fixed = np.ones(n)
        cov = pd.DataFrame(index=range(n))
        rank_columns = ()
    elif degeneracy == "explained":
        fixed = trend
        cov = pd.DataFrame({"c": trend})
        rank_columns = ("c",)
    else:
        fixed = trend
        # Intercept plus these indicator columns saturates the design matrix.
        cov = pd.DataFrame(np.eye(n)[:, 1:])
        rank_columns = ()
    for x, y in [(fixed, other), (other, fixed)]:
        observed_n, rho, p = m.partial_rank_correlation(
            pd.Series(x), pd.Series(y), cov, rank_columns=rank_columns)
        assert observed_n == n
        assert np.isnan(rho) and np.isnan(p)


def test_partial_rank_ranks_after_dropping_missing_rows():
    """Ranks must be recomputed on the retained rows (the legacy run ranked first)."""
    rng = np.random.default_rng(5)
    n = 120
    c = rng.normal(size=n)
    x = c + rng.normal(size=n)
    y = pd.Series(c + rng.normal(size=n))
    y.iloc[:40] = np.nan
    cov = pd.DataFrame({"c": c})
    _, rho_new, _ = m.partial_rank_correlation(pd.Series(x), y, cov, rank_columns=("c",))
    keep = y.notna().to_numpy()
    rx, ry, rc = (stats.rankdata(v[keep]) for v in (x, y.to_numpy(), c))
    assert rho_new == pytest.approx(np.corrcoef(_residual(rx, rc), _residual(ry, rc))[0, 1])
    # pre-ranked covariate reproduces the legacy value, which is different
    _, rho_legacy, _ = m.partial_rank_correlation(pd.Series(x), y, pd.DataFrame({"c": stats.rankdata(c)}))
    assert rho_legacy != pytest.approx(rho_new, abs=1e-9)


# --------------------------------------------------------------------------- #
# Clinical classification
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("ihc,fish,expected", [
    ("Positive", "Negative", "Discordant"),
    ("Negative", "Positive", "Discordant"),
    ("Positive", "Positive", "Positive"),
    ("Equivocal", "Negative", "Negative"),
    ("Positive", None, "Positive"),
    ("Positive", float("nan"), "Positive"),
    ("Equivocal", "Indeterminate", "Unresolved"),
    (None, None, "Unresolved"),
])
def test_classify_her2(ihc, fish, expected):
    assert m.classify_her2(ihc, fish) == expected


@pytest.mark.parametrize("er,pr,her2,expected", [
    ("Negative", "Negative", "Negative", "TNBC"),
    ("Positive", "Negative", "Negative", "non-TNBC"),
    ("Negative", "Negative", "Positive", "non-TNBC"),
    ("Indeterminate", "Negative", "Positive", "non-TNBC"),
    ("Negative", "Negative", "Discordant", "Unclassified"),
    ("Negative", "Indeterminate", "Negative", "Unclassified"),
    (None, "Negative", "Negative", "Unclassified"),
])
def test_classify_tnbc(er, pr, her2, expected):
    assert m.classify_tnbc(er, pr, her2) == expected


@pytest.mark.parametrize("subtype,expected", [
    ("BRCA_Basal", "Basal"), ("BRCA_LumA", "Other PAM50"), ("BRCA_Normal", "Other PAM50"),
    (None, "Unclassified"), (float("nan"), "Unclassified"),
])
def test_classify_basal(subtype, expected):
    assert m.classify_basal(subtype) == expected


@pytest.mark.parametrize("label,expected", [
    ("Stage IV", 1.0), ("Stage IIIA", 1.0), ("Stage IIIC", 1.0), ("stage iii", 1.0),
    ("Stage I", 0.0), ("Stage IA", 0.0), ("Stage IB", 0.0), ("Stage II", 0.0), ("Stage IIB", 0.0),
    ("Stage X", np.nan), ("[Discrepancy]", np.nan), ("[Not Available]", np.nan), (None, np.nan),
    (float("nan"), np.nan), ("", np.nan),
])
def test_stage_advanced(label, expected):
    out = m.stage_advanced(label)
    assert (np.isnan(out) and np.isnan(expected)) or out == expected


# --------------------------------------------------------------------------- #
# Symbol resolution and coverage
# --------------------------------------------------------------------------- #
def test_resolve_symbol_prefers_direct_then_alias():
    assert m.resolve_symbol("LAMTOR2", {"LAMTOR2", "ROBLD3"}) == "LAMTOR2"
    assert m.resolve_symbol("LAMTOR2", {"ROBLD3"}) == "ROBLD3"
    assert m.resolve_symbol("LINC03007", {"LOC646588"}) == "LOC646588"
    assert m.resolve_symbol("LINC03007", {"ABCB1"}) is None
    assert m.resolve_symbol("WDR24", {"ROBLD3"}) is None


def test_gene_coverage_flags_absent_and_duplicated_symbols():
    rna = pd.DataFrame({"Hugo_Symbol": ["ABCB1", "ROBLD3", "WT1", "WT1"],
                        "Entrez_Gene_Id": [5243, 28956, 7490, 7490], "S1": [1, 2, 3, 4]})
    alias, maps = m.gene_coverage(rna, genes=["ABCB1", "LAMTOR2", "WT1", "LINC03007"])
    by_gene = {x["gene"]: x for x in maps}
    assert alias == {"ABCB1": "ABCB1", "LAMTOR2": "ROBLD3", "WT1": "WT1", "LINC03007": None}
    assert by_gene["ABCB1"]["rna_coverage"] == "available"
    assert by_gene["LAMTOR2"]["entrez_gene_id"] == "28956"
    assert by_gene["WT1"] == {"gene": "WT1", "source_symbol": "WT1", "rna_rows": 2,
                              "entrez_gene_id": "", "rna_coverage": "unavailable_or_ambiguous"}
    assert by_gene["LINC03007"]["rna_rows"] == 0


def test_entrez_int_accepts_integral_strings_and_floats():
    assert m._entrez_int("28956") == 28956
    assert m._entrez_int("28956.0") == 28956
    assert m._entrez_int(28956.0) == 28956
    assert m._entrez_int("") is None and m._entrez_int(None) is None
    assert m._entrez_int("28956.5") is None and m._entrez_int(float("nan")) is None


def test_resolve_cna_aliases_direct_then_unique_entrez():
    cna = pd.DataFrame({"Hugo_Symbol": ["ABCB1", "ROBLD3", "A", "B"],
                        "Entrez_Gene_Id": [5243, 28956, 7, 7], "S1": [0, 0, 0, 0]})
    maps = [{"gene": "ABCB1", "entrez_gene_id": "5243"}, {"gene": "LAMTOR2", "entrez_gene_id": "28956"},
            {"gene": "WT1", "entrez_gene_id": "7"}, {"gene": "LINC03007", "entrez_gene_id": ""}]
    assert m.resolve_cna_aliases(cna, maps) == {"ABCB1": "ABCB1", "LAMTOR2": "ROBLD3",
                                                "WT1": None, "LINC03007": None}


# --------------------------------------------------------------------------- #
# Hallmark scores
# --------------------------------------------------------------------------- #
def test_hallmark_scores_are_mean_gene_z_with_constant_genes_dropped():
    logrna = pd.DataFrame({"S1": [1.0, 2.0, 5.0, 3.0], "S2": [3.0, 2.0, 5.0, 1.0], "S3": [2.0, 2.0, 5.0, 2.0]},
                          index=["G1", "G2", "G3", "G4"])
    sets = {"P": ["G1", "G3", "G4", "MISSING"], "E2F_TARGETS": ["G2"]}
    zrna, scores, coverage = m.hallmark_scores(logrna, sets, pathways=["P", "E2F_TARGETS"])
    assert np.isnan(zrna.loc["G2"]).all() and np.isnan(zrna.loc["G3"]).all()  # zero-variance rows
    expected = (zrna.loc["G1"] + zrna.loc["G4"]) / 2  # G3 is NaN and skipped by the mean
    assert np.allclose(scores["P"], expected)
    assert np.isnan(scores["E2F_TARGETS"]).all()
    assert coverage.set_index("pathway").loc["P", "rna_matched_genes"] == 3
    assert coverage.set_index("pathway").loc["P", "source_genes"] == 4


# --------------------------------------------------------------------------- #
# Cox wrapper
# --------------------------------------------------------------------------- #
def test_fit_cox_recovers_planted_hazard_and_records_warnings_field():
    rng = np.random.default_rng(6)
    n = 400
    z = rng.normal(size=n)
    time = rng.exponential(1 / (0.01 * np.exp(0.5 * z)))
    censor = rng.uniform(20, 200, n)
    d = pd.DataFrame({"time_days": np.minimum(time, censor), "event": (time <= censor).astype(float),
                      "gene_z": z})
    out = m.fit_cox(d)
    assert out["status"] == "fitted"
    assert 1.3 < out["HR_per_expression_SD"] < 2.1
    assert out["CI95_lower"] < out["HR_per_expression_SD"] < out["CI95_upper"]
    assert out["p_value"] < 1e-6 and 0 <= out["PH_rank_p_value"] <= 1
    assert isinstance(out["fit_warnings"], str)


def test_fit_cox_keeps_lifelines_warning_text():
    rng = np.random.default_rng(7)
    n = 60
    d = pd.DataFrame({"time_days": rng.exponential(100, n), "event": rng.integers(0, 2, n).astype(float),
                      "gene_z": rng.normal(size=n)})
    d["age_per10"] = 5 + rng.normal(size=n) * 1e-7  # near-constant covariate
    out = m.fit_cox(d)
    assert out["status"] == "fitted"
    assert "very low variance" in out["fit_warnings"]


def test_survival_models_record_fit_failure_instead_of_raising():
    rng = np.random.default_rng(8)
    n = 80
    index = pd.Index([f"S{i}" for i in range(n)], name="SAMPLE_ID")
    z = rng.normal(size=n)
    clinical = pd.DataFrame({"TNBC_STATUS": "non-TNBC", "BASAL_STATUS": "Other PAM50",
                             "CDR_AGE": 20 * z, "STAGE_ADVANCED": 0.0}, index=index)
    for ep in m.ENDPOINTS:
        clinical[ep + ".time"] = rng.exponential(300, n)
        clinical[ep] = rng.integers(0, 2, n).astype(float)
    expr = pd.DataFrame({"WT1": z}, index=index)
    surv = m.survival_models(expr, clinical)
    adjusted = surv.loc[surv.model.eq("age_stage_adjusted") & surv.cohort.eq("BRCA")]
    assert adjusted.status.str.startswith("fit_failed").all()  # age is collinear with gene_z
    assert surv.loc[surv.model.eq("univariable") & surv.cohort.eq("BRCA"), "status"].eq("fitted").all()
    assert surv.loc[surv.cohort.ne("BRCA"), "status"].eq("insufficient_events").all()
