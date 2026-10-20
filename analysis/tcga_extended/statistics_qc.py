"""Explicit RNA identifier resolution and complete-case TCGA statistics."""
import numpy as np
import pandas as pd
from scipy import stats

REDACTION_POLICY = ('Conservative TCGA-CDR flag-based QC: exclude Redaction=Redacted before all integrated analyses. '
                    'Historical flags can include corrected aliquots whose replacement identity cannot be established '
                    'from this sample-level matrix; flags do not imply all cases are intrinsically invalid or consent-withdrawn. Aggregate counts only.')


def exclude_cdr_redactions(clinical):
    """Apply the declared conservative CDR flag policy, returning data and mask.

    A historical flag need not mean every aliquot is invalid. This policy avoids
    assuming that a sample-level matrix contains the documented replacement.
    The caller reports aggregate counts rather than excluded case identifiers.
    """
    if 'Redaction' not in clinical:
        raise ValueError('TCGA-CDR Redaction is required for the declared QC policy.')
    excluded = clinical.Redaction.fillna('').astype(str).str.strip().str.casefold().eq('redacted')
    return clinical.loc[~excluded].copy(), excluded


def resolve_required_rna_symbols(rna):
    """Resolve ambiguous required symbols by an explicitly reviewed stable ID.

    Input must already be restricted to candidate/signature symbols. Historical
    Entrez 200058 (FLJ23867) was merged into protein-coding QSOX1/5768 in 2018:
    https://www.ncbi.nlm.nih.gov/gene/200058 . The two source rows have different
    measurements, so silently keeping whichever occurs first is inappropriate.
    """
    resolved = rna.copy()
    canonical_ids = {'QSOX1': 5768}
    audit = []
    duplicate_rows = resolved.loc[resolved.Hugo_Symbol.duplicated(keep=False)]
    for symbol, group in duplicate_rows.groupby('Hugo_Symbol', sort=True):
        canonical = canonical_ids.get(symbol)
        selected = (group.loc[pd.to_numeric(group.Entrez_Gene_Id, errors='coerce').eq(canonical)]
                    if canonical is not None else group.iloc[0:0])
        if len(selected) != 1:
            raise ValueError(f'Ambiguous required RNA symbol {symbol}: resolve its stable gene identifiers explicitly before analysis.')
        audit.append({'gene': symbol, 'source_entrez_ids': ';'.join(group.Entrez_Gene_Id.astype(str)),
                      'selected_entrez_id': canonical,
                      'resolution': 'current protein-coding Entrez identifier; retired duplicate excluded'})
        resolved = resolved.drop(index=group.index.difference(selected.index))
    return resolved, pd.DataFrame(audit, columns=['gene', 'source_entrez_ids', 'selected_entrez_id', 'resolution'])


def partialrank(x, y, cov, rank_columns=None, min_n=30):
    """Complete-case partial rank correlation; undefined fits return NaN."""
    d = pd.concat([x.rename('x'), y.rename('y'), cov], axis=1).dropna()
    n = len(d)
    if n < min_n:
        return n, np.nan, np.nan
    # None preserves the audited all-covariate convention. Named columns allow
    # continuous covariates to be ranked while binary subtype indicators remain
    # unchanged; the two binary encodings span the same adjustment space.
    design = d.iloc[:, 2:].astype(float).copy()
    columns = design.columns if rank_columns is None else rank_columns
    for column in columns:
        if column in design:
            design[column] = stats.rankdata(design[column])
    a = np.column_stack([np.ones(n), design.to_numpy(float)])
    k = np.linalg.matrix_rank(a) - 1
    df = n - k - 2
    if df <= 0:
        return n, np.nan, np.nan
    xr, yr = stats.rankdata(d.x), stats.rankdata(d.y)
    # Exact linear explanation can leave floating-point residual noise. Its
    # correlation is undefined, not evidence of a partial association.
    residuals = []
    for ranked in (xr, yr):
        residual = ranked - a @ np.linalg.lstsq(a, ranked, rcond=None)[0]
        tolerance = np.finfo(float).eps * max(a.shape) * max(np.linalg.norm(ranked), 1.0)
        if np.linalg.norm(residual - residual.mean()) <= tolerance:
            return n, np.nan, np.nan
        residuals.append(residual)
    xr, yr = residuals
    r = np.corrcoef(xr, yr)[0, 1]
    p = 2 * stats.t.sf(abs(r * np.sqrt(df / max(1 - r * r, 1e-15))), df)
    return n, r, p
