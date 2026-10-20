"""Regression checks for scientific QC decisions, without clinical source data."""
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis' / 'tcga_extended'))
from statistics_qc import exclude_cdr_redactions, resolve_required_rna_symbols, partialrank


class TcgaQcTests(unittest.TestCase):
    def test_qsox1_choice_is_independent_of_source_row_order(self):
        source = pd.DataFrame({'Hugo_Symbol': ['QSOX1', 'QSOX1', 'PARG'],
                               'Entrez_Gene_Id': [200058, 5768, 8505],
                               'tumour': [900., 3600., 200.]})
        for frame in [source, source.iloc[::-1]]:
            selected, audit = resolve_required_rna_symbols(frame)
            qsox = selected.loc[selected.Hugo_Symbol.eq('QSOX1')]
            self.assertEqual(qsox.Entrez_Gene_Id.tolist(), [5768])
            self.assertEqual(qsox.tumour.tolist(), [3600.])
            self.assertEqual(audit.selected_entrez_id.tolist(), [5768])
            self.assertEqual(len(selected), 2)

    def test_unreviewed_ambiguity_fails_instead_of_choosing_first(self):
        source = pd.DataFrame({'Hugo_Symbol': ['UNKNOWN', 'UNKNOWN'],
                               'Entrez_Gene_Id': [1, 2], 'tumour': [10., 20.]})
        with self.assertRaisesRegex(ValueError, 'Ambiguous required RNA symbol UNKNOWN'):
            resolve_required_rna_symbols(source)

    def test_incomplete_observations_cannot_change_complete_case_ranks(self):
        rng = np.random.default_rng(171)
        n = 50
        cov = pd.DataFrame({'age': np.arange(n, dtype=float), 'purity': rng.uniform(size=n)})
        x = pd.Series(.1 * cov.age + rng.normal(size=n))
        y = pd.Series(.3 * cov.age + rng.normal(size=n))
        expected = partialrank(x, y, cov)
        # Interspersed covariates would distort ranks if ranked before dropna.
        extra = pd.DataFrame({'age': np.repeat(10.5, 25), 'purity': np.repeat(.25, 25)}, index=np.arange(n, n + 25))
        extended_cov = pd.concat([cov, extra])
        extended_x = x.reindex(extended_cov.index)
        extended_y = pd.concat([y, pd.Series(np.arange(25), index=extra.index)])
        observed = partialrank(extended_x, extended_y, extended_cov)
        self.assertEqual(observed[0], n)
        np.testing.assert_allclose(observed, expected, rtol=0, atol=1e-12)

    def test_cdr_policy_removes_only_explicit_flag_and_preserves_input(self):
        source = pd.DataFrame({'Redaction': [np.nan, 'Redacted', ' redacted ', ''],
                               'tumour_expression': [1., 2., 3., 4.]})
        eligible, excluded = exclude_cdr_redactions(source)
        self.assertEqual(excluded.tolist(), [False, True, True, False])
        self.assertEqual(eligible.tumour_expression.tolist(), [1., 4.])
        self.assertEqual(len(source), 4)
        with self.assertRaisesRegex(ValueError, 'Redaction is required'):
            exclude_cdr_redactions(source.drop(columns='Redaction'))


if __name__ == '__main__':
    unittest.main()
