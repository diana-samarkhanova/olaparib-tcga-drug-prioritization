"""Regression checks for the 2026-10-20 pharmacology audit."""
import ast
import importlib.util
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('prism_summary', ROOT / 'analysis/prism/summarize_curves.py')
summary = importlib.util.module_from_spec(spec)
spec.loader.exec_module(summary)


class DoseContextRegression(unittest.TestCase):
    def records(self):
        # The two archived selinexor/468 source records that exposed pooling.
        return pd.DataFrame({
            'ccle_name': ['MDAMB468_BREAST'] * 2,
            'screen_id': ['HTS002', 'MTS006'],
            'broad_id': ['BRD-K21361524-001-01-1'] * 2,
            'dose_min_uM': [10/16384, 2.5/16384],
            'dose_max_uM': [10.0, 2.5],
            'auc': [0.700781, 0.523278],
            'r2': [0.658758, 0.871472],
            'ic50': [0.509477, 0.020076],
            'ic50_nonfinite_or_missing': [False, False],
            'ic50_outside_tested_range': [False, False],
            'r2_nonpositive': [False, False],
            'auc_outside_0_1': [False, False],
        })

    def test_different_windows_do_not_produce_pooled_potency_or_auc(self):
        result = summary.summarize(self.records())
        self.assertFalse(result.pooled_dose_window_compatible)
        self.assertEqual(result.n_distinct_dose_windows, 2)
        for field in ['auc_median', 'r2_median', 'ic50_uM_median_within_tested_range']:
            self.assertTrue(np.isnan(result[field]), field)
        self.assertIn('MTS006', result.selected_curve_details_JSON)
        self.assertIn('HTS002', result.selected_curve_details_JSON)

    def test_each_context_retains_its_own_ic50(self):
        records = self.records()
        for index, expected in enumerate([0.509477, 0.020076]):
            result = summary.summarize(records.iloc[[index]])
            self.assertTrue(result.pooled_dose_window_compatible)
            self.assertAlmostEqual(result.ic50_uM_median_within_tested_range, expected)

    def test_source_rounding_does_not_create_false_window_difference(self):
        records = self.records()
        records['dose_max_uM'] = [10.0, 10.0]
        records['dose_min_uM'] = [0.00061035, 0.0006103515625]
        self.assertEqual(summary.dose_window_count(records), 1)

    def test_missing_window_cannot_be_pooled(self):
        records = self.records()
        records.loc[0, 'dose_max_uM'] = np.nan
        self.assertFalse(summary.summarize(records).pooled_dose_window_compatible)

    def test_selinexor_alias_survives_a_fresh_builder_run(self):
        # The archived table had this reviewed pair, but the builder omitted it.
        tree = ast.parse((ROOT / 'analysis/drugs_extended/curate_aliases.py').read_text())
        defs = next(ast.literal_eval(node.value) for node in tree.body
                    if isinstance(node, ast.Assign)
                    and any(isinstance(target, ast.Name) and target.id == 'defs'
                            for target in node.targets))
        self.assertIn(('B1464', 'selinexor', 'KPT-330'), [row[:3] for row in defs])


if __name__ == '__main__':
    unittest.main()
