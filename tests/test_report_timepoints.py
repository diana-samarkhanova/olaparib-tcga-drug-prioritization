"""Regression: shuffled input rows must not swap early and late estimates."""
import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis/report_extended"))
from report_helpers import select_screen_timepoint


class ScreenTimepointTests(unittest.TestCase):
    def setUp(self):
        self.panel = pd.DataFrame([
            {"gene": "A", "library": "KO", "timepoint": "T24", "log2FC": 24.},
            {"gene": "B", "library": "SAM", "timepoint": "T14", "log2FC": 14.},
            {"gene": "A", "library": "KO", "timepoint": "T16", "log2FC": 16.},
            {"gene": "B", "library": "SAM", "timepoint": "T23", "log2FC": 23.},
        ])

    def test_shuffled_timepoints_keep_their_biological_labels(self):
        early=select_screen_timepoint(self.panel,"early").set_index("gene")
        late=select_screen_timepoint(self.panel,"late").set_index("gene")
        self.assertEqual(early.log2FC.to_dict(), {"B":14.,"A":16.})
        self.assertEqual(late.log2FC.to_dict(), {"A":24.,"B":23.})

    def test_missing_comparison_does_not_silently_drop_a_gene(self):
        with self.assertRaisesRegex(ValueError,"Missing early"):
            select_screen_timepoint(self.panel.drop(index=2),"early")

    def test_duplicate_comparison_is_rejected(self):
        with self.assertRaisesRegex(ValueError,"Duplicate"):
            select_screen_timepoint(pd.concat([self.panel,self.panel.iloc[[2]]]),"early")


if __name__ == "__main__":
    unittest.main()
