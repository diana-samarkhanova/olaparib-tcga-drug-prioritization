"""Ensure context integration never silently uses the stale archived screen panel."""
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import pandas as pd


SOURCE = Path(__file__).resolve().parents[1] / 'analysis/astra_extended/integrate_context.py'


class ContextPanelSourceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.analysis = self.root / 'analysis/astra_extended'
        self.analysis.mkdir(parents=True)
        self.script = self.analysis / SOURCE.name
        shutil.copy2(SOURCE, self.script)
        pd.DataFrame({'gene': ['ABCB1'], 'scientific_priority': ['fixed curation']}).to_csv(
            self.analysis / 'final_gene_priority_matrix.csv', index=False)
        archive = self.root / 'tables'
        archive.mkdir()
        # The archival input deliberately disagrees with the generated panel.
        self.archived = pd.DataFrame({
            'gene': ['ABCB1', 'ABCB1'], 'timepoint': ['T14', 'T23'],
            'positive_FDR': [0.001, 0.001], 'log2FC': [2.0, 3.0],
        })
        self.archived.to_csv(archive / 'Selected_CRISPR_panel_20_genes.csv', index=False)

    def run_script(self):
        return subprocess.run([sys.executable, str(self.script)], capture_output=True, text=True)

    def test_generated_panel_wins_over_conflicting_archival_copy(self):
        target = self.root / 'outputs/TCGA_Drug_Analysis_EN_20261013/tables'
        target.mkdir(parents=True)
        current = self.archived.copy()
        current.loc[0, 'log2FC'] = -1.0
        current.to_csv(target / 'Selected_CRISPR_panel_20_genes.csv', index=False)
        result = self.run_script()
        self.assertEqual(result.returncode, 0, result.stderr)
        row = pd.read_csv(self.analysis / 'final_gene_priority_matrix.csv').iloc[0]
        self.assertEqual(row.screen_strict_positive_timepoints, 'T23')
        self.assertFalse(row.screen_all_displayed_timepoints_strict_positive)
        self.assertEqual(row.scientific_priority, 'fixed curation')

    def test_missing_generated_panel_does_not_fall_back_to_archive(self):
        result = self.run_script()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('plot_project_evidence.py', result.stderr)
        self.assertIn('not used as a fallback', result.stderr)


if __name__ == '__main__':
    unittest.main()
