"""QC decisions must survive the functional refactor and its synthetic tests."""
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from analysis.tcga_extended import analyze_tcga as m
import synthetic_tcga


def test_import_does_not_read_analysis_inputs_or_create_outputs(tmp_path):
    module_path = Path(m.__file__).resolve().parent
    program = """
import pathlib
import sys
from unittest.mock import patch
import pandas as pd

sys.path.insert(0, sys.argv[1])
def unexpected(*args, **kwargs):
    raise AssertionError('Analysis I/O occurred during import')
with patch.object(pd, 'read_csv', unexpected), patch.object(pd, 'read_excel', unexpected), \\
     patch.object(pathlib.Path, 'mkdir', unexpected), patch.object(pathlib.Path, 'write_text', unexpected):
    import analyze_tcga
assert callable(analyze_tcga.run)
"""
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    result = subprocess.run([sys.executable, "-c", program, str(module_path)],
                            cwd=tmp_path, env=env, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    assert list(tmp_path.iterdir()) == []


def test_clinical_builder_applies_redaction_qc_before_expression(tmp_path):
    raw = tmp_path / "raw"
    synthetic_tcga.make_dataset(raw, n=60, seed=18)
    path = raw / "TCGA-CDR-SupplementalTableS1.xlsx"
    cdr = pd.read_excel(path, sheet_name="TCGA-CDR")
    cdr["Redaction"] = cdr.Redaction.astype(object)
    cdr.loc[0, "Redaction"] = "Redacted"
    cdr.loc[1, "Redaction"] = " redacted "
    excluded = set(cdr.loc[[0, 1], "bcr_patient_barcode"])
    cdr.to_excel(path, sheet_name="TCGA-CDR", index=False)
    clinical, n_before_qc = m.build_clinical_cohort(raw)
    assert n_before_qc == 60
    assert len(clinical) == 58
    assert not excluded.intersection(clinical.PATIENT_ID)
    assert clinical.PATIENT_ID.is_unique


@pytest.mark.parametrize("reverse_rows", [False, True])
def test_expression_preparation_uses_reviewed_qsox1_for_pathway_scores(reverse_rows):
    samples = ["TCGA-AA-0001-01", "TCGA-AA-0002-01", "TCGA-AA-0003-01"]
    rna = pd.DataFrame({
        "Hugo_Symbol": ["QSOX1", "QSOX1", "PARG", "OTHER"],
        "Entrez_Gene_Id": [200058, 5768, 8505, 7],
        samples[0]: [1023.0, 1.0, 3.0, 15.0],
        samples[1]: [7.0, 31.0, 15.0, 7.0],
        samples[2]: [0.0, 127.0, 63.0, 1.0],
    })
    if reverse_rows:
        rna = rna.iloc[::-1]
    clinical = pd.DataFrame({"SAMPLE_ID": samples})
    alias, maps = m.gene_coverage(rna, genes=["PARG"])
    _, logrna, expr = m.prepare_expression(rna, clinical, alias, maps)
    np.testing.assert_allclose(logrna.loc["QSOX1"], np.log2([2.0, 32.0, 128.0]))
    assert expr.columns.tolist() == ["PARG"]
    sets = {"GLYCOLYSIS": ["QSOX1", "OTHER"]}
    _, scores, _ = m.hallmark_scores(logrna, sets, pathways=["GLYCOLYSIS"])
    direct = np.log2(np.array([[2.0, 32.0, 128.0], [16.0, 8.0, 2.0]]))
    expected = ((direct - direct.mean(axis=1, keepdims=True)) /
                direct.std(axis=1, ddof=1, keepdims=True)).mean(axis=0)
    np.testing.assert_allclose(scores.GLYCOLYSIS, expected)
