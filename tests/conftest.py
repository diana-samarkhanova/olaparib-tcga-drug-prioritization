"""Shared fixtures: import path for the analysis modules and a synthetic TCGA snapshot."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
for module_dir in [REPO, REPO / "analysis" / "tcga_extended", REPO / "tests"]:
    if str(module_dir) not in sys.path:
        sys.path.insert(0, str(module_dir))

import synthetic_tcga  # noqa: E402


@pytest.fixture(scope="session")
def synthetic_raw(tmp_path_factory) -> tuple[Path, dict]:
    """Synthetic raw snapshot shared by the pipeline tests (built once per session)."""
    raw = tmp_path_factory.mktemp("synthetic") / "raw"
    truth = synthetic_tcga.make_dataset(raw, n=320, seed=0)
    return raw, truth


@pytest.fixture(scope="session")
def pipeline_outputs(synthetic_raw, tmp_path_factory):
    """Run the full analysis once (without figures) and expose root plus returned tables."""
    from analysis.tcga_extended import analyze_tcga

    raw, truth = synthetic_raw
    root = tmp_path_factory.mktemp("run")
    tables = analyze_tcga.run(root, raw, figures=False)
    return root, tables, truth
