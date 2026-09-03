"""Smoke test: the pipeline runs end to end and its output has the expected shape.

Run:  .venv/bin/python -m pytest tests -q
"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssc_coh.clean import build_processed          # noqa: E402
from ssc_coh.config import TABLES                  # noqa: E402
from ssc_coh.features import build_features        # noqa: E402
from ssc_coh.raw import load_raw                   # noqa: E402

# frames that are not one row per measurement, so they carry no subject_id
NO_SUBJECT_ID = {"issues", "subject_id_map"}


@pytest.fixture(scope="module")
def processed():
    """One pipeline run shared by every test in this file."""
    frames = build_processed()
    frames["features"] = build_features(
        frames["subjects"], frames["ssc_subtype"], frames["vitals"],
        frames["labs"], frames["pft"], frames["mrss"], frames["medications"],
        frames["antibodies"], frames["bal"], frames["biopsies"],
        frames["libraries"],
    )
    return frames


def test_load_raw_returns_the_11_source_tables():
    raw_tables = load_raw()
    assert len(raw_tables) == 11
    assert sorted(raw_tables) == sorted(TABLES)
    assert all(len(table) for table in raw_tables.values())


def test_the_pipeline_returns_the_19_frames_written_to_data_processed(processed):
    assert len(processed) == 19          # 18 from build_processed, plus features
    assert "features" in processed
    assert "issues" in processed
    assert all(len(frame.columns) for frame in processed.values())


def test_issue_log_is_populated(processed):
    issues = processed["issues"]
    assert list(issues.columns) == ["table", "issue", "decision", "n_affected", "detail"]
    assert len(issues) > 0
    assert issues["table"].notna().all()
    assert issues["decision"].notna().all()


def test_every_measurement_frame_has_subject_id(processed):
    for name, frame in processed.items():
        if name in NO_SUBJECT_ID:
            continue
        assert "subject_id" in frame.columns, f"{name} has no subject_id column"
