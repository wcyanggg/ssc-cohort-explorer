"""Rule-level tests: the cleaning decisions and the shared statistics helpers.

Run:  .venv/bin/python -m pytest tests -q
"""
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssc_coh import clean                            # noqa: E402
from ssc_coh.stats import association_test, prevalence_summary   # noqa: E402


@pytest.fixture(autouse=True)
def fresh_issue_log():
    """Each test gets its own issue log and remap, so the rows it reads are its own."""
    clean._ISSUES = clean.IssueLog()
    clean._ID_REMAP = {}
    yield
    clean._ISSUES = clean.IssueLog()
    clean._ID_REMAP = {}


def issue_log() -> pd.DataFrame:
    return clean._ISSUES.frame()


# ------------------------------------------------------------ prevalence bounds
def test_prevalence_with_nothing_missing_puts_both_bounds_on_the_complete_case_rate():
    summary = prevalence_summary(pd.Series([True, True, False, False], dtype="boolean"))
    assert (summary.n_all, summary.n_recorded, summary.n_missing) == (4, 4, 0)
    assert summary.complete_case_rate == 0.5
    assert summary.lower_bound == 0.5
    assert summary.upper_bound == 0.5


def test_prevalence_with_everything_missing_has_no_complete_case_rate_and_full_bounds():
    summary = prevalence_summary(pd.Series([pd.NA, pd.NA, pd.NA], dtype="boolean"))
    assert (summary.n_all, summary.n_recorded, summary.n_missing) == (3, 0, 3)
    assert pd.isna(summary.complete_case_rate)
    assert summary.lower_bound == 0.0
    assert summary.upper_bound == 1.0


def test_prevalence_with_partial_missing_separates_the_three_rates():
    # 6 positives, 2 negatives, 2 unrecorded
    flag_values = pd.Series([True] * 6 + [False] * 2 + [pd.NA] * 2, dtype="boolean")
    summary = prevalence_summary(flag_values)
    assert (summary.n_all, summary.n_recorded, summary.n_missing) == (10, 8, 2)
    assert summary.complete_case_rate == 0.75
    assert summary.lower_bound == 0.6
    assert summary.upper_bound == 0.8
    assert summary.lower_bound <= summary.complete_case_rate <= summary.upper_bound


def test_prevalence_of_the_cohort_matches_the_rates_reported_in_the_app():
    features = pd.read_parquet(ROOT / "data" / "processed" / "features.parquet")
    summary = prevalence_summary(features["dx_ild"])
    assert (summary.n_all, summary.n_recorded, summary.n_missing) == (1498, 943, 555)
    assert round(summary.complete_case_rate * 100, 1) == 67.1
    assert round(summary.lower_bound * 100, 1) == 42.3
    assert round(summary.upper_bound * 100, 1) == 79.3


# ------------------------------------------------------- association test gate
def test_association_test_uses_chi_square_when_every_expected_count_is_large():
    result = association_test(pd.DataFrame([[50, 50], [40, 60]]))
    assert result.test_name == "chi-square"
    assert result.smallest_expected_count >= 5


def test_association_test_falls_back_to_fisher_when_an_expected_count_is_small():
    # observed cells look usable but one expected count is far below 5
    result = association_test(pd.DataFrame([[1, 1], [1, 20]]))
    assert result.test_name == "Fisher's exact"
    assert result.smallest_expected_count < 5
    assert 0 <= result.p_value <= 1


def test_association_test_reports_descriptives_only_for_a_small_table_larger_than_two_by_two():
    result = association_test(pd.DataFrame([[1, 1, 1], [1, 1, 20], [2, 2, 2]]))
    assert result.test_name == "none"
    assert "does not apply" in result.reason


def test_association_test_needs_two_groups_and_two_outcome_values():
    result = association_test(pd.DataFrame([[5, 5]]))
    assert result.test_name == "none"
    assert "only one group" in result.reason


# ------------------------------------------------------------ onset order flag
def subtype_frame(raynaud_dates, nonraynaud_dates) -> pd.DataFrame:
    return pd.DataFrame({
        "study_code": [f"subject_{index}" for index in range(len(raynaud_dates))],
        "raynaud_date": raynaud_dates,
        "nonraynaud_date": nonraynaud_dates,
        "diagnosis_date": ["2020-06-01"] * len(raynaud_dates),
        "ssc_subtype": ["lcSSc"] * len(raynaud_dates),
        "other_dx": ["ILD"] * len(raynaud_dates),
    })


def test_onset_order_flag_is_unknown_when_either_onset_date_is_missing():
    cleaned = clean.clean_ssc_subtype(subtype_frame(
        raynaud_dates=["2010-01-01", "2015-01-01", None, "2012-01-01"],
        nonraynaud_dates=["2012-01-01", "2010-01-01", "2011-01-01", None]))
    assert list(cleaned["onset_order_flag"]) == [
        "ok", "raynaud_after_nonraynaud", "unknown", "unknown"]
    logged = issue_log()
    assert logged["issue"].eq("raynaud_date or nonraynaud_date missing").any()


# ------------------------------------------------- numeric coercion is logged
def pft_frame(values) -> pd.DataFrame:
    return pd.DataFrame({
        "case_number": [f"subject_{index}" for index in range(len(values))],
        "PFT_dts": ["2020-01-01"] * len(values),
        "NAME": ["FVC"] * len(values),
        "ORD_VALUE": values,
    })


def test_clean_pft_logs_non_numeric_values_and_keeps_blanks_out_of_the_count():
    cleaned = clean.clean_pft(pft_frame(["88.0", None, "not done"]))
    assert cleaned["value"].tolist()[0] == 88.0
    assert cleaned["value"].isna().sum() == 2
    coercion_rows = issue_log().query("issue == 'non-numeric PFT value'")
    assert len(coercion_rows) == 1
    assert coercion_rows["n_affected"].iloc[0] == 1          # the blank is not a coercion


def demographics_frame(heights, weights) -> pd.DataFrame:
    return pd.DataFrame({
        "case number": [f"subject_{index}" for index in range(len(heights))],
        "birth date": ["1970-01-01"] * len(heights),
        "first name": ["Ada"] * len(heights),
        "last name": [f"Lovelace{index}" for index in range(len(heights))],
        "height": heights,
        "weight": weights,
        "ethnicity": ["not hispanic"] * len(heights),
        "gender": ["female"] * len(heights),
        "races": ["white"] * len(heights),
        "diagnosis": ["systemic sclerosis"] * len(heights),
        "state": ["IL"] * len(heights),
    })


def test_clean_demographics_logs_non_numeric_height_and_weight():
    cleaned = clean.clean_demographics(
        demographics_frame(heights=["170", "unknown", None],
                           weights=["70", "70", "see chart"]),
        kg_patients=[])
    assert cleaned["height_cm"].notna().sum() == 1
    logged = issue_log()
    assert logged.query("issue == 'non-numeric height'")["n_affected"].tolist() == [1]
    assert logged.query("issue == 'non-numeric weight'")["n_affected"].tolist() == [1]


# ------------------------------------------------------ duplicate registrations
def raw_tables_for(demographics: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Empty longitudinal tables carrying only the id column each one is read by."""
    return {name: pd.DataFrame({clean.ID_COLS[name]: []}) for name in clean.LONG_TABLES}


def test_an_incomplete_identity_key_is_never_merged_automatically():
    demographics = demographics_frame(heights=["170"] * 3, weights=["70"] * 3)
    demographics["last name"] = ["Lovelace", "Lovelace", "Lovelace"]
    demographics["birth date"] = ["1970-01-01", "1970-01-01", None]
    id_map = clean.find_duplicate_registrations(demographics, raw_tables_for(demographics))
    # the first two share a complete identity key and merge; the third has no birth date
    assert list(id_map["original_id"]) == ["subject_1"]
    incomplete = issue_log().query("issue.str.startswith('incomplete identity key')")
    assert incomplete["n_affected"].tolist() == [1]
    assert "needs adjudication" in incomplete["detail"].iloc[0]
    assert "subject_2" in incomplete["detail"].iloc[0]


def test_a_unique_patient_produces_no_merge_and_no_incomplete_key_row():
    demographics = demographics_frame(heights=["170", "165"], weights=["70", "60"])
    id_map = clean.find_duplicate_registrations(demographics, raw_tables_for(demographics))
    assert id_map.empty
    assert issue_log().empty


# ------------------------------------------------------------ registry linkage
def test_build_processed_fails_fast_when_the_two_registry_tables_cover_different_patients():
    real_raw = clean.load_raw()

    def raw_with_one_extra_id_on_each_side() -> dict[str, pd.DataFrame]:
        tables = {name: frame.copy() for name, frame in real_raw.items()}
        extra_demographics = tables["demographics"].iloc[[0]].copy()
        extra_demographics["case number"] = "subject_demographics_only"
        extra_demographics["first name"] = "Solo"
        tables["demographics"] = pd.concat([tables["demographics"], extra_demographics],
                                           ignore_index=True)
        extra_subtype = tables["ssc_subtype"].iloc[[0]].copy()
        extra_subtype["study_code"] = "subject_subtype_only"
        tables["ssc_subtype"] = pd.concat([tables["ssc_subtype"], extra_subtype],
                                          ignore_index=True)
        return tables

    original_load_raw = clean.load_raw
    clean.load_raw = raw_with_one_extra_id_on_each_side
    try:
        with pytest.raises(ValueError) as raised:
            clean.build_processed()
    finally:
        clean.load_raw = original_load_raw
    message = str(raised.value)
    assert "subject_demographics_only" in message
    assert "subject_subtype_only" in message


# --------------------------------------------------------------- quarantine
def test_quarantine_controls_returns_both_frames_on_a_continuous_index():
    table = pd.DataFrame({
        "subject_id": ["subject_1", "SSC_NORM_0101", "subject_2", "SSC_NORM_0102",
                       "subject_not_registered"],
        "value": [1, 2, 3, 4, 5],
    })
    keep, quarantined = clean.quarantine_controls(table, "vitals", {"subject_1", "subject_2"})
    assert list(keep.index) == list(range(len(keep)))
    assert list(quarantined.index) == list(range(len(quarantined)))
    assert keep["subject_id"].tolist() == ["subject_1", "subject_2"]
    assert len(quarantined) == 3


# --------------------------------------------------------------- FVC slope
def test_the_fvc_slope_is_computed_only_for_three_or_more_distinct_measurement_dates():
    features = pd.read_parquet(ROOT / "data" / "processed" / "features.parquet")
    pft = pd.read_parquet(ROOT / "data" / "processed" / "pft.parquet")
    fvc_dates = pft[pft["measure"].eq("FVC")].groupby("subject_id")["date"]
    eligible_ids = set(fvc_dates.nunique()[fvc_dates.nunique() >= 3].index)
    eligible_ids &= set((fvc_dates.max() - fvc_dates.min())[
        (fvc_dates.max() - fvc_dates.min()) > pd.Timedelta(0)].index)
    with_a_slope = set(features.loc[features["fvc_slope_pct_yr"].notna(), "subject_id"])
    assert with_a_slope == eligible_ids
    # every patient counted here contributes at least three separate time points
    assert fvc_dates.nunique().loc[sorted(with_a_slope)].min() >= 3
