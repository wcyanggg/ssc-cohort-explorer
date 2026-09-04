"""Rule-level tests: the cleaning decisions and the shared statistics helpers.

Run:  .venv/bin/python -m pytest tests -q
"""
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssc_coh import clean                            # noqa: E402
from ssc_coh.stats import (association_test, prevalence_by_group,   # noqa: E402
                           prevalence_summary)


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


# ----------------------------------------------- per-group bounds and their overlap
# the three comparisons the report and the app quote, as (lower bound, upper bound) in percent
ILD_BOUNDS_BY_GROUP = {
    "Scl-70 positive": (64.7, 87.7),
    "Scl-70 negative": (32.3, 75.8),
    "anti-centromere positive": (18.1, 69.3),
    "anti-centromere negative": (50.4, 82.3),
    "dcSSc": (63.2, 89.8),
    "lcSSc": (28.3, 72.3),
}
ILD_COMPARISONS = [("Scl-70 positive", "Scl-70 negative"),
                   ("anti-centromere positive", "anti-centromere negative"),
                   ("dcSSc", "lcSSc")]


def ild_bounds_by_group() -> dict[str, tuple[float, float]]:
    """The ILD lower and upper bound of each group, in percent, from the committed features."""
    features = pd.read_parquet(ROOT / "data" / "processed" / "features.parquet")
    group_masks = {
        "Scl-70 positive": features["ab_scl70"].eq("positive"),
        "Scl-70 negative": features["ab_scl70"].eq("negative"),
        "anti-centromere positive": features["ab_aca"].eq("positive"),
        "anti-centromere negative": features["ab_aca"].eq("negative"),
        "dcSSc": features["ssc_subtype"].eq("dcSSc"),
        "lcSSc": features["ssc_subtype"].eq("lcSSc"),
    }
    bounds = {}
    for group_name, group_mask in group_masks.items():
        summary = prevalence_summary(features.loc[group_mask, "dx_ild"])
        bounds[group_name] = (round(summary.lower_bound * 100, 1),
                              round(summary.upper_bound * 100, 1))
    return bounds


def test_each_group_carries_the_ild_bounds_the_report_quotes():
    assert ild_bounds_by_group() == ILD_BOUNDS_BY_GROUP


def test_the_ild_bounds_of_the_two_groups_overlap_in_every_comparison():
    bounds = ild_bounds_by_group()
    for first_group, second_group in ILD_COMPARISONS:
        first_low, first_high = bounds[first_group]
        second_low, second_high = bounds[second_group]
        assert first_low <= second_high and second_low <= first_high, (
            f"{first_group} and {second_group} no longer overlap")


def test_the_ild_bounds_leave_room_for_the_ordering_to_reverse():
    """Overlapping bounds mean an assignment exists that flips the complete-case ordering.

    The group with the higher complete-case rate can be pushed to its lower bound while the
    other is pushed to its upper bound, so no claim that the direction is robust to
    differential missingness can be made from these bounds.
    """
    bounds = ild_bounds_by_group()
    for higher_group, lower_group in [("Scl-70 positive", "Scl-70 negative"),
                                      ("anti-centromere negative", "anti-centromere positive"),
                                      ("dcSSc", "lcSSc")]:
        assert bounds[higher_group][0] < bounds[lower_group][1], (
            f"{higher_group} can no longer fall below {lower_group}")


def test_opposite_assignments_in_the_two_groups_reverse_the_complete_case_ordering():
    # group_high: 6 of 8 recorded positive (75.0%) with 12 unrecorded
    # group_low: 5 of 10 recorded positive (50.0%) with 2 unrecorded
    frame = pd.DataFrame({
        "group": ["group_high"] * 20 + ["group_low"] * 12,
        "outcome": pd.array([True] * 6 + [False] * 2 + [pd.NA] * 12
                            + [True] * 5 + [False] * 5 + [pd.NA] * 2, dtype="boolean"),
    })
    by_group = prevalence_by_group(frame, "group", "outcome").set_index("group")
    assert by_group.loc["group_high", "complete_case_rate"] == 0.75
    assert by_group.loc["group_low", "complete_case_rate"] == 0.5
    # read group_high's unrecorded rows as negatives and group_low's as positives
    assert by_group.loc["group_high", "lower_bound"] == 0.3
    assert by_group.loc["group_low", "upper_bound"] == 7 / 12
    assert by_group.loc["group_high", "lower_bound"] < by_group.loc["group_low", "upper_bound"]


def test_prevalence_by_group_keeps_every_group_on_its_own_denominator():
    # group_sparse is unrecorded for 6 of 10; group_complete is recorded for all 10
    frame = pd.DataFrame({
        "group": ["group_sparse"] * 10 + ["group_complete"] * 10,
        "outcome": pd.array([True] * 2 + [False] * 2 + [pd.NA] * 6
                            + [True] * 5 + [False] * 5, dtype="boolean"),
    })
    by_group = prevalence_by_group(frame, "group", "outcome").set_index("group")
    cohort = prevalence_summary(frame["outcome"])
    assert (by_group.loc["group_sparse", "lower_bound"],
            by_group.loc["group_sparse", "upper_bound"]) == (0.2, 0.8)
    # a group with nothing unrecorded has no width at all, which cohort bounds would hide
    assert (by_group.loc["group_complete", "lower_bound"],
            by_group.loc["group_complete", "upper_bound"]) == (0.5, 0.5)
    assert (cohort.lower_bound, cohort.upper_bound) == (0.35, 0.65)
    for group_name in ["group_sparse", "group_complete"]:
        assert (by_group.loc[group_name, "lower_bound"],
                by_group.loc[group_name, "upper_bound"]) != (cohort.lower_bound,
                                                             cohort.upper_bound)


# ------------------------------------------------------- association test gate
def test_association_test_uses_chi_square_when_every_expected_count_is_large():
    result = association_test(pd.DataFrame([[50, 50], [40, 60]]))
    assert result.test_name == "chi-square"
    assert result.smallest_expected_count >= 5
    assert result.statistic_name == "chi-square"
    assert result.statistic > 0


def test_association_test_falls_back_to_fisher_when_an_expected_count_is_small():
    # observed cells look usable but one expected count is far below 5
    result = association_test(pd.DataFrame([[1, 1], [1, 20]]))
    assert result.test_name == "Fisher's exact"
    assert result.smallest_expected_count < 5
    assert result.statistic_name == "odds ratio"
    assert result.statistic > 0
    assert 0 <= result.p_value <= 1


def test_association_test_reports_descriptives_only_for_a_small_table_larger_than_two_by_two():
    result = association_test(pd.DataFrame([[1, 1, 1], [1, 1, 20], [2, 2, 2]]))
    assert result.test_name == "none"
    assert "does not apply" in result.reason


def test_association_test_needs_two_groups_and_two_outcome_values():
    result = association_test(pd.DataFrame([[5, 5]]))
    assert result.test_name == "none"
    assert "only one group" in result.reason


def test_the_notebook_uses_the_shared_association_test():
    """One association gate, not two: the notebook calls into ssc_coh.stats like the app."""
    notebook = json.loads((ROOT / "notebooks" / "01_disease_patterns.ipynb").read_text())
    notebook_code = "\n".join("".join(cell["source"]) for cell in notebook["cells"]
                              if cell["cell_type"] == "code")
    assert "from ssc_coh.stats import association_test" in notebook_code
    assert "chi2_contingency" not in notebook_code


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
