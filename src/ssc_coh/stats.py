"""Shared statistics: prevalence under partial recording, and the test gate.

The comorbidity flags are three-state: true, false, and unknown where the source
field was empty. One function computes every rate the notebook and the app show,
so the three numbers are defined once. A second function decides which
association test a contingency table can carry.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact


@dataclass(frozen=True)
class PrevalenceSummary:
    """Counts and the three rates for one nullable-boolean column.

    complete_case_rate uses the recorded patients as the denominator. lower_bound
    reads every unrecorded patient as a negative and upper_bound reads every one
    of them as a positive, so the two are the widest pair of rates the data allow
    without an assumption about why the field is empty.
    """

    n_all: int
    n_recorded: int
    n_missing: int
    n_positive: int
    complete_case_rate: float
    lower_bound: float
    upper_bound: float

    def as_dict(self) -> dict:
        return asdict(self)


def prevalence_summary(flag_values: pd.Series) -> PrevalenceSummary:
    """Prevalence of a nullable-boolean flag under three denominators.

    complete_case_rate = positives / recorded
    lower_bound        = positives / all patients
    upper_bound        = (positives + unrecorded) / all patients

    With nothing recorded the complete-case rate is undefined and comes back as
    NaN; the bounds are still 0 and 1. With an empty input every rate is NaN.
    """
    n_all = int(len(flag_values))
    n_recorded = int(flag_values.notna().sum())
    n_missing = n_all - n_recorded
    n_positive = int((flag_values == True).sum())        # noqa: E712 - pd.NA is not True
    complete_case_rate = n_positive / n_recorded if n_recorded else float("nan")
    lower_bound = n_positive / n_all if n_all else float("nan")
    upper_bound = (n_positive + n_missing) / n_all if n_all else float("nan")
    return PrevalenceSummary(n_all=n_all, n_recorded=n_recorded, n_missing=n_missing,
                             n_positive=n_positive, complete_case_rate=complete_case_rate,
                             lower_bound=lower_bound, upper_bound=upper_bound)


MIN_EXPECTED_COUNT = 5           # the smallest expected cell a chi-square approximation needs


@dataclass(frozen=True)
class AssociationTest:
    """Which test a contingency table supports, and the result if it supports one."""

    test_name: str
    p_value: float
    smallest_expected_count: float
    reason: str


def association_test(contingency: pd.DataFrame) -> AssociationTest:
    """Pick the association test the table can carry, judged on expected counts.

    The chi-square approximation depends on the expected frequencies, not on the
    observed cells, so the expected table decides. Every expected count at or above
    MIN_EXPECTED_COUNT gives a chi-square. A 2x2 table below it falls back to
    Fisher's exact test. A larger table below it gets no test and the reason why.
    """
    if contingency.shape[0] < 2 or contingency.shape[1] < 2:
        return AssociationTest("none", float("nan"), float("nan"),
                               "only one group or one outcome value is present")

    chi_square, chi_square_p, _, expected_counts = chi2_contingency(contingency)
    smallest_expected = float(expected_counts.min())
    if smallest_expected >= MIN_EXPECTED_COUNT:
        return AssociationTest("chi-square", float(chi_square_p), smallest_expected,
                               f"smallest expected count {smallest_expected:.1f}")
    if contingency.shape == (2, 2):
        _, fisher_p = fisher_exact(contingency)
        return AssociationTest("Fisher's exact", float(fisher_p), smallest_expected,
                               f"smallest expected count {smallest_expected:.1f}, below the "
                               f"{MIN_EXPECTED_COUNT} a chi-square needs")
    return AssociationTest("none", float("nan"), smallest_expected,
                           f"smallest expected count {smallest_expected:.1f}, below the "
                           f"{MIN_EXPECTED_COUNT} a chi-square needs, and the table is "
                           f"{contingency.shape[0]} by {contingency.shape[1]}, so Fisher's "
                           "exact test does not apply either")
