"""Subject-level feature matrix for cohort analysis and the discovery page.

Long tables become one row per subject: latest and mean and baseline
measurements, visit counts, antibody status, and event indicators.

A `_latest` column holds the latest valid value, not the value on the latest
visit: the aggregation skips a missing reading and reads back to the most
recent visit that carries a number.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _antibody_status(antibodies: pd.DataFrame) -> pd.DataFrame:
    """Per subject per test: latest result, ever-positive, and whether the test was run.

    A subject appears in a test's frame only if that test has a row for them, so
    presence is what ab_*_tested records. Subjects with no row for a test at all
    come out missing here and are set to False in build_features.
    """
    per_test = []
    for test_short, group in antibodies.groupby("test_short"):
        group = group.sort_values("date")
        status = pd.DataFrame({
            f"ab_{test_short}": group.groupby("subject_id")["value"].last(),
            f"ab_{test_short}_ever_pos": (group.assign(is_positive=group["value"].eq("positive"))
                                          .groupby("subject_id")["is_positive"].any()),
        })
        status[f"ab_{test_short}_tested"] = True
        per_test.append(status)
    out = per_test[0] if per_test else pd.DataFrame()
    for status in per_test[1:]:
        out = out.join(status, how="outer")
    return out.reset_index().rename(columns={"index": "subject_id"})


def _summarize_long(df: pd.DataFrame, value_col: str, prefix: str,
                    cols: list[str] | None = None) -> pd.DataFrame:
    """latest / mean / n per subject for one long table."""
    if cols is not None:
        df = df[df["measure"].isin(cols)]
    if df.empty:
        return pd.DataFrame()
    return (df.sort_values("date")
            .groupby("subject_id")[value_col]
            .agg(**{f"{prefix}_count": "count", f"{prefix}_mean": "mean",
                    f"{prefix}_latest": "last", f"{prefix}_baseline": "first"})
            .reset_index())


def build_features(subjects: pd.DataFrame, ssc: pd.DataFrame,
                   vitals: pd.DataFrame, labs: pd.DataFrame,
                   pft: pd.DataFrame, mrss: pd.DataFrame,
                   meds: pd.DataFrame, antibodies: pd.DataFrame,
                   bal: pd.DataFrame, biopsies: pd.DataFrame,
                   libraries: pd.DataFrame) -> pd.DataFrame:
    """Wide subject-level matrix used by EDA and the app's discovery page."""
    feat = subjects.merge(
        ssc[["subject_id", "ssc_subtype", "dx_ild", "dx_gerd", "dx_pah",
             "onset_order_flag", "disease_duration_years"]],
        on="subject_id", how="left")

    antibody_status = _antibody_status(antibodies)
    feat = feat.merge(antibody_status, on="subject_id", how="left")
    tested_columns = [c for c in antibody_status.columns if c.endswith("_tested")]
    feat[tested_columns] = feat[tested_columns].fillna(False).astype(bool)

    # antibody result flips (SSc antibodies are clinically stable; flips are
    # a data-quality signal surfaced to the app, latest result used above).
    # clean_antibodies() already flagged every row of a flipping patient-test pair.
    flip_any = antibodies.groupby("subject_id")["result_flip_flag"].any().rename("ab_result_flip")
    feat = feat.merge(flip_any.reset_index(), on="subject_id", how="left")
    feat["ab_result_flip"] = feat["ab_result_flip"].fillna(False).astype(bool)

    # longitudinal summaries. "last" reads the last row of each group, so the sort is explicit
    vit_wide = (vitals.sort_values("date")
                .pivot_table(index="subject_id", columns="measure",
                             values="value", aggfunc="last")
                .rename(columns=lambda c: f"vit_{c.lower().replace(' ', '_')}"))
    feat = feat.merge(vit_wide.reset_index(), on="subject_id", how="left")

    # rows without a date carry no position in the series, so they are dropped before "last"
    lab_wide = (labs[labs["date"].notna()].sort_values("date")
                .pivot_table(index="subject_id", columns="component",
                             values="value_num", aggfunc="last")
                .rename(columns=str.lower))
    feat = feat.merge(lab_wide.reset_index(), on="subject_id", how="left")

    pft_latest = (pft.sort_values("date").groupby(["subject_id", "measure"])["value"]
                  .last().unstack().rename(columns=str.lower))
    pft_first = (pft.sort_values("date").groupby(["subject_id", "measure"])["value"]
                 .first().unstack().rename(columns=lambda c: f"{c.lower()}_baseline"))
    feat = feat.merge(pft_latest.reset_index(), on="subject_id", how="left")
    feat = feat.merge(pft_first.reset_index(), on="subject_id", how="left")

    mrss_s = _summarize_long(mrss, "score", "mrss")
    feat = feat.merge(mrss_s, on="subject_id", how="left")

    # visit / event counts
    counts = pd.DataFrame({
        "n_vitals": vitals.groupby("subject_id").size(),
        "n_labs": labs.groupby("subject_id").size(),
        "n_pft": pft.groupby("subject_id").size(),
        "n_mrss": mrss.groupby("subject_id").size(),
        "n_meds": meds.groupby("subject_id").size(),
        "n_antibodies": antibodies.groupby("subject_id").size(),
        "n_bal": bal.groupby("subject_id").size(),
        "n_biopsies": biopsies.groupby("subject_id").size(),
        "n_libraries": libraries.groupby("subject_id").size(),
    }).fillna(0).astype(int).reset_index()
    feat = feat.merge(counts, on="subject_id", how="left")
    feat[counts.columns[1:]] = feat[counts.columns[1:]].fillna(0).astype(int)

    # Exploratory FVC slope (% predicted per year) for patients with at least three
    # distinct FVC measurement dates. Repeat tests on one date add no time points, so
    # the rule counts unique dates rather than rows. Everyone else keeps a missing
    # slope; no annual rate is extrapolated for them.
    slopes = []
    for subject_id, fvc_rows in pft[pft["measure"] == "FVC"].groupby("subject_id"):
        fvc_rows = fvc_rows.sort_values("date")
        if fvc_rows["date"].nunique() >= 3:
            years_from_first_test = (fvc_rows["date"] - fvc_rows["date"].min()).dt.days / 365.25
            if years_from_first_test.max() > 0:
                slopes.append({"subject_id": subject_id,
                               "fvc_slope_pct_yr": np.polyfit(years_from_first_test,
                                                              fvc_rows["value"], 1)[0]})
    if slopes:
        feat = feat.merge(pd.DataFrame(slopes), on="subject_id", how="left")
    else:
        feat["fvc_slope_pct_yr"] = np.nan

    return feat
