"""Level 3: user-defined group comparisons with descriptive statistics."""
import pandas as pd
import plotly.express as px
import streamlit as st
from scipy.stats import kruskal, mannwhitneyu

from common import (BLUE, GROUP_COLORS, GROUP_VARS, NUMERIC_VARS,
                    RATE_OUTCOMES, data_ready, footer, load, page_setup,
                    yes_no)
from ssc_coh.stats import association_test, prevalence_by_group

page_setup("Compare Groups")
st.title("Compare groups")
if not data_ready():
    st.stop()

feat = load("features")

c1, c2, c3 = st.columns([1.2, 1.2, 1])
group_col = c1.selectbox("Group patients by", list(GROUP_VARS),
                         format_func=GROUP_VARS.get)
mode = c2.radio("Compare", ["a numeric measure", "an outcome rate"],
                horizontal=True)
hide_amb = c3.checkbox("Hide borderline / indeterminate", value=True)

df = feat.copy()
df[group_col] = yes_no(df[group_col])          # dx_* are nullable boolean
if hide_amb:
    df = df[~df[group_col].isin(["borderline", "indeterminate"])]
df = df.dropna(subset=[group_col])

if mode == "a numeric measure":
    metric = st.selectbox("Measure", list(NUMERIC_VARS),
                          format_func=NUMERIC_VARS.get, index=3)
    numeric_subset = df.dropna(subset=[metric])
    groups = sorted(numeric_subset[group_col].astype(str).unique())
    samples = [numeric_subset.loc[numeric_subset[group_col].astype(str) == g, metric]
               for g in groups]

    fig = px.box(numeric_subset, x=group_col, y=metric,
                 color=numeric_subset[group_col].astype(str),
                 color_discrete_map=GROUP_COLORS, points="outliers")
    fig.update_layout(height=420, showlegend=False,
                      xaxis_title=GROUP_VARS[group_col],
                      yaxis_title=NUMERIC_VARS[metric],
                      margin=dict(l=10, r=10, t=20, b=10))
    st.plotly_chart(fig, width="stretch")

    stats = pd.DataFrame({
        "group": groups,
        "n": [len(s) for s in samples],
        "median": [round(s.median(), 2) for s in samples],
        "IQR": [f"{s.quantile(.25):.1f} to {s.quantile(.75):.1f}" for s in samples],
    })
    st.dataframe(stats, hide_index=True)

    if len(groups) == 2 and all(len(s) >= 5 for s in samples):
        u, p = mannwhitneyu(samples[0], samples[1])
        st.markdown(f"Mann-Whitney U ({groups[0]} vs {groups[1]}): "
                    f"**p = {p:.2e}**")
    elif len(groups) > 2 and all(len(s) >= 5 for s in samples):
        h, p = kruskal(*samples)
        st.markdown(f"Kruskal-Wallis across {len(groups)} groups: "
                    f"**p = {p:.2e}**")
    else:
        st.caption("Groups too small for a distribution test; showing "
                   "descriptives only.")
    st.caption("Rank-based tests throughout: several measures are skewed and "
               "a few groups are small, so no normality assumption is made. "
               "p-values are descriptive here (no multiplicity correction); "
               "treat them as screening, not confirmation.")
else:
    # the grouping variable cannot also be the outcome
    outcome = st.selectbox("Outcome", [o for o in RATE_OUTCOMES if o != group_col],
                           format_func=RATE_OUTCOMES.get)
    # the prevalence of each group is computed on the filtered group universe, before the
    # unrecorded outcomes are dropped, so every group carries its own missingness and its
    # own bounds rather than the cohort's
    group_prevalence = prevalence_by_group(df, group_col, outcome)
    group_prevalence["bar_label"] = [f"{recorded} of {total}" for recorded, total
                                     in zip(group_prevalence["n_recorded"],
                                            group_prevalence["n_all"])]
    outcome_subset = df.dropna(subset=[outcome])
    fig = px.bar(group_prevalence, x="group", y="complete_case_rate", text="bar_label",
                 color_discrete_sequence=[BLUE])
    fig.update_layout(height=420, yaxis_tickformat=".0%",
                      xaxis_title=GROUP_VARS[group_col],
                      yaxis_title=f"{RATE_OUTCOMES[outcome]}, complete-case "
                                  f"(bar label = recorded of total)",
                      margin=dict(l=10, r=10, t=20, b=10))
    st.plotly_chart(fig, width="stretch")
    rate_columns = ["complete_case_rate", "lower_bound", "upper_bound"]
    st.dataframe(group_prevalence[["group", "n_all", "n_recorded", "n_missing"] + rate_columns]
                 .assign(**{column: (group_prevalence[column] * 100).round(1)
                            for column in rate_columns}),
                 hide_index=True)

    # the chi-square approximation is judged on expected counts, not on the observed
    # cells; a 2x2 table that fails the check falls back to Fisher's exact test
    contingency = pd.crosstab(outcome_subset[group_col], outcome_subset[outcome])
    association = association_test(contingency)
    if association.test_name == "none":
        st.caption(f"No association test is shown: {association.reason}. The rates and "
                   "counts above are the result.")
    else:
        st.markdown(f"{association.test_name} association: "
                    f"**p = {association.p_value:.2e}** ({association.reason})")
    unrecorded_total = int(group_prevalence["n_missing"].sum())
    if unrecorded_total:
        own_missingness = (f"lower_bound reads every unrecorded patient in the group as a negative "
                           f"and upper_bound reads every one as a positive, so each group is bounded "
                           f"by its own missingness: the n_missing column shows how the "
                           f"{unrecorded_total} unrecorded patients of the "
                           f"{int(group_prevalence['n_all'].sum())} in this selection fall across the "
                           f"groups. Those bounds are a sensitivity range, not a confidence interval, "
                           f"and two groups whose complete-case rates differ can still have "
                           f"overlapping bounds.")
    else:
        own_missingness = ("The outcome is recorded for every patient in this selection, so each "
                           "group's bounds sit on its complete-case rate.")
    st.caption(f"Each bar is the complete-case rate of its group: the patients with "
               f"{RATE_OUTCOMES[outcome]} divided by the patients in that group whose field was "
               f"recorded. {own_missingness} The association test above uses only the "
               f"{len(outcome_subset)} patients whose outcome was recorded.")

footer()
