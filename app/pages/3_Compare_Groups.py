"""Level 3: user-defined group comparisons with descriptive statistics."""
import pandas as pd
import plotly.express as px
import streamlit as st
from scipy.stats import chi2_contingency, kruskal, mannwhitneyu

from common import (BLUE, GROUP_COLORS, GROUP_VARS, NUMERIC_VARS,
                    RATE_OUTCOMES, data_ready, footer, load, page_setup,
                    yes_no)

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
    outcome_subset = df.dropna(subset=[outcome])
    rate = (outcome_subset.groupby(group_col, observed=True)[outcome]
            .agg(rate="mean", n="size").reset_index())
    fig = px.bar(rate, x=group_col, y="rate", text="n",
                 color_discrete_sequence=[BLUE])
    fig.update_layout(height=420, yaxis_tickformat=".0%",
                      xaxis_title=GROUP_VARS[group_col],
                      yaxis_title=f"{RATE_OUTCOMES[outcome]} (bar label = n)",
                      margin=dict(l=10, r=10, t=20, b=10))
    st.plotly_chart(fig, width="stretch")
    st.dataframe(rate.assign(rate=(rate["rate"] * 100).round(1)),
                 hide_index=True)

    ct = pd.crosstab(outcome_subset[group_col], outcome_subset[outcome])
    if ct.shape[0] >= 2 and ct.shape[1] == 2 and (ct.values >= 5).all():
        chi2, p, _, _ = chi2_contingency(ct)
        st.markdown(f"Chi-square association: **p = {p:.2e}**")
    else:
        st.caption("Cells too small for a chi-square test.")
    if outcome.startswith("dx_"):
        st.caption("Comorbidity fields come from free text missing for 37% "
                   "of patients: rates are lower bounds on prevalence.")

footer()
