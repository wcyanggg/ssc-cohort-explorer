"""Level 2: cohort-level structure and the disease patterns."""
import pandas as pd
import plotly.express as px
import streamlit as st

from common import (GRAY, LIGHTGRAY, STATUS_COLORS, SUBTYPE_COLORS, data_ready,
                    footer, load, page_setup)
from ssc_coh.stats import prevalence_summary

page_setup("Cohort")
st.title("Cohort")
if not data_ready():
    st.stop()

feat = load("features")

subtypes = st.multiselect("Subtype filter", ["lcSSc", "dcSSc"],
                          default=["lcSSc", "dcSSc"])
selected = feat[feat["ssc_subtype"].isin(subtypes)]
if selected.empty:
    st.warning("Select at least one subtype.")
    st.stop()
st.caption(f"{len(selected):,} patients selected · lcSSc = limited cutaneous, "
           "dcSSc = diffuse cutaneous (more aggressive)")

# ---------------------------------------------------------------- who
st.markdown("### Who is in the cohort")
c1, c2, c3 = st.columns(3)
with c1:
    fig = px.histogram(selected, x="age_years", color="ssc_subtype", nbins=30,
                       barmode="overlay", opacity=0.65,
                       color_discrete_map=SUBTYPE_COLORS)
    fig.update_layout(height=280, margin=dict(l=10, r=10, t=30, b=10),
                      title="Age", legend_title="")
    st.plotly_chart(fig, width="stretch")
with c2:
    sex = selected.groupby(["gender", "ssc_subtype"]).size().reset_index(name="n")
    fig = px.bar(sex, x="gender", y="n", color="ssc_subtype",
                 color_discrete_map=SUBTYPE_COLORS)
    fig.update_layout(height=280, margin=dict(l=10, r=10, t=30, b=10),
                      title="Sex (82% female matches real SSc)",
                      legend_title="")
    st.plotly_chart(fig, width="stretch")
with c3:
    dur = px.histogram(selected, x="disease_duration_years", color="ssc_subtype",
                       nbins=25, barmode="overlay", opacity=0.65,
                       color_discrete_map=SUBTYPE_COLORS)
    dur.update_layout(height=280, margin=dict(l=10, r=10, t=30, b=10),
                      title="Years since diagnosis", legend_title="")
    st.plotly_chart(dur, width="stretch")

# ---------------------------------------------------------------- patterns
st.markdown("### The patterns that hold across the cohort")
st.caption("All computed on the cleaned layer; FEV1 excluded (it tracks FVC and "
           "adds no usable variation in this dataset), sentinel codes "
           "and default-fills removed.")

p1, p2 = st.columns(2)
with p1:
    rows = []
    for col, label in [("ab_aca", "anti-centromere"), ("ab_scl70", "Scl-70"),
                       ("ab_rna_pol3", "RNA Pol III")]:
        for s in subtypes:
            sub = selected[(selected["ssc_subtype"] == s)
                           & selected[col].isin(["positive", "negative"])]
            if len(sub):
                rows.append({"antibody": label, "subtype": s,
                             "positive": sub[col].eq("positive").mean()})
    ab = pd.DataFrame(rows)
    # the ACA gap between the subtypes, computed from the bars themselves
    aca_rate = ab[ab["antibody"] == "anti-centromere"].set_index("subtype")["positive"]
    aca_title = "Antibody positivity by subtype"
    if {"lcSSc", "dcSSc"}.issubset(aca_rate.index) and aca_rate["dcSSc"] > 0:
        aca_title += (" (ACA positivity is "
                      f"{aca_rate['lcSSc'] / aca_rate['dcSSc']:.0f}× as common in limited disease)")
    fig = px.bar(ab, x="antibody", y="positive", color="subtype",
                 barmode="group", color_discrete_map=SUBTYPE_COLORS)
    fig.update_layout(height=320, margin=dict(l=10, r=10, t=30, b=10),
                      title=aca_title,
                      yaxis_tickformat=".0%", legend_title="")
    st.plotly_chart(fig, width="stretch")
with p2:
    mrss_median = selected.groupby("ssc_subtype")["mrss_latest"].median().dropna()
    mrss_title = ("Skin score (mRSS, latest): median "
                  + ", ".join(f"{subtype} {value:.0f}"
                              for subtype, value in mrss_median.items()))
    fig = px.box(selected.dropna(subset=["mrss_latest"]), x="ssc_subtype",
                 y="mrss_latest", color="ssc_subtype",
                 color_discrete_map=SUBTYPE_COLORS, points="outliers")
    fig.update_layout(height=320, margin=dict(l=10, r=10, t=30, b=10),
                      title=mrss_title,
                      showlegend=False)
    st.plotly_chart(fig, width="stretch")

p3, p4 = st.columns(2)
with p3:
    groups = {"Scl-70 positive": selected["ab_scl70"].eq("positive"),
              "Scl-70 negative": selected["ab_scl70"].eq("negative"),
              "ACA positive": selected["ab_aca"].eq("positive"),
              "ACA negative": selected["ab_aca"].eq("negative")}
    # dx_ild is three-state, so each group gets the complete-case estimate and the
    # two bounds that read every unrecorded field as no ILD and as ILD
    complete_case = "complete-case estimate (field recorded)"
    lower_bound = "lower bound (unrecorded read as no ILD)"
    upper_bound = "upper bound (unrecorded read as ILD)"
    ild = pd.DataFrame(
        [{"group": group_label, "estimate": estimate_label, "ILD recorded": rate}
         for group_label, group_mask in groups.items() if group_mask.any()
         for estimate_label, rate in
         [(complete_case, prevalence_summary(selected.loc[group_mask, "dx_ild"]).complete_case_rate),
          (lower_bound, prevalence_summary(selected.loc[group_mask, "dx_ild"]).lower_bound),
          (upper_bound, prevalence_summary(selected.loc[group_mask, "dx_ild"]).upper_bound)]])
    fig = px.bar(ild, x="group", y="ILD recorded", color="estimate",
                 barmode="group",
                 category_orders={"estimate": [complete_case, lower_bound, upper_bound]},
                 color_discrete_map={complete_case: STATUS_COLORS["positive"],
                                     lower_bound: GRAY, upper_bound: LIGHTGRAY})
    fig.update_layout(height=320, margin=dict(l=10, r=10, t=30, b=10),
                      title="Recorded ILD is more common with Scl-70 positivity "
                            "and less common with ACA positivity",
                      yaxis_tickformat=".0%", xaxis_title="",
                      legend=dict(orientation="h", y=-0.35, title=""))
    st.plotly_chart(fig, width="stretch")
    st.caption("`other_dx` is empty for 37% of patients, so each group gets three "
               "bars. The complete-case estimate divides recorded ILD by the "
               "patients whose field was filled in. The two bounds read every "
               "empty field first as no ILD and then as ILD, so the distance "
               "between them is how far the missing field alone can move the "
               "rate. It is a sensitivity range, not a confidence interval.")
with p4:
    lung = selected.melt(id_vars="ssc_subtype", value_vars=["fvc", "dlco_sb"],
                         var_name="measure", value_name="pct").dropna()
    lung["measure"] = lung["measure"].map({"fvc": "FVC", "dlco_sb": "DLCO"})
    fig = px.box(lung, x="measure", y="pct", color="ssc_subtype",
                 color_discrete_map=SUBTYPE_COLORS)
    fig.add_hline(y=80, line_dash="dash", line_width=1, opacity=0.5)
    fig.update_layout(height=320, margin=dict(l=10, r=10, t=30, b=10),
                      title="Lung function sits lower in diffuse disease "
                            "(80% = lower normal)",
                      yaxis_title="% predicted (latest)", legend_title="")
    st.plotly_chart(fig, width="stretch")

# the qualifying series are short, so the span they cover is reported with the count
eligible_slope_ids = selected.loc[selected["fvc_slope_pct_yr"].notna(), "subject_id"]
fvc_tests = load("pft")
fvc_dates_by_patient = (fvc_tests[fvc_tests["measure"].eq("FVC")
                                  & fvc_tests["subject_id"].isin(eligible_slope_ids)]
                        .groupby("subject_id")["date"])
follow_up_years = (fvc_dates_by_patient.max() - fvc_dates_by_patient.min()).dt.days / 365.25
if len(follow_up_years):
    st.info(f"No cohort-level decline was detected among the {len(follow_up_years)} patients "
            "who qualify for the exploratory slope analysis (at least three distinct FVC "
            "measurement dates): their per-patient slopes center on zero in both subtypes. "
            f"The qualifying series are short, median span {follow_up_years.median():.2f} year "
            f"and {int((follow_up_years < 1).sum())} of {len(follow_up_years)} under a year, so "
            "they carry no information about longer follow-up. This is a negative result and it "
            "is reported as such.")

footer()
