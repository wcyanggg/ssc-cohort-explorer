"""Level 2: cohort-level structure and the disease patterns."""
import pandas as pd
import plotly.express as px
import streamlit as st

from common import (GRAY, STATUS_COLORS, SUBTYPE_COLORS, data_ready, footer,
                    load, page_setup)

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
st.caption("All computed on the cleaned layer; FEV1 excluded (it tracks FVC "
           "as pure noise and carries no independent signal), sentinel codes "
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
        aca_title += (" (ACA marks limited disease, "
                      f"{aca_rate['lcSSc'] / aca_rate['dcSSc']:.0f}×)")
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
    # dx_ild is three-state, so the rate is shown under both readings of an
    # unrecorded other_dx: field recorded only, and everyone with missing = no
    recorded_only = "among patients with a recorded comorbidity field"
    all_patients = "among all patients, missing read as none"
    ild = pd.DataFrame(
        [{"group": k, "denominator": recorded_only,
          "ILD recorded": selected.loc[m, "dx_ild"].mean()}
         for k, m in groups.items() if m.any()]
        + [{"group": k, "denominator": all_patients,
            "ILD recorded": selected.loc[m, "dx_ild"].fillna(False).mean()}
           for k, m in groups.items() if m.any()])
    fig = px.bar(ild, x="group", y="ILD recorded", color="denominator",
                 barmode="group",
                 color_discrete_map={recorded_only: STATUS_COLORS["positive"],
                                     all_patients: GRAY})
    fig.update_layout(height=320, margin=dict(l=10, r=10, t=30, b=10),
                      title="Lung fibrosis follows Scl-70; ACA is protective",
                      yaxis_tickformat=".0%", xaxis_title="",
                      legend=dict(orientation="h", y=-0.25, title=""))
    st.plotly_chart(fig, width="stretch")
    st.caption("`other_dx` is missing for 37% of patients, so each group gets "
               "two bars: the rate among patients whose comorbidity field was "
               "recorded, and the rate among all patients with a missing "
               "field read as none. The true rate lies between the two.")
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

n_slope_patients = int(selected["fvc_slope_pct_yr"].notna().sum())
st.info("Per-patient FVC *slopes* center on zero for both subtypes "
        f"({n_slope_patients} patients with three or more tests). The dataset "
        "carries cross-sectional lung differences and no progressive decline. "
        "This is a negative result and it is reported as such.")

footer()
