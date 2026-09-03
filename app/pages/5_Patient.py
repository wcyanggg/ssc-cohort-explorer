"""Level 5: one patient's longitudinal trajectory, down to research samples."""
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots

from common import (MEASURE_COLORS, MILESTONES, data_ready, footer, load,
                    page_setup)

page_setup("Patient")
st.title("Patient trajectory")
if not data_ready():
    st.stop()

subjects = load("subjects")
ssc = load("ssc_subtype")

sid = st.selectbox("Patient", subjects["subject_id"].sort_values(),
                   help="Type to search, e.g. 2005")

person = subjects.set_index("subject_id").loc[sid]
disease = ssc.set_index("subject_id").loc[sid]

c1, c2, c3, c4, c5 = st.columns(5)
c1.metric("Age", f"{person.age_years:.0f}")
c2.metric("Sex", person.gender)
c3.metric("Subtype", disease.ssc_subtype)
c4.metric("Years since diagnosis", f"{disease.disease_duration_years:.1f}")
comorbidities = [("ILD", disease.dx_ild), ("GERD", disease.dx_gerd),
                 ("PAH", disease.dx_pah)]
recorded = [name for name, value in comorbidities
            if not pd.isna(value) and bool(value)]
if recorded:
    comorbidity_text = ", ".join(recorded)
elif all(not pd.isna(value) for _, value in comorbidities):
    comorbidity_text = "none"          # other_dx was filled in and named none of the three
else:
    comorbidity_text = "not recorded"  # other_dx was empty, so the flags are unknown
c5.metric("Recorded comorbidities", comorbidity_text)

ab = load("antibodies")
ab_p = ab[ab.subject_id == sid].sort_values("date")
latest_ab = ab_p.groupby("test_short")["value"].last()
st.caption(
    "Latest antibodies: "
    + " · ".join(f"**{t}**: {latest_ab.get(t, 'not tested')}"
                 for t in ("aca", "scl70", "rna_pol3"))
    + (f" · onset-order flag: {disease.onset_order_flag}"
       if disease.onset_order_flag != "ok" else "")
)

# ---------------------------------------------------------------- timeline
mrss = load("mrss")
pft = load("pft")
vit = load("vitals")
mrss_p = mrss[mrss.subject_id == sid]
pft_p = pft[pft.subject_id == sid]
vit_p = vit[vit.subject_id == sid]

panels = []
if len(mrss_p):
    panels.append(("Skin score (mRSS, 0-51)", [
        go.Scatter(x=mrss_p.date, y=mrss_p.score, mode="lines+markers",
                   name="mRSS", line=dict(color=MEASURE_COLORS["mRSS"], width=2),
                   marker=dict(size=8))]))
if len(pft_p):
    traces = []
    for m, label in (("FVC", "FVC"), ("DLCO_SB", "DLCO")):
        d = pft_p[pft_p.measure == m].dropna(subset=["value"])
        if len(d):
            traces.append(go.Scatter(
                x=d.date, y=d.value, mode="lines+markers", name=label,
                line=dict(color=MEASURE_COLORS[label], width=2),
                marker=dict(size=8)))
    if traces:
        panels.append(("Lung function (% predicted; FEV1 excluded as "
                       "uninformative, it tracks FVC)", traces))
weights = vit_p[vit_p.measure == "WEIGHT IN POUND"].dropna(subset=["value"])
if len(weights):
    panels.append(("Weight (lb)", [
        go.Scatter(x=weights.date, y=weights.value, mode="lines+markers", name="Weight",
                   line=dict(color=MEASURE_COLORS["Weight"], width=2),
                   marker=dict(size=7))]))
systolic = vit_p[vit_p.measure == "BP SYSTOLIC"].dropna(subset=["value"])
diastolic = vit_p[vit_p.measure == "BP DIASTOLIC"].dropna(subset=["value"])
if len(systolic) or len(diastolic):
    traces = []
    if len(systolic):
        traces.append(go.Scatter(x=systolic.date, y=systolic.value, mode="lines+markers",
                                 name="systolic", marker=dict(size=7),
                                 line=dict(color=MEASURE_COLORS["BP systolic"], width=2)))
    if len(diastolic):
        traces.append(go.Scatter(x=diastolic.date, y=diastolic.value, mode="lines+markers",
                                 name="diastolic", marker=dict(size=7),
                                 line=dict(color=MEASURE_COLORS["BP diastolic"], width=2)))
    panels.append(("Blood pressure", traces))

if panels:
    fig = make_subplots(rows=len(panels), cols=1, shared_xaxes=True,
                        subplot_titles=[t for t, _ in panels],
                        vertical_spacing=0.07)
    for i, (_, traces) in enumerate(panels, start=1):
        for tr in traces:
            fig.add_trace(tr, row=i, col=1)
    for col, label in MILESTONES:
        d = disease[col]
        if pd.notna(d):
            fig.add_vline(x=d, line_dash="dot", line_width=1,
                          line_color="#64748B")
            fig.add_annotation(x=d, y=1.02, yref="paper", text=label,
                               showarrow=False, font=dict(size=10, color="#64748B"),
                               textangle=-25)
    fig.update_layout(height=210 * len(panels) + 80,
                      margin=dict(l=10, r=10, t=60, b=10),
                      legend=dict(orientation="h", y=-0.05))
    st.plotly_chart(fig, width="stretch")
    st.caption("Dotted lines: disease milestones, in order Raynaud onset, "
               "first non-Raynaud symptom, diagnosis.")
else:
    st.info("No longitudinal measurements recorded for this patient.")

# ---------------------------------------------------------------- detail tabs
meds = load("medications")
labs = load("labs")
bal = load("bal")
biopsies = load("biopsies")
libraries = load("libraries")

t_med, t_ab, t_lab, t_res = st.tabs(
    ["Medications", "Antibody history", "Labs", "Research samples"])

with t_med:
    m = meds[meds.subject_id == sid].sort_values("date")
    if len(m):
        st.dataframe(m.drop(columns="subject_id"), hide_index=True)
    else:
        st.caption("No medication records.")

with t_ab:
    if len(ab_p):
        st.dataframe(ab_p.drop(columns="subject_id"), hide_index=True)
    else:
        st.caption("No antibody tests.")

with t_lab:
    lp = labs[labs.subject_id == sid]
    if len(lp):
        comp = st.selectbox("Component", sorted(lp.component.unique()))
        d = lp[lp.component == comp].sort_values("date")
        fig = go.Figure(go.Scatter(x=d.date, y=d.value_num,
                                   mode="lines+markers", marker=dict(size=7),
                                   line=dict(color=MEASURE_COLORS["FVC"], width=2)))
        fig.update_layout(height=300, margin=dict(l=10, r=10, t=10, b=10),
                          yaxis_title=comp)
        st.plotly_chart(fig, width="stretch")
    else:
        st.caption("No lab results.")

with t_res:
    bal_p = bal[bal.subject_id == sid].sort_values("date")
    st.markdown(f"**BAL procedures**: {len(bal_p)}")
    if len(bal_p):
        st.dataframe(bal_p.drop(columns="subject_id"), hide_index=True)
    biopsies_p = biopsies[biopsies.subject_id == sid].sort_values("date")
    st.markdown(f"**Skin biopsies**: {len(biopsies_p)}")
    if len(biopsies_p):
        st.dataframe(biopsies_p.drop(columns="subject_id"), hide_index=True)
    libraries_p = libraries[libraries.subject_id == sid]
    st.markdown(f"**RNA-seq library records**: {len(libraries_p)}")
    if len(libraries_p):
        st.dataframe(libraries_p.drop(columns="subject_id"), hide_index=True)

footer()
