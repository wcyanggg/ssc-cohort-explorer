"""Level 4: open-ended structure discovery in the multidimensional data."""
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from common import (GROUP_COLORS, GROUP_VARS, NUMERIC_VARS, data_ready,
                    footer, load, page_setup, yes_no)

page_setup("Discover Structure")
st.title("Discover structure")
st.caption("Not a fixed set of summary charts: correlation structure, a PCA "
           "projection, and free X-Y exploration over the patient-level "
           "feature matrix.")
if not data_ready():
    st.stop()

feat = load("features")

tab_corr, tab_pca, tab_xy = st.tabs(
    ["Correlation matrix", "PCA projection", "X-Y explorer"])

with tab_corr:
    method = st.radio("Correlation", ["spearman", "pearson"], horizontal=True)
    cols = st.multiselect("Variables", list(NUMERIC_VARS),
                          default=list(NUMERIC_VARS),
                          format_func=NUMERIC_VARS.get)
    if len(cols) >= 2:
        corr = feat[cols].corr(method=method)
        labels = [NUMERIC_VARS[c] for c in corr.columns]
        fig = px.imshow(corr.values, x=labels, y=labels, zmin=-1, zmax=1,
                        color_continuous_scale="RdBu_r", aspect="auto",
                        text_auto=".2f")
        fig.update_layout(height=560, margin=dict(l=10, r=10, t=10, b=10))
        st.plotly_chart(fig, width="stretch")
        st.caption("Diverging scale centered at 0. Pairwise complete "
                   "observations; DLCO and slope columns have substantial "
                   "missingness, so their cells rest on fewer patients.")

with tab_pca:
    color_by = st.selectbox("Color points by", list(GROUP_VARS),
                            format_func=GROUP_VARS.get)
    cols = [c for c in NUMERIC_VARS if c != "fvc_slope_pct_yr"]
    features_raw = feat[cols]
    # median imputation is applied and documented: PCA needs complete rows and
    # dropping every row with any gap would discard most of the cohort (the latest
    # DLCO alone is missing for most patients; the share is computed below)
    dlco_missing_share = feat["dlco_sb"].isna().mean()
    imputed = features_raw.fillna(features_raw.median(numeric_only=True))
    scaled = StandardScaler().fit_transform(imputed)
    pca = PCA(n_components=2, random_state=0)
    projection = pca.fit_transform(scaled)
    plot = pd.DataFrame({"PC1": projection[:, 0], "PC2": projection[:, 1],
                         "subject": feat["subject_id"],
                         "group": yes_no(feat[color_by]).astype(str)})
    fig = px.scatter(plot, x="PC1", y="PC2", color="group",
                     color_discrete_map=GROUP_COLORS, hover_name="subject",
                     opacity=0.7)
    fig.update_traces(marker=dict(size=6))
    fig.update_layout(height=520, margin=dict(l=10, r=10, t=10, b=10),
                      legend_title=GROUP_VARS[color_by])
    st.plotly_chart(fig, width="stretch")
    ev = pca.explained_variance_ratio_
    st.caption(f"Explained variance: PC1 {ev[0]:.0%}, PC2 {ev[1]:.0%}. "
               "Standardized features, median-imputed. The latest "
               f"DLCO is missing for {dlco_missing_share:.0%} of patients, so for them it "
               "is the cohort median, and imputation pulls incomplete patients toward "
               "the center.")
    load_tbl = (pd.DataFrame(pca.components_.T, index=[NUMERIC_VARS[c] for c in cols],
                             columns=["PC1", "PC2"])
                .assign(strongest=lambda d: d.abs().max(axis=1))
                .sort_values("strongest", ascending=False).drop(columns="strongest")
                .round(2).head(8))
    st.markdown("Top loadings (what the axes mean):")
    st.dataframe(load_tbl)

with tab_xy:
    c1, c2, c3 = st.columns(3)
    x = c1.selectbox("X", list(NUMERIC_VARS), index=3,
                     format_func=NUMERIC_VARS.get)
    y = c2.selectbox("Y", list(NUMERIC_VARS), index=5,
                     format_func=NUMERIC_VARS.get)
    color_by = c3.selectbox("Color by", list(GROUP_VARS),
                            format_func=GROUP_VARS.get, key="xy_color")
    sub = feat.dropna(subset=[x, y])
    fig = px.scatter(sub, x=x, y=y, color=yes_no(sub[color_by]).astype(str),
                     color_discrete_map=GROUP_COLORS,
                     hover_name="subject_id", opacity=0.7,
                     labels={x: NUMERIC_VARS[x], y: NUMERIC_VARS[y]})
    fig.update_traces(marker=dict(size=6))
    fig.update_layout(height=520, margin=dict(l=10, r=10, t=10, b=10),
                      legend_title=GROUP_VARS[color_by])
    st.plotly_chart(fig, width="stretch")
    if len(sub) > 10:
        rho = sub[[x, y]].corr(method="spearman").iloc[0, 1]
        st.caption(f"Spearman ρ = {rho:.2f} over {len(sub):,} patients with "
                   "both values present.")

footer()
