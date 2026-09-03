"""Level 1: raw dimensions, data dictionary, and the quality-decision log."""
import pandas as pd
import plotly.express as px
import streamlit as st

from common import (BLUE, MAIN_TABLES, SOURCE_TABLES, TABLE_DESCRIPTIONS,
                    data_ready, footer, load, page_setup)

page_setup("Data & Quality")
st.title("Data & quality")
if not data_ready():
    st.stop()

tab_dict, tab_cov, tab_issues, tab_quar = st.tabs(
    ["Dictionary", "Coverage & linkage", "Issue log", "Quarantine"])

with tab_dict:
    name = st.selectbox("Table", MAIN_TABLES)
    df = load(name)
    st.markdown(f"**{TABLE_DESCRIPTIONS[name]}**")
    st.caption(f"{len(df):,} rows × {df.shape[1]} columns")

    info = pd.DataFrame({
        "column": df.columns,
        "dtype": [str(t) for t in df.dtypes],
        "non-null %": [(df[c].notna().mean() * 100).round(1) for c in df.columns],
        "unique values": [df[c].nunique() for c in df.columns],
    })
    left, right = st.columns([1, 2])
    with left:
        st.dataframe(info, hide_index=True, height=420)
    with right:
        st.markdown("Sample rows")
        st.dataframe(df.head(10), hide_index=True, height=420)

with tab_cov:
    n_registered = load("subjects")["subject_id"].nunique()
    st.markdown(
        f"Unique patients per table, out of **{n_registered:,}** in the cleaned "
        "registry. Research tables cover nested subsets by design (only part of "
        "the cohort undergoes bronchoscopy or biopsy), which is structure "
        "rather than missingness."
    )
    rows = [{"table": name, "patients": load(name)["subject_id"].nunique()}
            for name in SOURCE_TABLES]
    cov = pd.DataFrame(rows).sort_values("patients")
    fig = px.bar(cov, x="patients", y="table", orientation="h",
                 color_discrete_sequence=[BLUE], text="patients")
    fig.update_layout(height=420, margin=dict(l=10, r=10, t=10, b=10),
                      xaxis_title="unique patients", yaxis_title="")
    st.plotly_chart(fig, width="stretch")
    st.markdown(
        "- Subject ids arrive under **four different column names** "
        "(`case_number`, `reg_id`, `case number`, `study_code`); the cleaning "
        "layer normalizes all of them to `subject_id`.\n"
        "- **4 healthy-control subjects** (`SSC_NORM_*`) from a companion "
        "study leaked into the export; their rows are quarantined, not "
        "deleted (see the Quarantine tab).\n"
        "- **2 duplicate registrations** (same name + birth date under two "
        "ids) were merged; all longitudinal rows follow the canonical id."
    )

with tab_issues:
    issues = load("issues")
    st.markdown(
        "Every rule in the cleaning pipeline that changed or flagged data "
        "wrote a row here. Clear errors are set to missing and logged, "
        "structural artifacts are kept and flagged and left out of the "
        "analyses they would distort, and ambiguous values are flagged "
        "only, never silently altered."
    )
    pick = st.multiselect("Filter by table", sorted(issues["table"].unique()))
    view = issues[issues["table"].isin(pick)] if pick else issues
    st.dataframe(view, hide_index=True, height=480)
    st.download_button("Download issue log (CSV)",
                       issues.to_csv(index=False).encode(),
                       "issues.csv", "text/csv")

with tab_quar:
    st.markdown(
        "Rows recorded under the four `SSC_NORM_*` control subjects. They do "
        "not belong to the SSc registry, but deleting them would hide a "
        "linkage problem, so they are kept visible here."
    )
    for t in ("controls_vitals", "controls_pft", "controls_mrss",
              "controls_libraries"):
        q = load(t)
        st.markdown(f"**{t}**: {len(q)} rows")
        st.dataframe(q, hide_index=True)

footer()
