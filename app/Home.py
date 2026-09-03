"""Entry page for the SSc Cohort Explorer."""
import streamlit as st

from common import (LONG_TABLES, RESEARCH_TABLES, SOURCE_TABLES,
                    TABLE_DESCRIPTIONS, data_ready, footer, load, page_setup)

page_setup("Home")

st.title("SSc Cohort Explorer")
st.markdown(
    "An interactive browser for a **synthetic systemic-sclerosis (SSc) registry**: "
    "1,500 registered patients, 11 source tables, clinical follow-up plus a "
    "research sample chain (lavage, then RNA-seq libraries). Built for a reader "
    "with no prior knowledge of the disease or the dataset."
)

if not data_ready():
    st.stop()

subjects = load("subjects")
issues = load("issues")
long_rows = sum(len(load(name)) for name in LONG_TABLES)
research_rows = sum(len(load(name)) for name in RESEARCH_TABLES)

c1, c2, c3, c4 = st.columns(4)
c1.metric("Patients (after cleaning)", f"{len(subjects):,}")
c2.metric("Longitudinal records", f"{long_rows:,}")
c3.metric("Research-sample records", f"{research_rows:,}")
c4.metric("Quality decisions logged", len(issues))

st.markdown("### How the tables connect")
st.graphviz_chart("""
digraph {
  rankdir=LR; node [shape=box, style=rounded, fontsize=11];
  subgraph cluster_reg { label="registry (one row per patient)"; color=gray;
    demographics; ssc_subtype; }
  subgraph cluster_clin { label="clinical follow-up (longitudinal)"; color=gray;
    vitals; labs; pft; mrss; medications; antibodies; }
  subgraph cluster_res { label="research samples"; color=gray;
    bal; biopsies; libraries; }
  demographics -> ssc_subtype [dir=none, label="1:1"];
  demographics -> {vitals labs pft mrss medications antibodies} [dir=none];
  demographics -> bal [dir=none];
  demographics -> biopsies [dir=none];
  bal -> libraries [label="sample"];
}
""")

st.markdown("### Five levels of the data, five pages")
st.markdown(f"""
| Page | Level | What you can do |
|---|---|---|
| **Data & Quality** | raw dimensions & dictionary | every table, column, and the full cleaning-decision log |
| **Cohort** | cohort-level structure | composition and the disease patterns that hold across {len(subjects):,} patients |
| **Compare Groups** | group / variable comparisons | pick any grouping and any measure, with descriptive statistics |
| **Discover Structure** | multidimensional relationships | correlation matrix, PCA, free X-Y exploration |
| **Patient** | one patient's trajectory | full longitudinal timeline down to individual research samples |
""")

with st.expander("What each source table is (plain-language)"):
    for name in SOURCE_TABLES:
        st.markdown(f"**{name}**: {TABLE_DESCRIPTIONS[name]}")

footer()
