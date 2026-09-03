"""Shared data access, labels, and palette for the Streamlit app.

Every page reads the cleaned parquet layer produced by scripts/build.py and
uses the same fixed color assignments. The blue and the red that carry meaning
stay distinguishable under simulated protanopia, deuteranopia and tritanopia;
grays are reserved for borderline and indeterminate.
"""
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssc_coh.config import PROC_DIR, REFERENCE_DATE  # noqa: E402

# ---------------------------------------------------------------- palette
BLUE, RED, TEAL, ORANGE = "#2E6FE0", "#E0483B", "#0E9D8C", "#E88712"
GRAY, LIGHTGRAY = "#94A3B8", "#CBD5E1"

SUBTYPE_COLORS = {"lcSSc": BLUE, "dcSSc": RED}
STATUS_COLORS = {"positive": ORANGE, "negative": BLUE,
                 "borderline": GRAY, "indeterminate": LIGHTGRAY}
MEASURE_COLORS = {"mRSS": RED, "FVC": BLUE, "DLCO": TEAL, "Weight": ORANGE,
                  "BP systolic": BLUE, "BP diastolic": TEAL, "Pulse": ORANGE}

# one color per group label, shared by the Compare and Discover pages
GROUP_COLORS = {**SUBTYPE_COLORS, **STATUS_COLORS,
                "yes": STATUS_COLORS["positive"], "no": STATUS_COLORS["negative"]}

# ---------------------------------------------------------------- table groups
LONG_TABLES = ["vitals", "labs", "pft", "mrss", "medications", "antibodies"]
RESEARCH_TABLES = ["bal", "biopsies", "libraries"]
SOURCE_TABLES = ["subjects", "ssc_subtype"] + LONG_TABLES + RESEARCH_TABLES
MAIN_TABLES = SOURCE_TABLES + ["features"]

# ---------------------------------------------------------------- metadata
TABLE_DESCRIPTIONS = {
    "subjects": "One row per patient: demographics, standardized height/weight, derived age and BMI.",
    "ssc_subtype": "One row per patient: disease subtype, onset milestone dates, recorded comorbidities.",
    "vitals": "Longitudinal vitals (long format): BP, pulse, weight, BMI per visit.",
    "labs": "Longitudinal CBC lab results (long format), sentinel codes removed.",
    "lab_differential_type": "Lab method metadata split out of the numeric lab table.",
    "pft": "Longitudinal lung function (% predicted): FVC, FEV1 (flagged), DLCO.",
    "mrss": "Longitudinal skin score (0-51) with the scoring clinician.",
    "medications": "Prescriptions with brand names mapped to generic.",
    "antibodies": "SSc-specific autoantibody results over time.",
    "bal": "Bronchoalveolar lavage procedures: site, instilled/recovered volume.",
    "biopsies": "Skin biopsies with pathology image references (format mismatches flagged).",
    "libraries": "RNA-seq library prep records (LIMS export, columns normalized).",
    "features": "Derived one-row-per-patient matrix used by Compare and Discover pages.",
    "issues": "The pipeline's decision log: every quality rule that changed or flagged data.",
    "subject_id_map": "Duplicate-registration ids mapped to their canonical id.",
    "controls_vitals": "Quarantined healthy-control rows (SSC_NORM_*).",
    "controls_pft": "Quarantined healthy-control rows (SSC_NORM_*).",
    "controls_mrss": "Quarantined healthy-control rows (SSC_NORM_*).",
    "controls_libraries": "Quarantined healthy-control rows (SSC_NORM_*).",
}

NUMERIC_VARS = {
    "age_years": "Age (years)",
    "bmi_calc": "BMI (derived)",
    "disease_duration_years": "Disease duration (years)",
    "mrss_latest": "mRSS, latest",
    "mrss_mean": "mRSS, mean across visits",
    "fvc": "FVC % predicted, latest",
    "fvc_baseline": "FVC % predicted, first",
    "dlco_sb": "DLCO % predicted, latest",
    "fvc_slope_pct_yr": "FVC slope (% predicted / year)",
    "hemoglobin": "Hemoglobin, latest",
    "wbc": "WBC, latest",
    "platelet count": "Platelets, latest",
    "vit_bp_systolic": "Systolic BP, latest",
    "vit_pulse": "Pulse, latest",
    "weight_kg": "Weight (kg)",
}

GROUP_VARS = {
    "ssc_subtype": "Disease subtype",
    "gender": "Sex",
    "ab_aca": "Anti-centromere (latest)",
    "ab_scl70": "Scl-70 (latest)",
    "ab_rna_pol3": "RNA Pol III (latest)",
    "dx_ild": "ILD recorded",
    "dx_gerd": "GERD recorded",
    "dx_pah": "PAH recorded",
    "ab_result_flip": "Antibody result flipped",
    "onset_order_flag": "Onset-order flag",
}

RATE_OUTCOMES = {
    "dx_ild": "ILD recorded",
    "dx_gerd": "GERD recorded",
    "dx_pah": "PAH recorded",
    "ab_result_flip": "Antibody result flipped",
}

MILESTONES = [("raynaud_date", "Raynaud onset"),
              ("nonraynaud_date", "First non-Raynaud symptom"),
              ("diagnosis_date", "SSc diagnosis")]


# ---------------------------------------------------------------- data access
@st.cache_data(show_spinner=False)
def load(name: str) -> pd.DataFrame:
    return pd.read_parquet(PROC_DIR / f"{name}.parquet")


def page_setup(title: str) -> None:
    st.set_page_config(page_title=f"{title} · SSc Cohort Explorer",
                       page_icon="🫁", layout="wide")


def data_ready() -> bool:
    # features.parquet is written last by scripts/build.py, so its presence means
    # the whole cleaned layer is there, not just an empty directory
    if (PROC_DIR / "features.parquet").exists():
        return True
    st.error("Processed data not found. Run `python scripts/build.py` first "
             "(see README).")
    return False


def yes_no(values: pd.Series) -> pd.Series:
    """Boolean columns read as yes / no on every page. Missing stays missing,
    because a comorbidity that was never recorded is not a no."""
    if pd.api.types.is_bool_dtype(values):
        return values.map({True: "yes", False: "no"})
    return values


def footer() -> None:
    st.caption(f"Fully synthetic SSc cohort (no PHI) · cleaned layer from "
               f"`scripts/build.py` · reference date {REFERENCE_DATE} · "
               f"every cleaning decision is logged on the Data & Quality page.")
