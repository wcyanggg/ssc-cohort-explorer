"""Cleaning layer: every data-quality decision as explicit, logged code.

Input: the 11 raw tables from raw.py. Output: the cleaned long tables (same
shape as the source, one row per measurement or visit), one row per patient
for the two registry tables, the quarantined control rows, and the issue log.

The rules come from the register in REPORT.md (section 3), which groups them
by treatment:

  clear errors          -> missing + a row in the issue log
                           (sentinel codes, the 160.6 default weight, impossible
                           doses, one weight impossible in any unit)
  structural artifacts  -> kept + flagged, excluded from specific analyses
                           (FEV1, blood-pressure defaults, clipped PFT values)
  ambiguous values      -> a flag column, never altered
                           (onset order, antibody flips, same-day scores,
                           timeline paradoxes, lab status contradictions)
  identity and structure -> deterministic fixes, logged
                           (four id column names, duplicate registrations,
                           control subjects from a companion study)

Nothing here modifies the raw files.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .config import (
    ANTIBODY_TESTS,
    CONTROL_PREFIX,
    DATE_COLS,
    ID_COLS,
    LIBRARIES_RENAME,
    MEDICATION_ALIASES,
    REFERENCE_DATE,
)
from .raw import load_raw

REFERENCE = pd.Timestamp(REFERENCE_DATE)

LONG_TABLES = ("vitals", "lab_report", "pft", "mrss", "medications",
               "antibodies", "bal", "skin_biopsies", "libraries")

# registry fields that describe the disease itself. Two registrations of one person
# that disagree on any of these need a human decision, not a row count.
DISEASE_FIELDS = ("ssc_subtype", "other_dx", "diagnosis", "diagnosis_date",
                  "nonraynaud_sx", "raynaud_date", "nonraynaud_date")

# Fixed decisions, numbered as in REPORT.md (section 3)
HEIGHT_CM_THRESHOLD = 100            # issue 1: a height below this is inches, at or above it is cm
BMI_PLAUSIBLE = (12, 70)             # issue 2: generous bounds, so genuinely obese patients are not "cleaned" away
KG_IF_RECORDED_BMI_BELOW = 15        # issue 11: a recorded BMI below this only happens when the weights are already kg
DEFAULT_WEIGHT_LB = 160.6            # issue 10: the template default that replaces whole weight histories
BP_DEFAULTS = {"BP SYSTOLIC": 124.0, "BP DIASTOLIC": 77.0, "PULSE": 76.0}   # issue 12: template defaults
LAB_MISSING_CODES = (999.0, 9999.0)  # issue 15: classic placeholder codes for a missing result
PFT_BOUNDS = (40.0, 130.0)           # issue 26: every PFT value sits inside these; a pile-up at a bound is clipping
MAX_TABLETS_PER_DOSE = 10            # issue 29: impossible doses
MAX_GRAMS_PER_DOSE = 10
MAX_MG_PER_DOSE = 5000
SSC_SPECIFIC_DRUGS = ("mycophenolate mofetil", "rituximab", "tocilizumab")  # issue 30: started only after the diagnosis
RNA_CONCENTRATION_OUTLIER_FACTOR = 10  # issue 39: a concentration this many times the median is a suspected unit slip

LB_PER_KG = 2.20462
CM_PER_INCH = 2.54


@dataclass
class IssueLog:
    """Append-only log of every decision applied by the pipeline."""

    rows: list[dict] = field(default_factory=list)

    def add(self, table: str, issue: str, decision: str, n_affected: int = 0,
            detail: str = "") -> None:
        self.rows.append(dict(table=table, issue=issue, decision=decision,
                              n_affected=int(n_affected), detail=str(detail)))

    def frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            self.rows,
            columns=["table", "issue", "decision", "n_affected", "detail"],
        )


# module state for one pipeline run: the log and the duplicate-registration remap
_ISSUES = IssueLog()
_ID_REMAP: dict[str, str] = {}


# --------------------------------------------------------------- helpers
def subject_ids(table: pd.DataFrame, table_name: str) -> pd.Series:
    """The patient id column of a raw table, as text, under the canonical id.

    The id of a duplicate registration is rewritten to the canonical id here,
    so every step after this one sees one id per person.
    """
    ids = table[ID_COLS[table_name]].astype(str).str.strip()
    is_remapped = ids.isin(_ID_REMAP)
    if is_remapped.any():
        ids = ids.replace(_ID_REMAP)
        _ISSUES.add(table_name, "rows recorded under a duplicate-registration id",
                    "subject_id rewritten to the canonical id", int(is_remapped.sum()))
    return ids.rename("subject_id")


def parse_dates(table: pd.DataFrame, table_name: str,
                drop_source: bool = False) -> pd.DataFrame:
    """Parse every date column listed in DATE_COLS for this table.

    A value that is present but does not parse becomes missing and is logged.
    With drop_source the text column is dropped once it has been parsed, which
    is what the libraries table needs (its date columns are renamed as well).
    """
    parsed = table.copy()
    for source_column, target_column in DATE_COLS[table_name].items():
        dates = pd.to_datetime(parsed[source_column], errors="coerce")
        n_unparseable = int((parsed[source_column].notna() & dates.isna()).sum())
        if n_unparseable:
            _ISSUES.add(table_name, f"unparseable dates in {source_column}",
                        "set to missing", n_unparseable)
        if drop_source and source_column != target_column:
            parsed = parsed.drop(columns=[source_column])
        parsed[target_column] = dates
    return parsed


def open_long_table(table_name: str, table: pd.DataFrame,
                    keep_columns: list[str]) -> pd.DataFrame:
    """First pass shared by every longitudinal table.

    Canonical id, parsed dates, exact duplicate rows dropped, the columns the
    pipeline needs, sorted by patient and date.
    """
    opened = table.copy()
    opened["subject_id"] = subject_ids(opened, table_name)
    opened = parse_dates(opened, table_name)
    n_exact_duplicates = int(opened.duplicated().sum())
    if n_exact_duplicates:
        _ISSUES.add(table_name, "exact duplicate rows", "dropped", n_exact_duplicates)
        opened = opened.drop_duplicates()
    return opened[keep_columns].sort_values(["subject_id", "date"]).reset_index(drop=True)


def flag_before_diagnosis(table: pd.DataFrame, table_name: str, what: str,
                          diagnosis_date_by_patient: pd.Series) -> pd.Series:
    """Rows dated before the patient's diagnosis date (issues 30, 31, 42): flagged, never moved."""
    diagnosis_date = table["subject_id"].map(diagnosis_date_by_patient)
    before = (table["date"] < diagnosis_date).fillna(False)
    if before.any():
        days_before = (diagnosis_date - table["date"]).dt.days[before]
        _ISSUES.add(table_name, f"{what} dated before the diagnosis", "flagged, dates not changed",
                    int(before.sum()),
                    detail=f"{table.loc[before, 'subject_id'].nunique()} patients; "
                           f"{int(days_before.min())} to {int(days_before.max())} days before, "
                           f"median {int(days_before.median())}")
    return before.astype(bool)


def drop_merged_second_rows(table: pd.DataFrame, table_name: str, original_ids: pd.Series) -> pd.DataFrame:
    """After the id remap a merged patient has two registry rows; keep the one that
    originally carried the canonical id (its clinical rows are the majority)."""
    is_canonical_row = original_ids.astype(str).str.strip().eq(table["subject_id"])
    ordered = table.assign(_canonical=is_canonical_row).sort_values("_canonical", ascending=False, kind="stable")
    is_second_row = ordered.duplicated(subset=["subject_id"])
    if is_second_row.any():
        # which fields the two registrations disagree on, so the choice is visible in the log
        disagreements = []
        for patient_id, rows in ordered[ordered["subject_id"].isin(ordered.loc[is_second_row, "subject_id"])].groupby("subject_id"):
            compared = rows.drop(columns=["_canonical", "subject_id"]).astype(str)
            differing = [column for column in compared.columns if compared[column].nunique() > 1]
            note = " needs adjudication" if set(differing) & set(DISEASE_FIELDS) else ""
            disagreements.append(f"{patient_id}: {differing if differing else 'identical'}{note}")
        _ISSUES.add(table_name, "second registry row of a merged patient",
                    "dropped; the row of the canonical id is kept", int(is_second_row.sum()),
                    detail="fields on which the two rows disagree: " + "; ".join(disagreements))
        ordered = ordered[~is_second_row]
    ordered["merged_registration_flag"] = ordered["subject_id"].isin(set(_ID_REMAP.values()))   # kept row of a merged patient
    return ordered.drop(columns="_canonical").sort_index()


# ------------------------------------------------------ duplicate registrations
def find_duplicate_registrations(demographics: pd.DataFrame,
                                 raw_tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Same first name, last name and birth date under different ids (issue 4).

    One id per person is kept as canonical: the one with more rows across the
    longitudinal tables, the smaller number on a tie. The other id is remapped.

    A row with a missing first name, last name or birth date is not matched
    automatically: two people can share the fields that remain, so the evidence
    is not enough to merge on. Those rows are logged for manual review instead.
    """
    identity_key = ["first name", "last name", "birth date"]
    identity_key_incomplete = demographics[identity_key].isna().any(axis=1)
    if identity_key_incomplete.any():
        _ISSUES.add("demographics", "incomplete identity key (first name, last name or birth date missing)",
                    "left unmerged; automatic duplicate matching needs all three fields",
                    int(identity_key_incomplete.sum()),
                    detail="needs adjudication: "
                           + str(demographics.loc[identity_key_incomplete, "case number"].tolist()))
    complete_identity = demographics[~identity_key_incomplete]

    ids_by_person = complete_identity.groupby(identity_key, dropna=False)["case number"].apply(list)
    remap: dict[str, str] = {}
    for ids in ids_by_person[ids_by_person.str.len() > 1]:
        rows_by_id = {
            patient_id: sum(int((raw_tables[name][ID_COLS[name]] == patient_id).sum())
                            for name in LONG_TABLES)
            for patient_id in ids
        }
        canonical = max(ids, key=lambda patient_id: (rows_by_id[patient_id],
                                                     -int(re.sub(r"\D", "", patient_id) or 0)))
        for patient_id in ids:
            if patient_id != canonical:
                remap[patient_id] = canonical
                _ISSUES.add("demographics", "duplicate registration (same name and birth date)",
                            f"remapped {patient_id} to {canonical} (the id with more clinical rows)",
                            1, detail=f"1 registration; longitudinal rows: {rows_by_id}")
    return pd.DataFrame([{"original_id": original, "canonical_id": canonical}
                         for original, canonical in remap.items()],
                        columns=["original_id", "canonical_id"])


# --------------------------------------------------------------- vitals
def clean_vitals(vitals_raw: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """Vitals: the default weight, the kg patients, blood-pressure defaults and inversions.

    Returns the cleaned table and the list of patients whose weights are kg,
    which the demographics rules need.
    """
    vitals = open_long_table("vitals", vitals_raw,
                             ["subject_id", "date", "vital_type_name_category", "vital_value"])
    vitals = vitals.rename(columns={"vital_type_name_category": "measure"})
    vitals["value"] = pd.to_numeric(vitals["vital_value"], errors="coerce")
    n_not_numeric = int((vitals["vital_value"].notna() & vitals["value"].isna()).sum())
    if n_not_numeric:
        _ISSUES.add("vitals", "non-numeric vital_value", "set to missing", n_not_numeric)

    is_weight = vitals["measure"].eq("WEIGHT IN POUND")
    is_bmi = vitals["measure"].eq("BMI")

    # issue 10: the template default replaces whole weight histories
    is_default_weight = is_weight & vitals["value"].eq(DEFAULT_WEIGHT_LB)
    weights_by_patient = (vitals[is_weight].assign(is_default=is_default_weight[is_weight])
                          .groupby("subject_id")["is_default"].agg(n_weights="size", n_default="sum"))
    patients_with_default = weights_by_patient[weights_by_patient["n_default"] > 0]
    n_patients_all_weights_default = int((patients_with_default["n_default"] == patients_with_default["n_weights"]).sum())
    vitals.loc[is_default_weight, "value"] = np.nan
    _ISSUES.add("vitals", f"weight equals the template default {DEFAULT_WEIGHT_LB}",
                "set to missing", int(is_default_weight.sum()),
                detail=f"{int(is_default_weight.sum())} weight rows over {len(patients_with_default)} patients, "
                       f"{n_patients_all_weights_default} of them with no other weight on record")

    # issue 11: a recorded BMI below 15 only happens when the weights are already kg, and the
    # stored BMI was computed from that kg number with the pound formula
    kg_patients = sorted(vitals.loc[is_bmi & (vitals["value"] < KG_IF_RECORDED_BMI_BELOW), "subject_id"].unique())
    weight_is_kg = is_weight & vitals["subject_id"].isin(kg_patients)
    corrupted_bmi = is_bmi & vitals["subject_id"].isin(kg_patients)
    vitals.loc[weight_is_kg, "value"] = vitals.loc[weight_is_kg, "value"] * LB_PER_KG
    vitals.loc[corrupted_bmi, "value"] = np.nan
    _ISSUES.add("vitals", "kg weights recorded in the pound column",
                "converted kg to lb; the derived BMI rows set to missing (the patient-level BMI is recomputed in demographics from the standardized weight)",
                int(weight_is_kg.sum()),
                detail=f"{len(kg_patients)} patients (recorded BMI < {KG_IF_RECORDED_BMI_BELOW}): {kg_patients}")

    # issue 12: blood pressure and pulse sit on template defaults and barely vary; reported, not changed
    by_visit = vitals.pivot_table(index=["subject_id", "date"], columns="measure",
                                  values="value", aggfunc="first")
    default_shares = {measure: float(vitals.loc[vitals["measure"].eq(measure), "value"].eq(default).mean())
                      for measure, default in BP_DEFAULTS.items()}
    systolic_by_patient = by_visit.groupby(level="subject_id")["BP SYSTOLIC"].agg(n="count", distinct="nunique")
    followed = systolic_by_patient[systolic_by_patient["n"] >= 4]
    n_near_constant = int((followed["distinct"] <= 2).sum())
    n_default_rows = int(sum(vitals.loc[vitals["measure"].eq(measure), "value"].eq(default).sum()
                             for measure, default in BP_DEFAULTS.items()))
    _ISSUES.add("vitals", "blood pressure and pulse on template defaults (124 / 77 / 76) and near-constant series",
                "reported only; not usable as a longitudinal signal, nothing changed", n_default_rows,
                detail=", ".join(f"{measure} = {default:.0f} in {share:.0%}" for (measure, default), share
                                 in zip(BP_DEFAULTS.items(), default_shares.values()))
                       + f"; {n_near_constant} of {len(followed)} patients with >= 4 systolic readings "
                         f"have <= 2 distinct values")

    # issue 13: systolic below diastolic on the same visit
    inverted_visits = by_visit.index[by_visit["BP SYSTOLIC"] < by_visit["BP DIASTOLIC"]]
    visit_key = pd.MultiIndex.from_frame(vitals[["subject_id", "date"]])
    vitals["bp_inverted_flag"] = (visit_key.isin(inverted_visits)
                                  & vitals["measure"].isin(["BP SYSTOLIC", "BP DIASTOLIC"]))
    n_inverted_rows = int(vitals["bp_inverted_flag"].sum())
    _ISSUES.add("vitals", "systolic below diastolic on the same visit", "flagged, possible field swap",
                n_inverted_rows,
                detail=f"{len(inverted_visits)} visits, {n_inverted_rows} BP rows; "
                       f"patients: {sorted(set(inverted_visits.get_level_values('subject_id')))}")

    return vitals[["subject_id", "date", "measure", "value", "bp_inverted_flag"]], kg_patients


# ---------------------------------------------------------- demographics
def clean_demographics(demographics_raw: pd.DataFrame, kg_patients: list[str]) -> pd.DataFrame:
    """Demographics: canonical id, height to cm, weight to kg, BMI and age."""
    demographics = parse_dates(demographics_raw, "demographics")          # adds birth_date
    demographics["subject_id"] = subject_ids(demographics, "demographics")
    demographics = drop_merged_second_rows(demographics, "demographics", demographics["case number"])
    demographics = demographics.drop(columns=["birth date", "case number"])

    # issue 1: height in two units
    height = pd.to_numeric(demographics["height"], errors="coerce")
    height_not_numeric = int((demographics["height"].notna() & height.isna()).sum())
    if height_not_numeric:
        _ISSUES.add("demographics", "non-numeric height", "set to missing", height_not_numeric)
    height_is_inches = height < HEIGHT_CM_THRESHOLD
    demographics["height_cm"] = np.where(height_is_inches, height * CM_PER_INCH, height)
    _ISSUES.add("demographics", "height in two units (inches and cm)",
                f"standardized to cm; values below {HEIGHT_CM_THRESHOLD} read as inches",
                int(height_is_inches.sum()))

    # issues 2, 3, 11: weight is pounds unless (a) the patient is on the kg list from
    # vitals or (b) the pound reading gives an impossible BMI while kg is plausible
    weight = pd.to_numeric(demographics["weight"], errors="coerce")
    weight_not_numeric = int((demographics["weight"].notna() & weight.isna()).sum())
    if weight_not_numeric:
        _ISSUES.add("demographics", "non-numeric weight", "set to missing", weight_not_numeric)
    height_m2 = (demographics["height_cm"] / 100) ** 2
    bmi_if_lb = (weight / LB_PER_KG) / height_m2
    bmi_if_kg = weight / height_m2
    lb_plausible = bmi_if_lb.between(*BMI_PLAUSIBLE)
    kg_plausible = bmi_if_kg.between(*BMI_PLAUSIBLE)
    on_kg_list = demographics["subject_id"].isin(kg_patients)
    has_weight = weight.notna()
    weight_is_kg = on_kg_list | (~lb_plausible & kg_plausible)
    weight_impossible = has_weight & ~weight_is_kg & ~lb_plausible   # a missing weight is not impossible
    demographics["weight_kg"] = np.select(
        [weight_is_kg, ~weight_is_kg & lb_plausible],
        [weight, weight / LB_PER_KG],
        default=np.nan,
    )
    _ISSUES.add("demographics", "kg weights among pounds",
                "standardized to kg (vitals cross-check first, BMI plausibility second)",
                int(weight_is_kg.sum()),
                detail=f"{int(on_kg_list.sum())} from the vitals kg list, "
                       f"{int((weight_is_kg & ~on_kg_list).sum())} from the BMI rule alone")
    if weight_impossible.any():
        _ISSUES.add("demographics", "weight impossible in both lb and kg", "set to missing",
                    int(weight_impossible.sum()),
                    detail=str(demographics.loc[weight_impossible, "subject_id"].tolist()))

    demographics["bmi_calc"] = demographics["weight_kg"] / height_m2
    demographics["age_years"] = (REFERENCE - demographics["birth_date"]).dt.days / 365.25

    keep = ["subject_id", "first name", "last name", "birth_date", "age_years",
            "ethnicity", "gender", "races", "diagnosis", "state",
            "height_cm", "weight_kg", "bmi_calc", "merged_registration_flag"]
    return demographics[keep].sort_values("subject_id").reset_index(drop=True)


# ----------------------------------------------------------- ssc_subtype
def clean_ssc_subtype(ssc_subtype_raw: pd.DataFrame) -> pd.DataFrame:
    """Subtype table: milestone dates, onset-order flags, three-state comorbidity flags."""
    ssc_subtype = parse_dates(ssc_subtype_raw, "ssc_subtype")
    ssc_subtype["subject_id"] = subject_ids(ssc_subtype, "ssc_subtype")
    ssc_subtype = drop_merged_second_rows(ssc_subtype, "ssc_subtype", ssc_subtype["study_code"])
    ssc_subtype = ssc_subtype.drop(columns=["study_code"])

    # issue 7: Raynaud's dated after the first non-Raynaud symptom, at the same rate in both
    # subtypes: an unconstrained generator draw, flagged and never swapped
    # a patient missing either date has no onset order to compare, so the flag is
    # "unknown" there rather than "ok"
    onset_dates_known = ssc_subtype["raynaud_date"].notna() & ssc_subtype["nonraynaud_date"].notna()
    raynaud_after = onset_dates_known & (ssc_subtype["raynaud_date"] > ssc_subtype["nonraynaud_date"])
    ssc_subtype["onset_order_flag"] = np.select(
        [~onset_dates_known, raynaud_after],
        ["unknown", "raynaud_after_nonraynaud"],
        default="ok",
    )
    _ISSUES.add("ssc_subtype", "raynaud_date after nonraynaud_date",
                "flagged, dates not changed", int(raynaud_after.sum()),
                detail="rate by subtype: " + str(raynaud_after.groupby(ssc_subtype["ssc_subtype"]).mean().round(3).to_dict()))
    if (~onset_dates_known).any():
        _ISSUES.add("ssc_subtype", "raynaud_date or nonraynaud_date missing",
                    "onset_order_flag set to unknown", int((~onset_dates_known).sum()))

    # issue 8: diagnosed before the first non-Raynaud symptom (possible; flagged only)
    ssc_subtype["diagnosis_before_symptom_flag"] = (
        ssc_subtype["diagnosis_date"] < ssc_subtype["nonraynaud_date"]).fillna(False).astype(bool)
    _ISSUES.add("ssc_subtype", "diagnosis_date before nonraynaud_date",
                "flagged only (early diagnosis on Raynaud's plus antibodies is possible)",
                int(ssc_subtype["diagnosis_before_symptom_flag"].sum()))

    # issue 6: comorbidities come from a semicolon-joined free-text field that is missing
    # for 37% of patients; missing means "not recorded", so the flags stay unknown there
    other_dx = ssc_subtype["other_dx"]
    for condition in ("ILD", "GERD", "PAH"):
        recorded = other_dx.str.contains(condition, regex=False).astype("boolean")
        ssc_subtype[f"dx_{condition.lower()}"] = recorded.mask(other_dx.isna(), pd.NA)   # unknown, not False
    _ISSUES.add("ssc_subtype", "other_dx missing", "comorbidity flags left unknown (not false) where other_dx is missing",
                int(other_dx.isna().sum()),
                detail=f"{other_dx.isna().mean():.0%} of patients")

    ssc_subtype["disease_duration_years"] = (REFERENCE - ssc_subtype["diagnosis_date"]).dt.days / 365.25
    ssc_subtype["onset_to_dx_years"] = (ssc_subtype["diagnosis_date"] - ssc_subtype["nonraynaud_date"]).dt.days / 365.25

    keep = ["subject_id", "ssc_subtype", "other_dx",
            "raynaud_date", "nonraynaud_date", "diagnosis_date",
            "onset_order_flag", "diagnosis_before_symptom_flag",
            "dx_ild", "dx_gerd", "dx_pah",
            "disease_duration_years", "onset_to_dx_years", "merged_registration_flag"]
    return ssc_subtype[keep].sort_values("subject_id").reset_index(drop=True)


# ------------------------------------------------------------ lab_report
def clean_labs(lab_report_raw: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Lab results: the method row split out, sentinel and impossible values to missing.

    Returns the numeric results and the separate differential-method table.
    """
    lab_report = open_long_table("lab_report", lab_report_raw,
                                 ["subject_id", "date", "component_name", "value"])
    lab_report = lab_report.rename(columns={"component_name": "component"})

    # issue 14: DIFFERENTIAL TYPE carries the text "Automated", the method of the count, not a result
    is_method_row = lab_report["component"].eq("DIFFERENTIAL TYPE")
    differential_method = (lab_report[is_method_row]
                           .rename(columns={"value": "differential_method"})
                           [["subject_id", "date", "differential_method"]])
    lab_report = lab_report[~is_method_row].copy()
    _ISSUES.add("lab_report", "DIFFERENTIAL TYPE is the method of the count, not a measurement",
                "split into its own table", int(is_method_row.sum()))

    lab_report["value_num"] = pd.to_numeric(lab_report["value"], errors="coerce")
    n_not_numeric = int((lab_report["value"].notna() & lab_report["value_num"].isna()).sum())
    if n_not_numeric:
        _ISSUES.add("lab_report", "non-numeric lab value", "set to missing", n_not_numeric)

    # issue 15: placeholder codes and impossible values. A red-cell count of 0 is impossible;
    # a basophil or eosinophil count of 0 is legitimate and stays.
    is_missing_code = lab_report["value_num"].isin(LAB_MISSING_CODES)
    is_negative = lab_report["value_num"] < 0
    is_rbc_zero = lab_report["component"].eq("RBC") & lab_report["value_num"].eq(0)
    is_sentinel = is_missing_code | is_negative | is_rbc_zero
    detail = (lab_report.loc[is_sentinel].groupby(["component", "value_num"]).size().to_dict())
    lab_report.loc[is_sentinel, "value_num"] = np.nan
    _ISSUES.add("lab_report", "sentinel and impossible values (999, 9999, negatives, RBC = 0)",
                "set to missing", int(is_sentinel.sum()), detail=str(detail))

    return (lab_report[["subject_id", "date", "component", "value_num"]],
            differential_method)


# ------------------------------------------------------------------- pft
def clean_pft(pft_raw: pd.DataFrame) -> pd.DataFrame:
    """Lung function: DLCO missingness logged, FEV1 flagged, clipping logged."""
    pft = open_long_table("pft", pft_raw, ["subject_id", "date", "NAME", "ORD_VALUE"])
    pft = pft.rename(columns={"NAME": "measure", "ORD_VALUE": "value"})
    pft_value_num = pd.to_numeric(pft["value"], errors="coerce")
    n_not_numeric = int((pft["value"].notna() & pft_value_num.isna()).sum())
    if n_not_numeric:
        _ISSUES.add("pft", "non-numeric PFT value", "set to missing", n_not_numeric)
    pft["value"] = pft_value_num

    # issue 24: DLCO is measured less often; genuine absence, no imputation
    is_dlco = pft["measure"].eq("DLCO_SB")
    _ISSUES.add("pft", "DLCO_SB missing", "kept missing, no imputation",
                int(pft.loc[is_dlco, "value"].isna().sum()),
                detail=f"{pft.loc[is_dlco, 'value'].isna().mean():.0%} of DLCO rows")

    # issue 25: in this dataset FEV1 tracks FVC and adds no usable variation. In % predicted,
    # FEV1 above FVC is not impossible; what stands out here is a ratio sitting at 1.00 with
    # no obstructive tail and no association with subtype or with recorded fibrosis. The pairs
    # are flagged and FEV1 is left out of the analyses, which is a statement about this export.
    by_visit = pft.pivot_table(index=["subject_id", "date"], columns="measure",
                               values="value", aggfunc="first")
    if {"FEV1", "FVC"}.issubset(by_visit.columns):
        fev1_above_fvc_visits = by_visit.index[by_visit["FEV1"] > by_visit["FVC"]]
        visit_key = pd.MultiIndex.from_frame(pft[["subject_id", "date"]])
        pft["fev1_gt_fvc_flag"] = visit_key.isin(fev1_above_fvc_visits)
        n_flagged_rows = int(pft["fev1_gt_fvc_flag"].sum())
        ratio = by_visit["FEV1"] / by_visit["FVC"]
        _ISSUES.add("pft", "FEV1 tracks FVC and adds no usable variation in this dataset",
                    "kept + flagged; FEV1 and FEV1/FVC not used in analyses",
                    n_flagged_rows,
                    detail=f"{len(fev1_above_fvc_visits)} visits, {n_flagged_rows} rows (every PFT row of "
                           f"those visits); ratio median {ratio.median():.3f}, "
                           f"share below 0.90 {(ratio < 0.9).mean():.1%}")

    # issue 26: values clipped at exactly 40 and 130
    at_bound = pft["value"].isin(PFT_BOUNDS)
    _ISSUES.add("pft", f"values at the generator bounds {PFT_BOUNDS}", "kept, documented",
                int(at_bound.sum()),
                detail=str(pft.loc[at_bound].groupby(["measure", "value"]).size().to_dict()))
    return pft


# ------------------------------------------------------------------ mrss
def clean_mrss(mrss_raw: pd.DataFrame) -> pd.DataFrame:
    """Skin score: numeric score, rater, same-day conflicts flagged."""
    mrss = open_long_table("mrss", mrss_raw, ["subject_id", "date", "ENTRY_USER_NAME", "mrss_score"])
    mrss = mrss.rename(columns={"ENTRY_USER_NAME": "rater"})
    mrss["score"] = pd.to_numeric(mrss["mrss_score"], errors="coerce")
    n_not_numeric = int((mrss["mrss_score"].notna() & mrss["score"].isna()).sum())
    if n_not_numeric:
        _ISSUES.add("mrss", "non-numeric mrss_score", "set to missing", n_not_numeric)

    # issue 21: two scores for one patient on one day cannot both be right
    mrss["same_day_conflict"] = mrss.duplicated(subset=["subject_id", "date"], keep=False)
    if mrss["same_day_conflict"].any():
        _ISSUES.add("mrss", "two scores on the same visit date", "kept + flagged",
                    int(mrss["same_day_conflict"].sum()),
                    detail=str(mrss.loc[mrss["same_day_conflict"], "subject_id"].unique().tolist()))
    return mrss[["subject_id", "date", "rater", "score", "same_day_conflict"]]


# ----------------------------------------------------------- medications
DOSE_PATTERN = re.compile(r"^\s*(-?\d+(?:\.\d+)?)\s*([A-Za-z]+)")


def dose_is_impossible(dose) -> bool:
    """Issue 29: a dose that cannot be real (zero, negative, absurd mass, mixed units).

    Doses without a leading number and unit (for example "PRN") are left alone.
    """
    if not isinstance(dose, str):
        return False
    match = DOSE_PATTERN.match(dose)
    if match is None:
        return False
    amount, unit = float(match.group(1)), match.group(2).lower()
    if amount <= 0:
        return True
    if unit == "g" and amount >= MAX_GRAMS_PER_DOSE:
        return True
    if unit == "mg" and amount >= MAX_MG_PER_DOSE:
        return True
    if unit.startswith("tablet") and amount >= MAX_TABLETS_PER_DOSE:
        return True
    if "every hour" in dose.lower():
        return True
    if unit == "ml" and "tablet" in dose.lower():
        return True
    return False


def clean_medications(medications_raw: pd.DataFrame,
                      diagnosis_date_by_patient: pd.Series) -> pd.DataFrame:
    """Prescriptions: one name per drug, impossible doses to missing, timeline flags."""
    medications = open_long_table("medications", medications_raw,
                                  ["subject_id", "date", "medication", "dose"])

    # issue 27: one drug under three names
    medications["medication_generic"] = medications["medication"].replace(MEDICATION_ALIASES).str.strip()
    n_renamed = int((medications["medication"].ne(medications["medication_generic"])
                     & medications["medication"].notna()).sum())
    _ISSUES.add("medications", "brand name and abbreviation variants of one drug",
                "mapped to the generic name", n_renamed, detail=str(MEDICATION_ALIASES))

    # issue 28: rows without a drug name or a dose are kept, out of drug-level counts
    n_no_name = int(medications["medication"].isna().sum())
    n_no_dose = int(medications["dose"].isna().sum())
    _ISSUES.add("medications", "rows without a drug name", "kept, excluded from drug counts", n_no_name)
    _ISSUES.add("medications", "rows without a dose", "kept as missing", n_no_dose)

    # issue 29: impossible doses
    impossible = medications["dose"].map(dose_is_impossible).astype(bool)
    detail = medications.loc[impossible, "dose"].value_counts().to_dict()
    medications.loc[impossible, "dose"] = pd.NA
    _ISSUES.add("medications", "impossible doses (zero, negative, absurd mass, mixed units)",
                "dose set to missing", int(impossible.sum()), detail=str(detail))

    # issue 30: general drugs precede the diagnosis routinely (nifedipine for Raynaud's);
    # SSc-specific immunosuppressants should not
    medications["before_diagnosis_flag"] = flag_before_diagnosis(
        medications, "medications", "prescriptions", diagnosis_date_by_patient)
    medications["ssc_specific_before_diagnosis_flag"] = (
        medications["before_diagnosis_flag"] & medications["medication_generic"].isin(SSC_SPECIFIC_DRUGS))
    _ISSUES.add("medications", "SSc-specific immunosuppressants prescribed before the diagnosis",
                "flagged only (likely the diagnosis date was entered after treatment started)",
                int(medications["ssc_specific_before_diagnosis_flag"].sum()))

    return medications[["subject_id", "date", "medication_generic", "dose",
                        "before_diagnosis_flag", "ssc_specific_before_diagnosis_flag"]]


# ------------------------------------------------------------ antibodies
def clean_antibodies(antibodies_raw: pd.DataFrame) -> pd.DataFrame:
    """Antibody results: short test names, flips between positive and negative flagged."""
    antibodies = open_long_table("antibodies", antibodies_raw, ["subject_id", "date", "test", "value"])
    antibodies["test_short"] = antibodies["test"].map(ANTIBODY_TESTS)
    unmapped = antibodies["test"].notna() & antibodies["test_short"].isna()
    if unmapped.any():
        _ISSUES.add("antibodies", "test name not listed in ANTIBODY_TESTS",
                    "short name left missing; add the name to the config map", int(unmapped.sum()),
                    detail=str(sorted(antibodies.loc[unmapped, "test"].unique())))

    # issue 18: these antibodies are lifelong markers; a patient switching between positive
    # and negative on repeat testing is flagged, both rows kept, the latest result used downstream
    definite = antibodies[antibodies["value"].isin(["positive", "negative"])]
    n_results = definite.groupby(["subject_id", "test_short"])["value"].nunique()
    flipping_pairs = n_results[n_results > 1].index
    pair_key = pd.MultiIndex.from_frame(antibodies[["subject_id", "test_short"]])
    antibodies["result_flip_flag"] = pair_key.isin(flipping_pairs)
    n_flagged_rows = int(antibodies["result_flip_flag"].sum())
    _ISSUES.add("antibodies", "positive and negative results for the same patient and test",
                "kept + flagged; analyses use the latest result per test",
                n_flagged_rows,
                detail=f"{len(flipping_pairs)} patient-test pairs, {n_flagged_rows} rows, "
                       f"{flipping_pairs.get_level_values('subject_id').nunique()} patients "
                       f"(counted after duplicate registrations are merged)")

    return antibodies[["subject_id", "date", "test_short", "value", "result_flip_flag"]]


# ------------------------------------------------------------------- bal
def clean_bal(bal_raw: pd.DataFrame, diagnosis_date_by_patient: pd.Series) -> pd.DataFrame:
    """Lavage procedures: the volume rule checked, procedures before the diagnosis flagged."""
    bal = open_long_table("bal", bal_raw, ["subject_id", "date", "procedure_site",
                                           "volume_instilled_ml", "volume_recovered_ml", "bal_comment"])
    for column in ["volume_instilled_ml", "volume_recovered_ml"]:
        as_number = pd.to_numeric(bal[column], errors="coerce")
        n_not_numeric = int((bal[column].notna() & as_number.isna()).sum())
        if n_not_numeric:
            _ISSUES.add("bal", f"non-numeric {column}", "set to missing", n_not_numeric)
        bal[column] = as_number
    recovered_more = bal["volume_recovered_ml"] > bal["volume_instilled_ml"]
    bal["volume_flag"] = np.where(recovered_more, "recovered>instilled", "ok")
    _ISSUES.add("bal", "recovered volume above instilled volume", "checked; kept + flagged if any",
                int(recovered_more.sum()))
    bal["before_diagnosis_flag"] = flag_before_diagnosis(bal, "bal", "lavage procedures",
                                                         diagnosis_date_by_patient)
    return bal


# --------------------------------------------------------- skin_biopsies
def clean_biopsies(skin_biopsies_raw: pd.DataFrame,
                   diagnosis_date_by_patient: pd.Series) -> pd.DataFrame:
    """Biopsies: file extension against the declared format, biopsies before the diagnosis flagged."""
    skin_biopsies = open_long_table("skin_biopsies", skin_biopsies_raw,
                                    ["subject_id", "date", "ENTRY_USER_NAME", "biopsy_site",
                                     "clinical_indication", "specimen_accession",
                                     "image_file_path", "image_format"])
    skin_biopsies = skin_biopsies.rename(columns={"ENTRY_USER_NAME": "clinician"})

    # issue 41: the extension in the path disagrees with image_format in 61% of rows;
    # which field is right needs the files, so both are kept and the row is flagged
    extension = (skin_biopsies["image_file_path"].str.extract(r"\.(\w+)$")[0]
                 .str.upper().replace({"TIF": "TIFF"}))
    skin_biopsies["path_extension"] = extension
    skin_biopsies["format_mismatch"] = extension.ne(skin_biopsies["image_format"].str.upper())
    _ISSUES.add("skin_biopsies", "file extension disagrees with image_format",
                "both fields kept, flagged for pathology adjudication",
                int(skin_biopsies["format_mismatch"].sum()),
                detail="extensions: " + str(extension.value_counts().to_dict()))

    skin_biopsies["before_diagnosis_flag"] = flag_before_diagnosis(
        skin_biopsies, "skin_biopsies", "biopsies", diagnosis_date_by_patient)
    return skin_biopsies


# ------------------------------------------------------------- libraries
def clean_libraries(libraries_raw: pd.DataFrame) -> pd.DataFrame:
    """RNA-seq library log: snake_case columns, numeric QC fields, status and identifier flags."""
    libraries = libraries_raw.copy()
    libraries["subject_id"] = subject_ids(libraries, "libraries")
    libraries = parse_dates(libraries, "libraries", drop_source=True)
    libraries = libraries.drop(columns=["reg_id"]).rename(columns=LIBRARIES_RENAME)   # issue 34
    n_exact_duplicates = int(libraries.duplicated().sum())
    if n_exact_duplicates:
        _ISSUES.add("libraries", "exact duplicate rows", "dropped", n_exact_duplicates)
        libraries = libraries.drop_duplicates()
    for column in ["cell_viability_pct", "rin", "rna_concentration_pg_ul", "elution_vol_ul", "ul_for_250pg"]:
        as_number = pd.to_numeric(libraries[column], errors="coerce")
        n_not_numeric = int((libraries[column].notna() & as_number.isna()).sum())
        if n_not_numeric:
            _ISSUES.add("libraries", f"non-numeric {column}", "set to missing", n_not_numeric)
        libraries[column] = as_number

    # issue 35: marked "not for sequencing" but with library work done
    libraries["status_contradiction"] = (libraries["sequence_flag"].eq("No")
                                         & libraries["complete_flag"].isin(["Complete", "Partial"]))
    _ISSUES.add("libraries", "Sequence? = No but library work recorded (Complete or Partial)",
                "both flags kept, contradiction flagged", int(libraries["status_contradiction"].sum()))

    # issue 36: sequenced samples without a RIN although the QC was run
    is_sequenced = libraries["sequence_flag"].eq("Yes")
    libraries["rin_missing_for_sequenced"] = is_sequenced & libraries["rin"].isna()
    _ISSUES.add("libraries", "sequenced samples without a RIN", "kept missing; flagged; a question for the wet lab",
                int(libraries["rin_missing_for_sequenced"].sum()),
                detail=f"of {int(is_sequenced.sum())} sequenced samples; TapeStation comments: "
                       + str(libraries.loc[libraries["rin_missing_for_sequenced"], "ts_comment"]
                             .fillna("(empty)").value_counts().to_dict()))

    # issue 37: a sample id used twice
    has_sample_id = libraries["sample_id"].notna()
    libraries["sample_id_reused"] = has_sample_id & libraries["sample_id"].duplicated(keep=False)
    _ISSUES.add("libraries", "sample id used more than once", "flagged for lab adjudication",
                int(libraries["sample_id_reused"].sum()),
                detail=str(sorted(libraries.loc[libraries["sample_id_reused"], "sample_id"].unique())))

    # issue 38: a sequencing batch of one sample cannot be batch-corrected
    batch_sizes = libraries["rnaseq_batch"].value_counts()
    libraries["single_sample_batch"] = libraries["rnaseq_batch"].map(batch_sizes).eq(1).fillna(False).astype(bool)
    _ISSUES.add("libraries", "sequencing batches with a single sample", "flagged",
                int(libraries["single_sample_batch"].sum()),
                detail=str(sorted(batch_sizes[batch_sizes == 1].index)))

    # issue 39: an RNA concentration far above the rest is a suspected unit slip
    concentration = libraries["rna_concentration_pg_ul"]
    libraries["rna_concentration_outlier"] = (concentration
                                              > RNA_CONCENTRATION_OUTLIER_FACTOR * concentration.median()).fillna(False).astype(bool)
    _ISSUES.add("libraries", "RNA concentration far above the median (suspected unit slip)",
                "kept + flagged", int(libraries["rna_concentration_outlier"].sum()),
                detail=f"median {concentration.median():.0f} pg/ul, max {concentration.max():.0f}")

    columns = ["subject_id", "processing_date", "cell_viability_pct",
               "rna_isolation_kit", "kit_lot_number", "elution_vol_ul",
               "macrophage_rna_tube_id", "macrophage_rna_tube_box_id",
               "rna_qc_date", "rna_volume_qc_ul", "rin",
               "rna_concentration_pg_ul", "comment", "ts_comment", "techcore_comment",
               "sequence_flag", "tapestation_assay_type", "ul_for_250pg",
               "complete_flag", "library_prep_date", "sample_id", "rna_tube_id",
               "library_prep_plate", "rnaseq_batch",
               "status_contradiction", "rin_missing_for_sequenced", "sample_id_reused",
               "single_sample_batch", "rna_concentration_outlier"]
    return libraries[columns].sort_values(["subject_id", "processing_date"]).reset_index(drop=True)


# -------------------------------------------------------------- pipeline
def quarantine_controls(table: pd.DataFrame, table_name: str,
                        registered_ids: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Issue 9: control subjects from a companion study, and any row with no registry match,
    move to a separate frame. Nothing is deleted.

    build_processed() runs this on every longitudinal table, so the same rule decides
    every table; today only vitals, pft, mrss and libraries have rows to quarantine."""
    is_control = table["subject_id"].str.startswith(CONTROL_PREFIX)
    is_orphan = ~is_control & ~table["subject_id"].isin(registered_ids)
    if is_control.any():
        _ISSUES.add(table_name, f"control subjects ({CONTROL_PREFIX}*) in the export",
                    "quarantined to a controls frame", int(is_control.sum()),
                    detail=str(sorted(table.loc[is_control, "subject_id"].unique())))
    if is_orphan.any():
        _ISSUES.add(table_name, "rows with no registry match", "quarantined to the controls frame",
                    int(is_orphan.sum()))
    # both frames are reindexed so the in-memory result matches what comes back
    # from the parquet files, which are written with index=False
    keep = table[~(is_control | is_orphan)].reset_index(drop=True)
    quarantined = table[is_control | is_orphan].reset_index(drop=True)
    return keep, quarantined


def build_processed() -> dict[str, pd.DataFrame]:
    """Run the full cleaning pipeline. Returns the cleaned frames keyed by name."""
    global _ISSUES, _ID_REMAP
    _ISSUES = IssueLog()
    raw_tables = load_raw()

    # 1. identity: duplicate registrations become one canonical id per person; every
    #    subject_ids() call after this line applies the remap
    id_map = find_duplicate_registrations(raw_tables["demographics"], raw_tables)
    _ID_REMAP = dict(zip(id_map["original_id"], id_map["canonical_id"]))

    # 2. registry tables. vitals first: its kg list feeds the demographics weight rule
    vitals, kg_patients = clean_vitals(raw_tables["vitals"])
    demographics = clean_demographics(raw_tables["demographics"], kg_patients)
    ssc_subtype = clean_ssc_subtype(raw_tables["ssc_subtype"])
    registered_ids = set(demographics["subject_id"])
    subtype_ids = set(ssc_subtype["subject_id"])
    if registered_ids != subtype_ids:
        # the two registry tables define the cohort together, so a mismatch is a linkage
        # failure to resolve at the source rather than something to quietly narrow down
        raise ValueError(
            "demographics and ssc_subtype cover different patients. "
            f"In demographics only ({len(registered_ids - subtype_ids)}): "
            f"{sorted(registered_ids - subtype_ids)}. "
            f"In ssc_subtype only ({len(subtype_ids - registered_ids)}): "
            f"{sorted(subtype_ids - registered_ids)}."
        )
    diagnosis_date_by_patient = ssc_subtype.set_index("subject_id")["diagnosis_date"]

    # 3. the longitudinal tables
    lab_report, differential_method = clean_labs(raw_tables["lab_report"])
    pft = clean_pft(raw_tables["pft"])
    mrss = clean_mrss(raw_tables["mrss"])
    medications = clean_medications(raw_tables["medications"], diagnosis_date_by_patient)
    antibodies = clean_antibodies(raw_tables["antibodies"])
    bal = clean_bal(raw_tables["bal"], diagnosis_date_by_patient)
    skin_biopsies = clean_biopsies(raw_tables["skin_biopsies"], diagnosis_date_by_patient)
    libraries = clean_libraries(raw_tables["libraries"])

    processed = {
        "subjects": demographics,
        "ssc_subtype": ssc_subtype,
        "vitals": vitals,
        "labs": lab_report,
        "lab_differential_type": differential_method,
        "pft": pft,
        "mrss": mrss,
        "medications": medications,
        "antibodies": antibodies,
        "bal": bal,
        "biopsies": skin_biopsies,
        "libraries": libraries,
        "subject_id_map": id_map,
    }

    # 4. every longitudinal table goes through the same quarantine rule; a controls frame
    #    is written only where something was quarantined (four tables carry control rows)
    for name in ("vitals", "labs", "pft", "mrss", "medications", "antibodies",
                 "bal", "biopsies", "libraries"):
        processed[name], quarantined = quarantine_controls(processed[name], name, registered_ids)
        if len(quarantined):
            processed[f"controls_{name}"] = quarantined

    # 5. the log is materialized last so every step above is in it
    processed["issues"] = _ISSUES.frame()
    return processed
