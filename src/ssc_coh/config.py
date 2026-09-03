"""Central configuration: paths, schema maps, and fixed decisions.

The maps from source conventions to canonical names live here, apart from
the cleaning rules that use them.
"""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "data" / "raw"
PROC_DIR = ROOT / "data" / "processed"
APP_DIR = ROOT / "app"

# Reference date for age calculations (audit date; deterministic).
REFERENCE_DATE = "2026-08-30"

TABLES = [
    "antibodies",
    "bal",
    "demographics",
    "lab_report",
    "libraries",
    "medications",
    "mrss",
    "pft",
    "skin_biopsies",
    "ssc_subtype",
    "vitals",
]

# canonical subject-id column per source table (4 different conventions)
ID_COLS = {
    "antibodies": "case_number",
    "bal": "reg_id",
    "demographics": "case number",
    "lab_report": "reg_id",
    "libraries": "reg_id",
    "medications": "reg_id",
    "mrss": "reg_id",
    "pft": "case_number",
    "skin_biopsies": "case_number",
    "ssc_subtype": "study_code",
    "vitals": "reg_id",
}

# observation-date columns per table, mapped to a canonical name per table
DATE_COLS = {
    "antibodies": {"dts": "date"},
    "bal": {"procedure_date": "date"},
    "demographics": {"birth date": "birth_date"},
    "lab_report": {"order_date": "date"},
    "libraries": {
        "processing_date": "processing_date",
        "Date of RNA QC": "rna_qc_date",
        "Library prep date": "library_prep_date",
    },
    "medications": {"date": "date"},
    "mrss": {"date": "date"},
    "pft": {"PFT_dts": "date"},
    "skin_biopsies": {"biopsy_date": "date"},
    "ssc_subtype": {
        "raynaud_date": "raynaud_date",
        "nonraynaud_date": "nonraynaud_date",
        "diagnosis_date": "diagnosis_date",
    },
    "vitals": {"date": "date"},
}

# healthy-control subjects from a companion study that leaked into the export
CONTROL_PREFIX = "SSC_NORM_"

# same drug under three names -> canonical generic name
MEDICATION_ALIASES = {
    "CellCept": "mycophenolate mofetil",
    "MMF": "mycophenolate mofetil",
}

# LIMS column renames (spaces / ? / units in names)
LIBRARIES_RENAME = {
    "cell_viability": "cell_viability_pct",
    "processing_date": "processing_date",
    "comment": "comment",
    "RNA isolation kit": "rna_isolation_kit",
    "Kit lot number": "kit_lot_number",
    "Elution vol (ul)": "elution_vol_ul",
    "Macrophage RNA Tube ID": "macrophage_rna_tube_id",
    "Macrophage RNA Tube Location Box ID": "macrophage_rna_tube_box_id",
    "Date of RNA QC": "rna_qc_date",
    "RNA volume for QC": "rna_volume_qc_ul",
    "RIN": "rin",
    "RNA concentration (pg/ul)": "rna_concentration_pg_ul",
    "TS comment": "ts_comment",
    "TechCore comment": "techcore_comment",
    "Sequence?": "sequence_flag",
    "TapeStation assay type": "tapestation_assay_type",
    "ul for 250 pg": "ul_for_250pg",
    "Complete?": "complete_flag",
    "Library prep date": "library_prep_date",
    "Sample": "sample_id",
    "RNA Tube ID": "rna_tube_id",
    "Library prep plate": "library_prep_plate",
    "RNAseq batch": "rnaseq_batch",
}

ANTIBODY_TESTS = {
    "scl70": "scl70",
    "anti-centromere antibodies": "aca",
    "rna polymerase III": "rna_pol3",
}
