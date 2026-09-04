# SSc cohort explorer

An audit, cleaning pipeline and interactive browser for a fully synthetic systemic sclerosis
registry of 11 CSV source tables, covering clinical follow-up and the research sample chain from
bronchoalveolar lavage to RNA-seq libraries. No PHI anywhere in the project.

Live app: https://ssc-cohort-explorer-vrvc7q36vhdm865v7rdx7p.streamlit.app/ (the first load takes
about twenty seconds while the instance wakes up).

## Quickstart

Python 3.12.

```bash
python3.12 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/streamlit run app/Home.py
```

The app needs `data/processed/` and nothing else. That cleaned parquet layer is committed, so
the app runs without the raw files. Notebooks additionally need
`pip install -r requirements-dev.txt`.

## Rebuilding the processed layer

The 11 source CSVs are committed in `data/raw/` (synthetic data, no PHI). To rebuild, run:

```bash
.venv/bin/python scripts/build.py
```

That rewrites every parquet frame in `data/processed/` plus `issues.csv`, the pipeline's decision
log. Nothing in the pipeline modifies a raw file.

## Layout

```
data/raw/*.csv                      11 source tables, committed (synthetic, no PHI)
   │
   │   src/ssc_coh/raw.py           reads them, changes nothing
   ▼
src/ssc_coh/                        the pipeline
   config.py                        paths, id and date column per table, aliases, fixed thresholds
   clean.py                         every quality rule, each one written to the issue log
   features.py                      one row per patient
   │
   ▼
scripts/build.py                    runs clean, then features, writes everything below
   │
   ▼
data/processed/                     19 parquet frames plus issues.csv, committed
   │
   ├── app/                         Streamlit: common.py, Home.py, pages/ (five pages); reads parquet only
   ├── notebooks/01_disease_patterns.ipynb   disease patterns tested on the cleaned layer
   └── REPORT.md                    the write-up; its section 3 is the register the issue log follows

notebooks/00_first_look.ipynb       reads data/raw directly: the table-by-table first look
                                    whose findings became the rules in config.py and clean.py
```

Notebook 00 sits beside the raw data. Everything else sits downstream of `scripts/build.py` and
never touches the raw files.

### Directories

```
data/raw/         11 source CSVs, committed (synthetic, no PHI)
data/processed/   19 parquet frames plus issues.csv, committed
src/ssc_coh/      config.py, raw.py, clean.py, features.py
scripts/          build.py
app/              common.py, Home.py, pages/ (five pages)
notebooks/        00_first_look.ipynb (raw), 01_disease_patterns.ipynb (processed)
tests/            test_smoke.py
REPORT.md         the write-up
```

## Findings and decisions

Everything the audit found, what was done about each issue and why, the disease patterns that
came out of the cleaned data, and the AI workflow behind it are in [REPORT.md](REPORT.md).
