# SSc cohort: what the data is, what is wrong with it, and what I did about it

A synthetic systemic sclerosis registry of 1,500 patients in 11 source tables,
profiled in `notebooks/00_first_look.ipynb`, cleaned by `src/ssc_coh` and `scripts/build.py`,
tested in `notebooks/01_disease_patterns.ipynb`, and browsable in the Streamlit app under `app/`.
Every number below sits under an output in those notebooks or in `data/processed/issues.csv`.

## 1. What the data is

Systemic sclerosis is a rare autoimmune disease in which the immune system attacks the body's own
connective tissue and the repair response lays down too much collagen. Skin thickens and hardens,
and the same scarring reaches internal organs, most dangerously the lungs. There is no cure, so
patients are followed for life, which is why nearly every table here is longitudinal. Two
subtypes divide the disease by how far the skin thickening spreads. Limited cutaneous disease
(lcSSc) stays on the hands, forearms and face and moves slowly. Diffuse cutaneous disease (dcSSc)
reaches the trunk and thighs, moves faster, and takes the organs with it. The cohort is 900
limited and 598 diffuse after cleaning.

Three blood tests carry most of the prognosis. Anti-centromere antibody marks limited disease and
a low risk of lung scarring. Anti-topoisomerase I, usually called Scl-70, marks diffuse disease
and interstitial lung disease (ILD), which is fibrosis of the lung tissue itself and the leading
cause of death in this population. Anti-RNA polymerase III marks diffuse disease and kidney
crisis. Skin severity is scored with the modified Rodnan skin score (mRSS): a clinician pinches
the skin at 17 sites, scores each 0 to 3, and adds them up, so the total runs 0 to 51 and higher
is worse. Lung function comes as three percentages of what a healthy person of the same age, sex
and height would produce. FVC is the total volume the patient can blow out, FEV1 the volume in
the first second, and DLCO how well gas crosses from air into blood. Fibrosis stiffens the lung,
so FVC and DLCO fall; 100 means exactly as expected and 80 is the lower limit of normal.

The tables open in the order a patient moves through the clinic.


| table           | one row per        | rows   | patients                             |
| --------------- | ------------------ | ------ | ------------------------------------ |
| `demographics`  | registered patient | 1,500  | 1,500                                |
| `ssc_subtype`   | registered patient | 1,500  | 1,500                                |
| `vitals`        | measurement        | 40,250 | 1,503 (1,500 registered, 3 controls) |
| `lab_report`    | test result        | 30,941 | 1,313                                |
| `antibodies`    | test result        | 3,511  | 1,431                                |
| `mrss`          | skin-score visit   | 3,113  | 824 (821 registered, 3 controls)     |
| `pft`           | measurement        | 5,069  | 1,174 (1,171 registered, 3 controls) |
| `medications`   | prescription       | 3,273  | 1,327                                |
| `bal`           | lavage procedure   | 1,170  | 820                                  |
| `libraries`     | sample             | 918    | 702 (698 registered, 4 controls)     |
| `skin_biopsies` | biopsy             | 215    | 150                                  |

The research tables cover nested subsets by design: only part of the cohort goes to bronchoscopy,
and a smaller part again to biopsy. In `bal`, bronchoalveolar lavage, saline is washed into one
lung segment and suctioned back out with the immune cells living there; those cells become the
RNA-seq samples in `libraries`.

## 2. How I looked

I ran the same routine on all 11 tables before drawing any conclusion: shape and `head()` to see
what a row is, `describe()` for ranges, `value_counts()` on every categorical and on every
measure in the long tables, and a histogram of each numeric column. `value_counts()` on the
continuous variables is what found the template defaults. A measured quantity should not have a
most frequent value: weight 160.6 appears in 870 rows against 39 for the runner-up, and systolic
124 fills 25% of the blood-pressure readings.

Single tables cannot settle unit or timeline questions, so four cross-table checks carried the
audit. Registration weight against the BMI the nurses recorded in `vitals` decided which patients
were weighed in kilograms, including nine whose registration weights of 84 to 104 look like
ordinary pounds. Prescription dates against `diagnosis_date` separated the drugs that routinely
precede a diagnosis from the immunosuppressants that cannot, and the same check ran against
lavages and biopsies. The link from `bal` to `libraries` confirmed the research chain: 367 of 372
processed samples have a lavage on the same patient 1 to 4 days earlier, and the five without are
the controls, who have no lavage rows.

Every claim in the notebooks sits directly under the output that prints its number. That rule
is why the audit could be re-checked line by line later, and why six wrong numbers surfaced
instead of being carried forward.

## 3. What I found and what I did about it

Forty-two issues across the ten stations, grouped here by treatment rather than by table,
because the cleaning code is written by treatment.

### Identity and structure, fixed deterministically

Each of these becomes a rule in `clean.py`. The change is mechanical and reversible.


| #  | table        | what                                                                                                                                           | action                                                                 |
| -- | ------------ | ---------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------- |
| 5  | all tables   | patient id under four column names, one of them with a space                                                                                   | rename to one id column                                                |
| 34 | libraries    | 24 column names carrying spaces, units and `?`                                                                                                  | rename to snake_case before anything else                              |
| 4  | demographics | two duplicate registrations (same name and birth date, different ids), so one person is counted twice and their longitudinal records are split | one canonical id per person, chosen by which id has more clinical rows |
| 9  | vitals       | 3 `SSC_NORM_*` ids not in the registry, which must not be counted as SSc patients                                                               | quarantine                                                             |
| 20 | mrss         | the same 3 `SSC_NORM_*` control ids                                                                                                             | quarantine                                                             |
| 23 | pft          | the same 3 `SSC_NORM_*` control ids (9 rows)                                                                                                    | quarantine                                                             |
| 33 | libraries    | 4 `SSC_NORM_*` control ids (5 rows); `SSC_NORM_0104` appears nowhere else                                                                       | quarantine                                                             |
| 14 | lab_report   | `DIFFERENTIAL TYPE` puts the text "Automated" in the value column (232 rows), so the column holds metadata that cannot be averaged             | split into its own field                                               |
| 27 | medications  | one drug under three names (`mycophenolate mofetil`, `CellCept`, `MMF`; 592 rows under the three names, of which 381 are alias rows), which splits any drug count three ways | map the 381 alias rows to the generic name                             |
| 17 | antibodies   | 4 rows are exact duplicates (2 pairs), which double counts them                                                                                | drop one copy of each                                                  |

### Units and sentinels, fixed only with evidence and logged

Each writes a row to `issues.csv` with the number of values changed.


| #  | table                | what                                                                                                                                                  | action                                                                 |
| -- | -------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------- |
| 1  | demographics         | height in two units, inches and cm, split at 100; every BMI and every plot on height is wrong until this is fixed                                     | convert to one unit                                                    |
| 2  | demographics         | a few weights read as kg among pounds (rows below 80)                                                                                                 | convert; the exact list comes from station 3                           |
| 11 | vitals, demographics | 14 patients' weights are kg (recorded BMI below 15), and their stored BMI was computed from the kg number with the imperial formula                   | convert the weights, blank the derived BMI                             |
| 3  | demographics         | one weight (22.5) impossible in any unit, so it cannot be repaired                                                                                    | missing and logged                                                     |
| 10 | vitals               | weight 160.6 is a template default filling 870 rows, the entire weight history of 194 patients, so any weight mean or trend is wrong until it goes    | missing and logged                                                     |
| 15 | lab_report           | placeholder and impossible values: 999, 9999, negatives, RBC = 0 (12 rows)                                                                            | missing and logged; legitimate zeros in eosinophils and basophils kept |
| 29 | medications          | 45 rows with impossible doses (zero, negative, 999 tablets, 250 g, 10000 mg every hour, mL tablets) where the pipeline caught only `999 tablets daily` | dose set to missing and logged; widen the rule                         |

### Structural artifacts, kept, flagged, excluded from the analyses they would distort

Nothing is deleted; the app states the exclusion.


| #  | table  | what                                                                                                                                                       | action                                                                                        |
| -- | ------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------- |
| 25 | pft    | FEV1 is FVC plus noise: ratio pinned at 1.00, no obstructive tail, no link to subtype or fibrosis, and the pattern drifts by year                          | exclude FEV1 and FEV1/FVC as uninformative; flag the 904 visits (2,655 PFT rows) where FEV1 came out above FVC |
| 12 | vitals | SBP 124, DBP 77 and pulse 76 fill 25%, 25% and 16% of readings, and 726 of the 1,061 patients with at least four readings have near-constant series, so blood pressure is not a longitudinal signal | report, no silent fix                                                                         |
| 26 | pft    | values bounded at exactly 40 and 130, with DLCO piling up at 40                                                                                            | keep as generator clipping, note it in the dictionary                                         |

### Ambiguous and timeline issues, flagged only, never altered

The data cannot say which value is right, so each gets a flag column the app surfaces.


| #  | table         | what                                                                                                                                                           | action                                                    |
| -- | ------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------- |
| 7  | ssc_subtype   | Raynaud date after non-Raynaud date for 308 patients at the same rate in both subtypes, so any "time from onset" variable is unreliable for 20% of the cohort  | flag; never swap the dates                                |
| 8  | ssc_subtype   | diagnosis before the first non-Raynaud symptom, 51 patients, which is clinically possible                                                                      | flag only                                                 |
| 13 | vitals        | systolic below diastolic, `subject_4632`, 2 visits and the 4 BP rows on them, reading as a field swap                                                           | flag only                                                 |
| 18 | antibodies    | 13 patients flip between positive and negative on repeat testing, although antibody status is treated as fixed in SSc                                          | keep and flag; use the latest result and report the count |
| 21 | mrss          | two scores for one patient on one day (`subject_2545`: 20 and 14), with no basis to pick                                                                       | keep both and flag; never average them silently           |
| 30 | medications   | 14 SSc-specific immunosuppressant prescriptions dated 5 to 80 days before the diagnosis, likely because the diagnosis date was entered after treatment started | flag only                                                 |
| 31 | bal           | 104 procedures (9%) dated before the diagnosis, up to 3.3 years earlier, which is too far to be an entry lag                                                   | flag; no fix possible                                     |
| 35 | libraries     | status contradiction: 13 samples marked `Sequence? = No` but `Complete? = Complete` (23 "No" samples with library work)                                         | keep both values and flag                                 |
| 36 | libraries     | 88 of 188 sequenced samples have no RIN although the QC date and concentration exist, and the RIN column is partly hand-edited                                 | keep missing and flag; a question for the wet lab         |
| 37 | libraries     | 4 sample ids used twice, on different processing dates or batches                                                                                              | flag for lab adjudication                                 |
| 38 | libraries     | 4 sequencing batches hold a single sample (09, 10, 12, 15), so they cannot be batch-corrected                                                                  | flag                                                      |
| 39 | libraries     | RNA concentration maximum 49,785 pg/ul against a median of 738, a possible unit slip                                                                           | flag and check                                            |
| 41 | skin_biopsies | the file extension disagrees with `image_format` in 132 of 215 rows, and `.svs` is never declared, so no image pipeline can trust the format column             | keep both fields and flag; pathology adjudicates          |
| 42 | skin_biopsies | 17 biopsies dated before the diagnosis                                                                                                                         | flag                                                      |

### Not errors, documented


| #  | table       | what                                                                                                     | action                                                        |
| -- | ----------- | -------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------- |
| 6  | ssc_subtype | `other_dx` is free text and missing for 37% of patients, so a missing value means unrecorded, not absent | parse into per-condition flags; missing stays unknown         |
| 16 | lab_report  | test names from a lab interface (`MONOCYTES, BODY FLUID`, `NM BKR ABSOLUTE NEUTROPHIL`)                  | keep as is, document in the dictionary                        |
| 19 | antibodies  | 78 borderline and indeterminate results, which are legitimate assay outcomes                             | keep; out of positive-against-negative comparisons by default |
| 24 | pft         | DLCO missing in 49% of its rows (794 of 1,623) while FVC and FEV1 are complete                           | keep missing, no imputation; display the missingness          |
| 28 | medications | 20 rows with no drug name, 85 with no dose, but the dates still matter                                   | keep the rows; exclude from drug-level counts                 |
| 32 | bal         | the comment column is a 6-item pick-list and site counts are flat across 5 sites                         | document                                                      |

### Corrections to earlier work


| #  | what                                                                                                                                                         | action                                                                                  |
| -- | ------------------------------------------------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------------------------- |
| 22 | the mRSS rater spread of 7.7 to 25.8 is case mix, not a rater effect: within one subtype it collapses to 7.1 to 9.5 for limited and 22.5 to 29.5 for diffuse | correct the claim; my earlier audit document reported a rater problem that is not there |
| 40 | the libraries `comment` is the lab's own note, identical to `bal_comment` in only 7 of 121 cases                                                              | correct the data-dictionary note; keep the column in the cleaned layer                  |

### Three principles

Clear errors become missing and get a row in the log. I never repair one by guessing, which is
why the 22.5 weight, the 9999 lab values and the 250 g doses are blank rather than plausible.

Structural artifacts stay in the data, carry a flag, and are excluded from the analyses they
would distort. FEV1 and the blood-pressure trend are real values that mean nothing, so deleting
them would hide the problem and using them would launder it.

Ambiguous values are flagged, surfaced in the app, and never altered: onset order, antibody
flips, the library status contradictions, the timeline paradoxes. The data cannot say which
reading is right, and neither can I.

### What the table-by-table walk changed compared with my first audit

The 160.6 weight does not patch missing visits. It replaces the whole weight history of 194
patients, who therefore have no weight at all.

The Raynaud-after-symptom gap is the left tail of one smooth bell, identical in both subtypes,
with only 32 patients within 30 days of zero: an unconstrained generator draw, not swapped
dates.

FEV1 is excluded as uninformative, because it is FVC plus noise with no physiology in the ratio.
My first draft called it impossible, which was wrong.

The mRSS rater spread is case mix. Dr. Fischer at 7.7 scored no diffuse patients and Dr. Rossi at
25.8 scored 88% diffuse; inside one subtype the spread is a few points.

Impossible doses are 45 rows in seven patterns, not the 7 rows of `999 tablets daily` that the
pipeline caught.

The libraries `comment` is the lab's own note, not a copy of the lavage comment, so the cleaned
`libraries` frame keeps it beside the two TapeStation comment columns.

## 4. The cleaning pipeline

`src/ssc_coh` holds four files. `config.py` carries the fixed decisions: the id column per table,
the date columns, the control prefix, the medication aliases and the LIMS column renames.
`raw.py` reads the 11 CSVs and does nothing else, so the raw files are never modified. `clean.py`
implements every rule above, numbered as in section 3. `features.py` collapses the long tables
into one row per patient. `scripts/build.py` runs all of it and writes `data/processed/`.

The decision log in `data/processed/issues.csv` now has 49 entries, one per rule that changed or
flagged data, each with the table, the issue, the decision, an `n_affected` count and a detail
field naming the patients or values involved. `n_affected` counts rows of the table named in the
row. Where a rule works at the level of a visit or a patient-test pair, the flag still lands on
every row of that visit or pair, so the detail field states both counts, for example
`2 visits, 4 BP rows`.

Nothing is repaired in place without a flag. These columns mark what stayed.


| flag column                          | table                      | what it marks                                                            | flagged                                                                    |
| ------------------------------------ | -------------------------- | ------------------------------------------------------------------------ | -------------------------------------------------------------------------- |
| `merged_registration_flag`           | subjects, ssc_subtype      | the kept row of a patient who was registered twice                       | 2 rows in each table, 2 patients                                           |
| `diagnosis_before_symptom_flag`      | ssc_subtype                | diagnosis recorded before the first non-Raynaud symptom                  | 51 rows, one row per patient                                               |
| `onset_order_flag`                   | ssc_subtype                | Raynaud date after the non-Raynaud date                                  | 307 rows (308 before the merge)                                            |
| `bp_inverted_flag`                   | vitals                     | the BP rows of a visit whose systolic reading is below the diastolic     | 4 rows over 2 visits                                                       |
| `result_flip_flag`                   | antibodies                 | every row of a patient-test pair holding both positive and negative      | 31 rows over 15 patient-test pairs, 15 patients (13 before the merge)      |
| `fev1_gt_fvc_flag`                   | pft                        | every PFT row of a visit where FEV1 came out above FVC                   | 2,655 rows over 904 visits                                                 |
| `same_day_conflict`                  | mrss                       | two skin scores for one patient on one day                               | 2 rows, 1 patient                                                          |
| `before_diagnosis_flag`              | medications, bal, biopsies | an event dated before the diagnosis                                      | 70, 104, 17 rows                                                           |
| `ssc_specific_before_diagnosis_flag` | medications                | an SSc-specific immunosuppressant among those                            | 14 rows                                                                    |
| `volume_flag`                        | bal                        | recovered volume above instilled volume                                  | 0 rows; the rule holds on all 1,170 procedures                             |
| `format_mismatch`                    | biopsies                   | the file extension disagrees with `image_format`                         | 132 of 215 rows                                                            |
| `status_contradiction`               | libraries                  | `Sequence? = No` with library work recorded                              | 23 rows                                                                    |
| `rin_missing_for_sequenced`          | libraries                  | a sequenced sample with no RIN although QC ran                           | 88 rows, of the 184 sequenced samples left after quarantine                |
| `sample_id_reused`                   | libraries                  | a sample id that appears more than once                                  | 8 rows, 4 ids                                                              |
| `single_sample_batch`                | libraries                  | a sequencing batch holding one sample                                    | 4 rows in the export, all four of them control samples, so 0 rows remain in the cleaned frame |
| `rna_concentration_outlier`          | libraries                  | a concentration more than 10 times the median                            | 41 rows                                                                    |

The dose rule reads the leading number and unit off the free-text dose string and rejects what
cannot be a prescription: any amount at or below zero, 10 g or more, 5,000 mg or more, 10 tablets
or more, anything dosed every hour, and milliliters of tablets. A dose with no leading number,
such as `PRN`, is left alone. It catches 45 rows across seven strings where the first version
caught 7.

The comorbidity flags `dx_ild`, `dx_gerd` and `dx_pah` are parsed out of the free-text `other_dx`
field and are three-state. Where `other_dx` is empty the flag is unknown rather than false, which
is 555 of the 1,498 patients. Every rate built on them therefore has two readings: among the 943
patients whose field was recorded, and among all 1,498 with an empty field read as no
comorbidity. Both appear side by side in the notebook and in the app.

Merging the two duplicate registrations was the one identity fix that cost information. Each pair
shares a name and a birth date but disagrees elsewhere. In `demographics`, `subject_7539`
disagrees with its twin on `case number` and `state`, and `subject_8994` on `case number` alone.
In `ssc_subtype`, `subject_7539` disagrees on `study_code`, `raynaud_date`, `nonraynaud_date`,
`nonraynaud_sx` and `diagnosis_date`, and `subject_8994` on those five plus `ssc_subtype` and
`other_dx`, so one of these patients is limited under one id and diffuse under the other. I keep
the id with more clinical rows, remap the longitudinal rows to it, drop the second registry row,
and log the disagreeing fields by name so the choice can be reversed.

## 5. Findings

I tested four patterns against the cleaned layer. The cohort lands where SSc epidemiology
predicts before any of them: 1,498 patients, 82.4% female, 900 limited against 598 diffuse, ages
36 to 82.

Antibodies track subtype in the directions the literature reports. Anti-centromere is positive in
47.1% of limited patients against 5.0% of diffuse, a ninefold gap over 914 tested patients
(chi-square 174.2, p = 8.8e-40). Scl-70 leans the other way, 52.3% diffuse against 36.9% limited,
and RNA polymerase III leans the same way weakly, 13.1% against 8.4% (p = 0.039). All three
directions coming out right is the evidence that the antibody table is joined to the right
patients.

Skin severity separates the subtypes almost completely. Latest mRSS has a median of 24 in diffuse
disease against 8 in limited, quartiles 21 to 28 against 5 to 12, so the two middle halves do not
overlap, and the mean across visits gives the same separation, 23.4 against 8.0. 821 of the 1,498
patients have a skin score at all, across 3,110 visits.

The lung triangle holds on all three edges. Among patients whose comorbidity field was recorded,
ILD is noted for 84.0% of Scl-70 positive patients against 57.2% of Scl-70 negative, for 37.1% of
anti-centromere positive against 74.0% of anti-centromere negative, and for 86.1% of diffuse
against 50.6% of limited. Reading an empty field as no comorbidity lowers every rate and leaves
the ordering untouched: 64.7% against 32.3% by Scl-70, 18.1% against 50.4% by anti-centromere,
63.2% against 28.3% by subtype. The chi-square is decisive under both denominators. The true rate
lies between the two columns, because an empty field can mean there was nothing to record or that
nobody recorded it, and no analysis here picks one. Lung function agrees: latest FVC median 78.5
in diffuse against 88.0 in limited, a gap of 9.5 points, and DLCO 74.0 against 77.9.

The fourth pattern is a negative result. Annual FVC decline is the standard endpoint in SSc
trials, so I fit a per-patient slope for the 108 patients with at least three FVC tests. The
slopes center just above zero, median 0.60 and mean 0.18% predicted per year, with 47 patients
(43.5%) falling rather than rising. The subtypes are indistinguishable (p = 0.96) and the whole
set does not differ from zero (p = 0.31). This cohort carries cross-sectional lung differences
and no progressive decline.

## 6. The app

Five pages, each a level further into the data, all reading the same cleaned parquet layer.

Data & Quality carries the dictionary for every table with dtypes and fill rates, the coverage
and linkage report, the issue log as a filterable download, and the quarantined control rows. The
design decision is that quarantine is a tab rather than a delete, because removing the four
`SSC_NORM_*` subjects would hide the linkage problem that found them.

Cohort shows who is in the registry (age, sex, disease duration) and the patterns that hold
across it. The design decision is that the ILD chart plots both denominators side by side rather
than picking one, since `other_dx` is unrecorded for 37% of patients.

Compare Groups crosses any grouping variable with any numeric measure or outcome rate, with a box
plot, a median and IQR table, and a Mann-Whitney, Kruskal-Wallis or chi-square test. The design
decision is that the tests are rank-based and the caption calls the p-values screening rather
than confirmation, since there is no multiplicity correction.

Discover Structure is the open-ended page: a correlation matrix in Spearman or Pearson, a PCA
projection colored by any grouping with its top loadings printed, and a free X-Y explorer. The
design decision is that PCA imputes missing values at the column median and says so, because
dropping incomplete rows would discard most of the cohort.

Patient draws one person's trajectory on a shared time axis: skin score, lung function, weight
and blood pressure in aligned panels with the three disease milestones as vertical lines across
all of them, then tabs down to medications, antibody history, labs and research samples. The
design decision is that shared axis, because the question the page exists to answer is whether
the skin and the lungs move together.

## 7. How I worked with AI tools, and how I checked the output

A first-pass weight-unit rule blanked 20 demographics weights where my audit had predicted 1.
Reconciling the two numbers showed the rule was working from BMI plausibility alone and erasing
genuinely obese patients. I redesigned it around the vitals BMI cross-check; the final rule
converts 14 patients from the vitals evidence, 0 from BMI plausibility alone, and sets 1 weight
to missing.

I had an independent sub-agent recompute every number in the two audit documents straight from
the raw files. 50 of 56 matched and 6 were corrected, among them the metric-subject counts, the
direction of the unit conversion, and DLCO missingness by subtype, which the draft had the wrong
way round.

I asked why FEV1 above FVC should be impossible. It is not. Both columns are % predicted with
different denominators, and in restrictive disease FEV1 % predicted at or above FVC % predicted
is common. The finding was reframed from impossible to uninformative, FVC plus noise.

Re-reading `clean.py` line by line showed the issue-log frame was materialized before the
quarantine step, so four entries never reached `issues.csv`. Fixing that exposed a second error:
control rows were also being counted as orphans, a double count.

The table-by-table walk in notebook 00 found 45 impossible dose rows in seven patterns where the
pipeline caught 7, showed the mRSS rater spread to be case mix, and showed the libraries
`comment` to be the lab's own note rather than a copy of the lavage comment.

The rebuilt layer showed that the two merged registrations disagree on their disease fields,
which is now logged with the field names and flagged on the kept row.

The rule that came out of this is that generated code gets checked against an independently
derived expectation before I trust it, and a mismatch is investigated rather than settled by
adjusting the expectation. The second is that every number in a document has to be traceable to
an output, which is what made the six corrections findable at all.

## 8. Limitations, tradeoffs and what I would do next

The PCA on the Discover page fills every missing feature with the column median. The latest DLCO
is missing for 63% of patients, so those patients sit at the cohort median on that axis and the
projection pulls them toward the center. Keeping only complete rows would leave 44 patients of
1,498, which is why I impute, but any tight group in the middle of that plot is partly an artifact
of the choice, and the caption on the page says so.

The group comparisons are rank tests with no multiplicity correction. The Compare page offers 150
grouping-by-measure combinations, so a p-value there is a screening number and not a confirmation.
I show the group sizes, medians and IQRs beside it so the size of the difference is visible and
not only its p-value.

Merging the two duplicate registrations cost information. Each pair shares a name and a birth date
and disagrees elsewhere, on subtype for one of the two. I keep the registry row of the id with more
clinical rows, drop the other, and log the disagreeing field names so the merge can be reversed,
but the analysis runs on one of two readings of those two patients.

There is no schema validation, and the only automated check is `tests/test_smoke.py`: the raw
loader returns 11 tables, the build returns 19 frames, the issue log has its five columns and is
not empty, and every measurement frame carries `subject_id`. Nothing tests dtypes, value ranges or
the cleaning rules themselves.

What I would do next, in this order. Declare each raw table's schema with pandera and fail the
build on a violation, so a new export cannot change a column type without being noticed. Require a
minimum follow-up span before fitting an FVC slope: the 108 patients I fit have a median of 3 tests
over a median span of 1.00 year, and 72 of them span less than a year, which is too short to call a
trend. Carry the RNA-seq batch into any analysis of the libraries, since it is a source of
technical variation that nothing here accounts for, and four batches held a single sample before
quarantine. Deploy the app to Streamlit Cloud so a reviewer does not have to install anything.
Version the issue log with the build that produced it, so two runs can be compared rule by rule
instead of file by file.
